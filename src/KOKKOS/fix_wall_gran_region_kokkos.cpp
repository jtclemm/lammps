/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_wall_gran_region_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory_kokkos.h"
#include "region_block_kokkos.h"
#include "region_sphere_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

enum { NOSTYLE = -1, XPLANE = 0, YPLANE = 1, ZPLANE = 2, ZCYLINDER, REGION };

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixWallGranRegionKokkos<DeviceType>::FixWallGranRegionKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixWallGranOld(lmp, narg, arg), region(nullptr), ncontact(nullptr), walls(nullptr),
  history_many(nullptr)
{
  if (wallstyle != REGION)
    error->all(FLERR, "Fix {} requires the region keyword", style);
  if (pairstyle != HOOKE_HISTORY)
    error->all(FLERR, "Fix {} currently supports only the hooke/history interaction", style);

  // the per-atom array of the "contacts" keyword is not maintained here

  if (peratom_flag)
    error->all(FLERR, "Fix {} does not yet support the contacts keyword", style);

  restart_global = 1;
  motion_resetflag = 0;

  region = domain->get_region_by_id(idregion);
  if (!region) error->all(FLERR, "Region {} for fix {} does not exist", idregion, style);
  tmax = region->tmax;

  kokkosable = 1;
  exchange_comm_device = sort_device = 1;
  stride = 1 + tmax * (1 + size_history);
  maxexchange = stride;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | OMEGA_MASK | TORQUE_MASK | RADIUS_MASK |
    RMASS_MASK | MASK_MASK;
  datamask_modify = F_MASK | TORQUE_MASK;

  // replace the single-contact history of the base class
  // the base class already registered with the Atom class

  memory->destroy(history_one);
  history_one = nullptr;

  grow_arrays(atom->nmax);

  // initialize history as if particle is not touching region

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) ncontact[i] = 0;
  k_ncontact.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixWallGranRegionKokkos<DeviceType>::~FixWallGranRegionKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_ncontact, ncontact);
  memoryKK->destroy_kokkos(k_walls, walls);
  memoryKK->destroy_kokkos(k_history_many, history_many);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixWallGranRegionKokkos<DeviceType>::region_supported() const
{
  if (RegBlockKokkos<DeviceType>::cast(region)) return 1;
  if (RegSphereKokkos<DeviceType>::cast(region)) return 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranRegionKokkos<DeviceType>::init()
{
  FixWallGranOld::init();

  if (fix_rigid) error->all(FLERR, "Fix {} is not yet compatible with rigid bodies", style);

  auto *newregion = domain->get_region_by_id(idregion);
  if (!newregion) error->all(FLERR, "Region {} for fix {} does not exist", idregion, style);

  // check if region properties changed between runs
  // reset if restart info was inconsistent

  if (newregion != region) {
    if (newregion->tmax != tmax)
      error->all(FLERR, "Region {} for fix {} was replaced by a region with a different "
                 "number of walls", idregion, style);
    region = newregion;
    if (comm->me == 0)
      error->warning(FLERR, "Region properties for region {} changed between runs, "
                     "resetting its motion", idregion);
    region->reset_vel();
  }

  if (motion_resetflag) {
    if (comm->me == 0)
      error->warning(FLERR, "Region properties for region {} are inconsistent with restart "
                     "file, resetting its motion", idregion);
    region->reset_vel();
  }

  if (!region_supported())
    error->all(FLERR, "Fix {} requires region style block/kk or sphere/kk", style);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranRegionKokkos<DeviceType>::post_force(int /*vflag*/)
{
  // do not update shear history during setup

  history_update = 1;
  if (update->setupflag) history_update = 0;

  coeffs.kn = static_cast<KK_FLOAT>(kn);
  coeffs.kt = static_cast<KK_FLOAT>(kt);
  coeffs.gamman = static_cast<KK_FLOAT>(gamman);
  coeffs.gammat = static_cast<KK_FLOAT>(gammat);
  coeffs.xmu = static_cast<KK_FLOAT>(xmu);
  coeffs.dt = static_cast<KK_FLOAT>(dt);
  dt_region = update->dt;

  // the KOKKOS regions store the periodic box in prematch(), so it is called
  // every step; set_velocity() also updates the previous region position

  regiondynamic = region->dynamic_check();
  region->prematch();
  if (regiondynamic) region->set_velocity();

  atomKK->sync(execution_space, datamask_read);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  d_omega = atomKK->k_omega.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  d_radius = atomKK->k_radius.view<DeviceType>();
  int nlocal = atom->nlocal;

  k_ncontact.template sync<DeviceType>();
  k_walls.template sync<DeviceType>();
  k_history_many.template sync<DeviceType>();
  d_ncontact = k_ncontact.template view<DeviceType>();
  d_walls = k_walls.template view<DeviceType>();
  d_history_many = k_history_many.template view<DeviceType>();

  // nerror = # of particles touching more walls than the history can hold

  int nerror = 0;

  // the functor copies this fix and the region; their copies must not free memory

  copymode = 1;
  region->copymode = 1;

  if (auto *regionKK = RegBlockKokkos<DeviceType>::cast(region)) {
    FixWallGranRegionKokkosFunctor<DeviceType, RegBlockKokkos<DeviceType>> functor(this, regionKK);
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0, nlocal), functor, nerror);
  } else if (auto *regionKK = RegSphereKokkos<DeviceType>::cast(region)) {
    FixWallGranRegionKokkosFunctor<DeviceType, RegSphereKokkos<DeviceType>> functor(this, regionKK);
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0, nlocal), functor, nerror);
  }

  region->copymode = 0;
  copymode = 0;

  atomKK->modified(execution_space, datamask_modify);
  k_ncontact.template modify<DeviceType>();
  k_walls.template modify<DeviceType>();
  k_history_many.template modify<DeviceType>();

  if (nerror) error->one(FLERR, "Too many wall/gran/region contacts for one particle");
}

/* ----------------------------------------------------------------------
   velocity of the region wall at the contact point
   mirrors Region::velocity_contact(), see there
------------------------------------------------------------------------- */

template<class DeviceType>
template<class RegionKK>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallGranRegionKokkos<DeviceType>::velocity_contact(const RegionKK &regionKK,
                                                           const double *xi,
                                                           const Region::Contact &c,
                                                           double *vwall) const
{
  double xc[3];
  xc[0] = xi[0] - c.delx;
  xc[1] = xi[1] - c.dely;
  xc[2] = xi[2] - c.delz;

  vwall[0] = vwall[1] = vwall[2] = 0.0;

  if (regionKK.moveflag) {
    vwall[0] = regionKK.v[0];
    vwall[1] = regionKK.v[1];
    vwall[2] = regionKK.v[2];
  }
  if (regionKK.rotateflag) {
    const double *w = regionKK.omega;
    const double *rp = regionKK.rpoint;
    vwall[0] += w[1] * (xc[2] - rp[2]) - w[2] * (xc[1] - rp[1]);
    vwall[1] += w[2] * (xc[0] - rp[0]) - w[0] * (xc[2] - rp[2]);
    vwall[2] += w[0] * (xc[1] - rp[1]) - w[1] * (xc[0] - rp[0]);
  }

  if (regionKK.varshape && c.varflag)
    regionKK.velocity_contact_shape_kokkos(vwall, xc, dt_region);
}

/* ----------------------------------------------------------------------
   interaction of particle I with all walls of the region it touches
   history is matched to contacts via the id of the wall (iwall)
------------------------------------------------------------------------- */

template<class DeviceType>
template<class RegionKK>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallGranRegionKokkos<DeviceType>::wall_particle(const RegionKK &regionKK, const int i,
                                                        int &nerror) const
{
  if (!(mask(i) & groupbit)) return;

  const double xi[3] = {static_cast<double>(x(i,0)), static_cast<double>(x(i,1)),
                        static_cast<double>(x(i,2))};
  if (!regionKK.match_kokkos(xi[0], xi[1], xi[2])) return;

  const KK_FLOAT radius = d_radius(i);

  // the contact list lives on the stack, so that concurrently running
  // threads do not overwrite each other's contacts

  Region::Contact contact[RegionKK::MAXCONTACT];
  const int nc = regionKK.surface_kokkos(xi[0], xi[1], xi[2], static_cast<double>(radius),
                                         contact);
  if (nc > tmax) {
    nerror++;
    return;
  }

  // history maintenance
  // update ncontact, walls, history_many for particle I
  //   to reflect new and persistent history values
  // slot[ic] = index into the history of particle I for contact ic

  int slot[RegionKK::MAXCONTACT];
  int n = d_ncontact(i);

  if (nc == 0) {
    d_ncontact(i) = 0;
    return;
  }

  if ((nc == 1) && (n == 0)) {
    slot[0] = 0;
    d_walls(i,0) = contact[0].iwall;
    for (int k = 0; k < size_history; k++) d_history_many(i,0,k) = 0.0;
    n = 1;
  } else if ((nc == 1) && (n == 1) && (contact[0].iwall == d_walls(i,0))) {
    slot[0] = 0;
  } else {

    // delete old contacts that are not in the new contact list
    //   by copying the last contact over them

    int iold = 0;
    while (iold < n) {
      int m;
      for (m = 0; m < nc; m++)
        if (contact[m].iwall == d_walls(i,iold)) break;
      if (m >= nc) {
        const int ilast = n - 1;
        for (int k = 0; k < size_history; k++)
          d_history_many(i,iold,k) = d_history_many(i,ilast,k);
        d_walls(i,iold) = d_walls(i,ilast);
        n--;
      } else iold++;
    }

    // match new contacts to the remaining old ones,
    //   add the others with zeroed history

    const int nold = n;
    for (int inew = 0; inew < nc; inew++) {
      const int iwall = contact[inew].iwall;
      int m;
      for (m = 0; m < nold; m++)
        if (d_walls(i,m) == iwall) break;
      if (m < nold) slot[inew] = m;
      else {
        slot[inew] = n;
        for (int k = 0; k < size_history; k++) d_history_many(i,n,k) = 0.0;
        d_walls(i,n) = iwall;
        n++;
      }
    }
  }

  d_ncontact(i) = n;

  // process current contacts

  const KK_FLOAT vi[3] = {v(i,0), v(i,1), v(i,2)};
  const KK_FLOAT wi[3] = {d_omega(i,0), d_omega(i,1), d_omega(i,2)};
  const KK_FLOAT meff = rmass(i);

  for (int ic = 0; ic < nc; ic++) {
    const KK_FLOAT rsq = static_cast<KK_FLOAT>(contact[ic].r * contact[ic].r);
    const KK_FLOAT dx = static_cast<KK_FLOAT>(contact[ic].delx);
    const KK_FLOAT dy = static_cast<KK_FLOAT>(contact[ic].dely);
    const KK_FLOAT dz = static_cast<KK_FLOAT>(contact[ic].delz);

    double vwall[3] = {0.0, 0.0, 0.0};
    if (regiondynamic) velocity_contact(regionKK, xi, contact[ic], vwall);
    const KK_FLOAT vw[3] = {static_cast<KK_FLOAT>(vwall[0]), static_cast<KK_FLOAT>(vwall[1]),
                            static_cast<KK_FLOAT>(vwall[2])};

    const int s = slot[ic];
    KK_FLOAT fi[3], ti[3];
    WallGranKokkos::hooke_history(dx, dy, dz, rsq, vw, vi, wi, radius, meff,
                                  d_history_many(i,s,0), d_history_many(i,s,1),
                                  d_history_many(i,s,2), history_update, coeffs, fi, ti);

    f(i,0) += static_cast<KK_ACC_FLOAT>(fi[0]);
    f(i,1) += static_cast<KK_ACC_FLOAT>(fi[1]);
    f(i,2) += static_cast<KK_ACC_FLOAT>(fi[2]);
    torque(i,0) += static_cast<KK_ACC_FLOAT>(ti[0]);
    torque(i,1) += static_cast<KK_ACC_FLOAT>(ti[1]);
    torque(i,2) += static_cast<KK_ACC_FLOAT>(ti[2]);
  }
}

/* ----------------------------------------------------------------------
   make the host side of the history current, or mark it as modified
------------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranRegionKokkos<DeviceType>::sync_history_host()
{
  k_ncontact.sync_host();
  k_walls.sync_host();
  k_history_many.sync_host();
}

template<class DeviceType>
void FixWallGranRegionKokkos<DeviceType>::modify_history_host()
{
  k_ncontact.modify_host();
  k_walls.modify_host();
  k_history_many.modify_host();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double FixWallGranRegionKokkos<DeviceType>::memory_usage()
{
  double nmax = atom->nmax;
  double bytes = 0.0;
  bytes += nmax * sizeof(int);                                  // ncontact
  bytes += nmax * tmax * sizeof(int);                           // walls
  bytes += nmax * tmax * size_history * sizeof(double);        // history_many
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranRegionKokkos<DeviceType>::grow_arrays(int nmax)
{
  // force reallocation on host

  sync_history_host();
  modify_history_host();

  memoryKK->grow_kokkos(k_ncontact, ncontact, nmax, "wall/gran/region/kk:ncontact");
  memoryKK->grow_kokkos(k_walls, walls, nmax, tmax, "wall/gran/region/kk:walls");
  memoryKK->grow_kokkos(k_history_many, history_many, nmax, tmax, size_history,
                        "wall/gran/region/kk:history_many");

  d_ncontact = k_ncontact.template view<DeviceType>();
  d_walls = k_walls.template view<DeviceType>();
  d_history_many = k_history_many.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranRegionKokkos<DeviceType>::copy_arrays(int i, int j, int /*delflag*/)
{
  sync_history_host();

  int n = ncontact[i];
  for (int iwall = 0; iwall < n; iwall++) {
    walls[j][iwall] = walls[i][iwall];
    for (int m = 0; m < size_history; m++) history_many[j][iwall][m] = history_many[i][iwall][m];
  }
  ncontact[j] = ncontact[i];

  modify_history_host();
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranRegionKokkos<DeviceType>::set_arrays(int i)
{
  k_ncontact.sync_host();
  ncontact[i] = 0;
  k_ncontact.modify_host();
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixWallGranRegionKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  sync_history_host();

  int n = 0;
  int count = ncontact[i];
  buf[n++] = ubuf(count).d;
  for (int iwall = 0; iwall < count; iwall++) {
    buf[n++] = ubuf(walls[i][iwall]).d;
    for (int m = 0; m < size_history; m++) buf[n++] = history_many[i][iwall][m];
  }
  return n;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange
------------------------------------------------------------------------- */

template<class DeviceType>
int FixWallGranRegionKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  sync_history_host();

  int n = 0;
  int count = ncontact[nlocal] = (int) ubuf(buf[n++]).i;
  for (int iwall = 0; iwall < count; iwall++) {
    walls[nlocal][iwall] = (int) ubuf(buf[n++]).i;
    for (int m = 0; m < size_history; m++) history_many[nlocal][iwall][m] = buf[n++];
  }

  modify_history_host();
  return n;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

template<class DeviceType>
int FixWallGranRegionKokkos<DeviceType>::pack_restart(int i, double *buf)
{
  sync_history_host();

  int n = 1;
  int count = ncontact[i];

  buf[n++] = ubuf(count).d;
  for (int iwall = 0; iwall < count; iwall++) {
    buf[n++] = ubuf(walls[i][iwall]).d;
    for (int m = 0; m < size_history; m++) buf[n++] = history_many[i][iwall][m];
  }

  // pack buf[0] this way because other fixes unpack it

  buf[0] = n;
  return n;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranRegionKokkos<DeviceType>::unpack_restart(int nlocal, int nth)
{
  sync_history_host();

  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int>(extra[nlocal][m]);
  m++;

  int count = ncontact[nlocal] = (int) ubuf(extra[nlocal][m++]).i;
  for (int iwall = 0; iwall < count; iwall++) {
    walls[nlocal][iwall] = (int) ubuf(extra[nlocal][m++]).i;
    for (int k = 0; k < size_history; k++) history_many[nlocal][iwall][k] = extra[nlocal][m++];
  }

  modify_history_host();
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

template<class DeviceType>
int FixWallGranRegionKokkos<DeviceType>::maxsize_restart()
{
  return 2 + tmax * (size_history + 1);
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

template<class DeviceType>
int FixWallGranRegionKokkos<DeviceType>::size_restart(int nlocal)
{
  k_ncontact.sync_host();
  return 2 + ncontact[nlocal] * (size_history + 1);
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranRegionKokkos<DeviceType>::write_restart(FILE *fp)
{
  if (comm->me) return;
  int len = 0;
  region->length_restart_string(len);
  fwrite(&len, sizeof(int), 1, fp);
  region->write_restart(fp);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranRegionKokkos<DeviceType>::restart(char *buf)
{
  int n = 0;
  if (!region->restart(buf, n)) motion_resetflag = 1;
}

/* ----------------------------------------------------------------------
   sort local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranRegionKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  // always sort on the device

  k_ncontact.sync_device();
  k_walls.sync_device();
  k_history_many.sync_device();

  Sorter.sort(LMPDeviceType(), k_ncontact.view_device());
  Sorter.sort(LMPDeviceType(), k_walls.view_device());
  Sorter.sort(LMPDeviceType(), k_history_many.view_device());

  k_ncontact.modify_device();
  k_walls.modify_device();
  k_history_many.modify_device();
}

/* ----------------------------------------------------------------------
   the device exchange buffer has a fixed number of values per atom:
     ncontact, then for each of the tmax possible contacts
     the id of the wall and its history
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallGranRegionKokkos<DeviceType>::operator()(TagFixWallGranRegionPackExchange,
                                                     const int &mysend) const
{
  const int i = d_sendlist(mysend);
  int m = mysend * stride;
  d_buf(m++) = d_ubuf(d_ncontact(i)).d;
  for (int c = 0; c < tmax; c++) {
    d_buf(m++) = d_ubuf(d_walls(i,c)).d;
    for (int k = 0; k < size_history; k++)
      d_buf(m++) = static_cast<double>(d_history_many(i,c,k));
  }

  // move the last atom into the hole left by the sent one

  const int j = d_copylist(mysend);
  if (j > -1) {
    d_ncontact(i) = d_ncontact(j);
    for (int c = 0; c < tmax; c++) {
      d_walls(i,c) = d_walls(j,c);
      for (int k = 0; k < size_history; k++)
        d_history_many(i,c,k) = d_history_many(j,c,k);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixWallGranRegionKokkos<DeviceType>::pack_exchange_kokkos(
   const int &nsend, DAT::tdual_double_2d_lr &k_buf,
   DAT::tdual_int_1d k_sendlist, DAT::tdual_int_1d k_copylist,
   ExecutionSpace space)
{
  k_ncontact.template sync<DeviceType>();
  k_walls.template sync<DeviceType>();
  k_history_many.template sync<DeviceType>();
  d_ncontact = k_ncontact.template view<DeviceType>();
  d_walls = k_walls.template view<DeviceType>();
  d_history_many = k_history_many.template view<DeviceType>();

  k_buf.sync<DeviceType>();
  k_sendlist.sync<DeviceType>();
  k_copylist.sync<DeviceType>();

  d_sendlist = k_sendlist.view<DeviceType>();
  d_copylist = k_copylist.view<DeviceType>();

  d_buf = typename AT::t_double_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixWallGranRegionPackExchange>(0,nsend),
                       *this);

  copymode = 0;

  // the kernel also moves the history of the last atoms into the holes

  k_ncontact.template modify<DeviceType>();
  k_walls.template modify<DeviceType>();
  k_history_many.template modify<DeviceType>();

  // MPI sends the buffer from the exchange space, so make it current there

  k_buf.modify<DeviceType>();
  if (space == HostKK) k_buf.sync_host();
  else k_buf.sync_device();

  return nsend * stride;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallGranRegionKokkos<DeviceType>::operator()(TagFixWallGranRegionUnpackExchange,
                                                     const int &i) const
{
  const int index = d_indices(i);
  if (index > -1) {
    int m = i * stride;
    d_ncontact(index) = (int) d_ubuf(d_buf(m++)).i;
    for (int c = 0; c < tmax; c++) {
      d_walls(index,c) = (int) d_ubuf(d_buf(m++)).i;
      for (int k = 0; k < size_history; k++)
        d_history_many(index,c,k) = static_cast<KK_FLOAT>(d_buf(m++));
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranRegionKokkos<DeviceType>::unpack_exchange_kokkos(
  DAT::tdual_double_2d_lr &k_buf, DAT::tdual_int_1d &k_indices, int nrecv,
  int /*nrecv1*/, int /*nextrarecv1*/,
  ExecutionSpace /*space*/)
{
  k_buf.template sync<DeviceType>();
  k_indices.template sync<DeviceType>();

  d_buf = typename AT::t_double_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  d_indices = k_indices.view<DeviceType>();

  // the kernel below writes only the rows of the atoms that arrived, so the
  // rest have to be current on the device first

  k_ncontact.template sync<DeviceType>();
  k_walls.template sync<DeviceType>();
  k_history_many.template sync<DeviceType>();
  d_ncontact = k_ncontact.template view<DeviceType>();
  d_walls = k_walls.template view<DeviceType>();
  d_history_many = k_history_many.template view<DeviceType>();

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixWallGranRegionUnpackExchange>(0,nrecv),
                       *this);

  copymode = 0;

  k_ncontact.template modify<DeviceType>();
  k_walls.template modify<DeviceType>();
  k_history_many.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixWallGranRegionKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixWallGranRegionKokkos<LMPHostType>;
#endif
}
