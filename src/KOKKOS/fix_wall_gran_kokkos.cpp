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

#include "fix_wall_gran_kokkos.h"
#include "atom_kokkos.h"
#include "error.h"
#include "memory_kokkos.h"
#include "atom_masks.h"
#include "update.h"

using namespace LAMMPS_NS;

enum{XPLANE=0,YPLANE=1,ZPLANE=2,ZCYLINDER,REGION};
enum{HOOKE,HOOKE_HISTORY,HERTZ_HISTORY,BONDED_HISTORY};
enum{NONE,CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixWallGranKokkos<DeviceType>::FixWallGranKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixWallGranOld(lmp, narg, arg)
{
  // the per-atom array of the "contacts" keyword is not maintained here

  if (peratom_flag)
    error->all(FLERR, "Fix wall/gran/kk does not yet support the contacts keyword");

  kokkosable = 1;
  exchange_comm_device = sort_device = 1;
  maxexchange = size_history;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | OMEGA_MASK | TORQUE_MASK | RADIUS_MASK | RMASS_MASK | MASK_MASK;
  datamask_modify = F_MASK | TORQUE_MASK;

  memory->destroy(history_one);
  history_one = nullptr;
  grow_arrays(atom->nmax);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixWallGranKokkos<DeviceType>::~FixWallGranKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_history_one, history_one);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranKokkos<DeviceType>::init()
{
  FixWallGranOld::init();

  if (fix_rigid)
    error->all(FLERR, "Fix wall/gran/kk not yet compatible with rigid bodies");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranKokkos<DeviceType>::post_force(int /*vflag*/)
{
  // do not update shear history during setup

  history_update = 1;
  if (update->setupflag) history_update = 0;

  // contact model coefficients for the kernel

  coeffs.kn = static_cast<KK_FLOAT>(kn);
  coeffs.kt = static_cast<KK_FLOAT>(kt);
  coeffs.gamman = static_cast<KK_FLOAT>(gamman);
  coeffs.gammat = static_cast<KK_FLOAT>(gammat);
  coeffs.xmu = static_cast<KK_FLOAT>(xmu);
  coeffs.dt = static_cast<KK_FLOAT>(dt);
  coeffs.limit_damping = limit_damping;

  // set position of wall to initial settings and velocity to 0.0
  // if wiggle or shear, set wall position and velocity accordingly

  wlo = static_cast<KK_FLOAT>(lo);
  whi = static_cast<KK_FLOAT>(hi);
  vwall[0] = vwall[1] = vwall[2] = 0.0;
  if (wiggle) {
    double arg = omega * (update->ntimestep - time_origin) * dt;
    if (wallstyle == axis) {
      wlo = static_cast<KK_FLOAT>(lo + amplitude - amplitude*cos(arg));
      whi = static_cast<KK_FLOAT>(hi + amplitude - amplitude*cos(arg));
    }
    vwall[axis] = static_cast<KK_FLOAT>(amplitude*omega*sin(arg));
  } else if (wshear) vwall[axis] = static_cast<KK_FLOAT>(vshear);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  d_omega = atomKK->k_omega.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  d_radius = atomKK->k_radius.view<DeviceType>();
  int nlocal = atom->nlocal;

  atomKK->sync(execution_space,datamask_read);

  if (use_history) {
    k_history_one.template sync<DeviceType>();
    d_history_one = k_history_one.template view<DeviceType>();
  }

  copymode = 1;

  if (pairstyle == HOOKE)
    error->all(FLERR, "Fix wall/gran/kk doesn't yet support hooke style");
  else if (pairstyle == HOOKE_HISTORY)
    launch<0>(nlocal);
  else if (pairstyle == HERTZ_HISTORY)
    launch<1>(nlocal);
  else if (pairstyle == GRANULAR)
    error->all(FLERR, "Fix wall/gran/kk doesn't yet support granular style");

  atomKK->modified(execution_space,datamask_modify);

  if (use_history) k_history_one.template modify<DeviceType>();

  copymode = 0;
}

/* ----------------------------------------------------------------------
   run the kernel for the wall style, Hertz = 1 for hertz/history
------------------------------------------------------------------------- */

template<class DeviceType>
template<int Hertz>
void FixWallGranKokkos<DeviceType>::launch(int nlocal)
{
  if (wallstyle == XPLANE)
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType, TagFixWallGranHistory<XPLANE, Hertz>>(0, nlocal), *this);
  else if (wallstyle == YPLANE)
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType, TagFixWallGranHistory<YPLANE, Hertz>>(0, nlocal), *this);
  else if (wallstyle == ZPLANE)
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType, TagFixWallGranHistory<ZPLANE, Hertz>>(0, nlocal), *this);
  else if (wallstyle == ZCYLINDER)
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType, TagFixWallGranHistory<ZCYLINDER, Hertz>>(0, nlocal), *this);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int WallStyle, int Hertz>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallGranKokkos<DeviceType>::operator()(TagFixWallGranHistory<WallStyle, Hertz>, const int &i) const
{
  KK_FLOAT vwall_[3];
  vwall_[0] = vwall[0];
  vwall_[1] = vwall[1];
  vwall_[2] = vwall[2];

  const KK_FLOAT cylradius_kk = static_cast<KK_FLOAT>(cylradius);
  const KK_FLOAT vshear_kk = static_cast<KK_FLOAT>(vshear);

  if (mask[i] & groupbit) {
    KK_FLOAT radius = d_radius(i);

    KK_FLOAT dx = 0.0;
    KK_FLOAT dy = 0.0;
    KK_FLOAT dz = 0.0;

    // rwall = radius of curvature of the wall, only used by hertz/history

    KK_FLOAT rwall = 0.0;

    if (WallStyle == XPLANE) {
      KK_FLOAT del1 = x(i,0) - wlo;
      KK_FLOAT del2 = whi - x(i,0);
      if (del1 < del2) dx = del1;
      else dx = -del2;
    } else if (WallStyle == YPLANE) {
      KK_FLOAT del1 = x(i,1) - wlo;
      KK_FLOAT del2 = whi - x(i,1);
      if (del1 < del2) dy = del1;
      else dy = -del2;
    } else if (WallStyle == ZPLANE) {
      KK_FLOAT del1 = x(i,2) - wlo;
      KK_FLOAT del2 = whi - x(i,2);
      if (del1 < del2) dz = del1;
      else dz = -del2;
    } else if (WallStyle == ZCYLINDER) {
      KK_FLOAT delxy = Kokkos::sqrt(x(i,0)*x(i,0) + x(i,1)*x(i,1));
      KK_FLOAT delr = cylradius_kk - delxy;
      if (delr > radius) {
        dz = cylradius_kk;
      } else {
        dx = -delr/delxy * x(i,0);
        dy = -delr/delxy * x(i,1);
        // rwall = -2r_c if inside cylinder, 2r_c outside
        rwall = (delxy < cylradius_kk) ? -2*cylradius_kk : 2*cylradius_kk;
       if (wshear && axis != 2) {
          vwall_[0] += vshear_kk * x(i,1)/delxy;
          vwall_[1] += -vshear_kk * x(i,0)/delxy;
          vwall_[2] = 0.0;
        }
      }
    }

    KK_FLOAT rsq = dx*dx + dy*dy + dz*dz;

    if (rsq > radius*radius) {
      if (use_history)
        for (int j = 0; j < 3; j++)
          d_history_one(i,j) = 0.0;
    } else {
      // meff = effective mass of sphere = rmass(i)

      const KK_FLOAT vi[3] = {v(i,0), v(i,1), v(i,2)};
      const KK_FLOAT wi[3] = {d_omega(i,0), d_omega(i,1), d_omega(i,2)};
      KK_FLOAT fi[3], ti[3];

      WallGranKokkos::classic_history<Hertz>(dx, dy, dz, rsq, vwall_, vi, wi, rwall, radius,
                                             rmass(i), d_history_one(i,0), d_history_one(i,1),
                                             d_history_one(i,2), history_update, coeffs, fi, ti);

      f(i,0) += static_cast<KK_ACC_FLOAT>(fi[0]);
      f(i,1) += static_cast<KK_ACC_FLOAT>(fi[1]);
      f(i,2) += static_cast<KK_ACC_FLOAT>(fi[2]);
      torque(i,0) += static_cast<KK_ACC_FLOAT>(ti[0]);
      torque(i,1) += static_cast<KK_ACC_FLOAT>(ti[1]);
      torque(i,2) += static_cast<KK_ACC_FLOAT>(ti[2]);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranKokkos<DeviceType>::grow_arrays(int nmax)
{
  if (use_history) {
    k_history_one.sync_host(); // force reallocation on host
    k_history_one.modify_host();
    memoryKK->grow_kokkos(k_history_one,history_one,nmax,size_history,"wall/gran/kk:history_one");
    d_history_one = k_history_one.template view<DeviceType>();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranKokkos<DeviceType>::copy_arrays(int i, int j, int delflag)
{
  if (use_history) {
    k_history_one.sync_host();
    FixWallGranOld::copy_arrays(i,j,delflag);
    k_history_one.modify_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranKokkos<DeviceType>::set_arrays(int i)
{
  if (use_history) {
    k_history_one.sync_host();
    FixWallGranOld::set_arrays(i);
    k_history_one.modify_host();
  }
}

/* ----------------------------------------------------------------------
   sort local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  // always sort on the device

  k_history_one.sync_device();

  Sorter.sort(LMPDeviceType(), k_history_one.view_device());

  k_history_one.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixWallGranKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_history_one.sync_host();

  return FixWallGranOld::pack_exchange(i,buf);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixWallGranKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  int n = FixWallGranOld::unpack_exchange(nlocal,buf);

  k_history_one.modify_host();

  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixWallGranKokkos<DeviceType>::pack_restart(int i, double *buf)
{
  if (!use_history) return 0;

  k_history_one.sync_host();

  return FixWallGranOld::pack_restart(i,buf);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranKokkos<DeviceType>::unpack_restart(int nlocal, int nth)
{
  if (!use_history) return;

  FixWallGranOld::unpack_restart(nlocal,nth);

  k_history_one.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallGranKokkos<DeviceType>::operator()(TagFixWallGranPackExchange, const int &mysend) const
{
  const int i = d_sendlist(mysend);
  int m = mysend*size_history;
  for (int v = 0; v < size_history; v++)
    d_buf(m++) = static_cast<double>(d_history_one(i,v));

  const int j = d_copylist(mysend);
  if (j > -1) {
    for (int v = 0; v < size_history; v++) {
      d_history_one(i,v) = d_history_one(j,v);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixWallGranKokkos<DeviceType>::pack_exchange_kokkos(
   const int &nsend, DAT::tdual_double_2d_lr &k_buf,
   DAT::tdual_int_1d k_sendlist, DAT::tdual_int_1d k_copylist,
   ExecutionSpace space)
{
  k_history_one.template sync<DeviceType>();

  k_buf.sync<DeviceType>();
  k_sendlist.sync<DeviceType>();
  k_copylist.sync<DeviceType>();

  d_sendlist = k_sendlist.view<DeviceType>();
  d_copylist = k_copylist.view<DeviceType>();

  d_buf = typename AT::t_double_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixWallGranPackExchange>(0,nsend),*this);

  copymode = 0;

  // the kernel also moves the history of the last atoms into the holes

  k_history_one.template modify<DeviceType>();

  // MPI sends the buffer from the exchange space, so make it current there

  k_buf.modify<DeviceType>();
  if (space == HostKK) k_buf.sync_host();
  else k_buf.sync_device();

  return nsend*size_history;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixWallGranKokkos<DeviceType>::operator()(TagFixWallGranUnpackExchange, const int &i) const
{
  int index = d_indices(i);
  if (index > -1) {
    int m = i*size_history;
    for (int v = 0; v < size_history; v++)
      d_history_one(index,v) = static_cast<KK_FLOAT>(d_buf(m++));
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallGranKokkos<DeviceType>::unpack_exchange_kokkos(
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
  // rest have to be current on the device first.  syncing here also retires
  // any outstanding host claim, which the modify<DeviceType>() at the end
  // would otherwise hit as a concurrent modification

  k_history_one.template sync<DeviceType>();

  d_history_one = k_history_one.template view<DeviceType>();

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixWallGranUnpackExchange>(0,nrecv),*this);

  copymode = 0;

  k_history_one.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixWallGranKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixWallGranKokkos<LMPHostType>;
#endif
}
