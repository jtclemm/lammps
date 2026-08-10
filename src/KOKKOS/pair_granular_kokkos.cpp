// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "pair_granular_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "fix_neigh_history_kokkos.h"
#include "force.h"
#include "gran_sub_mod.h"
#include "granular_model.h"
#include "kokkos.h"
#include "modify.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace Granular_NS;

namespace {
enum {NORMAL_HOOKE, NORMAL_HERTZ, NORMAL_HERTZ_MATERIAL};
enum {DAMP_NONE, DAMP_VELOCITY, DAMP_MASS_VELOCITY, DAMP_VISCOELASTIC, DAMP_TSUJI, DAMP_COR};
enum {TANG_NONE, TANG_LINEAR_NOHISTORY, TANG_LINEAR_HISTORY, TANG_MINDLIN, TANG_MINDLIN_FORCE};
enum {ROLL_NONE, ROLL_SDS};
enum {TWIST_NONE, TWIST_SDS, TWIST_MARSHALL};
enum {HEAT_NONE, HEAT_RADIUS, HEAT_AREA};
}

template<class DeviceType>
PairGranularKokkos<DeviceType>::PairGranularKokkos(LAMMPS *lmp) : PairGranular(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  no_virial_fdotr_compute = 0;
  datamask_read = X_MASK | V_MASK | OMEGA_MASK | F_MASK | TORQUE_MASK | TYPE_MASK | MASK_MASK |
      RMASS_MASK | RADIUS_MASK;
  datamask_modify = F_MASK | TORQUE_MASK;
  fix_historyKK = nullptr;
}

template<class DeviceType>
PairGranularKokkos<DeviceType>::~PairGranularKokkos()
{
  if (copymode) return;
}

template<class DeviceType>
const char *PairGranularKokkos<DeviceType>::history_fix_style() const
{
  return execution_space == Device ? "NEIGH_HISTORY/KK/DEVICE" : "NEIGH_HISTORY/KK/HOST";
}

template<class DeviceType>
void PairGranularKokkos<DeviceType>::validate_models() const
{
  const char *categories[] = {"normal", "damping", "tangential", "rolling", "twisting", "heat"};
  const char *supported[][6] = {
      {"hooke", "hertz", "hertz/material", nullptr, nullptr, nullptr},
      {"none", "velocity", "mass_velocity", "viscoelastic", "tsuji", "coeff_restitution"},
      {"none", "linear_nohistory", "linear_history", "mindlin", "mindlin/force", nullptr},
      {"none", "sds", nullptr, nullptr, nullptr, nullptr},
      {"none", "sds", "marshall", nullptr, nullptr, nullptr},
      {"none", "radius", "area", nullptr, nullptr, nullptr}};

  for (int n = 0; n < nmodels; ++n) {
    auto *model = models_list[n];
    for (int category = 0; category < NSUBMODELS; ++category) {
      const auto &name = model->sub_models[category]->name;
      bool found = false;
      for (int i = 0; supported[category][i]; ++i)
        if (name == supported[category][i]) found = true;
      if (!found) {
        if ((name == "mindlin_rescale") || (name == "mindlin_rescale/force"))
          error->all(FLERR, "Pair granular/kk does not support tangential model {} because it requires nondefault history transfer", name);
        error->all(FLERR, "Pair granular/kk does not support {} model {}", categories[category], name);
      }
    }
    if (model->nondefault_history_transfer)
      error->all(FLERR, "Pair granular/kk does not support a model requiring nondefault history transfer");
  }
}

template<class DeviceType>
void PairGranularKokkos<DeviceType>::init_style()
{
  validate_models();
  PairGranular::init_style();

  if (use_history) {
    fix_historyKK = dynamic_cast<FixNeighHistoryKokkos<DeviceType> *>(fix_history);
    if (!fix_historyKK) error->all(FLERR, "Pair granular/kk requires NEIGH_HISTORY/KK");
  }

  const int neighflag = lmp->kokkos->neighflag;
  if (neighflag == FULL)
    error->all(FLERR, "Pair granular/kk requires a half neighbor list");
  auto *request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
}

template<class DeviceType>
double PairGranularKokkos<DeviceType>::init_one(int i, int j)
{
  const double cutoff = PairGranular::init_one(i,j);
  build_models();
  return cutoff;
}

template<class DeviceType>
void PairGranularKokkos<DeviceType>::build_models()
{
  const int ntypes = atom->ntypes;
  k_models = Kokkos::DualView<Model **,Kokkos::LayoutRight,DeviceType>("granular:models",
                                                                         ntypes + 1, ntypes + 1);
  auto h_models = k_models.view_host();
  for (int i = 1; i <= ntypes; ++i) {
    for (int j = 1; j <= ntypes; ++j) {
      auto *gm = models_list[types_indices[i][j]];
      Model p{};
      auto *normal = gm->sub_models[NORMAL];
      auto *damping = gm->sub_models[DAMPING];
      auto *tangential = gm->sub_models[TANGENTIAL];
      auto *rolling = gm->sub_models[ROLLING];
      auto *twisting = gm->sub_models[TWISTING];
      auto *heat = gm->sub_models[HEAT];

      if (normal->name == "hooke") { p.normal = NORMAL_HOOKE; p.normal_k = normal->coeffs[0]; }
      else if (normal->name == "hertz") { p.normal = NORMAL_HERTZ; p.normal_k = normal->coeffs[0]; }
      else {
        p.normal = NORMAL_HERTZ_MATERIAL;
        const auto e = normal->coeffs[0], nu = normal->coeffs[2];
        p.normal_k = static_cast<KK_FLOAT>(2.0 * e / (3.0 * (1.0 - nu * nu)));
      }
      const double normal_damp = normal->coeffs[1];
      if (damping->name == "none") { p.damping = DAMP_NONE; p.damping_coeff = 0; }
      else if (damping->name == "velocity") { p.damping = DAMP_VELOCITY; p.damping_coeff = normal_damp; }
      else if (damping->name == "mass_velocity") { p.damping = DAMP_MASS_VELOCITY; p.damping_coeff = normal_damp; }
      else if (damping->name == "viscoelastic") { p.damping = DAMP_VISCOELASTIC; p.damping_coeff = normal_damp; }
      else {
        p.damping = damping->name == "tsuji" ? DAMP_TSUJI : DAMP_COR;
        if (p.damping == DAMP_TSUJI) {
          const double q = normal_damp;
          p.damping_coeff = 1.2728 - 4.2783*q + 11.087*q*q - 22.348*q*q*q +
              27.467*q*q*q*q - 18.022*q*q*q*q*q + 4.8218*q*q*q*q*q*q;
        } else {
          const double l = std::log(normal_damp);
          constexpr double pi = 3.14159265358979323846;
          p.damping_coeff = normal->name == "hooke" ?
              -2.0*l/std::sqrt(pi*pi+l*l) :
              -1.8257418583505538*l/std::sqrt(pi*pi+l*l);
        }
      }
      if (tangential->name == "none") p.tangential = TANG_NONE;
      else if (tangential->name == "linear_nohistory") p.tangential = TANG_LINEAR_NOHISTORY;
      else if (tangential->name == "linear_history") p.tangential = TANG_LINEAR_HISTORY;
      else if (tangential->name == "mindlin") p.tangential = TANG_MINDLIN;
      else p.tangential = TANG_MINDLIN_FORCE;
      p.tangential_k = tangential->num_coeffs ? tangential->coeffs[0] : 0;
      if ((p.tangential == TANG_MINDLIN) || (p.tangential == TANG_MINDLIN_FORCE)) {
        if (p.tangential_k < 0) {
          const double e = normal->coeffs[0], nu = normal->coeffs[2];
          p.tangential_k = 2.0 * e / ((2.0-nu)*(1.0+nu));
        }
      }
      p.tangential_xt = tangential->num_coeffs > 1 ? tangential->coeffs[1] : 0;
      p.tangential_mu = tangential->num_coeffs > 1 ?
          tangential->coeffs[tangential->num_coeffs - 1] : 0;
      p.tangential_history = tangential->history_index;

      p.rolling = rolling->name == "sds" ? ROLL_SDS : ROLL_NONE;
      p.rolling_k = rolling->num_coeffs ? rolling->coeffs[0] : 0;
      p.rolling_damp = rolling->num_coeffs ? rolling->coeffs[1] : 0;
      p.rolling_mu = rolling->num_coeffs ? rolling->coeffs[2] : 0;
      p.rolling_history = rolling->history_index;
      p.twisting = twisting->name == "sds" ? TWIST_SDS :
          (twisting->name == "marshall" ? TWIST_MARSHALL : TWIST_NONE);
      p.twisting_k = twisting->name == "marshall" ? p.tangential_k : (twisting->num_coeffs ? twisting->coeffs[0] : 0);
      p.twisting_damp = twisting->name == "marshall" ? p.tangential_xt : (twisting->num_coeffs ? twisting->coeffs[1] : 0);
      p.twisting_mu = twisting->name == "marshall" ? p.tangential_mu : (twisting->num_coeffs ? twisting->coeffs[2] : 0);
      p.twisting_history = twisting->history_index;
      p.heat = heat->name == "radius" ? HEAT_RADIUS : (heat->name == "area" ? HEAT_AREA : HEAT_NONE);
      p.heat_coeff = heat->num_coeffs ? heat->coeffs[0] : 0;
      p.limit_damping = gm->limit_damping;
      p.synchronized_verlet = gm->synchronized_verlet;
      h_models(i,j) = p;
    }
  }
  k_models.modify_host();
  models = k_models.template view<DeviceType>();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGranularKokkos<DeviceType>::evaluate(const Model &p, Contact &c, Result &out,
                                               const KK_FLOAT dt, const bool update)
{
  for (int k = 0; k < 3; ++k) out.force[k] = out.torquei[k] = out.torquej[k] =
      out.fs[k] = out.fr[k] = 0;
  out.heat = out.twist = 0;
  const KK_FLOAT rsq = c.dx[0]*c.dx[0] + c.dx[1]*c.dx[1] + c.dx[2]*c.dx[2];
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT rinv = 1.0/r;
  const KK_FLOAT delta = c.radi + c.radj - r;
  const KK_FLOAT reff = c.radi*c.radj/(c.radi+c.radj);
  const KK_FLOAT a = Kokkos::sqrt(delta*reff);
  KK_FLOAT n[3], vr[3], vt[3], vtr[3], wr[3];
  for (int k = 0; k < 3; ++k) { n[k] = c.dx[k]*rinv; vr[k] = c.vi[k]-c.vj[k]; }
  const KK_FLOAT vnnr = vr[0]*n[0]+vr[1]*n[1]+vr[2]*n[2];
  for (int k = 0; k < 3; ++k) vt[k] = vr[k]-vnnr*n[k];
  for (int k = 0; k < 3; ++k) wr[k] = c.radi*c.wi[k]+c.radj*c.wj[k];
  vtr[0] = vt[0] - (wr[1]*n[2]-wr[2]*n[1]);
  vtr[1] = vt[1] - (wr[2]*n[0]-wr[0]*n[2]);
  vtr[2] = vt[2] - (wr[0]*n[1]-wr[1]*n[0]);
  const KK_FLOAT vrel = Kokkos::sqrt(vtr[0]*vtr[0]+vtr[1]*vtr[1]+vtr[2]*vtr[2]);
  const KK_FLOAT fnormal = p.normal == NORMAL_HOOKE ? p.normal_k*delta : p.normal_k*a*delta;
  KK_FLOAT dpref = 0;
  if (p.damping == DAMP_VELOCITY) dpref = p.damping_coeff;
  else if (p.damping == DAMP_MASS_VELOCITY) dpref = p.damping_coeff*c.meff;
  else if (p.damping == DAMP_VISCOELASTIC) dpref = p.damping_coeff*c.meff*a;
  else if ((p.damping == DAMP_TSUJI) || (p.damping == DAMP_COR)) {
    const KK_FLOAT q = delta > 0 ? c.meff*fnormal/delta : 0;
    dpref = p.damping_coeff*Kokkos::sqrt(q > 0 ? q : 0);
  }
  KK_FLOAT fntot = fnormal-dpref*vnnr;
  if (p.limit_damping && fntot < 0) fntot = 0;
  const KK_FLOAT fncrit = Kokkos::fabs(fntot);

  if (p.tangential == TANG_LINEAR_NOHISTORY) {
    const KK_FLOAT mag = p.tangential_xt*dpref*vrel;
    const KK_FLOAT scale = vrel > 0 ? (mag < p.tangential_mu*fncrit ? mag : p.tangential_mu*fncrit)/vrel : 0;
    for (int k = 0; k < 3; ++k) out.fs[k] = -scale*vtr[k];
  } else if (p.tangential != TANG_NONE) {
    KK_FLOAT *h = c.history + p.tangential_history;
    KK_FLOAT shear[3] = {h[0],h[1],h[2]};
    if (update) {
      const KK_FLOAT dot = shear[0]*n[0]+shear[1]*n[1]+shear[2]*n[2];
      for (int k = 0; k < 3; ++k) shear[k] -= dot*n[k];
      const KK_FLOAT ks = p.tangential == TANG_LINEAR_HISTORY ? p.tangential_k : p.tangential_k*a;
      for (int k = 0; k < 3; ++k) shear[k] += (p.tangential == TANG_MINDLIN_FORCE ? -ks*dt : dt)*vtr[k];
    }
    const KK_FLOAT ks = p.tangential == TANG_LINEAR_HISTORY ? p.tangential_k : p.tangential_k*a;
    for (int k = 0; k < 3; ++k)
      out.fs[k] = p.tangential == TANG_MINDLIN_FORCE ? shear[k]-p.tangential_xt*dpref*vtr[k] :
          -ks*shear[k]-p.tangential_xt*dpref*vtr[k];
    const KK_FLOAT mag = Kokkos::sqrt(out.fs[0]*out.fs[0]+out.fs[1]*out.fs[1]+out.fs[2]*out.fs[2]);
    const KK_FLOAT lim = p.tangential_mu*fncrit;
    if (mag > lim && mag > 0) for (int k = 0; k < 3; ++k) out.fs[k] *= lim/mag;
    if (update) for (int k = 0; k < 3; ++k) h[k] = shear[k];
  }
  for (int k = 0; k < 3; ++k) out.force[k] = fntot*n[k]+out.fs[k];
  const KK_FLOAT di = c.radi-0.5*delta, dj = c.radj-0.5*delta;
  const KK_FLOAT cross[3] = {n[1]*out.fs[2]-n[2]*out.fs[1], n[2]*out.fs[0]-n[0]*out.fs[2],
                              n[0]*out.fs[1]-n[1]*out.fs[0]};
  for (int k = 0; k < 3; ++k) { out.torquei[k] = -di*cross[k]; out.torquej[k] = -dj*cross[k]; }

  if (p.rolling == ROLL_SDS) {
    KK_FLOAT *h = c.history+p.rolling_history, rel[3], vrl[3], fr[3];
    for (int k = 0; k < 3; ++k) rel[k] = c.wi[k]-c.wj[k];
    vrl[0] = reff*(rel[1]*n[2]-rel[2]*n[1]); vrl[1] = reff*(rel[2]*n[0]-rel[0]*n[2]);
    vrl[2] = reff*(rel[0]*n[1]-rel[1]*n[0]);
    if (update) for (int k = 0; k < 3; ++k) h[k] += dt*vrl[k];
    for (int k = 0; k < 3; ++k) fr[k] = -p.rolling_k*h[k]-p.rolling_damp*vrl[k];
    const KK_FLOAT mag = Kokkos::sqrt(fr[0]*fr[0]+fr[1]*fr[1]+fr[2]*fr[2]), lim = p.rolling_mu*fncrit;
    if (mag > lim && mag > 0) for (int k = 0; k < 3; ++k) fr[k] *= lim/mag;
    const KK_FLOAT tr[3] = {reff*(n[1]*fr[2]-n[2]*fr[1]), reff*(n[2]*fr[0]-n[0]*fr[2]),
                            reff*(n[0]*fr[1]-n[1]*fr[0])};
    for (int k = 0; k < 3; ++k) { out.fr[k] = fr[k]; out.torquei[k] += tr[k]; out.torquej[k] -= tr[k]; }
  }
  if (p.twisting != TWIST_NONE) {
    KK_FLOAT *h = c.history+p.twisting_history;
    const KK_FLOAT twist = (c.wi[0]-c.wj[0])*n[0]+(c.wi[1]-c.wj[1])*n[1]+(c.wi[2]-c.wj[2])*n[2];
    const KK_FLOAT k = p.twisting == TWIST_MARSHALL ? 0.5*p.twisting_k*a*a : p.twisting_k;
    const KK_FLOAT damp = p.twisting == TWIST_MARSHALL ? 0.5*p.twisting_damp*dpref*a*a : p.twisting_damp;
    const KK_FLOAT lim = (p.twisting == TWIST_MARSHALL ? 2.0*a/3.0*p.twisting_mu : p.twisting_mu)*fncrit;
    if (update) h[0] += twist*dt;
    KK_FLOAT mt = -k*h[0]-damp*twist;
    if (Kokkos::fabs(mt) > lim) mt = mt > 0 ? lim : -lim;
    out.twist = mt;
    for (int q = 0; q < 3; ++q) { out.torquei[q] += mt*n[q]; out.torquej[q] -= mt*n[q]; }
  }
  if (p.heat == HEAT_RADIUS) out.heat = 2*p.heat_coeff*a*(c.tj-c.ti);
  else if (p.heat == HEAT_AREA) out.heat = p.heat_coeff*3.141592653589793*a*a*(c.tj-c.ti);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGranularKokkos<DeviceType>::operator()(TagPairGranularKokkosCompute, const int ii) const
{
  const int i = d_ilist[ii];
  const int jnum = d_numneigh[i];
  for (int jj = 0; jj < jnum; ++jj) {
    int j = d_neighbors(i,jj);
    const KK_FLOAT factor = special_lj[sbmask(j)];
    j &= NEIGHMASK;
    if (factor == 0) continue;
    Contact c{};
    for (int k = 0; k < 3; ++k) { c.dx[k] = x(i,k)-x(j,k); c.vi[k] = v(i,k); c.vj[k] = v(j,k);
      c.wi[k] = omega(i,k); c.wj[k] = omega(j,k); }
    c.radi = radius[i]; c.radj = radius[j];
    const KK_FLOAT rsq = c.dx[0]*c.dx[0]+c.dx[1]*c.dx[1]+c.dx[2]*c.dx[2];
    if (rsq >= (c.radi+c.radj)*(c.radi+c.radj)) {
      if (use_history) { d_firsttouch(i,jj) = 0; for (int k = 0; k < size_history; ++k) d_firsthistory(i,size_history*jj+k) = 0; }
      continue;
    }
    const KK_FLOAT mi = rmass[i], mj = rmass[j];
    c.meff = mi*mj/(mi+mj);
    if (mask[i] & freeze_group_bit) c.meff = mj;
    if (mask[j] & freeze_group_bit) c.meff = mi;
    c.ti = heat_flag ? d_temperature[i] : 0; c.tj = heat_flag ? d_temperature[j] : 0;
    c.history = use_history ? &d_firsthistory(i,size_history*jj) : nullptr;
    if (use_history) d_firsttouch(i,jj) = 1;
    Result out{};
    evaluate(models(type[i],type[j]),c,out,dt_kk,history_update);
    for (int k = 0; k < 3; ++k) {
      Kokkos::atomic_add(&f(i,k),factor*out.force[k]); Kokkos::atomic_add(&torque(i,k),factor*out.torquei[k]);
      if (newton_pair || j < nlocal) { Kokkos::atomic_add(&f(j,k),-factor*out.force[k]);
        Kokkos::atomic_add(&torque(j,k),factor*out.torquej[k]); }
    }
    if (heat_flag) { Kokkos::atomic_add(&d_heatflow[i],out.heat);
      if (newton_pair || j < nlocal) Kokkos::atomic_add(&d_heatflow[j],-out.heat); }
  }
}

template<class DeviceType>
void PairGranularKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  ev_init(eflag_in,vflag_in,0);
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);
  k_models.template sync<DeviceType>();
  x = atomKK->k_x.view<DeviceType>(); v = atomKK->k_v.view<DeviceType>();
  omega = atomKK->k_omega.view<DeviceType>(); f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>(); type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>(); rmass = atomKK->k_rmass.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  nlocal = atom->nlocal; newton_pair = force->newton_pair; history_update = update->setupflag == 0; dt_kk = update->dt;
  for (int k = 0; k < 4; ++k) special_lj[k] = force->special_lj[k];
  auto *k_list = static_cast<NeighListKokkos<DeviceType> *>(list);
  d_neighbors = k_list->d_neighbors; d_ilist = k_list->d_ilist; d_numneigh = k_list->d_numneigh;
  if (use_history) {
    fix_historyKK->k_firstflag.template sync<DeviceType>(); fix_historyKK->k_firstvalue.template sync<DeviceType>();
    d_firsttouch = fix_historyKK->k_firstflag.template view<DeviceType>();
    d_firsthistory = fix_historyKK->k_firstvalue.template view<DeviceType>();
  }
  if (heat_flag) {
    const int nall = atom->nlocal + atom->nghost;
    k_temperature = Kokkos::DualView<KK_FLOAT *,DeviceType>("granular:temperature",nall);
    k_heatflow = Kokkos::DualView<KK_FLOAT *,DeviceType>("granular:heatflow",nall);
    auto ht = k_temperature.view_host(), hh = k_heatflow.view_host();
    for (int i = 0; i < nall; ++i) { ht(i) = atom->temperature[i]; hh(i) = 0; }
    k_temperature.modify_host(); k_heatflow.modify_host(); k_temperature.template sync<DeviceType>();
    k_heatflow.template sync<DeviceType>(); d_temperature = k_temperature.template view<DeviceType>(); d_heatflow = k_heatflow.template view<DeviceType>();
  }
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairGranularKokkosCompute>(0,list->inum),*this);
  copymode = 0;
  if (use_history) { fix_historyKK->k_firstflag.template modify<DeviceType>(); fix_historyKK->k_firstvalue.template modify<DeviceType>(); }
  if (heat_flag) {
    k_heatflow.template modify<DeviceType>(); k_heatflow.sync_host();
    auto hh = k_heatflow.view_host();
    for (int i = 0; i < nlocal; ++i) atom->heatflow[i] += hh(i);
  }
  if (vflag_fdotr) pair_virial_fdotr_compute(this);
}

template<class DeviceType>
double PairGranularKokkos<DeviceType>::single(int i, int j, int itype, int jtype,
                                              double, double, double, double &fforce)
{
  Model p = k_models.view_host()(itype,jtype);
  Contact c{};
  for (int k = 0; k < 3; ++k) { c.dx[k] = atom->x[i][k]-atom->x[j][k]; c.vi[k] = atom->v[i][k];
    c.vj[k] = atom->v[j][k]; c.wi[k] = atom->omega[i][k]; c.wj[k] = atom->omega[j][k]; }
  c.radi = atom->radius[i]; c.radj = atom->radius[j]; c.meff = atom->rmass[i]*atom->rmass[j]/(atom->rmass[i]+atom->rmass[j]);
  KK_FLOAT history[9] = {}; c.history = history;
  Result out{}; evaluate(p,c,out,update->dt,false);
  const double r = std::sqrt(c.dx[0]*c.dx[0]+c.dx[1]*c.dx[1]+c.dx[2]*c.dx[2]);
  fforce = (out.force[0]*c.dx[0]+out.force[1]*c.dx[1]+out.force[2]*c.dx[2])/(r*r);
  for (int k = 0; k < single_extra; ++k) svector[k] = 0;
  for (int k = 0; k < 3; ++k) { svector[k] = out.fs[k]; svector[4+k] = out.fr[k]; svector[9+k] = c.dx[k]; }
  svector[3] = std::sqrt(out.fs[0]*out.fs[0]+out.fs[1]*out.fs[1]+out.fs[2]*out.fs[2]);
  svector[7] = std::sqrt(out.fr[0]*out.fr[0]+out.fr[1]*out.fr[1]+out.fr[2]*out.fr[2]); svector[8] = out.twist;
  return 0.0;
}

template class PairGranularKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairGranularKokkos<LMPHostType>;
#endif
