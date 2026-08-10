/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(granular/kk,PairGranularKokkos<LMPDeviceType>);
PairStyle(granular/kk/device,PairGranularKokkos<LMPDeviceType>);
PairStyle(granular/kk/host,PairGranularKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_GRANULAR_KOKKOS_H
#define LMP_PAIR_GRANULAR_KOKKOS_H

#include "kokkos_type.h"
#include "pair_granular.h"
#include "pair_kokkos.h"

namespace LAMMPS_NS {

template <class DeviceType> class FixNeighHistoryKokkos;
struct TagPairGranularKokkosCompute {};

template <class DeviceType>
class PairGranularKokkos : public PairGranular {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairGranularKokkos(class LAMMPS *);
  ~PairGranularKokkos() override;
  void compute(int, int) override;
  void init_style() override;
  double init_one(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;

  struct Model {
    int normal, damping, tangential, rolling, twisting, heat;
    int tangential_history, rolling_history, twisting_history;
    int limit_damping, synchronized_verlet;
    KK_FLOAT normal_k, damping_coeff;
    KK_FLOAT tangential_k, tangential_xt, tangential_mu;
    KK_FLOAT rolling_k, rolling_damp, rolling_mu;
    KK_FLOAT twisting_k, twisting_damp, twisting_mu;
    KK_FLOAT heat_coeff;
  };

 protected:
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_randomread v;
  typename AT::t_kkfloat_1d_3_randomread omega;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread radius;
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;
  typename AT::t_int_2d d_firsttouch;
  typename AT::t_kkfloat_2d d_firsthistory;

  Kokkos::DualView<Model **,Kokkos::LayoutRight,DeviceType> k_models;
  typename Kokkos::DualView<Model **,Kokkos::LayoutRight,DeviceType>::t_dev_const models;
  Kokkos::DualView<KK_FLOAT *,DeviceType> k_temperature, k_heatflow;
  typename AT::t_kkfloat_1d d_temperature, d_heatflow;

  FixNeighHistoryKokkos<DeviceType> *fix_historyKK;
  int nlocal, newton_pair, history_update;
  KK_FLOAT dt_kk;
  KK_FLOAT special_lj[4];

  const char *history_fix_style() const override;
  void build_models();
  void validate_models() const;

  struct Contact {
    KK_FLOAT dx[3], vi[3], vj[3], wi[3], wj[3];
    KK_FLOAT radi, radj, meff, ti, tj;
    KK_FLOAT *history;
  };
  struct Result {
    KK_FLOAT force[3], torquei[3], torquej[3], heat;
    KK_FLOAT fs[3], fr[3], twist;
  };

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static void evaluate(const Model &, Contact &, Result &, const KK_FLOAT, const bool);

 public:
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGranularKokkosCompute, const int) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static int sbmask(const int j) { return j >> SBBITS & 3; }

  friend void pair_virial_fdotr_compute<PairGranularKokkos>(PairGranularKokkos *);
};

}    // namespace LAMMPS_NS

#endif
#endif
