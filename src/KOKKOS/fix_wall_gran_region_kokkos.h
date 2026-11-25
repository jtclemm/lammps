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

#ifdef FIX_CLASS
// clang-format off
FixStyle(wall/gran/region/kk,FixWallGranRegionKokkos<LMPDeviceType>);
FixStyle(wall/gran/region/kk/device,FixWallGranRegionKokkos<LMPDeviceType>);
FixStyle(wall/gran/region/kk/host,FixWallGranRegionKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_WALL_GRAN_REGION_KOKKOS_H
#define LMP_FIX_WALL_GRAN_REGION_KOKKOS_H

#include "fix_wall_gran_old.h"

#include "fix_wall_gran_models_kokkos.h"
#include "kokkos_base.h"
#include "kokkos_type.h"
#include "region.h"

namespace LAMMPS_NS {

struct TagFixWallGranRegionPackExchange{};
struct TagFixWallGranRegionUnpackExchange{};

template<class DeviceType>
class FixWallGranRegionKokkos : public FixWallGranOld, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef int value_type;

  FixWallGranRegionKokkos(class LAMMPS *, int, char **);
  ~FixWallGranRegionKokkos() override;
  void init() override;
  void post_force(int) override;
  void write_restart(FILE *) override;
  void restart(char *) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;
  int pack_exchange_kokkos(const int &nsend, DAT::tdual_double_2d_lr &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space) override;
  void unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf,
                              DAT::tdual_int_1d &indices, int nrecv,
                              int nrecv1, int nrecv1extra,
                              ExecutionSpace space) override;

  // interaction of one particle with all its contacts with the region
  // RegionKK is the KOKKOS region class, Hertz = 1 for hertz/history,
  // error is set if a particle touches more walls than the history can hold

  template<class RegionKK, int Hertz>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void wall_particle(const RegionKK &, const int, int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallGranRegionPackExchange, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallGranRegionUnpackExchange, const int &) const;

 private:
  class Region *region;
  int tmax;                  // max # of region walls one particle can touch
  int motion_resetflag;      // set by restart() if the region motion must be reset
  int regiondynamic;         // 1 if the region moves or changes shape this step
  int stride;                // # of values per atom in the device exchange buffer

  // contact history for multiple contacts per particle
  // host side arrays alias the host side of the dual views

  int *ncontact;             // # of contacts per particle
  int **walls;               // which wall each contact is with
  double ***history_many;    // history per particle per contact

  DAT::tdual_int_1d k_ncontact;
  DAT::tdual_int_2d k_walls;
  DAT::ttransform_kkfloat_3d k_history_many;
  typename AT::t_int_1d d_ncontact;
  typename AT::t_int_2d d_walls;
  typename AT::t_kkfloat_3d d_history_many;

  WallGranKokkos::Coeffs coeffs;
  double dt_region;

  typename AT::t_kkfloat_1d_3_lr x;
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_kkfloat_1d_3 d_omega;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d mask;
  typename AT::t_kkfloat_1d rmass;
  typename AT::t_kkfloat_1d d_radius;

  typename AT::t_int_1d d_sendlist;
  typename AT::t_double_1d d_buf;
  typename AT::t_int_1d d_copylist;
  typename AT::t_int_1d d_indices;

  int region_supported() const;
  template<class RegionKK>
  int launch(RegionKK *, int);
  void sync_history_host();
  void modify_history_host();

  template<class RegionKK>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void velocity_contact(const RegionKK &, const double *, const Region::Contact &,
                        double *) const;
};

template <class DeviceType, class RegionKK, int Hertz>
struct FixWallGranRegionKokkosFunctor {
  typedef DeviceType device_type;
  typedef int value_type;
  FixWallGranRegionKokkos<DeviceType> c;

  // held by value so device kernels can use it

  RegionKK regionKK;

  FixWallGranRegionKokkosFunctor(FixWallGranRegionKokkos<DeviceType> *c_ptr,
                                 RegionKK *region_ptr) :
    c(*c_ptr), regionKK(*region_ptr) {}

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, int &error) const {
    c.template wall_particle<RegionKK, Hertz>(regionKK, i, error);
  }
};

}

#endif
#endif
