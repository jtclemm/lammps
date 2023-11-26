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
FixStyle(def/uef/kk,FixDefUEFKokkos);
FixStyle(def/uef/kk/device,FixDefUEFKokkos);
FixStyle(def/uef/kk/host,FixDefUEFKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_DEF_UEF_KOKKOS_H
#define LMP_FIX_DEF_UEF_KOKKOS_H

#include "fix_def_uef.h"

namespace LAMMPS_NS {

class FixDefUEFKokkos : public FixDefUEF {
 public:
  FixDefUEFKokkos(class LAMMPS *, int, char **);

  void pre_exchange() override;
  void end_of_step() override;

  KOKKOS_INLINE_FUNCTION
  void rotate_x_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void inv_rotate_x_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void rotate_v_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void inv_rotate_v_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void rotate_f_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void inv_rotate_f_item(int) const;

 private:
  class DomainKokkos *domainKK;
  double *boxlo;

};

}

#endif
#endif

