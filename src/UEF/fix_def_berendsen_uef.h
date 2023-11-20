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
FixStyle(def/berendsen/uef,FixDefBerendsenUEF)
// clang-format on
#else

#ifndef LMP_FIX_DEF_BERENDSEN_UEF_H
#define LMP_FIX_DEF_BERENDSEN_UEF_H

#include "fix.h"

namespace LAMMPS_NS {
  // forward declaration
  namespace UEF_utils {
    class UEFBox;
  };

class FixDefBerendsenUEF : public Fix {
 public:
  int remapflag;                   // whether x,v are remapped across PBC
  int dimflag[6];                  // which dims are deformed
  int flip;

  FixDefBerendsenUEF(class LAMMPS *, int, char **);
  virtual ~FixDefBerendsenUEF();
  int setmask();
  void init();
  virtual void setup(int);
  virtual void pre_exchange();
  virtual void end_of_step();
  virtual void write_restart(FILE *);
  virtual void restart(char *buf);
  void get_rot(double[3][3]);
  void get_box(double[3][3]);
  double compute_scalar();
  double compute_vector(int);
  int modify_param(int, char **);

 protected:
  double rate[2],strain[3];
  double p_start, p_stop, p_period, p_target, p_current, dilation;

  char *id_temp, *id_press;
  class Compute *temperature, *pressure;
  int tflag, pflag;

  void rotate_x(double [3][3]);
  void inv_rotate_x(double[3][3]);
  void rotate_v(double[3][3]);
  void inv_rotate_v(double[3][3]);
  void rotate_f(double[3][3]);
  void inv_rotate_f(double[3][3]);

  class Irregular *irregular;
  UEF_utils::UEFBox *uefbox;
  double rot[3][3];

  bool nearly_equal(double,double,double);
};

}

#endif
#endif
