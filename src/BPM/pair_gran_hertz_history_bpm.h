/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(gran/hertz/history/bpm,PairGranHertzHistoryBPM)

#else

#ifndef LMP_PAIR_GRAN_HERTZ_HISTORY_BPM_H
#define LMP_PAIR_GRAN_HERTZ_HISTORY_BPM_H

#include "pair_gran_hooke_history.h"

namespace LAMMPS_NS {

class PairGranHertzHistoryBPM : public PairGranHookeHistory {
 public:
  PairGranHertzHistoryBPM(class LAMMPS *);
  virtual void compute(int, int);
  void settings(int, char **);
  void init_style();
  double single(int, int, int, int, double, double, double, double &);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

  private:
  double C_exp;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
