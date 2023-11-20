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
PairStyle(lj/cut/coul/cut/dpd/tstat,PairLJCutCoulCutDPDTstat);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_COUL_CUT_DPD_TSTAT_H
#define LMP_PAIR_LJ_CUT_COUL_CUT_DPD_TSTAT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutCoulCutDPDTstat : public Pair {
 public:
  PairLJCutCoulCutDPDTstat(class LAMMPS *);
  virtual ~PairLJCutCoulCutDPDTstat();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  void born_matrix(int, int, int, int, double, double, double, double &, double &) override;
  void *extract(const char *, int &);

 protected:
  double temperature;
  int seed;
  double **cut_lj, **cut_ljsq;
  double **cut_coul, **cut_coulsq;
  double **cut_dpd, **cut_dpdsq, **cut_dpdinv;
  double **epsilon, **sigma, **sigma2, **gamma;
  double **lj1, **lj2, **lj3, **lj4, **offset;

  double t_start, t_stop;
  class RanMars *random;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
