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
PairStyle(lj/cut/coul/long/dpd/tstat,PairLJCutCoulLongDPDTstat);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_COUL_LONG_DPD_TSTAT_H
#define LMP_PAIR_LJ_CUT_COUL_LONG_DPD_TSTAT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutCoulLongDPDTstat : public Pair {

 public:
  PairLJCutCoulLongDPDTstat(class LAMMPS *);
  ~PairLJCutCoulLongDPDTstat() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;

  void compute_inner() override;
  void compute_middle() override;
  void compute_outer(int, int) override;
  void *extract(const char *, int &) override;

 protected:
  double temperature;
  int seed;
  double **cut_lj, **cut_ljsq;
  double cut_coul, cut_coulsq;
  double **cut_dpd, **cut_dpdsq, **cut_dpdinv;
  double **epsilon, **sigma, **sigma2, **gamma;
  double **lj1, **lj2, **lj3, **lj4, **offset;
  double *cut_respa;
  double qdist;    // TIP4P distance from O site to negative charge
  double g_ewald;

  double t_start, t_stop;
  class RanMars *random;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
