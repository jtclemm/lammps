/* -*- c++ -*- ----------------------------------------------
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
PairStyle(granular/accelerated,PairGranularAccelerated);
// clang-format on
#else

#ifndef LMP_PAIR_GRANULAR_ACCELERATED_H
#define LMP_PAIR_GRANULAR_ACCELERATED_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGranularAccelerated : public Pair {
 public:
  PairGranularAccelerated(class LAMMPS *);
  ~PairGranularAccelerated() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void reset_dt() override;
  double single(int, int, int, int, double, double, double, double &) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double memory_usage() override;
  double atom2cut(int) override;
  double radii2cut(double, double) override;

 protected:
  int freeze_group_bit;
  int use_history;

  int neighprev;
  double *onerad_dynamic, *onerad_frozen;
  double *maxrad_dynamic, *maxrad_frozen;
  double **cut;

  class FixDummy *fix_dummy;
  class FixNeighHistory *fix_history;

  // storage of rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, null pointer if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  void allocate();
  void transfer_history(double *, double *, int, int) override;

 private:

  int heat_flag;

  // optional user-specified global cutoff, per-type user-specified cutoffs

  double **cutoff_type;
  double cutoff_global;

  void prep_model(int, int);

  // GranSubMod variables

  std::string *model_name;
  int *model_size_history;
  int *model_history_index;
  int *model_num_coeffs;
  double **model_transfer_history_factor;
  double *coeffs, ****model_coeffs;
  int max_num_coeffs;

  // GranularModel variables

  int contact_radius_flag, beyond_contact, limit_damping, history_update;
  int size_history, history_index, nondefault_history_transfer;
  double *history;
  double *transfer_history_factor;

  double Fnormal, forces[3], torquesi[3], torquesj[3], dq;

  double radi, radj, meff, dt, Ti, Tj, contact_radius;
  double Fntot, magtortwist;

  double *xi, *xj, *vi, *vj, *omegai, *omegaj;
  double fs[3], fr[3], ft[3];

  double dx[3], nx[3], r, rsq, rinv, Reff, radsum, delta, dR;
  double vr[3], vn[3], vnnr, vt[3], wr[3], vtr[3], vrl[3], relrot[3], vrel;
  double magtwist;

  int rolling_defined, twisting_defined, heat_defined; // Flag optional sub models

  // GranSubModNormal variables

  double k_norm, cohesion;
  double F_pulloff, Fne;
  double damp_norm;
  double Emix, Emod, poiss;
  double Fncrit;
  int material_properties, cohesive_flag;

  // GranSubModDamping variables

  double damp_prefactor, damp_damp;

  // GranSubModTangential variables

  double k_tang, damp_tang, mu_tang;
  int mindlin_rescale, mindlin_force;
  double xt;

  // GranSubModRolling variables

  double k_roll, mu_roll, gamma;

  // GranSubModTwisting variables

  double k_twist, mu_twist, damp_twist;

  // GranSubModHeat variables

  double conductivity, heat_transfer_coeff;

  // GranSubMod calculations

  double calculate_forces_normal();
  double calculate_forces_damping();
  void calculate_forces_tangential();
  void calculate_forces_rolling();
  void calculate_forces_twisting();
  double calculate_heat();

  void coeffs_to_local_normal();
  void coeffs_to_local_damping();
  void coeffs_to_local_tangential();
  void coeffs_to_local_rolling();
  void coeffs_to_local_twisting();
  void coeffs_to_local_heat();

  void init_normal();
  void init_damping();
  void init_tangential();
  void init_rolling();
  void init_twisting();
  void init_heat();

  bool touch_normal();
  double pulloff_distance_normal(double, double);
  double calculate_contact_radius_normal();
  void set_fncrit_normal();

  double mix_stiffnessE(double, double, double, double);
  double mix_stiffnessG(double, double, double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif
