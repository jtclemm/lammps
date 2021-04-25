/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS

BondStyle(bpm/beam,BondBPMBeam)

#else

#ifndef LMP_BOND_BPM_BEAM_H
#define LMP_BOND_BPM_BEAM_H

#include "bond.h"

namespace LAMMPS_NS {

class BondBPMBeam : public Bond {
 public:
  BondBPMBeam(class LAMMPS *);
  virtual ~BondBPMBeam();
  virtual void compute(int, int);
  virtual void settings(int, char **);  
  void coeff(int, char **);
  void init_style();
  double equilibrium_distance(int);  
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, double, int, int, double &);
 
  class FixBondStore *fix_bond_store;
 
 protected:
  double *E,*G, *R,*gamma,*gammaw,*eps_th,*theta_th,*C_exp;
  double r0_max_estimate;

  void get_perp_vec(double*, double*);
  void get_rot_quat( double*, double*, double*);
  void calc_theta( double*,  double*,  double*, double*, double*, double*);
  double bound_pm_one(double);
  double round_if_zero(double);

  int break_in_comp_flag; // 1 => can break in compression
  int overlay_pair_flag; // 1 => dont' subtract pair forces

  class FixBrokenBonds *fix_broken_bonds;
  void allocate();
  void store_data();  
  double store_bond(int, int, int);  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style does not support bond_style beam

The pair style does not have a single() function, so it can
not be invoked by bond_style beam.

E: Bond style beam cannot be used with 3,4-body interactions

No angle, dihedral, or improper styles can be defined when using
bond style beam.

E: Bond style beam cannot be used with atom style template

This bond style can change the bond topology which is not
allowed with this atom style.

E: Bond style bpm/beam requires gran pairstyle without overlay

Only the gran pairstyles are supportted unless pair forces are overlaid

E: Special bonds must be turned off for bond style beam

Special bond are not used by BPM bonds

W: Bond style bpm/beam not intended for 2d use, may be inefficient

This bond style will perform a lot of unnecessary calculations in 2d

*/
