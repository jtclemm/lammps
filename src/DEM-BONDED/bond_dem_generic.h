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

BondStyle(dem/generic,BondDEMGeneric)

#else

#ifndef LMP_BOND_DEM_GENERIC_H
#define LMP_BOND_DEM_GENERIC_H

#include "bond.h"

namespace LAMMPS_NS {

class BondDEMGeneric : public Bond {
 public:
  BondDEMGeneric(class LAMMPS *);
  virtual ~BondDEMGeneric();
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
  double *Kr, *Ks, *Kt, *Kb, *gamma, *gammaw, *Fcr, *Fcs, *Gct, *Gcb, *C_exp;
  double max_r0;
int printflag;
  int overlay_pair_flag; // 1 => dont' subtract pair forces
  void calc_forces(int, double, double, double*, double*, double*,double*, double*, double*, double*, double &, double &, double &, double &);
  double acos_limit(double);

  class FixBrokenBonds *fix_broken_bonds;
  void allocate();
  void store_data();  
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

E: Bond style beam requires special_bonds = 1,1,1

This is a restriction of the current bond beam implementation.

*/
