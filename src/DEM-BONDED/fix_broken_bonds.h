/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(bonds/broken,FixBrokenBonds)

#else

#ifndef LMP_FIX_BROKEN_BONDS_H
#define LMP_FIX_BROKEN_BONDS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBrokenBonds : public Fix {
 public:
  FixBrokenBonds(class LAMMPS *, int, char **);
  ~FixBrokenBonds();
  int setmask();  
  void post_constructor();
  void init();
  void post_force(int);
  double memory_usage();
  void add_bond(int, int);

 private:
  int nvalues, nevery;
  int nmax;
  int store_flag;
  int index_i, index_j;
  
  double *vector;
  double **array;
  
  double lx;
  double ly;
  double lz;    

  int ncount;

  void reallocate(int);

  typedef void (FixBrokenBonds::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions
  
  void pack_id1(int);
  void pack_id2(int);
  
  void pack_time(int);  
  
  void pack_x(int);  
  void pack_y(int);  
  void pack_z(int);  
  
  void pack_xstore(int);  
  void pack_ystore(int);  
  void pack_zstore(int);   
  
   
  char *id_fix;
  int index_x, index_y, index_z;  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix property/local cannot use these inputs together

Only inputs that generate the same number of datums can be used
togther.  E.g. bond and angle quantities cannot be mixed.

E: Invalid keyword in compute property/local command

Self-explanatory.

E: Fix property/local does not (yet) work with atom_style template

Self-explanatory.

E: Fix property/local for property that isn't allocated

Self-explanatory.

E: No pair style is defined for compute property/local

Self-explanatory.

E: Pair style does not support compute property/local

The pair style does not have a single() function, so it can
not be invoked by fix bond/swap.

*/
