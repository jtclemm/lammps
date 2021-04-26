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

#ifdef FIX_CLASS

FixStyle(BOND_STORE,FixBondStore)

#else

#ifndef LMP_FIX_BOND_STORE_H
#define LMP_FIX_BOND_STORE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondStore : public Fix {
 public:
 
  FixBondStore(class LAMMPS *, int, char **);
  ~FixBondStore();
  int setmask();
  void post_constructor();
  void init();
  void setup_post_neighbor();
  void post_neighbor();
  void pre_exchange();
  double memory_usage();
  void write_restart(FILE *fp);
  void restart(char *buf);
  void set_arrays(int);
  
  void update_atom_value(int, int, int, double);
  double get_atom_value(int, int, int);
  double **bondstore;   
  int stored_flag;

 protected:
 
  void allocate();
 
  int update_flag;
  int nbond, maxbond, ndata;
  int index;
  char *new_fix_id;
  char *array_id;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

E: Index exceeded in fix bond store

Bond requested non-existant data

E: Cannot store bond variables without any bonds

Atoms must have a nonzero number of bonds to store data



*/

