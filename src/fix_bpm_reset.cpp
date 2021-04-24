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

#include "fix_bpm_reset.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "neighbor.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
//#include "bond_dem_simple.h"
#include "bond_dem_beam.h"
#include "bond_dem_generic.h"
#include "fix_bond_store.h"

using namespace LAMMPS_NS;
using namespace FixConst;
enum{SIMPLE,GENERIC,BEAM};

/* ---------------------------------------------------------------------- */

FixBPMReset::FixBPMReset(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),fix_bond_store(NULL)
{
  if (narg != 4) error->all(FLERR,"Illegal fix bpm/reset command");

  fix_bond_store = NULL;
  max_mol = -1;

  // required args
  zero_flag = utils::inumeric(FLERR,arg[3],false,lmp);
}

/* ---------------------------------------------------------------------- */

FixBPMReset::~FixBPMReset()
{

}

/* ---------------------------------------------------------------------- */

int FixBPMReset::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBPMReset::init()
{
  //class BondDEMSimple bondDEMSimple = (BondDEMSimple *) force->bond_match("dem/simple");
  class BondDEMBeam *bondDEMBeam = (BondDEMBeam *) force->bond_match("dem/beam");
  class BondDEMGeneric *bondDEMGeneric = (BondDEMGeneric *) force->bond_match("dem/generic");

  //if (bondDEMSimple != NULL){
  //  bond_model = SIMPLE;
  //  fix_bond_store = (FixBondStore *) bondDEMSimple->fix_bond_store;
  if (bondDEMBeam != NULL) {
    bond_model = BEAM;
    fix_bond_store = (FixBondStore *) bondDEMBeam->fix_bond_store;
  } else if (bondDEMGeneric != NULL) {
    bond_model = GENERIC;
    fix_bond_store = (FixBondStore *) bondDEMGeneric->fix_bond_store;
  } else {
    error->all(FLERR,"Need dem bond style to use fix bpm/reset");
  }
}

/* ----------------------------------------------------------------------
   perform particle insertion
------------------------------------------------------------------------- */

void FixBPMReset::setup_pre_force(int /*vflag*/)
{
  int i,j,m,type;
  double delx,dely,delz,r,rsq;

  int *mol = atom->molecule;
  double **x = atom->x;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  double rinv;

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  
  int new_max_mol = -1;
  double **bondstore = fix_bond_store->bondstore;

  for (i = 0; i < atom->nlocal; i++) {
    for (m = 0; m < atom->num_bond[i]; m++) {
      type = atom->bond_type[i][m];
      //Skip if bond was turned off
      if(type < 0)
          continue;    

      j = atom->map(atom->bond_atom[i][m]);          
          
      if(j == -1) error->all(FLERR, "Atom missing in DEM bond");
      
      if(mol[i] <= max_mol or mol[j] <= max_mol) continue;
      if(mol[i] > new_max_mol) new_max_mol = mol[i];
      if(mol[j] > new_max_mol) new_max_mol = mol[j];      
      
      // Save orientation as pointing towards small tag
      if(tag[i] < tag[j]){
        delx = x[i][0] - x[j][0]; 
        dely = x[i][1] - x[j][1]; 
        delz = x[i][2] - x[j][2]; 
      } else {
        delx = x[j][0] - x[i][0]; 
        dely = x[j][1] - x[i][1]; 
        delz = x[j][2] - x[i][2];
      }
      
      r = sqrt(delx*delx + dely*dely + delz*delz);
      rinv = 1.0/r;
      fix_bond_store->update_atom_value(i, m, 0, r);
      fix_bond_store->update_atom_value(i, m, 1, delx*rinv); 
      fix_bond_store->update_atom_value(i, m, 2, dely*rinv); 
      fix_bond_store->update_atom_value(i, m, 3, delz*rinv);     
    }
  }
        
  if(max_mol < new_max_mol)
    max_mol = new_max_mol;  
}

/* ----------------------------------------------------------------------
   perform particle insertion
------------------------------------------------------------------------- */

void FixBPMReset::pre_force(int /*vflag*/)
{
  int i1,i2,n,m;
  double delx,dely,delz,r,r0,rsq;

  int *mol = atom->molecule;
  double **x = atom->x;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  double rinv;

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  
  int new_max_mol = -1;
  double **bondstore = fix_bond_store->bondstore;
    
  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    r0 = bondstore[n][0];
    
    if (not mask[i1] & groupbit) continue;
    if (not mask[i2] & groupbit) continue;
    if (bondlist[n][2] <= 0) continue;    
    //if(zero_flag and r0 != 0) continue; - note this doesn't actually work (anymore?) since fix property atom does not zero arrays after atoms are created, can load old memory
    if(mol[i1] <= max_mol or mol[i2] <= max_mol) continue;
    if(mol[i1] > new_max_mol) new_max_mol = mol[i1];
    if(mol[i2] > new_max_mol) new_max_mol = mol[i2];
    
    if(tag[i1] < tag[i2]){
      delx = x[i1][0] - x[i2][0];
      dely = x[i1][1] - x[i2][1];
      delz = x[i1][2] - x[i2][2];
    } else {
      delx = x[i2][0] - x[i1][0];
      dely = x[i2][1] - x[i1][1];
      delz = x[i2][2] - x[i1][2];
    }
    rsq = delx*delx + dely*dely + delz*delz; 
    r = sqrt(rsq);
    rinv = 1.0/r;
                    
    bondstore[n][0] = r;
    if(bond_model != SIMPLE) {
      bondstore[n][1] = delx*rinv;
      bondstore[n][2] = dely*rinv;
      bondstore[n][3] = delz*rinv;
    }
    
    if (i1 < atom->nlocal){
      for (m = 0; m < atom->num_bond[i1]; m++){
        if (atom->bond_atom[i1][m] == atom->tag[i2]){
          fix_bond_store->update_atom_value(i1, m, 0, r);
          if(bond_model != SIMPLE) {
            fix_bond_store->update_atom_value(i1, m, 1, delx*rinv);
            fix_bond_store->update_atom_value(i1, m, 2, dely*rinv);
            fix_bond_store->update_atom_value(i1, m, 3, delz*rinv);
          }
        }
      }
    }
        
    if (i2 < atom->nlocal){
      for (m = 0; m < atom->num_bond[i2]; m++){
        if (atom->bond_atom[i2][m] == atom->tag[i1]){
          fix_bond_store->update_atom_value(i2, m, 0, r);
          if(bond_model != SIMPLE) {          
            fix_bond_store->update_atom_value(i2, m, 1, delx*rinv);
            fix_bond_store->update_atom_value(i2, m, 2, dely*rinv);
            fix_bond_store->update_atom_value(i2, m, 3, delz*rinv);            
          }
        }
      }
    }
        
  }
   

  if(max_mol < new_max_mol)
    max_mol = new_max_mol;  
}

