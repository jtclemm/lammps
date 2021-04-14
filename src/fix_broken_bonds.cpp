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

#include <string.h>
#include "fix_broken_bonds.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "group.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define DELTA 100

/* ---------------------------------------------------------------------- */

FixBrokenBonds::FixBrokenBonds(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nvalues(0),
  array(NULL), vector(NULL), pack_choice(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal fix broken bonds command");
  store_flag = 0;
  local_flag = 1;
  nvalues = narg - 4; 
  
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix broken bonds command");

  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  pack_choice = new FnPtrPack[nvalues];
  
  int i;
  for (int iarg = 4; iarg < narg; iarg++) {
    i = iarg -4;
    
    if (strcmp(arg[iarg],"id1") == 0) {
      pack_choice[i] = &FixBrokenBonds::pack_id1;
    } else if (strcmp(arg[iarg],"id2") == 0) {
      pack_choice[i] = &FixBrokenBonds::pack_id2;
           
    } else if (strcmp(arg[iarg],"time") == 0) {
      pack_choice[i] = &FixBrokenBonds::pack_time;      
      
    } else if (strcmp(arg[iarg],"x") == 0) {
      pack_choice[i] = &FixBrokenBonds::pack_x;    
    } else if (strcmp(arg[iarg],"y") == 0) {
      pack_choice[i] = &FixBrokenBonds::pack_y;    
    } else if (strcmp(arg[iarg],"z") == 0) {
      pack_choice[i] = &FixBrokenBonds::pack_z;    
      
    } else if (strcmp(arg[iarg],"xstore") == 0) {
      pack_choice[i] = &FixBrokenBonds::pack_xstore;    
      store_flag = 1;
    } else if (strcmp(arg[iarg],"ystore") == 0) {
      pack_choice[i] = &FixBrokenBonds::pack_ystore;    
      store_flag = 1;
    } else if (strcmp(arg[iarg],"zstore") == 0) {
      pack_choice[i] = &FixBrokenBonds::pack_zstore;      
      store_flag = 1;      
   
    } else error->all(FLERR, "Invalid keyword in fix broken bonds command");
  }

  nmax = 0;
  ncount = 0;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

FixBrokenBonds::~FixBrokenBonds()
{
  if (modify->nfix & store_flag == 1) {
    modify->delete_fix(id_fix);
    delete [] id_fix;
  }
    
  delete [] pack_choice;

  memory->destroy(vector);  
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

int FixBrokenBonds::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBrokenBonds::post_constructor()
{
  //If use stored x,y,z values, use store property (can transfer to ghost atoms) to store positions
  
  if(store_flag == 1){
      
    lx = domain->xprd;
    ly = domain->yprd;
    lz = domain->zprd;        

    int nn = strlen(id) + strlen("_FIX_PROP_ATOM") + 1;
    id_fix = new char[nn];
    strcpy(id_fix,id);
    strcat(id_fix,"_FIX_PROP_ATOM");
    
    int ifix = modify->find_fix(id_fix);
    if (ifix < 0) {
    
      int n_x = strlen(id) + 4;
      
      char * lab1 = new char[n_x];
      strcpy(lab1, "d_x");
      strcat(lab1, id);
      char * lab2 = new char[n_x];
      strcpy(lab2, "d_y");
      strcat(lab2, id);
      char * lab3 = new char[n_x];
      strcpy(lab3, "d_z");
      strcat(lab3, id);
        
      char **newarg = new char*[8];
      newarg[0] = id_fix;
      newarg[1] = group->names[igroup];
      newarg[2] = (char *) "property/atom"; 
      newarg[3] = (char *) lab1;         
      newarg[4] = (char *) lab2;         
      newarg[5] = (char *) lab3; 
      newarg[6] = (char *) "ghost";
      newarg[7] = (char *) "yes";

      modify->add_fix(8,newarg); 
      //Needs ghost atoms to calculate CoM

      int type_flag;
      int col_flag;

      strcpy(lab1, "x");
      strcat(lab1, id);
      strcpy(lab2, "y");
      strcat(lab2, id);
      strcpy(lab3, "z");
      strcat(lab3, id);

      index_x = atom->find_custom(lab1, type_flag, col_flag);
      index_y = atom->find_custom(lab2, type_flag, col_flag);
      index_z = atom->find_custom(lab3, type_flag, col_flag);
      delete [] newarg;    
      delete [] lab1;
      delete [] lab2;
      delete [] lab3;
    } 
      
    ifix = modify->find_fix(id_fix);
    if (ifix < 0) error->all(FLERR,"Could not find fix broken bonds fix ID");
    if (modify->fix[ifix]->restart_reset) {
        modify->fix[ifix]->restart_reset = 0;
    } else {

      double *xi = atom->dvector[index_x];
      double *yi = atom->dvector[index_y];
      double *zi = atom->dvector[index_z];
      
      double **xs = atom->x;
      int nlocal = atom->nlocal;
      
      for (int i = 0; i < nlocal; i++){
        xi[i] = xs[i][0];
        yi[i] = xs[i][1];
        zi[i] = xs[i][2];   
      }
    }    
  }  
}

/* ---------------------------------------------------------------------- */

void FixBrokenBonds::init()
{
  // Set size of array/vector
  ncount = 0;

  if (ncount >= nmax) {
      reallocate(ncount);}
  size_local_rows = ncount;
 
}

/* ---------------------------------------------------------------------- */

void FixBrokenBonds::add_bond(int i, int j)
{    
    if(ncount == nmax) reallocate(ncount);
    
    index_i = i;
    index_j = j;
    
    // fill vector or array with local values
    if (nvalues == 1) {
      (this->*pack_choice[0])(0);
    } else {
      for (int n = 0; n < nvalues; n++)
        (this->*pack_choice[n])(n); 
    }    
//printf("%d: Adding bond at ncount %d (%d-%d tags %d %d)\n", update->ntimestep, ncount, index_i, index_j, atom->tag[index_i], atom->tag[index_j]);
    ncount += 1;    
}

/* ---------------------------------------------------------------------- */

void FixBrokenBonds::post_force(int /*vflag*/) 
{   
  if(update->ntimestep % nevery == 0){
    size_local_rows = ncount;
    ncount = 0;  
  //  printf("%d: resetting ncount\n", update->ntimestep);
  }
}


/* ---------------------------------------------------------------------- */

void FixBrokenBonds::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax <= n) nmax += DELTA;
  
  if (nvalues == 1) {
    memory->grow(vector,nmax,"fix_broken_bonds:vector");
    vector_local = vector;
  } else {
    memory->grow(array,nmax,nvalues,"fix_broken_bonds:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double FixBrokenBonds::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  bytes += nmax*2 * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword fix broken/bond can output
   the atom property is packed into array or vector
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void FixBrokenBonds::pack_id1(int n) 
{
  int i;
  tagint *tag = atom->tag;
  
  if(nvalues == 1)
    vector[ncount] = tag[index_i];
  else
    array[ncount][n] = tag[index_i];
}

/* ---------------------------------------------------------------------- */

void FixBrokenBonds::pack_id2(int n)
{
  int i;
  tagint *tag = atom->tag;

  if(nvalues == 1)
    vector[ncount] = tag[index_j];
  else
    array[ncount][n] = tag[index_j];
}


/* ---------------------------------------------------------------------- */

void FixBrokenBonds::pack_time(int n)
{
  bigint time = update->ntimestep;

  if(nvalues == 1)
    vector[ncount] = time;
  else
    array[ncount][n] = time;
}

/* ---------------------------------------------------------------------- */

void FixBrokenBonds::pack_x(int n)
{
  double lx_new = domain->xprd;
  double **x = atom->x; 
  double xj = x[index_j][0];
    
  if(x[index_i][0] - xj > lx_new/2){
      xj  += -lx_new;
  } else if (x[index_i][0] - xj < -lx_new/2){
      xj  += lx_new;
  }    
  
  if(nvalues == 1)
    vector[ncount] = (x[index_i][0] + xj)/2;
  else
    array[ncount][n] = (x[index_i][0] + xj)/2;
}


/* ---------------------------------------------------------------------- */

void FixBrokenBonds::pack_y(int n)
{
  double lx_new = domain->yprd;
  double **x = atom->x; 
  double xj = x[index_j][1];
    
  if(x[index_i][1] - xj > lx_new/2){
      xj  += -lx_new;
  } else if (x[index_i][1] - xj < -lx_new/2){
      xj  += lx_new;
  }    
  
  if(nvalues == 1)
    vector[ncount] = (x[index_i][1] + xj)/2;
  else
    array[ncount][n] = (x[index_i][1] + xj)/2;
}


/* ---------------------------------------------------------------------- */

void FixBrokenBonds::pack_z(int n)
{
  double lx_new = domain->zprd;
  double **x = atom->x; 
  double xj = x[index_j][2];
    
  if(x[index_i][2] - xj > lx_new/2){
      xj  += -lx_new;
  } else if (x[index_i][2] - xj < -lx_new/2){
      xj  += lx_new;
  }    
  
  if(nvalues == 1)
    vector[ncount] = (x[index_i][2] + xj)/2;
  else
    array[ncount][n] = (x[index_i][2] + xj)/2;
}


/* ---------------------------------------------------------------------- */

void FixBrokenBonds::pack_xstore(int n)
{
  double *arrayx = atom->dvector[index_x];
  double xj = arrayx[index_j];
        
  if(arrayx[index_i] - xj > lx/2){
      xj  += -lx;
  } else if (arrayx[index_i] - xj < -lx/2){
      xj  += lx;
  }     
  
  if(nvalues == 1)
    vector[ncount] = (arrayx[index_i] + xj)/2;
  else
    array[ncount][n] = (arrayx[index_i] + xj)/2;
}


/* ---------------------------------------------------------------------- */

void FixBrokenBonds::pack_ystore(int n)
{
  double *arrayx = atom->dvector[index_y];
  double xj = arrayx[index_j];
        
  if(arrayx[index_i] - xj > ly/2){
      xj  += -ly;
  } else if (arrayx[index_i] - xj < -ly/2){
      xj  += ly;
  }     
  
  if(nvalues == 1)
    vector[ncount] = (arrayx[index_i] + xj)/2;
  else
    array[ncount][n] = (arrayx[index_i] + xj)/2;
}


/* ---------------------------------------------------------------------- */

void FixBrokenBonds::pack_zstore(int n)
{
  double *arrayx = atom->dvector[index_z];
  double xj = arrayx[index_j];
        
  if(arrayx[index_i] - xj > lz/2){
      xj  += -lz;
  } else if (arrayx[index_i] - xj < -lz/2){
      xj  += lz;
  }     
  
  if(nvalues == 1)
    vector[ncount] = (arrayx[index_i] + xj)/2;
  else
    array[ncount][n] = (arrayx[index_i] + xj)/2;
}

