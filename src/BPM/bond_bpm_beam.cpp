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

/* ----------------------------------------------------------------------
   Contributing authors: Chris Lorenz and Mark Stevens (SNL)
------------------------------------------------------------------------- */

#include "bond_bpm_beam.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "memory.h"
#include "error.h"
#include "modify.h"
#include "domain.h"
#include "math_extra.h"
#include "math_const.h"
#include "fix_bond_store.h"
#include "fix_broken_bonds.h"
#include "update.h"
#include "citeme.h"

#define TOLERANCE 1e-15
#define EPSILON 1e-10

using namespace LAMMPS_NS;
using namespace MathExtra;
using namespace MathConst;

// todo write single
// add support for granular? Redo pair overlay...
// Move elongational force to fbend
// Damp tangental velocity?

static const char cite_bond_bpm_beam[] =
  "bond bpm/beam command:\n\n"
  "@Article{Carmona08,\n"
  " author =  {Carmona, H. A. and Wittel, F. K. and Kun, F. and Herrmann, H. J.},\n"
  " title =   {Fragmentation processes in impact of spheres},\n"
  " journal = {Physical Review E},\n"
  " year =    2008,\n"
  " number =  5,\n"
  " pages =   {051302}\n"
  "}\n"
  "@Article{Andre12,\n"
  " author =  {Andre, Damien and Iordanoff, Ivan and Charles, Jean Luc and Neauport, Jerome},\n"
  " title =   {Discrete element method to simulate continuous material by using the cohesive beam model},\n"
  " journal = {Computer Methods in Applied Mechanics and Engineering},\n"
  " year =    2012,\n"
  " volume =  {213--216},\n"
  " issue =   {1--3},\n"
  " pages =   {113--125}\n"  
  "}\n\n";

/* ---------------------------------------------------------------------- */

BondBPMBeam::BondBPMBeam(LAMMPS *lmp) : Bond(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_bond_bpm_beam);    
    
  partial_flag = 1;
  fix_broken_bonds = NULL;
  fix_bond_store = NULL;
  
  break_in_comp_flag = 0; 
  overlay_pair_flag = 0;
  r0_max_estimate = 0;
}

/* ---------------------------------------------------------------------- */

BondBPMBeam::~BondBPMBeam()
{
  if(fix_bond_store) modify->delete_fix("BOND_STORE_BPM_BEAM");    
    
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(E);
    memory->destroy(G);
    memory->destroy(R);
    memory->destroy(gamma);
    memory->destroy(gammaw);
    memory->destroy(eps_th);
    memory->destroy(theta_th);
    memory->destroy(C_exp);
  }
}

/* ---------------------------------------------------------------------- */
// abitrary perpendicular to v, ignore smallest component, don't get 2 zeros

void BondBPMBeam::get_perp_vec(double* v, double* perp)
{
  // Always use z axis in 2d
  if (domain->dimension == 2) {
    perp[0] = 0.0;    
    perp[1] = 0.0;
    perp[2] = 1.0;
    return;
  }
  double nx_abs = fabs(v[0]);
  double ny_abs = fabs(v[1]);
  double nz_abs = fabs(v[2]);
  if (nx_abs <= ny_abs && nx_abs <= nz_abs) {
    perp[0] = 0.0;
    perp[1] = -v[2];
    perp[2] = v[1];
  } else if (ny_abs <= nz_abs) {
    perp[0] = -v[2];
    perp[1] = 0.0;
    perp[2] = v[0];
  } else {
    perp[0] = -v[1];
    perp[1] = v[0];
    perp[2] = 0.0;
  }    
  
  MathExtra::norm3(perp);  
}

/* ---------------------------------------------------------------------- */
// Find quaternion to rotate u to v via the shortest distance
// https://stackoverflow.com/questions/1171849/finding-quaternion-representing-the-rotation-from-one-vector-to-another

void BondBPMBeam::get_rot_quat(double* u, double* v, double* rot_quat)
{
  double cross[3], dot;
  
  double dx_round = round_if_zero(u[0]-v[0]);
  double dy_round = round_if_zero(u[1]-v[1]);
  double dz_round = round_if_zero(u[2]-v[2]);
  
  //If already aligned, return identity
  if (dx_round == 0.0 && dy_round == 0.0 && dz_round == 0.0) {
    rot_quat[0] = 1.0;
    rot_quat[1] = 0.0;
    rot_quat[2] = 0.0;
    rot_quat[3] = 0.0;
    return;
  }
  
  //If it's anti parallel, return 180 rotation around some orthogonal vector
  if (abs(dx_round) > 0.99 && dy_round == 0.0 && dz_round == 0.0) {
    double perp[3];
    get_perp_vec(u, perp);
    
    rot_quat[0] = 0;
    rot_quat[1] = perp[0];
    rot_quat[2] = perp[1];
    rot_quat[3] = perp[2];
    return;
  }    
   
  // get cross product (u x v) and dot
  MathExtra::cross3(u, v, cross);
  // Replace cross product in 2d
  if (domain->dimension == 2) {
    cross[0] = 0.0;    
    cross[1] = 0.0;
    cross[2] = (u[0]*v[1]-u[1]*v[0]);
  } 
  dot = dot3(u, v);
  
  // Construct quaternion
  rot_quat[0] = sqrt(lensq3(u)*lensq3(v));
  rot_quat[0] += dot;
  rot_quat[1] = cross[0];
  rot_quat[2] = cross[1];
  rot_quat[3] = cross[2];
  
  // Normalize final quaternion
  MathExtra::qnormalize(rot_quat);
}

/* ---------------------------------------------------------------------- */

double BondBPMBeam::round_if_zero(double value)
{
  double abs_value = fabs(value);
  
  if (abs_value < TOLERANCE) return 0.0;
  return value;
}

/* ---------------------------------------------------------------------- */

double BondBPMBeam::bound_pm_one(double value){
    double result = value;
    if (value < -1) {
      if (value+1 < -TOLERANCE) error->one(FLERR,"Exceeded tolerance for -1");
      result = -1.0;
    } else if (value > 1) {
      if (value-1 > TOLERANCE) error->one(FLERR,"Exceeded tolerance for +1");
      result = 1.0;
    }
    
    return result;
}

/* ---------------------------------------------------------------------- */

void BondBPMBeam::calc_theta(double* nij, double* nb, double* q_nij_ex, double* q_ex_nij, double* q_atom, double* theta)
{
  double q_ij_b[4], q_temp[4], q_tot[4], q_tot_rot[4], nb_current[3], dot;
  
  //Find quaternion, rotate nij to nb, could do once...
  get_rot_quat(nij, nb, q_ij_b);
  
  // Find quaternion to rotate nij to current nb
  // nb_current = qatom*nb*qatom^-1
  // nb = qij2b nij * qijb^-1
  // => nb_current = qatom*qij2b* nij (qatom*qij2b)^-1
  MathExtra::quatquat(q_atom, q_ij_b, q_tot);
  MathExtra::qnormalize(q_tot);
  
  // Rotate q_tot into FoR where nij -> ex for simplicity
  // qtot' = q_nij_to_ex * q_tot * q_ex_to_nij^-1
  MathExtra::quatquat(q_nij_ex, q_tot, q_temp);
  MathExtra::quatquat(q_temp, q_ex_nij, q_tot_rot);
  MathExtra::qnormalize(q_tot_rot);
  
  // Calculate Euler angles
  theta[0] = asin(q_tot_rot[1])*2.0;
  theta[1] = asin(q_tot_rot[2])*2.0;
  theta[2] = asin(q_tot_rot[3])*2.0;
  
  //Find absolute angle between current nb and nij for breaking criteria
  MathExtra::quatrotvec(q_atom, nb, nb_current);
  MathExtra::norm3(nb_current);    
  
  dot = dot3(nb_current, nij);
  dot = bound_pm_one(dot);
  theta[3] = acos(dot);    
}



/* ---------------------------------------------------------------------- */

double BondBPMBeam::store_bond(int n,int i,int j)
{
  int m,k;
  double delx, dely, delz, r, rinv;
  double **x = atom->x; 
  tagint *tag = atom->tag;
  double **bondstore = fix_bond_store->bondstore;

  if (tag[i] < tag[j]) {
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
  
  bondstore[n][0] = r;
  bondstore[n][1] = delx*rinv;
  bondstore[n][2] = dely*rinv;
  bondstore[n][3] = delz*rinv;

  for (m = 0; m < atom->num_bond[i]; m ++) {
    if (atom->bond_atom[i][m] == tag[j]) {
      fix_bond_store->update_atom_value(i, m, 0, r);
      fix_bond_store->update_atom_value(i, m, 1, delx*rinv); 
      fix_bond_store->update_atom_value(i, m, 2, dely*rinv); 
      fix_bond_store->update_atom_value(i, m, 3, delz*rinv); 
    }
  }
  
  for (m = 0; m < atom->num_bond[j]; m ++) {
    if (atom->bond_atom[j][m] == tag[i]) {
      fix_bond_store->update_atom_value(j, m, 0, r);
      fix_bond_store->update_atom_value(j, m, 1, delx*rinv); 
      fix_bond_store->update_atom_value(j, m, 2, dely*rinv); 
      fix_bond_store->update_atom_value(j, m, 3, delz*rinv); 
    }
  }
  
  return r;
}

/* ---------------------------------------------------------------------- */

void BondBPMBeam::store_data()
{        
  int j, type;
  double delx, dely, delz, r, rinv;
  double **x = atom->x; 
  int **bond_type = atom->bond_type;    
  tagint *tag = atom->tag;
  
  for (int i = 0; i < atom->nlocal; i ++) {
    for (int m = 0; m < atom->num_bond[i]; m ++) {
      type = bond_type[i][m];
              
      //Skip if bond was turned off
      if(type < 0)
        continue;                
              
      // map to find index n for tag
      j = atom->map(atom->bond_atom[i][m]);          
      if (j == -1) error->one(FLERR, "Atom missing in BPM bond");
      
      // Save orientation as pointing towards small tag
      if (tag[i] < tag[j]) {
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

  fix_bond_store->post_neighbor();    
}

/* ---------------------------------------------------------------------- */

void BondBPMBeam::compute(int eflag, int vflag)
{
  if (! fix_bond_store->stored_flag) {
    fix_bond_store->stored_flag = true;
    store_data();   
  }      
    
  int i1,i2,m,n,type,itype,jtype;
  double delx,dely,delz,ebond,fbond,felong,evdwl,fpair;
  double forces[3], force_rot[3], tori[3], tori_rot[3], torj[3], torj_rot[3];
  double rsq,r,r2inv,r0,r0inv,r0sqinv, rinv,dr,drsq,e,esq, ethinv;
  double fdrag, dot, delvx, delvy, delvz, wd, wd2;
  double stretch,bend,twist,breaking;
  double thetai[4], thetaj[4], nij[3], nb[3], qi[4], qj[4];
  double q_rot[4], q_inv_rot[4];   
  double fs1, fs2, fs3, radi1, radi2, tor1, tor2, tor3;
  double wr[3], wrn[3], wrt[3], tdamp;
  double ex[3] = {0};
  double I, Itor;
  ex[0] = 1.0;  
  
  ebond = 0.0;
  ev_init(eflag,vflag);

  if (vflag_global == 2)
    force->pair->vflag_either = force->pair->vflag_global = 1;

  double **cutsq = force->pair->cutsq;
  double **x = atom->x;
  double **v = atom->v;                                                             
  double **omega = atom->omega; 
  double **f = atom->f;
  double **torque = atom->torque;
  double **quat = atom->quat;
  double *radius = atom->radius; 
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  
  double **bondstore = fix_bond_store->bondstore;

  for (n = 0; n < nbondlist; n++) {

    // skip bond if already broken

    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    r0 = bondstore[n][0];
    
    // If bond hasn't been set (always initialized to zero?)
    if (r0 < EPSILON || isnan(r0))
      r0 = store_bond(n,i1,i2);  
    
    nb[0] = bondstore[n][1];
    nb[1] = bondstore[n][2];
    nb[2] = bondstore[n][3];
    
    //Stored nb points small tag to large tag
    if(tag[i1] > tag[i2]){
      nb[0] = -nb[0];    
      nb[1] = -nb[1];
      nb[2] = -nb[2];      
    }    
    
    qj[0] = quat[i1][0]; // vector rij points twoards i1 (j) atom
    qj[1] = quat[i1][1]; // use notation from Herrmann
    qj[2] = quat[i1][2];
    qj[3] = quat[i1][3];
    
    qi[0] = quat[i2][0];
    qi[1] = quat[i2][1];
    qi[2] = quat[i2][2];
    qi[3] = quat[i2][3];
    
    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);

    r0inv = 1.0/r0;
    r0sqinv = r0inv*r0inv;
    rinv = 1.0/r;
    
    nij[0] = delx*rinv;
    nij[1] = dely*rinv;
    nij[2] = delz*rinv;    

    get_rot_quat(nij, ex, q_rot);
    MathExtra::qconjugate(q_rot, q_inv_rot);
    
    calc_theta(nij, nb, q_rot, q_inv_rot, qi, thetai);
    calc_theta(nij, nb, q_rot, q_inv_rot, qj, thetaj);    

    e = r*r0inv - 1.0;
    esq = e*e;
    
    // check break criteria
    // Can still break in compression if beam twists
    ethinv = 1.0/eps_th[type];
    stretch = esq*ethinv*ethinv;
    if (! break_in_comp_flag && e < 0.0) stretch = 0.0;
    bend = std::max(thetai[3], thetaj[3])/theta_th[type];
    twist = fabs(thetai[0] - thetaj[0])/theta_th[type];
    breaking = stretch+bend+twist;

    if (breaking >= 1.0) {
      bondlist[n][2] = 0;
      for (m = 0; m < atom->num_bond[i1]; m++)
        if (atom->bond_atom[i1][m] == atom->tag[i2])
          atom->bond_type[i1][m] = 0;
      if (i2 < atom->nlocal)
        for (m = 0; m < atom->num_bond[i2]; m++)
          if (atom->bond_atom[i2][m] == atom->tag[i1])
            atom->bond_type[i2][m] = 0;
        
      if (fix_broken_bonds != NULL) fix_broken_bonds->add_bond(i1, i2);          
        
      continue;
    }

    wd = 1.0 - breaking;     
    wd2 = wd*wd;

    // Elongational and bend forces
    felong = -E[type]*MY_PI*R[type]*R[type]*e;
    if (! break_in_comp_flag && e < 0.0) felong *= exp(C_exp[type]*(r0*rinv-1.0));
    
    I = MY_PI*0.25*pow(R[type], 4);
    force_rot[0] = felong;
    force_rot[1] = +3.0*E[type]*I*(thetai[2]+thetaj[2])*r0sqinv; //Q on j, "Qzb"
    force_rot[2] = -3.0*E[type]*I*(thetai[1]+thetaj[1])*r0sqinv; // -1 for RHR, "Qyb"
    //Rotate image st ey is in ez direction, then ez points down
    
    // Torques
    Itor = I*2;
    tori_rot[0] = -G[type]*Itor*(thetai[0] - thetaj[0])*r0inv;
    torj_rot[0] = -G[type]*Itor*(thetaj[0] - thetai[0])*r0inv;

    tori_rot[1] = E[type]*I*(thetaj[1] + 2*thetai[1])*r0inv;
    torj_rot[1] = E[type]*I*(thetai[1] + 2*thetaj[1])*r0inv;

    tori_rot[2] = E[type]*I*(thetaj[2] + 2*thetai[2])*r0inv;
    torj_rot[2] = E[type]*I*(thetai[2] + 2*thetaj[2])*r0inv;

    if (eflag) ebond = -wd; //Not really, but it's just to give a sense

    // apply force to each of 2 atoms
    MathExtra::quatrotvec(q_inv_rot, force_rot, forces);
    MathExtra::quatrotvec(q_inv_rot, torj_rot, torj);
    MathExtra::quatrotvec(q_inv_rot, tori_rot, tori);
    
    // Damping relative tangental rotation
    wr[0] = (omega[i1][0] + omega[i2][0]);
    wr[1] = (omega[i1][1] + omega[i2][1]);
    wr[2] = (omega[i1][2] + omega[i2][2]);    
    
    // Subtract out normal component
    dot = MathExtra::dot3(wr, nij);
    wrt[0] = wr[0] - dot*nij[0];
    wrt[1] = wr[1] - dot*nij[1];
    wrt[2] = wr[2] - dot*nij[2];
    
    tori[0] -= gammaw[type]*wrt[0];
    tori[1] -= gammaw[type]*wrt[1];
    tori[2] -= gammaw[type]*wrt[2];
    
    torj[0] -= gammaw[type]*wrt[0];
    torj[1] -= gammaw[type]*wrt[1];
    torj[2] -= gammaw[type]*wrt[2];
    
    // Damping relative normal component, twist
    wr[0] = (omega[i1][0] - omega[i2][0]);
    wr[1] = (omega[i1][1] - omega[i2][1]);
    wr[2] = (omega[i1][2] - omega[i2][2]);   
    
    dot = MathExtra::dot3(wr, nij);
    wrn[0] = dot*nij[0];
    wrn[1] = dot*nij[1];
    wrn[2] = dot*nij[2];    
    
    tori[0] += gammaw[type]*wrn[0];
    tori[1] += gammaw[type]*wrn[1];
    tori[2] += gammaw[type]*wrn[2];
            
    torj[0] -= gammaw[type]*wrn[0];
    torj[1] -= gammaw[type]*wrn[1];
    torj[2] -= gammaw[type]*wrn[2];   

    // DPD Damping force
    r2inv = 1.0/rsq;    
    delvx = v[i1][0] - v[i2][0];                                       
    delvy = v[i1][1] - v[i2][1];                                       
    delvz = v[i1][2] - v[i2][2];     
    dot = delx*delvx + dely*delvy + delz*delvz;
    //No factor of wd b/c added later
    //fdrag = -gamma[type]; 
    fdrag = -gamma[type]*dot*r2inv;

    forces[0] += delx*fdrag;
    forces[1] += dely*fdrag;
    forces[2] += delz*fdrag;

    if (newton_bond || i1 < nlocal) {

      f[i1][0] += forces[0]*wd2;
      f[i1][1] += forces[1]*wd2;
      f[i1][2] += forces[2]*wd2;   
      
      torque[i1][0] += torj[0]*wd2;
      torque[i1][1] += torj[1]*wd2; 
      torque[i1][2] += torj[2]*wd2;
    }

    if (newton_bond || i2 < nlocal) {
      
      f[i2][0] -= forces[0]*wd2;
      f[i2][1] -= forces[1]*wd2;
      f[i2][2] -= forces[2]*wd2; 

      torque[i2][0] += tori[0]*wd2;
      torque[i2][1] += tori[1]*wd2; 
      torque[i2][2] += tori[2]*wd2;      
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,felong,delx,dely,delz);

    // subtract out pairwise contribution from 2 atoms via pair->single()

    // Skip if overlaying
    if (overlay_pair_flag) continue;

    itype = atom->type[i1];
    jtype = atom->type[i2];
    
    if (rsq < cutsq[itype][jtype]) {
      evdwl = -force->pair->single(i1,i2,itype,jtype,rsq,1.0,1.0,fpair);
      fpair = -fpair;
      
      fs1 = -force->pair->svector[0];
      fs2 = -force->pair->svector[1];
      fs3 = -force->pair->svector[2];
      radi1 = radius[i1];
      radi2 = radius[i2];
      tor1 = rinv * (dely*fs3 - delz*fs2);
      tor2 = rinv * (delz*fs1 - delx*fs3);
      tor3 = rinv * (delx*fs2 - dely*fs1);
        
      //Note only works with my modified pair granular for now      
      if (newton_bond || i1 < nlocal) {
        f[i1][0] += delx*fpair + fs1;
        f[i1][1] += dely*fpair + fs2;
        f[i1][2] += delz*fpair + fs3;          
        torque[i1][0] -= radi1*tor1;
        torque[i1][1] -= radi1*tor2;
        torque[i1][2] -= radi1*tor3;
      }
      if (newton_bond || i2 < nlocal) {
        f[i2][0] -= delx*fpair + fs1;
        f[i2][1] -= dely*fpair + fs2;
        f[i2][2] -= delz*fpair + fs3;          
        torque[i2][0] -= radi2*tor1;
        torque[i2][1] -= radi2*tor2;
        torque[i2][2] -= radi2*tor3;
      }
    
      if (evflag) force->pair->ev_tally(i1,i2,nlocal,newton_bond,
                                        evdwl,0.0,fpair,delx,dely,delz);
    }
  }
}

/* ---------------------------------------------------------------------- */

void BondBPMBeam::settings(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "break/compress") == 0) {
      break_in_comp_flag = 1;
    } else if (strcmp(arg[iarg], "overlay/pair") == 0) {
      overlay_pair_flag = 1;
    } else error->all(FLERR,"Illegal pair_style command");
    iarg++;   
  }
}

/* ---------------------------------------------------------------------- */

void BondBPMBeam::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(E,n+1,"bond:E");
  memory->create(G,n+1,"bond:G");
  memory->create(R,n+1,"bond:R");
  memory->create(gamma,n+1,"bond:gamma");
  memory->create(gammaw,n+1,"bond:gammaw");
  memory->create(eps_th,n+1,"bond:eps_th");
  memory->create(theta_th,n+1,"bond:theta_th");
  memory->create(C_exp,n+1,"bond:C_exp");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;

}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondBPMBeam::coeff(int narg, char **arg)
{
  if (narg != 9) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nbondtypes,ilo,ihi,error);

  double E_one = utils::numeric(FLERR,arg[1],false,lmp);
  double G_one = utils::numeric(FLERR,arg[2],false,lmp);
  double R_one = utils::numeric(FLERR,arg[3],false,lmp);
  double eps_th_one = utils::numeric(FLERR,arg[4],false,lmp);
  double theta_th_one = utils::numeric(FLERR,arg[5],false,lmp);
  double gamma_one = utils::numeric(FLERR,arg[6],false,lmp);
  double gammaw_one = utils::numeric(FLERR,arg[7],false,lmp);
  double Cexp_one = utils::numeric(FLERR,arg[8],false,lmp);
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    E[i] = E_one;
    G[i] = G_one;
    R[i] = R_one;
    eps_th[i] = eps_th_one;
    theta_th[i] = theta_th_one;
    gamma[i] = gamma_one;
    gammaw[i] = gammaw_one;
    C_exp[i] = Cexp_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");

}

/* ----------------------------------------------------------------------
   check if pair defined and special_bond settings are valid
------------------------------------------------------------------------- */

void BondBPMBeam::init_style()
{
  if (!atom->quat_flag)
    error->all(FLERR,"Bond beam requires atom attributes quaternion");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Bond beam requires ghost atoms store velocity");

  if (force->pair == NULL || force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support bond_style bpm/beam");

  if (force->angle || force->dihedral || force->improper)
    error->all(FLERR,
               "Bond style beam cannot be used with 3,4-body interactions");
  if (atom->molecular == 2)
    error->all(FLERR,
               "Bond style beam cannot be used with atom style template");

  if (atom->special_flag)
    error->all(FLERR, "Special bonds must be turned off for bond style beam");

  if (domain->dimension == 2)
    error->warning(FLERR, "Bond style beam not intended for 2d use, may be inefficient");

  // Determine if correct pair style is used
  if (not overlay_pair_flag) {
    int correct_pair = 0;
    if (force->pair_match("gran/hooke",0)) correct_pair = 1;
    if (force->pair_match("gran/hooke/history",0)) correct_pair = 1;
    if (force->pair_match("gran/hertz",0)) correct_pair = 1;
    if (force->pair_match("gran/hertz/history",0)) correct_pair = 1;
    if (! correct_pair) 
      error->all(FLERR, "Bond style bpm/beam requires gran pairstyle without overlay");
  }    
    
  //Define bond store
  if (fix_bond_store == NULL) {
    char **fixarg = new char*[5];
    fixarg[0] = (char *) "BOND_STORE_BPM_BEAM";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "BOND_STORE";
    fixarg[3] = (char *) "0";
    fixarg[4] = (char *) "4";
    modify->add_fix(5,fixarg,1);
    delete [] fixarg;
    int ifix = modify->find_fix("BOND_STORE_BPM_BEAM");
    fix_bond_store = (FixBondStore *) modify->fix[ifix];
    //Note don't use most recent nfix b/c fix bond store creates a fix property atom    
  }
  
  int ifix = modify->find_fix_by_style("bonds/broken");
  if (ifix != -1) fix_broken_bonds = (FixBrokenBonds *) modify->fix[ifix];      
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length - not perfect, estimates based on local-local only
------------------------------------------------------------------------- */

double BondBPMBeam::equilibrium_distance(int i)
{
  // Ghost atoms not yet communicated, so some will be skipped
  if (r0_max_estimate == 0) {
    int type, j;
    double delx, dely, delz, r;
    double **x = atom->x;
    for (int i = 0; i < atom->nlocal; i ++) {
      for (int m = 0; m < atom->num_bond[i]; m ++) {
        type = atom->bond_type[i][m];
        if (type == 0)
            continue;                
                
        j = atom->map(atom->bond_atom[i][m]);        
        if(j == -1) continue;
        
        delx = x[i][0] - x[j][0]; 
        dely = x[i][1] - x[j][1]; 
        delz = x[i][2] - x[j][2]; 
        
        r = sqrt(delx*delx + dely*dely + delz*delz);
        if(r > r0_max_estimate) r0_max_estimate = r;
      }
    }
    
    double temp;
    MPI_Allreduce(&r0_max_estimate,&temp,1,MPI_DOUBLE,MPI_MAX,world);
    r0_max_estimate = temp;  
    
    if (comm->me == 0)
      utils::logmesg(lmp,fmt::format("Estimating longest bond = {}\n",r0_max_estimate));
  }
  
  double r_break = r0_max_estimate*(1+eps_th[i]);
  
  // Divide out heuristic prefactor added in comm class
  return r_break/1.5;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondBPMBeam::write_restart(FILE *fp)
{
  fwrite(&E[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&G[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&R[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&eps_th[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&theta_th[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&gamma[1],sizeof(double),atom->nbondtypes,fp);  
  fwrite(&gammaw[1],sizeof(double),atom->nbondtypes,fp);  
  fwrite(&C_exp[1],sizeof(double),atom->nbondtypes,fp);  
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondBPMBeam::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&E[1],sizeof(double),atom->nbondtypes,fp);
    fread(&G[1],sizeof(double),atom->nbondtypes,fp);
    fread(&R[1],sizeof(double),atom->nbondtypes,fp);
    fread(&eps_th[1],sizeof(double),atom->nbondtypes,fp);
    fread(&theta_th[1],sizeof(double),atom->nbondtypes,fp);
    fread(&gamma[1],sizeof(double),atom->nbondtypes,fp);    
    fread(&gammaw[1],sizeof(double),atom->nbondtypes,fp);    
    fread(&C_exp[1],sizeof(double),atom->nbondtypes,fp);    
  }
  MPI_Bcast(&E[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&G[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&R[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&eps_th[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&theta_th[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&gamma[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&gammaw[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&C_exp[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  
  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondBPMBeam::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
     fprintf(fp,"%d %g %g %g %g %g %g %g %g\n",i,E[i],G[i],R[i],eps_th[i],theta_th[i],gamma[i],gammaw[i], C_exp[i]);
}

/* ---------------------------------------------------------------------- */

double BondBPMBeam::single(int type, double rsq, int i, int j,
                           double &fforce)
{
  // Incomplete
  if (type <= 0) return 0.0;
  double r0; 
  
  //Access stored values - inefficient call rarely (I think only compute local)
  for(int n = 0; n < atom->num_bond[i]; n ++){
    if(atom->bond_atom[i][n] == atom->tag[j]){
      r0 = fix_bond_store->get_atom_value(i, n, 0);
    }
  }  

  double r = sqrt(rsq);
  double dr = r - r0;
  double drsq = dr*dr;

  fforce = 0;
  if (r > 0.0) fforce = -E[type]*MY_PI*R[type]*R[type] * dr/r;
  return E[type]*MY_PI*R[type]*R[type] * drsq;
}
