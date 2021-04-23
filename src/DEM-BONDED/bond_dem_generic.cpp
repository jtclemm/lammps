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

#include "bond_dem_generic.h"
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

#define TOLERANCE 1e-15

using namespace LAMMPS_NS;
using namespace MathExtra;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

BondDEMGeneric::BondDEMGeneric(LAMMPS *lmp) : Bond(lmp)
{
  partial_flag = 1;
  fix_broken_bonds = NULL;
  fix_bond_store = NULL;
  
  overlay_pair_flag = 0;
  max_r0 = 1.0; // Placeholder
}

/* ---------------------------------------------------------------------- */

BondDEMGeneric::~BondDEMGeneric()
{
  if(fix_bond_store) modify->delete_fix("BOND_STORE_DEM_GENERIC");    
    
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(Kr);
    memory->destroy(Ks);
    memory->destroy(Kt);
    memory->destroy(Kb);
    memory->destroy(Fcr);
    memory->destroy(Fcs);
    memory->destroy(Gct);
    memory->destroy(Gcb);
    memory->destroy(gamma);
    memory->destroy(gammaw);
    memory->destroy(C_exp);
  }
}

/* ---------------------------------------------------------------------- */
// Calculate acos after bounding quantity in range [-1, 1]

double BondDEMGeneric::acos_limit(double c)
{
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;
  return acos(c);    
}

/* ---------------------------------------------------------------------- */

void BondDEMGeneric::store_data()
{        
  int j, type;
  double delx, dely, delz, r, rinv;
  double **x = atom->x; 
  int **bond_type = atom->bond_type;    
  tagint *tag = atom->tag;
  
  max_r0 = 0.0;
  for(int i = 0; i < atom->nlocal; i ++){
    for(int m = 0; m < atom->num_bond[i]; m ++){
      type = bond_type[i][m];
              
      //Skip if bond was turned off
      if(type < 0)
          continue;                
              
      // map to find index n for tag
      j = atom->map(atom->bond_atom[i][m]);          
      if(j == -1) error->all(FLERR, "Atom missing in DEM bond");
      
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
      domain->minimum_image(delx,dely,delz);        
      r = sqrt(delx*delx + dely*dely + delz*delz);
      rinv = 1.0/r;
      fix_bond_store->update_atom_value(i, m, 0, r);
      fix_bond_store->update_atom_value(i, m, 1, delx); 
      fix_bond_store->update_atom_value(i, m, 2, dely); 
      fix_bond_store->update_atom_value(i, m, 3, delz); 
      if(r > max_r0) max_r0 = r;
    }
  }
  
  double temp;
  MPI_Allreduce(&max_r0,&temp,1,MPI_DOUBLE,MPI_MAX,world);
  max_r0 = temp;
  fix_bond_store->post_neighbor();    
}

/* ---------------------------------------------------------------------- */

void BondDEMGeneric::calc_forces(int type, double r_mag, double r0_mag, double *q1, double *q2, double *r0, double *r, double *force1on2, double *torque1on2, double *torque2on1, double &Fs_mag, double &Fr_mag, double &Tb_mag, double &Tt_mag)  
{
  double q2inv[4], rb[3], rb_x_r0[3], s[3], t[3], Fs[3], q21[4], qp21[4], Tbp[3], Ttp[3];
  double Tsp[3], Fsp[3], m[4], minv[4], Ttmp[3], Ftmp[3], Tt[3], Tb[3], Ts[3], F_rot[3], T_rot[3], qtmp[4];
  double r0_dot_rb, gamma, c, psi, theta, sin_phi, cos_phi, temp, mag_in_plane, mag_out_plane;
  double r_mag_inv = 1.0/r_mag;
  
  // Calculate normal forces, rb = bond vector in particle 1's frame
  MathExtra::qconjugate(q2, q2inv);
  MathExtra::quatrotvec(q2inv, r, rb);
  
  Fr_mag = Kr[type]*(r_mag - r0_mag);
  
  if (r_mag < r0_mag) 
    Fr_mag *= exp(C_exp[type]*(r0_mag*r_mag_inv-1.0));
  
  F_rot[0] = Fr_mag*rb[0]*r_mag_inv;
  F_rot[1] = Fr_mag*rb[1]*r_mag_inv;
  F_rot[2] = Fr_mag*rb[2]*r_mag_inv; 
    
  // Calculate forces due to tangential displacements (no rotation)
  r0_dot_rb = dot3(r0, rb); 
  c = r0_dot_rb*r_mag_inv/r0_mag;
  gamma = acos_limit(c);
  
  MathExtra::cross3(rb, r0, rb_x_r0);
  MathExtra::cross3(rb, rb_x_r0, s);
  MathExtra::norm3(s);
  
  Fs[0] = Ks[type]*r_mag*gamma*s[0];
  Fs[1] = Ks[type]*r_mag*gamma*s[1];
  Fs[2] = Ks[type]*r_mag*gamma*s[2];
    
  // Calculate torque due to tangential displacements
  MathExtra::cross3(r0, rb, t);
  MathExtra::norm3(t);
  
  Ts[0] = 0.5*r_mag*Ks[type]*r_mag*gamma*t[0];
  Ts[1] = 0.5*r_mag*Ks[type]*r_mag*gamma*t[1];
  Ts[2] = 0.5*r_mag*Ks[type]*r_mag*gamma*t[2];
  
  // Relative rotation force/torque
  // Use representation of X'Y'Z' rotations from Wang, Mora 2009
  
  temp = r_mag + rb[2];
  if (temp < 0.0) temp = 0.0;
  m[0] = sqrt(2)*0.5*sqrt(temp*r_mag_inv);

  temp = sqrt(rb[0]*rb[0]+rb[1]*rb[1]);
  if (temp != 0.0) {
    //m[1] = -sqrt(2)*0.5*sqrt((r_mag - rb[2])*r_mag_inv)/temp;
    m[1] = -sqrt(2)*0.5/temp;
    temp = r_mag - rb[2];
    if (temp < 0.0) temp = 0.0;
    m[1] *= sqrt(temp*r_mag_inv);    
    m[2] = -m[1];
    m[1] *= rb[1];
    m[2] *= rb[0];    
  } else {
    // If aligned along z axis, x,y terms zero (r_mag-rb[2] = 0)
    m[1] = 0.0;
    m[2] = 0.0;   
  }
  m[3] = 0.0;  

  // qp21 = opposite of r^\circ_21 in Wang
  // q21 = opposite of r_21 in Wang
  MathExtra::quatquat(q2inv, q1, qp21);
  MathExtra::qconjugate(m, minv);
  MathExtra::quatquat(minv,qp21,qtmp);
  MathExtra::quatquat(qtmp,m,q21);
  
  temp = sqrt(q21[0]*q21[0] + q21[3]*q21[3]);
  if (temp != 0.0) {
    c = q21[0]/temp;
    psi = 2.0*acos_limit(c);
  } else {
    c = 0.0;
    psi = 0.0;
  }
  
  // Map negative rotations
  if(q21[3] < 0.0) // sin = q21[3]/temp
    psi = -psi;
    
  if(q21[3] == 0.0)
    psi = 0.0;
  
  c = q21[0]*q21[0] - q21[1]*q21[1] - q21[2]*q21[2] + q21[3]*q21[3];
  theta = acos_limit(c);
  
  // Separately calculte magnitude of quaternion in x-y and out of x-y planes
  // to avoid dividing by zero
  mag_out_plane = (q21[0]*q21[0] + q21[3]*q21[3]);
  mag_in_plane = (q21[1]*q21[1] + q21[2]*q21[2]);
    
  if (mag_in_plane == 0.0) {
    // No rotation => no bending/shear torque or extra shear force
    // achieve by setting cos/sin = 0
    cos_phi = 0.0;
    sin_phi = 0.0;
  } else if (mag_out_plane == 0.0) {
    // Calculate angle in plane
    cos_phi =  q21[2]/sqrt(mag_in_plane);
    sin_phi = -q21[1]/sqrt(mag_in_plane);
  } else {
    // Default equations in Mora, Wang 2009
    cos_phi = q21[1]*q21[3] + q21[0]*q21[2];
    sin_phi = q21[2]*q21[3] - q21[0]*q21[1];
    
    cos_phi /= sqrt(mag_out_plane*mag_in_plane);
    sin_phi /= sqrt(mag_out_plane*mag_in_plane);
  }
  
  Tbp[0] = -Kb[type]*theta*sin_phi;
  Tbp[1] = Kb[type]*theta*cos_phi;
  Tbp[2] = 0.0;
  
  Ttp[0] = 0.0;
  Ttp[1] = 0.0;
  Ttp[2] = Kt[type]*psi;
  
  Fsp[0] = -0.5*Ks[type]*r_mag*theta*cos_phi;
  Fsp[1] = -0.5*Ks[type]*r_mag*theta*sin_phi;
  Fsp[2] = 0.0;
  
  Tsp[0] = 0.25*Ks[type]*r_mag*theta*sin_phi;
  Tsp[1] = -0.25*Ks[type]*r_mag*r_mag*theta*cos_phi;
  Tsp[2] = 0.0;
  
  // Rotate forces/torques back to 1st particle's frame
  
  MathExtra::quatrotvec(m, Fsp, Ftmp);
    
  Fs[0] += Ftmp[0];
  Fs[1] += Ftmp[1];
  Fs[2] += Ftmp[2];

  MathExtra::quatrotvec(m, Tbp, Ttmp);
  Tb[0] = Ttmp[0];
  Tb[1] = Ttmp[1];
  Tb[2] = Ttmp[2];
  
  MathExtra::quatrotvec(m, Ttp, Ttmp);
  Tt[0] = Ttmp[0];
  Tt[1] = Ttmp[1];
  Tt[2] = Ttmp[2];

  MathExtra::quatrotvec(m, Tsp, Ttmp);
  Ts[0] += Ttmp[0];
  Ts[1] += Ttmp[1];
  Ts[2] += Ttmp[2];

  // Sum forces and calculate magnitudes
  F_rot[0] += Fs[0];
  F_rot[1] += Fs[1];
  F_rot[2] += Fs[2];
  
  MathExtra::quatrotvec(q2, F_rot, force1on2);  
  
  T_rot[0] = Ts[0] + Tt[0] + Tb[0];;
  T_rot[1] = Ts[1] + Tt[1] + Tb[1];;
  T_rot[2] = Ts[2] + Tt[2] + Tb[2];;

  MathExtra::quatrotvec(q2, T_rot, torque1on2);  

  T_rot[0] = Ts[0] - Tt[0] - Tb[0];;
  T_rot[1] = Ts[1] - Tt[1] - Tb[1];;
  T_rot[2] = Ts[2] - Tt[2] - Tb[2];;

  MathExtra::quatrotvec(q2, T_rot, torque2on1);  

  Fs_mag = MathExtra::len3(Fs);
  Tt_mag = MathExtra::len3(Tt);
  Tb_mag = MathExtra::len3(Tb);    
}

/* ---------------------------------------------------------------------- */

void BondDEMGeneric::compute(int eflag, int vflag)
{
    
  if(not fix_bond_store->stored_flag){
    fix_bond_store->stored_flag = true;
    store_data();   
  }      
    
  int i1,i2,m,n,type,itype,jtype;
  double evdwl,fpair,rsq,ebond;
  double q1[4], q2[4], r[3], r0[3];
  double r0_mag, r_mag, r_mag_inv, Fr_mag, Fs_mag;
  double Tt_mag, Tb_mag;
  double force1on2[3], torque1on2[3], torque2on1[3];
  double breaking, wd, wd2;
  double r2inv, wr[3], wrn[3], wrt[3], tdamp, dot;
  double fdrag, delvx, delvy, delvz;
  double radi1, radi2;
  double tor1, tor2, tor3, fs1, fs2, fs3;
  
  ev_init(eflag,vflag);

  if (vflag_global == 2)
    force->pair->vflag_either = force->pair->vflag_global = 1;

  double **cutsq = force->pair->cutsq;
  double **x = atom->x;
  double **v = atom->v;                                                             
  double **omega = atom->omega; 
  double **f = atom->f;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double **quat = atom->quat;
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
    r0_mag = bondstore[n][0];
    r0[0] = bondstore[n][1];
    r0[1] = bondstore[n][2];
    r0[2] = bondstore[n][3];
    
    //Stored nb points towards small tag
    //Reexpress so points towards atom i1
    if(tag[i2] < tag[i1]){
      r0[0] = -r0[0];    
      r0[1] = -r0[1];
      r0[2] = -r0[2];      
    }    
    
    q1[0] = quat[i1][0];
    q1[1] = quat[i1][1]; 
    q1[2] = quat[i1][2];
    q1[3] = quat[i1][3];
    
    q2[0] = quat[i2][0];
    q2[1] = quat[i2][1];
    q2[2] = quat[i2][2];
    q2[3] = quat[i2][3];
    
    // Note this is the reverse of Mora & Wang
    r[0] = x[i1][0] - x[i2][0];
    r[1] = x[i1][1] - x[i2][1];
    r[2] = x[i1][2] - x[i2][2];

    rsq = MathExtra::lensq3(r);
    r_mag = sqrt(rsq);
    
    // Calculate forces from Wang 2009
    calc_forces(type, r_mag, r0_mag, q1, q2, r0, r, force1on2, torque1on2, torque2on1, Fs_mag, Fr_mag, Tb_mag, Tt_mag);
    
    breaking = Fr_mag/Fcr[type] + Fs_mag/Fcs[type] + Tb_mag/Gcb[type] + Tt_mag/Gct[type];

    if (breaking >= 1.0) {
      bondlist[n][2] = 0;
      for (m = 0; m < atom->num_bond[i1]; m++)
        if (atom->bond_atom[i1][m] == atom->tag[i2])
          atom->bond_type[i1][m] = 0;
      if (i2 < atom->nlocal)
        for (m = 0; m < atom->num_bond[i2]; m++)
          if (atom->bond_atom[i2][m] == atom->tag[i1])
            atom->bond_type[i2][m] = 0;
        
      if(fix_broken_bonds != NULL) fix_broken_bonds->add_bond(i1, i2);          
        
      continue;
    }

    wd = 1.0 - breaking;     
    wd2 = wd*wd;

    if (eflag) ebond = -wd; //Not really, but it's just to give a sense

    // Damping relative tangental rotation
    wr[0] = (omega[i1][0] + omega[i2][0]);
    wr[1] = (omega[i1][1] + omega[i2][1]);
    wr[2] = (omega[i1][2] + omega[i2][2]);    
    
    // Subtract out normal component
    dot = MathExtra::dot3(wr, r);
    dot /= r_mag;
    
    wrt[0] = wr[0] - dot*r[0]; 
    wrt[1] = wr[1] - dot*r[1];
    wrt[2] = wr[2] - dot*r[2];
    
    torque2on1[0] -= gammaw[type]*wrt[0];
    torque2on1[1] -= gammaw[type]*wrt[1];
    torque2on1[2] -= gammaw[type]*wrt[2];
    
    torque1on2[0] -= gammaw[type]*wrt[0];
    torque1on2[1] -= gammaw[type]*wrt[1];
    torque1on2[2] -= gammaw[type]*wrt[2];
    
    // Damping relative normal component, twist
    wr[0] = (omega[i1][0] - omega[i2][0]);
    wr[1] = (omega[i1][1] - omega[i2][1]);
    wr[2] = (omega[i1][2] - omega[i2][2]);   
    
    dot = MathExtra::dot3(wr, r);
    dot /= r_mag;
    wrn[0] = dot*r[0]; 
    wrn[1] = dot*r[1];
    wrn[2] = dot*r[2];    
    
    torque1on2[0] += gammaw[type]*wrn[0];
    torque1on2[1] += gammaw[type]*wrn[1];
    torque1on2[2] += gammaw[type]*wrn[2];
            
    torque2on1[0] -= gammaw[type]*wrn[0];
    torque2on1[1] -= gammaw[type]*wrn[1];
    torque2on1[2] -= gammaw[type]*wrn[2];   

    // DPD Damping force
    r2inv = 1.0/rsq;    
    delvx = v[i1][0] - v[i2][0];                                       
    delvy = v[i1][1] - v[i2][1];                                       
    delvz = v[i1][2] - v[i2][2];     
    dot = r[0]*delvx + r[1]*delvy + r[2]*delvz;
    fdrag = -gamma[type]*dot*r2inv;

    force1on2[0] -= r[0]*fdrag;
    force1on2[1] -= r[1]*fdrag;
    force1on2[2] -= r[2]*fdrag;

    if (newton_bond || i1 < nlocal) {

      f[i1][0] -= force1on2[0]*wd2;
      f[i1][1] -= force1on2[1]*wd2;
      f[i1][2] -= force1on2[2]*wd2;   
      
      torque[i1][0] += torque2on1[0]*wd2;
      torque[i1][1] += torque2on1[1]*wd2; 
      torque[i1][2] += torque2on1[2]*wd2;
    }

    if (newton_bond || i2 < nlocal) {
      
      f[i2][0] += force1on2[0]*wd2;
      f[i2][1] += force1on2[1]*wd2;
      f[i2][2] += force1on2[2]*wd2; 

      torque[i2][0] += torque1on2[0]*wd2;
      torque[i2][1] += torque1on2[1]*wd2; 
      torque[i2][2] += torque1on2[2]*wd2;      
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,Fr_mag,r[0],r[1],r[2]);

    // subtract out pairwise contribution from 2 atoms via pair->single()

    // Skip if overlaying
    if(overlay_pair_flag) continue;

    itype = atom->type[i1];
    jtype = atom->type[i2];
    
    if (rsq < cutsq[itype][jtype]) {
      evdwl = -force->pair->single(i1,i2,itype,jtype,rsq,1.0,1.0,fpair);
      fpair = -fpair;
      
      r_mag_inv = 1/r_mag;
      fs1 = -force->pair->svector[0];
      fs2 = -force->pair->svector[1];
      fs3 = -force->pair->svector[2];
      radi1 = radius[i1];
      radi2 = radius[i2];
      tor1 = r_mag_inv * (r[1]*fs3 - r[2]*fs2);
      tor2 = r_mag_inv * (r[2]*fs1 - r[0]*fs3);
      tor3 = r_mag_inv * (r[0]*fs2 - r[1]*fs1);
        
      //Note only works with my modified pair granular for now      
      if (newton_bond || i1 < nlocal) {
        f[i1][0] += r[0]*fpair + fs1;
        f[i1][1] += r[1]*fpair + fs2;
        f[i1][2] += r[2]*fpair + fs3;          
        torque[i1][0] -= radi1*tor1;
        torque[i1][1] -= radi1*tor2;
        torque[i1][2] -= radi1*tor3;
      }
      if (newton_bond || i2 < nlocal) {
        f[i2][0] -= r[0]*fpair + fs1;
        f[i2][1] -= r[1]*fpair + fs2;
        f[i2][2] -= r[2]*fpair + fs3;          
        torque[i2][0] -= radi2*tor1;
        torque[i2][1] -= radi2*tor2;
        torque[i2][2] -= radi2*tor3;
      }
      
      if (evflag) force->pair->ev_tally(i1,i2,nlocal,newton_bond,
                                        evdwl,0.0,fpair,r[0],r[1],r[2]);
    }
  }
}

/* ---------------------------------------------------------------------- */

void BondDEMGeneric::settings(int narg, char **arg)
{
  int iarg = 0;
  while(iarg < narg){
    if (strcmp(arg[iarg], "overlay/pair") == 0){
      overlay_pair_flag = 1;
    } else error->all(FLERR,"Illegal pair_style command");
    iarg++;   
  }
}

/* ---------------------------------------------------------------------- */

void BondDEMGeneric::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(Kr,n+1,"bond:Kr");
  memory->create(Ks,n+1,"bond:Ks");
  memory->create(Kt,n+1,"bond:Kt");
  memory->create(Kb,n+1,"bond:Kb");
  memory->create(Fcr,n+1,"bond:Fcr");
  memory->create(Fcs,n+1,"bond:Fcs");
  memory->create(Gct,n+1,"bond:Gct");
  memory->create(Gcb,n+1,"bond:Gcb");
  memory->create(gamma,n+1,"bond:gamma");
  memory->create(gammaw,n+1,"bond:gammaw");
  memory->create(C_exp,n+1,"bond:C_exp");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;

}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondDEMGeneric::coeff(int narg, char **arg)
{
  if (narg != 12) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nbondtypes,ilo,ihi,error);

  double Kr_one = utils::numeric(FLERR,arg[1],false,lmp);
  double Ks_one = utils::numeric(FLERR,arg[2],false,lmp);
  double Kt_one = utils::numeric(FLERR,arg[3],false,lmp);
  double Kb_one = utils::numeric(FLERR,arg[4],false,lmp);
  double Fcr_one = utils::numeric(FLERR,arg[5],false,lmp);
  double Fcs_one = utils::numeric(FLERR,arg[6],false,lmp);
  double Gct_one = utils::numeric(FLERR,arg[7],false,lmp);
  double Gcb_one = utils::numeric(FLERR,arg[8],false,lmp);
  double gamma_one = utils::numeric(FLERR,arg[9],false,lmp);
  double gammaw_one = utils::numeric(FLERR,arg[10],false,lmp);
  double Cexp_one = utils::numeric(FLERR,arg[11],false,lmp);
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    Kr[i] = Kr_one;
    Ks[i] = Ks_one;
    Kt[i] = Kt_one;
    Kb[i] = Kb_one;
    Fcr[i] = Fcr_one;
    Fcs[i] = Fcs_one;
    Gct[i] = Gct_one;
    Gcb[i] = Gcb_one;
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

void BondDEMGeneric::init_style()
{
  if (!atom->quat_flag)
    error->all(FLERR,"Bond generic requires atom attributes quaternion");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Bond generic requires ghost atoms store velocity");

  if (force->pair == NULL || force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support bond_style generic");

  if (force->angle || force->dihedral || force->improper)
    error->all(FLERR,
               "Bond style generic cannot be used with 3,4-body interactions");
  if (atom->molecular == 2)
    error->all(FLERR,
               "Bond style generic cannot be used with atom style template");

  if(domain->dimension == 2)
    error->warning(FLERR, "Bond style generic not intended for 2d use, may be inefficient");
  // Determine if correct pair style is used
  if(not overlay_pair_flag) {
    int correct_pair = 0;
    if(force->pair_match("gran/hooke",0)) correct_pair = 1;
    if(force->pair_match("gran/hooke/history",0)) correct_pair = 1;
    if(force->pair_match("gran/hertz",0)) correct_pair = 1;
    if(force->pair_match("gran/hertz/history",0)) correct_pair = 1;
    if (! correct_pair) 
      error->all(FLERR, "Bond style DEM generic requires gran pairstyle");
  }    
    
  //Define bond store
  if(fix_bond_store == NULL){
    char **fixarg = new char*[5];
    fixarg[0] = (char *) "BOND_STORE_DEM_GENERIC";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "BOND_STORE";
    fixarg[3] = (char *) "0";
    fixarg[4] = (char *) "4";
    modify->add_fix(5,fixarg,1);
    delete [] fixarg;
    int ifix = modify->find_fix("BOND_STORE_DEM_GENERIC");
    fix_bond_store = (FixBondStore *) modify->fix[ifix];
    //Note don't use most recent nfix b/c fix bond store creates a fix property atom    
  }
  
  int ifix = modify->find_fix_by_style("bonds/broken");
  if(ifix != -1) fix_broken_bonds = (FixBrokenBonds *) modify->fix[ifix];      
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondDEMGeneric::equilibrium_distance(int i)
{
  return max_r0;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondDEMGeneric::write_restart(FILE *fp)
{
  fwrite(&Kr[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&Ks[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&Kt[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&Kb[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&Fcr[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&Fcs[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&Gct[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&Gcb[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&gamma[1],sizeof(double),atom->nbondtypes,fp);  
  fwrite(&gammaw[1],sizeof(double),atom->nbondtypes,fp);  
  fwrite(&C_exp[1],sizeof(double),atom->nbondtypes,fp);  
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondDEMGeneric::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&Kr[1],sizeof(double),atom->nbondtypes,fp);
    fread(&Ks[1],sizeof(double),atom->nbondtypes,fp);
    fread(&Kt[1],sizeof(double),atom->nbondtypes,fp);
    fread(&Kb[1],sizeof(double),atom->nbondtypes,fp);
    fread(&Fcr[1],sizeof(double),atom->nbondtypes,fp);
    fread(&Fcs[1],sizeof(double),atom->nbondtypes,fp);
    fread(&Gct[1],sizeof(double),atom->nbondtypes,fp);
    fread(&Gcb[1],sizeof(double),atom->nbondtypes,fp);
    fread(&gamma[1],sizeof(double),atom->nbondtypes,fp);    
    fread(&gammaw[1],sizeof(double),atom->nbondtypes,fp);    
    fread(&C_exp[1],sizeof(double),atom->nbondtypes,fp);    
  }
  MPI_Bcast(&Kr[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Ks[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Kt[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Kb[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Fcr[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Fcs[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Gct[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Gcb[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&gamma[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&gammaw[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&C_exp[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  
  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondDEMGeneric::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
     fprintf(fp,"%d %g %g %g %g %g %g %g %g %g %g %g\n",i,Kr[i],Ks[i],Kt[i],Kb[i],Fcr[i], Fcs[i], Gct[i], Gcb[i], gamma[i],gammaw[i], C_exp[i]);
}

/* ---------------------------------------------------------------------- */
// Not implemented

double BondDEMGeneric::single(int type, double rsq, int i, int j,
                           double &fforce)
{
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
  return 0.0;
}
