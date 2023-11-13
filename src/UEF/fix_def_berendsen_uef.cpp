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
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fix_def_berendsen_uef.h"
#include "compute_pressure_def_uef.h"
#include "compute_temp_def_uef.h"
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "compute.h"
#include "irregular.h"
#include "domain.h"
#include "lattice.h"
#include "force.h"
#include "modify.h"
#include "math_const.h"
#include "input.h"
#include "error.h"
#include "output.h"
#include "timer.h"
#include "neighbor.h"
#include "uef_utils.h"
#include "neigh_list.h"
#include "pair.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;



FixDefBerendsenUEF::FixDefBerendsenUEF(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
irregular(NULL), uefbox(NULL)
{
  if (narg != 8 && narg != 10) error->all(FLERR,"Illegal fix uef deform command");

  no_change_box = 1;    //Defined in deform and in nh
                        //1 => cannot swap ortho <-> triclinic using change_box while this fix is running
  restart_global = 1;   //Defined in deform and in nh
                        //1 => fix saves global state -not entirely clear used in modify.cpp
  pre_exchange_migrate = 1;
  vector_flag = 1;
  size_vector = 12;

  rate[0] = utils::numeric(FLERR,arg[3],false,lmp);
  rate[1] = utils::numeric(FLERR,arg[4],false,lmp);

  p_start = utils::numeric(FLERR,arg[5],false,lmp);
  p_stop = utils::numeric(FLERR,arg[6],false,lmp);
  p_period = utils::numeric(FLERR,arg[7],false,lmp);

  // initialize the UEFBox class which computes the box at each step
  uefbox = new UEF_utils::UEFBox();

  if (narg == 10) {
    strain[0] = utils::numeric(FLERR,arg[8],false,lmp);
    strain[1] = utils::numeric(FLERR,arg[9],false,lmp);
    strain[2] = -strain[0] - strain[1];
    uefbox->set_strain(strain[0],strain[1]);
  } else {
    strain[0] = 0;
    strain[1] = 0;
    strain[2] = 0;
  }

  box_change |= BOX_CHANGE_X;
  box_change |= BOX_CHANGE_Y;
  box_change |= BOX_CHANGE_Z;
  box_change |= BOX_CHANGE_YZ;
  box_change |= BOX_CHANGE_XZ;
  box_change |= BOX_CHANGE_XY;


  irregular = new Irregular(lmp);
  force_reneighbor = 1;
  next_reneighbor = -1;
  nevery = 1;

  scalar_flag = 1;
  extscalar = 0;
  flip = 0;

  // create a new compute temp style

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} all temp/def/uef",id_temp));
  tflag = 1;

  // create a new compute pressure style

  id_press = utils::strdup(std::string(id) + "_press");
  modify->add_compute(fmt::format("{} all pressure/def/uef {}",id_press, id_temp));
  pflag = 1;

}

/* ---------------------------------------------------------------------- */

FixDefBerendsenUEF::~FixDefBerendsenUEF()
{
  delete irregular;
  delete uefbox;
  if (tflag) modify->delete_compute(id_temp);
  if (pflag) modify->delete_compute(id_press);
  delete [] id_temp;
  delete [] id_press;
}

/* ---------------------------------------------------------------------- */

int FixDefBerendsenUEF::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDefBerendsenUEF::init()
{
  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"deform") == 0) count++;
  if (count > 1) error->all(FLERR,"More than one fix deform");

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix press/berendsen does not exist");
  temperature = modify->compute[icompute];

  icompute = modify->find_compute(id_press);
  if (icompute < 0)
    error->all(FLERR,"Pressure ID for fix press/berendsen does not exist");
  pressure = modify->compute[icompute];
}

/* ----------------------------------------------------------------------
 * Run setup() make sure the box is OK and set the rotation matrix
 * for the first step
 * ---------------------------------------------------------------------- */

void FixDefBerendsenUEF::setup(int j)
{
  double box[3][3];
  double vol = domain->xprd * domain->yprd * domain->zprd;
  uefbox->get_box(box,vol);
  double tol = 1e-4; // too many box resizes may cause drift
  // ensure the box is ok for uef
  bool isok = true;
  isok &= nearly_equal(domain->h[0],box[0][0],tol);
  isok &= nearly_equal(domain->h[1],box[1][1],tol);
  isok &= nearly_equal(domain->h[2],box[2][2],tol);
  isok &= nearly_equal(domain->xy,box[0][1],tol);
  isok &= nearly_equal(domain->yz,box[1][2],tol);
  isok &= nearly_equal(domain->xz,box[0][2],tol);

  if (!isok)
    error->all(FLERR,"Initial box is not close enough to the expected uef box");

  uefbox->get_rot(rot);

  // trigger virial computation on next timestep

  ((ComputeTempDefUef*) temperature)->yes_rot();
  ((ComputePressureDefUef*) pressure)->in_fix = true;
  ((ComputePressureDefUef*) pressure)->update_rot();

  pressure->addstep(update->ntimestep+1);
}

/* ----------------------------------------------------------------------
  box flipped on previous step
  perform irregular on atoms in lamda coords to migrate atoms to new procs
------------------------------------------------------------------------- */

void FixDefBerendsenUEF::pre_exchange()
{
  if (flip != 0) {

    // go to lab frame
    inv_rotate_x(rot);
    inv_rotate_v(rot);
    inv_rotate_f(rot);


    // get & set the new box and rotation matrix
    double vol = domain->xprd * domain->yprd * domain->zprd;
    double box[3][3];
    uefbox->get_box(box,vol);
    domain->boxhi[0] = domain->boxlo[0]+box[0][0];
    domain->boxhi[1] = domain->boxlo[1]+box[1][1];
    domain->boxhi[2] = domain->boxlo[2]+box[2][2];
    domain->xy = box[0][1];
    domain->xz = box[0][2];
    domain->yz = box[1][2];
    domain->set_global_box();
    domain->set_local_box();
    uefbox->get_rot(rot);

    // rotate to the new upper triangular frame
    rotate_v(rot);
    rotate_x(rot);
    rotate_f(rot);

    // put all atoms in the new box
    double **x = atom->x;
    imageint *image = atom->image;
    int nlocal = atom->nlocal;
    for (int i=0; i<nlocal; i++) domain->remap(x[i],image[i]);

    // move atoms to the right processors
    domain->x2lamda(atom->nlocal);
    irregular->migrate_atoms();
    domain->lamda2x(atom->nlocal);

    flip = 0;
  }
}

/* ----------------------------------------------------------------------
 * The following are routines to rotate between the lab and upper triangular
 * (UT) frames. For most of the time the simulation is in the UT frame.
 * To get to the lab frame, apply the inv_rotate_[..](rot) and to
 * get back to the UT frame apply rotate_[..](rot).
 *
 * Note: the rotate_x() functions also apply a shift to/from the fixedpoint
 * to make the integration a little simpler.
 * ---------------------------------------------------------------------- */

void FixDefBerendsenUEF::rotate_x(double r[3][3])
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double xn[3];
  for (int i=0;i<nlocal;i++)
  {
    if (mask[i] & groupbit)
    {
      xn[0]=r[0][0]*x[i][0]+r[0][1]*x[i][1]+r[0][2]*x[i][2];
      xn[1]=r[1][0]*x[i][0]+r[1][1]*x[i][1]+r[1][2]*x[i][2];
      xn[2]=r[2][0]*x[i][0]+r[2][1]*x[i][1]+r[2][2]*x[i][2];
      x[i][0]=xn[0]+domain->boxlo[0];
      x[i][1]=xn[1]+domain->boxlo[1];
      x[i][2]=xn[2]+domain->boxlo[2];
    }
  }
}

void FixDefBerendsenUEF::inv_rotate_x(double r[3][3])
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double xn[3];
  for (int i=0;i<nlocal;i++)
  {
    if (mask[i] & groupbit)
    {
      x[i][0] -= domain->boxlo[0];
      x[i][1] -= domain->boxlo[1];
      x[i][2] -= domain->boxlo[2];
      xn[0]=r[0][0]*x[i][0]+r[1][0]*x[i][1]+r[2][0]*x[i][2];
      xn[1]=r[0][1]*x[i][0]+r[1][1]*x[i][1]+r[2][1]*x[i][2];
      xn[2]=r[0][2]*x[i][0]+r[1][2]*x[i][1]+r[2][2]*x[i][2];
      x[i][0]=xn[0];
      x[i][1]=xn[1];
      x[i][2]=xn[2];
    }
  }
}

void FixDefBerendsenUEF::rotate_v(double r[3][3])
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double vn[3];
  for (int i=0;i<nlocal;i++)
  {
    if (mask[i] & groupbit)
    {
      vn[0]=r[0][0]*v[i][0]+r[0][1]*v[i][1]+r[0][2]*v[i][2];
      vn[1]=r[1][0]*v[i][0]+r[1][1]*v[i][1]+r[1][2]*v[i][2];
      vn[2]=r[2][0]*v[i][0]+r[2][1]*v[i][1]+r[2][2]*v[i][2];
      v[i][0]=vn[0]; v[i][1]=vn[1]; v[i][2]=vn[2];
    }
  }
}

void FixDefBerendsenUEF::inv_rotate_v(double r[3][3])
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double vn[3];
  for (int i=0;i<nlocal;i++)
  {
    if (mask[i] & groupbit)
    {
      vn[0]=r[0][0]*v[i][0]+r[1][0]*v[i][1]+r[2][0]*v[i][2];
      vn[1]=r[0][1]*v[i][0]+r[1][1]*v[i][1]+r[2][1]*v[i][2];
      vn[2]=r[0][2]*v[i][0]+r[1][2]*v[i][1]+r[2][2]*v[i][2];
      v[i][0]=vn[0]; v[i][1]=vn[1]; v[i][2]=vn[2];
    }
  }
}

void FixDefBerendsenUEF::rotate_f(double r[3][3])
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double fn[3];
  for (int i=0;i<nlocal;i++)
  {
    if (mask[i] & groupbit)
    {
      fn[0]=r[0][0]*f[i][0]+r[0][1]*f[i][1]+r[0][2]*f[i][2];
      fn[1]=r[1][0]*f[i][0]+r[1][1]*f[i][1]+r[1][2]*f[i][2];
      fn[2]=r[2][0]*f[i][0]+r[2][1]*f[i][1]+r[2][2]*f[i][2];
      f[i][0]=fn[0]; f[i][1]=fn[1]; f[i][2]=fn[2];
    }
  }
}

void FixDefBerendsenUEF::inv_rotate_f(double r[3][3])
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  double fn[3];
  for (int i=0;i<nlocal;i++)
  {
    if (mask[i] & groupbit)
    {
      fn[0]=r[0][0]*f[i][0]+r[1][0]*f[i][1]+r[2][0]*f[i][2];
      fn[1]=r[0][1]*f[i][0]+r[1][1]*f[i][1]+r[2][1]*f[i][2];
      fn[2]=r[0][2]*f[i][0]+r[1][2]*f[i][1]+r[2][2]*f[i][2];
      f[i][0]=fn[0]; f[i][1]=fn[1]; f[i][2]=fn[2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixDefBerendsenUEF::end_of_step()
{
  double iv = domain->xprd*domain->yprd*domain->zprd;
  double dtv = update->dt;
  double ex = rate[0]*dtv;
  double ey = rate[1]*dtv;
  strain[0] += ex;
  strain[1] += ey;
  strain[2] += -ex-ey;

  int i;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      domain->x2lamda(x[i],x[i]);

  //Run after NVE_X (in initial_integrate) can move to post_integrate
  uefbox->step_deform(ex,ey);
  double box[3][3];
  double vol = domain->xprd * domain->yprd * domain->zprd;
  uefbox->get_box(box,vol);
  domain->boxhi[0] = domain->boxlo[0]+box[0][0];
  domain->boxhi[1] = domain->boxlo[1]+box[1][1];
  domain->boxhi[2] = domain->boxlo[2]+box[2][2];
  domain->xy = box[0][1];
  domain->xz = box[0][2];
  domain->yz = box[1][2];


 // compute new T,P
  ((ComputePressureDefUef*) pressure)->update_rot();
  temperature->compute_vector();
  pressure->compute_vector();
  double *tensor = pressure->vector;

  p_current = 0.5 * (tensor[0] + tensor[1]);
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  p_target = p_start + delta * (p_stop-p_start);
  dilation = pow(1.0 - update->dt/p_period * (p_target-p_current),1.0/3.0);
  double oldlo,oldhi,ctr;

  for (i = 0; i < 2; i++) {
      oldlo = domain->boxlo[i];
      oldhi = domain->boxhi[i];
      ctr = 0.5 * (oldlo + oldhi);
      domain->boxlo[i] = (oldlo-ctr)*dilation + ctr;
      domain->boxhi[i] = (oldhi-ctr)*dilation + ctr;
  }

  domain->set_global_box();
  domain->set_local_box();
  uefbox->get_rot(rot);


  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      domain->lamda2x(x[i],x[i]);

  if(uefbox->reduce()) flip = 1;
  if(flip) next_reneighbor = update->ntimestep + 1;

  // trigger virial computation on next timestep
  pressure->addstep(update->ntimestep+1);
}


/* ----------------------------------------------------------------------
   write Set data to restart file
------------------------------------------------------------------------- */

void FixDefBerendsenUEF::write_restart(FILE *fp)
{
  if (comm->me == 0) {
    int size = 3*sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(strain,sizeof(double),3,fp);
  }
}

/* ----------------------------------------------------------------------
   use selected state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixDefBerendsenUEF::restart(char *buf)
{

  double *strain_restart = (double *) buf;
  uefbox->set_strain(strain_restart[0],strain_restart[1]);
  strain[0] = strain_restart[0];
  strain[1] = strain_restart[1];
}

/* ----------------------------------------------------------------------
 * public read for rotation matrix
 * ---------------------------------------------------------------------- */

void FixDefBerendsenUEF::get_rot(double r[3][3])
{
  r[0][0] = rot[0][0];
  r[0][1] = rot[0][1];
  r[0][2] = rot[0][2];
  r[1][0] = rot[1][0];
  r[1][1] = rot[1][1];
  r[1][2] = rot[1][2];
  r[2][0] = rot[2][0];
  r[2][1] = rot[2][1];
  r[2][2] = rot[2][2];
}

/* ----------------------------------------------------------------------
 * public read for simulation box
 * ---------------------------------------------------------------------- */
void FixDefBerendsenUEF::get_box(double b[3][3])
{
  double box[3][3];
  double vol = domain->xprd * domain->yprd * domain->zprd;
  uefbox->get_box(box,vol);
  b[0][0] = box[0][0];
  b[0][1] = box[0][1];
  b[0][2] = box[0][2];
  b[1][0] = box[1][0];
  b[1][1] = box[1][1];
  b[1][2] = box[1][2];
  b[2][0] = box[2][0];
  b[2][1] = box[2][1];
  b[2][2] = box[2][2];
}

/* ----------------------------------------------------------------------
 * comparing floats
 * it's imperfect, but should work provided no infinities
 * ---------------------------------------------------------------------- */
bool FixDefBerendsenUEF::nearly_equal(double a, double b, double epsilon)
{
  double absa = fabs(a);
  double absb = fabs(b);
  double diff = fabs(a-b);
  if (a == b) return true;
  else if ( (absa+absb) < epsilon)
    return diff < epsilon*epsilon;
  else
    return diff/(absa+absb) < epsilon;
}

double FixDefBerendsenUEF::compute_scalar()
{
  return strain[0];
}

double FixDefBerendsenUEF::compute_vector(int n)
{
  uefbox->get_rot(rot);

  if(n == 0) return rot[0][0];
  if(n == 1) return rot[0][1];
  if(n == 2) return rot[0][2];

  if(n == 3) return rot[1][0];
  if(n == 4) return rot[1][1];
  if(n == 5) return rot[1][2];

  if(n == 6) return rot[2][0];
  if(n == 7) return rot[2][1];
  if(n == 8) return rot[2][2];

  if(n == 9)  return strain[0];
  if(n == 10) return strain[1];
  if(n == 11) return strain[2];

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int FixDefBerendsenUEF::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    id_temp = utils::strdup(arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,"Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning(FLERR,"Temperature for NPT is not for group all");

    // reset id_temp of pressure to new temperature ID

    icompute = modify->find_compute(id_press);
    if (icompute < 0)
      error->all(FLERR,"Pressure ID for fix press/berendsen does not exist");
    modify->compute[icompute]->reset_extra_compute_fix(id_temp);

    return 2;

  } else if (strcmp(arg[0],"press") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (pflag) {
      modify->delete_compute(id_press);
      pflag = 0;
    }
    delete [] id_press;
    id_press = utils::strdup(arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify pressure ID");
    pressure = modify->compute[icompute];

    if (pressure->pressflag == 0)
      error->all(FLERR,"Fix_modify pressure ID does not compute pressure");
    return 2;
  }
  return 0;
}
