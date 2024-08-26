// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Dan Bolintineanu (SNL), Joel Clemmer (SNL), Ishan Srivastava (SNL),
   Jeremy Lechman(SNL), Leo Silbert (SNL), Gary Grest (SNL)
----------------------------------------------------------------------- */

// TODO:
//       move init and coeff to local to function
//       single()
//       doc files
//       test

#include "pair_granular_accelerated.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "fix_dummy.h"
#include "fix_neigh_history.h"
#include "force.h"
#include "granular_model.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace MathExtra;

using MathConst::MY_2PI;
using MathConst::MY_PI;

static constexpr double PI27SQ = 266.47931882941264802866;      // 27*PI**2
static constexpr double THREEROOT3 = 5.19615242270663202362;    // 3*sqrt(3)
static constexpr double SIXROOT6 = 14.69693845669906728801;     // 6*sqrt(6)
static constexpr double INVROOT6 = 0.40824829046386307274;      // 1/sqrt(6)
static constexpr double FOURTHIRDS = (4.0 / 3.0);               // 4/3
static constexpr double JKRPREFIX = 1.2277228507842888;         // cbrt(3*PI**2/16)
static constexpr double EPSILON = 1e-10;

/* ---------------------------------------------------------------------- */

PairGranularAccelerated::PairGranularAccelerated(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  no_virial_fdotr_compute = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  finitecutflag = 1;

  single_extra = 12;
  svector = new double[single_extra];

  neighprev = 0;
  nmax = 0;
  mass_rigid = nullptr;

  onerad_dynamic = nullptr;
  onerad_frozen = nullptr;
  maxrad_dynamic = nullptr;
  maxrad_frozen = nullptr;

  // --------------------------------------------------------------- //
  //  CUSTOMIZE: variables from constructors of GranSubMod classes
  //   Normal, damping, and tangential are mandatory
  //   Rolling, twisting, and heat are optional (otherwise use "none")
  //   Must manually confirm submodels are compatable, i.e.:
  //     normal submodel specifies material properties, if needed
  //     some submodels don't allow cohesive normal models
  // --------------------------------------------------------------- //

  contact_radius_flag = 1;            // If set in any submodel
  nondefault_history_transfer = 1;    // If set in any submodel

  model_name[NORMAL] = "hertz/material";
  model_size_history[NORMAL] = 0;
  model_num_coeffs[NORMAL] = 3;
  beyond_contact = 0;

  model_name[DAMPING] = "tsuji";
  model_size_history[DAMPING] = 0;
  model_num_coeffs[DAMPING] = 0;

  model_name[TANGENTIAL] = "mindlin_rescale/force";
  model_size_history[TANGENTIAL] = 4;
  model_num_coeffs[TANGENTIAL] = 3;
  mindlin_force = 1;                  // Only defined in mindlin classes
  mindlin_rescale = 1;                // Only defined in mindlin classes

  model_name[ROLLING] = "sds";
  model_size_history[ROLLING] = 3;
  model_num_coeffs[ROLLING] = 3;

  model_name[TWISTING] = "sds";
  model_size_history[TWISTING] = 3;
  model_num_coeffs[TWISTING] = 3;

  model_name[HEAT] = "radius";
  model_size_history[HEAT] = 0;
  model_num_coeffs[HEAT] = 1;

  // -------------------------------------------------------------- //

  if (model_name[ROLLING] == "none")
    rolling_defined = 0;
  else
    rolling_defined = 1;

  if (model_name[TWISTING] == "none")
    twisting_defined = 0;
  else
    twisting_defined = 1;

  if (model_name[HEAT] == "none")
    heat_defined = 0;
  else
    heat_defined = 1;

  size_history = 0;
  for (int i = 0; i < NSUBMODELS; i++)
    size_history += model_size_history[i];

  transfer_history_factor = nullptr;
  if (nondefault_history_transfer) {
    transfer_history_factor = new double[size_history];
    for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1;

    double *transfer_factor_normal = transfer_history_factor;
    double *transfer_factor_damping = transfer_factor_normal + model_size_history[0];
    double *transfer_factor_tangential = transfer_factor_damping + model_size_history[1];
    double *transfer_factor_rolling = transfer_factor_tangential + model_size_history[2];
    double *transfer_factor_twisting = transfer_factor_rolling + model_size_history[3];
    double *transfer_factor_heat = transfer_factor_twisting + model_size_history[4];

    // --------------------------------------------------------------- //
    //  CUSTOMIZE: add nondefault history transfer factors from
    //    constructor, if neeeded. Use custom variables, demonstrated
    //    below, w/ same index in submodel
    // --------------------------------------------------------------- //

    //transfer_factor_normal[] = ;
    //transfer_factor_damping[] = ;
    transfer_factor_tangential[3] = +1;
    //transfer_factor_rolling[] = ;
    //transfer_factor_twisting[] = ;
    //transfer_factor_heat[] = ;

    // --------------------------------------------------------------- //
  }

  // set comm size needed by this Pair if used with fix rigid

  comm_forward = 1;

  // create dummy fix as placeholder for FixNeighHistory
  // this is so final order of Modify:fix will conform to input script

  fix_history = nullptr;
  fix_dummy = dynamic_cast<FixDummy *>(modify->add_fix("NEIGH_HISTORY_GRANULAR_DUMMY all DUMMY"));
}

/* ---------------------------------------------------------------------- */

PairGranularAccelerated::~PairGranularAccelerated()
{
  delete[] svector;

  if (!fix_history) modify->delete_fix("NEIGH_HISTORY_GRANULAR_DUMMY");
  else modify->delete_fix("NEIGH_HISTORY_GRANULAR");

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutoff_type);

    delete [] onerad_dynamic;
    delete [] onerad_frozen;
    delete [] maxrad_dynamic;
    delete [] maxrad_frozen;

    memory->destroy(model_coeffs);
  }

  memory->destroy(mass_rigid);
  delete [] transfer_history_factor;
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::compute(int eflag, int vflag)
{
  int i, j, k, ii, jj, inum, jnum, itype, jtype;
  double factor_lj, mi, mj;

  int *ilist, *jlist, *numneigh, **firstneigh;
  int *touch, **firsttouch;
  double *history, *allhistory, **firsthistory;
  int prior_itype = -1;
  int prior_jtype = -1;

  bool touchflag = false;
  const bool history_update = update->setupflag == 0;

  ev_init(eflag,vflag);

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    auto mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm(this);
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  double *heatflow, *temperature;
  if (heat_defined) {
    heatflow = atom->heatflow;
    temperature = atom->temperature;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  if (use_history) {
    firsttouch = fix_history->firstflag;
    firsthistory = fix_history->firstvalue;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    xi = x[i];
    radi = radius[i];
    mi = rmass[i];
    vi = v[i];
    omegai = omega[i];
    if (heat_defined)
      Ti = temperature[i];

    if (fix_rigid)
      if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
    if (use_history) {
      touch = firsttouch[i];
      allhistory = firsthistory[i];
    }
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      if (factor_lj == 0) continue;

      jtype = type[j];

      if (itype != prior_itype || jtype != prior_jtype) {
        prep_model(itype, jtype);
        prior_itype = itype;
        prior_jtype = jtype;
      }

      // Reset model and copy initial geometric data
      xj = x[j];
      radj = radius[j];
      if (use_history) touch = touch[jj];

      // GranularModel->check_contact()
      sub3(xi, xj, dx);
      rsq = lensq3(dx);
      radsum = radi + radj;
      Reff = radi * radj / radsum;

      touchflag = touch();

      if (!touchflag) {
        // unset non-touching neighbors
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history * jj];
          for (k = 0; k < size_history; k++) history[k] = 0.0;
        }
        continue;
      }

      // if any history is needed
      if (use_history) touch[jj] = 1;

      // meff = effective mass of pair of particles
      // if I or J part of rigid body, use body mass
      // if I or J is frozen, meff is other particle
      mj = rmass[j];
      if (fix_rigid)
        if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
      meff = mi * mj / (mi + mj);
      if (mask[i] & freeze_group_bit) meff = mj;
      if (mask[j] & freeze_group_bit) meff = mi;

      // Copy additional information and prepare force calculations
      meff = meff;
      vj = v[j];
      omegaj = omega[j];
      if (use_history)
        history = &allhistory[size_history * jj];

      if (heat_defined)
        Tj = temperature[j];

      // GranularModel->calculate_forces();

      rinv = 1.0 / r;
      delta = radsum - r;
      dR = delta * Reff;
      scale3(rinv, dx, nx);

      // relative translational velocity
      sub3(vi, vj, vr);

      // normal component
      vnnr = dot3(vr, nx);
      scale3(vnnr, nx, vn);

      // tangential component
      sub3(vr, vn, vt);

      // relative rotational velocity
      scaleadd3(radi, omegai, radj, omegaj, wr);

      // relative tangential velocities
      double temp[3];
      cross3(wr, nx, temp);
      sub3(vt, temp, vtr);
      vrel = len3(vtr);

      // calculate forces/torques
      double Fdamp, dist_to_contact;
      history_index = model_history_index[NORMAL];
      if (contact_radius_flag)
        contact_radius = calculate_contact_radius();
      Fnormal = calculate_forces_normal();

      history_index = model_history_index[DAMPING];
      Fdamp = calculate_forces_damping();
      Fntot = Fnormal + Fdamp;
      if (limit_damping && Fntot < 0.0) Fntot = 0.0;

      set_fncrit(); // Needed for tangential, rolling, twisting
      history_index = model_history_index[TANGENTIAL];
      calculate_forces_tangential();

      // sum normal + tangential contributions

      scale3(Fntot, nx, forces);
      add3(forces, fs, forces);

      // May need to eventually rethink tris..
      cross3(nx, fs, torquesi);
      scale3(-1, torquesi);

      copy3(torquesi, torquesj);

      dist_to_contact = radi - 0.5 * delta;
      scale3(dist_to_contact, torquesi);
      dist_to_contact = radj - 0.5 * delta;
      scale3(dist_to_contact, torquesj);


      // Extra modes

      if (rolling_defined || twisting_defined)
        sub3(omegai, omegaj, relrot);

      if (rolling_defined) {
        // rolling velocity, see eq. 31 of Wang et al, Particuology v 23, p 49 (2015)
        // this is different from the Marshall papers, which use the Bagi/Kuhn formulation
        // for rolling velocity (see Wang et al for why the latter is wrong)
        vrl[0] = Reff * (relrot[1] * nx[2] - relrot[2] * nx[1]);
        vrl[1] = Reff * (relrot[2] * nx[0] - relrot[0] * nx[2]);
        vrl[2] = Reff * (relrot[0] * nx[1] - relrot[1] * nx[0]);

        history_index = model_history_index[ROLLING];
        calculate_forces_rolling();

        double torroll[3];
        cross3(nx, fr, torroll);
        scale3(Reff, torroll);
        add3(torquesi, torroll, torquesi);
        if (contact_type == PAIR) sub3(torquesj, torroll, torquesj);
      }

      if (twisting_defined) {
        // omega_T (eq 29 of Marshall)
        magtwist = dot3(relrot, nx);

        history_index = model_history_index[TWISTING];
        calculate_forces_twisting();

        double tortwist[3];
        scale3(magtortwist, nx, tortwist);
        add3(torquesi, tortwist, torquesi);
        if (contact_type == PAIR) sub3(torquesj, tortwist, torquesj);
      }

      if (heat_defined) {
        history_index = model_history_index[HEAT];
        dq = calculate_heat();
      }

      // apply forces & torques
      scale3(factor_lj, forces);
      add3(f[i], forces, f[i]);

      scale3(factor_lj, torquesi);
      add3(torque[i], torquesi, torque[i]);

      if (force->newton_pair || j < nlocal) {
        sub3(f[j], forces, f[j]);
        scale3(factor_lj, torquesj);
        add3(torque[j], torquesj, torque[j]);
      }

      if (heat_defined) {
        heatflow[i] += dq;
        if (force->newton_pair || j < nlocal) heatflow[j] -= dq;
      }

      if (evflag) {
        ev_tally_xyz(i, j, nlocal, force->newton_pair,
          0.0, 0.0, forces[0], forces[1], forces[2], dx[0], dx[1], dx[2]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGranularAccelerated::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutoff_type,n+1,n+1,"pair:cutoff_type");
  memory->create(types_indices,n+1,n+1,"pair:types_indices");

  onerad_dynamic = new double[n+1];
  onerad_frozen = new double[n+1];
  maxrad_dynamic = new double[n+1];
  maxrad_frozen = new double[n+1];

  max_num_coeffs = 0;
  for (int i = 0; i < NSUBMODELS; i++)
    max_num_coeffs = MAX(max_num_coeffs, model_num_coeffs[i]);
  memory->create(model_coeffs, n+1, n+1, NSUBMODELS, max_num_coeffs);
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranularAccelerated::settings(int narg, char **arg)
{
  if (narg == 1) {
    cutoff_global = utils::numeric(FLERR,arg[0],false,lmp);
  } else {
    cutoff_global = -1; // will be set based on particle sizes, model choice
  }
}

/* ----------------------------------------------------------------------
  set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGranularAccelerated::coeff(int narg, char **arg)
{
  double cutoff_one = -1;

  if (narg < 3)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  model_defined[NSUBMODELS];
  for (int i = 0; i < NSUBMODELS; i++) model_defined[i] = 0;

  //Parse mandatory specification

  int iarg = 2;
  if (arg[iarg++] != model_names[NORMAL])
    error->all(FLERR, "Normal model does not match");

  if (iarg + model_num_coeffs[NORMAL] > narg)
    error->all(FLERR, "Insufficient arguments provided for {} model", model_names[NORMAL]);
  for (int i = 0; i < model_num_coeffs[NORMAL]; i++) {
    // A few parameters (e.g. kt for tangential mindlin) allow null
    if (strcmp(arg[iarg + i], "NULL") == 0)
      coeffs[i] = -1;
    else
      coeffs[i] = utils::numeric(FLERR, arg[iarg + i], false, lmp);
  }

  iarg += model_num_coeffs[NORMAL];
  model_defined[NORMAL] = 1;

  //Parse optional arguments

  int model_index;
  while (iarg < narg) {

    model_index = -1;
    if (strcmp(arg[iarg], "tangential") == 0) {
      model_index = TANGENTIAL;
    } else if (strcmp(arg[iarg], "damping") == 0) {
      model_index = DAMPING;
    } else if (strcmp(arg[iarg], "rolling") == 0) {
      model_index = ROLLING;
    } else if (strcmp(arg[iarg], "twisting") == 0) {
      model_index = TWISTING;
    } else if (strcmp(arg[iarg], "heat") == 0) {
      model_index = HEAT;
    } else if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal pair_coeff command, not enough parameters for cutoff keyword");
      cutoff_one = utils::numeric(FLERR,arg[iarg + 1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "limit_damping") == 0) {
      limit_damping = 1;
      iarg += 1;
    } else error->all(FLERR, "Illegal pair_coeff command {}", arg[iarg]);

    if (model_index != -1) {
      if (arg[iarg++] != model_names[model_index])
      error->all(FLERR, "Submodel model does not match");

      if (iarg + model_num_coeffs[model_index] > narg)
        error->all(FLERR, "Insufficient arguments provided for {} model", model_names[model_index]);
      for (int n = 0; n < model_num_coeffs[model_index]; i++) {
        // A few parameters (e.g. kt for tangential mindlin) allow null
        for (int i = ilo; i <= ihi; i++) {
          for (int j = MAX(jlo,i); j <= jhi; j++) {
            if (strcmp(arg[iarg + i], "NULL") == 0)
              model_coeffs[i][j][model_index][n] = -1;
            else
              model_coeffs[i][j][model_index][n] = utils::numeric(FLERR, arg[iarg + i], false, lmp);
          }
        }
      }

      iarg += model_num_coeffs[model_index];
      model_defined[model_index] = 1;
    }
  }

  for (int i = 0; i < NSUBMODELS; i++) {
    if (!model_defined[i] && (model_name[i] != "none"))
      error->all(FLERR, "Submodel missing from pair_coeff command");
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cutoff_type[i][j] = cutoff_type[j][i] = cutoff_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
  init specific to this pair style
------------------------------------------------------------------------- */

void PairGranularAccelerated::init_style()
{
  // error and warning checks

  if (!atom->radius_flag || !atom->rmass_flag || !atom->omega_flag)
    error->all(FLERR,"Pair granular requires atom attributes radius, rmass, omega");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair granular requires ghost atoms store velocity");

  if (heat_defined) {
    if (!atom->temperature_flag)
      error->all(FLERR,"Heat conduction in pair granular requires atom style with temperature property");
    if (!atom->heatflow_flag)
      error->all(FLERR,"Heat conduction in pair granular requires atom style with heatflow property");
  }

  // Calculate history variables

  if (size_history == 0)
    use_history = 0;
  else
    use_history = 1;

  model_history_index[0] = 0;
  for (int i = 1; i < NSUBMODELS; i++)
    model_history_index[i] = mdoel_history_index[i - 1] + model_size_history[i - 1];

  if (use_history) neighbor->add_request(this, NeighConst::REQ_SIZE | NeighConst::REQ_HISTORY);
  else neighbor->add_request(this, NeighConst::REQ_SIZE);

  // if history is stored and first init, create Fix to store history
  // it replaces FixDummy, created in the constructor
  // this is so its order in the fix list is preserved

  if (use_history && fix_history == nullptr) {
    fix_history = dynamic_cast<FixNeighHistory *>(modify->replace_fix("NEIGH_HISTORY_GRANULAR_DUMMY",
                                                          "NEIGH_HISTORY_GRANULAR"
                                                          " all NEIGH_HISTORY "
                                                          + std::to_string(size_history),1));
    fix_history->pair = this;
  } else if (use_history) {
    fix_history = dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id("NEIGH_HISTORY_GRANULAR"));
    if (!fix_history) error->all(FLERR,"Could not find pair fix neigh history ID");
  }

  // check for FixFreeze and set freeze_group_bit

  auto fixlist = modify->get_fix_by_style("^freeze");
  if (fixlist.size() == 0)
    freeze_group_bit = 0;
  else if (fixlist.size() > 1)
    error->all(FLERR, "Only one fix freeze command at a time allowed");
  else
    freeze_group_bit = fixlist.front()->groupbit;

  // check for FixRigid so can extract rigid body masses

  fix_rigid = nullptr;
  for (const auto &ifix : modify->get_fix_list()) {
    if (ifix->rigid_flag) {
      if (fix_rigid)
        error->all(FLERR, "Only one fix rigid command at a time allowed");
      else fix_rigid = ifix;
    }
  }

  // check for FixPour and FixDeposit so can extract particle radii

  auto pours = modify->get_fix_by_style("^pour");
  auto deps = modify->get_fix_by_style("^deposit");

  // set maxrad_dynamic and maxrad_frozen for each type
  // include future FixPour and FixDeposit particles as dynamic

  int itype;
  for (int i = 1; i <= atom->ntypes; i++) {
    onerad_dynamic[i] = onerad_frozen[i] = 0.0;
    for (auto &ipour : pours) {
      itype = i;
      double maxrad = *((double *) ipour->extract("radius", itype));
      if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
    }
    for (auto &idep : deps) {
      itype = i;
      double maxrad = *((double *) idep->extract("radius", itype));
      if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
    }
  }

  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & freeze_group_bit)
      onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]], radius[i]);
    else
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);
  }

  MPI_Allreduce(&onerad_dynamic[1],&maxrad_dynamic[1],atom->ntypes,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&onerad_frozen[1],&maxrad_frozen[1],atom->ntypes,MPI_DOUBLE,MPI_MAX,world);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGranularAccelerated::init_one(int i, int j)
{
  double cutoff = 0.0;

  if (setflag[i][j] == 0)
    error->all(FLERR, "Pair style granular/accelerated does not support mixing coefficients");

  // It is possible that cut[i][j] at this point is still 0.0.
  // This can happen when
  // there is a future fix_pour after the current run. A cut[i][j] = 0.0 creates
  // problems because neighbor.cpp uses min(cut[i][j]) to decide on the bin size
  // To avoid this issue, for cases involving  cut[i][j] = 0.0 (possible only
  // if there is no current information about radius/cutoff of type i and j).
  // we assign cutoff = max(cut[i][j]) for i,j such that cut[i][j] > 0.0.

  double pulloff;
  if (cutoff_type[i][j] < 0 && cutoff_global < 0) {
    if (((maxrad_dynamic[i] > 0.0) && (maxrad_dynamic[j] > 0.0)) ||
        ((maxrad_dynamic[i] > 0.0) &&  (maxrad_frozen[j] > 0.0)) ||
        // radius info about both i and j exist
        ((maxrad_frozen[i] > 0.0)  && (maxrad_dynamic[j] > 0.0))) {
      cutoff = maxrad_dynamic[i] + maxrad_dynamic[j];
      pulloff = 0.0;
      if (beyond_contact) {
        pulloff = pulloff_distance(maxrad_dynamic[i], maxrad_dynamic[j]);
        cutoff += pulloff;

        pulloff = pulloff_distance(maxrad_frozen[i], maxrad_dynamic[j]);
        cutoff = MAX(cutoff, maxrad_frozen[i] + maxrad_dynamic[j] + pulloff);

        pulloff = pulloff_distance(maxrad_dynamic[i], maxrad_frozen[j]);
        cutoff = MAX(cutoff,maxrad_dynamic[i] + maxrad_frozen[j] + pulloff);
      }
    } else {
      // radius info about either i or j does not exist
      // (i.e. not present and not about to get poured;
      // set to largest value to not interfere with neighbor list)

      double cutmax = 0.0;
      for (int k = 1; k <= atom->ntypes; k++) {
        cutmax = MAX(cutmax,2.0*maxrad_dynamic[k]);
        cutmax = MAX(cutmax,2.0*maxrad_frozen[k]);
      }
      cutoff = cutmax;
    }
  } else if (cutoff_type[i][j] > 0) {
    cutoff = cutoff_type[i][j];
  } else if (cutoff_global > 0) {
    cutoff = cutoff_global;
  }

  dt = update->dt;
  for (int n = 0; n < NSUBMODELS; n++)
    for (int m = 0; m < max_num_coeffs; m++)
      model_coeffs[j][i][n][m] = model_coeffs[i][j][n][m];

  return cutoff;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGranularAccelerated::write_restart(FILE *fp)
{
  int i, j, n, m;

  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&cutoff_type[i][j], sizeof(double), 1, fp);
        for (n = 0; n < NSUBMODELS; n++)
          for (m = 0; m < max_num_coeffs; m++)
            fwrite(&model_coeffs[i][j][n][m], sizeof(double), 1, fp);
      }
    }
  }

  fwrite(&limit_damping, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGranularAccelerated::read_restart(FILE *fp)
{
  allocate();
  int i, j, n, m;
  int me = comm->me;

  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0)
          utils::sfread(FLERR, &cutoff_type[i][j], sizeof(double), 1, fp, nullptr, error);
        MPI_Bcast(&cutoff_type[i][j], 1, MPI_DOUBLE, 0, world);
        for (n = 0; n < NSUBMODELS; n++) {
          for (m = 0; m < max_num_coeffs; m++) {
            if (me == 0)
              utils::sfread(FLERR,&model_coeffs[i][j][n][m], sizeof(double), 1, fp, nullptr, error);
            MPI_Bcast(&model_coeffs[i][j][n][m], 1, MPI_DOUBLE, 0, world);
          }
        }
      }
    }
  }

  if (comm->me == 0)
    utils::sfread(FLERR, &limit_damping, sizeof(int), 1, fp, nullptr, error);
  MPI_Bcast(&limit_damping, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

double PairGranularAccelerated::single(int i, int j, int itype, int jtype,
                            double /*rsq*/, double /* factor_coul */,
                            double factor_lj, double &fforce)
{
  if (factor_lj == 0) {
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }

  int nall = atom->nlocal + atom->nghost;
  if ((i >= nall) || (j >= nall))
    error->all(FLERR,"Not enough atoms for pair granular single function");

  class GranularModel* model = models_list[types_indices[itype][jtype]];

  // Reset model and copy initial geometric data
  double **x = atom->x;
  double *radius = atom->radius;

  xi = x[i];
  xj = x[j];
  radi = radius[i];
  radj = radius[j];
  history_update = 0; // Don't update history

  // If history is needed
  double *history,*allhistory;
  int jnum = list->numneigh[i];
  int *jlist = list->firstneigh[i];
  if (use_history) {
    if ((fix_history == nullptr) || (fix_history->firstvalue == nullptr))
      error->one(FLERR,"Pair granular single computation needs history");
    allhistory = fix_history->firstvalue[i];
    for (int jj = 0; jj < jnum; jj++) {
      neighprev++;
      if (neighprev >= jnum) neighprev = 0;
      if (jlist[neighprev] == j) break;
    }
    history = &allhistory[size_history * neighprev];
    touch = fix_history->firstflag[i][neighprev];
  }

  int touchflag = check_contact();

  if (!touchflag) {
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }

  // meff = effective mass of pair of particles
  // if I or J part of rigid body, use body mass
  // if I or J is frozen, meff is other particle
  double *rmass = atom->rmass;
  int *mask = atom->mask;

  double mi = rmass[i];
  double mj = rmass[j];
  if (fix_rigid) {
    if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
    if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
  }
  double meff = mi * mj / (mi + mj);
  if (mask[i] & freeze_group_bit) meff = mj;
  if (mask[j] & freeze_group_bit) meff = mi;

  // Copy additional information and calculate forces
  double **v = atom->v;
  double **omega = atom->omega;

  meff = meff;
  vi = v[i];
  vj = v[j];
  omegai = omega[i];
  omegaj = omega[j];
  history = history;

  calculate_forces();

  // apply forces & torques
  // Calculate normal component, normalized by r
  fforce = Fntot * rinv;

  // set single_extra quantities
  svector[0] = fs[0];
  svector[1] = fs[1];
  svector[2] = fs[2];
  svector[3] = MathExtra::len3(fs);
  svector[4] = fr[0];
  svector[5] = fr[1];
  svector[6] = fr[2];
  svector[7] = MathExtra::len3(fr);
  svector[8] = magtortwist;
  svector[9] = dx[0];
  svector[10] = dx[1];
  svector[11] = dx[2];

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairGranularAccelerated::pack_forward_comm(int n, int *list, double *buf,
                                    int /* pbc_flag */, int * /* pbc */)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = mass_rigid[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    mass_rigid[i] = buf[m++];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairGranularAccelerated::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   transfer history during fix/neigh/history exchange
   only needed if any history entries i-j are not just negative of j-i entries
------------------------------------------------------------------------- */

void PairGranularAccelerated::transfer_history(double* source, double* target, int itype, int jtype)
{
  if (nondefault_history_transfer) {
    for (int i = 0; i < size_history; i++) {
      target[i] = transfer_history_factor[i] * source[i];
    }
  } else {
    for (int i = 0; i < size_history; i++) {
      target[i] = -source[i];
    }
  }
}

/* ----------------------------------------------------------------------
   self-interaction range of particle
------------------------------------------------------------------------- */

double PairGranularAccelerated::atom2cut(int i)
{
  double cut;

  cut = atom->radius[i] * 2;
  if (beyond_contact) {
    cut += pulloff_distance(cut, cut);
  }

  return cut;
}

/* ----------------------------------------------------------------------
   maximum interaction range for two finite particles
------------------------------------------------------------------------- */

double PairGranularAccelerated::radii2cut(double r1, double r2)
{
  double cut = 0.0;

  if (beyond_contact) {
    int n = atom->ntypes;
    double temp = pulloff_distance(r1, r2);
    if (temp > cut) cut = temp;
  }

  cut += r1 + r2;

  return cut;
}

/* ----------------------------------------------------------------------
   map coeffs and initialize
------------------------------------------------------------------------- */

void PairGranularAccelerated::prep_model(int itype, int jtype)
{
  coeffs = &model_coeffs[itype][jtype][NORMAL][0];
  coeffs_to_local_normal();
  coeffs = &model_coeffs[itype][jtype][DAMPING][0];
  coeffs_to_local_damping();
  coeffs = &model_coeffs[itype][jtype][TANGENTIAL][0];
  coeffs_to_local_tangential();
  if (rolling_defined) {
    coeffs = &model_coeffs[itype][jtype][ROLLING][0];
    coeffs_to_local_rolling();
  }
  if (twisting_defined) {
    coeffs = &model_coeffs[itype][jtype][TWISTING][0];
    coeffs_to_local_twisting();
  }
  if (heat_defined) {
    coeffs = &model_coeffs[itype][jtype][HEAT][0];
    coeffs_to_local_heat();
  }

  init_normal();
  init_damping();
  init_tangential();
  if (rolling_defined) init_rolling();
  if (twisting_defined) init_twisting();
  if (heat_defined) init_heat();
}

/* ----------------------------------------------------------------------
   mixing of Young's modulus (E)
------------------------------------------------------------------------- */

double PairGranularAccelerated::mix_stiffnessE(double E1, double E2, double poiss1, double poiss2)
{
  double factor1 = (1 - poiss1 * poiss1) / E1;
  double factor2 = (1 - poiss2 * poiss2) / E2;
  return 1 / (factor1 + factor2);
}

/* ----------------------------------------------------------------------
   mixing of shear modulus (G)
------------------------------------------------------------------------ */

double PairGranularAccelerated::mix_stiffnessG(double E1, double E2, double poiss1, double poiss2)
{
  double factor1 = 2 * (2 - poiss1) * (1 + poiss1) / E1;
  double factor2 = 2 * (2 - poiss2) * (1 + poiss2) / E2;
  return 1 / (factor1 + factor2);
}

/* ----------------------------------------------------------------------
      CUSTOMIZE: Methods from Normal SubModel
        Delete all instances of ""
------------------------------------------------------------------------- */

bool PairGranularAccelerated::touch()
{
  bool touchflag = (rsq < radsum * radsum);
  return touchflag;
}

/* ---------------------------------------------------------------------- */

double PairGranularAccelerated::pulloff_distance(double /*radi*/, double /*radj*/)
{
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double PairGranularAccelerated::calculate_contact_radius()
{
  return sqrt(dR);
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::set_fncrit()
{
  Fncrit = fabs(Fntot);
}

/* --------------------------------------------------------------------------
      CUSTOMIZE: Methods from each SubModel, necessary changes include:
        Delete all instances of ""
        Contact_type is always PAIR
        Replace unnecessary accessors directly w/ variable
          normal_get_emod() -> Emod
          normal_get_poiss() -> poiss
          normal_get_fncrit() -> Fncrit
          normal_get_damp() -> damp_norm
          damping_get_damp_prefactor() -> damp_prefactor
        Delete normal_get_material_properties() checks, ensure manually
----------------------------------------------------------------------------- */

void PairGranularAccelerated::coeffs_to_local_normal()
{
  Emod = coeffs[0];
  damp_norm = coeffs[1];
  poiss = coeffs[2];
  if (!mixed_coefficients)
    k_norm = FOURTHIRDS * mix_stiffnessE(Emod, Emod, poiss, poiss);

  if (Emod < 0.0 || damp_norm < 0.0) error->all(FLERR, "Illegal Hertz material normal model");
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::init_normal()
{
}

/* ---------------------------------------------------------------------- */

double PairGranularAccelerated::calculate_forces_normal()
{
  return k_norm * contact_radius * delta;
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::coeffs_to_local_damping()
{
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::init_damping()
{
  double tmp = damp_norm;
  damp_damp = 1.2728 - 4.2783 * tmp + 11.087 * square(tmp);
  damp_damp += -22.348 * cube(tmp) + 27.467 * powint(tmp, 4);
  damp_damp += -18.022 * powint(tmp, 5) + 4.8218 * powint(tmp, 6);
}

/* ---------------------------------------------------------------------- */

double PairGranularAccelerated::calculate_forces_damping()
{
  // in case argument <= 0 due to precision issues
  double sqrt1;
  if (delta > 0.0)
    sqrt1 = MAX(0.0, meff * Fnormal / delta);
  else
    sqrt1 = 0.0;
  damp_prefactor = damp_damp * sqrt(sqrt1);
  return -damp_prefactor * vnnr;
}


/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::coeffs_to_local_tangential()
{
  k_tang = coeffs[0];
  xt = coeffs[1];
  mu_tang = coeffs[2];

  if (k_tang == -1)
    k_tang = 8.0 * mix_stiffnessG(Emod, Emod, poiss, poiss);

  if (k_tang < 0.0 || xt < 0.0 || mu_tang < 0.0) error->all(FLERR, "Illegal Mindlin tangential model");
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::init_tangential()
{
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::calculate_forces_tangential()
{
  double k_scaled, magfs, magfs_inv, rsht, shrmag, prjmag, temp_dbl;
  double temp_array[3];
  int frame_update = 0;

  damp_tang = xt * damp_prefactor;

  double *history = &history[history_index];
  double Fscrit = Fncrit * mu_tang;

  k_scaled = k_tang * contact_radius;

  // on unloading, rescale the shear displacements/force
  if (mindlin_rescale)
    if (contact_radius < history[3]) scale3(contact_radius / history[3], history);

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (history_update) {
    rsht = dot3(history, nx);
    if (mindlin_force) {
      frame_update = fabs(rsht) > (EPSILON * Fscrit);
    } else {
      frame_update = (fabs(rsht) * k_scaled) > (EPSILON * Fscrit);
    }

    if (frame_update) {
      shrmag = len3(history);
      // projection
      scale3(rsht, nx, temp_array);
      sub3(history, temp_array, history);
      // also rescale to preserve magnitude
      prjmag = len3(history);
      if (prjmag > 0)
        temp_dbl = shrmag / prjmag;
      else
        temp_dbl = 0;
      scale3(temp_dbl, history);
    }

    // update history
    if (mindlin_force) {
      // tangential force
      // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
      scale3(-k_scaled * dt, vtr, temp_array);
    } else {
      scale3(dt, vtr, temp_array);
    }
    add3(history, temp_array, history);

    if (mindlin_rescale) history[3] = contact_radius;
  }

  // tangential forces = history + tangential velocity damping
  scale3(-damp_tang, vtr, fs);

  if (!mindlin_force) {
    scale3(k_scaled, history, temp_array);
    sub3(fs, temp_array, fs);
  } else {
    add3(fs, history, fs);
  }

  // rescale frictional displacements and forces if needed
  magfs = len3(fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, fs, history);
      scale3(damp_tang, vtr, temp_array);
      add3(history, temp_array, history);

      if (!mindlin_force) scale3(-1.0 / k_scaled, history);

      scale3(Fscrit * magfs_inv, fs);
    } else {
      zero3(fs);
    }
  }
}


/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::coeffs_to_local_rolling()
{
  k_roll= coeffs[0];
  gamma = coeffs[1];
  mu_roll = coeffs[2];

  if (k_roll< 0.0 || mu_roll < 0.0 || gamma < 0.0) error->all(FLERR, "Illegal SDS rolling model");
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::init_rolling()
{
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::calculate_forces_rolling()
{
  int rhist0, rhist1, rhist2, frameupdate;
  double Frcrit, rolldotn, rollmag, prjmag, magfr, hist_temp[3], scalefac, temp_array[3];
  double k_inv, magfr_inv;

  rhist0 = history_index;
  rhist1 = rhist0 + 1;
  rhist2 = rhist1 + 1;

  Frcrit = mu_roll * Fncrit;

  if (history_update) {
    hist_temp[0] = history[rhist0];
    hist_temp[1] = history[rhist1];
    hist_temp[2] = history[rhist2];
    rolldotn = dot3(hist_temp, nx);

    frameupdate = (fabs(rolldotn) * k_roll) > (EPSILON * Frcrit);
    if (frameupdate) {    // rotate into tangential plane
      rollmag = len3(hist_temp);
      // projection
      scale3(rolldotn, nx, temp_array);
      sub3(hist_temp, temp_array, hist_temp);

      // also rescale to preserve magnitude
      prjmag = len3(hist_temp);
      if (prjmag > 0)
        scalefac = rollmag / prjmag;
      else
        scalefac = 0;
      scale3(scalefac, hist_temp);
    }
    scale3(dt, vrl, temp_array);
    add3(hist_temp, temp_array, hist_temp);
  }

  scaleadd3(-k_roll, hist_temp, -gamma, vrl, fr);

  // rescale frictional displacements and forces if needed
  magfr = len3(fr);
  if (magfr > Frcrit) {
    rollmag = len3(hist_temp);
    if (rollmag != 0.0) {
      k_inv = 1.0 / k_roll;
      magfr_inv = 1.0 / magfr;
      scale3(-Frcrit * k_inv * magfr_inv, fr, hist_temp);
      scale3(-gamma * k_inv, vrl, temp_array);
      add3(hist_temp, temp_array, hist_temp);

      scale3(Frcrit * magfr_inv, fr);
    } else {
      zero3(fr);
    }
  }

  if (history_update) {
    history[rhist0] = hist_temp[0];
    history[rhist1] = hist_temp[1];
    history[rhist2] = hist_temp[2];
  }
}


/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::coeffs_to_local_twisting()
{
  k_twist = coeffs[0];
  damp_twist = coeffs[1];
  mu_twist = coeffs[2];

  if (k_twist < 0.0 || mu_twist < 0.0 || damp_twist < 0.0) error->all(FLERR, "Illegal SDS twisting model");
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::init_twisting()
{
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::calculate_forces_twisting()
{
  double signtwist, Mtcrit;

  if (history_update) { history[history_index] += magtwist * dt; }

  // M_t torque (eq 30)
  magtortwist = -k_twist * history[history_index] - damp_twist * magtwist;
  signtwist = (magtwist > 0) - (magtwist < 0);
  Mtcrit = mu_twist * Fncrit;    // critical torque (eq 44)

  if (fabs(magtortwist) > Mtcrit) {
    history[history_index] = (Mtcrit * signtwist - damp_twist * magtwist) / k_twist;
    magtortwist = -Mtcrit * signtwist;    // eq 34
  }
}


/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::coeffs_to_local_heat()
{
  conductivity = coeffs[0];

  if (conductivity < 0.0) error->all(FLERR, "Illegal radius heat model");
}

/* ---------------------------------------------------------------------- */

void PairGranularAccelerated::init_heat()
{
}

/* ---------------------------------------------------------------------- */

double PairGranularAccelerated::calculate_heat()
{
  return 2 * conductivity * contact_radius * (Tj - Ti);
}
