/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "gran_sub_mod_twisting.h"

#include "error.h"
#include "gran_sub_mod_normal.h"
#include "gran_sub_mod_tangential.h"
#include "granular_model.h"
#include "math_const.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace Granular_NS;

using MathConst::TWOTHIRDS;

/* ----------------------------------------------------------------------
   Default twisting model
------------------------------------------------------------------------- */

GranSubModTwisting::GranSubModTwisting(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp) {}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModTwistingNone::GranSubModTwistingNone(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTwisting(gm, lmp)
{
}

/* ----------------------------------------------------------------------
   Marshall twisting model
------------------------------------------------------------------------- */

GranSubModTwistingMarshall::GranSubModTwistingMarshall(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTwisting(gm, lmp)
{
  num_coeffs = 0;
  size_history = 3;
  contact_radius_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModTwistingMarshall::init()
{
  k_tang = gm->tangential_model->get_k();
  mu_tang = gm->tangential_model->get_mu();
}

/* ---------------------------------------------------------------------- */

void GranSubModTwistingMarshall::calculate_forces()
{
  double signtwist, Mtcrit;

  // Calculate twist coefficients from tangential model & contact geometry
  // eq 32 of Marshall paper
  double k_twist = 0.5 * k_tang * gm->contact_radius * gm->contact_radius;
  double damp_twist = 0.5 * gm->tangential_model->get_damp() * gm->contact_radius * gm->contact_radius;
  double mu_twist = TWOTHIRDS * mu_tang * gm->contact_radius;

  if (gm->history_update) { gm->history[history_index] += gm->magtwist * gm->dt; }

  // M_t torque (eq 30)
  gm->magtortwist = -k_twist * gm->history[history_index] - damp_twist * gm->magtwist;
  signtwist = (gm->magtwist > 0) - (gm->magtwist < 0);
  Mtcrit = mu_twist * gm->normal_model->get_fncrit();    // critical torque (eq 44)

  if (fabs(gm->magtortwist) > Mtcrit) {
    gm->history[history_index] = (Mtcrit * signtwist - damp_twist * gm->magtwist) / k_twist;
    gm->magtortwist = -Mtcrit * signtwist;    // eq 34
  }
}

/* ----------------------------------------------------------------------
   SDS twisting model
------------------------------------------------------------------------- */

GranSubModTwistingSDS::GranSubModTwistingSDS(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTwisting(gm, lmp)
{
  num_coeffs = 3;
  size_history = 3;
}

/* ---------------------------------------------------------------------- */

void GranSubModTwistingSDS::coeffs_to_local()
{
  k_twist = coeffs[0];
  damp_twist = coeffs[1];
  mu_twist = coeffs[2];

  if (k_twist < 0.0 || mu_twist < 0.0 || damp_twist < 0.0) error->all(FLERR, "Illegal SDS twisting model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTwistingSDS::calculate_forces()
{
  double signtwist, Mtcrit;

  if (gm->history_update) { gm->history[history_index] += gm->magtwist * gm->dt; }

  // M_t torque (eq 30)
  gm->magtortwist = -k_twist * gm->history[history_index] - damp_twist * gm->magtwist;
  signtwist = (gm->magtwist > 0) - (gm->magtwist < 0);
  Mtcrit = mu_twist * gm->normal_model->get_fncrit();    // critical torque (eq 44)

  if (fabs(gm->magtortwist) > Mtcrit) {
    gm->history[history_index] = (Mtcrit * signtwist - damp_twist * gm->magtwist) / k_twist;
    gm->magtortwist = -Mtcrit * signtwist;    // eq 34
  }
}
