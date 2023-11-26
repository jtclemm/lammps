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
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "fix_def_uef_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "irregular.h"
#include "kspace.h"
#include "math_const.h"
#include "modify.h"
#include "uef_utils.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixDefUEFKokkos::FixDefUEFKokkos(LAMMPS *lmp, int narg, char **arg) : FixDeform(lmp, narg, arg)
{
  kokkosable = 1;
  domainKK = (DomainKokkos *) domain;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ----------------------------------------------------------------------
  box flipped on previous step
  reset box tilts for flipped config and create new box in domain
  image_flip() adjusts image flags due to box shape change induced by flip
  remap() puts atoms outside the new box back into the new box
  perform irregular on atoms in lamda coords to migrate atoms to new procs
  important that image_flip comes before remap, since remap may change
    image flags to new values, making eqs in doc of Domain:image_flip incorrect
------------------------------------------------------------------------- */

void FixDefUEFKokkos::pre_exchange()
{
  if (flip == 0) return;

  // go to lab frame
  FixNVEKokkosInitialIntegrateFunctor<DeviceType,1> functor(this);
    Kokkos::parallel_for(nlocal,functor);

  boxlo[0] = domain->boxlo[0];
  boxlo[1] = domain->boxlo[1];
  boxlo[2] = domain->boxlo[2];

  inv_rotate_x();
  inv_rotate_v();
  inv_rotate_f();

  // get & set the new box and rotation matrix
  double vol = domainK->xprd * domain->yprd * domain->zprd;
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
  rotate_v();
  rotate_x();
  rotate_f();

  // put all atoms in the new box
  domainKK->remap_all();

  domainKK->x2lamda(atom->nlocal);
  atomKK->sync(Host,ALL_MASK);
  irregular->migrate_atoms();
  atomKK->modified(Host,ALL_MASK);
  domainKK->lamda2x(atom->nlocal);

  flip = 0;
}

/* ---------------------------------------------------------------------- */

void FixDefUEFKokkos::end_of_step()
{
  double iv = domain->xprd*domain->yprd*domain->zprd;
  double dtv = update->dt;
  double ex = rate[0] * dtv;
  double ey = rate[1] * dtv;
  strain[0] += ex;
  strain[1] += ey;
  strain[2] += -ex - ey;

  int i;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // convert atoms and rigid bodies to lamda coords

  domainKK->x2lamda(nlocal);

  if (rfix.size() > 0) {
    atomKK->sync(Host,ALL_MASK);
    for (auto &ifix : rfix)
      ifix->deform(0);
    atomKK->modified(Host,ALL_MASK);
  }

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
  domain->set_global_box();
  domain->set_local_box();
  uefbox->get_rot(rot);

  // convert atoms and rigid bodies back to box coords

  domainKK->lamda2x(nlocal);

  if (rfix.size() > 0) {
    atomKK->sync(Host,ALL_MASK);
    for (auto &ifix : rfix)
      ifix->deform(1);
    atomKK->modified(Host,ALL_MASK);
  }

  // redo KSpace coeffs since box has changed

  if (kspace_flag) force->kspace->setup();
  if(uefbox->reduce()) flip = 1;
  if(flip) next_reneighbor = update->ntimestep + 1;
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
void FixDefUEFKokkos::rotate_x()
{
  atomKK->sync(Device,X_MASK);

  x = atomKK->k_x.view<LMPDeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType, TagFixDefUEF_rotate_x>(0,n),*this);
  copymode = 0;

  atomKK->modified(Device,X_MASK);
}

KOKKOS_INLINE_FUNCTION
void FixDefUEFKokkos::operator()(TagFixDefUEF_rotate_x, const int &i) const {
  const double xi0 = x(i,0);
  const double xi1 = x(i,1);
  const double xi2 = x(i,2);

  x(i,0) = boxlo[0] + rot[0][0] * xi0 + rot[0][1] * xi1 + rot[0][2] * xi2;
  x(i,1) = boxlo[1] + rot[1][0] * xi0 + rot[1][1] * xi1 + rot[1][2] * xi2;
  x(i,2) = boxlo[2] + rot[2][0] * xi0 + rot[2][1] * xi1 + rot[2][2] * xi2;
}

void FixDefUEFKokkos::inv_rotate_x()
{
  atomKK->sync(Device,X_MASK);

  x = atomKK->k_x.view<LMPDeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType, TagFixDefUEF_inv_rotate_x>(0,n),*this);
  copymode = 0;

  atomKK->modified(Device,X_MASK);
}

KOKKOS_INLINE_FUNCTION
void FixDefUEFKokkos::operator()(TagFixDefUEF_inv_rotate_x, const int &i) const {
  const double xi0 = x(i,0) - boxlo[0];
  const double xi1 = x(i,1) - boxlo[1];
  const double xi2 = x(i,2) - boxlo[2];

  x(i,0) = rot[0][0] * xi0 + rot[1][0] * xi1 + rot[2][0] * xi2;
  x(i,1) = rot[0][1] * xi0 + rot[1][1] * xi1 + rot[2][1] * xi2;
  x(i,2) = rot[0][2] * xi0 + rot[1][2] * xi1 + rot[2][2] * xi2;
}

void FixDefUEFKokkos::rotate_v()
{
  atomKK->sync(Device,V_MASK);

  v = atomKK->k_v.view<LMPDeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType, TagFixDefUEF_rotate_v>(0,n),*this);
  copymode = 0;

  atomKK->modified(Device,V_MASK);
}

KOKKOS_INLINE_FUNCTION
void FixDefUEFKokkos::operator()(TagFixDefUEF_rotate_v, const int &i) const {
  const double vi0 = v(i,0);
  const double vi1 = v(i,1);
  const double vi2 = v(i,2);

  v(i,0) = rot[0][0] * vi0 + rot[0][1] * vi1 + rot[0][2] * vi2;
  v(i,1) = rot[1][0] * vi0 + rot[1][1] * vi1 + rot[1][2] * vi2;
  v(i,2) = rot[2][0] * vi0 + rot[2][1] * vi1 + rot[2][2] * vi2;
}

void FixDefUEFKokkos::inv_rotate_v()
{
  atomKK->sync(Device,V_MASK);

  v = atomKK->k_v.view<LMPDeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType, TagFixDefUEF_inv_rotate_v>(0,n),*this);
  copymode = 0;

  atomKK->modified(Device,V_MASK);
}

KOKKOS_INLINE_FUNCTION
void FixDefUEFKokkos::operator()(TagFixDefUEF_inv_rotate_v, const int &i) const {
  const double vi0 = v(i,0);
  const double vi1 = v(i,1);
  const double vi2 = v(i,2);

  v(i,0) = rot[0][0] * vi0 + rot[1][0] * vi1 + rot[2][0] * vi2;
  v(i,1) = rot[0][1] * vi0 + rot[1][1] * vi1 + rot[2][1] * vi2;
  v(i,2) = rot[0][2] * vi0 + rot[1][2] * vi1 + rot[2][2] * vi2;
}

void FixDefUEFKokkos::rotate_f()
{
  atomKK->sync(Device,F_MASK);

  f = atomKK->k_f.view<LMPDeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType, TagFixDefUEF_rotate_f>(0,n),*this);
  copymode = 0;

  atomKK->modified(Device,F_MASK);
}

KOKKOS_INLINE_FUNCTION
void FixDefUEFKokkos::operator()(TagFixDefUEF_rotate_f, const int &i) const {
  const double fi0 = f(i,0);
  const double fi1 = f(i,1);
  const double fi2 = f(i,2);

  f(i,0) = rot[0][0] * fi0 + rot[0][1] * fi1 + rot[0][2] * fi2;
  f(i,1) = rot[1][0] * fi0 + rot[1][1] * fi1 + rot[1][2] * fi2;
  f(i,2) = rot[2][0] * fi0 + rot[2][1] * fi1 + rot[2][2] * fi2;
}

void FixDefUEFKokkos::inv_rotate_f()
{
  atomKK->sync(Device,F_MASK);

  f = atomKK->k_f.view<LMPDeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType, TagFixDefUEF_inv_rotate_f>(0,n),*this);
  copymode = 0;

  atomKK->modified(Device,F_MASK);
}

KOKKOS_INLINE_FUNCTION
void FixDefUEFKokkos::operator()(TagFixDefUEF_inv_rotate_f, const int &i) const {
  const double fi0 = f(i,0);
  const double fi1 = f(i,1);
  const double fi2 = f(i,2);

  f(i,0) = rot[0][0] * fi0 + rot[1][0] * fi1 + rot[2][0] * fi2;
  f(i,1) = rot[0][1] * fi0 + rot[1][1] * fi1 + rot[2][1] * fi2;
  f(i,2) = rot[0][2] * fi0 + rot[1][2] * fi1 + rot[2][2] * fi2;
}
