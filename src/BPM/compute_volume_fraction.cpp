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

#include "compute_volume_fraction.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "random_park.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeVolumeFraction::ComputeVolumeFraction(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg),
    random(nullptr)
{
  if (narg != 5) error->all(FLERR, "Illegal compute volume/fraction command");

  ntest = utils::inumeric(FLERR, arg[3], false, lmp);
  seed = utils::inumeric(FLERR, arg[4], false, lmp);

  if (seed <= 0 || ntest <= 0) error->all(FLERR, "Illegal compute volume/fraction command");

  scalar_flag = 1;

  // random number generator, same for all procs
  // warm up the generator 30x to avoid correlations in first-particle
  // positions if runs are repeated with consecutive seeds

  random = new RanPark(lmp, seed);
  for (int ii = 0; ii < 30; ii++) random->uniform();
}

/* ---------------------------------------------------------------------- */

ComputeVolumeFraction::~ComputeVolumeFraction()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

void ComputeVolumeFraction::init()
{
  auto request = neighbor->add_request(this, NeighConst::REQ_OCCASIONAL | NeighConst::SKIP);

  // Skip all types, don't need npair class
  int ntypes = atom->ntypes;
  int *iskip = new int[ntypes + 1];
  int **ijskip;
  memory->create(ijskip, ntypes + 1, ntypes + 1, "compute_volume_fraction:ijskip");

  for (int itype = 1; itype <= ntypes; itype++) {
    iskip[itype] = 1;
    for (int jtype = 1; jtype <= ntypes; jtype++) {
      ijskip[itype][jtype] = 1;
    }
  }

  request->set_skip(iskip, ijskip);
  memory->destroy(ijskip);
}

/* ---------------------------------------------------------------------- */

void ComputeVolumeFraction::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeVolumeFraction::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  // Trigger a neighbor list build to update stencil
  // NPair skip won't build a neighbor list
  neighbor->build_one(list)

  // could grab list->np->find_neighbors()
  //  virtual int find_neighbors(double *, double);
  // how do you define a one atom cut size?
  // work with multi vs standard? finite radius vs by-type cut?
  // purely error with JKR

  double xlo, ylo, zlo, xhi, yhi, zhi, zmid;
  double delx, dely, delz, distsq, odistsq;
  double lamda[3], xone[3], *coord;
  double *boxlo, *boxhi;
  int triclinic = domain->triclinic;

  if (triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
    zmid = zlo + 0.5 * (zhi - zlo);
  } else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
    zmid = zlo + 0.5 * (zhi - zlo);
    boxlo = domain->boxlo_lamda;
    boxhi = domain->boxhi_lamda;
  }

  int collisions = 0;
  for (int i = 0; i < ntest; i++) {
    // Using same recipe from create_atoms.cpp
    xone[0] = xlo + random->uniform() * (xhi - xlo);
    xone[1] = ylo + random->uniform() * (yhi - ylo);
    xone[2] = zlo + random->uniform() * (zhi - zlo);
    if (domain->dimension == 2) xone[2] = zmid;

    if (triclinic) {
      domain->x2lamda(xone, lamda);
      coord = lamda;
      if (coord[0] < boxlo[0] || coord[0] >= boxhi[0] || coord[1] < boxlo[1] ||
          coord[1] >= boxhi[1] || coord[2] < boxlo[2] || coord[2] >= boxhi[2])
        continue;
    } else {
      coord = xone;
    }

    //if (list->check_collision(coord)) collision += 1;
  }

  MPI_Allreduce(&collisions, &scalar, 1, MPI_INT, MPI_SUM, world);

  scalar = scalar / ntest;
  return scalar;
}
