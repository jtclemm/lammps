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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "pair_lj_cut_coul_cut_dpd_tstat.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairLJCutCoulCutDPDTstat::PairLJCutCoulCutDPDTstat(LAMMPS *lmp) : Pair(lmp)
{
  born_matrix_enable = 1;
  single_enable = 0;
  writedata = 1;
  random = nullptr;
}

/* ---------------------------------------------------------------------- */

PairLJCutCoulCutDPDTstat::~PairLJCutCoulCutDPDTstat()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(cut_coul);
    memory->destroy(cut_coulsq);
    memory->destroy(cut_dpd);
    memory->destroy(cut_dpdsq);
    memory->destroy(cut_dpdinv);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(sigma2);
    memory->destroy(gamma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
  }

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulCutDPDTstat::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul,fpair;
  double vxtmp, vytmp, vztmp, delvx, delvy, delvz;
  double rsq, r2inv, r6inv, forcecoul, forcelj, force_dpd, factor_coul, factor_lj;
  double r, rinv, dot, wd, randnum, factor_dpd;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  // adjust sigma if target T is changing

  if (t_start != t_stop) {
    double delta = update->ntimestep - update->beginstep;
    if (delta != 0.0) delta /= update->endstep - update->beginstep;
    temperature = t_start + delta * (t_stop - t_start);
    double boltz = force->boltz;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        sigma2[i][j] = sigma2[j][i] = sqrt(2.0 * boltz * temperature * gamma[i][j]);
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0 / sqrt(update->dt);
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        rinv = 1.0 / r;
        r2inv = rinv * rinv;

        if (rsq < cut_coulsq[itype][jtype])
          forcecoul = qqrd2e * qtmp * q[j] * rinv * r2inv;
        else
          forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv * r2inv * r2inv;
          forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]) * r2inv;
        } else
          forcelj = 0.0;

        if (rsq < cut_dpdsq[itype][jtype]) {
          delvx = vxtmp - v[j][0];
          delvy = vytmp - v[j][1];
          delvz = vztmp - v[j][2];
          dot = delx * delvx + dely * delvy + delz * delvz;
          wd = 1.0 - r * cut_dpdinv[itype][jtype];
          randnum = random->gaussian();

          // drag force = -gamma * wd^2 * (delx dot delv) / r
          // random force = sigma * wd * rnd * dtinvsqrt;

          force_dpd = -gamma[itype][jtype] * wd * wd * dot * r2inv;
          force_dpd += sigma2[itype][jtype] * wd * randnum * dtinvsqrt * rinv;
        }

        fpair = (factor_coul * forcecoul + factor_lj * (forcelj + force_dpd));

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) {
          if (rsq < cut_coulsq[itype][jtype])
            ecoul = factor_coul * qqrd2e * qtmp * q[j] * sqrt(r2inv);
          else
            ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
            evdwl *= factor_lj;
          } else
            evdwl = 0.0;
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, ecoul, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCutCoulCutDPDTstat::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n, n, "pair:cutsq");

  memory->create(cut_lj, n, n, "pair:cut_lj");
  memory->create(cut_ljsq, n, n, "pair:cut_ljsq");
  memory->create(cut_coul, n, n, "pair:cut_coul");
  memory->create(cut_coulsq, n, n, "pair:cut_coulsq");
  memory->create(cut_dpd, n, n, "pair:cut_dpd");
  memory->create(cut_dpdsq, n, n, "pair:cut_dpdsq");
  memory->create(cut_dpdinv, n, n, "pair:cut_dpdinv");
  memory->create(epsilon, n, n, "pair:epsilon");
  memory->create(sigma, n, n, "pair:sigma");
  memory->create(sigma2, n, n, "pair:sigma2");
  memory->create(gamma, n, n, "pair:gamma");
  memory->create(lj1, n, n, "pair:lj1");
  memory->create(lj2, n, n, "pair:lj2");
  memory->create(lj3, n, n, "pair:lj3");
  memory->create(lj4, n, n, "pair:lj4");
  memory->create(offset, n, n, "pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJCutCoulCutDPDTstat::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR, "Illegal pair_style command");

  t_start = utils::numeric(FLERR,arg[0],false,lmp);
  t_stop = utils::numeric(FLERR,arg[1],false,lmp);
  seed = utils::inumeric(FLERR,arg[2],false,lmp);

  temperature = t_start;

  // initialize Marsaglia RNG with processor-unique seed

  if (seed <= 0) error->all(FLERR,"Illegal pair_style command");
  delete random;
  random = new RanMars(lmp, seed + comm->me);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCutCoulCutDPDTstat::coeff(int narg, char **arg)
{
  if (narg != 8) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);
  double gamma_one = utils::numeric(FLERR,arg[4],false,lmp);

  double cut_lj_one = utils::numeric(FLERR, arg[5], false, lmp);
  double cut_coul_one = utils::numeric(FLERR, arg[6], false, lmp);
  double cut_dpd_one = utils::numeric(FLERR, arg[7], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      gamma[i][j] = gamma_one;
      cut_lj[i][j] = cut_lj_one;
      cut_coul[i][j] = cut_coul_one;
      cut_dpd[i][j] = cut_dpd_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutCoulCutDPDTstat::init_style()
{
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair lj/cut/coul/cut/dpd/tstat requires ghost atoms store velocity");

  if (!atom->q_flag) error->all(FLERR, "Pair style lj/cut/coul/cut/dpd/tstat requires atom attribute q");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0) error->warning(FLERR,
      "Pair lj/cut/coul/cut/dpd/tstat needs newton pair on for momentum conservation");

  neighbor->request(this, instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutCoulCutDPDTstat::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], sigma[i][i], sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i], sigma[j][j]);
    cut_lj[i][j] = mix_distance(cut_lj[i][i], cut_lj[j][j]);
    cut_coul[i][j] = mix_distance(cut_coul[i][i], cut_coul[j][j]);
    cut_dpd[i][j] = mix_distance(cut_dpd[i][i], cut_dpd[j][j]);
  }

  double cut = MAX(MAX(cut_lj[i][j], cut_coul[i][j]), cut_dpd[i][j]);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];
  cut_coulsq[i][j] = cut_coul[i][j] * cut_coul[i][j];
  cut_dpdsq[i][j] = cut_dpd[i][j] * cut_dpd[i][j];
  cut_dpdinv[i][j] = 1.0 / cut_dpd[i][j];

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);

  if (offset_flag && (cut_lj[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut_lj[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 12.0) - pow(ratio, 6.0));
  } else
    offset[i][j] = 0.0;

  sigma2[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);

  cut_ljsq[j][i] = cut_ljsq[i][j];
  cut_coulsq[j][i] = cut_coulsq[i][j];
  cut_dpdsq[j][i] = cut_dpdsq[i][j];
  cut_dpdinv[j][i] = cut_dpdinv[i][j];

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];
  sigma2[j][i] = sigma2[i][j];
  gamma[j][i] = gamma[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2], all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count, all, 2, MPI_DOUBLE, MPI_SUM, world);

    double sig2 = sigma[i][j] * sigma[i][j];
    double sig6 = sig2 * sig2 * sig2;
    double rc3 = cut_lj[i][j] * cut_lj[i][j] * cut_lj[i][j];
    double rc6 = rc3 * rc3;
    double rc9 = rc3 * rc6;
    double prefactor = 8.0 * MY_PI * all[0] * all[1] * epsilon[i][j] * sig6 / (9.0 * rc9);
    etail_ij = prefactor * (sig6 - 3.0 * rc6);
    ptail_ij = 2.0 * prefactor * (2.0 * sig6 - 3.0 * rc6);
  }

  return cut;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutCoulCutDPDTstat::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&sigma2[i][j], sizeof(double), 1, fp);
        fwrite(&gamma[i][j], sizeof(double), 1, fp);
        fwrite(&cut_lj[i][j], sizeof(double), 1, fp);
        fwrite(&cut_coul[i][j], sizeof(double), 1, fp);
        fwrite(&cut_dpd[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutCoulCutDPDTstat::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &sigma2[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &gamma[i][j], sizeof(double), 1, fp,nullptr, error);
          utils::sfread(FLERR, &cut_lj[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut_coul[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut_dpd[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma2[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&gamma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut_lj[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut_coul[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut_dpd[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutCoulCutDPDTstat::write_restart_settings(FILE *fp)
{
  fwrite(&t_start, sizeof(double), 1, fp);
  fwrite(&t_stop, sizeof(double), 1, fp);
  fwrite(&seed, sizeof(int), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  fwrite(&tail_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutCoulCutDPDTstat::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &t_start,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR, &t_stop,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR, &seed, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&t_start, 1 ,MPI_DOUBLE, 0, world);
  MPI_Bcast(&t_stop, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&seed,1, MPI_INT, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);

  temperature = t_start;

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJCutCoulCutDPDTstat::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g %g %g %g %g %g\n", i, epsilon[i][i], sigma[i][i], gamma[i][i], cut_lj[i][i], cut_coul[i][i], cut_dpd[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJCutCoulCutDPDTstat::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g %g\n", i, j, epsilon[i][j], sigma[i][j], gamma[i][j], cut_lj[i][j], cut_coul[i][j], cut_dpd[i][j]);
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulCutDPDTstat::born_matrix(int i, int j, int itype, int jtype, double rsq,
                                   double factor_coul, double factor_lj, double &dupair,
                                   double &du2pair)
{
  double rinv, r2inv, r3inv, r6inv;
  double du_lj, du2_lj, du_coul, du2_coul;

  double *q = atom->q;
  double qqrd2e = force->qqrd2e;

  r2inv = 1.0 / rsq;
  rinv = sqrt(r2inv);
  r3inv = r2inv * rinv;
  r6inv = r2inv * r2inv * r2inv;

  // Reminder: lj1 = 48*e*s^12, lj2 = 24*e*s^6

  du_lj = r6inv * rinv * (lj2[itype][jtype] - lj1[itype][jtype] * r6inv);
  du2_lj = r6inv * r2inv * (13 * lj1[itype][jtype] * r6inv - 7 * lj2[itype][jtype]);

  du_coul = -qqrd2e * q[i] * q[j] * r2inv;
  du2_coul = 2.0 * qqrd2e * q[i] * q[j] * r3inv;

  dupair = factor_lj * du_lj + factor_coul * du_coul;
  du2pair = factor_lj * du2_lj + factor_coul * du2_coul;
}

/* ---------------------------------------------------------------------- */

void *PairLJCutCoulCutDPDTstat::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "cut_coul") == 0) return (void *) cut_coul;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;
  return nullptr;
}
