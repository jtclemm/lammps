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

// clang-format off
#ifndef LMP_FIX_WALL_GRAN_MODELS_KOKKOS_H
#define LMP_FIX_WALL_GRAN_MODELS_KOKKOS_H

#include "kokkos_type.h"

namespace LAMMPS_NS {
namespace WallGranKokkos {

// coefficients of the classic granular wall interactions,
// converted to the Kokkos floating point type

struct Coeffs {
  KK_FLOAT kn, kt, gamman, gammat, xmu, dt;
};

/* ----------------------------------------------------------------------
   hooke/history interaction of one sphere with one wall contact
   shared by fix wall/gran/kk and fix wall/gran/region/kk

   dx,dy,dz = vector from the contact point on the wall to the sphere center
   rsq = dx*dx + dy*dy + dz*dz, the sphere must overlap the wall
   vwall = velocity of the wall at the contact point
   v, w = velocity and angular velocity of the sphere
   h0,h1,h2 = tangential (shear) history of this contact, updated in place
     if history_update is set
   on return, fout = force and tout = torque to add to the sphere
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void hooke_history(const KK_FLOAT dx, const KK_FLOAT dy, const KK_FLOAT dz, const KK_FLOAT rsq,
                   const KK_FLOAT *vwall, const KK_FLOAT *v, const KK_FLOAT *w,
                   const KK_FLOAT radius, const KK_FLOAT meff,
                   KK_FLOAT &h0, KK_FLOAT &h1, KK_FLOAT &h2,
                   const int history_update, const Coeffs &p,
                   KK_FLOAT *fout, KK_FLOAT *tout)
{
  KK_FLOAT r = Kokkos::sqrt(rsq);
  KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0)/r;
  KK_FLOAT rsqinv = static_cast<KK_FLOAT>(1.0)/rsq;

  // relative translational velocity

  KK_FLOAT vr1 = v[0] - vwall[0];
  KK_FLOAT vr2 = v[1] - vwall[1];
  KK_FLOAT vr3 = v[2] - vwall[2];

  // normal component

  KK_FLOAT vnnr = vr1*dx + vr2*dy + vr3*dz;
  KK_FLOAT vn1 = dx*vnnr * rsqinv;
  KK_FLOAT vn2 = dy*vnnr * rsqinv;
  KK_FLOAT vn3 = dz*vnnr * rsqinv;

  // tangential component

  KK_FLOAT vt1 = vr1 - vn1;
  KK_FLOAT vt2 = vr2 - vn2;
  KK_FLOAT vt3 = vr3 - vn3;

  // relative rotational velocity

  KK_FLOAT wr1 = radius*w[0] * rinv;
  KK_FLOAT wr2 = radius*w[1] * rinv;
  KK_FLOAT wr3 = radius*w[2] * rinv;

  // normal forces = Hookian contact + normal velocity damping

  KK_FLOAT damp = meff*p.gamman*vnnr*rsqinv;
  KK_FLOAT ccel = p.kn*(radius-r)*rinv - damp;

  // relative velocities

  KK_FLOAT vtr1 = vt1 - (dz*wr2-dy*wr3);
  KK_FLOAT vtr2 = vt2 - (dx*wr3-dz*wr1);
  KK_FLOAT vtr3 = vt3 - (dy*wr1-dx*wr2);

  // shear history effects

  if (history_update) {
    h0 += vtr1*p.dt;
    h1 += vtr2*p.dt;
    h2 += vtr3*p.dt;
  }
  KK_FLOAT shrmag = Kokkos::sqrt(h0*h0 + h1*h1 + h2*h2);

  // rotate shear displacements

  KK_FLOAT rsht = h0*dx + h1*dy + h2*dz;
  rsht = rsht*rsqinv;
  if (history_update) {
    h0 -= rsht*dx;
    h1 -= rsht*dy;
    h2 -= rsht*dz;
  }

  // tangential forces = shear + tangential velocity damping

  KK_FLOAT fs1 = - (p.kt*h0 + meff*p.gammat*vtr1);
  KK_FLOAT fs2 = - (p.kt*h1 + meff*p.gammat*vtr2);
  KK_FLOAT fs3 = - (p.kt*h2 + meff*p.gammat*vtr3);

  // rescale frictional displacements and forces if needed

  KK_FLOAT fs = Kokkos::sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  KK_FLOAT fn = p.xmu * Kokkos::fabs(ccel*r);

  if (fs > fn) {
    if (shrmag != static_cast<KK_FLOAT>(0.0)) {
      h0 = (fn/fs) * (h0 + meff*p.gammat*vtr1/p.kt) - meff*p.gammat*vtr1/p.kt;
      h1 = (fn/fs) * (h1 + meff*p.gammat*vtr2/p.kt) - meff*p.gammat*vtr2/p.kt;
      h2 = (fn/fs) * (h2 + meff*p.gammat*vtr3/p.kt) - meff*p.gammat*vtr3/p.kt;
      fs1 *= fn/fs;
      fs2 *= fn/fs;
      fs3 *= fn/fs;
    } else fs1 = fs2 = fs3 = 0.0;
  }

  // forces & torques

  fout[0] = dx*ccel + fs1;
  fout[1] = dy*ccel + fs2;
  fout[2] = dz*ccel + fs3;

  KK_FLOAT tor1 = rinv * (dy*fs3 - dz*fs2);
  KK_FLOAT tor2 = rinv * (dz*fs1 - dx*fs3);
  KK_FLOAT tor3 = rinv * (dx*fs2 - dy*fs1);
  tout[0] = -(radius*tor1);
  tout[1] = -(radius*tor2);
  tout[2] = -(radius*tor3);
}

}    // namespace WallGranKokkos
}    // namespace LAMMPS_NS

#endif
