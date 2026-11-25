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

#ifdef REGION_CLASS
// clang-format off
RegionStyle(cylinder/kk,RegCylinderKokkos<LMPDeviceType>);
RegionStyle(cylinder/kk/device,RegCylinderKokkos<LMPDeviceType>);
RegionStyle(cylinder/kk/host,RegCylinderKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_REGION_CYLINDER_KOKKOS_H
#define LMP_REGION_CYLINDER_KOKKOS_H

#include "region_cylinder.h"

#include "kokkos_base.h"
#include "kokkos_type.h"
#include "region_remap_kokkos.h"

namespace LAMMPS_NS {

struct TagRegCylinderMatchAll{};

template<class DeviceType>
class RegCylinderKokkos : public RegCylinder, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  // maximum number of contacts a particle can have with this region, i.e. the
  // length of the contact list the caller has to provide.  matches RegCylinder::cmax

  static constexpr int MAXCONTACT = 3;

  RegCylinderKokkos(class LAMMPS *, int, char **);
  // checked conversion from a base class pointer, see the .cpp file
  static RegCylinderKokkos<DeviceType> *cast(class Region *);

  void match_all_kokkos(int, DAT::tdual_int_1d) override;

  void prematch() override
  {
    Region::prematch();
    boxremap.capture(domain);
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRegCylinderMatchAll, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int match_kokkos(double x, double y, double z) const
  {
    boxremap.remap(x,y,z);
    if (dynamic) inverse_transform(x,y,z);
    if (openflag) return 1;
    return !(k_inside(x,y,z) ^ interior);
  }

  // device version of Region::surface()

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int surface_kokkos(double x, double y, double z, double cutoff, Contact *contact) const
  {
    int ncontact;
    double xs, ys, zs;
    double xnear[3], xorig[3];

    boxremap.remap(x, y, z);

    xorig[0] = x; xorig[1] = y; xorig[2] = z;
    if (dynamic)
      inverse_transform(x, y, z);

    xnear[0] = x; xnear[1] = y; xnear[2] = z;

    if (!openflag) {
      if (interior)
        ncontact = surface_interior_kokkos(xnear, cutoff, contact);
      else
        ncontact = surface_exterior_kokkos(xnear, cutoff, contact);
    } else {
      // most of the time, one of surface_int/ext() will return 0
      //   however, when exactly on top of a periodic boundary
      //   both could return 1, so run exterior then interior
      ncontact = surface_exterior_kokkos(xnear, cutoff, contact);
      if (ncontact == 0) ncontact = surface_interior_kokkos(xnear, cutoff, contact);
    }

    if (rotateflag && ncontact) {
      for (int i = 0; i < ncontact; i++) {
        xs = xnear[0] - contact[i].delx;
        ys = xnear[1] - contact[i].dely;
        zs = xnear[2] - contact[i].delz;
        forward_transform(xs, ys, zs);
        contact[i].delx = xorig[0] - xs;
        contact[i].dely = xorig[1] - ys;
        contact[i].delz = xorig[2] - zs;
      }
    }

    return ncontact;
  }

  // device version of RegCylinder::velocity_contact_shape()
  // add velocity due to shape change to wall velocity

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void velocity_contact_shape_kokkos(double *vwall, const double *xc, double dt) const
  {
    double delx, dely, delz;    // Displacement of contact point in x,y,z
    if (axis == 'x') {
      delx = 0;
      dely = (xc[1] - xcenter[1]) * (1 - rprev / radius);
      delz = (xc[2] - xcenter[2]) * (1 - rprev / radius);
    } else if (axis == 'y') {
      delx = (xc[0] - xcenter[0]) * (1 - rprev / radius);
      dely = 0;
      delz = (xc[2] - xcenter[2]) * (1 - rprev / radius);
    } else {
      delx = (xc[0] - xcenter[0]) * (1 - rprev / radius);
      dely = (xc[1] - xcenter[1]) * (1 - rprev / radius);
      delz = 0;
    }
    vwall[0] += delx / dt;
    vwall[1] += dely / dt;
    vwall[2] += delz / dt;
  }

 private:
  int groupbit;
  RegionRemapKokkos boxremap;
  typename AT::t_int_1d d_match;
  typename AT::t_kkfloat_1d_3_lr_randomread d_x;
  typename AT::t_int_1d_randomread d_mask;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int k_inside(double x, double y, double z) const
  {
    double del1, del2, dist;

    if (axis == 'x') {
      del1 = y - c1;
      del2 = z - c2;
      dist = sqrt(del1 * del1 + del2 * del2);
      if (dist <= radius && x >= lo && x <= hi) return 1;
    } else if (axis == 'y') {
      del1 = x - c1;
      del2 = z - c2;
      dist = sqrt(del1 * del1 + del2 * del2);
      if (dist <= radius && y >= lo && y <= hi) return 1;
    } else {
      del1 = x - c1;
      del2 = y - c2;
      dist = sqrt(del1 * del1 + del2 * del2);
      if (dist <= radius && z >= lo && z <= hi) return 1;
    }
    return 0;
  }

  static constexpr double BIG = 1.0e200;

  // device versions of the CPU surface_interior() and surface_exterior()
  // of RegCylinder, see there

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int surface_interior_kokkos(double *x, double cutoff, Contact *contact) const
  {
    double del1, del2, r, delta;

    int n = 0;

    if (axis == 'x') {
      del1 = x[1] - c1;
      del2 = x[2] - c2;
      r = Kokkos::sqrt(del1 * del1 + del2 * del2);

      // x is exterior to cylinder

      if (r > radius || x[0] < lo || x[0] > hi) return 0;

      // x is interior to cylinder or on its surface

      delta = radius - r;
      if (delta < cutoff && r > 0.0 && !open_faces[2]) {
        contact[n].r = delta;
        contact[n].delx = 0.0;
        contact[n].dely = del1 * (1.0 - radius / r);
        contact[n].delz = del2 * (1.0 - radius / r);
        contact[n].radius = -2.0 * radius;
        contact[n].iwall = 2;
        contact[n].varflag = 1;
        n++;
      }
      delta = x[0] - lo;
      if (delta < cutoff && !open_faces[0]) {
        contact[n].r = delta;
        contact[n].delx = delta;
        contact[n].dely = contact[n].delz = 0.0;
        contact[n].radius = 0;
        contact[n].iwall = 0;
        contact[n].varflag = 0;
        n++;
      }
      delta = hi - x[0];
      if (delta < cutoff && !open_faces[1]) {
        contact[n].r = delta;
        contact[n].delx = -delta;
        contact[n].dely = contact[n].delz = 0.0;
        contact[n].radius = 0;
        contact[n].iwall = 1;
        contact[n].varflag = 0;
        n++;
      }

    } else if (axis == 'y') {
      del1 = x[0] - c1;
      del2 = x[2] - c2;
      r = Kokkos::sqrt(del1 * del1 + del2 * del2);

      // y is exterior to cylinder

      if (r > radius || x[1] < lo || x[1] > hi) return 0;

      // y is interior to cylinder or on its surface

      delta = radius - r;
      if (delta < cutoff && r > 0.0 && !open_faces[2]) {
        contact[n].r = delta;
        contact[n].delx = del1 * (1.0 - radius / r);
        contact[n].dely = 0.0;
        contact[n].delz = del2 * (1.0 - radius / r);
        contact[n].radius = -2.0 * radius;
        contact[n].iwall = 2;
        contact[n].varflag = 1;
        n++;
      }
      delta = x[1] - lo;
      if (delta < cutoff && !open_faces[0]) {
        contact[n].r = delta;
        contact[n].dely = delta;
        contact[n].delx = contact[n].delz = 0.0;
        contact[n].radius = 0;
        contact[n].iwall = 0;
        contact[n].varflag = 0;
        n++;
      }
      delta = hi - x[1];
      if (delta < cutoff && !open_faces[1]) {
        contact[n].r = delta;
        contact[n].dely = -delta;
        contact[n].delx = contact[n].delz = 0.0;
        contact[n].radius = 0;
        contact[n].iwall = 1;
        contact[n].varflag = 0;
        n++;
      }

    } else {
      del1 = x[0] - c1;
      del2 = x[1] - c2;
      r = Kokkos::sqrt(del1 * del1 + del2 * del2);

      // z is exterior to cylinder

      if (r > radius || x[2] < lo || x[2] > hi) return 0;

      // z is interior to cylinder or on its surface

      delta = radius - r;
      if (delta < cutoff && r > 0.0 && !open_faces[2]) {
        contact[n].r = delta;
        contact[n].delx = del1 * (1.0 - radius / r);
        contact[n].dely = del2 * (1.0 - radius / r);
        contact[n].delz = 0.0;
        contact[n].radius = -2.0 * radius;
        contact[n].iwall = 2;
        contact[n].varflag = 1;
        n++;
      }
      delta = x[2] - lo;
      if (delta < cutoff && !open_faces[0]) {
        contact[n].r = delta;
        contact[n].delz = delta;
        contact[n].delx = contact[n].dely = 0.0;
        contact[n].radius = 0;
        contact[n].iwall = 0;
        contact[n].varflag = 0;
        n++;
      }
      delta = hi - x[2];
      if (delta < cutoff && !open_faces[1]) {
        contact[n].r = delta;
        contact[n].delz = -delta;
        contact[n].delx = contact[n].dely = 0.0;
        contact[n].radius = 0;
        contact[n].iwall = 1;
        contact[n].varflag = 0;
        n++;
      }
    }

    return n;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int surface_exterior_kokkos(double *x, double cutoff, Contact *contact) const
  {
    double del1, del2, r;
    double xp, yp, zp;
    double dx, dr, dr2, d2, d2prev;

    // with all three faces open there is no surface left to contact

    if (open_faces[0] && open_faces[1] && open_faces[2]) return 0;

    // radius of curvature for granular
    // 0 for flat surfaces (infinite case), 2*radius for curved portion

    double crad = 0.0;
    int varflag = 0;

    if (axis == 'x') {
      del1 = x[1] - c1;
      del2 = x[2] - c2;
      r = Kokkos::sqrt(del1 * del1 + del2 * del2);

      // x is far enough from cylinder that there is no contact
      // x is interior to cylinder

      if (r >= radius + cutoff || x[0] <= lo - cutoff || x[0] >= hi + cutoff) return 0;
      if (r < radius && x[0] > lo && x[0] < hi) return 0;

      // x is exterior to cylinder or on its surface
      // xp,yp,zp = point on surface of cylinder that x is closest to
      //            could be edge of cylinder
      // do not add contact point if r >= cutoff

      d2prev = BIG;
      if (!openflag) {
        if (r > radius) {
          yp = c1 + del1 * radius / r;
          zp = c2 + del2 * radius / r;
          crad = 2.0 * radius;
          varflag = 1;
        } else {
          yp = x[1];
          zp = x[2];
        }
        if (x[0] < lo)
          xp = lo;
        else if (x[0] > hi)
          xp = hi;
        else
          xp = x[0];

      } else {

        // closest point on curved surface

        dr = r - radius;
        dr2 = dr * dr;
        if (!open_faces[2]) {
          yp = c1 + del1 * radius / r;
          zp = c2 + del2 * radius / r;
          if (x[0] < lo) {
            dx = lo - x[0];
            xp = lo;
          } else if (x[0] > hi) {
            dx = x[0] - hi;
            xp = hi;
          } else {
            dx = 0;
            xp = x[0];
          }
          d2 = d2prev = dr2 + dx * dx;
          if (r > radius) {
            crad = 2.0 * radius;
            varflag = 1;
          }
        }

        // closest point on bottom cap

        if (!open_faces[0]) {
          dx = lo - x[0];
          if (r < radius)
            d2 = dx * dx;
          else
            d2 = dr2 + dx * dx;
          if (d2 < d2prev) {
            xp = lo;
            if (r < radius) {
              yp = x[1];
              zp = x[2];
            } else {
              yp = c1 + del1 * radius / r;
              zp = c2 + del2 * radius / r;
            }
            d2prev = d2;
          }
        }

        // closest point on top cap

        if (!open_faces[1]) {
          dx = hi - x[0];
          if (r < radius)
            d2 = dx * dx;
          else
            d2 = dr2 + dx * dx;
          if (d2 < d2prev) {
            xp = hi;
            if (r < radius) {
              yp = x[1];
              zp = x[2];
            } else {
              yp = c1 + del1 * radius / r;
              zp = c2 + del2 * radius / r;
            }
          }
        }
      }

      add_contact_kokkos(0, x, xp, yp, zp, contact);
      contact[0].radius = crad;
      contact[0].varflag = varflag;
      contact[0].iwall = 0;
      if (contact[0].r < cutoff) return 1;
      return 0;

    } else if (axis == 'y') {
      del1 = x[0] - c1;
      del2 = x[2] - c2;
      r = Kokkos::sqrt(del1 * del1 + del2 * del2);

      // y is far enough from cylinder that there is no contact
      // y is interior to cylinder

      if (r >= radius + cutoff || x[1] <= lo - cutoff || x[1] >= hi + cutoff) return 0;
      if (r < radius && x[1] > lo && x[1] < hi) return 0;

      // y is exterior to cylinder or on its surface
      // xp,yp,zp = point on surface of cylinder that x is closest to
      //            could be edge of cylinder
      // do not add contact point if r >= cutoff

      d2prev = BIG;
      if (!openflag) {
        if (r > radius) {
          xp = c1 + del1 * radius / r;
          zp = c2 + del2 * radius / r;
          crad = 2.0 * radius;
          varflag = 1;
        } else {
          xp = x[0];
          zp = x[2];
        }
        if (x[1] < lo)
          yp = lo;
        else if (x[1] > hi)
          yp = hi;
        else
          yp = x[1];

      } else {

        // closest point on curved surface

        dr = r - radius;
        dr2 = dr * dr;
        if (!open_faces[2]) {
          xp = c1 + del1 * radius / r;
          zp = c2 + del2 * radius / r;
          if (x[1] < lo) {
            dx = lo - x[1];
            yp = lo;
          } else if (x[1] > hi) {
            dx = x[1] - hi;
            yp = hi;
          } else {
            dx = 0;
            yp = x[1];
          }
          d2 = d2prev = dr2 + dx * dx;
          if (r > radius) {
            crad = 2.0 * radius;
            varflag = 1;
          }
        }

        // closest point on bottom cap

        if (!open_faces[0]) {
          dx = lo - x[1];
          if (r < radius)
            d2 = dx * dx;
          else
            d2 = dr2 + dx * dx;
          if (d2 < d2prev) {
            yp = lo;
            if (r < radius) {
              xp = x[0];
              zp = x[2];
            } else {
              xp = c1 + del1 * radius / r;
              zp = c2 + del2 * radius / r;
            }
            d2prev = d2;
          }
        }

        // closest point on top cap

        if (!open_faces[1]) {
          dx = hi - x[1];
          if (r < radius)
            d2 = dx * dx;
          else
            d2 = dr2 + dx * dx;
          if (d2 < d2prev) {
            yp = hi;
            if (r < radius) {
              xp = x[0];
              zp = x[2];
            } else {
              xp = c1 + del1 * radius / r;
              zp = c2 + del2 * radius / r;
            }
          }
        }
      }

      add_contact_kokkos(0, x, xp, yp, zp, contact);
      contact[0].radius = crad;
      contact[0].varflag = varflag;
      contact[0].iwall = 0;
      if (contact[0].r < cutoff) return 1;
      return 0;

    } else {
      del1 = x[0] - c1;
      del2 = x[1] - c2;
      r = Kokkos::sqrt(del1 * del1 + del2 * del2);

      // z is far enough from cylinder that there is no contact
      // z is interior to cylinder

      if (r >= radius + cutoff || x[2] <= lo - cutoff || x[2] >= hi + cutoff) return 0;
      if (r < radius && x[2] > lo && x[2] < hi) return 0;

      // z is exterior to cylinder or on its surface
      // xp,yp,zp = point on surface of cylinder that x is closest to
      //            could be edge of cylinder
      // do not add contact point if r >= cutoff

      d2prev = BIG;
      if (!openflag) {
        if (r > radius) {
          xp = c1 + del1 * radius / r;
          yp = c2 + del2 * radius / r;
          crad = 2.0 * radius;
          varflag = 1;
        } else {
          xp = x[0];
          yp = x[1];
        }
        if (x[2] < lo)
          zp = lo;
        else if (x[2] > hi)
          zp = hi;
        else
          zp = x[2];

      } else {

        // closest point on curved surface

        dr = r - radius;
        dr2 = dr * dr;
        if (!open_faces[2]) {
          xp = c1 + del1 * radius / r;
          yp = c2 + del2 * radius / r;
          if (x[2] < lo) {
            dx = lo - x[2];
            zp = lo;
          } else if (x[2] > hi) {
            dx = x[2] - hi;
            zp = hi;
          } else {
            dx = 0;
            zp = x[2];
          }
          d2prev = dr2 + dx * dx;
          if (r > radius) {
            crad = 2.0 * radius;
            varflag = 1;
          }
        }

        // closest point on bottom cap

        if (!open_faces[0]) {
          dx = lo - x[2];
          if (r < radius)
            d2 = dx * dx;
          else
            d2 = dr2 + dx * dx;
          if (d2 < d2prev) {
            zp = lo;
            if (r < radius) {
              xp = x[0];
              yp = x[1];
            } else {
              xp = c1 + del1 * radius / r;
              yp = c2 + del2 * radius / r;
            }
            d2prev = d2;
          }
        }

        // closest point on top cap

        if (!open_faces[1]) {
          dx = hi - x[2];
          if (r < radius)
            d2 = dx * dx;
          else
            d2 = dr2 + dx * dx;
          if (d2 < d2prev) {
            zp = hi;
            if (r < radius) {
              xp = x[0];
              yp = x[1];
            } else {
              xp = c1 + del1 * radius / r;
              yp = c2 + del2 * radius / r;
            }
          }
        }
      }

      add_contact_kokkos(0, x, xp, yp, zp, contact);
      contact[0].radius = crad;
      contact[0].varflag = varflag;
      contact[0].iwall = 0;
      if (contact[0].r < cutoff) return 1;
      return 0;
    }
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void add_contact_kokkos(int n, double *x, double xp, double yp, double zp,
                          Contact *contact) const
  {
    double delx = x[0] - xp;
    double dely = x[1] - yp;
    double delz = x[2] - zp;
    contact[n].r = Kokkos::sqrt(delx * delx + dely * dely + delz * delz);
    contact[n].radius = 0;
    contact[n].delx = delx;
    contact[n].dely = dely;
    contact[n].delz = delz;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void forward_transform(double &x, double &y, double &z) const
  {
    if (rotateflag) rotate(x, y, z, theta);
    if (moveflag) {
      x += dx;
      y += dy;
      z += dz;
    }
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void inverse_transform(double &x, double &y, double &z) const
  {
    if (moveflag) {
      x -= dx;
      y -= dy;
      z -= dz;
    }
    if (rotateflag) rotate(x,y,z,-theta);
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void rotate(double &x, double &y, double &z, double angle) const
  {
    double a[3],b[3],c[3],d[3],disp[3];

    double sine = sin(angle);
    double cosine = cos(angle);
    d[0] = x - point[0];
    d[1] = y - point[1];
    d[2] = z - point[2];
    double x0dotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
    c[0] = x0dotr * runit[0];
    c[1] = x0dotr * runit[1];
    c[2] = x0dotr * runit[2];
    a[0] = d[0] - c[0];
    a[1] = d[1] - c[1];
    a[2] = d[2] - c[2];
    b[0] = runit[1]*a[2] - runit[2]*a[1];
    b[1] = runit[2]*a[0] - runit[0]*a[2];
    b[2] = runit[0]*a[1] - runit[1]*a[0];
    disp[0] = a[0]*cosine  + b[0]*sine;
    disp[1] = a[1]*cosine  + b[1]*sine;
    disp[2] = a[2]*cosine  + b[2]*sine;
    x = point[0] + c[0] + disp[0];
    y = point[1] + c[1] + disp[1];
    z = point[2] + c[2] + disp[2];
  }

};

}

#endif
#endif
