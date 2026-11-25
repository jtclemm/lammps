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
RegionStyle(cone/kk,RegConeKokkos<LMPDeviceType>);
RegionStyle(cone/kk/device,RegConeKokkos<LMPDeviceType>);
RegionStyle(cone/kk/host,RegConeKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_REGION_CONE_KOKKOS_H
#define LMP_REGION_CONE_KOKKOS_H

#include "region_cone.h"

#include "kokkos_base.h"
#include "kokkos_type.h"
#include "region_remap_kokkos.h"

namespace LAMMPS_NS {

struct TagRegConeMatchAll{};

template<class DeviceType>
class RegConeKokkos : public RegCone, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  // maximum number of contacts a particle can have with this region, i.e. the
  // length of the contact list the caller has to provide.  matches RegCone::cmax

  static constexpr int MAXCONTACT = 3;

  RegConeKokkos(class LAMMPS *, int, char **);
  // checked conversion from a base class pointer, see the .cpp file
  static RegConeKokkos<DeviceType> *cast(class Region *);

  void match_all_kokkos(int, DAT::tdual_int_1d) override;

  void prematch() override
  {
    Region::prematch();
    boxremap.capture(domain);
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRegConeMatchAll, const int&) const;

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

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void velocity_contact_shape_kokkos(double * /*vwall*/, const double * /*xc*/,
                                     double /*dt*/) const
  {
    // like RegCone, no wall velocity is added for a changing shape
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
    double del1, del2, dist, currentradius;

    if (axis == 'x') {
      del1 = y - c1;
      del2 = z - c2;
      dist = sqrt(del1 * del1 + del2 * del2);
      currentradius = radiuslo + (x - lo) * (radiushi - radiuslo) / (hi - lo);
      if (dist <= currentradius && x >= lo && x <= hi) return 1;
    } else if (axis == 'y') {
      del1 = x - c1;
      del2 = z - c2;
      dist = sqrt(del1 * del1 + del2 * del2);
      currentradius = radiuslo + (y - lo) * (radiushi - radiuslo) / (hi - lo);
      if (dist <= currentradius && y >= lo && y <= hi) return 1;
    } else {
      del1 = x - c1;
      del2 = y - c2;
      dist = sqrt(del1 * del1 + del2 * del2);
      currentradius = radiuslo + (z - lo) * (radiushi - radiuslo) / (hi - lo);
      if (dist <= currentradius && z >= lo && z <= hi) return 1;
    }
    return 0;
  }

  static constexpr double BIG = 1.0e200;

  // device versions of the CPU surface_interior() and surface_exterior()
  // of RegCone, see there

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int surface_interior_kokkos(double *x, double cutoff, Contact *contact) const
  {
    double del1, del2, r, currentradius, delx, dely, delz, dist, delta;
    double surflo[3], surfhi[3], xs[3];

    int n = 0;

    if (axis == 'x') {
      del1 = x[1] - c1;
      del2 = x[2] - c2;
      r = Kokkos::sqrt(del1 * del1 + del2 * del2);
      currentradius = radiuslo + (x[0] - lo) * (radiushi - radiuslo) / (hi - lo);

      // x is exterior to cone

      if (r > currentradius || x[0] < lo || x[0] > hi) return 0;

      // x is interior to cone or on its surface
      // surflo = pt on outer circle of bottom end plane, same dir as x vs axis
      // surfhi = pt on outer circle of top end plane, same dir as x vs axis

      if (r > 0.0 && !open_faces[2]) {
        surflo[0] = lo;
        surflo[1] = c1 + del1 * radiuslo / r;
        surflo[2] = c2 + del2 * radiuslo / r;
        surfhi[0] = hi;
        surfhi[1] = c1 + del1 * radiushi / r;
        surfhi[2] = c2 + del2 * radiushi / r;
        point_on_line_segment_kokkos(surflo, surfhi, x, xs);
        delx = x[0] - xs[0];
        dely = x[1] - xs[1];
        delz = x[2] - xs[2];
        dist = Kokkos::sqrt(delx * delx + dely * dely + delz * delz);
        if (dist < cutoff) {
          contact[n].r = dist;
          contact[n].delx = delx;
          contact[n].dely = dely;
          contact[n].delz = delz;
          contact[n].radius = -2.0 * (radiuslo + (xs[0] - lo) * (radiushi - radiuslo) / (hi - lo));
          contact[n].iwall = 2;
          n++;
        }
      }

      delta = x[0] - lo;
      if (delta < cutoff && !open_faces[0]) {
        contact[n].r = delta;
        contact[n].delx = delta;
        contact[n].dely = contact[n].delz = 0.0;
        contact[n].radius = 0;
        contact[n].iwall = 0;
        n++;
      }
      delta = hi - x[0];
      if (delta < cutoff && !open_faces[1]) {
        contact[n].r = delta;
        contact[n].delx = -delta;
        contact[n].dely = contact[n].delz = 0.0;
        contact[n].radius = 0;
        contact[n].iwall = 1;
        n++;
      }

    } else if (axis == 'y') {
      del1 = x[0] - c1;
      del2 = x[2] - c2;
      r = Kokkos::sqrt(del1 * del1 + del2 * del2);
      currentradius = radiuslo + (x[1] - lo) * (radiushi - radiuslo) / (hi - lo);

      // y is exterior to cone

      if (r > currentradius || x[1] < lo || x[1] > hi) return 0;

      // y is interior to cone or on its surface
      // surflo = pt on outer circle of bottom end plane, same dir as y vs axis
      // surfhi = pt on outer circle of top end plane, same dir as y vs axis

      if (r > 0.0 && !open_faces[2]) {
        surflo[0] = c1 + del1 * radiuslo / r;
        surflo[1] = lo;
        surflo[2] = c2 + del2 * radiuslo / r;
        surfhi[0] = c1 + del1 * radiushi / r;
        surfhi[1] = hi;
        surfhi[2] = c2 + del2 * radiushi / r;
        point_on_line_segment_kokkos(surflo, surfhi, x, xs);
        delx = x[0] - xs[0];
        dely = x[1] - xs[1];
        delz = x[2] - xs[2];
        dist = Kokkos::sqrt(delx * delx + dely * dely + delz * delz);
        if (dist < cutoff) {
          contact[n].r = dist;
          contact[n].delx = delx;
          contact[n].dely = dely;
          contact[n].delz = delz;
          contact[n].iwall = 2;
          contact[n].radius = -2.0 * (radiuslo + (xs[1] - lo) * (radiushi - radiuslo) / (hi - lo));
          n++;
        }
      }

      delta = x[1] - lo;
      if (delta < cutoff && !open_faces[0]) {
        contact[n].r = delta;
        contact[n].delz = delta;
        contact[n].delx = contact[n].dely = 0.0;
        contact[n].iwall = 0;
        contact[n].radius = 0;
        n++;
      }
      delta = hi - x[1];
      if (delta < cutoff && !open_faces[1]) {
        contact[n].r = delta;
        contact[n].delz = -delta;
        contact[n].delx = contact[n].dely = 0.0;
        contact[n].iwall = 1;
        contact[n].radius = 0;
        n++;
      }

    } else {
      del1 = x[0] - c1;
      del2 = x[1] - c2;
      r = Kokkos::sqrt(del1 * del1 + del2 * del2);
      currentradius = radiuslo + (x[2] - lo) * (radiushi - radiuslo) / (hi - lo);

      // z is exterior to cone

      if (r > currentradius || x[2] < lo || x[2] > hi) return 0;

      // z is interior to cone or on its surface
      // surflo = pt on outer circle of bottom end plane, same dir as z vs axis
      // surfhi = pt on outer circle of top end plane, same dir as z vs axis

      if (r > 0.0 && !open_faces[2]) {
        surflo[0] = c1 + del1 * radiuslo / r;
        surflo[1] = c2 + del2 * radiuslo / r;
        surflo[2] = lo;
        surfhi[0] = c1 + del1 * radiushi / r;
        surfhi[1] = c2 + del2 * radiushi / r;
        surfhi[2] = hi;
        point_on_line_segment_kokkos(surflo, surfhi, x, xs);
        delx = x[0] - xs[0];
        dely = x[1] - xs[1];
        delz = x[2] - xs[2];
        dist = Kokkos::sqrt(delx * delx + dely * dely + delz * delz);
        if (dist < cutoff) {
          contact[n].r = dist;
          contact[n].delx = delx;
          contact[n].dely = dely;
          contact[n].delz = delz;
          contact[n].iwall = 2;
          contact[n].radius = -2.0 * (radiuslo + (xs[2] - lo) * (radiushi - radiuslo) / (hi - lo));
          n++;
        }
      }

      delta = x[2] - lo;
      if (delta < cutoff && !open_faces[0]) {
        contact[n].r = delta;
        contact[n].delz = delta;
        contact[n].delx = contact[n].dely = 0.0;
        contact[n].iwall = 0;
        contact[n].radius = 0;
        n++;
      }
      delta = hi - x[2];
      if (delta < cutoff && !open_faces[1]) {
        contact[n].r = delta;
        contact[n].delz = -delta;
        contact[n].delx = contact[n].dely = 0.0;
        contact[n].iwall = 1;
        contact[n].radius = 0;
        n++;
      }
    }

    return n;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int surface_exterior_kokkos(double *x, double cutoff, Contact *contact) const
  {
    double del1, del2, r, currentradius, distsq, distsqprev, crad;
    double corner1[3], corner2[3], corner3[3], corner4[3], xp[3], nearest[3];

    if (axis == 'x') {
      del1 = x[1] - c1;
      del2 = x[2] - c2;
      r = Kokkos::sqrt(del1 * del1 + del2 * del2);
      currentradius = radiuslo + (x[0] - lo) * (radiushi - radiuslo) / (hi - lo);

      // radius of curvature, only used for granular walls

      crad = 0.0;

      // x is far enough from cone that there is no contact
      // x is interior to cone

      if (r >= maxradius + cutoff || x[0] <= lo - cutoff || x[0] >= hi + cutoff) return 0;
      if (r < currentradius && x[0] > lo && x[0] < hi) return 0;

      // x is exterior to cone or on its surface
      // corner1234 = 4 corner pts of half trapezoid = cone surf in plane of x
      // project x to 3 line segments in half trapezoid (4th is axis of cone)
      // nearest = point on surface of cone that x is closest to
      //           could be edge of cone
      // do not add contact point if r >= cutoff

      corner1[0] = lo;
      corner1[1] = c1 + del1 * radiuslo / r;
      corner1[2] = c2 + del2 * radiuslo / r;
      corner2[0] = hi;
      corner2[1] = c1 + del1 * radiushi / r;
      corner2[2] = c2 + del2 * radiushi / r;
      corner3[0] = lo;
      corner3[1] = c1;
      corner3[2] = c2;
      corner4[0] = hi;
      corner4[1] = c1;
      corner4[2] = c2;

      distsq = BIG;

      if (!open_faces[2]) {
        point_on_line_segment_kokkos(corner1, corner2, x, xp);
        distsq = closest_kokkos(x, xp, nearest, distsq);
        crad = -2.0 * (radiuslo + (nearest[0] - lo) * (radiushi - radiuslo) / (hi - lo));
      }

      if (!open_faces[0]) {
        point_on_line_segment_kokkos(corner1, corner3, x, xp);
        distsqprev = distsq;
        distsq = closest_kokkos(x, xp, nearest, distsq);
        if (distsq < distsqprev) crad = 0.0;
      }

      if (!open_faces[1]) {
        point_on_line_segment_kokkos(corner2, corner4, x, xp);
        distsqprev = distsq;
        distsq = closest_kokkos(x, xp, nearest, distsq);
        if (distsq < distsqprev) crad = 0.0;
      }

      if (distsq == BIG) return 0;
      add_contact_kokkos(0, x, nearest[0], nearest[1], nearest[2], contact);
      contact[0].radius = crad;
      contact[0].iwall = 0;
      if (contact[0].r < cutoff) return 1;
      return 0;

    } else if (axis == 'y') {
      del1 = x[0] - c1;
      del2 = x[2] - c2;
      r = Kokkos::sqrt(del1 * del1 + del2 * del2);
      currentradius = radiuslo + (x[1] - lo) * (radiushi - radiuslo) / (hi - lo);

      // radius of curvature, only used for granular walls

      crad = 0.0;

      // y is far enough from cone that there is no contact
      // y is interior to cone

      if (r >= maxradius + cutoff || x[1] <= lo - cutoff || x[1] >= hi + cutoff) return 0;
      if (r < currentradius && x[1] > lo && x[1] < hi) return 0;

      // y is exterior to cone or on its surface
      // corner1234 = 4 corner pts of half trapezoid = cone surf in plane of y
      // project x to 3 line segments in half trapezoid (4th is axis of cone)
      // nearest = point on surface of cone that y is closest to
      //           could be edge of cone
      // do not add contact point if r >= cutoff

      corner1[0] = c1 + del1 * radiuslo / r;
      corner1[1] = lo;
      corner1[2] = c2 + del2 * radiuslo / r;
      corner2[0] = c1 + del1 * radiushi / r;
      corner2[1] = hi;
      corner2[2] = c2 + del2 * radiushi / r;
      corner3[0] = c1;
      corner3[1] = lo;
      corner3[2] = c2;
      corner4[0] = c1;
      corner4[1] = hi;
      corner4[2] = c2;

      distsq = BIG;

      if (!open_faces[2]) {
        point_on_line_segment_kokkos(corner1, corner2, x, xp);
        distsq = closest_kokkos(x, xp, nearest, distsq);
        crad = -2.0 * (radiuslo + (nearest[1] - lo) * (radiushi - radiuslo) / (hi - lo));
      }

      if (!open_faces[0]) {
        point_on_line_segment_kokkos(corner1, corner3, x, xp);
        distsqprev = distsq;
        distsq = closest_kokkos(x, xp, nearest, distsq);
        if (distsq < distsqprev) crad = 0;
      }

      if (!open_faces[1]) {
        point_on_line_segment_kokkos(corner2, corner4, x, xp);
        distsqprev = distsq;
        distsq = closest_kokkos(x, xp, nearest, distsq);
        if (distsq < distsqprev) crad = 0;
      }

      if (distsq == BIG) return 0;
      add_contact_kokkos(0, x, nearest[0], nearest[1], nearest[2], contact);
      contact[0].radius = crad;
      contact[0].iwall = 0;
      if (contact[0].r < cutoff) return 1;
      return 0;

    } else {
      del1 = x[0] - c1;
      del2 = x[1] - c2;
      r = Kokkos::sqrt(del1 * del1 + del2 * del2);
      currentradius = radiuslo + (x[2] - lo) * (radiushi - radiuslo) / (hi - lo);

      // radius of curvature, only used for granular walls

      crad = 0.0;

      // z is far enough from cone that there is no contact
      // z is interior to cone

      if (r >= maxradius + cutoff || x[2] <= lo - cutoff || x[2] >= hi + cutoff) return 0;
      if (r < currentradius && x[2] > lo && x[2] < hi) return 0;

      // z is exterior to cone or on its surface
      // corner1234 = 4 corner pts of half trapezoid = cone surf in plane of z
      // project x to 3 line segments in half trapezoid (4th is axis of cone)
      // nearest = point on surface of cone that z is closest to
      //           could be edge of cone
      // do not add contact point if r >= cutoff

      corner1[0] = c1 + del1 * radiuslo / r;
      corner1[1] = c2 + del2 * radiuslo / r;
      corner1[2] = lo;
      corner2[0] = c1 + del1 * radiushi / r;
      corner2[1] = c2 + del2 * radiushi / r;
      corner2[2] = hi;
      corner3[0] = c1;
      corner3[1] = c2;
      corner3[2] = lo;
      corner4[0] = c1;
      corner4[1] = c2;
      corner4[2] = hi;

      distsq = BIG;

      if (!open_faces[2]) {
        point_on_line_segment_kokkos(corner1, corner2, x, xp);
        distsq = closest_kokkos(x, xp, nearest, distsq);
        crad = -2.0 * (radiuslo + (nearest[2] - lo) * (radiushi - radiuslo) / (hi - lo));
      }

      if (!open_faces[0]) {
        point_on_line_segment_kokkos(corner1, corner3, x, xp);
        distsqprev = distsq;
        distsq = closest_kokkos(x, xp, nearest, distsq);
        if (distsq < distsqprev) crad = 0;
      }

      if (!open_faces[1]) {
        point_on_line_segment_kokkos(corner2, corner4, x, xp);
        distsqprev = distsq;
        distsq = closest_kokkos(x, xp, nearest, distsq);
        if (distsq < distsqprev) crad = 0;
      }

      if (distsq == BIG) return 0;
      add_contact_kokkos(0, x, nearest[0], nearest[1], nearest[2], contact);
      contact[0].radius = crad;
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

  // device version of Region::point_on_line_segment()

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void point_on_line_segment_kokkos(const double *a, const double *b, const double *c,
                                    double *d) const
  {
    const double ba[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const double ca[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    double t = (ca[0] * ba[0] + ca[1] * ba[1] + ca[2] * ba[2]) /
      (ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2]);
    if (t <= 0.0) {
      d[0] = a[0];
      d[1] = a[1];
      d[2] = a[2];
    } else if (t >= 1.0) {
      d[0] = b[0];
      d[1] = b[1];
      d[2] = b[2];
    } else {
      d[0] = a[0] + t * ba[0];
      d[1] = a[1] + t * ba[1];
      d[2] = a[2] + t * ba[2];
    }
  }

  // device version of RegCone::closest()

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  double closest_kokkos(const double *x, const double *near, double *nearest, double dsq) const
  {
    double delx = x[0] - near[0];
    double dely = x[1] - near[1];
    double delz = x[2] - near[2];
    double rsq = delx * delx + dely * dely + delz * delz;
    if (rsq >= dsq) return dsq;

    nearest[0] = near[0];
    nearest[1] = near[1];
    nearest[2] = near[2];
    return rsq;
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
