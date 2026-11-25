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
RegionStyle(plane/kk,RegPlaneKokkos<LMPDeviceType>);
RegionStyle(plane/kk/device,RegPlaneKokkos<LMPDeviceType>);
RegionStyle(plane/kk/host,RegPlaneKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_REGION_PLANE_KOKKOS_H
#define LMP_REGION_PLANE_KOKKOS_H

#include "region_plane.h"

#include "kokkos_base.h"
#include "kokkos_type.h"
#include "region_remap_kokkos.h"

namespace LAMMPS_NS {

struct TagRegPlaneMatchAll{};

template<class DeviceType>
class RegPlaneKokkos : public RegPlane, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  // maximum number of contacts a particle can have with this region, i.e. the
  // length of the contact list the caller has to provide.  matches RegPlane::cmax

  static constexpr int MAXCONTACT = 1;

  RegPlaneKokkos(class LAMMPS *, int, char **);
  // checked conversion from a base class pointer, see the .cpp file
  static RegPlaneKokkos<DeviceType> *cast(class Region *);

  void match_all_kokkos(int, DAT::tdual_int_1d) override;

  void prematch() override
  {
    Region::prematch();
    boxremap.capture(domain);
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRegPlaneMatchAll, const int&) const;

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
    // like RegPlane, no wall velocity is added for a changing plane
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
    double dot = (x - xp) * normal[0] + (y - yp) * normal[1] + (z - zp) * normal[2];

    if (dot >= 0.0) return 1;
    return 0;
  }

  // device versions of RegPlane::surface_interior() and surface_exterior()
  // one contact if 0 <= x < cutoff from the normal (interior) or
  //   non-normal (exterior) side of the plane
  // delxyz = vector from nearest projected point on plane to x

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int surface_interior_kokkos(double *x, double cutoff, Contact *contact) const
  {
    double dot = (x[0] - xp) * normal[0] + (x[1] - yp) * normal[1] + (x[2] - zp) * normal[2];
    if (dot < cutoff && dot >= 0.0) {
      contact[0].r = dot;
      contact[0].delx = dot * normal[0];
      contact[0].dely = dot * normal[1];
      contact[0].delz = dot * normal[2];
      contact[0].radius = 0;
      contact[0].iwall = 0;
      return 1;
    }
    return 0;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int surface_exterior_kokkos(double *x, double cutoff, Contact *contact) const
  {
    double dot = (x[0] - xp) * normal[0] + (x[1] - yp) * normal[1] + (x[2] - zp) * normal[2];
    dot = -dot;
    if (dot < cutoff && dot >= 0.0) {
      contact[0].r = dot;
      contact[0].delx = -dot * normal[0];
      contact[0].dely = -dot * normal[1];
      contact[0].delz = -dot * normal[2];
      contact[0].radius = 0;
      contact[0].iwall = 0;
      return 1;
    }
    return 0;
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
