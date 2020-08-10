// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>
#include <cmath>

DEAL_II_NAMESPACE_OPEN

// ============================================================
// PolarManifold
// ============================================================

template <int dim, int spacedim>
PolarManifold<dim, spacedim>::PolarManifold(const Point<spacedim> center)
  : ChartManifold<dim, spacedim, spacedim>(
      PolarManifold<dim, spacedim>::get_periodicity())
  , center(center)
{}


template <int dim, int spacedim>
Point<spacedim>
PolarManifold<dim, spacedim>::get_periodicity()
{
  Point<spacedim> periodicity;
  // In two dimensions, theta is periodic.
  // In three dimensions things are a little more complicated, since the only
  // variable that is truly periodic is phi, while theta should be bounded
  // between 0 and pi. There is currently no way to enforce this, so here we
  // only fix periodicity for the last variable, corresponding to theta in 2d
  // and phi in 3d.
  periodicity[spacedim - 1] = 2 * numbers::PI;
  return periodicity;
}



template <int dim, int spacedim>
Point<spacedim>
PolarManifold<dim, spacedim>::push_forward(
  const Point<spacedim> &spherical_point) const
{
  Assert(spherical_point[0] >= 0.0,
         ExcMessage("Negative radius for given point."));
  const double rho   = spherical_point[0];
  const double theta = spherical_point[1];

  Point<spacedim> p;
  if (rho > 1e-10)
    switch (spacedim)
      {
        case 2:
          p[0] = rho * std::cos(theta);
          p[1] = rho * std::sin(theta);
          break;
        case 3:
          {
            const double phi = spherical_point[2];
            p[0]             = rho * std::sin(theta) * std::cos(phi);
            p[1]             = rho * std::sin(theta) * std::sin(phi);
            p[2]             = rho * std::cos(theta);
            break;
          }
        default:
          Assert(false, ExcNotImplemented());
      }
  return p + center;
}



template <int dim, int spacedim>
Point<spacedim>
PolarManifold<dim, spacedim>::pull_back(
  const Point<spacedim> &space_point) const
{
  const Tensor<1, spacedim> R   = space_point - center;
  const double              rho = R.norm();

  Point<spacedim> p;
  p[0] = rho;

  switch (spacedim)
    {
      case 2:
        {
          p[1] = std::atan2(R[1], R[0]);
          if (p[1] < 0)
            p[1] += 2 * numbers::PI;
          break;
        }

      case 3:
        {
          const double z = R[2];
          p[2]           = std::atan2(R[1], R[0]); // phi
          if (p[2] < 0)
            p[2] += 2 * numbers::PI; // phi is periodic
          p[1] = std::atan2(std::sqrt(R[0] * R[0] + R[1] * R[1]), z); // theta
          break;
        }

      default:
        Assert(false, ExcNotImplemented());
    }
  return p;
}



template <int dim, int spacedim>
DerivativeForm<1, spacedim, spacedim>
PolarManifold<dim, spacedim>::push_forward_gradient(
  const Point<spacedim> &spherical_point) const
{
  Assert(spherical_point[0] >= 0.0,
         ExcMessage("Negative radius for given point."));
  const double rho   = spherical_point[0];
  const double theta = spherical_point[1];

  DerivativeForm<1, spacedim, spacedim> DX;
  if (rho > 1e-10)
    switch (spacedim)
      {
        case 2:
          {
            DX[0][0] = std::cos(theta);
            DX[0][1] = -rho * std::sin(theta);
            DX[1][0] = std::sin(theta);
            DX[1][1] = rho * std::cos(theta);
            break;
          }

        case 3:
          {
            const double phi = spherical_point[2];
            DX[0][0]         = std::sin(theta) * std::cos(phi);
            DX[0][1]         = rho * std::cos(theta) * std::cos(phi);
            DX[0][2]         = -rho * std::sin(theta) * std::sin(phi);

            DX[1][0] = std::sin(theta) * std::sin(phi);
            DX[1][1] = rho * std::cos(theta) * std::sin(phi);
            DX[1][2] = rho * std::sin(theta) * std::cos(phi);

            DX[2][0] = std::cos(theta);
            DX[2][1] = -rho * std::sin(theta);
            DX[2][2] = 0;
            break;
          }

        default:
          Assert(false, ExcNotImplemented());
      }
  return DX;
}



namespace
{
  template <int dim, int spacedim>
  bool
  spherical_face_is_horizontal(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    const Point<spacedim> &                                     manifold_center)
  {
    // We test whether a face is horizontal by checking that the vertices
    // all have roughly the same distance from the center: If the
    // maximum deviation for the distances from the vertices to the
    // center is less than 1.e-5 of the distance between vertices (as
    // measured by the minimum distance from any of the other vertices
    // to the first vertex), then we call this a horizontal face.
    constexpr unsigned int n_vertices =
      GeometryInfo<spacedim>::vertices_per_face;
    std::array<double, n_vertices>     sqr_distances_to_center;
    std::array<double, n_vertices - 1> sqr_distances_to_first_vertex;
    sqr_distances_to_center[0] =
      (face->vertex(0) - manifold_center).norm_square();
    for (unsigned int i = 1; i < n_vertices; ++i)
      {
        sqr_distances_to_center[i] =
          (face->vertex(i) - manifold_center).norm_square();
        sqr_distances_to_first_vertex[i - 1] =
          (face->vertex(i) - face->vertex(0)).norm_square();
      }
    const auto minmax_sqr_distance =
      std::minmax_element(sqr_distances_to_center.begin(),
                          sqr_distances_to_center.end());
    const auto min_sqr_distance_to_first_vertex =
      std::min_element(sqr_distances_to_first_vertex.begin(),
                       sqr_distances_to_first_vertex.end());

    return (*minmax_sqr_distance.second - *minmax_sqr_distance.first <
            1.e-10 * *min_sqr_distance_to_first_vertex);
  }
} // namespace


/*
template <int dim, int spacedim>
Tensor<1, spacedim>
PolarManifold<dim, spacedim>::normal_vector(
  const typename Triangulation<dim, spacedim>::face_iterator &face,
  const Point<spacedim> &                                     p) const
{
  // Let us first test whether we are on a "horizontal" face
  // (tangential to the sphere).  In this case, the normal vector is
  // easy to compute since it is proportional to the vector from the
  // center to the point 'p'.
  if (spherical_face_is_horizontal<dim, spacedim>(face, center))
    {
      // So, if this is a "horizontal" face, then just compute the normal
      // vector as the one from the center to the point 'p', adequately
      // scaled.
      const Tensor<1, spacedim> unnormalized_spherical_normal = p - center;
      const Tensor<1, spacedim> normalized_spherical_normal =
        unnormalized_spherical_normal / unnormalized_spherical_normal.norm();
      return normalized_spherical_normal;
    }
  else
    // If it is not a horizontal face, just use the machinery of the
    // base class.
    return Manifold<dim, spacedim>::normal_vector(face, p);

  return Tensor<1, spacedim>();
}
*/

// ============================================================
// SphericalManifold
// ============================================================
template <int dim, int spacedim>
SphericalManifold<dim,spacedim>::SphericalManifold(const Point<spacedim> center):
  ChartManifold<dim,spacedim,spacedim>(SphericalManifold<dim,spacedim>::get_periodicity()),
  center(center)
{}



template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim,spacedim>::get_periodicity()
{
  Point<spacedim> periodicity;
  periodicity[spacedim-1] = 2*numbers::PI; // theta and phi period.
  return periodicity;
}


template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim,spacedim>::get_new_point(const Quadrature<spacedim> &quad) const
{
  if (spacedim == 2)
    return ChartManifold<dim,spacedim,spacedim>::get_new_point(quad);
  else
    {
      double rho_average = 0;
      Point<spacedim> mid_point;
      for (unsigned int i=0; i<quad.size(); ++i)
        {
          rho_average += quad.weight(i)*(quad.point(i)-center).norm();
          mid_point += quad.weight(i)*quad.point(i);
        }
      // Project the mid_point back to the right location
      Tensor<1,spacedim> R = mid_point-center;
      // Scale it to have radius rho_average
      R *= rho_average/R.norm();
      // And return it.
      return center+R;
    }
}



template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim,spacedim>::push_forward(const Point<spacedim> &spherical_point) const
{
  Assert(spherical_point[0] >=0.0,
         ExcMessage("Negative radius for given point."));
  const double rho = spherical_point[0];
  const double theta = spherical_point[1];

  Point<spacedim> p;
  if (rho > 1e-10)
    switch (spacedim)
      {
      case 2:
        p[0] = rho*cos(theta);
        p[1] = rho*sin(theta);
        break;
      case 3:
      {
        const double &phi= spherical_point[2];
        p[0] = rho*sin(theta)*cos(phi);
        p[1] = rho*sin(theta)*sin(phi);
        p[2] = rho*cos(theta);
      }
      break;
      default:
        Assert(false, ExcInternalError());
      }
  return p+center;
}

template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim,spacedim>::pull_back(const Point<spacedim> &space_point) const
{
  const Tensor<1,spacedim> R = space_point-center;
  const double rho = R.norm();

  Point<spacedim> p;
  p[0] = rho;

  switch (spacedim)
    {
    case 2:
      p[1] = atan2(R[1],R[0]);
      if (p[1] < 0)
        p[1] += 2*numbers::PI;
      break;
    case 3:
    {
      const double z = R[2];
      p[2] = atan2(R[1],R[0]); // phi
      if (p[2] < 0)
        p[2] += 2*numbers::PI; // phi is periodic
      p[1] = atan2(sqrt(R[0]*R[0]+R[1]*R[1]),z);  // theta
    }
    break;
    default:
      Assert(false, ExcInternalError());
    }
  return p;
}


// ============================================================
// CylindricalManifold
// ============================================================

template <int dim, int spacedim>
CylindricalManifold<dim,spacedim>::CylindricalManifold(const unsigned int axis,
                                                       const double tolerance) :
  direction (Point<spacedim>::unit_vector(axis)),
  point_on_axis (Point<spacedim>()),
  tolerance(tolerance)
{
  Assert(spacedim > 1, ExcImpossibleInDim(1));
}


template <int dim, int spacedim>
CylindricalManifold<dim,spacedim>::CylindricalManifold(const Point<spacedim> &direction,
                                                       const Point<spacedim> &point_on_axis,
                                                       const double tolerance) :
  direction (direction),
  point_on_axis (point_on_axis),
  tolerance(tolerance)
{
  Assert(spacedim > 2, ExcImpossibleInDim(spacedim));
}




template <int dim, int spacedim>
Point<spacedim>
CylindricalManifold<dim,spacedim>::
get_new_point (const Quadrature<spacedim> &quad) const
{
  const std::vector<Point<spacedim> > &surrounding_points = quad.get_points();
  const std::vector<double> &weights = quad.get_weights();

  // compute a proposed new point
  Point<spacedim> middle = flat_manifold.get_new_point(quad);

  double radius = 0;
  Tensor<1,spacedim> on_plane;

  for (unsigned int i=0; i<surrounding_points.size(); ++i)
    {
      on_plane = surrounding_points[i]-point_on_axis;
      on_plane = on_plane - (on_plane*direction) * direction;
      radius += weights[i]*on_plane.norm();
    }

  // we then have to project this point out to the given radius from
  // the axis. to this end, we have to take into account the offset
  // point_on_axis and the direction of the axis
  const Tensor<1,spacedim> vector_from_axis = (middle-point_on_axis) -
                                              ((middle-point_on_axis) * direction) * direction;

  // scale it to the desired length and put everything back together,
  // unless we have a point on the axis
  if (vector_from_axis.norm() <= tolerance * middle.norm())
    return middle;

  else
    return Point<spacedim>((vector_from_axis / vector_from_axis.norm() * radius +
                            ((middle-point_on_axis) * direction) * direction +
                            point_on_axis));
}


// ============================================================
// FunctionChartManifold
// ============================================================
template <int dim, int spacedim, int chartdim>
FunctionManifold<dim,spacedim,chartdim>::FunctionManifold
(const Function<chartdim> &push_forward_function,
 const Function<spacedim> &pull_back_function,
 const Point<chartdim> periodicity,
 const double tolerance):
  ChartManifold<dim,spacedim,chartdim>(periodicity),
  push_forward_function(&push_forward_function),
  pull_back_function(&pull_back_function),
  tolerance(tolerance),
  owns_pointers(false)
{
  AssertDimension(push_forward_function.n_components, spacedim);
  AssertDimension(pull_back_function.n_components, chartdim);
}

template <int dim, int spacedim, int chartdim>
FunctionManifold<dim,spacedim,chartdim>::FunctionManifold
(const std::string push_forward_expression,
 const std::string pull_back_expression,
 const Point<chartdim> periodicity,
 const typename FunctionParser<spacedim>::ConstMap const_map,
 const std::string chart_vars,
 const std::string space_vars,
 const double tolerance) :
  ChartManifold<dim,spacedim,chartdim>(periodicity),
  const_map(const_map),
  tolerance(tolerance),
  owns_pointers(true)
{
  FunctionParser<chartdim> *pf = new FunctionParser<chartdim>(spacedim);
  FunctionParser<spacedim> *pb = new FunctionParser<spacedim>(chartdim);
  pf->initialize(chart_vars, push_forward_expression, const_map);
  pb->initialize(space_vars, pull_back_expression, const_map);
  push_forward_function = pf;
  pull_back_function = pb;
}

template <int dim, int spacedim, int chartdim>
FunctionManifold<dim,spacedim,chartdim>::~FunctionManifold()
{
  if (owns_pointers == true)
    {
      const Function<chartdim> *pf = push_forward_function;
      push_forward_function = 0;
      delete pf;

      const Function<spacedim> *pb = pull_back_function;
      pull_back_function = 0;
      delete pb;
    }
}

template <int dim, int spacedim, int chartdim>
Point<spacedim>
FunctionManifold<dim,spacedim,chartdim>::push_forward(const Point<chartdim> &chart_point) const
{
  Vector<double> pf(spacedim);
  Point<spacedim> result;
  push_forward_function->vector_value(chart_point, pf);
  for (unsigned int i=0; i<spacedim; ++i)
    result[i] = pf[i];

#ifdef DEBUG
  Vector<double> pb(chartdim);
  pull_back_function->vector_value(result, pb);
  for (unsigned int i=0; i<chartdim; ++i)
    Assert((chart_point.norm() > tolerance &&
            (std::abs(pb[i]-chart_point[i]) < tolerance*chart_point.norm())) ||
           (std::abs(pb[i]-chart_point[i]) < tolerance),
           ExcMessage("The push forward is not the inverse of the pull back! Bailing out."));
#endif

  return result;
}


template <int dim, int spacedim, int chartdim>
Point<chartdim>
FunctionManifold<dim,spacedim,chartdim>::pull_back(const Point<spacedim> &space_point) const
{
  Vector<double> pb(chartdim);
  Point<chartdim> result;
  pull_back_function->vector_value(space_point, pb);
  for (unsigned int i=0; i<chartdim; ++i)
    result[i] = pb[i];
  return result;
}

// explicit instantiations
#include "manifold_lib.inst"

DEAL_II_NAMESPACE_CLOSE
