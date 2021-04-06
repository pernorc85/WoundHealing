/* Authors: Jean-Paul Pelteret, University of Cape Town,            */
/*          Andrew McBride, University of Erlangen-Nuremberg, 2010  */
/*                                                                  */
/*    Copyright (C) 2010, 2011, 2012 by the deal.II authors         */
/*                        & Jean-Paul Pelteret and Andrew McBride   */
/*                                                                  */
/*    This file is subject to QPL and may not be  distributed       */
/*    without copyright and license information. Please refer       */
/*    to the file deal.II/doc/license.html for the  text  and       */
/*    further information on this license.                          */

 
#include "step-44-wound-mechanics.h"

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition_selector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <fstream>

#include "nr.h"
#include "util.h"
#include "BMP.h"
#include "make_grid.hh"

extern const int xstep;
extern const int ystep;

extern Mat_DP *fibroblast_density_ptr;
extern Mat_DP *speedfield_ptr;
extern Mat_DP *thetafield_ptr;

extern DP current_time;

// We then stick everything that relates to this tutorial program into a
// namespace of its own, and import all the deal.II function and class names
// into it:
namespace Step44
{
  using namespace dealii;

// @sect3{Run-time parameters}
//
// There are several parameters that can be set in the code so we set up a
// ParameterHandler object to read in the choices at run-time.

//{collagen fiber direction}
template <int dim>
Direction<dim>::Direction ()
                :
                Function<dim> (dim)
{}

template <int dim>
Direction<dim>::Direction (const Mat_DP& collagen_direction_)
        :collagen_direction(collagen_direction_),
         Function<dim> (dim)
{}

template <int dim>
void Direction<dim>::Set (const Mat_DP& collagen_direction_) {
    collagen_direction = collagen_direction_;
}


template <int dim>
inline
void Direction<dim>::vector_value (const Point<dim> &p,
                                       Vector<double>   &values) const
{
  Assert (values.size() == dim,
          ExcDimensionMismatch (values.size(), dim));
  Assert (dim >= 2, ExcNotImplemented());

  double collagenv = collagen_direction[(int)(p(1)*0.1)][(int)(p(0)*0.1)];

  values(0) = cos(collagenv);
  values(1) = sin(collagenv);
}


template <int dim>
void Direction<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                            std::vector<Vector<double> >   &value_list) const
{
  Assert (value_list.size() == points.size(),
          ExcDimensionMismatch (value_list.size(), points.size()));

  const unsigned int n_points = points.size();
  for (unsigned int p=0; p<n_points; ++p){
    Direction<dim>::vector_value (points[p],
                                      value_list[p]);
  } 
}

//========================================
template <int dim>
Magnitude<dim>::Magnitude ()
                :
                Function<dim> (dim)
{}

template <int dim>
Magnitude<dim>::Magnitude (const Mat_DP& collagen_density_, double scale_)
        :collagen_density(collagen_density_),
         scale(scale_),
         Function<dim> (dim)
{}

template <int dim>
void Magnitude<dim>::Set (const Mat_DP& collagen_density_, double scale_) {
    collagen_density = collagen_density_;
    scale = scale_;
}

template <int dim>
inline
void Magnitude<dim>::vector_value (const Point<dim> &p,
                                       double   &value) const
{
    value = collagen_density[(int)(p(1)*scale)][(int)(p(0)*scale)];
}


template <int dim>
void Magnitude<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                            std::vector<double>   &value_list) const
{
    Assert (value_list.size() == points.size(),
          ExcDimensionMismatch (value_list.size(), points.size()));

    const unsigned int n_points = points.size();

    for (unsigned int p=0; p<n_points; ++p)
        Magnitude<dim>::vector_value (points[p],
                                      value_list[p]);
}

//=======================================================================
template <int dim>
RightHandSide<dim>::RightHandSide ()
		:
		Function<dim> (dim)
{}


template <int dim>
inline
void RightHandSide<dim>::vector_value (const Point<dim> &p,
				       Vector<double>   &values) const {
  Assert (values.size() == dim, 
	  ExcDimensionMismatch (values.size(), dim));
  Assert (dim >= 2, ExcNotImplemented());

  Mat_DP &speedfield = *speedfield_ptr;
  Mat_DP &thetafield = *thetafield_ptr;    
  double dx=1,dy=1;
  double value0,value1;


  double theta, speed;
  if(p(0) == 0)
      value0 = (fb_density(p(1),p(0)+dx)-fb_density(p(1),p(0)))/dx;
  else if(p(0) == xstep)
      value0 = (fb_density(p(1),p(0))-fb_density(p(1),p(0)-dx))/dx;
  else  value0 = (fb_density(p(1),p(0)+dx)-fb_density(p(1),p(0)-dx))/2/dx;
  
  if(p(1) == 0)    
      value1 = (fb_density(p(1)+dy,p(0))-fb_density(p(1),p(0)))/dy;
  else if(p(1) == ystep)
      value1 = (fb_density(p(1),p(0))-fb_density(p(1)-dy,p(0)))/dy;
  else value1 = (fb_density(p(1)+dy,p(0))-fb_density(p(1)-dy,p(0)))/2/dy;
  
  speed = speedfield[(int)p(1)][(int)p(0)]; 
  theta = thetafield[(int)p(1)][(int)p(0)]; 
  value0 = value0*50 + 3*fb_density(p(1),p(0))*speed/15*cos(theta);
  value1 = value1*50 + 3*fb_density(p(1),p(0))*speed/15*sin(theta);

  int type = 2;
  double s = 1.0;//default = 1.0
  double dist;
  switch(type){
      case 0:
          //skin tension field in last paper
          if(sqrt(pow(p(0)-xstep*0.2, 2.0)+pow(p(1)-ystep*0.2, 2.0))<500 && p(0)<xstep*0.5 && p(1)<ystep*0.5){
              dist = sqrt(pow(p(0)-xstep*0.2, 2.0)+pow(p(1)-ystep*0.2, 2.0));
              values(0) = -(p(0)-xstep*0.2)/dist;
              values(1) = -(p(1)-ystep*0.2)/dist;
          } else if(sqrt(pow(p(0)-xstep*0.2, 2.0) + pow(p(1)-ystep*0.8, 2.0))<500 && p(0)<xstep*0.5 && p(1)>ystep*0.5){
              dist = sqrt(pow(p(0)-xstep*0.2, 2.0)+pow(p(1)-ystep*0.8, 2.0));
              values(0) = -(p(0)-xstep*0.2)/dist;
              values(1) = -(p(1)-ystep*0.8)/dist;
          } else if(sqrt(pow(p(0)-xstep*0.8, 2.0) + pow(p(1)-ystep*0.2, 2.0))<500 && p(0)>xstep*0.5 && p(1)<ystep*0.5){
              dist = sqrt(pow(p(0)-xstep*0.8, 2.0)+pow(p(1)-ystep*0.2, 2.0));
              values(0) = -(p(0)-xstep*0.8)/dist;
              values(1) = -(p(1)-ystep*0.2)/dist;
          } else if(sqrt(pow(p(0)-xstep*0.8, 2.0) + pow(p(1)-ystep*0.8, 2.0))<500 && p(0)>xstep*0.5 && p(1)>ystep*0.5){
              dist = sqrt(pow(p(0)-xstep*0.8, 2.0)+pow(p(1)-ystep*0.8, 2.0));
              values(0) = -(p(0)-xstep*0.8)/dist;
              values(1) = -(p(1)-ystep*0.8)/dist;
          } else {
              values(0) = 0.0;
              values(1) = 0.0;
          } 
          values(0) *= 5 * s;
          values(1) *= 5 * s;
          //values(0) += value0;
          //values(1) += value1;
          break;

      case 1:{
          //skin tension field A
          int index = int(current_time / 15.0);
          //double sign = index%2 == 0 ? 1 : (-1);
          double sign = 1.0;
          double radius = 150.0;
         
          Point<dim> p1(xstep*0.2, ystep*0.75);
          Point<dim> p2(xstep*0.4, ystep*0.3);
          Point<dim> p3(xstep*0.6, ystep*0.75);
          Point<dim> p4(xstep*0.8, ystep*0.3);

          /*
          double p1x = 900, p1y = 2100;
          double p2x = 1300, p2y = 1300;
          double p3x = 1700, p3y = 2100;
          double p4x = 2100, p4y = 1300;
          */
          double dist1 = p.distance(p1);
          double dist2 = p.distance(p2);
          double dist3 = p.distance(p3);
          double dist4 = p.distance(p4);
          if (dist1 < radius){
              dist = dist1;
              if (dist < radius * 0.25) {
                  values(0) = 0.0;
                  values(1) = sign * s * 10 * 0.75;
              } else {
                  values(0) = 0.0;
                  values(1) = sign * s * 10 * (1 - dist/radius);
              }
          } else if (dist2 < radius){
              dist = dist2;
              if (dist < radius * 0.25) {
                  values(0) = 0.0;
                  values(1) = -sign * s * 10 * 0.75;
              } else {
                  values(0) = 0.0;
                  values(1) = -sign * s * 10 * (1 - dist/radius);
              }
          } else if (dist3 < radius){
              dist = dist3;
              if (dist < radius * 0.25) {
                  values(0) = 0.0;
                  values(1) = sign * s * 10 * 0.75;
              } else {
                  values(0) = 0.0;
                  values(1) = sign * s * 10 * (1 - dist/radius);
              }
          } else if (dist4 < radius){
              dist = dist4;
              if (dist < radius * 0.25) {
                  values(0) = 0.0;
                  values(1) = -sign * s * 10 * 0.75;
              } else {
                  values(0) = 0.0;
                  values(1) = -sign * s * 10 * (1 - dist/radius);
              }
          } else {
              values(0) = 0.0;
              values(1) = 0.0;
          }
          //values(0) += 5 * value0;
          //values(1) += 5 * value1;
          break;
      }

      case 2:{
          //skin tension field B 
          double sign = 1;//index%2 == 0 ? 1 : (-1);
          double radius = 150.0;
          Point<dim> p1(xstep*0.2, ystep*0.22);
          Point<dim> p2(xstep*0.6, ystep*0.22);
          Point<dim> p3(xstep*0.2, ystep*0.78);
          Point<dim> p4(xstep*0.6, ystep*0.78);

          double dist1 = p.distance(p1);
          double dist2 = p.distance(p2);
          double dist3 = p.distance(p3);
          double dist4 = p.distance(p4);

          if(dist1 < radius){
              dist = dist1;
              if (dist < radius * 0.25) {
                  values(0) = sign * (-s) * 7.5 * (0.75);
                  values(1) = sign * (-s) * 7.5 * (0.75);
              } else {
                  values(0) = sign * (-s) * 7.5 * (1 - dist/radius);
                  values(1) = sign * (-s) * 7.5 * (1 - dist/radius);
              }
          } else if (dist2 < radius){
              dist = dist2;
              if (dist < radius * 0.25) {
                  values(0) = sign * s * 7.5 * (0.75);
                  values(1) = sign* (-s) * 7.5 * (0.75);
              } else {
                  values(0) = sign * s * 7.5 * (1 - dist/radius);
                  values(1) = sign* (-s) * 7.5 * (1 - dist/radius);
              }
          } else if (dist3 < radius){
              dist = dist3;
              if (dist < radius * 0.25) {
                  values(0) = sign * (-s) * 7.5 * (0.75);
                  values(1) = sign * s * 7.5 * (0.75);
              } else {
                  values(0) = sign * (-s) * 7.5 * (1 - dist/radius);
                  values(1) = sign * s * 7.5 * (1 - dist/radius);
              }
          } else if (dist4 < radius){
              dist = dist4;
              if (dist < radius * 0.25) {
                  values(0) = sign * s * 7.5 * (0.75);
                  values(1) = sign * s * 7.5 * (0.75);
              } else {
                  values(0) = sign * s * 7.5 * (1 - dist/radius);
                  values(1) = sign * s * 7.5 * (1 - dist/radius);
              }
          } else {
              values(0) = 0.0;
              values(1) = 0.0;
          }
          //values(0) += value0;
          //values(1) += value1;
          break;
      }
      case 3:{
          values(0) = 0.0;
          values(1) = 0.0;
          double dist_to_center = p.distance(Point<dim>(xstep*0.5, ystep*0.5));
          if(dist_to_center > 700 and dist_to_center <750) {
              double dist_to_mid_ring = dist_to_center - 725;
              double theta = get_theta(p(0) - xstep*0.5, p(1) - ystep*0.5);
              while(theta > M_PI / 3.0) {
                  theta -= M_PI / 3.0;
              }
              if(theta < M_PI / 6.0) {
                  values(0) = s * (p(0) - xstep*0.5)/dist_to_center * (1 - 0.04 * dist_to_mid_ring);
                  values(1) = s * (p(1) - ystep*0.5)/dist_to_center * (1 - 0.04 * dist_to_mid_ring);
              }
          }
          //values(0) += value0;
          //values(1) += value1;
          break;
      }
      case 4:{
          //skin tension field 
          double sign = -1.0;
          double radius = 150.0;
         
          Point<dim> p1(xstep*0.83, ystep*0.75);
          Point<dim> p2(xstep*0.83, ystep*0.25);

          double dist1 = p.distance(p1);
          double dist2 = p.distance(p2);
          if (dist1 < radius){
              dist = dist1;
              if (dist < radius * 0.25) {
                  values(0) = sign * s * 10 * 0.75;
                  values(1) = 0.0;
              } else {
                  values(0) = sign * s * 10 * (1 - dist/radius);
                  values(1) = 0.0;
              }
          } else if (dist2 < radius){
              dist = dist2;
              if (dist < radius * 0.25) {
                  values(0) = sign * s * 10 * 0.75;
                  values(1) = 0.0;
              } else {
                  values(0) = sign * s * 10 * (1 - dist/radius);
                  values(1) = 0.0;
              }
          } 
          break;
      }

      default:
          values(0) = value0;
          values(1) = value1;
  } // end of switch

  //values(0) = d/dx*(10^-4N/cell*0.001cell/um^2*fb_density) = d/dx* fb_density * 100kPa
  //10^-4N/cell * 0.001cell/um^2 = 10^-7N/um^2 = 10^5Pa = 100kPa
  //that's why values are multiplied by 100
  //Elastic modulus of skin has unit kPa
}


template <int dim>
void RightHandSide<dim>::vector_value_list (const std::vector<Point<dim> > &points,
					    std::vector<Vector<double> >   &value_list) const 
{
  Assert (value_list.size() == points.size(), 
	  ExcDimensionMismatch (value_list.size(), points.size()));

  const unsigned int n_points = points.size();

  for (unsigned int p=0; p<n_points; ++p)
    RightHandSide<dim>::vector_value (points[p],
				      value_list[p]);
}


// @sect3{Quasi-static quasi-incompressible finite-strain solid}

// The Solid class is the central class in that it represents the problem at
// hand. It follows the usual scheme in that all it really has is a
// constructor, destructor and a <code>run()</code> function that dispatches
// all the work to private functions of this class:
// @sect3{Implementation of the <code>Solid</code> class}

// @sect4{Public interface}
// We initialise the Solid class using data extracted from the parameter file.
template <int dim>
Solid<dim>::Solid(const std::string &input_file, int xstep, int ystep)
    :
    mXstep(xstep), mYstep(ystep),
    parameters(input_file),
    triangulation(Triangulation<dim>::maximum_smoothing),
    time(parameters.end_time, parameters.delta_t),
    timer(std::cout,
          TimerOutput::summary,
          TimerOutput::wall_times),
    degree(parameters.poly_degree),
    // The Finite Element System is composed of dim continuous displacement
    // DOFs, and discontinuous pressure and dilatation DOFs. In an attempt to
    // satisfy the Babuska-Brezzi or LBB stability conditions (see Hughes
    // (2000)), we setup a $Q_n \times DGPM_{n-1} \times DGPM_{n-1}$
    // system. $Q_2 \times DGPM_1 \times DGPM_1$ elements satisfy this
    // condition, while $Q_1 \times DGPM_0 \times DGPM_0$ elements do
    // not. However, it has been shown that the latter demonstrate good
    // convergence characteristics nonetheless.
    fe(FE_Q<dim>(parameters.poly_degree), dim, // displacement
       FE_DGPMonomial<dim>(parameters.poly_degree - 1), 1, // pressure
       FE_DGPMonomial<dim>(parameters.poly_degree - 1), 1), // dilatation
    dof_handler_ref(triangulation),
    dofs_per_cell (fe.dofs_per_cell),
    u_fe(first_u_component),
    p_fe(p_component),
    J_fe(J_component),
    dofs_per_block(n_blocks),
    qf_cell(parameters.quad_order),
    qf_face(parameters.quad_order),
    n_q_points (qf_cell.size()),
    n_q_points_f (qf_face.size())
{
    determine_component_extractors();
}

// The class destructor simply clears the data held by the DOFHandler
template <int dim>
Solid<dim>::~Solid()
{
    dof_handler_ref.clear();
}


// In solving the quasi-static problem, the time becomes a loading parameter,
// i.e. we increasing the loading linearly with time, making the two concepts
// interchangeable. We choose to increment time linearly using a constant time
// step size.
//
// We start the function with preprocessing, setting the initial dilatation
// values, and then output the initial grid before starting the simulation
//  proper with the first time (and loading)
// increment.
//
// Care must be taken (or at least some thought given) when imposing the
// constraint $\widetilde{J}=1$ on the initial solution field. The constraint
// corresponds to the determinant of the deformation gradient in the undeformed
// configuration, which is the identity tensor.
// We use FE_DGPMonomial bases to interpolate the dilatation field, thus we can't
// simply set the corresponding dof to unity as they correspond to the
// monomial coefficients. Thus we use the VectorTools::project function to do
// the work for us. The VectorTools::project function requires an argument
// indicating the hanging node constraints. We have none in this program
// So we have to create a constraint object. In its original state, constraint
// objects are unsorted, and have to be sorted (using the ConstraintMatrix::close function)
// before they can be used. Have a look at step-21 for more information.
// We only need to enforce the initial condition on the dilatation.
// In order to do this, we make use of a ComponentSelectFunction which acts
// as a mask and sets the J_component of n_components to 1. This is exactly what
// we want. Have a look at its usage in step-20 for more information.
template <int dim>
void Solid<dim>::RunIncrementalNonlinear(Mat_DP &collagen_orientation, Mat_DP &collagen_dispersion, Mat_DP &collagen_density) {
    update_qph_for_material_property(collagen_orientation, collagen_dispersion, collagen_density);
    time.restart();
    time.increment();
    // We then declare the incremental solution update $\varDelta
    // \mathbf{\Xi}:= \{\varDelta \mathbf{u},\varDelta \widetilde{p},
    // \varDelta \widetilde{J} \}$ and start the loop over the time domain.
    //
    // At the beginning, we reset the solution update for this time step...
    BlockVector<double> solution_delta(dofs_per_block);
    //assemble_perturbation();
    bool first_round = true;
    while (time.current() <= time.end()) {
        solution_delta = 0.0;

        // ...solve the current time step and update total solution vector
        // $\mathbf{\Xi}_{\textrm{n}} = \mathbf{\Xi}_{\textrm{n-1}} +
        // \varDelta \mathbf{\Xi}$...
        //update_qph_for_growth();

        solve_nonlinear_timestep(solution_delta);
        solution_n += solution_delta;
        //if (first_round) {
            assemble_perturbation();
            solution_n += solution_perturbation;
            first_round = false;
        //}
        // ...and plot the results before moving on happily to the next time
        // step:
        time.increment();
    }
    output_deformation_profile("healing");
    //output_woundcontour();
    update_qph_restore();//solution is restored, but tissue_displacement are not
}

template <int dim>
void Solid<dim>::RunNonlinear(Mat_DP &collagen_orientation, Mat_DP &collagen_dispersion, Mat_DP &collagen_density){
    BlockVector<double> solution_delta(dofs_per_block);
    solution_delta = 0.0;

        // ...solve the current time step and update total solution vector
        // $\mathbf{\Xi}_{\textrm{n}} = \mathbf{\Xi}_{\textrm{n-1}} +
        // \varDelta \mathbf{\Xi}$...
        //update_qph_for_growth();
    update_qph_for_material_property(collagen_orientation, collagen_dispersion, collagen_density);

    solve_nonlinear_timestep(solution_delta);
    solution_n += solution_delta;
        // ...and plot the results before moving on happily to the next time
        // step:
    output_deformation_profile("healing");
    //output_woundcontour();
    update_qph_restore();
}

template <int dim>
void Solid<dim>::Setup(Mat_DP &collagen_orientation, Mat_DP &collagen_dispersion, 
                       Mat_DP &collagen_density)
{
    make_grid_1();
    //make_grid_3<dim>(triangulation, mXstep, mYstep, vol_current, vol_reference, parameters.global_refinement);
    mCollagenDirection.Set(collagen_orientation);
    mCollagenDispersion.Set(collagen_dispersion, 0.1);
    mCollagenDensity.Set(collagen_density, 0.2);
    system_setup();
    {
      ConstraintMatrix constraints;
      constraints.close();

      const ComponentSelectFunction<dim>
        J_mask (J_component, n_components);

      VectorTools::project (dof_handler_ref,
                            constraints,
                            QGauss<dim>(degree+2),
                            J_mask,
                            solution_n);
    }

    tissue_displacement_x = Mat_DP(mYstep, mXstep);
    tissue_displacement_y = Mat_DP(mYstep, mXstep);
    for (size_t i = 0; i < mYstep; i++) {
        for (size_t j = 0; j < mXstep; j++) {
            tissue_displacement_x[i][j] = 0.0;
            tissue_displacement_y[i][j] = 0.0;
        }
    }
}
// @sect3{Private interface}

// @sect4{Threading-building-blocks structures}

// The first group of private member functions is related to parallization.
// We use the Threading Building Blocks library (TBB) to perform as many
// computationally intensive distributed tasks as possible. In particular, we
// assemble the tangent matrix and right hand side vector, the static
// condensation contributions, and update data stored at the quadrature points
// using TBB. Our main tool for this is the WorkStream class (see the @ref
// threads module for more information).

// Firstly we deal with the tangent matrix assembly structures. The
// PerTaskData object stores local contributions.
  template <int dim>
  struct Solid<dim>::PerTaskData_K
  {
    FullMatrix<double>        cell_matrix;
    std::vector<types::global_dof_index> local_dof_indices;

    PerTaskData_K(const unsigned int dofs_per_cell)
      :
      cell_matrix(dofs_per_cell, dofs_per_cell),
      local_dof_indices(dofs_per_cell)
    {}

    void reset()
    {
      cell_matrix = 0.0;
    }
  };


// On the other hand, the ScratchData object stores the larger objects such as
// the shape-function values array (<code>Nx</code>) and a shape function
// gradient and symmetric gradient vector which we will use during the
// assembly.
  template <int dim>
  struct Solid<dim>::ScratchData_K
  {
    FEValues<dim> fe_values_ref;

    std::vector<std::vector<double> >                   Nx;
    std::vector<std::vector<Tensor<2, dim> > >          Grad_Nx;

    ScratchData_K(const FiniteElement<dim> &fe_cell,
                  const QGauss<dim> &qf_cell,
                  const UpdateFlags uf_cell)
      :
      fe_values_ref(fe_cell, qf_cell, uf_cell),
      Nx(qf_cell.size(),
         std::vector<double>(fe_cell.dofs_per_cell)),
      Grad_Nx(qf_cell.size(),
              std::vector<Tensor<2, dim> >(fe_cell.dofs_per_cell))
    {}

    ScratchData_K(const ScratchData_K &rhs)
      :
      fe_values_ref(rhs.fe_values_ref.get_fe(),
                    rhs.fe_values_ref.get_quadrature(),
                    rhs.fe_values_ref.get_update_flags()),
      Nx(rhs.Nx),
      Grad_Nx(rhs.Grad_Nx)
    {}

    void reset()
    {
      const unsigned int n_q_points = Nx.size();
      const unsigned int n_dofs_per_cell = Nx[0].size();
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          Assert( Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
          Assert( Grad_Nx[q_point].size() == n_dofs_per_cell,
                  ExcInternalError());
          for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
            {
              Nx[q_point][k] = 0.0;
              Grad_Nx[q_point][k] = 0.0;
            }
        }
    }

  };

// Next, the same approach is used for the right-hand side assembly.  The
// PerTaskData object again stores local contributions and the ScratchData
// object the shape function object and precomputed values vector:
  template <int dim>
  struct Solid<dim>::PerTaskData_RHS
  {
    Vector<double>            cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;

    PerTaskData_RHS(const unsigned int dofs_per_cell)
      :
      cell_rhs(dofs_per_cell),
      local_dof_indices(dofs_per_cell)
    {}

    void reset()
    {
      cell_rhs = 0.0;
    }
  };


  template <int dim>
  struct Solid<dim>::ScratchData_RHS
  {
    FEValues<dim>     fe_values_ref;
    FEFaceValues<dim> fe_face_values_ref;

    std::vector<std::vector<double> >                   Nx;
    std::vector<std::vector<Tensor<2, dim> > >          Grad_Nx;

    ScratchData_RHS(const FiniteElement<dim> & fe_cell,
                    const QGauss<dim> & qf_cell, const UpdateFlags uf_cell,
                    const QGauss<dim - 1> & qf_face, const UpdateFlags uf_face)
      :
      fe_values_ref(fe_cell, qf_cell, uf_cell),
      fe_face_values_ref(fe_cell, qf_face, uf_face),
      Nx(qf_cell.size(),
         std::vector<double>(fe_cell.dofs_per_cell)),
      Grad_Nx(qf_cell.size(),
                   std::vector<Tensor<2, dim> >
                   (fe_cell.dofs_per_cell))
    {}

    ScratchData_RHS(const ScratchData_RHS & rhs)
      :
      fe_values_ref(rhs.fe_values_ref.get_fe(),
                    rhs.fe_values_ref.get_quadrature(),
                    rhs.fe_values_ref.get_update_flags()),
      fe_face_values_ref(rhs.fe_face_values_ref.get_fe(),
                    rhs.fe_face_values_ref.get_quadrature(),
                    rhs.fe_face_values_ref.get_update_flags()),
      Nx(rhs.Nx),
      Grad_Nx(rhs.Grad_Nx)
    {}

    void reset()
    {
      const unsigned int n_q_points      = Nx.size();
      const unsigned int n_dofs_per_cell = Nx[0].size();
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          Assert( Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
          Assert( Grad_Nx[q_point].size() == n_dofs_per_cell,
                  ExcInternalError());
          for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
            {
              Nx[q_point][k] = 0.0;
              Grad_Nx[q_point][k] = 0.0;
            }
        }
    }

  };

// Then we define structures to assemble the statically condensed tangent
// matrix. Recall that we wish to solve for a displacement-based formulation.
// We do the condensation at the element level as the $\widetilde{p}$ and
// $\widetilde{J}$ fields are element-wise discontinuous.  As these operations
// are matrix-based, we need to setup a number of matrices to store the local
// contributions from a number of the tangent matrix sub-blocks.  We place
// these in the PerTaskData struct.
//
// We choose not to reset any data in the <code>reset()</code> function as the
// matrix extraction and replacement tools will take care of this
  template <int dim>
  struct Solid<dim>::PerTaskData_SC
  {
    FullMatrix<double>        cell_matrix;
    std::vector<types::global_dof_index> local_dof_indices;

    FullMatrix<double>        k_orig;
    FullMatrix<double>        k_pu;
    FullMatrix<double>        k_pJ;
    FullMatrix<double>        k_JJ;
    FullMatrix<double>        k_pJ_inv;
    FullMatrix<double>        k_bbar;
    FullMatrix<double>        A;
    FullMatrix<double>        B;
    FullMatrix<double>        C;

    PerTaskData_SC(const unsigned int dofs_per_cell,
                   const unsigned int n_u,
                   const unsigned int n_p,
                   const unsigned int n_J)
      :
      cell_matrix(dofs_per_cell, dofs_per_cell),
      local_dof_indices(dofs_per_cell),
      k_orig(dofs_per_cell, dofs_per_cell),
      k_pu(n_p, n_u),
      k_pJ(n_p, n_J),
      k_JJ(n_J, n_J),
      k_pJ_inv(n_p, n_J),
      k_bbar(n_u, n_u),
      A(n_J,n_u),
      B(n_J, n_u),
      C(n_p, n_u)
    {}

    void reset()
    {}
  };


// The ScratchData object for the operations we wish to perform here is empty
// since we need no temporary data, but it still needs to be defined for the
// current implementation of TBB in deal.II.  So we create a dummy struct for
// this purpose.
  template <int dim>
  struct Solid<dim>::ScratchData_SC
  {
    void reset()
    {}
  };


// And finally we define the structures to assist with updating the quadrature
// point information. Similar to the SC assembly process, we do not need the
// PerTaskData object (since there is nothing to store here) but must define
// one nonetheless. Note that this is because for the operation that we have
// here -- updating the data on quadrature points -- the operation is purely
// local: the things we do on every cell get consumed on every cell, without
// any global aggregation operation as is usually the case when using the
// WorkStream class. The fact that we still have to define a per-task data
// structure points to the fact that the WorkStream class may be ill-suited to
// this operation (we could, in principle simply create a new task using
// Threads::new_task for each cell) but there is not much harm done to doing
// it this way anyway.
// Furthermore, should there be different material models associated with a
// quadrature point, requiring varying levels of computational expense, then
// the method used here could be advantageous.
  template <int dim>
  struct Solid<dim>::PerTaskData_UQPH
  {
    void reset()
    {}
  };


// The ScratchData object will be used to store an alias for the solution
// vector so that we don't have to copy this large data structure. We then
// define a number of vectors to extract the solution values and gradients at
// the quadrature points.
  template <int dim>
  struct Solid<dim>::ScratchData_UQPH
  {
    const BlockVector<double>   &solution_total;

    std::vector<Tensor<2, dim> > solution_grads_u_total;
    std::vector<double>          solution_values_p_total;
    std::vector<double>          solution_values_J_total;

    FEValues<dim>                fe_values_ref;

    ScratchData_UQPH(const FiniteElement<dim> &fe_cell,
                     const QGauss<dim> &qf_cell,
                     const UpdateFlags uf_cell,
                     const BlockVector<double> &solution_total)
      :
      solution_total(solution_total),
      solution_grads_u_total(qf_cell.size()),
      solution_values_p_total(qf_cell.size()),
      solution_values_J_total(qf_cell.size()),
      fe_values_ref(fe_cell, qf_cell, uf_cell)
    {}

    ScratchData_UQPH(const ScratchData_UQPH &rhs)
      :
      solution_total(rhs.solution_total),
      solution_grads_u_total(rhs.solution_grads_u_total),
      solution_values_p_total(rhs.solution_values_p_total),
      solution_values_J_total(rhs.solution_values_J_total),
      fe_values_ref(rhs.fe_values_ref.get_fe(),
                    rhs.fe_values_ref.get_quadrature(),
                    rhs.fe_values_ref.get_update_flags())
    {}

    void reset()
    {
      const unsigned int n_q_points = solution_grads_u_total.size();
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          solution_grads_u_total[q] = 0.0;
          solution_values_p_total[q] = 0.0;
          solution_values_J_total[q] = 0.0;
        }
    }
  };


// @sect4{Solid::make_grid}

// On to the first of the private member functions. Here we create the
// triangulation of the domain, for which we choose the scaled cube with each
// face given a boundary ID number.  The grid must be refined at least once
// for the indentation problem.
//
// We then determine the volume of the reference configuration and print it
// for comparison:
template <int dim>
void Solid<dim>::make_grid_1() {        
    Point<dim> hole_center1(0.5*mXstep, 0.5*mYstep);
    const double inner_radius = 400;
    if (false) {
        //GridGenerator::hyper_cube (triangulation, 0, mXstep);
        const double outer_radius = 0.5 * mXstep;
        GridGenerator::hyper_cube_with_cylindrical_hole (triangulation, inner_radius, outer_radius);
        GridTools::shift(Point<dim>(0.5*mXstep, 0.5*mXstep), triangulation); 
        //const std::vector<unsigned int> holes = {1,1};
        //GridGenerator::cheese (triangulation, holes); 
        //GridTools::scale(2000.0/3.0, triangulation);
    } else {
        Point<dim> p0(0,0);
        Point<dim> p1(mXstep, mYstep);
        GridGenerator::hyper_rectangle (triangulation, p0, p1);
    }                           
    
    vol_reference = GridTools::volume(triangulation);
    vol_current = vol_reference;
    std::cout << "Grid:\n\t Reference volume: " << vol_reference << std::endl;
	    
    typename Triangulation<dim>::active_cell_iterator cell =
    triangulation.begin_active(), endc = triangulation.end();

    int boundary_count = 0;
    for (; cell != endc; ++cell)
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face){
          if (cell->face(face)->at_boundary()) {
              boundary_count++;
              Point<dim> face_center = cell->face(face)->center();
              double dist_to_center = face_center.distance(hole_center1);
              
              if (face_center(0) < 0.02 * mXstep) {
                  cell->face(face)->set_boundary_id (0);
              } else if (face_center(1) > 0.98 * mYstep or face_center(1) < 0.02 * mYstep) {
                  cell->face(face)->set_boundary_id (2);
                  cout << "boundary_id = 2" << endl;
              } else if (dist_to_center < inner_radius + 100) {
                  cout << "boundary_id = 4" << endl;
                  cell->face(face)->set_boundary_id (4);
              } else {
                  cell->face(face)->set_boundary_id (1);
              }
          }
      }
    cout << "all face boundary_id set" << std::endl;
    
    cout << "boundary_count = " << boundary_count << endl;

    static const PolarManifold<dim> myManifold(hole_center1); 
    triangulation.set_manifold(11, myManifold); 
    for (auto cell = triangulation.begin_active(); cell != triangulation.end(); ++cell) {
        for(unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face){
            cell->face(face)->set_all_manifold_ids(numbers::flat_manifold_id);
            if (cell->face(face)->at_boundary() and cell->face(face)->boundary_id() == 4) {
                cout << "cell face marked with curved manifold 11" << endl;
                cell->face(face)->set_all_manifold_ids(11);
            }
        }
    }

    triangulation.refine_global(std::max (1U, parameters.global_refinement));
/*
    for (auto cell = triangulation.begin_active(); cell != triangulation.end(); ++cell) {
        if (hole_center1.distance(cell->center()) < inner_radius + 150) {
            cell->set_refine_flag();
        }
    }
    triangulation.execute_coarsening_and_refinement();
*/
}

inline
void create_2d_grid(
    const std::vector<Point<2>> &vertices,
    const std::vector<
        std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>>
        &               vertex_indices,
    Triangulation<2> &coarse_grid)
{
  std::vector<CellData<2>> cells(vertex_indices.size());
  for (unsigned int i = 0; i < cells.size(); ++i)
    {
      for (unsigned int j = 0; j < GeometryInfo<2>::vertices_per_cell; ++j)
        cells[i].vertices[j] = vertex_indices[i][j];
    }
  coarse_grid.create_triangulation(vertices, cells, SubCellData());
}


template <int dim>
void Solid<dim>::make_grid_2() {        
    Point<dim> hole_center1(0.3 * mXstep, 0.5 * mYstep);
    Point<dim> hole_center2(0.7 * mXstep, 0.5 * mYstep);
    
    if (mXstep == mYstep) {
        GridGenerator::hyper_cube (triangulation, 0, mXstep);
        //const std::vector<unsigned int> holes = {1,1};
        //GridGenerator::cheese (triangulation, holes); 
        //GridTools::scale(2000.0/3.0, triangulation);
    } else {
        Point<dim> p0(0,0);
        Point<dim> p1(mXstep, mYstep);
        //GridGenerator::hyper_rectangle (triangulation, p0, p1);
        const double inner_radius = 400;
        std::vector<Point<dim>> vertices;
        //first push vertices on hole1
        double theta;
        for(int i = 0; i < 8; i++) {
            theta = 2 * M_PI * (double)i / 8.0;
            vertices.push_back(hole_center1 + Point<dim>(inner_radius * cos(theta), inner_radius * sin(theta)));
        }
        //then push vertices on hole2
        for(int i = 0; i < 8; i++) {
            theta = 2 * M_PI * (double)i / 8.0 + M_PI;
            vertices.push_back(hole_center2 + Point<dim>(inner_radius * cos(theta), inner_radius * sin(theta)));
        }
        //then push vertices on outer boundary
        int outer_start_index = 16;
        //index 16 - 22 will connect with index 1 - 7
        vertices.push_back(Point<dim>(hole_center1[0] + inner_radius, mYstep));
        vertices.push_back(Point<dim>(hole_center1[0], mYstep));
        vertices.push_back(Point<dim>(0.0, mYstep));
        vertices.push_back(Point<dim>(0.0, 0.5 * mYstep));
        vertices.push_back(Point<dim>(0.0, 0.0));
        vertices.push_back(Point<dim>(hole_center1[0], 0.0));
        vertices.push_back(Point<dim>(hole_center1[0] + inner_radius, 0.0));
        //index 23 - 29 will connect with index 9 - 15 
        vertices.push_back(Point<dim>(hole_center2[0] - inner_radius, 0.0));
        vertices.push_back(Point<dim>(hole_center2[0], 0.0));
        vertices.push_back(Point<dim>(mXstep, 0.0));
        vertices.push_back(Point<dim>(mXstep, 0.5 * mYstep));
        vertices.push_back(Point<dim>(mXstep, mYstep));
        vertices.push_back(Point<dim>(hole_center2[0], mYstep));
        vertices.push_back(Point<dim>(hole_center2[0] - inner_radius, mYstep));
        cout << "vertices number = " << vertices.size() << endl;

        std::vector<std::array<unsigned int, GeometryInfo<2>::vertices_per_cell>>
            cell_vertices = {};
        for(int i = 1; i < 7; i++) {
            int inner_i1 = i;
            int inner_i2 = i+1;
            int outer_i1 = i+15;
            int outer_i2 = i+16;
            cout << "push a cell composed of " << inner_i1 << ", " << outer_i1 << ", " << outer_i2 << ", " << inner_i2 << endl;
            cell_vertices.push_back({inner_i2, inner_i1, outer_i2, outer_i1});   
        }
        for(int i = 9; i < 15; i++) {
            int inner_i1 = i;
            int inner_i2 = i+1;  
            int outer_i1 = i+14;
            int outer_i2 = i+15;
            cell_vertices.push_back({inner_i2, inner_i1, outer_i2, outer_i1});
        }
        cout << "vert 0 (" << vertices[0][0] << ", " << vertices[0][1] << ")" << endl;
        cout << "vert 8 (" << vertices[8][0] << ", " << vertices[8][1] << ")" << endl;
        cout << "vert 7 (" << vertices[7][0] << ", " << vertices[7][1] << ")" << endl;
        cout << "vert 9 (" << vertices[9][0] << ", " << vertices[9][1] << ")" << endl;
        cout << "vert 22 (" << vertices[22][0] << ", " << vertices[22][1] << ")" << endl;
        cout << "vert 23 (" << vertices[23][0] << ", " << vertices[23][1] << ")" << endl;
        cell_vertices.push_back({1, 15, 16, 29});
        cell_vertices.push_back({0, 8, 1, 15});
        cell_vertices.push_back({7, 9, 0, 8});
        cell_vertices.push_back({7, 22, 9, 23});

        create_2d_grid(vertices, cell_vertices, triangulation);
    }                           
    
    vol_reference = GridTools::volume(triangulation);
    vol_current = vol_reference;
    std::cout << "Grid:\n\t Reference volume: " << vol_reference << std::endl;
	    
    typename Triangulation<dim>::active_cell_iterator cell =
    triangulation.begin_active(), endc = triangulation.end();

    int boundary_count = 0;
    for (; cell != endc; ++cell)
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face){
          if (cell->face(face)->at_boundary()) {
              boundary_count++;
              Point<dim> face_center = cell->face(face)->center();
              double dist_to_center1 = face_center.distance(hole_center1);
              double dist_to_center2 = face_center.distance(hole_center2);
              
              if (face_center(1) > 0.98 * mYstep or face_center(1) < 0.02 * mYstep) {
                  cell->face(face)->set_boundary_id (2);
                  cout << "boundary_id = 2" << endl;
              } else if (dist_to_center1 < 400) {
                  cout << "boundary_id = 4" << endl;
                  cell->face(face)->set_boundary_id (4);
              } else if (dist_to_center2 < 400) {
                  cout << "boundary_id = 5" << endl;
                  cell->face(face)->set_boundary_id (5);
              } else {
                  cell->face(face)->set_boundary_id (1);
              }
          }
      }
    cout << "all face boundary_id set" << std::endl;
    
    cout << "boundary_count = " << boundary_count << endl;

    static const PolarManifold<dim> myManifold1(hole_center1); 
    static const PolarManifold<dim> myManifold2(hole_center2); 
    triangulation.set_manifold(11, myManifold1); 
    triangulation.set_manifold(22, myManifold2); 
    for (auto cell = triangulation.begin_active(); cell != triangulation.end(); ++cell) {
        for(unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face){
            cell->face(face)->set_all_manifold_ids(numbers::flat_manifold_id);
            if (cell->face(face)->at_boundary() and cell->face(face)->boundary_id() == 4) {
                cout << "cell face marked with curved manifold 11" << endl;
                cell->face(face)->set_all_manifold_ids(11);
            } else if (cell->face(face)->at_boundary() and cell->face(face)->boundary_id() == 5) {
                cout << "cell face marked with curved manifold 22" << endl;
                cell->face(face)->set_all_manifold_ids(22);
            }
        }
    }


    triangulation.refine_global(std::max (1U, parameters.global_refinement));
  }


// @sect4{Solid::system_setup}

// Next we describe how the FE system is setup.  We first determine the number
// of components per block. Since the displacement is a vector component, the
// first dim components belong to it, while the next two describe scalar
// pressure and dilatation DOFs.
  template <int dim>
  void Solid<dim>::system_setup()
  {
    timer.enter_subsection("Setup system");

    std::vector<unsigned int> block_component(n_components, u_dof); // Displacement
    block_component[p_component] = p_dof; // Pressure
    block_component[J_component] = J_dof; // Dilatation

    // The DOF handler is then initialised and we renumber the grid in an
    // efficient manner. We also record the number of DOFs per block.
    dof_handler_ref.distribute_dofs(fe);
    DoFRenumbering::Cuthill_McKee(dof_handler_ref);
    DoFRenumbering::component_wise(dof_handler_ref, block_component);
    DoFTools::count_dofs_per_block(dof_handler_ref, dofs_per_block,
                                   block_component);

    std::cout << "Triangulation:"
              << "\n\t Number of active cells: " << triangulation.n_active_cells()
              << "\n\t Number of degrees of freedom: " << dof_handler_ref.n_dofs()
              << std::endl;

    // Setup the sparsity pattern and tangent matrix
    tangent_matrix.clear();
    {
      const types::global_dof_index n_dofs_u = dofs_per_block[u_dof];
      const types::global_dof_index n_dofs_p = dofs_per_block[p_dof];
      const types::global_dof_index n_dofs_J = dofs_per_block[J_dof];

      BlockDynamicSparsityPattern dsp(n_blocks, n_blocks);

      dsp.block(u_dof, u_dof).reinit(n_dofs_u, n_dofs_u);
      dsp.block(u_dof, p_dof).reinit(n_dofs_u, n_dofs_p);
      dsp.block(u_dof, J_dof).reinit(n_dofs_u, n_dofs_J);

      dsp.block(p_dof, u_dof).reinit(n_dofs_p, n_dofs_u);
      dsp.block(p_dof, p_dof).reinit(n_dofs_p, n_dofs_p);
      dsp.block(p_dof, J_dof).reinit(n_dofs_p, n_dofs_J);

      dsp.block(J_dof, u_dof).reinit(n_dofs_J, n_dofs_u);
      dsp.block(J_dof, p_dof).reinit(n_dofs_J, n_dofs_p);
      dsp.block(J_dof, J_dof).reinit(n_dofs_J, n_dofs_J);
      dsp.collect_sizes();

      // The global system matrix initially has the following structure
      // @f{align*}
      // \underbrace{\begin{bmatrix}
      //   \mathbf{\mathsf{K}}_{uu}  & \mathbf{\mathsf{K}}_{u\widetilde{p}} & \mathbf{0}
      //   \\ \mathbf{\mathsf{K}}_{\widetilde{p}u} & \mathbf{0} & \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}
      //   \\ \mathbf{0} & \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}} & \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
      // \end{bmatrix}}_{\mathbf{\mathsf{K}}(\mathbf{\Xi}_{\textrm{i}})}
      //      \underbrace{\begin{bmatrix}
      //          d \mathbf{\mathsf{u}}
      //      \\  d \widetilde{\mathbf{\mathsf{p}}}
      //      \\  d \widetilde{\mathbf{\mathsf{J}}}
      //      \end{bmatrix}}_{d \mathbf{\Xi}}
      // =
      // \underbrace{\begin{bmatrix}
      //  \mathbf{\mathsf{F}}_{u}(\mathbf{u}_{\textrm{i}})
      //  \\ \mathbf{\mathsf{F}}_{\widetilde{p}}(\widetilde{p}_{\textrm{i}})
      //  \\ \mathbf{\mathsf{F}}_{\widetilde{J}}(\widetilde{J}_{\textrm{i}})
      //\end{bmatrix}}_{ \mathbf{\mathsf{F}}(\mathbf{\Xi}_{\textrm{i}}) } \, .
      // @f}
      // We optimise the sparsity pattern to reflect this structure
      // and prevent unnecessary data creation for the right-diagonal
      // block components.
      Table<2, DoFTools::Coupling> coupling(n_components, n_components);
      for (unsigned int ii = 0; ii < n_components; ++ii)
        for (unsigned int jj = 0; jj < n_components; ++jj)
          if (((ii < p_component) && (jj == J_component))
              || ((ii == J_component) && (jj < p_component))
              || ((ii == p_component) && (jj == p_component)))
            coupling[ii][jj] = DoFTools::none;
          else
            coupling[ii][jj] = DoFTools::always;
      DoFTools::make_hanging_node_constraints(dof_handler_ref, constraints);
      constraints.close();
      DoFTools::make_sparsity_pattern(dof_handler_ref,
                                      coupling,
                                      dsp,
                                      constraints,
                                      false);
      sparsity_pattern.copy_from(dsp);
    }

    tangent_matrix.reinit(sparsity_pattern);

    // We then set up storage vectors
    system_rhs.reinit(dofs_per_block);
    system_rhs.collect_sizes();

    solution_n.reinit(dofs_per_block);
    solution_n.collect_sizes();

    solution_perturbation.reinit(dofs_per_block);
    solution_perturbation.collect_sizes();
    // ...and finally set up the quadrature
    // point history:
    setup_qph();

    timer.leave_subsection();
  }


// @sect4{Solid::determine_component_extractors}
// Next we compute some information from the FE system that describes which local
// element DOFs are attached to which block component.  This is used later to
// extract sub-blocks from the global matrix.
//
// In essence, all we need is for the FESystem object to indicate to which
// block component a DOF on the reference cell is attached to.  Currently, the
// interpolation fields are setup such that 0 indicates a displacement DOF, 1
// a pressure DOF and 2 a dilatation DOF.
  template <int dim>
  void
  Solid<dim>::determine_component_extractors()
  {
    element_indices_u.clear();
    element_indices_p.clear();
    element_indices_J.clear();

    for (unsigned int k = 0; k < fe.dofs_per_cell; ++k)
      {
        const unsigned int k_group = fe.system_to_base_index(k).first.first;
        if (k_group == u_dof)
          element_indices_u.push_back(k);
        else if (k_group == p_dof)
          element_indices_p.push_back(k);
        else if (k_group == J_dof)
          element_indices_J.push_back(k);
        else
          {
            Assert(k_group <= J_dof, ExcInternalError());
          }
      }
  }

// @sect4{Solid::setup_qph}
// The method used to store quadrature information is already described in
// step-18. Here we implement a similar setup for a SMP machine.
//
// Firstly the actual QPH data objects are created. This must be done only
// once the grid is refined to its finest level.
  template <int dim>
  void Solid<dim>::setup_qph()
  {
    std::cout << "    Setting up quadrature point data..." << std::endl;

    {
      triangulation.clear_user_data();
      {
        std::vector<PointHistory<dim> > tmp;
        tmp.swap(quadrature_point_history);
      }

      quadrature_point_history
      .resize(triangulation.n_active_cells() * n_q_points);

      unsigned int history_index = 0;
      for (typename Triangulation<dim>::active_cell_iterator cell =
             triangulation.begin_active(); cell != triangulation.end();
           ++cell)
        {
          cell->set_user_pointer(&quadrature_point_history[history_index]);
          history_index += n_q_points;
        }

      Assert(history_index == quadrature_point_history.size(),
             ExcInternalError());
    }

    // Next we setup the initial quadrature
    // point data:
    std::vector<Vector<double> >   dir_values (n_q_points, Vector<double>(dim));
    std::vector<double>            disp_values (n_q_points, 0.0);
    std::vector<double>            cde_values (n_q_points, 0.0);

    for (typename Triangulation<dim>::active_cell_iterator cell =
           triangulation.begin_active(); cell != triangulation.end(); ++cell)
      {
        PointHistory<dim>* lqph =
          reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

        FEValues<dim> fe_values (fe, qf_cell,//in step-8 here is "quadrature_formula"
			   update_values   | update_gradients | update_quadrature_points | update_JxW_values);
			   
        fe_values.reinit(cell);
                                              
        std::vector<Point<dim> > points = fe_values.get_quadrature_points();
        
        mCollagenDirection.vector_value_list (points, dir_values); 
        mCollagenDispersion.vector_value_list (points, disp_values); 
        mCollagenDensity.vector_value_list (points, cde_values);
                      
        Tensor<1,dim> M_par;
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point){
            //std::cout<<points[q_point](0)<<" "<<points[q_point](1)<<std::endl;
            M_par[0] = dir_values[q_point][0]; 
            M_par[1] = dir_values[q_point][1]; 
            //std::cout<<M_par[0]<<" "<<M_par[1]<<std::endl;    
            lqph[q_point].setup_lqp(parameters, cde_values[q_point], M_par, disp_values[q_point]);
        }
      }
  }

// @sect4{Solid::update_qph_incremental}
// As the update of QP information occurs frequently and involves a number of
// expensive operations, we define a multithreaded approach to distributing
// the task across a number of CPU cores.
//
// To start this, we first we need to obtain the total solution as it stands
// at this Newton increment and then create the initial copy of the scratch and
// copy data objects:
  template <int dim>
  void Solid<dim>::update_qph_incremental(const BlockVector<double> &solution_delta)
  {
    timer.enter_subsection("Update QPH data");
    std::cout << " UQPH " << std::flush;

    const BlockVector<double> solution_total(get_total_solution(solution_delta));

    const UpdateFlags uf_UQPH(update_values | update_gradients);
    PerTaskData_UQPH per_task_data_UQPH;
    ScratchData_UQPH scratch_data_UQPH(fe, qf_cell, uf_UQPH, solution_total);

    // We then pass them and the one-cell update function to the WorkStream to
    // be processed:                                    
    WorkStream::run(dof_handler_ref.begin_active(),
                    dof_handler_ref.end(),
                    *this,
                    &Solid::update_qph_one_cell,
                    &Solid::copy_local_to_global_UQPH,
                    scratch_data_UQPH,
                    per_task_data_UQPH);

    timer.leave_subsection();
  }

  template <int dim>
  void Solid<dim>::update_qph_restore()
  {
    timer.enter_subsection("Update QPH data restore");
    std::cout << " UQPH restore" << std::flush;

    BlockVector<double> solution_total(dofs_per_block);
    {
      ConstraintMatrix constraints;
      constraints.close();

      const ComponentSelectFunction<dim>
        J_mask (J_component, n_components);

      VectorTools::project (dof_handler_ref,
                            constraints,
                            QGauss<dim>(degree+2),
                            J_mask,
                            solution_total);
    }
    solution_n = solution_total;

    const UpdateFlags uf_UQPH(update_values | update_gradients);
    PerTaskData_UQPH per_task_data_UQPH;
    ScratchData_UQPH scratch_data_UQPH(fe, qf_cell, uf_UQPH, solution_total);

    // We then pass them and the one-cell update function to the WorkStream to
    // be processed:                                    
    WorkStream::run(dof_handler_ref.begin_active(),
                    dof_handler_ref.end(),
                    *this,
                    &Solid::update_qph_one_cell,
                    &Solid::copy_local_to_global_UQPH,
                    scratch_data_UQPH,
                    per_task_data_UQPH);

    timer.leave_subsection();
  }

// Now we describe how we extract data from the solution vector and pass it
// along to each QP storage object for processing.
  template <int dim>
  void
  Solid<dim>::update_qph_one_cell(const typename DoFHandler<dim>::active_cell_iterator & cell,
                                              ScratchData_UQPH & scratch,
                                              PerTaskData_UQPH & data)
  {
    PointHistory<dim>* lqph =
      reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
    Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

    Assert(scratch.solution_grads_u_total.size() == n_q_points,
           ExcInternalError());
    Assert(scratch.solution_values_p_total.size() == n_q_points,
           ExcInternalError());
    Assert(scratch.solution_values_J_total.size() == n_q_points,
           ExcInternalError());

    scratch.reset();

    // We first need to find the values and gradients at quadrature points 
    // inside the current cell and then we update each local QP using the 
    // displacement gradient and total pressure and dilatation solution 
    // values:
    scratch.fe_values_ref.reinit(cell);
    scratch.fe_values_ref[u_fe].get_function_gradients(scratch.solution_total,
                                                       scratch.solution_grads_u_total);
    scratch.fe_values_ref[p_fe].get_function_values(scratch.solution_total,
                                                    scratch.solution_values_p_total);
    scratch.fe_values_ref[J_fe].get_function_values(scratch.solution_total,
                                                    scratch.solution_values_J_total);

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      lqph[q_point].update_values(scratch.solution_grads_u_total[q_point],
                                  scratch.solution_values_p_total[q_point],
                                  scratch.solution_values_J_total[q_point]);
  }

  template <int dim>
  void
  Solid<dim>::update_qph_for_growth()
  {
    timer.enter_subsection("Update QPH data for growth");
    std::cout << " UQPHgrowth " << std::flush;
    FEValues<dim> fe_values(fe, qf_cell, update_values | update_gradients |
                                         update_quadrature_points);


    for(typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_ref.begin_active();
        cell != dof_handler_ref.end(); cell++){
      PointHistory<dim> *lqph = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

      fe_values.reinit(cell);
      double growth_rate = 1.0;
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
        if (cell->face(face)->at_boundary()
            && cell->face(face)->boundary_id() == 2)
          growth_rate = 1.0 + time.get_delta_t();

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point){
        Tensor<2, dim> F_g = growth_rate * Tensor<2, dim>(StandardTensors<dim>::I);
        lqph[q_point].update_values_for_growth(F_g);
      }
    }
    timer.leave_subsection();
  }

  template <int dim>
  void
  Solid<dim>::update_qph_for_material_property(const Mat_DP& collagen_direction,
                                               const Mat_DP& collagen_dispersion,
                                               const Mat_DP& collagen_density)
  {
    timer.enter_subsection("Update QPH data for collagen density and orientation change");
    std::cout << " UQPHmaterial_property " << std::flush;

    mCollagenDirection.Set(collagen_direction);
    mCollagenDispersion.Set(collagen_dispersion, 0.1);
    mCollagenDensity.Set(collagen_density, 0.2);

    FEValues<dim> fe_values(fe, qf_cell, update_values | update_gradients |
                                         update_quadrature_points);


    std::vector<Vector<double>> dir_values(n_q_points, Vector<double>(dim));
    std::vector<double> disp_values(n_q_points, 0.0);
    std::vector<double> cdensity_values(n_q_points, 0.0);

    for(typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_ref.begin_active();
            cell != dof_handler_ref.end(); cell++){
        PointHistory<dim> *lqph = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

        fe_values.reinit(cell);

        mCollagenDirection.vector_value_list(fe_values.get_quadrature_points(), dir_values);//mCollagenDirection is data member of Solid
        mCollagenDispersion.vector_value_list(fe_values.get_quadrature_points(), disp_values);//mCollagenDispersion is data member of Solid
        mCollagenDensity.vector_value_list(fe_values.get_quadrature_points(), cdensity_values);//mCollagenDensity is data member of Solid

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point){
            Tensor<1, dim> M_par;
            M_par[0] = dir_values[q_point][0];
            M_par[1] = dir_values[q_point][1];
            lqph[q_point].update_values_for_material_property(parameters,
                                                              M_par,
                                                              disp_values[q_point],
                                                              cdensity_values[q_point]);
        }
    }
    timer.leave_subsection();
  }
// @sect4{Solid::solve_nonlinear_timestep}

// The next function is the driver method for the Newton-Raphson scheme. At
// its top we create a new vector to store the current Newton update step,
// reset the error storage objects and print solver header.
  template <int dim>
  void
  Solid<dim>::solve_nonlinear_timestep(BlockVector<double> & solution_delta)
  {
    std::cout << std::endl << "Timestep " << time.get_timestep() << " @ "
              << time.current() << "s" << std::endl;

    BlockVector<double> newton_update(dofs_per_block);

    error_residual.reset();
    error_residual_0.reset();
    error_residual_norm.reset();
    error_update.reset();
    error_update_0.reset();
    error_update_norm.reset();

    print_conv_header();

                                     // We now perform a number of Newton
                                     // iterations to iteratively solve the
                                     // nonlinear problem.  Since the problem
                                     // is fully nonlinear and we are using a
                                     // full Newton method, the data stored in
                                     // the tangent matrix and right-hand side
                                     // vector is not reusable and must be
                                     // cleared at each Newton step.  We then
                                     // initially build the right-hand side
                                     // vector to check for convergence (and
                                     // store this value in the first
                                     // iteration).  The unconstrained DOFs
                                     // of the rhs vector hold the
                                     // out-of-balance forces. The building is
                                     // done before assembling the system
                                     // matrix as the latter is an expensive
                                     // operation and we can potentially avoid
                                     // an extra assembly process by not
                                     // assembling the tangent matrix when
                                     // convergence is attained.
    unsigned int newton_iteration = 0;
    double last_error_u = 0.0, 
           last_error_p = 0.0,
           last_error_J = 0.0;
    for (; newton_iteration < parameters.max_iterations_NR;
         ++newton_iteration)
      {
        std::cout << " " << std::setw(2) << newton_iteration << " " << std::flush;

        tangent_matrix = 0.0;
        system_rhs = 0.0;

        if (newton_iteration > 0) {
            if (time.end() - time.current() < time.get_delta_t()) {
                std::cout << " last timestep, hard to converge, we will accept the result from first newton iteration. " << std::endl; 
                print_conv_footer();
                break;
            }
        }

        assemble_system_rhs();
        get_error_residual(error_residual);

        if (newton_iteration == 0)
          error_residual_0 = error_residual;

                                         // We can now determine the
                                         // normalised residual error and
                                         // check for solution convergence:
        error_residual_norm = error_residual;
        error_residual_norm.normalise(error_residual_0);

        if (newton_iteration > 0 && error_update_norm.u < parameters.tol_u
            )//&& error_residual_norm.u <= parameters.tol_f) 
        {
            std::cout << " CONVERGED! " << std::endl;
            print_conv_footer();
            break;
        }
 
        if (newton_iteration > 1) {
            if ( error_residual_norm.u > 1.5 * last_error_u or 
                 error_residual_norm.p > 1.5 * last_error_p or
                 error_residual_norm.J > 1.5 * last_error_J ) {
                std::cout << "last_error_u = " << last_error_u << ", error_residual_norm.u = " << error_residual_norm.u << endl;
                std::cout << " cannot converge but we will accept the result " << std::endl;
                print_conv_footer();
                break;
            }
        }

        //record last error u
        last_error_u = error_residual_norm.u;
        last_error_p = error_residual_norm.p;
        last_error_J = error_residual_norm.J;
        /*
        if (newton_iteration == parameters.max_iterations_NR-1 
            && error_update_norm.u > 10e-2){
            time.dec_delta_t();
        }*/

                                         // If we have decided that we want to
                                         // continue with the iteration, we
                                         // assemble the tangent, make and
                                         // impose the Dirichlet constraints,
                                         // and do the solve of the linearised
                                         // system:
        assemble_system_tangent();
        make_constraints(newton_iteration);
        constraints.condense(tangent_matrix, system_rhs);

        const std::pair<unsigned int, double>
          lin_solver_output = solve_linear_system(newton_update);

        get_error_update(newton_update, error_update);
        if (newton_iteration == 0)
          error_update_0 = error_update;

                                         // We can now determine the
                                         // normalised Newton update error,
                                         // and perform the actual update of
                                         // the solution increment for the
                                         // current time step, update all
                                         // quadrature point information
                                         // pertaining to this new
                                         // displacement and stress state and
                                         // continue iterating:
        error_update_norm = error_update;
        error_update_norm.normalise(error_update_0);

        solution_delta += newton_update;
        //added 8/16/2020
        //constraints.distribute(solution_delta);
        update_qph_incremental(solution_delta);
        BlockVector<double> solution_temp(solution_n);
        solution_temp += solution_delta;
        //added 8/16/2020
        //constraints.distribute(solution_temp);

        output_results_temp(newton_iteration, solution_temp);


        std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(7)
                  << std::scientific << lin_solver_output.first << "  "
                  << lin_solver_output.second << "  " << error_residual_norm.norm
                  << "  " << error_residual_norm.u << "  "
                  << error_residual_norm.p << "  " << error_residual_norm.J
                  << "  " << error_update_norm.norm << "  " << error_update_norm.u
                  << "  " << error_update_norm.p << "  " << error_update_norm.J
                  << "  " << std::endl;
      }

                                     // At the end, if it turns out that we
                                     // have in fact done more iterations than
                                     // the parameter file allowed, we raise
                                     // an exception that can be caught in the
                                     // main() function. The call
                                     // <code>AssertThrow(condition,
                                     // exc_object)</code> is in essence
                                     // equivalent to <code>if (!cond) throw
                                     // exc_object;</code> but the former form
                                     // fills certain fields in the exception
                                     // object that identify the location
                                     // (filename and line number) where the
                                     // exception was raised to make it
                                     // simpler to identify where the problem
                                     // happened.
    AssertThrow (newton_iteration <= parameters.max_iterations_NR,
                 ExcMessage("No convergence in nonlinear solver!"));
  }


// @sect4{Solid::print_conv_header and Solid::print_conv_footer}

// This program prints out data in a nice table that is updated
// on a per-iteration basis. The next two functions set up the table
// header and footer:
  template <int dim>
  void Solid<dim>::print_conv_header()
  {
    static const unsigned int l_width = 155;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;

    std::cout << "                 SOLVER STEP                  "
              << " |  LIN_IT   LIN_RES    RES_NORM    "
              << " RES_U     RES_P      RES_J     NU_NORM     "
              << " NU_U       NU_P       NU_J " << std::endl;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;
  }



  template <int dim>
  void Solid<dim>::print_conv_footer()
  {
    static const unsigned int l_width = 155;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;

    const std::pair <double,double> error_dil = get_error_dilation();

    std::cout << "Relative errors:" << std::endl
              << "Displacement:\t" << error_update.u / error_update_0.u << std::endl
              << "Force: \t\t" << error_residual.u / error_residual_0.u << std::endl
              << "Dilatation:\t" << error_dil.first << std::endl
              << "v / V_0:\t" << vol_current << " / " << vol_reference
              << " = " << error_dil.second << std::endl;
  }


// @sect4{Solid::get_error_dilation}

// Calculate how well the dilatation $\widetilde{J}$ agrees with $J :=
// \textrm{det}\ \mathbf{F}$ from the $L^2$ error $ \bigl[ \int_{\Omega_0} {[ J
// - \widetilde{J}]}^{2}\textrm{d}V \bigr]^{1/2}$.
// We also return the ratio of the current volume of the
// domain to the reference volume. This is of interest for incompressible
// media where we want to check how well the isochoric constraint has been
// enforced.
  template <int dim>
  std::pair<double, double>
  Solid<dim>::get_error_dilation()
  {
    double dil_L2_error = 0.0;
    vol_current = 0.0;

    FEValues<dim> fe_values_ref(fe, qf_cell, update_JxW_values);

    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      {
        fe_values_ref.reinit(cell);

        PointHistory<dim> *lqph =
          reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

        Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
        Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const double det_F_qp = lqph[q_point].get_det_Fe();
            const double J_tilde_qp = lqph[q_point].get_J_tilde();
            const double the_error_qp_squared = std::pow((det_F_qp - J_tilde_qp),
                                                         2);
            const double JxW = fe_values_ref.JxW(q_point);

            dil_L2_error += the_error_qp_squared * JxW;
            vol_current += det_F_qp * JxW;
          }
        Assert(vol_current > 0, ExcInternalError());
      }

    std::pair<double, double> error_dil;
    error_dil.first = std::sqrt(dil_L2_error);
    error_dil.second = vol_current / vol_reference;

    return error_dil;
  }


// @sect4{Solid::get_error_residual}

// Determine the true residual error for the problem.  That is, determine the
// error in the residual for the unconstrained degrees of freedom.  Note that to
// do so, we need to ignore constrained DOFs by setting the residual in these
// vector components to zero.
  template <int dim>
  void Solid<dim>::get_error_residual(Errors &error_residual)
  {
    BlockVector<double> error_res(dofs_per_block);

    for (unsigned int i = 0; i < dof_handler_ref.n_dofs(); ++i)
      if (!constraints.is_constrained(i))
        error_res(i) = system_rhs(i);

 





    }

    error_residual.norm = error_res.l2_norm();
    error_residual.u = error_res.block(u_dof).l2_norm();
    error_residual.p = error_res.block(p_dof).l2_norm();
    error_residual.J = error_res.block(J_dof).l2_norm();
  }


// @sect4{Solid::get_error_udpate}

// Determine the true Newton update error for the problem
  template <int dim>
  void Solid<dim>::get_error_update(const BlockVector<double> &newton_update,
                                    Errors &error_update)
  {
    BlockVector<double> error_ud(dofs_per_block);
    for (unsigned int i = 0; i < dof_handler_ref.n_dofs(); ++i)
      if (!constraints.is_constrained(i))
        error_ud(i) = newton_update(i);

    error_update.norm = error_ud.l2_norm();
    error_update.u = error_ud.block(u_dof).l2_norm();
    error_update.p = error_ud.block(p_dof).l2_norm();
    error_update.J = error_ud.block(J_dof).l2_norm();
  }



// @sect4{Solid::get_total_solution}

// This function provides the total solution, which is valid at any Newton step.
// This is required as, to reduce computational error, the total solution is
// only updated at the end of the timestep.
  template <int dim>
  BlockVector<double>
  Solid<dim>::get_total_solution(const BlockVector<double> &solution_delta) const
  {
    BlockVector<double> solution_total(solution_n);
    solution_total += solution_delta;
    //added 8/16/2020
    //constraints.distribute(solution_total);
    return solution_total;
  }


// @sect4{Solid::assemble_system_tangent}

// Since we use TBB for assembly, we simply setup a copy of the
// data structures required for the process and pass them, along
// with the memory addresses of the assembly functions to the
// WorkStream object for processing. Note that we must ensure that
// the matrix is reset before any assembly operations can occur.
  template <int dim>
  void Solid<dim>::assemble_system_tangent()
  {
    timer.enter_subsection("Assemble tangent matrix");
    std::cout << " ASM_K " << std::flush;

    tangent_matrix = 0.0;

    const UpdateFlags uf_cell(update_values    |
                              update_gradients |
                              update_JxW_values);

    PerTaskData_K per_task_data(dofs_per_cell);
    ScratchData_K scratch_data(fe, qf_cell, uf_cell);

    WorkStream::run(dof_handler_ref.begin_active(),
                    dof_handler_ref.end(),
                    *this,
                    &Solid::assemble_system_tangent_one_cell,
                    &Solid::copy_local_to_global_K,
                    scratch_data,
                    per_task_data);

    timer.leave_subsection();
  }

// This function adds the local contribution to the system matrix.
// Note that we choose not to use the constraint matrix to do the
// job for us because the tangent matrix and residual processes have
// been split up into two separate functions.
  template <int dim>
  void Solid<dim>::copy_local_to_global_K(const PerTaskData_K &data)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        tangent_matrix.add(data.local_dof_indices[i],
                           data.local_dof_indices[j],
                           data.cell_matrix(i, j));
  }

// Of course, we still have to define how we assemble the tangent matrix
// contribution for a single cell.  We first need to reset and initialise some
// of the scratch data structures and retrieve some basic information
// regarding the DOF numbering on this cell.  We can precalculate the cell
// shape function values and gradients. Note that the shape function gradients
// are defined with regard to the current configuration.  That is
// $\textrm{grad}\ \boldsymbol{\varphi} = \textrm{Grad}\ \boldsymbol{\varphi}
// \ \mathbf{F}^{-1}$.
  template <int dim>
  void
  Solid<dim>::assemble_system_tangent_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                               ScratchData_K &scratch,
                                               PerTaskData_K &data)
  {
    data.reset();
    scratch.reset();
    scratch.fe_values_ref.reinit(cell);
    cell->get_dof_indices(data.local_dof_indices);
    PointHistory<dim> *lqph =
      reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
        for (unsigned int k = 0; k < dofs_per_cell; ++k) {
            const unsigned int k_group = fe.system_to_base_index(k).first.first;

            if (k_group == u_dof)
                scratch.Grad_Nx[q_point][k] = scratch.fe_values_ref[u_fe].gradient(k, q_point);
            else if (k_group == p_dof)
                scratch.Nx[q_point][k] = scratch.fe_values_ref[p_fe].value(k, q_point);
            else if (k_group == J_dof)
                scratch.Nx[q_point][k] = scratch.fe_values_ref[J_fe].value(k, q_point);
            else
              Assert(k_group <= J_dof, ExcInternalError());
        }
    }

    // Now we build the local cell stiffness matrix. Since the global and
    // local system matrices are symmetric, we can exploit this property by
    // building only the lower half of the local matrix and copying the values
    // to the upper half.  So we only assemble half of the
    // $\mathsf{\mathbf{k}}_{uu}$, $\mathsf{\mathbf{k}}_{\widetilde{p}
    // \widetilde{p}} = \mathbf{0}$, $\mathsf{\mathbf{k}}_{\widetilde{J}
    // \widetilde{J}}$ blocks, while the whole
    // $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$,
    // $\mathsf{\mathbf{k}}_{\mathbf{u} \widetilde{J}} = \mathbf{0}$,
    // $\mathsf{\mathbf{k}}_{\mathbf{u} \widetilde{p}}$ blocks are built.
    //
    // In doing so, we first extract some configuration dependent variables
    // from our QPH history objects for the current quadrature point.

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
        const Tensor<2, dim> S           = lqph[q_point].get_S();
        Tensor<4, dim> Cmaterial         = lqph[q_point].get_Cmaterial();
        const double d2Psi_vol_dJ2       = lqph[q_point].get_d2Psi_vol_dJ2();
        const double det_Fe              = lqph[q_point].get_det_Fe();
        const Tensor<2, dim> F           = lqph[q_point].get_F();
        const Tensor<2, dim> F_inv       = lqph[q_point].get_F_inv();
        // Next we define some aliases to make the assembly process easier to
        // follow        
        const std::vector<double>
        &N = scratch.Nx[q_point];
        const std::vector<Tensor<2, dim> >
        &Grad_Nx = scratch.Grad_Nx[q_point];
        const double JxW = scratch.fe_values_ref.JxW(q_point);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const unsigned int component_i = fe.system_to_component_index(i).first;
            const unsigned int i_group     = fe.system_to_base_index(i).first.first;

            for (unsigned int j = 0; j <= i; ++j)
              {
                const unsigned int component_j = fe.system_to_component_index(j).first;
                const unsigned int j_group     = fe.system_to_base_index(j).first.first;

                // This is the $\mathsf{\mathbf{k}}_{\mathbf{u} \mathbf{u}}$
                // contribution. It comprises a material contribution, and a
                // geometrical stress contribution which is only added along
                // the local matrix diagonals:
                if ((i_group == j_group) && (i_group == u_dof))
                  {
                    Tensor<2,dim> Grad_i_x_F = transpose(F) * Grad_Nx[i];
                    Tensor<2,dim> Grad_j_x_F = transpose(F) * Grad_Nx[j];
                    Tensor<2,dim> tmp = double_contract<0,0,1,1>(transpose(Grad_i_x_F), Cmaterial);
                    double tmp2 = double_contract<0,0,1,1>(tmp, Grad_j_x_F);
                    data.cell_matrix(i, j) += tmp2 // The material contribution:
                                               * JxW;
//                    data.cell_matrix(i, j) += symm_grad_Nx[i] * Jc // The material contribution:
//                                              * symm_grad_Nx[j] * JxW;
                    if (component_i == component_j) // geometrical stress contribution
                      data.cell_matrix(i, j) += Grad_Nx[i][component_i] * S 
                                                * Grad_Nx[j][component_j] * JxW;
                  }
                // Next is the $\mathsf{\mathbf{k}}_{ \widetilde{p} \mathbf{u}}$ contribution
                else if ((i_group == p_dof) && (j_group == u_dof))
                  {
                    data.cell_matrix(i, j) += N[i] * det_Fe
                                              * double_contract<0,0,1,1> (Grad_Nx[j] * F_inv,
                                                   Tensor<2, dim>(StandardTensors<dim>::I))
                                              * JxW;
                  }
                // and lastly the $\mathsf{\mathbf{k}}_{ \widetilde{J} \widetilde{p}}$
                // and $\mathsf{\mathbf{k}}_{ \widetilde{J} \widetilde{J}}$
                // contributions:
                else if ((i_group == J_dof) && (j_group == p_dof))
                  data.cell_matrix(i, j) -= N[i] * N[j] * JxW;
                else if ((i_group == j_group) && (i_group == J_dof))
                  data.cell_matrix(i, j) += N[i] * d2Psi_vol_dJ2 * N[j] * JxW;
                else
                  Assert((i_group <= J_dof) && (j_group <= J_dof),
                         ExcInternalError());
              }
          }
      }

    // Finally, we need to copy the lower half of the local matrix into the
    // upper half:
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
        data.cell_matrix(i, j) = data.cell_matrix(j, i);
  }

// @sect4{Solid::assemble_system_rhs}
// The assembly of the right-hand side process is similar to the
// tangent matrix, so we will not describe it in too much detail.
// Note that since we are describing a problem with Neumann BCs,
// we will need the face normals and so must specify this in the
// update flags.
  template <int dim>
  void Solid<dim>::assemble_system_rhs()
  {
    timer.enter_subsection("Assemble system right-hand side");
    std::cout << " ASM_R " << std::flush;

    system_rhs = 0.0;

    const UpdateFlags uf_cell(update_values |
                              update_gradients |
                              update_quadrature_points |
                              update_JxW_values);
    const UpdateFlags uf_face(update_values |
                              update_normal_vectors |
                              update_JxW_values);

    PerTaskData_RHS per_task_data(dofs_per_cell);
    ScratchData_RHS scratch_data(fe, qf_cell, uf_cell, qf_face, uf_face);

    WorkStream::run(dof_handler_ref.begin_active(),
                    dof_handler_ref.end(),
                    *this,
                    &Solid::assemble_system_rhs_one_cell,
                    &Solid::copy_local_to_global_rhs,
                    scratch_data,
                    per_task_data);

    timer.leave_subsection();
  }


  template <int dim>
  void Solid<dim>::copy_local_to_global_rhs(const PerTaskData_RHS &data)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      system_rhs(data.local_dof_indices[i]) += data.cell_rhs(i);
  }



  template <int dim>
  void
  Solid<dim>::assemble_system_rhs_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                           ScratchData_RHS &scratch,
                                           PerTaskData_RHS &data)
  {
    data.reset();
    scratch.reset();
    scratch.fe_values_ref.reinit(cell);
    cell->get_dof_indices(data.local_dof_indices);
    PointHistory<dim> *lqph =
      reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
        for (unsigned int k = 0; k < dofs_per_cell; ++k) {
            const unsigned int k_group = fe.system_to_base_index(k).first.first;

            if (k_group == u_dof) 
              scratch.Grad_Nx[q_point][k]
                = scratch.fe_values_ref[u_fe].gradient(k, q_point);
            else if (k_group == p_dof)
              scratch.Nx[q_point][k] = scratch.fe_values_ref[p_fe].value(k,
                                                                         q_point);
            else if (k_group == J_dof)
              scratch.Nx[q_point][k] = scratch.fe_values_ref[J_fe].value(k,
                                                                         q_point);
            else
              Assert(k_group <= J_dof, ExcInternalError());
        }
    }

    auto points = scratch.fe_values_ref.get_quadrature_points();
    
    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
        const Tensor<2, dim> S = lqph[q_point].get_S();
        const double det_Fe = lqph[q_point].get_det_Fe();
        const Tensor<2, dim> F = lqph[q_point].get_F();
        const double J_tilde = lqph[q_point].get_J_tilde();
        const double p_tilde = lqph[q_point].get_p_tilde();
        const double dPsi_vol_dJ = lqph[q_point].get_dPsi_vol_dJ();

        const std::vector<double>
        &Nx = scratch.Nx[q_point];
        const std::vector<Tensor<2, dim> >
        &Grad_Nx = scratch.Grad_Nx[q_point];
        const double JxW = scratch.fe_values_ref.JxW(q_point);

        // We first compute the contributions
        // from the internal forces.  Note, by
        // definition of the rhs as the negative
        // of the residual, these contributions
        // are subtracted.
        double tethering_coeff = 0.0;//0.004;
        Vector<double> deformation(dim+2);
        //const BlockVector<double> solution_total(get_total_solution(solution_delta));
        VectorTools::point_value(dof_handler_ref, solution_n, points[q_point], deformation);
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            const unsigned int i_group = fe.system_to_base_index(i).first.first;
            const unsigned int component_i = fe.system_to_component_index(i).first;
            if (i_group == u_dof) {
                data.cell_rhs(i) -= double_contract<0,0,1,1>(transpose(F)*Grad_Nx[i], S) * JxW;
                //if (deformation[0] > 50.0 or deformation[1] > 50.0) {
                //    cout << "deformation = " << deformation[0] << ", " << deformation[1] << endl;
                //}
                data.cell_rhs(i) -= tethering_coeff * deformation[component_i] * JxW;
            } else if (i_group == p_dof)
              data.cell_rhs(i) -= Nx[i] * (det_Fe - J_tilde) * JxW;
            else if (i_group == J_dof)
              data.cell_rhs(i) -= Nx[i] * (dPsi_vol_dJ - p_tilde) * JxW;
            else
              Assert(i_group <= J_dof, ExcInternalError());
        }
    }

    RightHandSide<dim> fibroblast_force;
    std::vector<Vector<double>> fibroblast_force_values(n_q_points, Vector<double>(dim));

    //auto points = scratch.fe_values_ref.get_quadrature_points();
    fibroblast_force.vector_value_list(points, fibroblast_force_values);
	//fibroblast_force should data member of Solid
   
    const double    time_ramp = time.current() / time.end(); 
    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
        Tensor<1, dim> traction;
        for(size_t i=0; i<dim; i++)
            traction[i] = time_ramp * fibroblast_force_values[q_point][i];

        //cout << "traction: " << traction[0] << ", " << traction[1] << endl;
		//this is the entry point for fibroblast-produced forces

        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            const unsigned int i_group = fe.system_to_base_index(i).first.first;

            if (i_group == u_dof) {
                const unsigned int component_i =
                    fe.system_to_component_index(i).first;
                const double Ni =
                    scratch.fe_values_ref.shape_value(i, q_point);
                const double JxW = scratch.fe_values_ref.JxW(q_point);

                data.cell_rhs(i) += (Ni * traction[component_i]) * JxW;
            }
        }
    }
  }

  template <int dim>
  void Solid<dim>::assemble_perturbation()
  {
    timer.enter_subsection("Assemble system perturbation");
    std::cout << " ASM_perturbation " << std::flush;

    solution_perturbation = 0.0;

    const UpdateFlags uf_cell(update_values |
                              update_gradients |
                              update_JxW_values);
    const UpdateFlags uf_face(update_values |
                              update_normal_vectors |
                              update_JxW_values);

    PerTaskData_RHS per_task_data(dofs_per_cell);
    ScratchData_RHS scratch_data(fe, qf_cell, uf_cell, qf_face, uf_face);

    WorkStream::run(dof_handler_ref.begin_active(),
                    dof_handler_ref.end(),
                    *this,
                    &Solid::assemble_system_perturb_one_cell,
                    &Solid::copy_local_to_global_perturb,
                    scratch_data,
                    per_task_data);

    timer.leave_subsection();
  }

  template <int dim>
  void Solid<dim>::copy_local_to_global_perturb(const PerTaskData_RHS &data)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      solution_perturbation(data.local_dof_indices[i]) += data.cell_rhs(i);
  }


  template <int dim>
  void
  Solid<dim>::assemble_system_perturb_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                           ScratchData_RHS &scratch,
                                           PerTaskData_RHS &data)
  {
    data.reset();
    scratch.reset();
    cell->get_dof_indices(data.local_dof_indices);


    srand(std::time(NULL));
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        const unsigned int i_group = fe.system_to_base_index(i).first.first;

        if (i_group == u_dof){
            data.cell_rhs(i) = double(rand())/double(RAND_MAX) * 0.001;
        } else if (i_group == p_dof){
            data.cell_rhs(i) = 0.0;
        } else if (i_group == J_dof){
            data.cell_rhs(i) = 0.0;
        } else {
            Assert(i_group <= J_dof, ExcInternalError());
        }
    }
  }




// The constraints for this problem are simple to describe.
// However, since we are dealing with an iterative Newton method,
// it should be noted that any displacement constraints should only
// be specified at the zeroth iteration and subsequently no
// additional contributions are to be made since the constraints
// are already exactly satisfied.
template <int dim>
void Solid<dim>::make_constraints(const int & it_nr)
{
    std::cout << " CST " << std::flush;

    // Since the constraints are different at different Newton iterations, we
    // need to clear the constraints matrix and completely rebuild
    // it. However, after the first iteration, the constraints remain the same
    // and we can simply skip the rebuilding step if we do not clear it.

    if (it_nr > 1)
      return;
    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler_ref, constraints);
    const bool apply_dirichlet_bc = (it_nr == 0);

    // The boundary conditions for the indentation problem are as follows: On
    // the -x, -y and -z faces (ID's 0,2,4) we set up a symmetry condition to
    // allow only planar movement while the +x and +y faces (ID's 1,3) are
    // traction free. In this contrived problem, part of the +z face (ID 5) is
    // set to have no motion in the x- and y-component. Finally, as described
    // earlier, the other part of the +z face has an the applied pressure but
    // is also constrained in the x- and y-directions.
    //
    // In the following, we will have to tell the function interpolation
    // boundary values which components of the solution vector should be
    // constrained (i.e., whether it's the x-, y-, z-displacements or
    // combinations thereof). This is done using ComponentMask objects (see
    // @ref GlossComponentMask) which we can get from the finite element if we
    // provide it with an extractor object for the component we wish to
    // select. To this end we first set up such extractor objects and later
    // use it when generating the relevant component masks:
    const FEValuesExtractors::Scalar x_displacement(0);
    const FEValuesExtractors::Scalar y_displacement(1);

    {
      const int boundary_id = 0;

      std::vector<bool> components(n_components, true);

      if (apply_dirichlet_bc == true)
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ConstantFunction<dim>(0.0, n_components),
                                                 constraints,
                                                 components);
                                                 //fe.component_mask(x_displacement));
      else
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 components);
                                                 //fe.component_mask(x_displacement));
    }
    {
      const int boundary_id = 1;

      std::vector<bool> components(n_components, true);

      if (apply_dirichlet_bc == true)
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 components);
                                                 //fe.component_mask(x_displacement));
      else
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 components);
                                                 //fe.component_mask(x_displacement));
    }
    {
      const int boundary_id = 2;

      std::vector<bool> components(n_components, true);

      if (apply_dirichlet_bc == true)
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 components);
                                                 //fe.component_mask(x_displacement));
      else
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 components);
                                                 //fe.component_mask(x_displacement));
    }
    {
      const int boundary_id = 4;

      std::vector<bool> components(n_components, true);

      if (apply_dirichlet_bc == true) {
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ConstantFunction<dim>(30.0, n_components),
                                                 constraints,
                                                 //components);
                                                 fe.component_mask(y_displacement));

        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ConstantFunction<dim>(0.0, n_components),
                                                 constraints,
                                                 //components);
                                                 fe.component_mask(y_displacement));
     
      } else
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 components);
    }
    {
      const int boundary_id = 5;

      std::vector<bool> components(n_components, true);

      if (apply_dirichlet_bc == true)
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 //components);
                                                 fe.component_mask(y_displacement));
      else
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 //components);
                                                 fe.component_mask(y_displacement));
    }


    bool extra_constraint = true;
    if (!extra_constraint) {
        constraints.close();
        return;
    } 
    std::vector<bool> dof_touched(dof_handler_ref.n_dofs(), false);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_ref.begin_active(),
                                                 endc = dof_handler_ref.end();

    //cout << "dofs_per_cell = " << fe.dofs_per_cell << endl;
    std::vector<unsigned int> local_dof_indices(fe.dofs_per_cell);
    for (; cell != endc; ++cell) {
        cell->get_dof_indices (local_dof_indices);
        Point<dim> cell_center = cell->center();

        for(unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; vertex++) {
            unsigned int dof_index_x = cell->vertex_dof_index(vertex,0);
            unsigned int dof_index_y = cell->vertex_dof_index(vertex,1);
            //cout << "vertex " << vertex << ", dof " << dof_index << endl;
            Point<dim> dof_location = cell->vertex(vertex);
            //cout << "location = " << dof_location(0) << ", " << dof_location(1) << endl;
            if (abs(dof_location(1) - 0.7 * ystep) < 80.0 and
                    abs(dof_location(0) - 0.5 * xstep) < 500.0 ) {
                constraints.add_line(dof_index_x);
                constraints.add_line(dof_index_y);
                constraints.set_inhomogeneity(dof_index_x, 0.0);//
                cout << "dof_index = " << dof_index_x << "constrint added" << endl;
                cout << "dof_index = " << dof_index_y << "constrint added" << endl;
            }
        }
    }
    constraints.close();
}

// @sect4{Solid::solve_linear_system}
// Solving the entire block system is a bit problematic as there are no
// contributions to the $\mathsf{\mathbf{k}}_{ \widetilde{J} \widetilde{J}}$
// block, rendering it non-invertible.
// Since the pressure and dilatation variables DOFs are discontinuous, we can
// condense them out to form a smaller displacement-only system which
// we will then solve and subsequently post-process to retrieve the
// pressure and dilatation solutions.
//
// At the top, we allocate two temporary vectors to help with the static
// condensation, and variables to store the number of linear solver iterations
// and the (hopefully converged) residual.
//
// For the following, recall that
// @f{align*}
//  \mathbf{\mathsf{K}}_{\textrm{store}}
//:=
//  \begin{bmatrix}
//      \mathbf{\mathsf{K}}_{\textrm{con}}      &       \mathbf{\mathsf{K}}_{u\widetilde{p}}    & \mathbf{0}
//  \\  \mathbf{\mathsf{K}}_{\widetilde{p}u}    &       \mathbf{0}      &       \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1}
//  \\  \mathbf{0}      &       \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}                & \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
//  \end{bmatrix} \, .
// @f}
// and
//  @f{align*}
//              d \widetilde{\mathbf{\mathsf{p}}}
//              & = \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \bigl[
//                       \mathbf{\mathsf{F}}_{\widetilde{J}}
//                       - \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}} d \widetilde{\mathbf{\mathsf{J}}} \bigr]
//              \\ d \widetilde{\mathbf{\mathsf{J}}}
//              & = \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1} \bigl[
//                      \mathbf{\mathsf{F}}_{\widetilde{p}}
//                      - \mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
//                      \bigr]
//               \\ \Rightarrow d \widetilde{\mathbf{\mathsf{p}}}
//              &=  \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{F}}_{\widetilde{J}}
//              - \underbrace{\bigl[\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
//              \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1}\bigr]}_{\overline{\mathbf{\mathsf{K}}}}\bigl[ \mathbf{\mathsf{F}}_{\widetilde{p}}
//              - \mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}} \bigr]
//  @f}
//  and thus
//  @f[
//              \underbrace{\bigl[ \mathbf{\mathsf{K}}_{uu} + \overline{\overline{\mathbf{\mathsf{K}}}}~ \bigr]
//              }_{\mathbf{\mathsf{K}}_{\textrm{con}}} d \mathbf{\mathsf{u}}
//              =
//          \underbrace{
//              \Bigl[
//              \mathbf{\mathsf{F}}_{u}
//                      - \mathbf{\mathsf{K}}_{u\widetilde{p}} \bigl[ \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{F}}_{\widetilde{J}}
//                      - \overline{\mathbf{\mathsf{K}}}\mathbf{\mathsf{F}}_{\widetilde{p}} \bigr]
//              \Bigr]}_{\mathbf{\mathsf{F}}_{\textrm{con}}}
//  @f]
//  where
//  @f[
//              \overline{\overline{\mathbf{\mathsf{K}}}} :=
//                      \mathbf{\mathsf{K}}_{u\widetilde{p}} \overline{\mathbf{\mathsf{K}}} \mathbf{\mathsf{K}}_{\widetilde{p}u} \, .
//  @f]
  template <int dim>
  std::pair<unsigned int, double>
  Solid<dim>::solve_linear_system(BlockVector<double> &newton_update)
  {
    BlockVector<double> A(dofs_per_block);
    BlockVector<double> B(dofs_per_block);

    unsigned int lin_it = 0;
    double lin_res = 0.0;

    // In the first step of this function, we solve for the incremental
    // displacement $d\mathbf{u}$.  To this end, we perform static
    // condensation to make
    //    $\mathbf{\mathsf{K}}_{\textrm{con}}
    //    = \bigl[ \mathbf{\mathsf{K}}_{uu} + \overline{\overline{\mathbf{\mathsf{K}}}}~ \bigr]$
    // and put
    // $\mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}$
    // in the original $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$ block.
    // That is, we make $\mathbf{\mathsf{K}}_{\textrm{store}}$.
    {
      assemble_sc();

      //              $
      //      \mathsf{\mathbf{A}}_{\widetilde{J}}
      //      =
      //              \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
      //              \mathsf{\mathbf{F}}_{\widetilde{p}}
      //              $
      tangent_matrix.block(p_dof, J_dof).vmult(A.block(J_dof),
                                               system_rhs.block(p_dof));
      //      $
      //      \mathsf{\mathbf{B}}_{\widetilde{J}}
      //      =
      //      \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
      //      \mathsf{\mathbf{F}}_{\widetilde{p}}
      //      $
      tangent_matrix.block(J_dof, J_dof).vmult(B.block(J_dof),
                                               A.block(J_dof));
      //      $
      //      \mathsf{\mathbf{A}}_{\widetilde{J}}
      //      =
      //      \mathsf{\mathbf{F}}_{\widetilde{J}}
      //      -
      //      \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
      //      \mathsf{\mathbf{F}}_{\widetilde{p}}
      //      $
      A.block(J_dof) = system_rhs.block(J_dof);
      A.block(J_dof) -= B.block(J_dof);
      //      $
      //      \mathsf{\mathbf{A}}_{\widetilde{J}}
      //      =
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{J} \widetilde{p}}
      //      [
      //      \mathsf{\mathbf{F}}_{\widetilde{J}}
      //      -
      //      \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
      //      \mathsf{\mathbf{F}}_{\widetilde{p}}
      //      ]
      //      $
      tangent_matrix.block(p_dof, J_dof).Tvmult(A.block(p_dof),
                                                A.block(J_dof));
      //      $
      //      \mathsf{\mathbf{A}}_{\mathbf{u}}
      //      =
      //      \mathsf{\mathbf{K}}_{\mathbf{u} \widetilde{p}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{J} \widetilde{p}}
      //      [
      //      \mathsf{\mathbf{F}}_{\widetilde{J}}
      //      -
      //      \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
      //      \mathsf{\mathbf{F}}_{\widetilde{p}}
      //      ]
      //      $
      tangent_matrix.block(u_dof, p_dof).vmult(A.block(u_dof),
                                               A.block(p_dof));
      //      $
      //      \mathsf{\mathbf{F}}_{\text{con}}
      //      =
      //      \mathsf{\mathbf{F}}_{\mathbf{u}}
      //      -
      //      \mathsf{\mathbf{K}}_{\mathbf{u} \widetilde{p}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{J} \widetilde{p}}
      //      [
      //      \mathsf{\mathbf{F}}_{\widetilde{J}}
      //      -
      //      \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
      //      \mathsf{\mathbf{K}}_{\widetilde{p}}
      //      ]
      //      $
      system_rhs.block(u_dof) -= A.block(u_dof);

      timer.enter_subsection("Linear solver");
      std::cout << " SLV " << std::flush;
      if (parameters.type_lin == "CG")
        {
          const int solver_its = tangent_matrix.block(u_dof, u_dof).m()
                                 * parameters.max_iterations_lin;
          const double tol_sol = parameters.tol_lin
                                 * system_rhs.block(u_dof).l2_norm();

          SolverControl solver_control(solver_its, tol_sol);

          GrowingVectorMemory<Vector<double> > GVM;
          SolverCG<Vector<double> > solver_CG(solver_control, GVM);

          // We've chosen by default a SSOR preconditioner as it appears to
          // provide the fastest solver convergence characteristics for this
          // problem on a single-thread machine.  However, this might not be
          // true for different problem sizes.
          PreconditionSelector<SparseMatrix<double>, Vector<double> >
          preconditioner (parameters.preconditioner_type,
                          parameters.preconditioner_relaxation);
          preconditioner.use_matrix(tangent_matrix.block(u_dof, u_dof));

          solver_CG.solve(tangent_matrix.block(u_dof, u_dof),
                          newton_update.block(u_dof),
                          system_rhs.block(u_dof),
                          preconditioner);

          lin_it = solver_control.last_step();
          lin_res = solver_control.last_value();
        }
      else if (parameters.type_lin == "Direct")
        {
          // Otherwise if the problem is small
          // enough, a direct solver can be
          // utilised.
          SparseDirectUMFPACK A_direct;
          A_direct.initialize(tangent_matrix.block(u_dof, u_dof));
          A_direct.vmult(newton_update.block(u_dof), system_rhs.block(u_dof));

          lin_it = 1;
          lin_res = 0.0;
        }
      else
        Assert (false, ExcMessage("Linear solver type not implemented"));

      timer.leave_subsection();
    }

    // Now that we have the displacement update, distribute the constraints
    // back to the Newton update:
    constraints.distribute(newton_update);

    timer.enter_subsection("Linear solver postprocessing");
    std::cout << " PP " << std::flush;

    // The next step after solving the displacement
    // problem is to post-process to get the
    // dilatation solution from the
    // substitution:
    //    $
    //     d \widetilde{\mathbf{\mathsf{J}}}
    //      = \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1} \bigl[
    //       \mathbf{\mathsf{F}}_{\widetilde{p}}
    //     - \mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
    //      \bigr]
    //    $
    {
      //      $
      //      \mathbf{\mathsf{A}}_{\widetilde{p}}
      //      =
      //      \mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
      //      $
      tangent_matrix.block(p_dof, u_dof).vmult(A.block(p_dof),
                                               newton_update.block(u_dof));
      //      $
      //      \mathbf{\mathsf{A}}_{\widetilde{p}}
      //      =
      //      -\mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
      //      $
      A.block(p_dof) *= -1.0;
      //      $
      //      \mathbf{\mathsf{A}}_{\widetilde{p}}
      //      =
      //      \mathbf{\mathsf{F}}_{\widetilde{p}}
      //      -\mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
      //      $
      A.block(p_dof) += system_rhs.block(p_dof);
      //      $
      //      d\mathbf{\mathsf{\widetilde{J}}}
      //      =
      //      \mathbf{\mathsf{K}}^{-1}_{\widetilde{p}\widetilde{J}}
      //      [
      //      \mathbf{\mathsf{F}}_{\widetilde{p}}
      //      -\mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
      //      ]
      //      $
      tangent_matrix.block(p_dof, J_dof).vmult(newton_update.block(J_dof),
                                               A.block(p_dof));
    }

    // we insure here that any Dirichlet constraints
    // are distributed on the updated solution:
    constraints.distribute(newton_update);

    // Finally we solve for the pressure
    // update with the substitution:
    //    $
    //    d \widetilde{\mathbf{\mathsf{p}}}
    //     =
    //    \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
    //    \bigl[
    //     \mathbf{\mathsf{F}}_{\widetilde{J}}
    //      - \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
    //    d \widetilde{\mathbf{\mathsf{J}}}
    //    \bigr]
    //    $
    {
      //      $
      //      \mathsf{\mathbf{A}}_{\widetilde{J}}
      //       =
      //      \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
      //      d \widetilde{\mathbf{\mathsf{J}}}
      //      $
      tangent_matrix.block(J_dof, J_dof).vmult(A.block(J_dof),
                                               newton_update.block(J_dof));
      //      $
      //      \mathsf{\mathbf{A}}_{\widetilde{J}}
      //       =
      //      -\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
      //      d \widetilde{\mathbf{\mathsf{J}}}
      //      $
      A.block(J_dof) *= -1.0;
	  //      $
      //      \mathsf{\mathbf{A}}_{\widetilde{J}}
      //       =
      //      \mathsf{\mathbf{F}}_{\widetilde{J}}
      //      -
      //      \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
      //      d \widetilde{\mathbf{\mathsf{J}}}
      //      $
      A.block(J_dof) += system_rhs.block(J_dof);
      // and finally....
      //    $
      //    d \widetilde{\mathbf{\mathsf{p}}}
      //     =
      //    \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
      //    \bigl[
      //     \mathbf{\mathsf{F}}_{\widetilde{J}}
      //      - \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
      //    d \widetilde{\mathbf{\mathsf{J}}}
      //    \bigr]
      //    $
      tangent_matrix.block(p_dof, J_dof).Tvmult(newton_update.block(p_dof),
                                                A.block(J_dof));
    }

    // We are now at the end, so we distribute all
    // constrained dofs back to the Newton
    // update:
    constraints.distribute(newton_update);

    timer.leave_subsection();

    return std::make_pair(lin_it, lin_res);
  }

// @sect4{Solid::assemble_system_SC}

// The static condensation process could be performed at a global level but we
// need the inverse of one of the blocks. However, since the pressure and
// dilatation variables are discontinuous, the static condensation (SC)
// operation can be done on a per-cell basis and we can produce the inverse of
// the block-diagonal $ \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}$
// block by inverting the local blocks. We can again
// use TBB to do this since each operation will be independent of one another.
//
// Using the TBB via the WorkStream class, we assemble the contributions to form
//  $
//  \mathbf{\mathsf{K}}_{\textrm{con}}
//  = \bigl[ \mathbf{\mathsf{K}}_{uu} + \overline{\overline{\mathbf{\mathsf{K}}}}~ \bigr]
//  $
// from each element's contributions. These
// contributions are then added to the global stiffness matrix. Given this
// description, the following two functions should be clear:
  template <int dim>
  void Solid<dim>::assemble_sc()
  {
    timer.enter_subsection("Perform static condensation");
    std::cout << " ASM_SC " << std::flush;

    PerTaskData_SC per_task_data(dofs_per_cell, element_indices_u.size(),
                                 element_indices_p.size(),
                                 element_indices_J.size());
    ScratchData_SC scratch_data;

    WorkStream::run(dof_handler_ref.begin_active(),
                    dof_handler_ref.end(),
                    *this,
                    &Solid::assemble_sc_one_cell,
                    &Solid::copy_local_to_global_sc,
                    scratch_data,
                    per_task_data);

    timer.leave_subsection();
  }


  template <int dim>
  void Solid<dim>::copy_local_to_global_sc(const PerTaskData_SC &data)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        tangent_matrix.add(data.local_dof_indices[i],
                           data.local_dof_indices[j],
                           data.cell_matrix(i, j));
  }


// Now we describe the static condensation process.  As per usual, we must
// first find out which global numbers the degrees of freedom on this cell
// have and reset some data structures:
  template <int dim>
  void
  Solid<dim>::assemble_sc_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                   ScratchData_SC &scratch,
                                   PerTaskData_SC &data)
  {
    data.reset();
    scratch.reset();
    cell->get_dof_indices(data.local_dof_indices);

    // We now extract the contribution of the dofs associated with the current
    // cell to the global stiffness matrix.  The discontinuous nature of the
    // $\widetilde{p}$ and $\widetilde{J}$ interpolations mean that their is
    // no coupling of the local contributions at the global level. This is not
    // the case with the u dof.  In other words,
    // $\mathsf{\mathbf{k}}_{\widetilde{J} \widetilde{p}}$,
    // $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{p}}$ and
    // $\mathsf{\mathbf{k}}_{\widetilde{J} \widetilde{p}}$, when extracted
    // from the global stiffness matrix are the element contributions.  This
    // is not the case for $\mathsf{\mathbf{k}}_{\mathbf{u} \mathbf{u}}$
    //
    // Note: A lower-case symbol is used to denote element stiffness matrices.

    // Currently the matrix corresponding to
    // the dof associated with the current element
    // (denoted somewhat loosely as $\mathsf{\mathbf{k}}$)
    // is of the form:
    // @f{align*}
    //    \begin{bmatrix}
    //       \mathbf{\mathsf{k}}_{uu}  &  \mathbf{\mathsf{k}}_{u\widetilde{p}}    & \mathbf{0}
    //    \\ \mathbf{\mathsf{k}}_{\widetilde{p}u} & \mathbf{0}  &  \mathbf{\mathsf{k}}_{\widetilde{p}\widetilde{J}}
    //    \\ \mathbf{0}  &  \mathbf{\mathsf{k}}_{\widetilde{J}\widetilde{p}}  & \mathbf{\mathsf{k}}_{\widetilde{J}\widetilde{J}}
    //    \end{bmatrix}
    // @f}
    //
    // We now need to modify it such that it appear as
    // @f{align*}
    //    \begin{bmatrix}
    //       \mathbf{\mathsf{k}}_{\textrm{con}}   & \mathbf{\mathsf{k}}_{u\widetilde{p}}    & \mathbf{0}
    //    \\ \mathbf{\mathsf{k}}_{\widetilde{p}u} & \mathbf{0} & \mathbf{\mathsf{k}}_{\widetilde{p}\widetilde{J}}^{-1}
    //    \\ \mathbf{0} & \mathbf{\mathsf{k}}_{\widetilde{J}\widetilde{p}} & \mathbf{\mathsf{k}}_{\widetilde{J}\widetilde{J}}
    //    \end{bmatrix}
    // @f}
    // with $\mathbf{\mathsf{k}}_{\textrm{con}} = \bigl[ \mathbf{\mathsf{k}}_{uu} +\overline{\overline{\mathbf{\mathsf{k}}}}~ \bigr]$
    // where
    // $               \overline{\overline{\mathbf{\mathsf{k}}}} :=
    // \mathbf{\mathsf{k}}_{u\widetilde{p}} \overline{\mathbf{\mathsf{k}}} \mathbf{\mathsf{k}}_{\widetilde{p}u}
    // $
    // and
    // $
    //    \overline{\mathbf{\mathsf{k}}} =
    //     \mathbf{\mathsf{k}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{k}}_{\widetilde{J}\widetilde{J}}
    //    \mathbf{\mathsf{k}}_{\widetilde{p}\widetilde{J}}^{-1}
    // $.
    //
    // At this point, we need to take note of
    // the fact that global data already exists
    // in the $\mathsf{\mathbf{K}}_{uu}$,
    // $\mathsf{\mathbf{K}}_{\widetilde{p} \widetilde{J}}$
    // and
    //  $\mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{p}}$
    // sub-blocks.  So if we are to modify them, we must account for the data
    // that is already there (i.e. simply add to it or remove it if
    // necessary).  Since the copy_local_to_global operation is a "+="
    // operation, we need to take this into account
    //
    // For the $\mathsf{\mathbf{K}}_{uu}$ block in particular, this means that
    // contributions have been added from the surrounding cells, so we need to
    // be careful when we manipulate this block.  We can't just erase the
    // sub-blocks.
    //
    // This is the strategy we will employ to get the sub-blocks we want:
    //
    // - $ {\mathbf{\mathsf{k}}}_{\textrm{store}}$:
    // Since we don't have access to $\mathsf{\mathbf{k}}_{uu}$,
    // but we know its contribution is added to
    // the global $\mathsf{\mathbf{K}}_{uu}$ matrix, we just want
    // to add the element wise
    // static-condensation $\overline{\overline{\mathbf{\mathsf{k}}}}$.
    //
    // - $\mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}$:
    //                      Similarly, $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$ exists in
    //          the subblock. Since the copy
    //          operation is a += operation, we
    //          need to subtract the existing
    //          $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$
    //                      submatrix in addition to
    //          "adding" that which we wish to
    //          replace it with.
    //
    // - $\mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{p}}$:
    //              Since the global matrix
    //          is symmetric, this block is the
    //          same as the one above and we
    //          can simply use
    //              $\mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}$
    //          as a substitute for this one.
    //
    // We first extract element data from the
    // system matrix. So first we get the
    // entire subblock for the cell, then
    // extract $\mathsf{\mathbf{k}}$
    // for the dofs associated with
    // the current element
    data.k_orig.extract_submatrix_from(tangent_matrix,
                                       data.local_dof_indices,
                                       data.local_dof_indices);
    // and next the local matrices for
	// $\mathsf{\mathbf{k}}_{ \widetilde{p} \mathbf{u}}$
	// $\mathsf{\mathbf{k}}_{ \widetilde{p} \widetilde{J}}$
	// and
	// $\mathsf{\mathbf{k}}_{ \widetilde{J} \widetilde{J}}$:
    data.k_pu.extract_submatrix_from(data.k_orig,
                                     element_indices_p,
                                     element_indices_u);
    data.k_pJ.extract_submatrix_from(data.k_orig,
                                     element_indices_p,
                                     element_indices_J);
    data.k_JJ.extract_submatrix_from(data.k_orig,
                                     element_indices_J,
                                     element_indices_J);

    // To get the inverse of $\mathsf{\mathbf{k}}_{\widetilde{p}
    // \widetilde{J}}$, we invert it directly.  This operation is relatively
    // inexpensive since $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$
    // since block-diagonal.

    data.k_pJ_inv.invert(data.k_pJ);

    // Now we can make condensation terms to
    // add to the $\mathsf{\mathbf{k}}_{\mathbf{u} \mathbf{u}}$
    // block and put them in
    // the cell local matrix
    //    $
    //    \mathsf{\mathbf{A}}
    //    =
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
    //    \mathsf{\mathbf{k}}_{\widetilde{p} \mathbf{u}}
    //    $:
    data.k_pJ_inv.mmult(data.A, data.k_pu);
    //      $
    //      \mathsf{\mathbf{B}}
    //      =
    //      \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{J}}
    //      \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
    //      \mathsf{\mathbf{k}}_{\widetilde{p} \mathbf{u}}
    //      $
    data.k_JJ.mmult(data.B, data.A);
    //    $
    //    \mathsf{\mathbf{C}}
    //    =
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{p}}
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{J}}
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
    //    \mathsf{\mathbf{k}}_{\widetilde{p} \mathbf{u}}
    //    $
    data.k_pJ_inv.Tmmult(data.C, data.B);
    //    $
    //    \overline{\overline{\mathsf{\mathbf{k}}}}
    //    =
    //    \mathsf{\mathbf{k}}_{\mathbf{u} \widetilde{p}}
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{p}}
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{J}}
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
    //    \mathsf{\mathbf{k}}_{\widetilde{p} \mathbf{u}}
    //    $
    data.k_pu.Tmmult(data.k_bbar, data.C);
    data.k_bbar.scatter_matrix_to(element_indices_u,
                                  element_indices_u,
                                  data.cell_matrix);

    // Next we place
    // $\mathsf{\mathbf{k}}^{-1}_{ \widetilde{p} \widetilde{J}}$
    // in the
    // $\mathsf{\mathbf{k}}_{ \widetilde{p} \widetilde{J}}$
    // block for post-processing.  Note again
    // that we need to remove the
    // contribution that already exists there.
    data.k_pJ_inv.add(-1.0, data.k_pJ);
    data.k_pJ_inv.scatter_matrix_to(element_indices_p,
                                    element_indices_J,
                                    data.cell_matrix);
  }

// @sect4{Solid::output_results}
// Here we present how the results are written to file to be viewed
// using ParaView or Visit. The method is similar to that shown in previous
// tutorials so will not be discussed in detail.
template <int dim>
void Solid<dim>::output_deformation_profile(std::string suffix)
{
    DataOut<dim> data_out;
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(dim,
                                  DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);

    std::vector<std::string> solution_name(dim, "displacement");
    solution_name.push_back("pressure");
    solution_name.push_back("dilatation");

    data_out.attach_dof_handler(dof_handler_ref);
    data_out.add_data_vector(solution_n,
                             solution_name,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);

    // Since we are dealing with a large deformation problem, it would be nice
    // to display the result on a displaced grid!  The MappingQEulerian class
    // linked with the DataOut class provides an interface through which this
    // can be achieved without physically moving the grid points in the
    // Triangulation object ourselves.  We first need to copy the solution to
    // a temporary vector and then create the Eulerian mapping. We also
    // specify the polynomial degree to the DataOut object in order to produce
    // a more refined output data set when higher order polynomials are used.
    Vector<double> soln(solution_n.size());
    for (unsigned int i = 0; i < soln.size(); ++i)
      soln(i) = solution_n(i);
    MappingQEulerian<dim> q_mapping(degree, soln, dof_handler_ref);
    data_out.build_patches(q_mapping, degree);

    std::ostringstream filename;
    static int out_i = 0;
    filename << "output/solution-" << out_i << suffix << ".vtk";
    out_i++;

    std::ofstream output(filename.str().c_str());
    data_out.write_vtk(output);

    //============================================================
    
    Point <dim> point;
    Vector<double> value(dim+2);

    for(int i=0; i<mYstep; i+=10) {
        for(int j=0; j<mXstep; j+=10) {
            point(0) = (double)j;
            point(1) = (double)i;
            try {
                VectorTools::point_value(dof_handler_ref, solution_n, point, value);
                tissue_displacement_x[i][j] = value(0);
                tissue_displacement_y[i][j] = value(1);
            //cout << "tissue_displacement" << value(0) << ", " << value(1) << ", "
            //                              << value(2) << ", " << value(3) << endl;
            } catch (...) {
                cout << "Caught exception from point_value" << endl; 
            }
        }
    }
 

    for(int i=0; i<mYstep; i++) {
        for(int j=0; j<mXstep; j++) {
            point(0) = (double)j;
            point(1) = (double)i;
            if (isOnWoundEdge(j, i, mXstep, mYstep)) {
                VectorTools::point_value(dof_handler_ref, solution_n, point, value);
                tissue_displacement_x[i][j] = value(0);
                tissue_displacement_y[i][j] = value(1);
            }
        }
    }

}

template <int dim>
void Solid<dim>::Output_woundcontour(){
    Mat_DP woundcontour(mYstep,mXstep);
    DP tdx,tdy;
    for(size_t i=0; i<mYstep; i++){
        for(size_t j=0; j<mXstep; j++){
            woundcontour[i][j] = 0;
        }
    }

    static int out_i = 0;
    for(size_t i=0; i<mYstep; i++){
        for(size_t j=0; j<mXstep; j++){
            if( isOnWoundEdge(j,i,mXstep,mYstep) ){
                tdx = tissue_displacement_x[i][j];
                tdy = tissue_displacement_y[i][j];
                if(i+(int)tdy>=0 && i+(int)tdy<mYstep && j+(int)tdx>=0 && j+(int)tdy<mXstep){
                    woundcontour[i+(int)tdy][j+(int)tdx] = 2.0;
                    if(out_i%10 == 0){
                        mWoundContourHistory[std::make_pair(i+(int)tdy, j+(int)tdx)] = out_i;
                    }
                }
            }
        }
    }

    char file_name[21];
    sprintf(file_name, "output/cont%05d.BMP", out_i++);
    BMP::output_BMP(file_name, 14, woundcontour, mXstep, mYstep);

    if(out_i%50 == 0){
        sprintf(file_name, "output/cont_aggregate.BMP");
        BMP::output_BMP3(file_name, 14, mWoundContourHistory, mXstep, mYstep);
    }
    return;
}

}//end of namespace Step44

// @sect3{Main function}
// Lastly we provide the main driver function which appears
// no different to the other tutorials.
/*
int main (int argc, char *argv[])
{
  using namespace dealii;
  using namespace Step44;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  Tensor<2,3,double> tensorA2, tensorB2;
  SymmetricTensor<2,3,double> s2_tensorA, s2_tensorB;
  s2_tensorA[0][0] = 1; s2_tensorB[0][0] = 2;
  s2_tensorA[0][1] = 2; s2_tensorB[0][1] = 3;
  s2_tensorA[0][2] = 3; s2_tensorB[0][2] = 4;
  s2_tensorA[1][1] = 4; s2_tensorB[1][0] = 5;
  s2_tensorA[1][2] = 5; s2_tensorB[1][1] = 6;
  s2_tensorA[2][2] = 6; s2_tensorB[1][2] = 7;
  tensorA2 = s2_tensorA; tensorB2 = s2_tensorB;
  std::cout<<"tensor converted from symmetric tensor..."<<std::endl;
  std::cout<<tensorA2[0][0]<<" "<<tensorA2[0][1]<<" "<<tensorA2[0][2]<<std::endl;
  std::cout<<tensorA2[1][0]<<" "<<tensorA2[1][1]<<" "<<tensorA2[1][2]<<std::endl;
  std::cout<<tensorA2[2][0]<<" "<<tensorA2[2][1]<<" "<<tensorA2[2][2]<<std::endl;
  SymmetricTensor<4,3,double> s4_tensor = outer_product(s2_tensorA,s2_tensorB);
  tensorA2[1][0] = 7;   tensorB2[1][0] = 8;
  tensorA2[2][0] = 8;   tensorB2[2][0] = 9;
  tensorA2[2][1] = 9;   tensorB2[2][1] = 10;
  Tensor<4,3,double> tensor4 = outer_product(tensorA2, tensorB2);
  double s_tensor_contraction = s2_tensorA*s4_tensor*s2_tensorA;
  Tensor<2,3,double> tensor_contraction1 = double_contract<2,0,3,1>(tensor4, tensorA2);
  double tensor_contraction2 = double_contract<0,0,1,1>(tensorA2,tensor_contraction1);
  std::cout<<"rank 4 symmetric tensor contracted with rank 2 symmetric tensor..."<<s_tensor_contraction<<std::endl;
//  std::cout<<s_tensor_contraction[2][0][0][0]<<" "<<s_tensor_contraction[2][0][0][1]<<" "<<s_tensor_contraction[2][0][0][2]<<std::endl;
//  std::cout<<s_tensor_contraction[2][0][1][0]<<" "<<s_tensor_contraction[2][0][1][1]<<" "<<s_tensor_contraction[2][0][1][2]<<std::endl;
//  std::cout<<s_tensor_contraction[2][0][2][0]<<" "<<s_tensor_contraction[2][0][2][1]<<" "<<s_tensor_contraction[2][0][2][2]<<std::endl;
  std::cout<<"rank 4 tensor contracted with rank 2 tensor..."<<tensor_contraction2<<std::endl;
  std::cout<<tensor_contraction1[0][0]<<" "<<tensor_contraction1[0][1]<<" "<<tensor_contraction1[0][2]<<std::endl;
  std::cout<<tensor_contraction1[1][0]<<" "<<tensor_contraction1[1][1]<<" "<<tensor_contraction1[1][2]<<std::endl;
  std::cout<<tensor_contraction1[2][0]<<" "<<tensor_contraction1[2][1]<<" "<<tensor_contraction1[2][2]<<std::endl;
  
  
   deallog.depth_console(0);

   Solid<3> solid_3d("parameters.prm", 2000, 2000);
   solid_3d.run();
   return 0;   
}
*/
