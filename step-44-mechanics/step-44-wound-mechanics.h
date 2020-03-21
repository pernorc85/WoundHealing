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

#pragma once 
#ifndef _STEP44_WOUND_MECHANICS_
#define _STEP44_WOUND_MECHANICS_

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
#include "parameters.h"
#include "Material_Transverse_Isotropic.h"
#include "Material_Transverse_Isotropic.hh"

extern const int xstep;
extern const int ystep;

extern Mat_DP *fibroblast_density_ptr;
extern Mat_DP *speedfield_ptr;
extern Mat_DP *thetafield_ptr;


// We then stick everything that relates to this tutorial program into a
// namespace of its own, and import all the deal.II function and class names
// into it:
namespace Step44
{
  using namespace dealii;


//{collagen fiber direction}
template <int dim>
class Direction :  public Function<dim>
{
public:
    Direction ();
    Direction (const Mat_DP& collagen_direction_);
    void Set (const Mat_DP& collagen_direction_);

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &values) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> >   &value_list) const;
protected:
    Mat_DP collagen_direction;
};

//========================================
template <int dim>
class Cdensity :  public Function<dim>
{
public:
    Cdensity ();
    Cdensity (const Mat_DP& collagen_density_);
    void Set(const Mat_DP& collagen_density_);

    virtual void vector_value (const Point<dim> &p,
                               double   &value) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<double>   &value_list) const;

protected:
    Mat_DP collagen_density;
};

//=======================================================================
template <int dim>
class RightHandSide : public Function<dim>
{
public:
    RightHandSide ();

    virtual void vector_value (const Point<dim> &p,
                               Vector<double> &value) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> > &value_list) const;
};


// @sect4{Finite Element system}   

// @sect3{Time class}
// A simple class to store time data. Its functioning is transparent so no
// discussion is necessary. For simplicity we assume a constant time step
// size.
  class Time
  {
  public:
    Time (const double time_end,
          const double delta_t)
      :
      timestep(0),
      time_current(0.0),
      time_end(time_end),
      delta_t(delta_t)
    {}

    virtual ~Time()
    {}

    double current() const
    {
      return time_current;
    }
    double end() const
    {
      return time_end;
    }
    double get_delta_t() const
    {
      return delta_t;
    }
    void dec_delta_t()
    {
      delta_t *= 0.25;
    }
    void inc_delta_t()
    {
      delta_t *= 1.1;
    }
    unsigned int get_timestep() const
    {
      return timestep;
    }
    void increment()
    {
      time_current += delta_t;
      ++timestep;
    }

  private:
    unsigned int timestep;
    double       time_current;
    const double time_end;
    double delta_t;
  };


// @sect3{Quadrature point history}

// As seen in step-18, the <code>
// PointHistory </code> class offers a method for storing data at the
// quadrature points.  Here each quadrature point holds a pointer to a
// material description.  Thus, different material models can be used in
// different regions of the domain.  Among other data, we choose to store the
// Kirchhoff stress $\boldsymbol{\tau}$ and the tangent $J\mathfrak{c}$ for
// the quadrature points.
  template <int dim>
  class PointHistory
  {
  public:
    PointHistory()
      :
      material(NULL),
      F(StandardTensors<dim>::I),
      F_inv(StandardTensors<dim>::I),
      F_g(StandardTensors<dim>::I),
      F_e(StandardTensors<dim>::I),
      S(Tensor<2, dim>()),
      d2Psi_vol_dJ2(0.0),
      dPsi_vol_dJ(0.0),
      Cmaterial(Tensor<4, dim>())
    {}

    virtual ~PointHistory()
    {
      delete material;
      material = NULL;
    }

    // The first function is used to create a material object and to
    // initialize all tensors correctly: The second one updates the stored
    // values and stresses based on the current deformation measure
    // $\textrm{Grad}\mathbf{u}_{\textrm{n}}$, pressure $\widetilde{p}$ and
    // dilation $\widetilde{J}$ field values.

    void setup_lqp (Parameters::AllParameters &parameters, double cdensity, Tensor<1,dim> M_par) {
      material = new Material_Compressible_TransIsotropic_Three_Field<dim>(parameters.mu,
                                                                           parameters.nu, 
                                                                           cdensity * parameters.k1, 
                                                                           parameters.k2,
                                                                         M_par);
                                                                         
      update_values(Tensor<2, dim>(), 0.0, 1.0);
    }

    // To this end, we calculate the deformation gradient $\mathbf{F}$ from
    // the displacement gradient $\textrm{Grad}\ \mathbf{u}$, i.e.
    // $\mathbf{F}(\mathbf{u}) = \mathbf{I} + \textrm{Grad}\ \mathbf{u}$ and
    // then let the material model associated with this quadrature point
    // update itself. When computing the deformation gradient, we have to take
    // care with which data types we compare the sum $\mathbf{I} +
    // \textrm{Grad}\ \mathbf{u}$: Since $I$ has data type SymmetricTensor,
    // just writing <code>I + Grad_u_n</code> would convert the second
    // argument to a symmetric tensor, perform the sum, and then cast the
    // result to a Tensor (i.e., the type of a possibly nonsymmetric
    // tensor). However, since <code>Grad_u_n</code> is nonsymmetric in
    // general, the conversion to SymmetricTensor will fail. We can avoid this
    // back and forth by converting $I$ to Tensor first, and then performing
    // the addition as between nonsymmetric tensors:


    void update_values (const Tensor<2, dim> &Grad_u_n,
                        const double p_tilde,
                        const double J_tilde)
    {
      F = Tensor<2, dim>(StandardTensors<dim>::I) + Grad_u_n;
      F_e = F * invert(F_g);
      material->update_material_data(F_e, p_tilde, J_tilde/determinant(F_g));

      // The material has been updated so we now calculate the Kirchhoff
      // stress $\mathbf{\tau}$, the tangent $J\mathfrak{c}$ and the first and
      // second derivatives of the volumetric free energy.
      //
      // We also store the inverse of the deformation gradient since we
      // frequently use it:
      F_inv = invert(F);
      S = material->get_S();
      Cmaterial = material->get_Cmaterial();
      dPsi_vol_dJ = material->get_dPsi_vol_dJ();
      d2Psi_vol_dJ2 = material->get_d2Psi_vol_dJ2();
    }

    void update_values_for_growth( const Tensor<2, dim> &F_g_)
    {
      F_g = F_g * F_g_;
      F_e = F * invert(F_g);//F is already stored in qph
      double J_tilde = get_J_tilde();
      double p_tilde = get_p_tilde();
      material->update_material_data(F_e, p_tilde, J_tilde/determinant(F_g));
    
      S = material->get_S();
      Cmaterial = material->get_Cmaterial();
      dPsi_vol_dJ = material->get_dPsi_vol_dJ();
      d2Psi_vol_dJ2 = material->get_d2Psi_vol_dJ2();
    }

    void update_values_for_material_property(const Tensor<1, dim> &collagen_orientation,
                                             double collagen_density)
    {
      material->update_material_data(collagen_orientation, collagen_density);

      S = material->get_S();
      Cmaterial = material->get_Cmaterial();
      dPsi_vol_dJ = material->get_dPsi_vol_dJ();
      d2Psi_vol_dJ2 = material->get_d2Psi_vol_dJ2();
    }

    // We offer an interface to retrieve certain data.  Here are the kinematic
    // variables:
    double get_J_tilde() const
    {
      return material->get_J_tilde();
    }

    double get_det_Fe() const
    {
      return material->get_det_Fe();
    }

    const Tensor<2, dim>& get_F() const
    {
      return F;
    }

    const Tensor<2, dim>& get_F_inv() const
    {
      return F_inv;
    }

    // ...and the kinetic variables.  These are used in the material and
    // global tangent matrix and residual assembly operations:
    double get_p_tilde() const
    {
      return material->get_p_tilde();
    }

    const Tensor<2, dim> &get_S() const
    {
      return S;
    }

    double get_dPsi_vol_dJ() const
    {
      return dPsi_vol_dJ;
    }

    double get_d2Psi_vol_dJ2() const
    {
      return d2Psi_vol_dJ2;
    }

    // And finally the tangent:
    const Tensor<4, dim> &get_Cmaterial() const
    {
      return Cmaterial;
    }
    // In terms of member functions, this class stores for the quadrature
    // point it represents a copy of a material type in case different
    // materials are used in different regions of the domain, as well as the
    // inverse of the deformation gradient...
  private:
    Material_Compressible_TransIsotropic_Three_Field<dim>* material;

    Tensor<2, dim> F;
    Tensor<2, dim> F_inv;
    Tensor<2, dim> F_g;
    Tensor<2, dim> F_e;

    // ... and stress-type variables along with the tangent $J\mathfrak{c}$:
    Tensor<2, dim> S;
    double                  d2Psi_vol_dJ2;
    double                  dPsi_vol_dJ;

    Tensor<4, dim> Cmaterial; 
  };


// @sect3{Quasi-static quasi-incompressible finite-strain solid}

// The Solid class is the central class in that it represents the problem at
// hand. It follows the usual scheme in that all it really has is a
// constructor, destructor and a <code>run()</code> function that dispatches
// all the work to private functions of this class:
template <int dim>
class Solid {
public:
    Solid(const std::string &input_file, int xstep, int ystep);

    virtual
    ~Solid();

    void
    run();

    void 
    Setup(Mat_DP &collagen_orientation, Mat_DP &collagen_density);

    void
    RunNonlinear(Mat_DP &collagen_orientation, Mat_DP &collagen_density);

public:
    Mat_DP tissue_displacement_x;
    Mat_DP tissue_displacement_y;

private:

    // In the private section of this class, we first forward declare a number
    // of objects that are used in parallelizing work using the WorkStream
    // object (see the @ref threads module for more information on this).
    //
    // We declare such structures for the computation of tangent (stiffness)
    // matrix, right hand side, static condensation, and for updating
    // quadrature points:

    struct PerTaskData_K;
    struct ScratchData_K;

    struct PerTaskData_RHS;
    struct ScratchData_RHS;

    struct PerTaskData_SC;
    struct ScratchData_SC;

    struct PerTaskData_UQPH;
    struct ScratchData_UQPH;

    // We start the collection of member functions with one that builds the
    // grid:
    void
    make_grid();

    // Set up the finite element system to be solved:
    void
    system_setup();

    void
    determine_component_extractors();

    // Several functions to assemble the system and right hand side matrices
    // using multithreading. Each of them comes as a wrapper function, one
    // that is executed to do the work in the WorkStream model on one cell,
    // and one that copies the work done on this one cell into the global
    // object that represents it:
    void
    assemble_system_tangent();

    void
    assemble_system_tangent_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                     ScratchData_K &scratch,
                                     PerTaskData_K &data);

    void
    copy_local_to_global_K(const PerTaskData_K &data);

    void
    assemble_system_rhs();

    void
    assemble_system_rhs_one_cell(const typename DoFHandler<dim>::active_cell_iterator & cell,
                                 ScratchData_RHS &scratch,
                                 PerTaskData_RHS &data);

    void
    copy_local_to_global_rhs(const PerTaskData_RHS &data);

    //********************** for perturbation
    void
    assemble_perturbation();

    void
    assemble_system_perturb_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                   ScratchData_RHS &scratch,
                                   PerTaskData_RHS &data);

    void
    copy_local_to_global_perturb(const PerTaskData_RHS &data);

//***********************************************************
    void
    assemble_sc();

    void
    assemble_sc_one_cell(const typename DoFHandler<dim>::active_cell_iterator & cell,
                         ScratchData_SC &scratch,
                         PerTaskData_SC &data);

    void
    copy_local_to_global_sc(const PerTaskData_SC &data);

    // Apply Dirichlet boundary conditions on the displacement field

    void
    make_constraints(const int &it_nr);

    // Create and update the quadrature points. Here, no data needs to be
    // copied into a global object, so the copy_local_to_global function is
    // empty:

    void
    setup_qph();

    void
    update_qph_incremental(const BlockVector<double> &solution_delta);

    void
    update_qph_incremental_one_cell(const typename DoFHandler<dim>::active_cell_iterator & cell,
                                    ScratchData_UQPH &scratch,
                                    PerTaskData_UQPH &data);



    void
    copy_local_to_global_UQPH(const PerTaskData_UQPH &data)
    {}
    void
    update_qph_for_growth();
    void
    update_qph_for_material_property(const Mat_DP& collagen_direction_in_rad,
                                     const Mat_DP& collagen_density);

    // Solve for the displacement using a Newton-Raphson method. We break this
    // function into the nonlinear loop and the function that solves the
    // linearized Newton-Raphson step:
    void
    solve_nonlinear_timestep(BlockVector<double> &solution_delta);

    std::pair<unsigned int, double>
    solve_linear_system(BlockVector<double> &newton_update);

    // Solution retrieval as well as post-processing and writing data to file:
    BlockVector<double>
    get_total_solution(const BlockVector<double> &solution_delta) const;

    void
    output_deformation_profile(std::string suffix);

    void
    output_results_temp(int iteration, BlockVector<double> &solution_temp) const
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
      data_out.add_data_vector(solution_temp,
                             solution_name,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);
      Vector<double> soln(solution_temp.size());
      for (unsigned int i = 0; i < soln.size(); ++i)
        soln(i) = solution_temp(i);
      MappingQEulerian<dim> q_mapping(degree, dof_handler_ref, soln);
      data_out.build_patches(q_mapping, degree);

      std::ostringstream filename;
      filename << "solution-" << time.get_timestep() << "-" << iteration << ".vtk";

      std::ofstream output(filename.str().c_str());
      data_out.write_vtk(output);
    }

    // Finally, some member variables that describe the current state: A
    // collection of the parameters used to describe the problem setup...
    Parameters::AllParameters        parameters;

    // ...the volume of the reference and current configurations...
    double                           vol_reference;
    double                           vol_current;

    // ...and description of the geometry on which the problem is solved:
    Triangulation<dim>               triangulation;

    // Also, keep track of the current time and the time spent evaluating
    // certain functions
    Time                             time;
    TimerOutput                      timer;

    // A storage object for quadrature point information.  See step-18 for
    // more on this:
    std::vector<PointHistory<dim> >  quadrature_point_history;

    // A description of the finite-element system including the displacement
    // polynomial degree, the degree-of-freedom handler, number of DoFs per
    // cell and the extractor objects used to retrieve information from the
    // solution vectors:
    const unsigned int               degree;
    const FESystem<dim>              fe;
    DoFHandler<dim>                  dof_handler_ref;
    const unsigned int               dofs_per_cell;
    const FEValuesExtractors::Vector u_fe;
    const FEValuesExtractors::Scalar p_fe;
    const FEValuesExtractors::Scalar J_fe;

    // Description of how the block-system is arranged. There are 3 blocks,
    // the first contains a vector DOF $\mathbf{u}$ while the other two
    // describe scalar DOFs, $\widetilde{p}$ and $\widetilde{J}$.
    static const unsigned int        n_blocks = 3;
    static const unsigned int        n_components = dim + 2;
    static const unsigned int        first_u_component = 0;
    static const unsigned int        p_component = dim;
    static const unsigned int        J_component = dim + 1;

    enum
    {
      u_dof = 0,
      p_dof = 1,
      J_dof = 2
    };

    std::vector<types::global_dof_index>        dofs_per_block;
    std::vector<types::global_dof_index>        element_indices_u;
    std::vector<types::global_dof_index>        element_indices_p;
    std::vector<types::global_dof_index>        element_indices_J;

    // Rules for Gauss-quadrature on both the cell and faces. The number of
    // quadrature points on both cells and faces is recorded.
    const QGauss<dim>                qf_cell;
    const QGauss<dim - 1>            qf_face;
    const unsigned int               n_q_points;
    const unsigned int               n_q_points_f;

    // Objects that store the converged solution and right-hand side vectors,
    // as well as the tangent matrix. There is a ConstraintMatrix object used
    // to keep track of constraints.  We make use of a sparsity pattern
    // designed for a block system.

    ConstraintMatrix                 constraints;
    BlockSparsityPattern             sparsity_pattern;
    BlockSparseMatrix<double>        tangent_matrix;
    BlockVector<double>              system_rhs;
    BlockVector<double>              solution_n;
    BlockVector<double>              solution_perturbation;
   
    int mXstep, mYstep; 
    Cdensity<dim>                    mCollagenDensity;
    Direction<dim>                   mCollagenDirection;
    RightHandSide<dim>               mFibroblastForce;

    // Then define a number of variables to store norms and update norms and
    // normalisation factors.
    struct Errors
    {
      Errors()
        :
        norm(1.0), u(1.0), p(1.0), J(1.0)
      {}

      void reset()
      {
        norm = 1.0;
        u = 1.0;
        p = 1.0;
        J = 1.0;
      }
      void normalise(const Errors &rhs)
      {
        if (rhs.norm != 0.0)
          norm /= rhs.norm;
        if (rhs.u != 0.0)
          u /= rhs.u;
        if (rhs.p != 0.0)
          p /= rhs.p;
        if (rhs.J != 0.0)
          J /= rhs.J;
      }

      double norm, u, p, J;
    };

    Errors error_residual, error_residual_0, error_residual_norm, error_update,
           error_update_0, error_update_norm;

    // Methods to calculate error measures
    void
    get_error_residual(Errors &error_residual);

    void
    get_error_update(const BlockVector<double> &newton_update,
                     Errors &error_update);

    std::pair<double, double>
    get_error_dilation();

    // Print information to screen in a pleasing way...
    static
    void
    print_conv_header();

    void
    print_conv_footer();
};

} //end of namespace Step44

#endif
