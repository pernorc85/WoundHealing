//step8_wound_mechanics_V6.cpp
//    VERSION: 6 (1 November 2015)
//    PURPOSE: Header file for FEM of force balance equation in wound contraction model.
//             Based on step-8.cc of deal.II package.
//             ($Id: step-8.cc 18448 2009-03-04 22:01:56Z kronbichler $ 
//             Author: Wolfgang Bangerth, University of Heidelberg, 2000)
//    INPUT FILES:
//             None.
//    OUTPUT FILES:
//             None.

//{Include files}
#pragma once
#ifndef _STEP8_WOUND_MECHANICS_
#define _STPE8_WOUND_MECHANICS_
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

                 //File for vector-valued finite elements.
#include <deal.II/fe/fe_system.h>
				 // We will compose the vector-valued
				 // finite elements from regular Q1
				 // elements which can be found here,
				 // as usual:
#include <deal.II/fe/fe_q.h>
				 // This again is C++:
#include <fstream>
#include <iostream>  
#include "nr.h"

using namespace dealii;




//{The ElasticProblem class template}

				 // The main class is, except for its
				 // name, almost unchanged with
				 // respect to the step-6 example.
				 //
				 // The only change is the use of a
				 // different class for the <code>fe</code>
				 // variable: Instead of a concrete
				 // finite element class such as
				 // <code>FE_Q</code>, we now use a more
				 // generic one, <code>FESystem</code>. In
				 // fact, <code>FESystem</code> is not really a
				 // finite element itself in that it
				 // does not implement shape functions
				 // of its own.  Rather, it is a class
				 // that can be used to stack several
				 // other elements together to form
				 // one vector-valued finite
				 // element. In our case, we will
				 // compose the vector-valued element
				 // of <code>FE_Q(1)</code> objects, as shown
				 // below in the constructor of this
				 // class.
template <int dim>
class ElasticProblem {
public:
    ElasticProblem (int xstep, int ystep);
    ~ElasticProblem ();
    void run (const Mat_DP&, const Mat_DP&, const Mat_DP&);
    void run_until_converge(const Mat_DP&, const Mat_DP&, const Mat_DP&);
public:
    Mat_DP tissue_displacement_x;
    Mat_DP tissue_displacement_y;    

private:
    void setup_system ();
    void assemble_system (const Mat_DP&, const Mat_DP&, const Mat_DP&);
    void make_constraints();
    void solve ();
    void refine_grid ();
    void output_deformation_profile ();
    void output_woundcontour();

    void elasticity_coefficient(Mat_DP F, Vec_DP M, DP c_density, Mat4D_DP& A);

    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;

    FESystem<dim>        fe;

    ConstraintMatrix     hanging_node_constraints;
    ConstraintMatrix     constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       solution_pre;
    Vector<double>       system_rhs;
    
    int first_iteration_tag;//state variable
    DP residue_local;//state variable
    DP residue_total;//state variable
    DP residue_max;//state variable
    int mXstep, mYstep;

    std::map<std::pair<int, int>, int> mWoundContourHistory;
};


//{Right hand side values for the force balance equation}

				 // Before going over to the
				 // implementation of the main class,
				 // we declare and define the class
				 // which describes the right hand
				 // side. This time, the right hand
				 // side is vector-valued, as is the
				 // solution, so we will describe the
				 // changes required for this in some
				 // more detail.
				 //
				 // The first thing is that
				 // vector-valued functions have to
				 // have a constructor, since they
				 // need to pass down to the base
				 // class of how many components the
				 // function consists. The default
				 // value in the constructor of the
				 // base class is one (i.e.: a scalar
				 // function), which is why we did not
				 // need not define a constructor for
				 // the scalar function used in
				 // previous programs.
template <int dim>
class RightHandSide :  public Function<dim> 
{
public:
    RightHandSide ();
    
    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &values) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
				    std::vector<Vector<double> >   &value_list) const;
				    
};

//{collagen fiber direction}
template <int dim>
class Direction :  public Function<dim> 
{
public:
    Direction ();
    Direction (const Mat_DP& collagen);
    
    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &values) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
				    std::vector<Vector<double> >   &value_list) const;
protected:
    Mat_DP collagen;  
};


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
//{collagen fiber density}
template <int dim>
class Cdensity :  public Function<dim> 
{
public:
    Cdensity ();
    Cdensity (const Mat_DP& collagen_density_);
    
    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &values) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
				    std::vector<Vector<double> >   &value_list) const;
				    
protected:
    Mat_DP collagen_density;
};

#endif
