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
#ifndef _STEP8_CC_
#define _STEP8_CC_
#include "step8_wound_mechanics.h"
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/grid_refinement.h>
#include <fstream>
#include <iostream>  
#include "nr.h"
#include "BMP.h"
#include "util.h"

using namespace dealii;
extern const int xstep;
extern const int ystep;

extern Mat_DP *fibroblast_density_ptr;
extern Mat_DP *speedfield_ptr;
extern Mat_DP *thetafield_ptr;


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
  
  values(0) = value0*50 + 3*fb_density(p(1),p(0))*speedfield[(int)p(1)][(int)p(0)]/15*cos(thetafield[(int)p(1)][(int)p(0)]);
  values(1) = value1*50 + 3*fb_density(p(1),p(0))*speedfield[(int)p(1)][(int)p(0)]/15*sin(thetafield[(int)p(1)][(int)p(0)]);

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

//================================================================

template <int dim>
Direction<dim>::Direction ()
		:
		Function<dim> (dim)
{}

template <int dim>
Direction<dim>::Direction (const Mat_DP& collagen_)
        :collagen(collagen_),
         Function<dim> (dim)
{}


template <int dim>
inline
void Direction<dim>::vector_value (const Point<dim> &p,
				       Vector<double>   &values) const 
{
  Assert (values.size() == dim, 
	  ExcDimensionMismatch (values.size(), dim));
  Assert (dim >= 2, ExcNotImplemented());
    
  double collagenv = collagen[(int)(p(1)*0.1)][(int)(p(0)*0.1)];
    
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

  for (unsigned int p=0; p<n_points; ++p)
    Direction<dim>::vector_value (points[p],
				      value_list[p]);
}

/////////////////////////////////////////////////////////////////////////////////////////

template <int dim>
Cdensity<dim>::Cdensity ()
		:
		Function<dim> (dim)
{}

template <int dim>
Cdensity<dim>::Cdensity (const Mat_DP& collagen_density_)
        :collagen_density(collagen_density_),
         Function<dim> (dim)
{}


template <int dim>
inline
void Cdensity<dim>::vector_value (const Point<dim> &p,
				       Vector<double>   &values) const 
{
  Assert (values.size() == dim, 
	  ExcDimensionMismatch (values.size(), dim));
  Assert (dim >= 2, ExcNotImplemented());
    
  values(0) = collagen_density[(int)(p(1)/5.0)][(int)(p(0)/5.0)];
  
}


template <int dim>
void Cdensity<dim>::vector_value_list (const std::vector<Point<dim> > &points,
					    std::vector<Vector<double> >   &value_list) const 
{
  Assert (value_list.size() == points.size(), 
	  ExcDimensionMismatch (value_list.size(), points.size()));

  const unsigned int n_points = points.size();

  for (unsigned int p=0; p<n_points; ++p)
    Cdensity<dim>::vector_value (points[p],
				      value_list[p]);
}

////////////////////////////////////////////////////////////////////////////////////////////



template <int dim>
ElasticProblem<dim>::ElasticProblem (int xstep, int ystep)
		:
                mXstep(xstep), mYstep(ystep),
		dof_handler (triangulation),
		fe (FE_Q<dim>(1), dim)
{}
				 // In fact, the <code>FESystem</code> class
				 // has several more constructors
				 // which can perform more complex
				 // operations than just stacking
				 // together several scalar finite
				 // elements of the same type into
				 // one; we will get to know these
				 // possibilities in later examples.



				 // The destructor, on the other hand,
				 // is exactly as in step-6:
template <int dim>
ElasticProblem<dim>::~ElasticProblem () 
{
  dof_handler.clear ();
}


//{ElasticProblem::setup_system}

				 // Setting up the system of equations
				 // is identitical to the function
				 // used in the step-6 example. The
				 // <code>DoFHandler</code> class and all other
				 // classes used here are fully aware
				 // that the finite element we want to
				 // use is vector-valued, and take
				 // care of the vector-valuedness of
				 // the finite element themselves. (In
				 // fact, they do not, but this does
				 // not need to bother you: since they
				 // only need to know how many degrees
				 // of freedom there are per vertex,
				 // line and cell, and they do not ask
				 // what they represent, i.e. whether
				 // the finite element under
				 // consideration is vector-valued or
				 // whether it is, for example, a
				 // scalar Hermite element with
				 // several degrees of freedom on each
				 // vertex).
template <int dim>
void ElasticProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  hanging_node_constraints.close ();
  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

  hanging_node_constraints.condense (sparsity_pattern);

  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  solution_pre.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());

  tissue_displacement_x = Mat_DP(mYstep, mXstep);
  tissue_displacement_y = Mat_DP(mYstep, mXstep);
}


//{ElasticProblem::assemble_system}

				 // The big changes in this program
				 // are in the creation of matrix and
				 // right hand side, since they are
				 // problem-dependent. We will go
				 // through that process step-by-step,
				 // since it is a bit more complicated
				 // than in previous examples.
				 //
				 // The first parts of this function
				 // are the same as before, however:
				 // setting up a suitable quadrature
				 // formula, initializing an
				 // <code>FEValues</code> object for the
				 // (vector-valued) finite element we
				 // use as well as the quadrature
				 // object, and declaring a number of
				 // auxiliary arrays. In addition, we
				 // declare the ever same two
				 // abbreviations: <code>n_q_points</code> and
				 // <code>dofs_per_cell</code>. The number of
				 // degrees of freedom per cell we now
				 // obviously ask from the composed
				 // finite element rather than from
				 // the underlying scalar Q1
				 // element. Here, it is <code>dim</code> times
				 // the number of degrees of freedom
				 // per cell of the Q1 element, though
				 // this is not explicit knowledge we
				 // need to care about:
template <int dim>
void ElasticProblem<dim>::assemble_system (const Mat_DP& collagen, const Mat_DP& collagen_density, const Mat_DP& fibroblast_density) 
{  
  QGauss<dim>  quadrature_formula(2);

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values   | update_gradients | update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

				   // As was shown in previous
				   // examples as well, we need a
				   // place where to store the values
				   // of the coefficients at all the
				   // quadrature points on a cell. In
				   // the present situation, we have
				   // two coefficients, lambda and mu.
  std::vector<double>     lambda_values (n_q_points);
  std::vector<double>     mu_values (n_q_points);

				   // Well, we could as well have
				   // omitted the above two arrays
				   // since we will use constant
				   // coefficients for both lambda and
				   // mu, which can be declared like
				   // this. They both represent
				   // functions always returning the
				   // constant value 1.0. Although we
				   // could omit the respective
				   // factors in the assemblage of the
				   // matrix, we use them here for
				   // purpose of demonstration.
  ConstantFunction<dim> lambda(1.), mu(1.);

				   // Then again, we need to have the
				   // same for the right hand
				   // side. This is exactly as before
				   // in previous examples. However,
				   // we now have a vector-valued
				   // right hand side, which is why
				   // the data type of the
				   // <code>rhs_values</code> array is
				   // changed. We initialize it by
				   // <code>n_q_points</code> elements, each of
				   // which is a <code>Vector@<double@></code>
				   // with <code>dim</code> elements.
  RightHandSide<dim>      right_hand_side;
  std::vector<Vector<double> > rhs_values (n_q_points, Vector<double>(dim));

  Direction<dim>      direction(collagen);
  std::vector<Vector<double> > dir_values (n_q_points, Vector<double>(dim));
  
  Cdensity<dim>       cdensity(collagen_density);
  std::vector<Vector<double> > cde_values (n_q_points, Vector<double>(dim));					   
  

				   // Now we can begin with the loop
				   // over all cells:
  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						 endc = dof_handler.end();

  if(first_iteration_tag == 3){
      residue_total = 0;
      residue_max = 0;
  }

  int counter = 0;						 
  for (; cell!=endc; ++cell){  
      counter++;
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);
      
				       // Next we get the values of
				       // the coefficients at the
				       // quadrature points. Likewise
				       // for the right hand side:
      lambda.value_list (fe_values.get_quadrature_points(), lambda_values);
      mu.value_list     (fe_values.get_quadrature_points(), mu_values);

      right_hand_side.vector_value_list (fe_values.get_quadrature_points(), rhs_values);
      cdensity.vector_value_list (fe_values.get_quadrature_points(), cde_values); 
      direction.vector_value_list (fe_values.get_quadrature_points(), dir_values);

				       // Then assemble the entries of
				       // the local stiffness matrix
				       // and right hand side
				       // vector. This follows almost
				       // one-to-one the pattern
				       // described in the
				       // introduction of this
				       // example.  One of the few
				       // comments in place is that we
				       // can compute the number
				       // <code>comp(i)</code>, i.e. the index
				       // of the only nonzero vector
				       // component of shape function
				       // <code>i</code> using the
				       // <code>fe.system_to_component_index(i).first</code>
				       // function call below.
				       //
				       // (By accessing the
				       // <code>first</code> variable of
				       // the return value of the
				       // <code>system_to_component_index</code>
				       // function, you might
				       // already have guessed
				       // that there is more in
				       // it. In fact, the
				       // function returns a
				       // <code>std::pair@<unsigned int,
				       // unsigned int@></code>, of
				       // which the first element
				       // is <code>comp(i)</code> and the
				       // second is the value
				       // <code>base(i)</code> also noted
				       // in the introduction, i.e.
				       // the index
				       // of this shape function
				       // within all the shape
				       // functions that are nonzero
				       // in this component,
				       // i.e. <code>base(i)</code> in the
				       // diction of the
				       // introduction. This is not a
				       // number that we are usually
				       // interested in, however.)
				       //
				       // With this knowledge, we can
				       // assemble the local matrix
				       // contributions:                                        
 
      for (unsigned int q_point=0; q_point<n_q_points;++q_point){         
          Mat_DP F(2,2);
          Mat_DP F_pre(2,2);
          F[0][0]=1;
          F[0][1]=0;
          F[1][0]=0;
          F[1][1]=1;
          F_pre[0][0] = 1;
          F_pre[0][1] = 0;
          F_pre[1][0] = 0;
          F_pre[1][1] = 1;
          for (unsigned int ii=0; ii<dofs_per_cell; ++ii){
              //cout<<"solution_local="<<solution(local_dof_indices[ii])<<endl;
              const unsigned int component_ii = fe.system_to_component_index(ii).first; 
              if(first_iteration_tag == 1){//if this is first iteration testing convergence
//                  F[0][0] += 0.5;
//                  F[0][1] += 0;
//                  F[1][0] += 0;
//                  F[1][1] += 0.5;
//                  Do nothing. F is initiated to (1,0,0,1);   
              }
              else if(first_iteration_tag == 2){
                   F[0][0] += 
//                (fe_values.shape_grad(1,q_point)[component_1]*solution(local_to_global(1)) 
//                +fe_values.shape_grad(3,q_point)[component_3]*solution(local_to_global(3))
//                +fe_values.shape_grad(5,q_point)[component_5]*solution(local_to_global(5))
//                +fe_values.shape_grad(7,q_point)[component_7]*solution(local_to_global(7));local_dof_indices[i]
                ((ii%2==1) ? fe_values.shape_grad(ii,q_point)[component_ii]*solution(local_dof_indices[ii]) : 0);
                   F[0][1] += 
//                (fe_values.shape_grad(1,q_point)[component_2]*solution(local_to_global(1)) 
//                +fe_values.shape_grad(3,q_point)[component_4]*solution(local_to_global(3))
//                +fe_values.shape_grad(5,q_point)[component_6]*solution(local_to_global(5))
//                +fe_values.shape_grad(7,q_point)[component_8]*solution(local_to_global(7));
                ((ii%2==1) ? fe_values.shape_grad(ii,q_point)[1-component_ii]*solution(local_dof_indices[ii]) : 0);
                   F[1][0] += 
//                (fe_values.shape_grad(2,q_point)[component_1]*solution(local_to_global(2)) 
//                +fe_values.shape_grad(4,q_point)[component_3]*solution(local_to_global(4))
//                +fe_values.shape_grad(6,q_point)[component_5]*solution(local_to_global(6))
//                +fe_values.shape_grad(8,q_point)[component_7]*solution(local_to_global(8));
                ((ii%2==0) ? fe_values.shape_grad(ii,q_point)[1-component_ii]*solution(local_dof_indices[ii]) : 0);
                   F[1][1] +=  
//                (fe_values.shape_grad(2,q_point)[component_2]*solution(local_to_global(2)) 
//                +fe_values.shape_grad(4,q_point)[component_4]*solution(local_to_global(4))
//                +fe_values.shape_grad(6,q_point)[component_6]*solution(local_to_global(6))
//                +fe_values.shape_grad(8,q_point)[component_8]*solution(local_to_global(8));
                ((ii%2==0) ? fe_values.shape_grad(ii,q_point)[component_ii]*solution(local_dof_indices[ii]) : 0);                           
              }
              else {
                  F[0][0] += 
                ((ii%2==1) ? fe_values.shape_grad(ii,q_point)[component_ii]*solution(local_dof_indices[ii]) : 0);
                  F_pre[0][0] +=
                ((ii%2==1) ? fe_values.shape_grad(ii,q_point)[component_ii]*solution_pre(local_dof_indices[ii]) : 0);
              //dX1/dx1                
                  F[0][1] += 
                ((ii%2==1) ? fe_values.shape_grad(ii,q_point)[1-component_ii]*solution(local_dof_indices[ii]) : 0);
                  F_pre[0][1] +=
                ((ii%2==1) ? fe_values.shape_grad(ii,q_point)[1-component_ii]*solution_pre(local_dof_indices[ii]) : 0);                          
              //dX1/dx2                    
                  F[1][0] += 
                ((ii%2==0) ? fe_values.shape_grad(ii,q_point)[1-component_ii]*solution(local_dof_indices[ii]) : 0);
                  F_pre[1][0] +=
                ((ii%2==0) ? fe_values.shape_grad(ii,q_point)[1-component_ii]*solution_pre(local_dof_indices[ii]) : 0);
              //dX2/dx1     
                  F[1][1] +=  
                ((ii%2==0) ? fe_values.shape_grad(ii,q_point)[component_ii]*solution(local_dof_indices[ii]) : 0);
                  F_pre[1][1] +=
                ((ii%2==0) ? fe_values.shape_grad(ii,q_point)[component_ii]*solution_pre(local_dof_indices[ii]) : 0);
              //dX2/dx2
              }
          }
           

        
          if(first_iteration_tag == 3){
              residue_local = sqrt(pow(F[0][0] - F_pre[0][0],2) + pow(F[0][1] - F_pre[0][1],2) 
                      + pow(F[1][0] - F_pre[1][0],2) + pow(F[1][1] - F_pre[1][1],2)) ;
              residue_total += residue_local;
              if(residue_local > residue_max) residue_max = residue_local;
              if(counter == 1000){
                  cout << "F_pre[0][0]=" << F_pre[0][0] << " F_pre[0][1] =" << F_pre[0][1] 
                       << " F_pre[1][0]=" << F_pre[1][0] << " F_pre[1][1] =" << F_pre[1][1] 
                       << "solution_pre=" << solution_pre(local_dof_indices[1]) << endl;
                  cout << "F[0][0]=" << F[0][0]     << " F[0][1]=" << F[0][1]     
                       << " F[1][0]=" << F[1][0]     << " F[1][1]=" << F[1][1]     
                       << "solution   =" << solution(local_dof_indices[1]) << endl;
              }
              F[0][0] = (F[0][0]+F_pre[0][0])/2;
              F[0][1] = (F[0][1]+F_pre[0][1])/2;
              F[1][0] = (F[1][0]+F_pre[1][0])/2;
              F[1][1] = (F[1][1]+F_pre[1][1])/2;
          }
          
           
          Vec_DP M(2); 
          M[0] = dir_values[q_point](0); 
          M[1] = dir_values[q_point](1);
          
          DP c_density = cde_values[q_point](0);
       
          Mat4D_DP A(2,2,2,2); 
/*
//isotropic
          for(int ii=0;ii<2;ii++)
              for(int jj=0;jj<2;jj++)
                  for(int aa=0;aa<2;aa++)
                      for(int bb=0;bb<2;bb++){
                          A[ii][jj][aa][bb] = ((ii==jj && aa==bb) ? lambda_values[q_point] : 0)
                                            + ((ii==aa && jj==bb) ? mu_values[q_point] : 0)
                                            + ((ii==bb && jj==aa) ? mu_values[q_point] : 0);
                      }
*/

/*
lambda+2mu  0  0  lambda
0           mu mu      0
0           mu mu      0
lambda      0  0  lambda+2mu
*/

/*
//transversely isotropic linear elasticity
          for(int ii=0;ii<2;ii++)
              for(int jj=0;jj<2;jj++)
                  for(int aa=0;aa<2;aa++)
                      for(int bb=0;bb<2;bb++){
                          A[ii][jj][aa][bb] = ((ii==jj && aa==bb) ? c12 : 0)
                                            + c66*(((ii==aa && jj==bb) ? 1 : 0) + (ii==bb && jj==aa) ? 1 : 0))
                                            + (c13-c12)*(((ii==jj) ? M[aa]*M[bb] : 0) + (aa==bb) ? M[ii]*M[jj] : 0))
                                            + (c44-c66)*(((ii==aa) ? M[jj]*M[bb] : 0) + (ii==bb) ? M[jj]*m[aa] : 0))+
                                                         ((jj==aa) ? M[ii]*M[bb] : 0) + (jj==bb) ? M[ii]*M[aa] : 0)))
                                            + (c11 + c33 - 2*c13 - 4*c44) * M[ii]*M[jj]*M[aa]*M[ll];
                      }                              
*/        
          elasticity_coefficient(F, M, c_density, A);
          
          if(counter == 100||counter == 200){
              cout<<"cell="<<cell<<" q_point="<<q_point<<endl;
              cout<<"A0000 = "<<A[0][0][0][0]<<" A0001 = "<<A[0][0][0][1]<<" A0010 = "<<A[0][0][1][0]<<" A0011 = "<<A[0][0][1][1]<<endl;
              cout<<"A0100 = "<<A[0][1][0][0]<<" A0101 = "<<A[0][1][0][1]<<" A0110 = "<<A[0][1][1][0]<<" A0111 = "<<A[0][1][1][1]<<endl;
              cout<<"A1000 = "<<A[1][0][0][0]<<" A1001 = "<<A[1][0][0][1]<<" A1010 = "<<A[1][0][1][0]<<" A1011 = "<<A[1][0][1][1]<<endl;
              cout<<"A1100 = "<<A[1][1][0][0]<<" A1101 = "<<A[1][1][0][1]<<" A1110 = "<<A[1][1][1][0]<<" A1111 = "<<A[1][1][1][1]<<endl;
          }
                       
          for (unsigned int i=0; i<dofs_per_cell; ++i){
	          const unsigned int component_i = fe.system_to_component_index(i).first;
/*	      cout<<"component_"<<i<<"="<<component_i<<".....................";
              cout<<"q_point="<<q_point<<"..."<<fe_values.shape_value(i,q_point)<<"	";
              cout<<fe_values.shape_grad(i,q_point)<<endl;   */
              //ËùÓÐfe_values.shape_value ºÍËùÓÐfe_values.shape_grad¶¼ÊÇÒ»ÑùµÄ£¬ÒòÎª°Ñ¾ßÌåÊýÖµ¶¼³é³öÈ¥ÁË°¡°¡°¡               
	      
              for (unsigned int j=0; j<dofs_per_cell; ++j){
              
                  const unsigned int component_j = fe.system_to_component_index(j).first;
                  cell_matrix(i,j) +=  //i=1~8;j=1~8
		                           //shape_value(1)=(1,0)*shape_base1;
		                           //shape_value(2)=(0,1)*shape_base1;
		                           //shape_value(3)=(1,0)*shape_base2;
		                           //shape_value(4)=(0,1)*shape_base2;
		                           //shape_value(5)=(1,0)*shape_base3;
		                           //shape_value(6)=(0,1)*shape_base3;
		                           //shape_value(7)=(1,0)*shape_base4;
		                           //shape_value(8)=(0,1)*shape_base4;                                    		                         
                  (
//                    (fe_values.shape_grad(1,q_point)[component_1==1] +
//        		       fe_values.shape_grad(3,q_point)[component_3==1] +
//                             fe_values.shape_grad(5,q_point)[component_5==1] +
//        		       fe_values.shape_grad(7,q_point)[component_7==1] +//dx1/dX1
//        		       )*A1i11+
//        		       
//                    (fe_values.shape_grad(1,q_point)[component_2==2] +
//        		       fe_values.shape_grad(3,q_point)[component_4==2] +
//                             fe_values.shape_grad(5,q_point)[component_6==2] +
//        		       fe_values.shape_grad(7,q_point)[component_8==2] +//dx1/dX2
//                             )*A1i21+
//        		           		       
//    	              (fe_values.shape_grad(2,q_point)[component_1==1] +
//        		       fe_values.shape_grad(4,q_point)[component_3==1] +
//        		       fe_values.shape_grad(6,q_point)[component_5==1] +
//        		       fe_values.shape_grad(8,q_point)[component_7==1] +//dx2/dX1
//        		       )*A1i12+
//        		       
//       	          (fe_values.shape_grad(2,q_point)[component_2==2] +
//        		       fe_values.shape_grad(4,q_point)[component_4==2] +
//        		       fe_values.shape_grad(6,q_point)[component_6==2] +
//        		       fe_values.shape_grad(8,q_point)[component_8==2] +//dx2/dX2
//                     )*A1i22 
                    (
                        fe_values.shape_grad(j,q_point)[component_j]
                        *((j%2==0)?A[0][0][0][0]:A[0][0][1][1]) +
                        fe_values.shape_grad(j,q_point)[1-component_j]
                        *((j%2==0)?A[0][0][1][0]:A[0][0][0][1])
                    )//A1i11*dx1/dX1+A1i21*dx1/dX2+A1i12*dx2/dX1+A1i22*dx2/dX2(i=1)
                    *((i%2==0)?fe_values.shape_grad(i,q_point)[component_i]:0)//dv1/dX1
    		       
                  +
                    (
                        fe_values.shape_grad(j,q_point)[component_j]
                        *((j%2==0)?A[1][0][0][0]:A[1][0][1][1]) +
                        fe_values.shape_grad(j,q_point)[1-component_j]
                        *((j%2==0)?A[1][0][1][0]:A[1][0][0][1])
                    )//A2i11*dx1/dX1+A2i21*dx1/dX2+A2i12*dx2/dX1+A2i22*dx2/dX2(i=1)
                    *((i%2==0)?fe_values.shape_grad(i,q_point)[1-component_i]:0)//dv1/dX2
    		       
    		      +
    		        (
                        fe_values.shape_grad(j,q_point)[component_j]
                        *((j%2==0)?A[0][1][0][0]:A[0][1][1][1]) +
                        fe_values.shape_grad(j,q_point)[1-component_j]
                        *((j%2==0)?A[0][1][1][0]:A[0][1][0][1])
                    )//A1i11*dx1/dX1+A1i21*dx1/dX2+A1i12*dx2/dX1+A1i22*dx2/dX2(i=2)
                    *((i%2==1)?fe_values.shape_grad(i,q_point)[1-component_i]:0)//dv2/dX1
   		          +
   		            (
                        fe_values.shape_grad(j,q_point)[component_j]
                        *((j%2==0)?A[1][1][0][0]:A[1][1][1][1]) +
                        fe_values.shape_grad(j,q_point)[1-component_j]
                        *((j%2==0)?A[1][1][1][0]:A[1][1][0][1])
                    )//A2i11*dx1/dX1+A2i21*dx1/dX2+A2i12*dx2/dX1+A2i22*dx2/dX2(i=2)
                    *((i%2==1)?fe_values.shape_grad(i,q_point)[component_i]:0)//dv2/
                  +
/*********************************tethering force********************************************************/                 
                    (
                        ((i%2==j%2)?0.0001*fe_values.shape_value(j,q_point)*fe_values.shape_value(i,q_point):0)
                    )//(X*V) == x1*v1+x2*v2
/********************************************************************************************************/ 
                  )*fe_values.JxW(q_point); 
              }
          }
      }

				       // Assembling the right hand
				       // side is also just as
				       // discussed in the
				       // introduction:
      for (unsigned int i=0; i<dofs_per_cell; ++i){
	      const unsigned int component_i = fe.system_to_component_index(i).first;
	  
	      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              cell_rhs(i) += fe_values.shape_value(i,q_point) *
			                 rhs_values[q_point](component_i) *
			                 fe_values.JxW(q_point);
      }

				       // The transfer from local
				       // degrees of freedom into the
				       // global matrix and right hand
				       // side vector does not depend
				       // on the equation under
				       // consideration, and is thus
				       // the same as in all previous
				       // examples. The same holds for
				       // the elimination of hanging
				       // nodes from the matrix and
				       // right hand side, once we are
				       // done with assembling the
				       // entire linear system:
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i){
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
	          system_matrix.add (local_dof_indices[i],local_dof_indices[j],cell_matrix(i,j));
	      system_rhs(local_dof_indices[i]) += cell_rhs(i);
	  }
  }//end of cell loop
  cout << "total # of cell is:" << counter << endl;
  cout << "residue_total=" << residue_total <<endl;
  cout << "residue_max=" << residue_max <<endl;
  
  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);

				   // The interpolation of the
				   // boundary values needs a small
				   // modification: since the solution
				   // function is vector-valued, so
				   // need to be the boundary
				   // values. The <code>ZeroFunction</code>
				   // constructor accepts a parameter
				   // that tells it that it shall
				   // represent a vector valued,
				   // constant zero function with that
				   // many components. By default,
				   // this parameter is equal to one,
				   // in which case the
				   // <code>ZeroFunction</code> object would
				   // represent a scalar
				   // function. Since the solution
				   // vector has <code>dim</code> components,
				   // we need to pass <code>dim</code> as
				   // number of components to the zero
				   // function as well.
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(dim),
					    boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs);
}



//{ElasticProblem::solve}

				 // The solver does not care about
				 // where the system of equations
				 // comes, as long as it stays
				 // positive definite and symmetric
				 // (which are the requirements for
				 // the use of the CG solver), which
				 // the system indeed is. Therefore,
				 // we need not change anything.
template <int dim>
void ElasticProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);
  
  cout << "in solve(),system_matrix(100,100)=" << system_matrix(100,100) << endl;
  
  cout << "solution pre=" << solution(1000);
  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);
  
  hanging_node_constraints.distribute (solution);
  cout << "  solution aft=" << solution(1000) << endl;
}


 
template <int dim>
void ElasticProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  typename FunctionMap<dim>::type neumann_boundary;
  KellyErrorEstimator<dim>::estimate (dof_handler,
				      QGauss<dim-1>(2),
				      neumann_boundary,
				      solution,
				      estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.3, 0.03);

  triangulation.execute_coarsening_and_refinement ();
}



//{ElasticProblem::output_results}

				 // The output happens mostly as has
				 // been shown in previous examples
				 // already. The only difference is
				 // that the solution function is
				 // vector valued. The <code>DataOut</code>
				 // class takes care of this
				 // automatically, but we have to give
				 // each component of the solution
				 // vector a different name.
template <int dim>
void ElasticProblem<dim>::output_deformation_profile () 
{
  static int out_i = 0;
  char filename[16];
  sprintf(filename, "output/solution-%d.vtk", out_i);
 
  cout << "output deformation profile..." << endl; 
  std::ofstream output (filename);

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
 
  std::vector<std::string> solution_names;
  switch (dim)
    {
      case 1:
	    solution_names.push_back ("displacement");
	    break;
      case 2:
	    solution_names.push_back ("x_displacement");	    
	    solution_names.push_back ("y_displacement");
	    std::cout<<"2D output"<<std::endl;
	    break;
      case 3:
	    solution_names.push_back ("x_displacement");	    
	    solution_names.push_back ("y_displacement");
	    solution_names.push_back ("z_displacement");
	    std::cout<<"3D output"<<std::endl;
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    }
	     
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(dim,DataComponentInterpretation::component_is_part_of_vector);
    //data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  if(data_component_interpretation.size() == 2)std::cout<<"data_dim=2"<<std::endl; 
  else std::cout<<"error in data dim"<<std::endl; 			    
  data_out.add_data_vector (solution, solution_names,
                            DataOut<dim>::type_dof_data,
                            data_component_interpretation);
  data_out.build_patches ();
  data_out.write_vtk (output);

  Point<dim> point; 
    //std::vector<Vector<double> > value_list (1000, Vector<double>(dim));
  Vector<double> value;
  for(int i=0;i<mYstep;i++){            
      for(int j=0;j<mXstep;j++){ 
          point(0)=(double)j; 
          point(1)=(double)i;
          if(isOnWoundEdge(j, i, mXstep, mYstep)) {
                VectorTools::point_value(dof_handler, solution, point, value);
                //ºÜÏÔÈ»£¬ÕâÀïµÄdof_handlerÀàÐÍ²»ÊÇ DoFHandler< dim, spacedim >,µ«Ò²¿ÉÒÔÕâÑùÓÃ
                //solution ÊÇ vector-valued FE function
                //ËùÒÔÕâÑùÌáÈ¡½á¹ûÊÇÍêÈ«ÕýÈ·µÄ 
                tissue_displacement_x[i][j] = value(0);
                tissue_displacement_y[i][j] = value(1);
                //cout<<"value(0)"<<value(0)<<"value(1)"<<value(1)<<endl; 
          }
      }  
  }

/*
template<int dim, class InVector , int spacedim> 
static void VectorTools::point_value  ( 
  const DoFHandler< dim, spacedim > &  dof,  
  const InVector &  fe_function,  
  const Point< spacedim > &  point,  
  Vector< double > &  value  
 ) 
*/
//see reference to step-7 for the use of VectorTools::___
  out_i++;
}

template <int dim>
void ElasticProblem<dim>::run_until_converge (const Mat_DP& collagen, const Mat_DP& collagen_density, const Mat_DP& fibroblast_density) 
{
    first_iteration_tag = 1;
    residue_total = 100000;
    int iter_counter = 0;
    int continue_tag;
    do{ //condition of convergence????????????????
        cout <<endl << "iteration:"<<first_iteration_tag<<endl;
        run(collagen, collagen_density, fibroblast_density); 
                              
        if(first_iteration_tag < 3){
            first_iteration_tag++; 
            residue_total = 100000;
        }else if(first_iteration_tag == 3){
            //first_iteration_tag = 1 : don't have solution
            //first_iteration_tag = 2 : have solution; don't have solution_pre
            //first_iteration_tag = 3 : have solution and solution_pre 
        }

        iter_counter++;
                    
        continue_tag = 1;
        if(iter_counter > 20 && iter_counter%10 == 0){
            printf("Continues?");
            scanf("%d", &continue_tag);
        }
        if(continue_tag == 0)break;     
    }while(residue_total > 20.0);

    output_deformation_profile();
    cout << "output wound contour...\n" << endl;
    output_woundcontour();
}

template <int dim>
void ElasticProblem<dim>::run(const Mat_DP& collagen, const Mat_DP& collagen_density, const Mat_DP& fibroblast_density) {
    static int cycle=0;
    if(cycle == 0){
        if (mXstep == mYstep){ 
            GridGenerator::hyper_cube (triangulation, 0, mXstep);
        } else {
            Point<dim> p0(0,0);
            Point<dim> p1(mXstep, mYstep);
            GridGenerator::hyper_rectangle (triangulation, p0, p1);
        }
        triangulation.refine_global (7);    
        std::cout << "   Number of active cells:       "
	              << triangulation.n_active_cells()
	              << std::endl;   
        setup_system ();    
        std::cout << "   Number of degrees of freedom: "
	              << dof_handler.n_dofs()
	              << std::endl;
        cycle++;
    }
    
//*************************************************
    system_matrix.reinit (sparsity_pattern);
    system_rhs.reinit (dof_handler.n_dofs());
//*********crucial steps***************************
    assemble_system (collagen, collagen_density, fibroblast_density);//test convergence using solution(new) and solution_pre(old)
    //cout << "system_matrix(100,100)=" << system_matrix(100,100) << endl;
    //cout << "system_rhs(100)=" << system_rhs(100) << endl;
    solution_pre = solution;//store the old solution to solution_pre
    solve ();//get the new solution
}

template <int dim>
void ElasticProblem<dim>::output_woundcontour(){
    Mat_DP woundcontour(mYstep,mXstep);
    DP tdx,tdy;
    for(size_t i=0; i<mYstep; i++){
        for(size_t j=0; j<mXstep; j++){
            woundcontour[i][j] = 0;
        }
    }

    for(size_t i=0; i<mYstep; i++){
        for(size_t j=0; j<mXstep; j++){
            if( isOnWoundEdge(j,i,mXstep,mYstep) ){
                tdx = tissue_displacement_x[i][j];
                tdy = tissue_displacement_y[i][j];
                if(i+(int)tdy>=0 && i+(int)tdy<mYstep && j+(int)tdx>=0 && j+(int)tdy<mXstep){
                    woundcontour[i+(int)tdy][j+(int)tdx] = 1;
                }
            }
        }
    }

    static int out_i = 0;
    char file_name[21];
    sprintf(file_name, "output/cont%05d.BMP", out_i++);
    BMP::output_BMP(file_name, 14, woundcontour, mXstep, mYstep);

    return;
}

//   EquationData::Pressure<dim> pres (t_0);
//   VectorTools::interpolate (dof_handler_pressure, pres, pres_n_minus_1);
//   pres.advance_time (dt);
//   VectorTools::interpolate (dof_handler_pressure, pres, pres_n);

//calculate the elasticity coefficients
template <int dim>
void ElasticProblem<dim>::elasticity_coefficient(Mat_DP F, Vec_DP M, DP c_density, Mat4D_DP& A) {   
    
    //from tissue displacement, get tissue deformation gradients in two different directions

    //calculate Aalpha_i_beta_j
    //how we get Falpha_i? solution of the last step
    //how we access solution at certain spatial point? 		            
    //solution(local_dof_indices[i])
    //how we access the two entries of the solution...difficult, 
    //by interpolate I guess,but anyway,there will be a way we can access it.
    //cell_rhs(i) += fe_values.shape_value(i,q_point) 
    //rhs_values[q_point](component_i) *
    //fe_values.JxW(q_point);
    Vec_DP m(2);
    Vec_DP dW_dI(5+1);

    Mat_DP F_inverse(2,2),C(2,2),B(2,2);
    Mat_DP ddW_dIdI(5+1,5+1);      
    Mat_DP dI1_dF(2,2), dI2_dF(2,2), dI3_dF(2,2), dI4_dF(2,2), dI5_dF(2,2);
    
    Mat4D_DP ddI1_dFdF(2,2,2,2), ddI2_dFdF(2,2,2,2), ddI3_dFdF(2,2,2,2), ddI4_dFdF(2,2,2,2), ddI5_dFdF(2,2,2,2);

      
//    W=c/2*(I1-3) + k1/2k2*(exp[k2*(I4-1)^2]-1) + k1/2k2*(exp[k2*(I4-1)^2]-1)
//    exp(x) ~ 1+x    
//    W=c/2*(I1-3) + k1/2*(I4-1)^2
//    exp(x) ~ 1+x+x^2/2
//    W=c/2*(I1-3) + k1/2*(I4-1)^2 + k1*k2/4*(I4-1)^4 
    //B=F*F_trans
    B[0][0]=F[0][0]*F[0][0]+F[0][1]*F[0][1];
    B[0][1]=F[0][0]*F[1][0]+F[0][1]*F[1][1];
    B[1][0]=F[1][0]*F[0][0]+F[1][1]*F[0][1];
    B[1][1]=F[1][0]*F[1][0]+F[1][1]*F[1][1];
    
    //C=F_trans*F
    C[0][0]=F[0][0]*F[0][0]+F[1][0]*F[1][0];
    C[0][1]=F[0][0]*F[0][1]+F[1][0]*F[1][1];
    C[1][0]=F[0][1]*F[0][0]+F[1][1]*F[1][0];
    C[1][1]=F[0][1]*F[0][1]+F[1][1]*F[1][1]; 
    
    double I1,I2,I3,I4;
    //I1 = tr(B)
    I1 = B[0][0]+B[1][1];                     //cout << "I1=" <<I1<<endl;
    //I2 = det(B)
    I2 = B[0][0]*B[1][1]-B[0][1]*B[1][0];     //cout << "I2=" <<I2<<endl;
    //I4=M*(CM)=M*(C[0][0]*M[0]+C[0][1]*M[1],C[1][0]*M[0]+C[1][1]*M[1])
    I3 = 0;
    I4 = M[0]*(C[0][0]*M[0]+C[0][1]*M[1])+M[1]*(C[1][0]*M[0]+C[1][1]*M[1]);
    //I5=M*(C^2M)
    //cout <<"I4=" <<I4<<endl;
             
       
    int n,ii,jj,aa,bb,cc;
    double c=50, k1=9*c_density;
    double k2=0.2;
    {        
        dW_dI[1] = c/2;
        dW_dI[2] = 0;
        dW_dI[3] = 0;
        dW_dI[4] = k1*(I4-1)+k1*k2*pow(I4-1,3);//k1*(I4-1); //k1*(I4-1)*exp(k2*pow(I4-1,2));
        //W=c/2*(I1-3) + k1/2*(I4-1)^2 + k1*k2/4*(I4-1)^4
        dW_dI[5] = 0;    
    }

        
    for(n=1;n<=5;n++){           
        ddW_dIdI[1][n] = 0;
        ddW_dIdI[2][n] = 0;
        ddW_dIdI[3][n] = 0;
        ddW_dIdI[4][n] = 0;
        ddW_dIdI[5][n] = 0;        
    }
    ddW_dIdI[1][1] = 0;//0.1;
    ddW_dIdI[4][4] = k1 + k1*k2*3*pow(I4-1,2);//k1;//(k1+2*k1*k2*pow(I4-1,2))*exp(k2*pow(I4-1,2)); 
    //cout<<"ddW_dIdI[4][4]="<<ddW_dIdI[4][4]<<endl;
    
    //m=F*M;
    m[0]=F[0][0]*M[0]+F[0][1]*M[1];
    m[1]=F[1][0]*M[0]+F[1][1]*M[1];
    
    Mat_DP M_tensormul_m(2,2);
    M_tensormul_m[0][0] = M[0]*m[0];
    M_tensormul_m[0][1] = M[0]*m[1];
    M_tensormul_m[1][0] = M[1]*m[0];
    M_tensormul_m[1][1] = M[1]*m[1];

/*  
    cout<< "M*m00=" << M_tensormul_m[0][0] << endl;
    cout<< "M*m01=" << M_tensormul_m[0][1] << endl;
    cout<< "M*m10=" << M_tensormul_m[1][0] << endl;
    cout<< "M*m11=" << M_tensormul_m[1][1] << endl;
*/
       
    matrix_inverse(F_inverse,F);
      
    for(ii=0;ii<2;ii++)
        for(jj=0;jj<2;jj++){
            //±êÁ¿¶Ô¾ØÕóÇóµ¼²»×ªÖÃ 
            dI1_dF[ii][jj] = 2*F[jj][ii];
            dI2_dF[ii][jj] = 2*(I1*F[jj][ii]-F[0][ii]*B[0][jj]-F[1][ii]*B[1][jj]);
            dI3_dF[ii][jj] = 2*I2*F_inverse[ii][jj];
            dI4_dF[ii][jj] = 2*M_tensormul_m[ii][jj];
            //ÔõÃ´¸Ä³É1Ö®ºó³ÌÐòºÃÏñ¶ÔÁË£¿
            dI5_dF[ii][jj] = 0;
            //dI5_dF[ii][jj] = 2*(M_tensormul_Bm[ii][jj]+CM_tensormul_m[ii][jj]);
        }   

    for(ii=0;ii<2;ii++)
        for(aa=0;aa<2;aa++)      
            for(jj=0;jj<2;jj++)                
                for(bb=0;bb<2;bb++){                   
                    ddI1_dFdF[ii][aa][jj][bb] = 
                        2*((aa==bb)?1:0) * ((ii==jj)?1:0);              
                    ddI2_dFdF[ii][aa][jj][bb] = 
                        2*I1*((aa==bb)?1:0) * ((ii==jj)?1:0) 
                        + 4*F[ii][aa]*F[jj][bb]
                        - 2*F[ii][bb]*F[jj][aa]
                        - 2*C[aa][bb]*((ii==jj)?1:0)
                        - 2*((aa==bb)?1:0)*B[ii][jj];                    
                    ddI3_dFdF[ii][aa][jj][bb] = 
                        4*I3*F_inverse[aa][ii]*F_inverse[bb][jj]
                        - 2*I3*F_inverse[aa][jj]*F_inverse[bb][ii];
                    ddI4_dFdF[ii][aa][jj][bb] = 
                        2*M[aa]*M[bb]*((ii==jj)?1:0);
                    
                    ddI5_dFdF[ii][aa][jj][bb] = 0;
                    for(cc=0;cc<2;cc++){
                        ddI5_dFdF[ii][aa][jj][bb] += 
                            M[aa]*M[bb]*B[ii][jj] 
                            + 2*M[aa]*C[bb][cc]*M[cc]*((ii==jj)?1:0)
                            + 2*C[aa][cc]*M[cc]*M[bb]*((ii==jj)?1:0)
                            + 2*((aa==bb)?1:0)*m[ii]*m[jj]
                            + 2*M[aa]*F[ii][bb]*m[jj]
                            + 2*M[bb]*F[jj][aa]*m[ii];   
                    }
                }
                 
    int i,alpha,j,beta;
    for(i=0;i<2;i++)
        for(alpha=0;alpha<2;alpha++)
            for(j=0;j<2;j++)
                for(beta=0;beta<2;beta++){                                        
                    A[alpha][i][beta][j] = 0;
                    A[alpha][i][beta][j] += dW_dI[1]*ddI1_dFdF[i][alpha][j][beta];
                    A[alpha][i][beta][j] += dW_dI[2]*ddI2_dFdF[i][alpha][j][beta];
                    A[alpha][i][beta][j] += dW_dI[3]*ddI3_dFdF[i][alpha][j][beta];
                    A[alpha][i][beta][j] += dW_dI[4]*ddI4_dFdF[i][alpha][j][beta];
//                    A[alpha][i][beta][j] += dW_dI[5]*ddI5_dFdF[i][alpha][j][beta];

                    //if(output == 1)cout << "dW_dI4=" << dW_dI[4] << "ddI4_dFdF=" <<ddI4_dFdF[i][alpha][j][beta] << endl;                    

                    A[alpha][i][beta][j] += ddW_dIdI[1][1]*dI1_dF[i][alpha]*dI1_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[1][2]*dI1_dF[i][alpha]*dI2_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[1][3]*dI1_dF[i][alpha]*dI3_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[1][4]*dI1_dF[i][alpha]*dI4_dF[j][beta];
//                    A[alpha][i][beta][j] += ddW_dIdI[1][5]*dI1_dF[i][alpha]*dI5_dF[j][beta];
                    
                    A[alpha][i][beta][j] += ddW_dIdI[2][1]*dI2_dF[i][alpha]*dI1_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[2][2]*dI2_dF[i][alpha]*dI2_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[2][3]*dI2_dF[i][alpha]*dI3_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[2][4]*dI2_dF[i][alpha]*dI4_dF[j][beta];
//                    A[alpha][i][beta][j] += ddW_dIdI[2][5]*dI2_dF[i][alpha]*dI5_dF[j][beta];
                    
                    A[alpha][i][beta][j] += ddW_dIdI[3][1]*dI3_dF[i][alpha]*dI1_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[3][2]*dI3_dF[i][alpha]*dI2_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[3][3]*dI3_dF[i][alpha]*dI3_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[3][4]*dI3_dF[i][alpha]*dI4_dF[j][beta];
//                    A[alpha][i][beta][j] += ddW_dIdI[3][5]*dI3_dF[i][alpha]*dI5_dF[j][beta];
                                                
                    A[alpha][i][beta][j] += ddW_dIdI[4][1]*dI4_dF[i][alpha]*dI1_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[4][2]*dI4_dF[i][alpha]*dI2_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[4][3]*dI4_dF[i][alpha]*dI3_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[4][4]*dI4_dF[i][alpha]*dI4_dF[j][beta];
//                    A[alpha][i][beta][j] += ddW_dIdI[4][5]*dI4_dF[i][alpha]*dI5_dF[j][beta];
                    
                    A[alpha][i][beta][j] += ddW_dIdI[5][1]*dI5_dF[i][alpha]*dI1_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[5][2]*dI5_dF[i][alpha]*dI2_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[5][3]*dI5_dF[i][alpha]*dI3_dF[j][beta];
                    A[alpha][i][beta][j] += ddW_dIdI[5][4]*dI5_dF[i][alpha]*dI4_dF[j][beta];
//                    A[alpha][i][beta][j] += ddW_dIdI[5][5]*dI5_dF[i][alpha]*dI5_dF[j][beta]; 

                } 
    return;
}


#endif
