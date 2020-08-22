#ifndef _MATERIAL_TRANSISOTROPIC_HGO_HH_
#define _MATERIAL_TRANSISOTROPIC_HGO_HH_

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

 
/*Reference:
 * 1) A viscoelastic model for fiber-reinforced composites at finite
 * strains: Continuum basis, computational aspects and applications.
 * Holzapfel (2001)
 * 2) Hyperelastic modeling of arterial layers with distributed collagen
 * fibre orientations. Gasser, Odgen, Holzapfel (2006)
 */

#include "Material_Transverse_Isotropic.h"
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
#include "util.h"


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

// @sect3{Some standard tensors}
// Now we define some frequently used second and fourth-order tensors:

  template <int dim>
  const SymmetricTensor<2, dim>
  StandardTensors<dim>::I = unit_symmetric_tensor<dim>();

  template <int dim>
  const SymmetricTensor<4, dim>
  StandardTensors<dim>::IxI = outer_product(I, I);

  template <int dim>
  const SymmetricTensor<4, dim>
  StandardTensors<dim>::II = identity_tensor<dim>();

  template <int dim>
  const SymmetricTensor<4, dim>
  StandardTensors<dim>::dev_P = deviator_tensor<dim>();


// @sect3{Compressible neo-Hookean material within a three-field formulation}

// As discussed in the Introduction, Neo-Hookean materials are a type of
// hyperelastic materials.  The entire domain is assumed to be composed of a
// compressible neo-Hookean material.  This class defines the behaviour of
// this material within a three-field formulation.  Compressible neo-Hookean
// materials can be described by a strain-energy function (SEF) $ \Psi =
// \Psi_{\text{iso}}(\overline{\mathbf{b}}) + \Psi_{\text{vol}}(\widetilde{J})
// $.
//
// The isochoric response is given by $
// \Psi_{\text{iso}}(\overline{\mathbf{b}}) = c_{1} [\overline{I}_{1} - 3] $
// where $ c_{1} = \frac{\mu}{2} $ and $\overline{I}_{1}$ is the first
// invariant of the left- or right-isochoric Cauchy-Green deformation tensors.
// That is $\overline{I}_1 :=\textrm{tr}(\overline{\mathbf{b}})$.  In this
// example the SEF that governs the volumetric response is defined as $
// \Psi_{\text{vol}}(\widetilde{J}) = \kappa \frac{1}{4} [ \widetilde{J}^2 - 1
// - 2\textrm{ln}\; \widetilde{J} ]$,  where $\kappa:= \lambda + 2/3 \mu$ is
// the <a href="http://en.wikipedia.org/wiki/Bulk_modulus">bulk modulus</a>
// and $\lambda$ is <a
// href="http://en.wikipedia.org/wiki/Lam%C3%A9_parameters">Lame's first
// parameter</a>.
//
// The following class will be used to characterize the material we work with,
// and provides a central point that one would need to modify if one were to
// implement a different material model. For it to work, we will store one
// object of this type per quadrature point, and in each of these objects
// store the current state (characterized by the values or measures  of the three fields)
// so that we can compute the elastic coefficients linearized around the
// current state.
//std::pow is not optimized for small integer
template <int dim>
Material_Compressible_TransIsotropic_Three_Field<dim>::Material_Compressible_TransIsotropic_Three_Field(const double mu,
                                               const double nu,
                                               const double dispersion_coeff_,
                                               const double k1_,
                                               const double k2_,
                                               Tensor<1,dim> M_)
      :
      kappa((2.0 * mu * (1.0 + nu)) / (3.0 * (1.0 - 2.0 * nu))),
      c_1(mu / 2.0),
      dispersion_coeff(dispersion_coeff_),
      k1(k1_),
      k2(k2_),
      M(M_),
      det_Fe(1.0),
      p_tilde(0.0),
      J_tilde(1.0),
      b_bar(StandardTensors<dim>::I)
      {
          Assert(kappa > 0, ExcInternalError());
          Tensor<2, dim> C = StandardTensors<dim>::I;
          I1 = trace(C);
          Tensor<1, dim> tmp = C * M;
          I4 = M * tmp;
      }

    // We update the material model with various deformation dependent data
    // based on $F$ and the pressure $\widetilde{p}$ and dilatation
    // $\widetilde{J}$, and at the end of the function include a physical
    // check for internal consistency:
template <int dim>
void Material_Compressible_TransIsotropic_Three_Field<dim>::update_material_data(const Tensor<2, dim> &Fe_par,
                              const double p_tilde_in,
                              const double J_tilde_in)
{
      Fe = Fe_par;
      det_Fe = determinant(Fe_par);
      b_bar = std::pow(det_Fe, -2.0 / 3.0) * symmetrize(Fe * transpose(Fe));
      C = symmetrize(transpose(Fe) * Fe);
      Cinv = symmetrize(invert(Fe) * invert(transpose(Fe)));
      p_tilde = p_tilde_in;
      J_tilde = J_tilde_in;

      I1 = trace(C);
      Tensor<1, dim> tmp = C * M;
      I4 = M * tmp;

      Assert(det_Fe > 0, ExcInternalError());
}

template <int dim>
void Material_Compressible_TransIsotropic_Three_Field<dim>
::update_material_data(const Tensor<1, dim> &collagen_orientation,
                              double collagen_dispersion_,
                              double mu,
                              double nu,
                              double k1_)
{
    M = collagen_orientation;
    dispersion_coeff = collagen_dispersion_;
    kappa = (2.0 * mu * (1.0 + nu)) / (3.0 * (1.0 - 2.0 * nu));
    k1 = k1_;  
 
    Tensor<1, dim> tmp = C * M;
    I4 = M * tmp;
}

/*
template <int dim> 
double Material_Compressible_TransIsotropic_Three_Field<dim>::get_energy() const
{
    return get_energy_iso() + get_energy_vol();
}

template <int dim>
double Material_Compressible_TransIsotropic_Three_Field<dim>::get_energy_iso() const
{
    double I1 = C[0][0] + C[1][1] + C[2][2];
    double I1_bar = pow(det_Fe, -2.0 / 3.0) * I1;
    double energy_iso1 =  c_1 * (I1_bar - 3);
    double energy_iso4 =  0.5 * k1/k2 *(exp( k2*dispersion_coeff*I1_bar + k2*(1-3*dispersion_coeff) * sq(I4-1) ) - 1);//I4 or I4bar???
    return energy_iso1 + energy_iso4;
}

template <int dim>
double Material_Compressible_TransIsotropic_Three_Field<dim>::get_energy_vol() const
{
    return 0.0;
}
*/
    // The second function determines the Kirchhoff stress $\boldsymbol{\tau}
    // = \boldsymbol{\tau}_{\textrm{iso}} + \boldsymbol{\tau}_{\textrm{vol}}$

template <int dim>
Tensor<2, dim> Material_Compressible_TransIsotropic_Three_Field<dim>::get_S()
{
    return get_S_iso1() + get_S_iso4() + get_S_vol();
}

                                       // The fourth-order elasticity tensor
                                       // in the spatial setting
                                       // $\mathfrak{c}$ is calculated from
                                       // the SEF $\Psi$ as $ J
                                       // \mathfrak{c}_{ijkl} = F_{iA} F_{jB}
                                       // \mathfrak{C}_{ABCD} F_{kC} F_{lD}$
                                       // where $ \mathfrak{C} = 4
                                       // \frac{\partial^2
                                       // \Psi(\mathbf{C})}{\partial
                                       // \mathbf{C} \partial \mathbf{C}}$

template <int dim>
Tensor<4, dim> Material_Compressible_TransIsotropic_Three_Field<dim>::get_Cmaterial() const
{
      return get_Cmaterial_vol() + get_Cmaterial_iso1() + get_Cmaterial_iso4();
}

// Derivative of the volumetric free energy with respect to
// $\widetilde{J}$ return $\frac{\partial
// \Psi_{\text{vol}}(\widetilde{J})}{\partial \widetilde{J}}$
template <int dim>
double Material_Compressible_TransIsotropic_Three_Field<dim>::get_dPsi_vol_dJ() const
{
      return (kappa / 2.0) * (J_tilde - 1.0 / J_tilde);
}

// Second derivative of the volumetric free energy wrt $\widetilde{J}$. We
// need the following computation explicitly in the tangent so we make it
// public.  We calculate $\frac{\partial^2
// \Psi_{\textrm{vol}}(\widetilde{J})}{\partial \widetilde{J} \partial
// \widetilde{J}}$
template <int dim>
double Material_Compressible_TransIsotropic_Three_Field<dim>::get_d2Psi_vol_dJ2() const
{
      return ( (kappa / 2.0) * (1.0 + 1.0 / (J_tilde * J_tilde)));
}

// The following functions are used internally in determining the result
// of some of the public functions above. The first one determines the
// volumetric Kirchhoff stress $\boldsymbol{\tau}_{\textrm{vol}}$:


/**********************************************************************/
/************* second Piola Kirchhoff *********************************/
/**********************************************************************/
template <int dim>
SymmetricTensor<2, dim> Material_Compressible_TransIsotropic_Three_Field<dim>::get_S_vol() const
{
      return p_tilde * det_Fe * Cinv;
}

template <int dim>
SymmetricTensor<2, dim> Material_Compressible_TransIsotropic_Three_Field<dim>::get_S_iso1() const
{
    //double energy_iso1 =  c_1 * (I1_bar - 3);
    //double energy_iso4 =  0.5 * k1/k2 *(exp( k2*dispersion_coeff*I1_bar + k2*(1-3*dispersion_coeff) * sq(I4-1)) - 1);
    //double energy_iso6 =  0.5 * k1/k2 *(exp( k2*dispersion_coeff*I1_bar + k2*(1-3*dispersion_coeff) * sq(I6-1)) - 1);
    
    double I1_bar = pow(det_Fe, -2.0 / 3.0) * I1;
    double exp_term4 = exp( k2*dispersion_coeff*I1_bar + k2*(1-3*dispersion_coeff) * sq(I4-1) );
    double dPhi_iso_dI1 = c_1 * std::pow(det_Fe, -2.0/3.0) 
              + 0.5 * k1 * dispersion_coeff * std::pow(det_Fe, -2.0/3.0) * exp_term4;
    Tensor<2, dim> tmp = StandardTensors<dim>::dev_P * C;
    return 2.0 * dPhi_iso_dI1 * symmetrize( tmp * Tensor<2, dim>(Cinv) );
}

template <int dim>
Tensor<2, dim> Material_Compressible_TransIsotropic_Three_Field<dim>::get_S_iso4() const
{
    double I1_bar = pow(det_Fe, -2.0 / 3.0) * I1;
    Tensor<2, dim> M_tensormul_M = outer_product(M, M);
    //if (I4 != 1.0) std::cout << "I4 = " << I4 << "det_Fe = " << det_Fe << std::endl;
    double exp_term4 = exp( k2*dispersion_coeff*I1_bar + k2*(1-3*dispersion_coeff) *sq (I4-1) );

    double dPhi_iso_dI4 = k1 * (1-3*dispersion_coeff) * (I4-1) * exp_term4;
    //if (I4 < 1) dPhi_iso_dI4 = 0.0;
    return 2.0 * dPhi_iso_dI4 * M_tensormul_M;
}

    // Calculate the isochoric part of the tangent $J
    // \mathfrak{c}_\textrm{iso}$:

                                       // Calculate the fictitious elasticity
                                       // tensor $\overline{\mathfrak{c}}$.
                                       // For the material model chosen this
                                       // is simply zero:
template <int dim>
SymmetricTensor<4, dim> Material_Compressible_TransIsotropic_Three_Field<dim>::get_c_bar() const
{
        return SymmetricTensor<4, dim>();//initialize every items as 0.
}

template <int dim>
SymmetricTensor<4, dim> Material_Compressible_TransIsotropic_Three_Field<dim>::get_c_bar1() const
{
        double dPhi_iso_dI1 = c_1 * std::pow(det_Fe, -2.0 / 3.0);
        double ddPhi_iso_dI1_dI1=0;
        double ddPhi_iso_dI1_dI4=0;
        SymmetricTensor<2, dim> dI1_dC=StandardTensors<dim>::I;//see Hozapfel paper D1
        SymmetricTensor<2, dim> dI4_dC = SymmetricTensor<2, dim>(outer_product(M,M));//see Holzapfel paper D4
        SymmetricTensor<4, dim> ddI1_dC_dC;
        return 4*(ddPhi_iso_dI1_dI1*outer_product(dI1_dC,dI1_dC)
                   +ddPhi_iso_dI1_dI4*outer_product(dI1_dC,dI4_dC)
                   +2*dPhi_iso_dI1*ddI1_dC_dC);
}


template <int dim>
SymmetricTensor<4, dim> 
Material_Compressible_TransIsotropic_Three_Field<dim>::get_Cmaterial_vol() const
{
      SymmetricTensor<4, dim> II_Cinv = SymmetricTensor<4, dim>();
      for(int i = 0; i < dim; i++){
          for(int j = 0; j < dim; j++){
              for(int k = 0; k < dim; k++){
                  for(int l = 0;l < dim; l++){
                      II_Cinv[i][j][k][l] = 1.0 / 2.0 * (Cinv[i][k] * Cinv[j][l] + Cinv[i][l] * Cinv[j][k]);

                  }
              }
          }
      }

      SymmetricTensor<4, dim> Cinv_x_Cinv = outer_product(Cinv, Cinv);

      return p_tilde * det_Fe * (Cinv_x_Cinv - 2*II_Cinv);
}

template <int dim>
SymmetricTensor<4, dim> 
Material_Compressible_TransIsotropic_Three_Field<dim>::get_Cmaterial_iso1() const
{
    SymmetricTensor<4, dim> Cinv_x_Cinv = outer_product(Cinv, Cinv);
    SymmetricTensor<4, dim> I_x_Cinv = outer_product(StandardTensors<dim>::I, Cinv);
    SymmetricTensor<4, dim> Cinv_x_I = outer_product(Cinv, StandardTensors<dim>::I);
    SymmetricTensor<4, dim> II_Cinv = SymmetricTensor<4, dim>();
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            for(int k = 0; k < dim; k++){
                for(int l = 0;l < dim; l++){
                    II_Cinv[i][j][k][l] = -1.0 / 2.0 * (Cinv[i][k] * Cinv[j][l] + Cinv[i][l] * Cinv[j][k]);
                }
            }
        }
    }
           
    auto result1 = (4.0 / 3.0) * c_1 * std::pow(det_Fe, -2.0 / 3.0) * (
              1.0 / 3.0 * trace(C) * Cinv_x_Cinv - I_x_Cinv - Cinv_x_I - trace(C)*II_Cinv
             );
    

    double I1_bar = pow(det_Fe, -2.0 / 3.0) * I1;
    double exp_term4 = exp( k2*dispersion_coeff*I1_bar + k2*(1-3*dispersion_coeff) * (I4-1)*(I4-1) );
    double dPhi_iso_dI1 = c_1 * std::pow(det_Fe, -2.0/3.0) 
                 + 0.5 * k1 * dispersion_coeff * std::pow(det_Fe, -2.0/3.0) * exp_term4;
    double ddPhi_iso_dI1_dI1 = 0.5 * k1 * k2 * pow(dispersion_coeff,2.0) * std::pow(det_Fe, -4.0/3.0) * exp_term4;
    double dPhi_iso_dI3 = (-1.0/3.0) * c_1 * I1 * pow(det_Fe, -8.0/3.0);
    double ddPhi_iso_dI3_dI3 = 4.0/9.0 * c_1 * I1 * pow(det_Fe, -14.0/3.0);
    double ddPhi_iso_dI1_dI3 = (-1.0/3.0) * c_1 * std::pow(det_Fe, -8.0/3.0);
	        
    SymmetricTensor<2, dim> dI1_dC = StandardTensors<dim>::I;
    SymmetricTensor<2, dim> dI3_dC = sq(det_Fe) * Cinv;
    SymmetricTensor<4, dim> ddI1_dC_dC = SymmetricTensor<4, dim>();
    SymmetricTensor<4, dim> ddI3_dC_dC = sq(det_Fe) * ( outer_product(Cinv, Cinv) + II_Cinv );

    return 4 * ddPhi_iso_dI1_dI1 * outer_product(dI1_dC, dI1_dC)
          + 4 * ddPhi_iso_dI1_dI3 * (outer_product(dI1_dC, dI3_dC) + outer_product(dI3_dC, dI1_dC))
          + 4 * ddPhi_iso_dI3_dI3 * outer_product(dI3_dC, dI3_dC)
          + 4 * dPhi_iso_dI1 * ddI1_dC_dC
          + 4 * dPhi_iso_dI3 * ddI3_dC_dC;
}

template <int dim>
Tensor<4, dim> 
Material_Compressible_TransIsotropic_Three_Field<dim>::get_Cmaterial_iso4() const
{
    double I1_bar = pow(det_Fe, -2.0 / 3.0) * I1;
    double exp_term4 = exp( k2*dispersion_coeff*I1_bar + k2*(1-3*dispersion_coeff) * sq (I4-1) ); 
    double dPhi_iso_dI4 = k1 * (1-3*dispersion_coeff) * (I4-1) * exp_term4;
    double aaa = 1-3*dispersion_coeff; 
    double ddPhi_iso_dI4_dI4 = k1 * (aaa + 2 * k2 * sq(aaa) * sq(I4-1)) * exp_term4;
    double ddPhi_iso_dI3_dI4 = 0;
    ddPhi_iso_dI3_dI4 *= exp_term4;
    
    //if (I4 < 1.0) {
    //    dPhi_iso_dI4 = 0.0;
    //    ddPhi_iso_dI4_dI4 = 0.0;
    //}

    SymmetricTensor<2, dim> dI4_dC = outer_product(M, M);
    SymmetricTensor<2, dim> dI3_dC = std::pow(det_Fe, 2.0) * Cinv;
    SymmetricTensor<4, dim> ddI4_dC_dC;
    return 4 * ddPhi_iso_dI3_dI4 * ( outer_product(dI3_dC, dI4_dC) + outer_product(dI4_dC, dI3_dC) )
            + 4 * ddPhi_iso_dI4_dI4 * outer_product(dI4_dC, dI4_dC)
            + 4 * dPhi_iso_dI4 * ddI4_dC_dC;
}

}//end namespace
#endif
