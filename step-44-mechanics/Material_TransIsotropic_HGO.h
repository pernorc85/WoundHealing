#ifndef _MATERIAL_TRANSISOTROPIC_H_
#define _MATERIAL_TRANSISOTROPIC_H_
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

 







// We start by including all the necessary deal.II header files and some C++
// related ones. They have been discussed in detail in previous tutorial
// programs, so you need only refer to past tutorials for details.
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

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


// We then stick everything that relates to this tutorial program into a
// namespace of its own, and import all the deal.II function and class names
// into it:
namespace Step44
{
  using namespace dealii;

// @sect3{Some standard tensors}
// Now we define some frequently used second and fourth-order tensors:
  template <int dim>
  class StandardTensors
  {
  public:

	// $\mathbf{I}$
    static const SymmetricTensor<2, dim> I;
	// $\mathbf{I} \otimes \mathbf{I}$
    static const SymmetricTensor<4, dim> IxI;
    // $\mathcal{S}$, note that as we only use this fourth-order unit tensor
    // to operate on symmetric second-order tensors.  To maintain notation
    // consistent with Holzapfel (2001) we name the tensor $\mathcal{I}$
    static const SymmetricTensor<4, dim> II;
	// Fourth-order deviatoric tensor such that
	// $\textrm{dev} \{ \bullet \} = \{ \bullet \} -
    //  [1/\textrm{dim}][ \{ \bullet\} :\mathbf{I}]\mathbf{I}$
    static const SymmetricTensor<4, dim> dev_P;
  };


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
  template <int dim>
  class Material_Compressible_TransIsotropic_Three_Field
  {
  public:
    Material_Compressible_TransIsotropic_Three_Field(const double mu,
                                               const double nu,
                                               const double dispersion_coeff_,
                                               const double k1_par,
                                               const double k2_par,
                                               Tensor<1,dim> M_par);

    ~Material_Compressible_TransIsotropic_Three_Field()
    {}

    // We update the material model with various deformation dependent data
    // based on $F$ and the pressure $\widetilde{p}$ and dilatation
    // $\widetilde{J}$, and at the end of the function include a physical
    // check for internal consistency:
    void update_material_data(const Tensor<2, dim> &Fe_par,
                              const double p_tilde_in,
                              const double J_tilde_in);
    
    void update_material_data(const Tensor<1, dim> &collagen_orientation,
                              double collagen_dispersion,
                              double mu,
                              double nu,
                              double k1_);


    Tensor<2, dim> get_S();

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

    Tensor<4, dim> get_Cmaterial() const;

    // Derivative of the volumetric free energy with respect to
    // $\widetilde{J}$ return $\frac{\partial
    // \Psi_{\text{vol}}(\widetilde{J})}{\partial \widetilde{J}}$
    double get_dPsi_vol_dJ() const;

    // Second derivative of the volumetric free energy wrt $\widetilde{J}$. We
    // need the following computation explicitly in the tangent so we make it
    // public.  We calculate $\frac{\partial^2
    // \Psi_{\textrm{vol}}(\widetilde{J})}{\partial \widetilde{J} \partial
    // \widetilde{J}}$

    double get_d2Psi_vol_dJ2() const;

    // The next few functions return various data that we choose to store with
    // the material:
    double get_det_Fe() const{
        return det_Fe;
    }

    double get_p_tilde() const{
        return p_tilde;
    }

    double get_J_tilde() const{
        return J_tilde;
    }

  protected:
    // Define constitutive model paramaters $\kappa$ (bulk modulus)and the
    // neo-Hookean model parameter $c_1$:
    double kappa;
    const double c_1;
    double k1;
    const double k2;
    double dispersion_coeff;
    Tensor<1,dim> M;
    Tensor<2,dim> Fe;

    // Model specific data that is convenient to store with the material:
    double det_Fe;
    double p_tilde;
    double J_tilde;
    SymmetricTensor<2, dim> b_bar;
    SymmetricTensor<2, dim> C;
    SymmetricTensor<2, dim> Cinv;
    double I1, I4;


    /**********************************************************************/
    /************* second Piola Kirchhoff *********************************/
    /**********************************************************************/
    SymmetricTensor<2, dim> get_S_vol() const;

    SymmetricTensor<2, dim> get_S_iso1() const;

    Tensor<2, dim> get_S_iso4() const;


    SymmetricTensor<4, dim> get_c_bar() const;

    SymmetricTensor<4, dim> get_c_bar1() const;

    SymmetricTensor<4, dim> get_Cmaterial_vol() const;

    SymmetricTensor<4, dim> get_Cmaterial_iso1() const;

    Tensor<4, dim> get_Cmaterial_iso4() const;

};

}//end namespace
#endif
