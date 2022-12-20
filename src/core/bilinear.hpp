// Here we assemble dirac_spinors together into Lorentz covariant bilinears
// e.g. ubar gamma_mu u
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef BILINEAR_HPP
#define BILINEAR_HPP

#include "constants.hpp"
#include "lorentz_tensor.hpp"
#include "dirac_matrix.hpp"
#include "dirac_spinor.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Covariant gamma matrix structures

    inline lorentz_tensor<dirac_matrix,1> gamma_vector()
    { return lorentz_vector<dirac_matrix>({gamma_0(), gamma_1(), gamma_2(), gamma_3()}); };

    // inline lorentz_tensor<dirac_matrix,2> sigma_tensor()
    // { return (tensor_product(gamma_vector(), gamma_vector())); };
    /* identity<dirac_matrix>()*metric_tensor() - */ 
};

#endif