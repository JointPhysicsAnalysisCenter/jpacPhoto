// Misc math functions that are useful to have
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _MATH_
#define _MATH_

#include "constants.hpp"
#include <iostream>
#include <complex>

namespace jpacPhoto
{
  template <typename T>
  T Kallen(T x, T y, T z)
  {
    return x*x + y*y + z*z - 2. * (x*y + x*z + y*z);
  };

  std::complex<double> cgamma(std::complex<double> z, int OPT = 0);

  // ---------------------------------------------------------------------------
  // Wigner d-function for half-integer spin
  double wigner_d_half(int j, int lam1, int lam2, double theta);

  // Wigner d-function for integer spin
  double wigner_d_int(int j, int lam1, int lam2, double theta);

  // Error message display function for the above
  double wigner_error(int j, int lam1, int lam2, bool half);
};

#endif
