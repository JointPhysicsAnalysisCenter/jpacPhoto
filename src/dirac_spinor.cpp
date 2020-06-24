// Class for dirac spinors for the spin-1/2 nucleon
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dirac_spinor.hpp"

// ---------------------------------------------------------------------------
// Energy part
std::complex<double> jpacPhoto::dirac_spinor::omega(int sign, double s)
{
  if (ANTI_PARTICLE)
  {
    sign *= -1;
  }

  std::complex<double> E = twobody.energy(particle, s);
  return sqrt(E + double(sign) * mass);
}

// ---------------------------------------------------------------------------
// Angular half angle factors
double jpacPhoto::dirac_spinor::xi(int lam, double theta)
{
  double result;
  (lam == 1) ? (result = cos(theta / 2.)) : (result = sin(theta / 2.));

  return result;
};

// ---------------------------------------------------------------------------
// Components for both the regular spinor or adjoint
// Assumed to be particle 2 but moving in the +z direction
std::complex<double> jpacPhoto::dirac_spinor::component(int i, int lambda, double s, double theta)
{
  if (abs(lambda) != 1)
  {
    std::cout << "\ndirac_spinor: Invalid helicity projection passed as argument. Quitting... \n";
    exit(0);
  }

  // theta convention
  switch (i)
  {
    case 0: return                  omega(+1, s) * xi( lambda, theta);
    case 1: return double(lambda) * omega(+1, s) * xi(-lambda, theta);
    case 2: return double(lambda) * omega(-1, s) * xi( lambda, theta);
    case 3: return                  omega(-1, s) * xi(-lambda, theta);
    default : std::cout << "dirac_spinor: Invalid component index " << i << " passed as argument. Quitting... \n";
              exit(0);
  }

};

std::complex<double> jpacPhoto::dirac_spinor::adjoint_component(int i, int lambda, double s, double theta)
{
  double phase;
  (i == 2 || i == 3) ? (phase = -1.) : (phase = 1.);

  return phase * component(i, lambda, s, theta);
};
