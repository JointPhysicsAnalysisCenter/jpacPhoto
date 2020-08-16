// Spin-3/2 exchange amplitude from perturbation theory
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _RARITA_SCHWINGER_
#define _RARITA_SCHWINGER_

#include "dirac_exchange.hpp"

namespace jpacPhoto
{
  class rarita_exchange : public dirac_exchange
  {
  public:
    // Constructor
    rarita_exchange(reaction_kinematics * xkinem, double mass, std::string name = "")
    : dirac_exchange(xkinem, mass, name)
    {};

    // Assemble the helicity amplitude by contracting the spinor indices
    std::complex<double> helicity_amplitude(std::vector<int> helicities, double xs, double xt);

  protected:

    // rank-2 traceless tensor
    std::complex<double> g_bar(int mu, int nu);

    // g_bar contracted with gamma^nu
    std::complex<double> slashed_g_bar(int mu, int i, int j);

    // Relative momentum entering or exiting the propagator
    std::complex<double> relative_momentum(int mu, std::string in_out);

    // Spin-3/2 propagator
    std::complex<double> rarita_propagator(int i, int j);
  };
};

#endif
