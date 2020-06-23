// Vector meson photoproduction dynamics proceeding through a pomeron exchange
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _POMERON_
#define _POMERON_

#include "amplitude.hpp"
#include "gamma_technology.hpp"
#include "regge_trajectory.hpp"

// ---------------------------------------------------------------------------
// The pomeron_exchange class describes the amplitude correspinding to
// vector meson photoproduction via a vector pomeron coupling.
//
// As a reggeon exchange model it is parameterized in terms three functions
// 1. top_vertex() coupling the vector meson to the incoming photon
// 2. bottom_vertex() coupling the two proton dirac spinors
// 3. regge_factor() the function describing the energy dependence of the amplitude
//
// To define a pomeron_exchange object all is needed is a pointer to
// a reaction_kinematics object containing the mass and name of the vector particle
// being produced
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class pomeron_exchange : public amplitude
  {
  public:
    // Constructor
    pomeron_exchange(reaction_kinematics * xkinem, regge_trajectory * alpha, std::string name = "")
    : amplitude(xkinem, name, 2), pomeron_traj(alpha)
    {};

    // Copy constructor
    pomeron_exchange(const pomeron_exchange & old)
    : amplitude(old),
      norm(old.norm), b0(old.b0)
    {};

    // Setting utility
    void set_params(std::vector<double> params)
    {
      check_Nparams(params);
      norm = params[0];
      b0 = params[1];
    };

    // Assemble the helicity amplitude by contracting the lorentz indices
    std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double t);

  private:

    double norm, b0; // Regge factor parameters: normalization and t-slope
    regge_trajectory * pomeron_traj;

    // Photon - Vector - Pomeron vertex
    std::complex<double> top_vertex(int mu, int lam_gam, int lam_vec, double s, double theta);

    // Nucleon - Nucleon - Pomeron vertex
    std::complex<double> bottom_vertex(int mu, int lam_targ, int lam_rec, double s, double theta);

    // Energy dependence from Pomeron propogator
    std::complex<double> regge_factor(double s, double t);
  };
};

#endif
