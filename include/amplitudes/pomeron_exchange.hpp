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
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class pomeron_exchange : public amplitude
  {
  public:

    // Constructor
    // need a pointer to kinematic object, pointer to trajectory.
    // Optional OLDMODEL if true will default to helicity conserving amplitude
    pomeron_exchange(reaction_kinematics * xkinem, regge_trajectory * alpha, bool OLDMODEL = false, std::string name = "")
    : amplitude(xkinem, name, 2), pomeron_traj(alpha), DELTA(OLDMODEL)
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

    bool DELTA = false; // Whether or not to use the helicity conserving model 
    double norm = 0., b0 = 0.; // Regge factor parameters: normalization and t-slope
    regge_trajectory * pomeron_traj;

    // Photon - Vector - Pomeron vertex
    std::complex<double> top_vertex(int mu, int lam_gam, int lam_vec);

    // Nucleon - Nucleon - Pomeron vertex
    std::complex<double> bottom_vertex(int mu, int lam_targ, int lam_rec);

    // Energy dependence from Pomeron propogator
    std::complex<double> regge_factor();
  };
};

#endif
