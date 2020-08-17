// Spin-1/2 exchange ampltiude from perturbation theory
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PROTON_
#define _PROTON_

#include <string>
#include <vector>
#include <iostream>

#include <iomanip>

#include "amplitude.hpp"
#include "gamma_technology.hpp"

namespace jpacPhoto
{
  class dirac_exchange : public amplitude
  {
  public:
    // constructor
    dirac_exchange(reaction_kinematics * xkinem, double mass, std::string name = "")
    : amplitude(xkinem, name, 2),
      mEx(mass), mEx2(mass*mass)
    {};

    // Setting utility
    void set_params(std::vector<double> params)
    {
      check_Nparams(params);
      gGam = params[0];
      gVec = params[1];
    };

    // Assemble the helicity amplitude by contracting the spinor indices
    std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double t);

    // debugging options to make either the photon or vector into scalars
    void set_debug(int i)
    {
      switch (i)
      {
        case 3: ScTOP = true; ScBOT = true; break;
        case 2: ScTOP = true; break;
        case 1: ScBOT = true; break;
      }
    }

  protected:
    // DEBUGGING PARAMS
    bool ScTOP = false, ScBOT = false;

    // Exchange nucleon mass
    double mEx, mEx2;

    // couplings
    double gGam = 0., gVec = 0.;

    // Should be exactly u_man(s, zs);
    double exchange_mass();

    // Four-momentum of the exhange (u - channel)
    std::complex<double> exchange_momentum(int mu);

    // Slashed momentumn
    std::complex<double> slashed_exchange_momentum(int i, int j);

    // Slashed polarization vectors
    std::complex<double> slashed_eps(int i, int j, double lam, polarization_vector * eps, bool STARRED, double s, double theta);

    // Photon - excNucleon - recNucleon vertex
    std::complex<double> top_vertex(int i, int lam_gam, int lam_rec);

    // excNucleon - recNucleon - Vector vertex
    std::complex<double> bottom_vertex(int j, int lam_vec, int lam_targ);

    // Spin-1/2 propagator
    std::complex<double> dirac_propagator(int i, int j);
  };
};
#endif
