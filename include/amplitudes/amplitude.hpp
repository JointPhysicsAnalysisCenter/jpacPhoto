// Abstract class for an amplitude. Used so we can easily build observables
// as the incoherent sum of amplitudes in s, t, and u channels.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _AMPLITUDE_
#define _AMPLITUDE_

// ---------------------------------------------------------------------------
// Abstract class to define helicity amplitudes. This will allow multiple different
// classes (for s, t, and u- channels but also multiple contibutions in each channel)
// to be added together and evaluated in observables.
//
// Any generic amplitude needs a reaction_kinematics object
// and a way to evaluate the helicity amplitude for given set of helicities,
// CoM energy and scattering angle.
// ---------------------------------------------------------------------------

#include "reaction_kinematics.hpp"

#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

#include <string>
#include <array>

namespace jpacPhoto
{
  class amplitude
  {
  public:
    // Constructor with only a kinematics object
    amplitude(reaction_kinematics * xkinem)
    : kinematics(xkinem)
    {};

    // Constructor with an amplitude id and number of parameters specified
    amplitude(reaction_kinematics * xkinem, std::string id, int N)
    : kinematics(xkinem), identifier(id), Nparams(N)
    {};

    // Kinematics object for thresholds and etc.
    reaction_kinematics * kinematics;

    // saved energies and angle 
    double s, t, theta;

    // Some saveable string by which to identify the amplitude
    std::string identifier;

    // How the calculate the helicity amplitude
    // Must be given a specific implementation in a user derived class
    virtual std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double t) = 0;

    // ---------------------------------------------------------------------------
    // Observables
    // Evaluatable in terms of s and t or an event object (see reaction_kinematics.hpp)

    // Modulus of the amplitude summed over all helicity combinations
    double probability_distribution(double s, double t);

    // Differential and total cross-section
    double differential_xsection(double s, double t);

    // integrated crossection
    double integrated_xsection(double s);

    // Spin asymmetries
    double A_LL(double s, double t); // Beam and target
    double K_LL(double s, double t); // Beam and recoil

    // Spin density matrix elements
    std::complex<double> SDME(int alpha, int lam, int lamp, double s, double t);

    // Beam Asymmetries
    double beam_asymmetry_y(double s, double t);    // Along the y direction
    double beam_asymmetry_4pi(double s, double t);  // integrated along phi

    // Parity asymmetry
    double parity_asymmetry(double s, double t);

    // ---------------------------------------------------------------------------
    // If helicity amplitudes have already been generated for a value of mV, s, t 
    // store them
    bool CACHED = false;
    double cached_mVec2 = 0., cached_s = 0., cached_t = 0.;
    std::array<std::complex<double>, 24> cached_helicity_amplitude;

    void check_cache(double _s, double _t);

    // ---------------------------------------------------------------------------
    // Nparams error message
    int Nparams = 0;
    void check_Nparams(std::vector<double> params)
    {
      if (params.size() != Nparams)
      {
        std::cout << "\nWarning! Invalid number of parameters (" << params.size() << ") passed to " << identifier << ".\n";
      }
    }

    // ---------------------------------------------------------------------------
    // Aliases for the above observables with option to change the produced meson mass
    inline double probability_distribution(double M2, double s, double t)
    {
      kinematics->set_mV2(M2);
      return probability_distribution(s, t);
    };

    inline double differential_xsection(double M2, double s, double t)
    {
      kinematics->set_mV2(M2);
      return differential_xsection(s, t);
    };

    inline double integrated_xsection(double M2, double s)
    {
      kinematics->set_mV2(M2);
      return integrated_xsection(s);
    };

    inline std::complex<double> SDME(int alpha, int lam, int lamp, double M2, double s, double t)
    {
      kinematics->set_mV2(M2);
      return SDME(alpha, lam, lamp, s, t);
    };

    inline double beam_asymmetry_y(double M2, double s, double t)
    {
      kinematics->set_mV2(M2);
      return beam_asymmetry_y(s, t);
    };

    inline double beam_asymmetry_4pi(double M2, double s, double t)
    {
      kinematics->set_mV2(M2);
      return beam_asymmetry_4pi(s, t);
    };

    inline double parity_asymmetry(double M2, double s, double t)
    {
      kinematics->set_mV2(M2);
      return parity_asymmetry(s, t);
    };
  };
};

#endif
