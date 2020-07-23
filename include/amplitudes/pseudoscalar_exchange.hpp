// Charged axial-vector meson photoproduction proceeding through a pseudoscalar (pion) exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:1503.02125 [hep-ph]
// ---------------------------------------------------------------------------

#ifndef _PSCALAR_
#define _PSCALAR_

#include "amplitude.hpp"
#include "regge_trajectory.hpp"
#include "gamma_technology.hpp"

// ---------------------------------------------------------------------------
// pseudoscalar_exchange class describes the amplitude for a spin-0 exchange
// in the t-channel. Derived in terms of simple feynman rules at tree level.
//
// Initialization required a reaction_kinematics object.
// Then either: the mass (in GeV) of the exchange (for fixed-spin exchange),
//          or: a pointer to a linear_trajectory object (for Reggeize exchange).
// and an optional string to identify the amplitude with.
//
// Evaluation requires two couplings:
// photon coupling, gGamma, and nucleon coupling, gNN respectively.
//
// Set couplings with amp.set_params({gGamma, gNN});
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class pseudoscalar_exchange : public amplitude
  {
  public:
    // constructor for fixed meson exchange
    pseudoscalar_exchange(reaction_kinematics * xkinem, double mass, std::string name = "")
    : amplitude(xkinem, name, 2), mEx2(mass*mass), REGGE(false)
    {};

    // constructors for regge exchange
    pseudoscalar_exchange(reaction_kinematics * xkinem, linear_trajectory * traj, std::string name = "")
    : amplitude(xkinem, name, 2), alpha(traj), REGGE(true)
    {};

    // Setting utility
    void set_params(std::vector<double> params)
    {
      check_Nparams(params);
      gGamma = params[0];
      gNN = params[1];
    };

    // Whether or not to include an exponential form factor (default false)
    void set_formfactor(bool FF, double bb = 0.)
    {
      IF_FF = FF;
      b = bb;
    }

    // Assemble the helicity amplitude by contracting the spinor indices
    std::complex<double> helicity_amplitude(std::vector<int> helicities, double xs, double xt);

  private:
    // Place to save fixed energies (and theta)
    double s, t, theta;

    // Whether to use fixed-spin propagator (false) or regge (true)
    bool REGGE;

    // Mass of the exchanged pseudo-scalar (if REGGE = false)
    // ignored otherwise
    double mEx2;

    // Regge trajectory for the pion (if REGGE = true)
    // ignored otherwise
    linear_trajectory * alpha;

    // Coupling constants
    double gGamma = 0.; // Gamma - Axial - Pseudoscalar coupling 
    double gNN = 0.;    // Pseudoscalar - Nucleon coupling

    bool IF_FF = false; // Whether to include the exponential form factor
    double b = 0.; // "t-slope" parameter in the FF

    // Photon - pseudoscalar - Axial vertex
    std::complex<double> top_vertex(double lam_gam, double lam_vec);

    // Pseudoscalar - Nucleon vertex
    std::complex<double> bottom_vertex(double lam_rec, double lam_targ);

    // Simple pole propagator
    std::complex<double> scalar_propagator();
  };
};

#endif
