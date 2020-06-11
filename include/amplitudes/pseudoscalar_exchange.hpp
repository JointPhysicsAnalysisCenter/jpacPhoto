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
// pseudoscalar_exchange class describes the amplitude for a fixed-spin-0 exchange
// in the t-channel. Derived in terms of simple feynman rules at tree level
//
// Initialization required a reaction_kinematics object, the mass of the exchange,
// and an optional string to identify the amplitude with.
//
//  Evaluation requires two couplings:
// photon coupling, gGamma, and nucleon coupling, gNN respectively.
//
// Set couplings with amp.set_params({gGamma, gNN});
// ---------------------------------------------------------------------------

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
    gPsi = params[0];
    gNN = params[1];
  };

  // Assemble the helicity amplitude by contracting the spinor indices
  std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs);

private:
  // Whether to use fixed-spin propagator or regge
  bool REGGE;

  // Mass of the exchanged pseudo-scalar
  double mEx2;

  // Regge trajectory for the pion
  linear_trajectory * alpha;

  // Coupling constants
  double gPsi = 0., gNN = 0.;

  // Pion form factors
  double LamPi = 0.7; // pi NN vertex cutoff parameter
  double form_factor(double m, double s, double zs);

  // VMD photon vertex
  std::complex<double> top_vertex(double lam_gam, double lam_vec, double s, double zs);

  // Nucleon vertex
  std::complex<double> bottom_vertex(double lam_rec, double lam_targ, double s, double zs);

  // Simple pole propagator
  std::complex<double> scalar_propagator(double s, double zs);
};

#endif
