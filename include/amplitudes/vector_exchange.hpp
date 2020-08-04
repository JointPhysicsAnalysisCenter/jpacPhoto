// Axial-vector meson photoproduction proceeding through a vector meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _AXIAL_
#define _AXIAL_

#include "amplitude.hpp"
#include "regge_trajectory.hpp"
#include "gamma_technology.hpp"

// ---------------------------------------------------------------------------
// vector_exchange class describes the amplitude for a fixed-spin-1 exchange
// in the t-channel. Derived in terms of simple feynman rules at tree level
//
// Initialization required a reaction_kinematics object, the mass of the exchange,
// and an optional string to identify the amplitude with.
//
//  Evaluation requires three couplings photon coupling, gGamma, and vector/tensor
// nucleon couplings, gV and gT respectively.
//
// Set couplings with amp.set_params({gGamma, gV, gT});
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class vector_exchange : public amplitude
  {
  public:
    // Constructor for fixed spin
    vector_exchange(reaction_kinematics * xkinem, double mass, std::string exchange = "")
    : amplitude(xkinem, exchange, 3), mEx2(mass*mass), REGGE(false)
    {};

    // Constructor for the reggized)
    vector_exchange(reaction_kinematics * xkinem, linear_trajectory * traj, std::string exchange = "")
    : amplitude(xkinem, exchange, 3), alpha(traj), REGGE(true)
    {};

    // Setting utility
    inline void set_params(std::vector<double> params)
    {
      check_Nparams(params); // make sure the right amout of params passed
      gGam = params[0];
      gV = params[1];
      gT = params[2];
    };

    // Whether or not to include an exponential form factor (default false)
    inline void set_formfactor(bool FF, double bb = 0.)
    {
      IF_FF = FF;
      b = bb;
    }

    // Special case of scalar photoproduction relevant for the X(6900)
    inline void set_scalarX(bool X)
    {
      IF_SCALAR_X = X;
      FOUR_VEC = true; // Currently scalar top vertex is not implemented in the faster analytic expression, need to do 4-vector algebra
    };

    // Assemble the helicity amplitude by contracting the lorentz indices
    std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double t);

  private:
    // Saved energies
    double s, t, theta;
    double zt;

    // if using reggeized propagator
    bool REGGE;
    // or the regge trajectory of the exchange
    linear_trajectory * alpha;

    // Whether using analytic or covariant expression
    bool FOUR_VEC = false;

    // Form factor parameters
    bool IF_FF = false;
    double b = 0.;

    // IF to treat the produced particle as a scalar
    bool IF_SCALAR_X = false;

    // Couplings to the axial-vector/photon and vector/tensor couplings to nucleon
    double gGam = 0., gpGam = 0., gV = 0., gT = 0.;

    // ---------------------------------------------------------------------------
    // Covariant evaluation

    // Mass of the exchange
    double mEx2;

    // Full covariant amplitude
    std::complex<double> covariant_amplitude(std::vector<int> helicities);
    
    // Four-momentum of the exhange
    std::complex<double> exchange_momenta(int mu);

    // Photon - Axial Vector - Vector vertex
    std::complex<double> top_vertex(int mu, int lam_gam, int lam_vec);

    // Nucleon - Nucleon - Vector vertex
    std::complex<double> bottom_vertex(int nu, int lam_targ, int lam_rec);

    // Vector propogator
    std::complex<double> vector_propagator(int mu, int nu);

    // ---------------------------------------------------------------------------
    // Analytic evaluation

    // Photon - Axial - Vector
    std::complex<double> top_residue(int lam_gam, int lam_vec);

    // Nucleon - Nucleon - Vector
    std::complex<double> bottom_residue(int lam_targ, int lam_rec);

    // Reggeon propagator
    std::complex<double> regge_propagator(int j, int lam, int lamp);

    // Half angle factors
    std::complex<double> half_angle_factor(int lam, int lamp);

    // Angular momentum barrier factor
    std::complex<double> barrier_factor(int j, int M);
  };
};

#endif
