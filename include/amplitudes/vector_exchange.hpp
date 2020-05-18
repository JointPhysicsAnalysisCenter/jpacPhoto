// Axial-vector meson photoproduction proceeding through a vector meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _AXIAL_
#define _AXIAL_

#include "amplitude.hpp"
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

class vector_exchange : public amplitude
{
public:
  // Constructor
  vector_exchange(reaction_kinematics * xkinem, double mass, std::string exchange = "", bool iffeyn = false)
  : amplitude(xkinem, exchange), mEx2(mass*mass), FEYN(iffeyn)
  {};

  // Copy constructor
  vector_exchange(const vector_exchange & old)
  : amplitude(old), mEx2(old.mEx2),
    gGam(old.gGam), gV(old.gV), gT(old.gT), FEYN(old.FEYN)
  {};

  // Setting utility
  void set_params(std::vector<double> params)
  {
    gGam = params[0];
    gV = params[1];
    gT = params[2];
  };

  // Calculation of the residues analytically
  std::complex<double> top_residue(int lam, double t);
  std::complex<double> bottom_residue(int lamp, double t);

  // Assemble the helicity amplitude by contracting the lorentz indices
  std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs);

private:
  // Whether to calculate with feynman rules or analytic residues
  // Give the same result but analytic is faster
  bool FEYN = false;

  // Mass of the exchange
  double mEx2;

  // Couplings to the axial-vector/photon and vector/tensor couplings to nucleon
  double gGam = 0., gpGam = 0., gV = 0., gT = 0.;

  // Four-momentum of the exhange
  std::complex<double> exchange_momenta(int mu, double s, double zs);

  // Photon - Axial Vector - Vector vertex
  std::complex<double> top_vertex(int mu, int lam_gam, int lam_vec, double s, double zs);

  // Nucleon - Nucleon - Vector vertex
  std::complex<double> bottom_vertex(int nu, int lam_targ, int lam_rec, double s, double zs);

  // Vector propogator
  std::complex<double> vector_propagator(int mu, int nu, double s, double zs);
};

#endif
