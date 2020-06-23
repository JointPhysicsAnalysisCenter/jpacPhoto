// Class to contain all relevant kinematic quantities. The kinematics of the reaction
// gamma p -> V p' is entirely determined by specifying the mass of the vector particle
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------


#ifndef _KINEMATICS_
#define _KINEMATICS_

#include "constants.hpp"
#include "misc_math.hpp"
#include "two_body_state.hpp"
#include "dirac_spinor.hpp"
#include "polarization_vector.hpp"

#include "TLorentzVector.h"

#include <vector>
#include <string>
#include <cmath>

namespace jpacPhoto
{
  // Simple struct to pass a full set of 4-vectors
  struct event
  {
    event(TLorentzVector pG, TLorentzVector pT, TLorentzVector qV, TLorentzVector qR)
    : pGam(pG), pTarg(pT), qVec(qV), qRec(qR)
    {};

    TLorentzVector pGam;
    TLorentzVector pTarg;
    TLorentzVector qVec;
    TLorentzVector qRec;

    // Mandelstam Variables
    double s_man()
    {
      return (pGam + pTarg).M2();
    }

    double t_man()
    {
      return (pGam - qVec).M2();
    }

    double u_man()
    {
      return (pGam - qRec).M2();
    }
  };

  // ---------------------------------------------------------------------------
  // The reaction kinematics object is intended to have all relevant kinematic quantities
  // forthe reaction. Here you'll find the momenta and energies of all particles,
  //  spinors for the nucleons and polarization vectors for the gamma and vector meson
  //
  // Additionally helicity combinations are stored for easier access
  // ---------------------------------------------------------------------------

  class reaction_kinematics
  {
  public:
    // Constructor
    reaction_kinematics(double vec_mass, std::string vec_name)
    : mVec(vec_mass), mVec2(vec_mass*vec_mass), vector_particle(vec_name),
      initial(0., mPro, "beam", "target"),
      final(vec_mass, mPro, vec_name, "recoil"),
      eps_vec(final, vec_name),
      eps_gamma(initial, "beam"),
      target(initial, "target"), recoil(final, "recoil")
    {};

  // Copy Constructor
    reaction_kinematics(const reaction_kinematics & old)
    : mVec(old.mVec), mVec2(old.mVec2), vector_particle(old.vector_particle),
      initial(old.initial), final(old.final),
      eps_vec(old.eps_vec), eps_gamma(old.eps_gamma),
      target(old.target), recoil(old.recoil)
    {};

    const std::string vector_particle;
    const double mVec, mVec2;

    // inital and final state kinematics
    const double sth = (mVec + mPro) * (mVec + mPro);
    two_body_state initial, final;
    polarization_vector eps_vec, eps_gamma;
    dirac_spinor target, recoil;

    // Get s, t, u from 4-vectors
    double s_man(event fvecs);
    double z_s(event fvecs);
    double z_s(double s, double t);

    // Invariant variables
    double t_man(double s, double theta);
    double u_man(double s, double theta);

    // Scattering angles
    double theta_s(double s, double t);
    std::complex<double> z_t(double s, double theta);
    std::complex<double> z_u(double s, double theta);

    // Cosine of crossing angles
    std::complex<double> crossing_angle(std::string particle, double s, double theta);

    // Helicity configurations
    // Photon [0], Incoming Proton [1], Vector meson [2], Outgoing Proton [3]
    std::vector< std::vector<int> > helicities =
    {
    //  {  γ,  p,  V,  p'}
        {  1, -1,  1, -1 }, // 0
        {  1, -1,  1,  1 }, // 1
        {  1, -1,  0, -1 }, // 2
        {  1, -1,  0,  1 }, // 3
        {  1, -1, -1, -1 }, // 4
        {  1, -1, -1,  1 }, // 5
        {  1,  1,  1, -1 }, // 6
        {  1,  1,  1,  1 }, // 7
        {  1,  1,  0, -1 }, // 8
        {  1,  1,  0,  1 }, // 9
        {  1,  1, -1, -1 }, // 10
        {  1,  1, -1,  1 }, // 11
        { -1, -1,  1, -1 }, // 12
        { -1, -1,  1,  1 }, // 13
        { -1, -1,  0, -1 }, // 14
        { -1, -1,  0,  1 }, // 15
        { -1, -1, -1, -1 }, // 16
        { -1, -1, -1,  1 }, // 17
        { -1,  1,  1, -1 }, // 18
        { -1,  1,  1,  1 }, // 19
        { -1,  1,  0, -1 }, // 20
        { -1,  1,  0,  1 }, // 21
        { -1,  1, -1, -1 }, // 22
        { -1,  1, -1,  1 }, // 23
    };
  };
};

#endif
