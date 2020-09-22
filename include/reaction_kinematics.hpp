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

#include "TMath.h"

#include <vector>
#include <string>
#include <cmath>

namespace jpacPhoto
{
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
    // Empty constructor,
    // defaults to compton scattering: gamma p -> gamma p
    reaction_kinematics()
    : mVec(0.), mVec2(0.), 
      mBar(mPro), mBar2(mPro2),
      Q2(0.)
    {
      initial   = new two_body_state(0., mPro2);
      eps_gamma = new polarization_vector(initial);
      target    = new dirac_spinor(initial);

      final     = new two_body_state(0., mPro2);
      eps_vec   = new polarization_vector(final);
      recoil    = new dirac_spinor(final);
    };

    // Constructor with a set vector mass mV2, 
    // defaults to proton as baryon and real photon
    reaction_kinematics(double mV2)
    : mVec(sqrt(mV2)), mVec2(mV2),
      mBar(mPro), mBar2(mPro2),
      Q2(0.)
    {
      initial   = new two_body_state(0., mPro2);
      eps_gamma = new polarization_vector(initial);
      target    = new dirac_spinor(initial);

      final     = new two_body_state(mV2, mPro2);
      eps_vec   = new polarization_vector(final);
      recoil    = new dirac_spinor(final);
    };

    // Constructor with a set mV2 and baryon mass mB2
    // defaults to real photon
    reaction_kinematics(double mV2, double mB2)
    : mVec(sqrt(mV2)), mVec2(mV2),
      mBar(sqrt(mB2)), mBar2(mB2),
      Q2(0.)
    {
      initial   = new two_body_state(0., mB2);
      eps_gamma = new polarization_vector(initial);
      target    = new dirac_spinor(initial);

      final     = new two_body_state(mV2, mB2);
      eps_vec   = new polarization_vector(final);
      recoil    = new dirac_spinor(final);
    };

    // Constructor with a set mV2 and baryon mass mB2
    // and q2 > 0 for a virtual photon
    reaction_kinematics(double mV2, double mB2, double q2)
    : mVec(sqrt(mV2)), mVec2(mV2),
      mBar(sqrt(mB2)), mBar2(mB2),
      Q2(q2)
    {
      initial   = new two_body_state(-q2, mB2);
      eps_gamma = new polarization_vector(initial);
      target    = new dirac_spinor(initial);

      final     = new two_body_state(mV2, mB2);
      eps_vec   = new polarization_vector(final);
      recoil    = new dirac_spinor(final);
    };

    // destructor
    ~reaction_kinematics()
    {
      delete initial, final;
      delete eps_gamma, eps_vec;
      delete target, recoil;
    }

    double mVec = 0., mVec2 = 0.; // mass and mass squared of the vector particle
    double mBar = 0., mBar2 = 0.; // mass nd mass squared of the baryon 
    double Q2 = 0.; // virtuality of the photon
    inline double Wth(){ return (mVec + mBar); }; // square root of the threshold
    inline double sth(){ return Wth() * Wth(); }; // final state threshold


    // Change the vector mass squared
    inline void set_mV2(double m2)
    {
      mVec  = sqrt(m2);
      mVec2 = m2;

      // also update the vector mass in two_body_state
      final->set_mV2(m2);
    };

    // Change the baryon mass squared
    inline void set_mB2(double m2)
    {
      mBar  = sqrt(m2);
      mBar2 = m2;

      // also update the baryon mass in two_body_state
      initial->set_mB2(m2);
      final->set_mB2(m2);
    };
    
    // Change virtuality of the photon
    // Q2 > 0
    inline void set_Q2(double q2)
    {
      if (q2 > 0) { std::cout << "Caution! set_Q2(x) requires x > 0! \n"; }
      Q2 = q2;
      initial->set_mV2(-q2);
    };
    
    // Helicity configurations
    // Photon [0], Incoming Proton [1], Vector meson [2], Outgoing Proton [3], hel_id [4]
    std::vector< std::vector<int> > helicities =
    {
    //  {  γ,  p,  V,  p'}
        {  1, -1,  1, -1 ,  0},
        {  1, -1,  1,  1 ,  1},
        {  1, -1,  0, -1 ,  2},
        {  1, -1,  0,  1 ,  3},
        {  1, -1, -1, -1 ,  4},
        {  1, -1, -1,  1 ,  5},
        {  1,  1,  1, -1 ,  6},
        {  1,  1,  1,  1 ,  7},
        {  1,  1,  0, -1 ,  8},
        {  1,  1,  0,  1 ,  9},
        {  1,  1, -1, -1 , 10},
        {  1,  1, -1,  1 , 11},
        { -1, -1,  1, -1 , 12},
        { -1, -1,  1,  1 , 13},
        { -1, -1,  0, -1 , 14},
        { -1, -1,  0,  1 , 15},
        { -1, -1, -1, -1 , 16},
        { -1, -1, -1,  1 , 17},
        { -1,  1,  1, -1 , 18},
        { -1,  1,  1,  1 , 19},
        { -1,  1,  0, -1 , 20},
        { -1,  1,  0,  1 , 21}, 
        { -1,  1, -1, -1 , 22},
        { -1,  1, -1,  1 , 23},
    };

    //--------------------------------------------------------------------------
    two_body_state * initial,  * final;
    polarization_vector * eps_vec, * eps_gamma;
    dirac_spinor * target, * recoil;

    // Get s-channel scattering angle from invariants
    inline double z_s(double s, double t)
    {
      std::complex<double> qdotqp = initial->momentum(s) * final->momentum(s);
      std::complex<double> E1E3   = initial->energy_V(s) * final->energy_V(s);

      double result = t - mVec2 + Q2 + 2.*abs(E1E3);
      result /= 2. * abs(qdotqp);

      return result;
    };

    // Scattering angle in the s-channel
    // Use TMath::ACos instead of std::acos because its safer at the end points
    inline double theta_s(double s, double t)
    {
      return TMath::ACos( z_s(s, t) );
    };

    // Invariant variables
    inline double t_man(double s, double theta)
    {
      std::complex<double> qdotqp = initial->momentum(s) * final->momentum(s);
      std::complex<double> E1E3 = initial->energy_V(s) * final->energy_V(s);

      return mVec*mVec - Q2 - 2. * abs(E1E3) + 2. * abs(qdotqp) * cos(theta);
    };

    inline double u_man(double s, double theta)
    { 
      return mVec2 - Q2 + 2. * mBar2 - s - t_man(s, theta);
    };

    // Scattering angles in t and u channel frames
    inline std::complex<double> z_t(double s, double theta)
    {
      double t = t_man(s, theta);
      std::complex<double> p_t = sqrt(xr * t - 4. * mBar2) / 2.;
      std::complex<double> q_t = sqrt(xr * Kallen(t, mVec2, -Q2)) / sqrt(xr * 4. * t);

      std::complex<double> result;
      result = 2. * s + t - 2. * mBar2 - mVec2 + Q2; // s - u
      result /= 4. * p_t * q_t;

      return result;
    };

    inline std::complex<double> z_u(double s, double theta)
    {
      // TODO: fix this
      std::cout << "z_u not implimented yet i fucked up here :p\n";
      std::cout << "if youre seeing this message email me and yell at me because i forgot\n";
      
      return xr;
    };
  };
};

#endif
