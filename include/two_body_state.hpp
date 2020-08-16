// Energy and momentum of a two-particle state in the center of mass scattering frame
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _TWO_BODY_
#define _TWO_BODY_

#include <string>
#include <complex>
#include <iostream>

#include "misc_math.hpp"

// ---------------------------------------------------------------------------
// The two_body_state is the base object for defining a reaction in the
// s-channel center of mass scatering frame.
//
// Two particles of mass mV and mB are defined with momenta opposite along the
// same axis such that the energy and momenta of both particles is entirely
// determined by the center-of-mass energy, s, and the cosing of the angle
// from the z-axis (define to be at theta = 0).
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class two_body_state
  {
  private:
    double mV2; // Vector mass
    double mB2; // Baryon mass (allowed to be float for N* or Δ)

  public:
    // Constructor
    two_body_state(double xm1, double xm2)
    : mV2(xm1 * xm1), mB2(xm2 * xm2)
    {};

    inline double get_mV() { return sqrt(mV2); };
    inline double get_mB() { return sqrt(mB2); };

    // set masses independently
    inline void set_mV(double _mV)
    {
      mV2 = _mV * _mV;
    };

    inline void set_mB(double _mB)
    {
      mB2 = _mB * _mB;
    };

    // Momenta
    // V is always particle 1 in + z direction, 
    inline std::complex<double> momentum(double s)
    {
      return sqrt( Kallen(xr * s, xr * mV2, xr * mB2)) / (2. * sqrt(xr * s));
    };

    // Energies
    inline std::complex<double> energy_V(double s)
    {
      return (s + mV2 - mB2) / (2. * sqrt(xr * s));
    };

    inline std::complex<double> energy_B(double s)
    {
      return (s - mV2 + mB2) / (2. * sqrt(xr * s));
    };

    // Full 4-momenta 
    std::complex<double> q(int mu, double s, double theta); // 4vector of vector, particle 1
    std::complex<double> p(int mu, double s, double theta); // 4vector of baryon, particle 2
  };
};

#endif
