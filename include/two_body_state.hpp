// Base class to work with four vectors in a scattering event
//
// All code can be easily modified to use TLorentzVector instead
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
// Two particles of mass m1 and m2 are defined with momenta opposite along the
// same axis such that the energy and momenta of both particles is entirely
// determined by the center-of-mass energy, s, and the cosing of the angle
// from the z-axis (define to be at theta = 0).
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class two_body_state
  {
  private:
    const double m1, m2;
    const std::string particle1, particle2;

  public:
    // Constructor
    two_body_state(double xm1, double xm2, std::string name1, std::string name2)
    : m1(xm1), m2(xm2), particle1(name1), particle2(name2)
    {};

    // Copy Constructor
    two_body_state(const two_body_state & old)
    : m1(old.m1), m2(old.m2),
      particle1(old.particle1), particle2(old.particle2)
    {};

    double get_mass(std::string name);

    std::complex<double> momentum(std::string name, double s);
    std::complex<double> energy(std::string name, double s);

    std::complex<double> component(int i, std::string name, double s, double theta);
  };
};

#endif
