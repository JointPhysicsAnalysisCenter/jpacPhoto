// Class for the polarization vector of vector particles
// coded up independently to not require ROOT to be installed
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _4VEC_
#define _4VEC_

#include <iostream>
#include "constants.hpp"
#include "two_body_state.hpp"

// ---------------------------------------------------------------------------
// Polarization vectors for vector particles
// in the s-channel center of mass frame
// Vector particles are always particle 1
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class polarization_vector
  {
  private:
      two_body_state state;
      const std::string particle;
      const double mass;

  public:
    // Constructor
    polarization_vector(two_body_state xstate, std::string name)
      : state(xstate), particle(name), mass(state.get_mass(name))
    {};

    // Copy Constructor
    polarization_vector(const polarization_vector & old)
      : state(old.state), particle(old.particle), mass(old.mass)
    {};

    // Destructor
    ~polarization_vector(){};

    // Components
    std::complex<double> component(int i, int lambda, double s, double theta);
    std::complex<double> conjugate_component(int i, int lambda, double s, double theta);
  };
};

#endif
