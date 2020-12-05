// Class to sum up any number of generic amplitudes and build observables.
// Amplitudes are loaded up in a vector and summed incoherently
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _SUM_
#define _SUM_

#include "amplitudes/amplitude.hpp"

// ---------------------------------------------------------------------------
// The amplitude_sum class can take a vector of the above amplitude objects
// and add them together to get observables!
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class amplitude_sum : public amplitude
  {
  private:
    // Store a vector of all the amplitudes you want to sum incoherently
    std::vector<amplitude*> amps;

  public:
    // Empty constructor
    amplitude_sum(reaction_kinematics * xkinem, std::string identifer = "")
    : amplitude(xkinem, identifer)
    {};

    // Constructor with a vector already set up
    amplitude_sum(reaction_kinematics * xkinem, std::vector<amplitude*> vec, std::string identifer = "")
    : amplitude(xkinem, identifer), amps(vec)
    {};

    // Add a new amplitude to the vector
    void add_amplitude(amplitude * new_amp)
    {
      amps.push_back(new_amp);
    };

    // Add all the members of an existing sum to a new sum
    void add_amplitude(amplitude_sum new_sum)
    {
      for (int i = 0; i < new_sum.amps.size(); i++)
      {
        amps.push_back(new_sum.amps[i]);
      }
    };

    // TODO: Add a set_params which timesi in one vector and allocates approriaten number of
    // params to each sub amplitude

    // Evaluate the sum for given set of helicites, energy, and cos
    std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double t);
  };
};

#endif
