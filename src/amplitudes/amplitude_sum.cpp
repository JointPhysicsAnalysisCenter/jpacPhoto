// Class to sum up any number of generic amplitudes and build observables.
// Amplitudes are loaded up in a vector and summed incoherently
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/amplitude_sum.hpp"

// Evaluate the sum for given set of helicites, and mandelstam invariant s and t
std::complex<double> jpacPhoto::amplitude_sum::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
    std::complex<double> result = 0.;
    for (int i = 0; i < _amps.size(); i++)
    {
        result += _amps[i]->helicity_amplitude(helicities, s, t);
    }

    return result;
};
