// Class for the polarization vector of vector particles
// coded up independently to not require ROOT to be installed
//
// Dependencies: None
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "polarization_vector.hpp"

// ---------------------------------------------------------------------------
// Components
// vectors are always particle 1
std::complex<double> jpacPhoto::polarization_vector::component(int i, int lambda, double s, double theta)
{
    // Check for massless photon
    if (lambda == 0 && abs(_state->get_mV()) < 0.01)
    {   
        return 0.;
    }

    int id = 10 * abs(lambda) + i;
    switch (id)
    {
        // Longitudinal
        case 0: return _state->momentum(s) / _state->get_mV();
        case 1: return _state->energy_V(s) * sin(theta) / _state->get_mV();
        case 2: return 0.;
        case 3: return _state->energy_V(s) * cos(theta) / _state->get_mV();

        // Transverse
        case 10: return 0.;
        case 11: return - double(lambda) * cos(theta) / sqrt(2.);
        case 12: return - XI / sqrt(2.);
        case 13: return double(lambda) * sin(theta) / sqrt(2.);

        default: 
        {
            std::cout << "polarization_vector: Invalid helicity! Quitting... \n";
            return 0.; 
        }
    };

};

std::complex<double> jpacPhoto::polarization_vector::conjugate_component(int i, int lambda, double s, double theta)
{
    return conj(component(i, lambda, s, theta));
};
