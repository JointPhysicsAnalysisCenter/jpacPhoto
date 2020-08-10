// ---------------------------------------------------------------------------
// Example of using the amplitudes as weight functions for an event generator
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/baryon_resonance.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "toy_monte_carlo.hpp"

#include <iostream>
#include <string>

using namespace jpacPhoto;

int main()
{
    // ---------------------------------------------------------------------------
    // AMPLITUDES
    // ---------------------------------------------------------------------------

    // Set up Kinematics
    reaction_kinematics * ptr = new reaction_kinematics(mJpsi, "jpsi");

    // ---------------------------------------------------------------------------
    // S - CHANNEL

    // Two different pentaquarks
    baryon_resonance P_c4450(ptr, 3, -1, 4.45, 0.040, "P_{c}(4450)");
    P_c4450.set_params({0.01, .7071}); // 1% branching fraction and equal photocouplings

    baryon_resonance P_c4380(ptr, 5, +1, 4.38, 0.205, "P_{c}(4380)");
    P_c4380.set_params({0.01, .7071}); // 1% branching fraction and equal photocouplings

    // ---------------------------------------------------------------------------
    // T - CHANNEL

    // Set up pomeron trajectory
    linear_trajectory alpha(+1, 0.941, 0.364, "pomeron");

    // Create amplitude with kinematics and trajectory
    pomeron_exchange background(ptr, &alpha, false, "Background");
    background.set_params({0.379, 0.12});

    // ---------------------------------------------------------------------------
    // SUM

    // Incoherent sum of the s and t channels
    amplitude_sum sum(ptr, {&background, &P_c4450, &P_c4380}, "10q Sum");

    // ---------------------------------------------------------------------------
    // MC GENERATOR
    // ---------------------------------------------------------------------------

    std::string filename = "5q_mc.root";

    // Make a toy mc object
    toy_monte_carlo mc(filename);

    // Pass the above preset amplitude
    mc.set_amplitude(&sum); 

    double beam_energy = 10.;
    int N = 220;

    mc.generate(beam_energy, N);

    delete ptr;
    return 1;
};