// ---------------------------------------------------------------------------
// Prediction for X(6900) photoproduction near threshold based on hypothetical
// omega exchange amplitudes.
// 
// Reproduces the plot in FIG 4 of [1].
// 
// USAGE:
// make X6900 && ./X6900
//
// OUTPUT:
// omega_exchange.pdf
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:2008.01001 [hep-ph]
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

    // ---------------------------------------------------------------------------
    // Amplitude
    // ---------------------------------------------------------------------------

    // X(6900)
    double mX = 6.900;
    reaction_kinematics * kX = new reaction_kinematics(mX*mX);

    double gV_psi = 1.6E-3, gT_psi = 0.;
    double gX_psi = 5.03;
    double gGam_psi = gX_psi * e * fJpsi / mJpsi;

    double gV_omega = 16., gT_omega = 0.;
    double gX_omega = 0.225;
    double gGam_omega = gX_omega * e * fJpsi / mJpsi;
    double bOmega = 0.68; // Cuttoff of LamOmega = 1.2 GeV

    vector_exchange X_psi(kX, mJpsi, "J/#psi exchange, BR = 100%");
    X_psi.set_params({gGam_psi, gV_psi, gT_psi});
    X_psi.set_scalarX(true);

    vector_exchange X_omega(kX, mOmega, "#it{X}(6900) with BR(#it{X #rightarrow #psi#omega}) = 1%");
    X_omega.set_params({gGam_omega, gV_omega, gT_omega});
    X_omega.set_formfactor(true, bOmega);
    X_omega.set_scalarX(true);
   
    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(&X_omega);

    int N = 30;

    double  xmin = 7.5;
    double  xmax = 15.;

    double  ymin = 1.E-2;
    double  ymax = 40.;

    std::string filename = "omega_exchange.pdf";
    std::string ylabel  = "#it{#sigma(#gamma p #rightarrow X p)}   [nb]";

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    for (int n = 0; n < amps.size(); n++)
    {
        std::cout << std::endl << "Printing amplitude: " << amps[n]->identifier << "\n";

        auto F = [&](double x)
        {
            return amps[n]->integrated_xsection(x*x);
        };

        std::array<std::vector<double>, 2> x_fx; 
        if (xmin < amps[n]->kinematics->Wth())
        {
            x_fx = vec_fill(N, F, amps[n]->kinematics->Wth() + EPS, 9., true);

            for (int j = 1; j <= 10; j++)
            {
                double Wi = 9. + double(j) * (xmax - 9.) / 10.;
                x_fx[0].push_back(Wi);
                x_fx[1].push_back(F(Wi));
            }
        }
        else
        {
            x_fx = vec_fill(N, F, xmin, xmax, true);
        }

        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->identifier);
    }

    plotter->SetXaxis(ROOT_italics("W_{#gammap}") + "  [GeV]", xmin, xmax);

    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);
    
    plotter->SetLegend(0.2, 0.75);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
    delete kX;

    return 1.;
};