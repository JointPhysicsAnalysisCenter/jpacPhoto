// ---------------------------------------------------------------------------
// 
// 
// 
// USAGE:
// make primakoff && ./primakoff
//
// OUTPUT:
// x_primakoff.pdf
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// REFERENCES:
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/primakoff_effect.hpp"
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
    double Q2 = 0.5;

    double mX = 3.872;
    reaction_kinematics * kX = new reaction_kinematics(mX);

    primakoff_effect U(kX, "^{238}U");
    U.set_params({92, 6.8054, 0.556, 3.2E-3});
    // U.set_Q2(Q2);

    primakoff_effect Yb(kX, "^{176}Yb");
    Yb.set_params({70, 6.3306, 0.486, 3.2E-3});
    // Yb.set_Q2(Q2);

    primakoff_effect Zn(kX, "^{70}Zn");
    Zn.set_params({30, 4.044, 0.583, 3.2E-3});
    // Zn.set_Q2(Q2);
   

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<primakoff_effect*> amps;
    amps.push_back(&Zn);
    amps.push_back(&Yb);
    amps.push_back(&U);

    int N = 50;

    double  xmin = 7.5;
    double  xmax = 30.;

    double  ymin = 2.E-2;
    double  ymax = 70.;

    // double  ymin = 1.E-2;
    // double  ymax = 10.;

    std::string filename = "primakoff_0.5.pdf";
    std::string ylabel  = "#it{#sigma(#gamma N #rightarrow X N)}   [nb]";


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
            x_fx = vec_fill(N, F, amps[n]->kinematics->Wth(), xmax, true);

        }
        else
        {
            x_fx = vec_fill(N, F, xmin, xmax, true);
        }

        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->identifier);
    }

    plotter->SetXaxis(ROOT_italics("W_{#gammaN}") + "  [GeV^{2}]", xmin, xmax);

      // Add a header to legend to specify the fixed energy
    std::ostringstream streamObj;
    streamObj << std::setprecision(4) << Q2;
    plotter->SetLegend(0.2, 0.76, "Q^{2} = " + streamObj.str() + " GeV^{2}");

    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYaxis(ylabel);
    

    plotter->SetYlogscale(true);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
    delete kX;

    return 1.;
};