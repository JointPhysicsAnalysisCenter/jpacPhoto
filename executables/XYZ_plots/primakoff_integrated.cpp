// --------------------------------------------------------------------------- 
// Prediction for X(3872) Primakoff production off a nuclear target
//
// USAGE:
// make primakoff_integrated && ../bin/primakoff_integrated
//
// OUTPUT:
// primakoff_integrated.pdf
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
    double Q2 = 0.1;
    double mX = 3.872;

    // Uranium
    double mU = 221.6977;
    reaction_kinematics * kU = new reaction_kinematics(mX*mX, mU*mU, Q2);

    primakoff_effect U(kU, "^{238}U");
    U.set_params({92, 34.48, 3.07, 3.2E-3});

    // Tin
    double mSn = 115.3924;
    reaction_kinematics * kSn = new reaction_kinematics(mX*mX, mSn*mSn, Q2);

    primakoff_effect Sn(kSn, "^{124}Sn");
    Sn.set_params({50, 27.56, 2.73, 3.2E-3});

    // Zinc
    double mZn = 65.1202;
    reaction_kinematics * kZn = new reaction_kinematics(mX*mX, mZn*mZn, Q2);

    primakoff_effect Zn(kZn, "^{70}Zn");
    Zn.set_params({30, 22.34, 2.954, 3.2E-3});

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<primakoff_effect*> amps;
    amps.push_back(&Zn);
    amps.push_back(&Sn);
    amps.push_back(&U);

    double xNs[3] = {70., 124., 238.}; // number of nucleons 

    int N = 50;
    std::string filename = "primakoff_integrated.pdf";

    // X axis params
    double  xmax =  5.;
    std::string xlabel  = "#it{W_{#gammaN}}    [GeV]";

    // y - axis params
    double  ymin = 2.E-1;
    double  ymax = 2.E5;
    std::string ylabel  = "#it{#sigma_{L} (#gamma* A #rightarrow X A)}     [nb]";

 

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
        
        double Wth = (amps[n]->kinematics->Wth() + EPS) / xNs[n];

        auto F = [&](double x)
        {
            x *= xNs[n];
            return amps[n]->integrated_xsection(x*x);
        };

        std::array<std::vector<double>, 2> x_fx; 
        x_fx = vec_fill(N, F, Wth, xmax, true);

        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->identifier);
    }

      // Add a header to legend to specify the fixed Q2
    std::ostringstream streamObj;
    streamObj << std::setprecision(4) << "Q^{2} = " << Q2 << " GeV^{2}";
    std::string header = streamObj.str();
    plotter->SetLegend(0.2, 0.66, header);

    // Set up axes
    plotter->SetXaxis(xlabel, kZn->Wth() / xNs[0], xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);

    // Output to file
    plotter->Plot(filename);

    // Cleanup
    delete plotter;
    delete kU, kSn, kZn;

    return 1.;
};