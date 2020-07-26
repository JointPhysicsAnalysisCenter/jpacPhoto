// ---------------------------------------------------------------------------
// Photoproduction of Y states by Pomeron exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
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
    // Preliminaries
    // ---------------------------------------------------------------------------

    // Chi_c1(1P)
    double mChi = 3.510;
    reaction_kinematics * kChi = new reaction_kinematics(mChi, "chi_c1");

    // X(3872)
    double mX = 3.87169;
    reaction_kinematics * kX = new reaction_kinematics(mX, "X(3872)");

    // Nucleon couplings 
    double gV_omega = 16., gT_omega = 0.;
    double bOmega = 0.68; // Cuttoff of LamOmega = 1.2 GeV

    double gV_rho = 2.4, gT_rho = 14.6;
    double bRho = 0.51; // Cuttoff of LamRho = 1.4 GeV

    double gV_phi = -6.2, gT_phi = 2.1;

    double gV_psi = 1.6E-3, gT_psi = 0.;

    // ---------------------------------------------------------------------------
    // Low-Energy Amplitudes
    // ---------------------------------------------------------------------------

    ////////////////////
    // Chi_c1(1P)
    double gChi_omega = 5.2E-4;
    vector_exchange Chi_omega(kChi, mOmega, "#omega");
    Chi_omega.set_params({gChi_omega, gV_omega, gT_omega});
    Chi_omega.set_formfactor(true, bOmega);

    double gChi_rho = 9.2E-4;
    vector_exchange Chi_rho(kChi, mRho, "#rho");
    Chi_rho.set_params({gChi_rho, gV_rho, gT_rho});
    Chi_rho.set_formfactor(true, bRho);

    double gChi_phi = 4.2E-4;
    vector_exchange Chi_phi(kChi, mPhi, "#phi");
    Chi_phi.set_params({gChi_phi, gV_phi, gT_phi});
    
    double gChi_psi = 1.;
    vector_exchange Chi_psi(kChi, mJpsi, "#psi");
    Chi_psi.set_params({gChi_psi, gV_psi, gT_psi});
    
    std::vector<amplitude*> chi_exchanges = {&Chi_omega, &Chi_rho, &Chi_phi, &Chi_psi};
    amplitude_sum chi(kChi, chi_exchanges, "#chi_{c1}(1P)");
    //////////////////

    //////////////////
    // X(3872)
    double gX_omega = 9.51E-3;
    vector_exchange X_omega(kX, mOmega, "#omega");
    X_omega.set_params({gX_omega, gV_omega, gT_omega});
    X_omega.set_formfactor(true, bOmega);

    double gX_rho = 3.81E-3;
    vector_exchange X_rho(kX, mRho, "#rho");
    X_rho.set_params({gX_rho, gV_rho, gT_rho});
    X_rho.set_formfactor(true, bRho);

    std::vector<amplitude*> X_exchanges = {&X_omega, &X_rho};
    amplitude_sum X(kX, X_exchanges, "#it{X}(3872)");
    /////////////////

    // ---------------------------------------------------------------------------
    // High-enegy Amplitudes
    // ---------------------------------------------------------------------------

    // Linear trajectory for the rho
    linear_trajectory alpha(-1, 0.5, 0.9, "#rho - #omega");

    // ////////////////////
    // Chi_c1(1P)
    vector_exchange Chi_omegaR(kChi, &alpha, "#omega");
    Chi_omegaR.set_params({gChi_omega, gV_omega, gT_omega});
    Chi_omegaR.set_formfactor(true, bOmega);

    vector_exchange Chi_rhoR(kChi, &alpha, "#rho");
    Chi_rhoR.set_params({gChi_rho, gV_rho, gT_rho});
    Chi_rhoR.set_formfactor(true, bRho);

    std::vector<amplitude*> chi_exchangesR = {&Chi_omegaR, &Chi_rhoR};
    amplitude_sum chiR(kChi, chi_exchangesR, "#chi_{c1}(1P)");
    
    //////////////////
    // X(3872)
    vector_exchange X_omegaR(kX, &alpha, "#omega");
    X_omegaR.set_params({gX_omega, gV_omega, gT_omega});
    X_omegaR.set_formfactor(true, bOmega);

    vector_exchange X_rhoR(kX, &alpha, "#rho");
    X_rhoR.set_params({gX_rho, gV_rho, gT_rho});
    X_rhoR.set_formfactor(true, bRho);   

    std::vector<amplitude*> X_exchangesR = {&X_omegaR, &X_rhoR};
    amplitude_sum XR(kX, X_exchangesR, "#it{X}(3872)");
    /////////////////
   
    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    ////////////////////////////////////////////

    // amps.push_back(&chi);
    // amps.push_back(&X);

    // int N = 10;

    // double  xmin = 4.;
    // double  xmax = 7.;

    // double  ymin = 2.E-3;
    // double  ymax = 800.;

    // std::string filename = "X_FS.pdf";

    ///////////////////////////////////////////

    amps.push_back(&chiR);
    amps.push_back(&XR);

    // Options
    int N = 20;

    double  xmin = 20.;
    double  xmax = 60.;

    double  ymin = 1.E-5;
    double  ymax = 1.;
    std::string filename = "X_regge.pdf";

    ///////////////////////////////////////////

    std::string ylabel  = ROOT_italics("#sigma(#gamma p #rightarrow X p)") + "   [nb]";

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
        if (xmin < amps[n]->kinematics->Wth)
        {
            double match =  amps[n]->kinematics->Wth + 3.;
            x_fx = vec_fill(20, F, amps[n]->kinematics->Wth + EPS, match, true);
            std::array<std::vector<double>, 2>  x_fx1 = vec_fill(N, F, match, xmax, true);

            for (int i = 0; i < x_fx1[0].size(); i++)
            {
            x_fx[0].push_back(x_fx1[0][i]);
            x_fx[1].push_back(x_fx1[1][i]);
            }

            x_fx[0].insert(x_fx[0].begin(), amps[n]->kinematics->Wth);
            x_fx[1].insert(x_fx[1].begin(), 0.);
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
    
    // plotter->SetLegend(0.23, 0.75);
    plotter->SetLegend(0.73, 0.65);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
    delete kChi, kX;

    return 1.;
};