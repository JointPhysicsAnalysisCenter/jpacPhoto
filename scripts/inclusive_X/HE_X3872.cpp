// Calculates the integrated cross section for both inclusive and exclusive 
// X(3872) via photon and vector meson exchanges
//
// OUTPUT: inclusive_X.pdf
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#include "plotter.hpp"
#include "semi_inclusive/vector_exchange.hpp"
#include "regge/vector_exchange.hpp"
#include "covariant/photon_exchange.hpp"

void HE_X3872()
{
    using namespace jpacPhoto;

    const int p = inclusive::vector_exchange::kProton;
    const int n = inclusive::vector_exchange::kNeutron;

    double Wth = M_X3872 + M_PROTON;
    
    //----------------------------------------------------------------------------
    // INPUTS

    // Couplings
    double gGamma  = 3.20E-3;
    double gRho    = 5.38E-3;
    double gOmega  = 3.54E-3;

    // VMD factors
    double etaRho   = 16.37;
    double etaOmega = 56.34;

    // Form factor cutoffs
    double lamRho   = 1.4;
    double lamOmega = 1.2;

    std::vector<double> parsGamma, parsRho, parsOmega;
    parsGamma = {gGamma, 1,        0       };
    parsRho   = {gRho,   etaRho,   lamRho  };
    parsOmega = {gOmega, etaOmega, lamOmega};

    // Nucleon couplings 
    double gV_omega = 16.,  gT_omega = 0.;
    double gV_rho = 2.4,    gT_rho = 14.6;
    double gV_phi = -6.2,   gT_phi = 2.1;
    double gV_psi = 1.6E-3, gT_psi = 0.;
    
    // Photon couplings
    double gX_omega     = 8.2E-3;
    double gX_rho       = 3.6E-3;
    
    // Reggeon trajectory
    double inter = 0.5;
    double slope = 0.9;

    kinematics kX = new_kinematics(M_X3872);
    kX->set_meson_JP(AXIALVECTOR);

    //----------------------------------------------------------------------------
    // Set up inclusive amplitudes

    semi_inclusive inc_gamma = new_semi_inclusive<inclusive::vector_exchange>(kX, 0, "Inclusive");
    inc_gamma->reggeized(true);
    inc_gamma->set_parameters(parsGamma);

    semi_inclusive inc_omega = new_semi_inclusive<inclusive::vector_exchange>(kX, M_OMEGA, "Inclusive");
    inc_omega->reggeized(true);
    inc_omega->set_parameters(parsOmega);

    semi_inclusive inc_rho   = new_semi_inclusive<inclusive::vector_exchange>(kX, M_RHO, "Inclusive");
    inc_rho->reggeized(true);
    inc_rho->set_parameters(parsRho);

    semi_inclusive inc_mesons = inc_omega + inc_rho;

    //----------------------------------------------------------------------------
    // Set up exclusives amplitudes

    amplitude exc_gamma = new_amplitude<covariant::photon_exchange>(kX, "#gamma^{*} exchange");
    exc_gamma->set_parameters(parsGamma);
    
    // X(3872)
    amplitude exc_omega = new_amplitude<regge::vector_exchange>(kX, "#omega exchange");
    exc_omega->set_parameters({gX_omega, gV_omega, gT_omega, lamOmega, inter, slope});

    amplitude exc_rho   = new_amplitude<regge::vector_exchange>(kX, "#rho exchange");
    exc_rho->set_parameters({gX_rho, gV_rho, gT_rho, lamRho, inter, slope});

    amplitude exc_mesons = exc_omega + exc_rho;
    exc_mesons->set_id("Exclusive");

    // --------------------------------------------------------------------------
    // Plot results

    // Bounds to plot
    std::array<double,2> HE = {20, 60};

    plotter plotter;
    plot p1 = plotter.new_plot();

    p1.set_curve_points(40);
    p1.set_logscale(false, true);
    p1.set_ranges(HE, {1E-3, 5});
    p1.set_legend(0.22, 0.75);
    p1.set_labels( "#it{W}_{#gamma#it{N}}  [GeV]", "#sigma  [nb]");

    p1.add_curve( HE, [&](double W){ return inc_gamma->integrated_xsection(W*W); }, "#gamma");
    // p1.add_dashed(HE, exc_primakoff_1E3);
    p1.add_curve( HE, [&](double W){ return inc_mesons->integrated_xsection(W*W); }, "#rho / #omega");
    p1.add_dashed(HE, [&](double W){ return exc_mesons->integrated_xsection(W*W); });

    p1.save("regge_X.pdf");
};