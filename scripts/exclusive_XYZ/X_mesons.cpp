// Photoproduction cross-sections for X(3872) and axial-vector charmonium mesons
// Reproduces fig. 3 of [1] in pdf form.
//
// OUTPUT: X_mesons.pdf
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:2008.01001 [hep-ph]
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "kinematics.hpp"
#include "plotter.hpp"

#include "analytic/vector_exchange.hpp"
#include "covariant/vector_exchange.hpp"
#include "regge/vector_exchange.hpp"

void X_mesons()
{
    using namespace jpacPhoto;

    // ---------------------------------------------------------------------------
    // Couplings and constants
    // ---------------------------------------------------------------------------'

    // Nucleon couplings 
    double gV_omega = 16., gT_omega = 0.;
    double gV_rho = 2.4,  gT_rho = 14.6;
    double gV_phi = -6.2, gT_phi = 2.1;
    double gV_psi = 1.6E-3, gT_psi = 0.;
    
    // Photon couplings
    double gChi_omega   = 5.2E-4;
    double gChi_rho     = 9.2E-4;
    double gChi_phi     = 4.2E-4;
    double gChi_psi     = 1.;
    double gX_omega     = 8.2E-3;
    double gX_rho       = 3.6E-3;
    
    // Form factor cutoffs
    double LamOmega = 1.2;
    double LamRho   = 1.4; 

    // Reggeon trajectory
    double inter = 0.5;
    double slope = 0.9;

    // ---------------------------------------------------------------------------
    // Kinematics
    // ---------------------------------------------------------------------------

    kinematics kChiC1  = new_kinematics( M_CHIC1 );
    kChiC1->set_meson_JP(AXIALVECTOR);

    kinematics kX      = new_kinematics( M_X3872 );
    kX->set_meson_JP(AXIALVECTOR);

    // ---------------------------------------------------------------------------
    // Low Energy amplitudes (Fixed-spin vector exchange)
    // ---------------------------------------------------------------------------

    // chi_c1
    amplitude ChiC1_omegaL = new_amplitude<analytic::vector_exchange>(kChiC1, M_OMEGA, "#omega exchange");
    ChiC1_omegaL->set_parameters({gChi_omega, gV_omega, gT_omega, LamOmega});

    amplitude ChiC1_rhoL = new_amplitude<analytic::vector_exchange>(kChiC1, M_RHO, "#rho exchange");
    ChiC1_rhoL->set_parameters({gChi_rho, gV_rho, gT_rho, LamRho});

    amplitude ChiC1_phiL = new_amplitude<analytic::vector_exchange>(kChiC1, M_PHI, "#phi exchange");
    ChiC1_phiL->set_option(analytic::vector_exchange::kNoFF);
    ChiC1_phiL->set_parameters({gChi_phi, gV_phi, gT_phi});

    amplitude ChiC1_psiL = new_amplitude<analytic::vector_exchange>(kChiC1, M_JPSI, "#it{J}/#psi exchange");
    ChiC1_psiL->set_option(analytic::vector_exchange::kNoFF);

    ChiC1_psiL->set_parameters({gChi_psi, gV_psi, gT_psi});
    
    amplitude ChiC1_L = ChiC1_omegaL + ChiC1_rhoL + ChiC1_phiL + ChiC1_psiL;
    ChiC1_L->set_id("#chi_{c1}");

    // X(3872)
    amplitude X_omegaL = new_amplitude<covariant::vector_exchange>(kX, M_OMEGA, "#omega exchange");
    X_omegaL->set_parameters({gX_omega, gV_omega, gT_omega, LamOmega});

    amplitude X_rhoL = new_amplitude<covariant::vector_exchange>(kX, M_RHO, "#rho exchange");
    X_rhoL->set_parameters({gX_rho, gV_rho, gT_rho, LamRho});
    
    // Total is the sum of the above exchanges
    amplitude X_L = X_omegaL + X_rhoL;
    X_L->set_id("#it{X}(3872)");

    // ---------------------------------------------------------------------------
    // High Energy amplitudes (Reggeized vector exchange)
    // ---------------------------------------------------------------------------

    // chi_c1
    amplitude ChiC1_omegaH = new_amplitude<regge::vector_exchange>(kChiC1, "#omega exchange");
    ChiC1_omegaH->set_parameters({gChi_omega, gV_omega, gT_omega, LamOmega, inter, slope});

    amplitude ChiC1_rhoH = new_amplitude<regge::vector_exchange>(kChiC1,  "#rho exchange");
    ChiC1_rhoH->set_parameters({gChi_rho, gV_rho, gT_rho, LamRho, inter, slope});

    amplitude ChiC1_H = ChiC1_omegaH + ChiC1_rhoH;
    ChiC1_H->set_id("#chi_{c1}");

    // X(3872)
    amplitude X_omegaH = new_amplitude<regge::vector_exchange>(kX, "#omega exchange");
    X_omegaH->set_parameters({gX_omega, gV_omega, gT_omega, LamOmega, inter, slope});

    amplitude X_rhoH = new_amplitude<regge::vector_exchange>(kX, "#rho exchange");
    X_rhoH->set_parameters({gX_rho, gV_rho, gT_rho, LamRho, inter, slope});
    
    // Total is the sum of the above exchanges
    amplitude X_H = X_omegaH + X_rhoH;
    X_H->set_id("#it{X}(3872)");

    // ---------------------------------------------------------------------------
    // Plot results
    // ---------------------------------------------------------------------------

    plotter plotter;

    // Low-energy plot
    plot p_low = plotter.new_plot();
    p_low.set_logscale(false, true);
    p_low.set_curve_points(100);
    p_low.set_legend(0.3, 0.73);
    p_low.set_ranges({4, 7}, {2E-3, 800});
    p_low.set_labels( "#it{W_{#gammap}}  [GeV]", "#sigma(#gamma#it{p} #rightarrow #it{Xp})  [nb]");
    p_low.add_curve(sigma_w, {ChiC1_L, X_L}, {4, 7});

    // High-energy plot
    plot p_high = plotter.new_plot();
    p_high.set_logscale(false, true);
    p_high.set_curve_points(100);
    p_high.set_legend(0.3, 0.73);
    p_high.set_ranges({20, 60}, {1E-5, 1});
    p_high.set_labels( "#it{W_{#gammap}}  [GeV]", "#sigma(#gamma#it{p} #rightarrow #it{Xp})  [nb]");
    p_high.add_curve(sigma_w, {ChiC1_H, X_H}, {20, 60});

    // Combine these into a single plot 
    plotter.combine({2,1}, {p_low, p_high}, "X_mesons.pdf");
};