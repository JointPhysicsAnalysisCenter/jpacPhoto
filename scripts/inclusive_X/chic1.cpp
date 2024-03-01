// Calculates the integrated cross section for both inclusive and exclusive 
// chi_c1 via photon and vector meson exchanges
//
// OUTPUT: inclusive_chic1.pdf
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#include "plotter.hpp"
#include "analytic/vector_exchange.hpp"
#include "covariant/photon_exchange.hpp"
#include "analytic/vector_exchange.hpp"

void chic1()
{
    using namespace jpacPhoto;

    double Wth = M_CHIC1 + M_PROTON;
    
    //----------------------------------------------------------------------------
    // INPUTS

    // Couplings
    double gRho    = 1.140E-3, gOmega  = 0.190E-3, gPhi    = 0.110E-3, gPsi    = 34.64E-3;

    // VMD factors
    double etaRho   = 16.37, etaOmega = 56.34, etaPhi   = 44.37, etaPsi   = 36.85;

    // Form factor cutoffs
    double lamRho   = 1.4, lamOmega = 1.2;

    std::vector<double> parsRho, parsOmega, parsPhi, parsPsi;
    parsRho   = {gRho,   etaRho,   lamRho  };
    parsOmega = {gOmega, etaOmega, lamOmega};
    parsPhi   = {gPhi,   etaPhi,   lamOmega};
    parsPsi   = {gPsi,   etaPsi,   0.};

    //----------------------------------------------------------------------------
    // Set up exclusives amplitudes

    kinematics kChi = new_kinematics(M_CHIC1);
    kChi->set_meson_JP(AXIALVECTOR);

    amplitude exc_omega = new_amplitude<covariant::photon_exchange>(kChi, M_OMEGA, "#omega");
    exc_omega->set_parameters(parsOmega);

    amplitude exc_rho   = new_amplitude<covariant::photon_exchange>(kChi, M_RHO, "#rho");
    exc_rho->set_parameters(parsRho);

    amplitude exc_phi   = new_amplitude<covariant::photon_exchange>(kChi, M_PHI, "#phi");
    exc_phi->set_parameters(parsPhi);

    amplitude exc_psi   = new_amplitude<covariant::photon_exchange>(kChi, M_JPSI, "J/#psi");
    exc_psi->set_parameters(parsPsi);

    // For the neutron target we need to flip the sign of the coupling for the rho
    amplitude exc_rho_m = new_amplitude<covariant::photon_exchange>(kChi, M_RHO, "#minus #rho Exhange");
    exc_rho->set_parameters({-gRho,   etaRho,  lamRho  });

    amplitude exc_mesons_p = exc_rho + exc_omega  + exc_phi/*  + exc_psi */;
    exc_mesons_p->set_id("Exclusive");
    amplitude exc_mesons_n = exc_rho_m + exc_omega + exc_phi/*  + exc_psi */;
    exc_mesons_n->set_id("Exclusive");

    //----------------------------------------------------------------------------
    // Compare with explicitly 

    // Nucleon couplings 
    double gV_omega = 16.,  gT_omega = 0.;
    double gV_rho = 2.4,    gT_rho = 14.6;
    double gV_phi = -6.2,   gT_phi = 2.1;
    double gV_psi = 1.6E-3, gT_psi = 0.;
    
    // Photon couplings
    double gChi_omega   = 5.2E-4;
    double gChi_rho     = 9.2E-4;
    double gChi_phi     = 4.2E-4;
    double gChi_psi     = 1.;
    double gX_omega     = 8.2E-3;
    double gX_rho       = 3.6E-3;

    // chi_c1
    amplitude ChiC1_omegaL = new_amplitude<analytic::vector_exchange>(kChi, M_OMEGA, "#omega exchange");
    ChiC1_omegaL->set_parameters({gChi_omega, gV_omega, gT_omega, lamOmega});

    amplitude ChiC1_rhoL = new_amplitude<analytic::vector_exchange>(kChi, M_RHO, "#rho exchange");
    ChiC1_rhoL->set_parameters({gChi_rho, gV_rho, gT_rho, lamRho});

    amplitude ChiC1_phiL = new_amplitude<analytic::vector_exchange>(kChi, M_PHI, "#phi exchange");
    ChiC1_phiL->set_option(analytic::vector_exchange::kNoFF);
    ChiC1_phiL->set_parameters({gChi_phi, gV_phi, gT_phi});

    amplitude ChiC1_psiL = new_amplitude<analytic::vector_exchange>(kChi, M_JPSI, "#it{J}/#psi exchange");
    ChiC1_psiL->set_option(analytic::vector_exchange::kNoFF);

    ChiC1_psiL->set_parameters({gChi_psi, gV_psi, gT_psi});
    
    amplitude ChiC1_L = ChiC1_rhoL + ChiC1_omegaL + ChiC1_phiL /* + ChiC1_psiL */;
    ChiC1_L->set_id("Sum");

    // --------------------------------------------------------------------------
    // Aux functions to help plotting easier

    // Bounds to plot
    std::array<double,2> NT = {Wth + EPS, 7};
    
    plotter plotter;
    plot p1 = plotter.new_plot();
    
    p1.set_curve_points(100);
    p1.set_logscale(false, true);
    p1.set_ranges({4.35, 7}, {1E-4, 1});
    p1.set_legend(0.65, 0.28);
    p1.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma(#gamma#it{p} #rightarrow #chi_{c1} #it{p})  [nb]");

    std::vector<amplitude> exc_exchanges = extract_subamplitudes(exc_mesons_p);
    std::vector<amplitude> hadron_exc    = extract_subamplitudes(ChiC1_L);
    for (int i = 0; i < exc_exchanges.size(); i++)
    {
        p1.add_curve(  NT, [&](double W){ return exc_exchanges[i]->integrated_xsection(W*W); }, exc_exchanges[i]->id());
        p1.add_dashed( NT, [&](double W){ return hadron_exc[i]->integrated_xsection(W*W); });
    };
    p1.save("chic1.pdf");
};