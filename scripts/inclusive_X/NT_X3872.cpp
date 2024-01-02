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
#include "inclusive/vector_exchange.hpp"
#include "analytic/vector_exchange.hpp"
#include "covariant/photon_exchange.hpp"

void NT_X3872()
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

    kinematics kX = new_kinematics(M_X3872);
    kX->set_meson_JP(AXIALVECTOR);

    //----------------------------------------------------------------------------
    // Set up inclusive amplitudes

    semi_inclusive inc_omega = new_semi_inclusive<inclusive::vector_exchange>(kX, M_OMEGA, "Inclusive");
    inc_omega->set_parameters(parsOmega);

    semi_inclusive inc_rho   = new_semi_inclusive<inclusive::vector_exchange>(kX, M_RHO, "Inclusive");
    inc_rho->set_parameters(parsRho);

    //----------------------------------------------------------------------------
    // Set up exclusives amplitudes

    amplitude exc_omega = new_amplitude<covariant::photon_exchange>(kX, M_OMEGA, "Omega Exchange");
    exc_omega->set_parameters(parsOmega);

    amplitude exc_rho   = new_amplitude<covariant::photon_exchange>(kX, M_RHO, "Rho Exchange");
    exc_rho->set_parameters(parsRho);

    // For the neutron target we need to flip the sign of the coupling for the rho
    amplitude exc_omega_n = new_amplitude<covariant::photon_exchange>(kX, M_OMEGA, "Omega Exchange");
    exc_omega_n->set_parameters(parsOmega);
    exc_omega_n->set_option(n);
    amplitude exc_rho_n = new_amplitude<covariant::photon_exchange>(kX, M_RHO, "#minus Rho Exhange");
    exc_rho_n->set_parameters({-gRho,   etaRho,   lamRho  });
    exc_rho_n->set_option(n);

    amplitude exc_mesons_p = exc_omega + exc_rho;
    exc_mesons_p->set_id("Exclusive");

    amplitude exc_mesons_n = exc_omega_n + exc_rho_n;
    exc_mesons_n->set_id("Exclusive");

    semi_inclusive inc_p = inc_omega + inc_rho + exc_mesons_p;
    semi_inclusive inc_n = inc_omega + inc_rho + exc_mesons_n;

    // --------------------------------------------------------------------------
    // Plot results

    // Bounds to plot
    std::array<double,2> NT = {Wth, 7.0};

    plotter plotter;
    plot p1 = plotter.new_plot();

    p1.set_curve_points(30);
    p1.set_logscale(false, true);
    p1.set_ranges({4.6, 7}, {10E-2, 2E3});
    p1.set_legend(0.27, 0.72);
    p1.set_labels( "#it{W}_{#gamma#it{N}}  [GeV]", "#sigma  [nb]");
    
    inc_rho->set_option(p); inc_omega->set_option(p); 
    p1.add_curve( NT, [&](double W){ return inc_p->integrated_xsection(W*W, 0.8); }, "#it{p}");
    p1.add_dashed(NT, [&](double W){ return exc_mesons_p->integrated_xsection(W*W); });

    inc_rho->set_option(n); inc_omega->set_option(n);
    p1.add_curve( NT, [&](double W){ return inc_n->integrated_xsection(W*W, 0.8); }, "#it{n}");
    p1.add_dashed(NT, [&](double W){ return exc_mesons_n->integrated_xsection(W*W); });

    p1.save("NT_X.pdf");
};