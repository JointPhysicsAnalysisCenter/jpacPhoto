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
#include "analytic/vector_exchange.hpp"
#include "covariant/photon_exchange.hpp"

void NT_X3872()
{
    using namespace jpacPhoto;
    using inclusive::vector_exchange;
    using covariant::photon_exchange;
    const int p = inclusive::vector_exchange::kProton;
    const int n = inclusive::vector_exchange::kNeutron;

    //----------------------------------------------------------------------------
    // INPUTS

    // Couplings
    double gXGG  = 3.20E-3;
    double gRho    = 1.140E-3, gOmega  = 0.190E-3, gPhi    = 0.110E-3;

    // VMD factors
    double etaRho   = 16.37, etaOmega = 56.34, etaPhi = 44.37;
    
    // Form factor cutoffs
    double lamRho   = 1.4, lamOmega = 1.2;

    //----------------------------------------------------------------------------
    // chi_c1

    kinematics kC = new_kinematics(M_CHIC1);
    kC->set_meson_JP(AXIALVECTOR);

    amplitude eC_omega = new_amplitude<photon_exchange>(kC, M_OMEGA, "Omega Exchange");
    eC_omega->set_parameters({gOmega, etaOmega, lamOmega});

    amplitude eC_rho   = new_amplitude<photon_exchange>(kC, M_RHO, "Rho Exchange");
    eC_rho->set_parameters({gRho, etaRho, lamRho});

    amplitude eC_phi   = new_amplitude<photon_exchange>(kC, M_PHI, "Phi Exchange");
    eC_phi->set_parameters({gPhi, etaPhi, 0.});

    semi_inclusive iC_omega = new_semi_inclusive<vector_exchange>(kC, M_OMEGA, "Inclusive");
    iC_omega->set_parameters({gOmega, etaOmega, lamOmega});

    semi_inclusive iC_rho   = new_semi_inclusive<vector_exchange>(kC, M_RHO, "Inclusive");
    iC_rho->set_parameters({gRho, etaRho, lamRho});
    
    semi_inclusive iC_phi   = new_semi_inclusive<vector_exchange>(kC, M_PHI, "Inclusive");
    iC_phi->set_parameters({gPhi, etaPhi, 0});

    amplitude      eC = eC_omega + eC_rho + eC_phi;
    semi_inclusive iC = iC_omega + iC_rho + iC_phi + eC;

    //----------------------------------------------------------------------------
    // X(3872)

    kinematics kX = new_kinematics(M_X3872);
    kX->set_meson_JP(AXIALVECTOR);

    amplitude eX_omega = new_amplitude<photon_exchange>(kX, M_OMEGA, "Omega Exchange");
    eX_omega->set_parameters({gXGG, etaOmega, lamOmega});

    amplitude eX_rho   = new_amplitude<photon_exchange>(kX, M_RHO, "Rho Exchange");
    eX_rho->set_parameters({gXGG, etaRho, lamRho});

    semi_inclusive iX_omega = new_semi_inclusive<vector_exchange>(kX, M_OMEGA, "Inclusive");
    iX_omega->set_parameters({gXGG, etaOmega, lamOmega});

    semi_inclusive iX_rho   = new_semi_inclusive<vector_exchange>(kX, M_RHO, "Inclusive");
    iX_rho->set_parameters({gXGG, etaRho, lamRho});

    amplitude      eX = eX_omega + eX_rho;
    semi_inclusive iX = iX_omega + iX_rho + eX;

    // --------------------------------------------------------------------------
    // Plot results

    // Bounds to plot
    std::array<double,2> X_NT = {kX->Wth() + EPS, 7.0};
    std::array<double,2> C_NT = {kC->Wth() + EPS, 7.0};

    plotter plotter;
    plot p1 = plotter.new_plot();

    p1.set_curve_points(30);
    p1.set_logscale(false, true);
    p1.set_ranges({4.6, 7}, {1E-2, 2E3});
    p1.set_legend(0.27, 0.72);
    p1.set_labels( "#it{W}_{#gamma#it{N}}  [GeV]", "#sigma  [nb]");
    p1.print_to_terminal(true);
    
    iX_rho->set_option(p); iX_omega->set_option(p); 
    iC_rho->set_option(p); iC_omega->set_option(p); iC_phi->set_option(p); 
    print("chic1 (p inc)"); divider(2);
    p1.add_curve( C_NT, [&](double W){ return iC->integrated_xsection(W*W, 0.8); }, "#chi_{#it{c}1} (#it{p})");
    print("chic1 (p exc)"); divider(2);
    p1.add_dashed(C_NT, [&](double W){ return eC->integrated_xsection(W*W); });

    iC_rho->set_option(n); iC_omega->set_option(n); iC_phi->set_option(n); 
    eC_rho->set_option(n); eC_omega->set_option(n); eC_phi->set_option(n);
    eC_rho->set_parameters({-gRho,   etaRho,   lamRho});

    print("chic1 (n inc)"); divider(2);
    p1.add_curve( C_NT, [&](double W){ return iC->integrated_xsection(W*W, 0.8); }, "#chi_{#it{c}1} (#it{n})");
    print("chic1 (n exc)"); divider(2);
    p1.add_dashed(C_NT, [&](double W){ return eC->integrated_xsection(W*W); });

    print("X (p inc)"); divider(2);
    p1.add_curve( X_NT, [&](double W){ return iX->integrated_xsection(W*W, 0.8); }, "#it{X}(3872) (#it{p})");
    print("X (p exc)"); divider(2);
    p1.add_dashed(X_NT, [&](double W){ return eX->integrated_xsection(W*W); });

    iX_rho->set_option(n); iX_omega->set_option(n);
    eX_rho->set_option(n); eX_omega->set_option(n);
    eX_rho->set_parameters({-gXGG,   etaRho,   lamRho});

    print("X (n inc)"); divider(2);
    p1.add_curve( X_NT, [&](double W){ return iX->integrated_xsection(W*W, 0.8); }, "#it{X}(3872) (#it{n})");
    print("X (n exc)"); divider(2);
    p1.add_dashed(X_NT, [&](double W){ return eX->integrated_xsection(W*W); });

    p1.shade_region({W_cm(22), 10}, {kBlack, 1001});

    p1.save("NT_X.pdf");
};