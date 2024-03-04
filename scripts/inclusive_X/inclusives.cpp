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

void inclusives()
{
    using namespace jpacPhoto;
    using inclusive::vector_exchange;
    using covariant::photon_exchange;

    //----------------------------------------------------------------------------
    // INPUTS

    // Couplings
    double gXGG   = 3.20E-3;
    double gRho   = 1.140E-3, gOmega   = 0.190E-3, gPhi   = 0.110E-3, gGam = (1.14+0.19+0.11+ 34.65)*1E-3;
    double etaRho = 16.37,    etaOmega = 56.34,    etaPhi = 44.37;
    
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

    amplitude eC_gam = new_amplitude<photon_exchange>(kC, 0., "#gamma exchange");
    eC_gam->set_parameters({gGam, 1, 0.});

    semi_inclusive iC_omega = new_semi_inclusive<vector_exchange>(kC, M_OMEGA, "Inclusive");
    iC_omega->set_parameters({gOmega, etaOmega, lamOmega});

    semi_inclusive iC_rho   = new_semi_inclusive<vector_exchange>(kC, M_RHO, "Inclusive");
    iC_rho->set_parameters({gRho, etaRho, lamRho});
    
    semi_inclusive iC_phi   = new_semi_inclusive<vector_exchange>(kC, M_PHI, "Inclusive");
    iC_phi->set_parameters({-gPhi, etaPhi, lamOmega});

    semi_inclusive iC_gam   = new_semi_inclusive<vector_exchange>(kC, 0, "Inclusive");
    iC_gam->reggeized(true);
    iC_gam->set_parameters({gGam, 1, 0.});

    amplitude      eC = eC_omega + eC_rho + eC_phi;
    semi_inclusive iC = iC_omega + iC_rho + iC_phi + iC_gam + eC;

    // Reggeized exclusive amplitude
    amplitude rC_omega = new_amplitude<regge::vector_exchange>(kC, "#omega exchange");
    rC_omega->set_parameters({5.2E-4, 16., 0., 1.2, 0.5, 0.9});

    amplitude rC_rho = new_amplitude<regge::vector_exchange>(kC, "#rho exchange");
    rC_rho->set_parameters({9.2E-4, 2.4, 14.6, 1.4, 0.5, 0.9});

    amplitude rC = rC_omega + rC_rho;
    semi_inclusive irC = iC_omega + iC_rho + iC_gam + rC;

    //---------------------------------------------------------------------------
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

    semi_inclusive iX_gam = new_semi_inclusive<inclusive::vector_exchange>(kX, 0, "Inclusive");
    iX_gam->reggeized(true);
    iX_gam->set_parameters({gXGG, 1, 0});

    amplitude      eX = eX_omega + eX_rho;
    semi_inclusive iX = iX_omega + iX_rho + eX;

    // Reggeized exclusive amplitude 
    amplitude rX_omega = new_amplitude<regge::vector_exchange>(kX, "#omega exchange");
    rX_omega->set_parameters({8.2E-3, 16., 0., lamOmega, 0.5, 0.9});

    amplitude rX_rho   = new_amplitude<regge::vector_exchange>(kX, "#rho exchange");
    rX_rho->set_parameters({3.6E-3, 2.4, 14.6, lamRho, 0.5, 0.9});

    amplitude rX = rX_omega + rX_rho;
    semi_inclusive irX = iX_omega + iX_rho + rX;

    // --------------------------------------------------------------------------
    // Plot results

    // Bounds to plot
    std::array<double,2> X_NT = {kX->Wth() + EPS, 7.0};
    std::array<double,2> C_NT = {kC->Wth() + EPS, 7.0};
    std::array<double,2> HE   = {20, 60};

    plotter plotter;

    plot p1 = plotter.new_plot();
    p1.set_curve_points(30);
    p1.set_logscale(false, true);
    p1.set_ranges({4.1, 7}, {1E-2, 2E3});
    p1.set_legend(0.27, 0.72);
    p1.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma  [nb]");
    p1.print_to_terminal(true);
    p1.shade_region({W_cm(22), 10}, {kBlack, 1001});
    print("chic1 (inclusive)"); divider(2);
    p1.add_curve( C_NT, [&](double W){ return iC->integrated_xsection(W*W, 0.8); }, "#chi_{#it{c}1}");
    print("chic1 (exclusive)"); divider(2);
    p1.add_dashed(C_NT, [&](double W){ return eC->integrated_xsection(W*W); });
    print("X(3872) (inclusive)"); divider(2);
    p1.add_curve( X_NT, [&](double W){ return iX->integrated_xsection(W*W, 0.8); }, "#it{X}(3872)");
    print("X(3872) (exclusive)"); divider(2);
    p1.add_dashed(X_NT, [&](double W){ return eX->integrated_xsection(W*W); });
    p1.save("NT.pdf");

    plot p2 = plotter.new_plot();
    p2.set_curve_points(30);
    p2.set_logscale(false, true);
    p2.set_ranges(HE, {5E-5, 0.5});
    p2.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma  [nb]");
    p2.set_legend(0.22, 0.30);
    iX_omega->reggeized(true); iX_rho->reggeized(true);
    p2.add_curve( HE, [&](double W){ return irX->integrated_xsection(W*W); },                   "Total #it{X}(3872) production");
    p2.add_curve( HE, [&](double W){ return (iX_omega+iX_rho)->integrated_xsection(W*W); },     "Inclusive #it{V} exchange");
    p2.add_curve( HE, [&](double W){ return rX->integrated_xsection(W*W); },                    "Exclusive #it{V} exchange");
    p2.add_curve( HE, [&](double W){ return iX_gam->integrated_xsection(W*W); },                "Inclusive #gamma exchange");

    plot p3 = plotter.new_plot();
    p3.set_curve_points(30);
    p3.set_logscale(false, true);
    p3.set_ranges(HE, {5E-4, 3E2});
    p3.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma  [pb]");
    p3.set_legend(0.20, 0.17);
    iC_omega->reggeized(true); iC_rho->reggeized(true); 
    p3.print_to_terminal(true);
    p3.add_curve( HE, [&](double W){ return (irC->integrated_xsection(W*W)+eC_gam->integrated_xsection(W*W)) * 1E3; }, "Total #chi_{#it{c}1}(1P) production");
    p3.add_curve(HE, [&](double W){  return (iC_rho+iC_omega)->integrated_xsection(W*W)  * 1E3; }, "Inclusive #it{V} exchange");
    p3.add_curve(HE, [&](double W){  return rC->integrated_xsection(W*W)  * 1E3; },                "Exclusive #it{V} exchange");
    p3.add_curve(HE, [&](double W){  return iC_gam->integrated_xsection(W*W)  * 1E3; },            "Inclusive #gamma exchange");
    p3.add_curve(HE, [&](double W){  return eC_gam->integrated_xsection(W*W)  * 1E3; },            "Exclusive #gamma exchange");
    plotter.combine({2,1}, {p3,p2}, "HE.pdf");

};