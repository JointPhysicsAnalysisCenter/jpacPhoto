// Exploration of the photon exchange amplitude for the X(3872)
// Produces the raw primakoff production rates for both proton and neutron targets
// Also uses VMD to relate this to the meson exchanges previously considered
//
// OUTPUT: primakoff_compare.pdf 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#include "plotter.hpp"
#include "covariant/photon_exchange.hpp"
#include "analytic/vector_exchange.hpp"

void X_compare()
{
    using namespace jpacPhoto;
    using namespace covariant;

    kinematics kX = new_kinematics(M_X3872);
    kX->set_meson_JP(AXIALVECTOR);

    // -----------------------------------------
    // Pure Primakoff

    amplitude primakoff = new_amplitude<photon_exchange>(kX, "#gamma^{*}");
    primakoff->set_parameters({3.2E-3, 1, 0});
    
    // cross section in femto
    auto sig_fb = [&] (double w)
    {
        return primakoff->integrated_xsection(w*w) * 1E6;
    };

    //-----------------------------------------
    // Rescaling the Primakoff with VMD to get the meson exchange

    amplitude rho_VMD   = new_amplitude<photon_exchange>(kX, M_RHO, "#rho");
    rho_VMD->set_parameters(  {3.2E-3, 16.37, 1.4});

    amplitude omega_VMD = new_amplitude<photon_exchange>(kX, M_OMEGA, "#omega");
    omega_VMD->set_parameters({3.2E-3, 56.34, 1.2});

    amplitude sum_VMD   = rho_VMD + omega_VMD;
    sum_VMD->set_id("Rescaled");

    //-----------------------------------------
    // Compare with the explicitly meson exchange amplitudes

    amplitude rho   = new_amplitude<analytic::vector_exchange>(kX, M_RHO, "#rho exchange");
    rho->set_parameters({3.6E-3, 2.4, 14.6, 1.4});

    amplitude omega = new_amplitude<analytic::vector_exchange>(kX, M_OMEGA, "#omega exchange");
    omega->set_parameters({8.2E-3, 16, 0, 1.2});

    amplitude sum = rho + omega;
    sum->set_id("Exclusive paper");

    // --------------------------------------------------------------------------
    // Plot results

    // double W = 6;
    // double costheta = 0.9;
    // double s = W*W;
    // double t = kX->t_man(s, acos(costheta));
    // amplitude amp = rho_VMD;

    // print("Inputs");
    // divider(2);
    // print("W", W);
    // print("cosTheta", costheta);
    // print("mX3872", M_X3872);
    // print("mProton", M_PROTON);
    // print("gXGamGam", 3.2E-3);
    // print("Lambda_omega", 1.2);
    // print("mOmega", M_OMEGA);
    // print("eta_omega", 56.34);
    // print("mRho", M_RHO);
    // print("Lambda_rho", 1.4);
    // print("eta_rho", 16.37);
    // line();
    // print("Outputs");
    // divider(2);
    // print("s ", s);
    // print("t ", t);
    // print("t'", t - kX->t_min(s));
    // print("dsig/dt", amp->differential_xsection(s,t));
    // print("sig",     amp->integrated_xsection(s));
    // exit(1);

    plotter plotter;
    
    // Plot of pure primakoff in fb
    plot p1 = plotter.new_plot();
    p1.set_curve_points(100);
    p1.set_logscale(false, true);
    p1.set_ranges({4.6, 7}, {5E-4, 4E1});
    p1.set_legend(0.27, 0.72);
    p1.set_labels( "#it{W}  [GeV]", "#sigma(#gamma#it{N} #rightarrow #it{X}#it{N})  [fb]");
    p1.add_header("#gamma^{*} exchange");

    p1.add_curve(  {kX->Wth(), 7}, sig_fb, "#it{p}");
    primakoff->set_option(  photon_exchange::kNeutron);
    p1.add_curve(  {kX->Wth(), 7}, sig_fb, "#it{n}");

    // Plot of vector meson exchanges in nb 
    plot p2 = plotter.new_plot();
    p2.color_offset(2);
    p2.set_curve_points(100);
    p2.set_logscale(false, true);
    p2.set_ranges({4.6, 7}, {6E-2, 4E2});
    p2.set_legend(0.7, 0.22);
    p2.set_labels( "#it{W}  [GeV]", "#sigma(#gamma#it{N} #rightarrow #it{X}#it{N})  [nb]");
    p2.add_header("#rho /#omega exchange");

    sum_VMD->set_id("#it{p}");
    p2.add_curve(  sigma_w, sum_VMD, {kX->Wth(), 7});
    p2.add_dashed( sigma_w, sum,     {kX->Wth(), 7});

    sum_VMD->set_option(  photon_exchange::kNeutron);
    sum_VMD->set_id("#it{n}");
    omega_VMD->set_parameters({-3.2E-3, 56.34, 1.2}); // For the neutron, the omega term gains an overall minus sign
    p2.add_curve(  sigma_w, sum_VMD, {kX->Wth(), 7});

    plotter.combine({2,1}, {p1, p2}, "primakoff_compare.pdf");
};