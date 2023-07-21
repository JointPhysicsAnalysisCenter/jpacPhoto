// Exploration of the photon exchange amplitude for the X(3872)
// Produces the raw primakoff production rates for both proton and neutron targets
// Also uses VMD to relate this to the meson exchanges previously considered
//
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

void elastic_primakoff()
{
    using namespace jpacPhoto;
    using namespace covariant;

    kinematics kX = new_kinematics(M_X3872);
    kX->set_meson_JP(AXIALVECTOR);

    // -----------------------------------------
    // Pure Primakoff
    amplitude primakoff = new_amplitude<photon_exchange>(kX, "#gamma^{*}");
    primakoff->set_parameters({3.2E-3, 1, 0});
    
    auto sig_pb = [&] (double w)
    {
        return primakoff->integrated_xsection(w*w) * 1E3;
    };

    //-----------------------------------------
    // Rescaling the Primakoff with VMD to get the meson exchange

    amplitude rho_VMD   = new_amplitude<photon_exchange>(kX, "#rho");
    rho_VMD->set_parameters(  {3.2E-3, 16.37, 0.77});

    amplitude omega_VMD = new_amplitude<photon_exchange>(kX, "#omega");
    omega_VMD->set_parameters({3.2E-3, 56.34, 0.78});

    amplitude sum_VMD   = rho_VMD + omega_VMD;
    sum_VMD->set_id("Rescaled");

    //-----------------------------------------
    // Compare with the explicitly meson exchange amplitudes

    amplitude omega = new_amplitude<analytic::vector_exchange>(kX, M_OMEGA, "#omega exchange");
    omega->set_parameters({8.2E-3, 16, 0, 1.2});

    amplitude rho   = new_amplitude<analytic::vector_exchange>(kX, M_RHO, "#rho exchange");
    rho->set_parameters({3.6E-3, 2.4, 14.6, 1.4});

    amplitude sum = rho + omega;
    sum->set_id("Explicit");

    // --------------------------------------------------------------------------
    // Plot results

    plotter plotter;

    plot p1 = plotter.new_plot();
    p1.set_curve_points(100);
    p1.set_logscale(false, true);
    p1.set_ranges({4.6, 7}, {1E-5, 0.2});
    p1.set_legend(0.27, 0.72);
    p1.set_labels( "#it{W}  [GeV]", "#sigma(#gamma#it{N} #rightarrow #it{X}#it{N})  [pb]");
    p1.add_header("#gamma^{*} exchange");

    p1.add_curve(  {kX->Wth(), 7}, sig_pb, "#it{p}");
    primakoff->set_option(  photon_exchange::kNeutron);
    p1.add_curve(  {kX->Wth(), 7}, sig_pb, "#it{n}");

    plot p2 = plotter.new_plot();
    p2.color_offset(2);
    p2.set_curve_points(100);
    p2.set_logscale(false, true);
    p2.set_ranges({4.6, 7}, {2E-1, 4E2});
    p2.set_legend(0.7, 0.22);
    p2.set_labels( "#it{W}  [GeV]", "#sigma(#gamma#it{N} #rightarrow #it{X}#it{N})  [nb]");
    p2.add_header("#rho /#omega exchange");

    sum_VMD->set_id("#it{p}");
    p2.add_curve(  sigma_w, sum_VMD, {kX->Wth(), 7});
    p2.add_dashed( sigma_w, sum,     {kX->Wth(), 7});

    sum_VMD->set_option(  photon_exchange::kNeutron);
    sum_VMD->set_id("#it{n}");
    p2.add_curve(  sigma_w, sum_VMD, {kX->Wth(), 7});

    plotter.combine({2,1}, {p1, p2}, "sigs.pdf");
};