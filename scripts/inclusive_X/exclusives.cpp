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

void exclusives()
{
    using namespace jpacPhoto;
    using covariant::photon_exchange;
    using analytic::vector_exchange;

    kinematics kX = new_kinematics(M_X3872);
    kX->set_meson_JP(AXIALVECTOR);

    kinematics kC = new_kinematics(M_CHIC1);
    kC->set_meson_JP(AXIALVECTOR);

    // -----------------------------------------
    // X(3872) amplitudes

    amplitude X_gam = new_amplitude<photon_exchange>(kX, "#gamma^{*}");
    X_gam->set_parameters({3.2E-3, 1, 0});

    // VMD rescaled amplitudes 

    amplitude X_rho_VMD   = new_amplitude<photon_exchange>(kX, M_RHO, "#rho");
    X_rho_VMD->set_parameters(  {3.2E-3, 16.37, 1.4});

    amplitude X_omega_VMD = new_amplitude<photon_exchange>(kX, M_OMEGA, "#omega");
    X_omega_VMD->set_parameters({3.2E-3, 56.34, 1.2});

    amplitude X_sum_VMD   = X_rho_VMD + X_omega_VMD;
    X_sum_VMD->set_id("Rescaled");

    // Compare with the explicitly meson exchange amplitudes

    amplitude X_rho   = new_amplitude<vector_exchange>(kX, M_RHO, "#rho exchange");
    X_rho->set_parameters({3.6E-3, 2.4, 14.6, 1.4});

    amplitude X_omega = new_amplitude<vector_exchange>(kX, M_OMEGA, "#omega exchange");
    X_omega->set_parameters({8.2E-3, 16, 0, 1.2});

    amplitude X_sum = X_rho + X_omega;
    X_sum->set_id("Exclusive paper");

    // -----------------------------------------
    // chi_c1 amplitudes

    amplitude C_omega_VMD = new_amplitude<photon_exchange>(kC, M_OMEGA, "#omega");
    C_omega_VMD->set_parameters({0.190E-3, 56.34, 1.2});

    amplitude C_rho_VMD   = new_amplitude<photon_exchange>(kC, M_RHO, "#rho");
    C_rho_VMD->set_parameters({1.140E-3, 16.37, 1.4});

    amplitude C_phi_VMD   = new_amplitude<photon_exchange>(kC, M_PHI, "#phi");
    C_phi_VMD->set_parameters({-0.110E-3, 44.37, 1.2});

    amplitude C_sum_VMD = C_omega_VMD + C_rho_VMD + C_phi_VMD;

    amplitude C_omega = new_amplitude<vector_exchange>(kC, M_OMEGA, "#omega exchange");
    C_omega->set_parameters({5.2E-4, 16., 0., 1.2});

    amplitude C_rho = new_amplitude<vector_exchange>(kC, M_RHO, "#rho exchange");
    C_rho->set_parameters({9.2E-4, 2.4, 14.6, 1.4});

    amplitude C_phi = new_amplitude<vector_exchange>(kC, M_PHI, "#phi exchange");
    C_phi->set_option(vector_exchange::kNoFF);
    C_phi->set_parameters({4.2E-4, -6.2, 2.1});

    amplitude C_sum = C_omega + C_rho + C_phi;

    // --------------------------------------------------------------------------
    // Plot results

    plotter plotter;
    
    // Compare chic1 exchanges
    plot p1 = plotter.new_plot();
    p1.set_curve_points(100);
    p1.set_logscale(false, true);
    p1.set_ranges({4.35, 7}, {5E-4, 1});
    p1.set_legend(0.65, 0.2);
    p1.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma(#gamma#it{p} #rightarrow #chi_{c1} #it{p})  [nb]");
    std::vector<amplitude> exc_exchanges = extract_subamplitudes(C_sum_VMD);
    std::vector<amplitude> hadron_exc    = extract_subamplitudes(C_sum);
    for (int i = 0; i < exc_exchanges.size(); i++)
    {
        p1.add_curve( {kC->Wth() + EPS, 7}, [&](double W){ return exc_exchanges[i]->integrated_xsection(W*W); }, exc_exchanges[i]->id());
        p1.add_dashed({kC->Wth() + EPS, 7}, [&](double W){ return hadron_exc[i]->integrated_xsection(W*W); });
    };
    p1.save("chic1_exchanges_compare.pdf");

    // Compare FF choices 
    plot p2 = plotter.new_plot();
    p2.set_curve_points(100);
    p2.set_logscale(false, true);
    p2.set_ranges({4.35, 7}, {1E-3, 5});
    p2.add_header("#omega exchange");
    p2.set_legend(0.22, 0.7);
    p2.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma(#gamma#it{p} #rightarrow #chi_{c1} #it{p})  [nb]");
    amplitude VMD = C_omega_VMD, EXP = C_omega;
    VMD->set_option(photon_exchange::kUseT);      EXP->set_option(vector_exchange::kUseT);   
    p2.add_curve(  {kC->Wth() + EPS, 7}, [&](double W){ return VMD->integrated_xsection(W*W); }, "#beta(#it{t})");
    p2.add_dashed( {kC->Wth() + EPS, 7}, [&](double W){ return EXP->integrated_xsection(W*W); });
    VMD->set_option(photon_exchange::kUseTprime); EXP->set_option(vector_exchange::kUseTprime);   
    p2.add_curve(  {kC->Wth() + EPS, 7}, [&](double W){ return VMD->integrated_xsection(W*W); }, "#beta(#it{t} - #it{t}_{min})");
    p2.add_dashed( {kC->Wth() + EPS, 7}, [&](double W){ return EXP->integrated_xsection(W*W); });
    p2.save("FF_compare.pdf");

    // X(3872) productions
    plot p3 = plotter.new_plot();
    p3.set_curve_points(50);
    p3.set_logscale(false, true);
    p3.set_ranges({4.35, 7}, {4E-3, 2E2});
    p3.set_legend(0.6, 0.22);
    p3.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma(#gamma#it{p} #rightarrow #it{X}#it{p})  [nb]");
    p3.add_curve( {kC->Wth(), 7}, [&](double W){ return C_sum_VMD->integrated_xsection(W*W); },  "#chi_{c1}(1P)");
    p3.add_dashed({kC->Wth(), 7}, [&](double W){ return C_sum->integrated_xsection(W*W);     });
    p3.add_curve( {kX->Wth(), 7}, [&](double W){ return X_sum_VMD->integrated_xsection(W*W); },  "#it{X}(3872)");
    p3.add_dashed({kX->Wth(), 7}, [&](double W){ return X_sum->integrated_xsection(W*W);     });
    p3.save("X_compare.pdf");
};