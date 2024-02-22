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

void chic1_FFs()
{
    using namespace jpacPhoto;
    using analytic::vector_exchange;
    using covariant::photon_exchange;

    double Wth = M_CHIC1 + M_PROTON;
    
    //----------------------------------------------------------------------------
    // INPUTS

    // Couplings
    double gOmega  = 0.190E-3;

    // VMD factors
    double etaOmega = 56.34;

    // Form factor cutoffs
    double lamOmega = 1.2;


    //----------------------------------------------------------------------------
    // Set up exclusives amplitudes

    kinematics kChi = new_kinematics(M_CHIC1);
    kChi->set_meson_JP(AXIALVECTOR);

    amplitude VMD = new_amplitude<photon_exchange>(kChi, M_OMEGA, "#omega");
    VMD->set_parameters({gOmega, etaOmega, lamOmega});


    //----------------------------------------------------------------------------
    // Compare with explicitly 

    // Nucleon couplings 
    double gV_omega = 16.,  gT_omega = 0.;
    
    // Photon couplings
    double gChi_omega   = 5.2E-4;

    // chi_c1
    amplitude EXP = new_amplitude<vector_exchange>(kChi, M_OMEGA, "#omega exchange");
    EXP->set_parameters({gChi_omega, gV_omega, gT_omega, lamOmega});


    // --------------------------------------------------------------------------
    // Aux functions to help plotting easier

    // Bounds to plot
    std::array<double,2> NT = {Wth + EPS, 7};
    
    plotter plotter;
    plot p1 = plotter.new_plot();
    
    p1.set_curve_points(100);
    p1.set_logscale(false, true);
    p1.set_ranges({4.35, 7}, {1E-3, 5});
    p1.add_header("#omega exchange");
    p1.set_legend(0.22, 0.7);
    p1.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma(#gamma#it{p} #rightarrow #chi_{c1} #it{p})  [nb]");
    // p1.print_to_terminal(true);

    VMD->set_option(photon_exchange::kUseT);      EXP->set_option(vector_exchange::kUseT);   
    p1.add_curve(  NT, [&](double W){ return VMD->integrated_xsection(W*W); }, "#beta(#it{t})");
    p1.add_dashed( NT, [&](double W){ return EXP->integrated_xsection(W*W); });
    VMD->set_option(photon_exchange::kUseTprime); EXP->set_option(vector_exchange::kUseTprime);   
    p1.add_curve(  NT, [&](double W){ return VMD->integrated_xsection(W*W); }, "#beta(#it{t} - #it{t}_{min})");
    p1.add_dashed( NT, [&](double W){ return EXP->integrated_xsection(W*W); });
    p1.save("FF_compare.pdf");
};