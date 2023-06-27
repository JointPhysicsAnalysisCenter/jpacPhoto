// Plots sigma_L and sigma_T at fixed Q2 as a function of W using the low-energy
// empirical formulae of [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] - https://arxiv.org/abs/0712.3731
// ------------------------------------------------------------------------------

#include "sigma_tot/CB_F.hpp"
#include "plotter.hpp"

void CB_sigmaLT()
{
    using namespace jpacPhoto;

    auto CB = dynamic_pointer_cast<CB_F>(new_total_xsection<CB_F>(1)); 

    // Pick some fixed Q2
    double Q2 = 0.5;

    // Lambdas for all the different curves we want to plot
    auto sigT_R = [&](double W)
    {
        return CB->sigmaT_R(W, Q2);
    };
    auto sigT_NR = [&](double W)
    {
        return CB->sigmaT_NR(W, Q2);
    };
    auto sigT = [&](double W)
    {
        return CB->sigma_T(W, Q2);
    };

    auto sigL_R = [&](double W)
    {
        return CB->sigmaL_R(W, Q2);
    };
    auto sigL_NR = [&](double W)
    {
        return CB->sigmaL_NR(W, Q2);
    };
    auto sigL = [&](double W)
    {
        return CB->sigma_L(W, Q2);
    };

    plotter plotter;

    // sigma_T plot
    plot p1 = plotter.new_plot();
    p1.set_curve_points(500);
    p1.set_ranges({1, 3}, {0, 400});
    p1.set_labels("#it{W} [GeV]", "#sigma_{#it{T}}(#it{W}, #it{Q}^{2})  [#mub]");
    p1.add_header("#it{Q}^{2} = 0.5 GeV^{2}");
    p1.set_legend(0.65,0.5);

    p1.add_curve({M_PROTON+M_PION, 3}, sigT_NR, "Non-resonant");
    p1.add_curve({M_PROTON+M_PION, 3}, sigT_R,  "Resonances");
    p1.add_curve({M_PROTON+M_PION, 3}, sigT,    "Total");

    // sigma_L
    plot p2 = plotter.new_plot();
    p2.set_curve_points(500);
    p2.set_ranges({1, 3}, {0, 35});
    p2.set_labels("#it{W} [GeV]", "#sigma_{#it{L}}(#it{W}, #it{Q}^{2})  [#mub]");
    p2.add_header("#it{Q}^{2} = 0.5 GeV^{2}");
    p2.set_legend(0.65,0.6);

    p2.add_curve({M_PROTON+M_PION, 3}, sigL_NR, "Non-resonant");
    p2.add_curve({M_PROTON+M_PION, 3}, sigL_R,  "Resonances");
    p2.add_curve({M_PROTON+M_PION, 3}, sigL,    "Total");

    plotter.combine({2,1}, {p1,p2}, "BC_sigmas.pdf");
};