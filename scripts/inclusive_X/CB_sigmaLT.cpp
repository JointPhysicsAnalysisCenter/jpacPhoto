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

#include "sigma_tots/CB_F.hpp"
#include "plotter.hpp"

void CB_sigmaLT()
{
    using namespace jpacPhoto;

    int p = CB_F::kProton, n = CB_F::kNeutron, N = CB_F::kAvgNucleon;
    auto CB = dynamic_pointer_cast<CB_F>(new_inclusive_function<CB_F>(1));

    // Pick some fixed -t 
    std:array<double,2> range = {M_PROTON+M_PION + 10.*EPS, 3.};

    // sigma_T plot
    plotter plotter;
    plot pT = plotter.new_plot();
    pT.set_curve_points(500);
    pT.set_ranges({1, 3}, {0, 400});
    pT.add_header("#it{t} = #minus 0.5 GeV^{2}");
    pT.set_labels("#it{M}_{#it{X}} [GeV]", "#sigma_{#it{T}}(#it{M}_{#it{X}}^{2}, #it{t})  [#mub]");
    pT.set_legend(0.65,0.5);

    pT.add_curve(  range, [&](double M){ return CB->sigma_T(  p, M, 0.5); }, "Total");
    pT.add_dashed( range, [&](double M){ return CB->sigma_T(  n, M, 0.5); });
    pT.add_curve(  range, [&](double M){ return CB->sigmaT_R( p, M, 0.5); }, "Resonant");
    pT.add_dashed( range, [&](double M){ return 2.*CB->sigmaT_R(N, M, 0.5) - CB->sigmaT_R(p, M, 0.5); });
    pT.add_curve(  range, [&](double M){ return CB->sigmaT_NR(p, M, 0.5); }, "Non-resonant");
    pT.add_dashed( range, [&](double M){ return 2.*CB->sigmaT_NR(N, M, 0.5) - CB->sigmaT_NR(p, M, 0.5); });

    plot pL = plotter.new_plot();
    pL.set_curve_points(500);
    pL.set_ranges({1, 3}, {0, 40});
    pL.add_header("#it{t} = #minus 0.5 GeV^{2}");
    pL.set_labels("#it{M}_{#it{X}} [GeV]", "#sigma_{#it{L}}(#it{M}_{#it{X}}^{2}, #it{t})  [#mub]");
    pL.set_legend(0.65,0.6);

    pL.add_curve(  range, [&](double M){ return CB->sigma_L(  p, M, 0.5); }, "Total");
    pL.add_dashed( range, [&](double M){ return CB->sigma_L(  n, M, 0.5); });
    pL.add_curve(  range, [&](double M){ return CB->sigmaL_R( p, M, 0.5); }, "Resonant");
    pL.add_curve(  range, [&](double M){ return CB->sigmaL_NR(p, M, 0.5); }, "Non-resonant");

    plotter.combine({2,1}, {pT, pL}, "sigmas.pdf");

};