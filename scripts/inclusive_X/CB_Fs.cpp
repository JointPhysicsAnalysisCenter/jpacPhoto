// Plots F_1 and F_2 at fixed Q2 as a function of W using the low-energy
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

#include "semi_inclusive/CB_F.hpp"
#include "plotter.hpp"

void CB_Fs()
{
    using namespace jpacPhoto;

    int p = CB_F::kProton, n = CB_F::kNeutron;

    auto pF1 = new_inclusive_function<CB_F>(1, p), pF2 = new_inclusive_function<CB_F>(2, p); 
    auto nF1 = new_inclusive_function<CB_F>(1, n), nF2 = new_inclusive_function<CB_F>(2, n); 


    std::array<double,2> range = {M_PROTON+M_PION, 3};
    // F1 plot
    plotter plotter;
    plot p1 = plotter.new_plot();
    p1.set_curve_points(1000);
    p1.set_ranges({1, 3}, {0, 4});
    p1.set_labels("#it{M}_{#it{X}} [GeV]", "#it{F}_{1}(#it{x}_{B}, #it{t})");
    p1.set_legend(0.25,0.75);

    p1.add_curve( range, [&](double w){ return pF1->evaluate(w*w, -0.1);}, "#it{t} = #minus 0.1 GeV^{2}");
    p1.add_dashed(range, [&](double w){ return nF1->evaluate(w*w, -0.1);});
    p1.add_curve( range, [&](double w){ return pF1->evaluate(w*w, -2.0);}, "#it{t} = #minus 2 GeV^{2}");
    p1.add_dashed(range, [&](double w){ return nF1->evaluate(w*w, -2.0);});
    p1.add_curve( range, [&](double w){ return pF1->evaluate(w*w, -10);}, "#it{t} = #minus 10 GeV^{2}");
    p1.add_dashed(range, [&](double w){ return nF1->evaluate(w*w, -10);});

    plot p2 = plotter.new_plot();
    p2.set_curve_points(1000);
    p2.set_ranges({1, 3}, {0, 0.5});
    p2.set_labels("#it{M}_{#it{X}} [GeV]", "#it{F}_{2}(#it{x}_{B}, #it{t})");
    p2.set_legend(0.25,0.75);

    p2.add_curve( range, [&](double w){ return pF2->evaluate(w*w, -0.1);}, "#it{t} = #minus 0.1 GeV^{2}");
    p2.add_dashed(range, [&](double w){ return nF2->evaluate(w*w, -0.1);});
    p2.add_curve( range, [&](double w){ return pF2->evaluate(w*w, -2.0);}, "#it{t} = #minus 2 GeV^{2}");
    p2.add_dashed(range, [&](double w){ return nF2->evaluate(w*w, -2.0);});
    p2.add_curve( range, [&](double w){ return pF2->evaluate(w*w, -10);}, "#it{t} = #minus 10 GeV^{2}");
    p2.add_dashed(range, [&](double w){ return nF2->evaluate(w*w, -10);});

    plotter.combine({2,1}, {p1,p2}, "CB_Fs.pdf");
};