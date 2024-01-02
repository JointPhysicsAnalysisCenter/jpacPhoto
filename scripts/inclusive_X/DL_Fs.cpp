// Plots inclusive unpolarized structure functions F_1 and F_2 at a variety of
// Q2 values as a function of Bjorken x.
// The high-energy parameterization of [1] is used.
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] - https://arxiv.org/abs/hep-ph/0402081
// ------------------------------------------------------------------------------

#include "sigma_tots/DL_F.hpp"
#include "plotter.hpp"

void DL_Fs()
{
    using namespace jpacPhoto;
    using complex = std::complex<double>;

    inclusive_function F1 = new_inclusive_function<DL_F>(1);
    inclusive_function F2 = new_inclusive_function<DL_F>(2);


    std::array<double,2> range = {M_PROTON+M_PION+EPS, 30};
    // F1 plot
    plotter plotter;
    plot p1 = plotter.new_plot();
    p1.set_curve_points(100);
    p1.set_logscale(false, true);
    p1.set_ranges({1, 30}, {0.2, 1000});
    p1.set_labels("#it{M}_{#it{X}} [GeV]", "#it{F}_{1}(#it{M}_{#it{X}}^{2}, #it{t})");
    p1.set_legend(0.25,0.75);

    p1.add_curve( range, [&](double w){ return F1->evaluate(w*w, -0.1);}, "#it{t} = #minus 0.1 GeV^{2}");
    p1.add_curve( range, [&](double w){ return F1->evaluate(w*w, -1.0);}, "#it{t} = #minus 1.0 GeV^{2}");
    p1.add_curve( range, [&](double w){ return F1->evaluate(w*w, -10);}, "#it{t} = #minus 10 GeV^{2}");

    plot p2 = plotter.new_plot();
    p2.set_curve_points(100);
    p2.set_logscale(false, true);
    p2.set_ranges({1, 30}, {1E-2, 10});
    p2.set_labels("#it{M}_{#it{X}} [GeV]", "#it{F}_{2}(#it{M}_{#it{X}}^{2}, #it{t})");
    p2.set_legend(0.25,0.75);

    p2.add_curve( range, [&](double w){ return F2->evaluate(w*w, -0.1);}, "#it{t} = #minus 0.1 GeV^{2}");
    p2.add_curve( range, [&](double w){ return F2->evaluate(w*w, -1.0);}, "#it{t} = #minus 1.0 GeV^{2}");
    p2.add_curve( range, [&](double w){ return F2->evaluate(w*w, -10);}, "#it{t} = #minus 10 GeV^{2}");

    plotter.combine({2,1}, {p1,p2}, "DL_Fs.pdf");
};