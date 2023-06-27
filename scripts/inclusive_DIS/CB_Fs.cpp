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

#include "sigma_tot/CB_F.hpp"
#include "plotter.hpp"

void CB_Fs()
{
    using namespace jpacPhoto;

    auto xF1 = new_total_xsection<CB_F>(1); 
    auto xF2 = new_total_xsection<CB_F>(2); 

    // Pick some fixed Q2
    double Q2 = 0.5;

    // Lambdas for all the different curves we want to plot
    auto F1 = [&](double W)
    {
        return xF1->evaluate(W*W, -Q2);
    };

    auto F2 = [&](double W)
    {
        return xF2->evaluate(W*W, -Q2);
    };

    plotter plotter;

    // F1 plot
    plot p1 = plotter.new_plot();
    p1.set_curve_points(1000);
    p1.set_ranges({1, 3}, {0, 4});
    p1.set_labels("#it{W} [GeV]", "#it{F}_{1}(#it{W}, #it{Q}^{2})");
    p1.set_legend(0.25,0.7);

    Q2 = 0.1;
    p1.add_curve({M_PROTON+M_PION, 3}, F1, "#it{Q}^{2} = 0.1 GeV^{2}");
    
    Q2 = 0.5;
    p1.add_curve({M_PROTON+M_PION, 3}, F1, "#it{Q}^{2} = 0.5 GeV^{2}");

    Q2 = 1.0;
    p1.add_curve({M_PROTON+M_PION, 3}, F1, "#it{Q}^{2} = 1.0 GeV^{2}");

    // plot p2 = plotter.new_plot();
    // p2.set_curve_points(1000);
    // p2.set_ranges({1, 3}, {0, 0.5});
    // p2.set_labels("#it{W} [GeV]", "#it{F}_{2}(#it{W}, #it{Q}^{2})");
    // p2.set_legend(0.65,0.35);

    // Q2 = 0.1;
    // p2.add_curve({M_PROTON+M_PION, 3}, F2, "#it{Q}^{2} = 0.1 GeV^{2}");
    
    // Q2 = 0.5;
    // p2.add_curve({M_PROTON+M_PION, 3}, F2, "#it{Q}^{2} = 0.5 GeV^{2}");

    // Q2 = 1.0;
    // p2.add_curve({M_PROTON+M_PION, 3}, F2, "#it{Q}^{2} = 1.0 GeV^{2}");


    plotter.combine({2,1}, {p1,p1}, "BC_Fs.pdf");
};