// Plots F_1 and F_2 at fixed Q2 as a function of W 
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
// [2] - https://arxiv.org/abs/hep-ph/0402081
// ------------------------------------------------------------------------------

#include "sigma_tot/ChristyBosted_F.hpp"
#include "sigma_tot/DonnachieLandshoff_F.hpp"
#include "sigma_tot/combined_F.hpp"

#include "plotter.hpp"

void combined_Fs()
{
    using namespace jpacPhoto;

    int p = ChristyBosted_F::kProton;
    double Wth = M_PROTON+M_PION;

    auto lF1 = new_inclusive_function<ChristyBosted_F>(1, p); 
    auto lF2 = new_inclusive_function<ChristyBosted_F>(2, p); 

    auto hF1 = new_inclusive_function<DonnachieLandshoff_F>(1); 
    auto hF2 = new_inclusive_function<DonnachieLandshoff_F>(2); 

    auto cF1 = new_inclusive_function<combined_F>(1, p); 
    auto cF2 = new_inclusive_function<combined_F>(2, p); 

    double Q2;     
    inclusive_function to_plot;

    // Lambdas for all the different curves we want to plot

    auto F = [&](double W)
    {
        return to_plot->evaluate(W*W, -Q2);
    };

    plotter plotter;

    // F1 plot
    plot p1 = plotter.new_plot();
    p1.set_curve_points(1000);
    p1.set_ranges({1, 10}, {0, 30});
    p1.set_labels("#it{W} [GeV]", "#it{F}_{1}(#it{W}, #it{Q}^{2})");
    p1.set_legend(0.25,0.7);

    Q2 = 0.1; 
    to_plot = cF1;
    p1.add_curve(  {Wth, 10}, F, "#it{Q}^{2} = 0.1 GeV^{2}");
    to_plot = lF1;
    p1.add_dashed( {Wth, 7}, F);
    to_plot = hF1;
    p1.add_dashed( {3, 10}, F);

    Q2 = 0.5; 
    to_plot = cF1;
    p1.add_curve(  {Wth, 10}, F, "#it{Q}^{2} = 0.5 GeV^{2}");
    to_plot = lF1;
    p1.add_dashed( {Wth, 7}, F);
    to_plot = hF1;
    p1.add_dashed( {3, 10}, F);

    Q2 = 1.0; 
    to_plot = cF1;
    p1.add_curve(  {Wth, 10}, F, "#it{Q}^{2} = 1.0 GeV^{2}");
    to_plot = lF1;
    p1.add_dashed( {Wth, 7}, F);
    to_plot = hF1;
    p1.add_dashed( {3, 10}, F);

    Q2 = 5.0; 
    to_plot = cF1;
    p1.add_curve(  {Wth, 10}, F, "#it{Q}^{2} = 5.0 GeV^{2}");
    to_plot = lF1;
    p1.add_dashed( {Wth, 7}, F);
    to_plot = hF1;
    p1.add_dashed( {3, 10}, F);

    // F2 plot
    plot p2 = plotter.new_plot();
    p2.set_curve_points(1000);
    p2.set_ranges({1, 10}, {0, 1});
    p2.set_labels("#it{W} [GeV]", "#it{F}_{2}(#it{W}, #it{Q}^{2})");
    p2.set_legend(0.5,0.7);

    Q2 = 0.1; 
    to_plot = cF2;
    p2.add_curve(  {Wth, 10}, F, "#it{Q}^{2} = 0.1 GeV^{2}");
    to_plot = lF2;
    p2.add_dashed( {Wth, 7}, F);
    to_plot = hF2;
    p2.add_dashed( {2, 10}, F);

    Q2 = 0.5; 
    to_plot = cF2;
    p2.add_curve(  {Wth, 10}, F, "#it{Q}^{2} = 0.5 GeV^{2}");
    to_plot = lF2;
    p2.add_dashed( {Wth, 7}, F);
    to_plot = hF2;
    p2.add_dashed( {2, 10}, F);

    Q2 = 1.0; 
    to_plot = cF2;
    p2.add_curve(  {Wth, 10}, F, "#it{Q}^{2} = 1.0 GeV^{2}");
    to_plot = lF2;
    p2.add_dashed( {Wth, 7}, F);
    to_plot = hF2;
    p2.add_dashed( {2, 10}, F);

    Q2 = 5.0; 
    to_plot = cF2;
    p2.add_curve(  {Wth, 10}, F, "#it{Q}^{2} = 5.0 GeV^{2}");
    to_plot = lF2;
    p2.add_dashed( {Wth, 7}, F);
    to_plot = hF2;
    p2.add_dashed( {2, 10}, F);

    plotter.combine({2,1}, {p1, p2}, "Fs.pdf");
};