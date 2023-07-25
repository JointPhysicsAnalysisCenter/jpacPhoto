// Comparison of b1 inclusive photoproduction predictions with differential data
// from the OmegaPHOTON collaboration 
// Reproduces fig. 5 in [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] 	arXiv:2209.05882 [hep-ph]
// ------------------------------------------------------------------------------

#include "inclusive_process.hpp"
#include "inclusive/pion_exchange.hpp"
#include "plotter.hpp"

void OmegaPhoton()
{
    using namespace jpacPhoto;
    using complex = std::complex<double>;
    
    plotter plotter;

    // ---------------------------------------------------------------------------
    // Data info
    // ---------------------------------------------------------------------------
    
    double s = 75.9421;

    std::vector<double> x, sig, dsig, dx;

    x    = {0.65, 0.75, 0.85, 0.95};
    dx   = {0.05, 0.05, 0.05, 0.05};
    sig  = {1.80957, 2.15690, 1.3661, 0.65901};
    dsig = {2.36188 - 1.80957, 2.47490 - 2.15690, 1.53345 - 1.36611, 0.76779 - 0.65901};

    // ---------------------------------------------------------------------------
    // Inclusive amplitude
    // ---------------------------------------------------------------------------

    inclusive_process b1p = new_inclusive_process<pion_exchange>(M_B1, +1, "b_{1}(1235)^{#plus}");
    b1p->reggeized(true);
    b1p->set_parameters({0.24});
    
    // ---------------------------------------------------------------------------
    // Make plot
    // ---------------------------------------------------------------------------

    int N = 1000;
    std::array<double,2> bounds = {0.7, 1};

    // b1 minus plot
    auto dsigdx = [&](double x)
    {
        return b1p->dsigma_dx(s, x) * 1E-3; // in mub
    };
    
    plot p1 = plotter.new_plot();
    p1.set_curve_points(N);
    p1.set_legend(0.3, 0.3);
    p1.set_ranges({0.7,1}, {0, 2.5});
    p1.set_labels("#it{W}_{#gamma#it{p}}  [GeV]", "d#sigma / d#it{x} [#mub]");

    // Add data
    p1.add_data(x, sig, {dx, dsig}, "Omega Photon");

    // Plot both the cross section with resonances 
    b1p->set_option(pion_exchange::kJPAC);
    p1.add_curve( bounds, dsigdx, "Inclusive #it{b}_{1}(1235)^{#plus}");
    // and without
    b1p->set_option(pion_exchange::kPDG);
    p1.add_dashed(bounds, dsigdx);


    p1.save("b1_OmegaPhoton.pdf");
};