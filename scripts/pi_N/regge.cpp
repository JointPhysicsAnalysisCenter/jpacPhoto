// Comparison of the Reggeized pion exchange with different values of interaction 
// radius, effective pole, and effective zero
// Reproduces Figure 2 of [1]
//
// OUTPUT: fig2.pdf
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Universitat Bonn (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] - arXiv:2407.19577 [hep-ph]
// ------------------------------------------------------------------------------

#include "constants.hpp"
#include "kinematics.hpp"
#include "plotter.hpp"

#include "pion/mgi_pion.hpp"
#include "pion/electric_nucleon.hpp"
#include "pion/magnetic_nucleon.hpp"
#include "piN/slac/data.hpp"

void regge()
{
    using namespace jpacPhoto;
    using namespace jpacPhoto::piN;

    // ---------------------------------------------------------------------------
    // Amplitude set up
    // ---------------------------------------------------------------------------


    kinematics kpi = new_kinematics(M_PION, M_PROTON);
    kpi->set_meson_JP(PSEUDOSCALAR);
    kpi->set_baryon_JP(HALFPLUS);

    amplitude pi_J0      = new_amplitude<mgi_pion>(kpi, mgi_pion::kFixedSpin, "J =0 ");
    amplitude pi_R_asymp = new_amplitude<mgi_pion>(kpi, mgi_pion::kHighEnergyLimit, "#it{s} #rightarrow #infty");
    pi_R_asymp->set_parameters({1/2.});

    amplitude pi_R_high = new_amplitude<mgi_pion>(kpi, mgi_pion::kResummed, "large jz");
    amplitude pi_R_low  = new_amplitude<mgi_pion>(kpi, mgi_pion::kResummed, "large jp");

    // ---------------------------------------------------------------------------
    // Plot results
    // ---------------------------------------------------------------------------

    double Egam = 16;
    double s = s_cm(Egam),  tmin = - kpi->t_min(s) + 0.003;

    plotter plotter;

    auto add_curves = [&](class plot & p, double R2, std::string label, jpacColor color)
    {

        pi_R_asymp->set_parameters({R2});
        p.color_offset(1);
        p.add_curve({tmin, 0.2}, [s,pi_R_asymp](double mt){ return pi_R_asymp->differential_xsection(s, -mt) * 1E-3; }, dashed(color, label)); 

        pi_R_high->set_parameters({R2, 1000,  0.5   });
        pi_R_low->set_parameters( {R2, 1,     1000.5});
        std::vector<double> x, yh, yl;
        for (int i = 0; i < 20; i++)
        {
            double tx = tmin + double(i)*(0.2 - tmin)/19.;
            x.push_back(tx);
            yh.push_back(pi_R_high->differential_xsection(s, -tx) * 1E-3);
            yl.push_back(pi_R_low->differential_xsection( s, -tx) * 1E-3);
        };
        p.add_band(x, {yh, yl});
    };

    auto plot_dsigma = [&](class plotter & pltr, std::array<double,2> yrange, bool legend = false)
    {          
        plot p = pltr.new_plot();
        p.set_legend(0.32, 0.20+0.07*!legend, 1.1);
        p.set_legend(legend);
        p.set_ranges({0, 0.2}, yrange);
        p.set_curve_points(20);
        p.set_logscale(false, true);
        p.set_labels("#minus #it{t} [GeV^{2}]", "d#sigma/d#it{t}  [#mub GeV^{-2}]");
        p.add_header(var_def("#it{E}_{#gamma}", Egam, "GeV"));

        double fm2gev = 5.068;

        p.add_curve({tmin, 0.2}, [&](double mt){ return pi_J0->differential_xsection(s, -mt) * 1E-3; }, "#it{J} = 0"); 

        add_curves(p, 1.,     "#it{R} = 1/#sqrt{#it{s}_{0}}",    jpacColor::Red);
        add_curves(p, pow(1*fm2gev, 2), "#it{R} = 1 fm",         jpacColor::Green);
        add_curves(p, pow(2*fm2gev, 2), "#it{R} = 2 fm",         jpacColor::Orange);

        return p;
    };

    plot p1 = plot_dsigma(plotter, {1E-2, 1}, true);

    p1.save("fig2.pdf");
};