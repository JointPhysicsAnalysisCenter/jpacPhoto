#include "constants.hpp"
#include "kinematics.hpp"
#include "plotter.hpp"

#include "pion/mgi_pion.hpp"
#include "pion/electric_nucleon.hpp"
#include "pion/magnetic_nucleon.hpp"
#include "piN/slac/data.hpp"

void regge_compare()
{
    using namespace jpacPhoto;
    using namespace jpacPhoto::piN;

    // ---------------------------------------------------------------------------
    // Amplitude set up
    // ---------------------------------------------------------------------------


    kinematics kpi = new_kinematics(M_PION, M_PROTON);
    kpi->set_meson_JP(PSEUDOSCALAR);
    kpi->set_baryon_JP(HALFPLUS);

    amplitude pi_R_asymp = new_amplitude<mgi_pion>(kpi, mgi_pion::kHighEnergyLimit, "#it{s} #rightarrow #infty");
    pi_R_asymp->set_parameters({1/2.});

    amplitude pi_R1 = new_amplitude<mgi_pion>(kpi, mgi_pion::kResummed, "#it{j}#kern[-0.15]{_{#it{p}}} = 1, #it{j}#kern[-0.2]{_{#it{z}}} = #frac{1}{2}");
    pi_R1->set_parameters({1/2., 1, 1./2});

    amplitude pi_R2 = new_amplitude<mgi_pion>(kpi, mgi_pion::kResummed, "#it{j}#kern[-0.15]{_{#it{p}}} = 1, #it{j}#kern[-0.2]{_{#it{z}}} = 1");
    pi_R2->set_parameters({1/2., 1, 1.});

    amplitude pi_R3 = new_amplitude<mgi_pion>(kpi, mgi_pion::kResummed, "#it{j}#kern[-0.15]{_{#it{p}}} = 2, #it{j}#kern[-0.2]{_{#it{z}}} = #frac{1}{2}");
    pi_R3->set_parameters({1/2., 2., 1./2.});


    // ---------------------------------------------------------------------------
    // Plot results
    // ---------------------------------------------------------------------------

    plotter plotter;

    std::vector<amplitude> to_plot = {pi_R_asymp, pi_R1, pi_R2, pi_R3};

    auto add_curves = [&](class plot & p, double Egam, double R2, std::string label, jpacColor color)
    {
        double s = s_cm(Egam),  tmin = - kpi->t_min(s) + 0.0001;
        pi_R_asymp->set_parameters({R2});
        pi_R1->set_parameters({R2, 1,  1./2});
        pi_R2->set_parameters({R2, 1,  1.});
        pi_R3->set_parameters({R2, 2., 1./2.});
        p.add_curve({tmin, 0.2}, [s,pi_R_asymp](double mt){ return pi_R_asymp->differential_xsection(s, -mt) * 1E-3; },  solid(color , label)); 
        p.add_curve({tmin, 0.2}, [s,pi_R1]     (double mt){ return pi_R1->differential_xsection(s, -mt)      * 1E-3; },  dashed(color)); 
        p.add_curve({tmin, 0.2}, [s,pi_R2]     (double mt){ return pi_R2->differential_xsection(s, -mt)      * 1E-3; },  dotted(color)); 
        p.add_curve({tmin, 0.2}, [s,pi_R3]     (double mt){ return pi_R3->differential_xsection(s, -mt)      * 1E-3; },  dash_dotted(color)); 
    };

    auto plot_dsigma = [&](class plotter & pltr, double Egam, std::array<double,2> yrange, bool legend = false)
    {  
        plot p = pltr.new_plot();
        p.set_legend(0.22, 0.65+0.07*!legend, 1.1);
        p.set_legend(legend);
        p.set_ranges({0, 0.2}, yrange);
        p.set_curve_points(40);
        p.set_logscale(false, true);
        p.set_labels("#minus #it{t} [GeV^{2}]", "d#sigma/d#it{t}  [#mub GeV^{-2}]");
        p.add_header(var_def("#it{E}_{#gamma}", Egam, "GeV"));
        if (legend)
        { 
            p.add_style_legend({  "#it{s} #rightarrow #infty", 
                              "#it{j}#kern[-0.15]{_{#it{p}}} = 1, #it{j}#kern[-0.2]{_{#it{z}}} = #frac{1}{2}", 
                              "#it{j}#kern[-0.15]{_{#it{p}}} = 1, #it{j}#kern[-0.2]{_{#it{z}}} = 1", 
                              "#it{j}#kern[-0.15]{_{#it{p}}} = 2, #it{j}#kern[-0.2]{_{#it{z}}} = #frac{1}{2}"} );
            p.set_style_legend(0.7, 0.65, 1.);
        };

        
        double fm2gev = 5.068;

        add_curves(p, Egam, 1/2.,     "#it{R}^{2} = (2#it{s}_{0})^{-1}", jpacColor::Blue);
        add_curves(p, Egam, 1*fm2gev, "#it{R}^{2} = 1 fm",               jpacColor::Red);
        add_curves(p, Egam, 2*fm2gev, "#it{R}^{2} = 2 fm",               jpacColor::Green);

        return p;
    };

    plot p1 = plot_dsigma(plotter, 16, {3E-2, 1}, true);
    plot p2 = plot_dsigma(plotter, 11, {1E-1, 1.2});
    plot p3 = plot_dsigma(plotter, 8 , {1E-1, 3});
    plot p4 = plot_dsigma(plotter, 5 , {4E-1, 10});

    plotter.combine({2,2}, {p1,p2,p3,p4}, "fig3.pdf");

};