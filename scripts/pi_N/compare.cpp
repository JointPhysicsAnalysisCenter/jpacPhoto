#include "constants.hpp"
#include "kinematics.hpp"
#include "plotter.hpp"

#include "pion/mgi_pion.hpp"
#include "pion/electric_nucleon.hpp"
#include "pion/magnetic_nucleon.hpp"
#include "pion/VGL.hpp"
#include "piN/slac/data.hpp"

void compare()
{
    using namespace jpacPhoto;
    using namespace jpacPhoto::piN;

    // ---------------------------------------------------------------------------
    // Amplitude set up
    // ---------------------------------------------------------------------------


    kinematics kpi = new_kinematics(M_PION, M_PROTON);
    kpi->set_meson_JP(PSEUDOSCALAR);
    kpi->set_baryon_JP(HALFPLUS);

    amplitude pi  = new_amplitude<mgi_pion>(kpi, mgi_pion::kFixedSpin, "Bare");

    amplitude pi_R = new_amplitude<mgi_pion>(kpi, mgi_pion::kHighEnergyLimit, "Regge");
    pi_R->set_parameters({1/2.});

    amplitude N_e = new_amplitude<electric_nucleon>(kpi);

    // Magnetic with and without FF
    amplitude N_m   = new_amplitude<magnetic_nucleon>(kpi);
    amplitude N_mFF = new_amplitude<magnetic_nucleon>(kpi);
    N_mFF->set_option(magnetic_nucleon::kAddFF);

    amplitude vgl   = new_amplitude<VGL>(kpi, "VGL");

    amplitude sum1 = pi + N_e + N_m;
    sum1->set_id("Elementary");

    amplitude sum2 = pi_R + N_e + N_m;
    sum2->set_id("Reggeized");

    amplitude sum3 = pi_R + N_e + N_mFF;
    sum3->set_id("Reggeized with FF");

    // ---------------------------------------------------------------------------
    // Data from SLAC (1968)
    // ---------------------------------------------------------------------------

    std::vector<data_set> boyarski = piN::differential();

    // ---------------------------------------------------------------------------
    // Plot results
    // ---------------------------------------------------------------------------

    plotter plotter;

    auto plot_dsigma = [&](class plotter & pltr, double Egam, std::array<double,2> yrange, bool legend = false)
    {  
        plot p = pltr.new_plot();
        p.print_to_terminal(true);
        p.set_legend(0.22, 0.65+0.07*!legend, 1.1);
        p.set_legend(legend);
        p.set_ranges({0,0.2}, yrange);
        p.set_labels("#minus #it{t} [GeV^{2}]", "d#sigma/d#it{t}  [#mub GeV^{-2}]");
        p.add_header(var_def("#it{E}_{#gamma}", Egam, "GeV"));
        if (legend)
        {
            p.add_style_legend({"#it{R}^{2} = (2#it{s}_{0})^{-1}", "#it{R}^{2} = 1 fm", "", ""});
            p.set_style_legend(0.70, 0.65, 0.5);
        };

        double s    = s_cm(Egam);
        double tmin = -kpi->t_min(s) + 0.001 ;

        p.add_curve( {tmin, 0.2}, [s, vgl] (double mt){ return vgl->differential_xsection(s, -mt) * 1E-3; },   vgl->id()); 
        p.add_curve( {tmin, 0.2}, [s, sum1](double mt){ return sum1->differential_xsection(s, -mt) * 1E-3; },  sum1->id()); 

        for (auto amp : {sum2, sum3}) 
        { 
            pi_R->set_parameters({1./2});
            p.add_curve( {tmin, 0.2}, [s, amp](double mt){ return amp->differential_xsection(s, -mt) * 1E-3; }, amp->id()); 
            pi_R->set_parameters({5.068});
            p.add_dashed({tmin, 0.2}, [s, amp](double mt){ return amp->differential_xsection(s, -mt) * 1E-3; }); 
        };
        return p;
    };

    plot p1 = plot_dsigma(plotter, 16, {0., 0.4}, true);
    plot p2 = plot_dsigma(plotter, 11, {0., 0.8});
    plot p3 = plot_dsigma(plotter, 8 , {0., 1.7});
    plot p4 = plot_dsigma(plotter, 5 , {0., 4.5});

    p1.add_data(boyarski[3]);
    p2.add_data(boyarski[2]);
    p3.add_data(boyarski[1]);
    p4.add_data(boyarski[0]);

    plotter.combine({2,2}, {p1,p2,p3,p4}, "fig5.pdf");

};