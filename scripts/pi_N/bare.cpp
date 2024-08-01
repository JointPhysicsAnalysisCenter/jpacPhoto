// Comparison of different components to the Born model of charged pi photoproduction
// with an unReggeized (i.e. bare) pion exchange
//
// OUTPUT: bare_pion.pdf
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Universitat Bonn (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------


#include "constants.hpp"
#include "kinematics.hpp"
#include "plotter.hpp"

#include "pion/mgi_pion.hpp"
#include "pion/electric_nucleon.hpp"
#include "pion/magnetic_nucleon.hpp"
#include "piN/slac/data.hpp"

void bare()
{
    using namespace jpacPhoto;
    using namespace jpacPhoto::piN;

    // ---------------------------------------------------------------------------
    // Amplitude set up
    // ---------------------------------------------------------------------------

    int pi_plus  = mgi_pion::kPiPlus;
    int pi_minus = mgi_pion::kPiMinus;

    kinematics kpi = new_kinematics(M_PION, M_PROTON);
    kpi->set_meson_JP(PSEUDOSCALAR);
    kpi->set_baryon_JP(HALFPLUS);

    amplitude pi  = new_amplitude<mgi_pion>(kpi, mgi_pion::kFixedSpin, "m.g.i. #pi");
    
    amplitude N_e = new_amplitude<electric_nucleon>(kpi);
    amplitude N_m = new_amplitude<magnetic_nucleon>(kpi);

    amplitude electric = pi + N_e;
    amplitude magnetic = N_m;

    amplitude total = electric + magnetic;
    total->set_id("Sum");

    // ---------------------------------------------------------------------------
    // Data from SLAC (1968)
    // ---------------------------------------------------------------------------

    std::vector<data_set> boyarski = piN::differential();

    // ---------------------------------------------------------------------------
    // Plot results
    // ---------------------------------------------------------------------------

    plotter plotter;

    amplitude to_plot = total;

    auto plot_dsigma = [&](class plotter & pltr, double Egam, std::array<double,2> yrange, bool legend = false)
    {  
        plot p = pltr.new_plot();
        p.print_to_terminal(true);
        p.set_legend(0.22, 0.65, 1.);
        p.set_legend(legend);
        p.set_ranges({0,0.2}, yrange);
        p.set_labels("#minus #it{t} [GeV^{2}]", "d#sigma/d#it{t}  [#mub GeV^{-2}]");
        p.add_header(var_def("#it{E}_{#gamma}", Egam, "GeV"));

        double s    = s_cm(Egam);
        double tmin = -to_plot->get_kinematics()->t_min(s) + 0.001 ;

        for (auto amp : expand(to_plot)) 
        { 
            amp->set_option(pi_plus);
            p.add_curve({tmin, 0.2}, [s, amp](double mt){ return amp->differential_xsection(s, -mt) * 1E-3; }, amp->id()); 
            amp->set_option(pi_minus);
            p.add_dashed({tmin, 0.2}, [s, amp](double mt){ return amp->differential_xsection(s, -mt) * 1E-3; }); 
        };
        return p;
    };

    plot p1 = plot_dsigma(plotter, 16, {0., 0.6}, true);
    plot p2 = plot_dsigma(plotter, 11, {0., 1.2});
    plot p3 = plot_dsigma(plotter, 8 , {0., 2.2});
    plot p4 = plot_dsigma(plotter, 5 , {0., 5.5});

    p1.add_data(boyarski[3]);
    p2.add_data(boyarski[2]);
    p3.add_data(boyarski[1]);
    p4.add_data(boyarski[0]);

    plotter.combine({2,2}, {p1,p2,p3,p4}, "bare_pion.pdf");
};