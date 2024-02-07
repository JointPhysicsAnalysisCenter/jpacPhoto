// Interface functions for the Boyarski and GlueX data of pi delta photoproduction
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Universitat Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------

#ifndef PIDELTA_PLOTS_HPP
#define PIDELTA_PLOTS_HPP

#include "data_set.hpp"
#include "plotter.hpp"
#include "piDelta/data.hpp"
#include "constants.hpp"

namespace jpacPhoto 
{ 
    namespace piDelta 
    {
        inline plot plot_differential(plotter& plotr)
        {
            plot p = plotr.new_plot();

            p.add_data(piDelta::differential());
            p.set_labels("#minus#it{t} [GeV]", "d#sigma/d#it{t}  [#mub]");
            p.set_logscale(true, true);
            p.add_header("#it{E}_{#gamma} = 8 GeV");
            p.set_legend(0.3, 0.4);
            p.set_ranges({1E-3, 2}, {5E-2, 7});
            return p;
        };

        inline plot plot_SDME(plotter& plotr, int a, int m, int mp)
        {
            plot p = plotr.new_plot();
            p.add_data(piDelta::SDME(a, m, mp));
            p.add_header("#it{E}_{#gamma} = 8.5 GeV");
            p.set_legend(0.7, 0.2);
            p.set_ranges({0, 1.2}, {-0.5, 0.5});
            p.set_labels("#minus #it{t}  [GeV^{2}]", 
                         "#rho^{" + std::to_string(a) + "}_{ " + std::to_string(m) + std::to_string(mp) + "}");
            return p;
        };

        inline std::vector<plot> plot_SDMEs(plotter & plotr)
        {
            return { plot_SDME(plotr, 0, 1, 1), plot_SDME(plotr, 0, 3,  1), plot_SDME(plotr, 0, 3, -1),
                     plot_SDME(plotr, 1, 1, 1), plot_SDME(plotr, 1, 3,  3), plot_SDME(plotr, 1, 3,  1), plot_SDME(plotr, 1, 3, -1),
                     plot_SDME(plotr, 2, 3, 1), plot_SDME(plotr, 2, 3, -1) };
        };

        inline plot plot_beam_asymmetry(plotter& plotr)
        {
            plot p = plotr.new_plot();
            p.add_data(piDelta::beam_asymmetry());
            p.add_header("#it{E}_{#gamma} = 8.5 GeV");
            p.set_legend(0.7, 0.2);
            p.set_ranges({0, 1.2}, {-1, 1});
            p.set_labels("#minus #it{t}  [GeV^{2}]", 
                            "#Sigma_{4#pi}");
            return p;
        };

    }; 
};

#endif