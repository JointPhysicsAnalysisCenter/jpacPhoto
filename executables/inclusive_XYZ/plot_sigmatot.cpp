#include "inclusive/total_xsection.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>
#include <memory>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    int N = 200;

    double  ymin = 1.;
    double  ymax = 1.E3;

    double xmin = (M_PION + M_PROTON) + EPS;   
    double xmax = 5.;

    std::string filename = "sigma.pdf";
    std::string xlabel   = "W  [GeV]";
    std::string ylabel   = "#sigma_{tot}^{#pip}   [mb] ";

    std::unique_ptr<total_xsection> sigma_tot_pimp( new PDG_parameterization(M_PION, M_PROTON, {-1., 1., 9.56, 1.767, 18.75}, "rpp2020-pimp_total.dat"));
    std::unique_ptr<total_xsection> sigma_tot_pipp( new PDG_parameterization(M_PION, M_PROTON, {+1., 1., 9.56, 1.767, 18.75}, "rpp2020-pipp_total.dat"));

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    auto plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the phase-space for each kinematics

    auto F = [&](double W)
    {
        double s = W*W;
        return sigma_tot_pipp->eval(s);
    };

    std::array<std::vector<double>, 2> x_fx; 
    x_fx = vec_fill(N, F, xmin, xmax, true);
    plotter->AddEntry(x_fx[0], x_fx[1], "#pi^{+} p");

    x_fx[0].clear(); x_fx[1].clear();

    auto G = [&](double W)
    {
        double s = W*W;
        return sigma_tot_pimp->eval(s);
    };

    x_fx = vec_fill(N, G, xmin, xmax, true);

    plotter->AddEntry(x_fx[0], x_fx[1], "#pi^{-} p");

    plotter->SetXaxis(xlabel, 0., xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);

    plotter->SetLegend(0.7, 0.3);
    plotter->SetLegendOffset(0.3, 0.075);

    // Output to file
    plotter->Plot(filename);

    return 1.;
};