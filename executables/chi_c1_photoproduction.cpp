// ---------------------------------------------------------------------------
// Analytic model for the photoproduction of chi_c1 near threshold at GlueX
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] 10.1103/PhysRevD.96.093008
// ---------------------------------------------------------------------------
// COMMAND LINE OPTIONS:
// -c double          # Change CM angle in degree (default: 0)
// -n int             # Number of points to plot (default: 25)
// -m double          # Maximum CM angle to plot (default: 10 GeV)
// -integ             # Plot integrated xsection (default: false)
// -y "[y1:y2]"       # Custom y bounds in output plot
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
  // ---------------------------------------------------------------------------
  // COMMAND LINE OPTIONS
  // ---------------------------------------------------------------------------

  // Default values
  int N = 25; // how many points to plot
  double theta = 0.;
  double max = 10.;
  bool INTEG = false;
  double y[2] = {0., 2.}; bool custom_y = false;
  std::string filename = "chi_c1_photoproduction.pdf";
  std::string ylabel = "d#sigma/dt  (nb / GeV^{2})";

  // ---------------------------------------------------------------------------
  // Parse command line arguments
  for (int i = 0; i < argc; i++)
  {
    // File name of output
    if (std::strcmp(argv[i],"-f")==0) filename = argv[i+1];
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-n")==0) N = atoi(argv[i+1]);
    if (std::strcmp(argv[i],"-m")==0) max = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-y")==0)
    {
      custom_y = true;
      y_range(argv[i+1], y);
    }
    if (std::strcmp(argv[i],"-integ")==0)
    {
       INTEG = true;
       ylabel = "#sigma  (nb)";
    }
  }

  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------

  // Set up kinematics for the chi_c1
  reaction_kinematics * ptr = new reaction_kinematics(3.510, "chi_c1");

  // Reaction proceeds through multiple vector exchanges.
  // Which we will sum incoherently
  std::vector<amplitude*> exchanges;

  vector_exchange omega(ptr, .780, "#omega");
  omega.set_params({5.2E-4, 16., 0.}); // hadronic and nucleon vector & tensor couplings
  exchanges.push_back(&omega); // Add to the sum vector

  // vector_exchange rho(ptr, .770, "#rho");
  // rho.set_params({9.2E-4, 2.4, 14.6});
  // exchanges.push_back(&rho);
  //
  // vector_exchange phi(ptr, 1.10, "#phi");
  // phi.set_params({4.2E-4, -6.2, 2.1});
  // exchanges.push_back(&phi);
  //
  // vector_exchange jpsi(ptr, 3.097, "J/#psi");
  // jpsi.set_params({1., 3.3E-3, 0.});
  // exchanges.push_back(&jpsi);
  //
  // // The total amplitude with all the above exchanges
  // amplitude_sum total(ptr, exchanges, "Sum");
  // exchanges.push_back(&total);


// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------
// Plotter object
jpacGraph1D* plotter = new jpacGraph1D();

// ---------------------------------------------------------------------------
// Print the desired observable for each amplitude
for (int n = 0; n < exchanges.size(); n++)
{
  auto F = [&](double W)
  {
    if (INTEG == false)
    {
      double t = ptr->t_man(W*W, theta * deg2rad);
      return exchanges[n]->differential_xsection(W*W, t);
    }
    else
    {
      return exchanges[n]->integrated_xsection(W*W);
    }
  };

  std::array<std::vector<double>, 2> x_fx = vec_fill(N, F, sqrt(ptr->sth), max, 1);
  plotter->AddEntry(x_fx[0], x_fx[1], exchanges[n]->identifier);
}

// Set up X-axis
plotter->SetXaxis("W  (GeV)", sqrt(ptr->sth), max);

// To change the range of the Y-axis or the position of the Legend change the arguments here
(custom_y == true) ? (plotter->SetYaxis(ylabel, y[0], y[1])) : (plotter->SetYaxis(ylabel));

// Position of the legend
plotter->SetLegend(0.2, .75);

// Output to file
plotter->Plot(filename);

delete ptr, plotter;
return 1.;
};
