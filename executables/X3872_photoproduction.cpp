// ---------------------------------------------------------------------------
// Prediction for the photoproduction of X(3872) or Chi_c1(3872) at the future
// EIC at low momentum transfer and high center of mass energy.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// COMMAND LINE OPTIONS:
// -c double          # Change CM angle in degree (default: 0)
// -n int             # Number of points to plot (default: 25)
// -m double          # Maximum energy to plot (default: 10 GeV)
// -integ             # Plot integrated xsection (default: false)
// -y "[y1:y2]"       # Custom y bounds in output plot
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"

#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/pomeron_exchange.hpp"
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
  int N = 50;
  double max = 10.;
  double theta = 0.;
  double y[2]; bool custom_y = false;
  bool INTEG = false; std::string ylabel = "d#sigma/dt  (nb GeV^{-2})";
  std::string filename = "X3872_photoproduction.pdf";

  // ---------------------------------------------------------------------------
  // Parse command line arguments
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-f")==0) filename = argv[i+1];
    if (std::strcmp(argv[i],"-y")==0) {y_range(argv[i+1], y); custom_y = true;}
    if (std::strcmp(argv[i],"-m")==0) max = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-n")==0) N = atoi(argv[i+1]);
    if (std::strcmp(argv[i],"-integ")==0)
    {
      INTEG = true;
      ylabel = "#sigma  (nb)";
    }
  }

  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------

  // Set up kinematics for the X(3872)
  reaction_kinematics * ptr = new reaction_kinematics(3.872, "X(3872)");

  // Sum individual contributions together incoherently
  amplitude_sum total(ptr, "Sum");

  // Linear trajectory for the rho
  linear_trajectory alpha(-1, 0.5, 0.9, "EXD_linear");

  // Initialize Reggeon amplitude with the above kinematics and regge_trajectory
  vector_exchange rho(ptr, &alpha, "#rho");
  rho.set_params({3.81E-3, 2.4, 14.6});
  total.add_amplitude(&rho);

  vector_exchange omega(ptr, &alpha, "#omega");
  omega.set_params({9.51E-3, 16, 0.});
  total.add_amplitude(&omega);

  // Vector of amplitudes to plot
  std::vector<amplitude*> exchanges;
  exchanges.push_back(&total);
  exchanges.push_back(&rho);
  exchanges.push_back(&omega);

  // ---------------------------------------------------------------------------
  // You shouldnt need to change anything below this line
  // ---------------------------------------------------------------------------

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
        double t = ptr->t_man(W*W, theta);
        return exchanges[n]->differential_xsection(W*W, t);
      }
      else
      {
        return exchanges[n]->integrated_xsection(W*W);
      }
    };

    std::array<std::vector<double>, 2> x_fx = vec_fill(N, F, sqrt(ptr->sth) + EPS, max);
    plotter->AddEntry(x_fx[0], x_fx[1], exchanges[n]->identifier);
  }

  // ---------------------------------------------------------------------------
  // Plotting Settings
  // ---------------------------------------------------------------------------

  // Tweak the axes
  (custom_y == true) ? (plotter->SetYaxis(ylabel, y[0], y[1])) : (plotter->SetYaxis(ylabel));

  plotter->SetXaxis("W   (GeV)", sqrt(ptr->sth), max);

  plotter->SetLegend(0.2, 0.7);
  plotter->Plot(filename);
}
