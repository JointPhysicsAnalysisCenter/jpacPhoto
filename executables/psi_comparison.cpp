// ---------------------------------------------------------------------------
// Comparing the photoproduction cross-sections of the Jpsi 1s and 2s states
// near threshold at GlueX
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] 10.1103/PhysRevD.100.034019
// ---------------------------------------------------------------------------
// COMMAND LINE OPTIONS:
// -c double          # Change CM angle in degree (default: 0)
// -n int             # Number of points to plot (default: 100)
// -m double          # Maximum energy to plot (default: 5.55 GeV)
// -ratio             # Plot ratio of 2S / 1S dxs (default: false)
// -y "[y1:y2]"       # Custom y bounds in output plot
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "regge_trajectory.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
  double theta = 0.;
  double max = 5.55;
  int N = 100; // how many points to plot
  bool RATIO = false;
  std::string ylabel = "d#sigma/dt  (nb GeV^{-2})";
  double y[2]; bool custom_y = false;
  std::string filename = "psi_photoproduction.pdf";
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-n")==0) N = atoi(argv[i+1]);
    if (std::strcmp(argv[i],"-m")==0) max = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-y")==0) {y_range(argv[i+1], y); custom_y = true;}
    if (std::strcmp(argv[i],"-ratio")==0)
    {
      RATIO = true;
      ylabel = "100 x #sigma(2S) / #sigma(1S)";
    }
  }


  // Set up pomeron trajectory
  // Here we use (real) linear trajectory with intercept and slope only free params
  // Best fit values from [1]
  linear_trajectory alpha(+1, 0.941, 0.364, "pomeron");

  // ---------------------------------------------------------------------------
  // PSI(1S)
  // ---------------------------------------------------------------------------

  // Set up kinematics, determined entirely by vector meson mass
  reaction_kinematics * ptr1S = new reaction_kinematics(3.097, "#psi(1S)");

  // Create amplitudes with kinematics and trajectory
  pomeron_exchange pomeron_1S(ptr1S, &alpha, false, "#psi(1S)");

  // Feed in other two parameters (normalization and t-slope)
  // Best fit from [1]
  pomeron_1S.set_params({0.379, 0.12});

  // ---------------------------------------------------------------------------
  // PSI(2S)
  // ---------------------------------------------------------------------------

  // higher mass kinematics
  reaction_kinematics * ptr2S = new reaction_kinematics(3.686, "#psi(2S)");

  // trajectory is the same but new kinematics for bigger phase-space
  pomeron_exchange pomeron_2S(ptr2S, &alpha, false, "100 x #psi(2S)");

  // Same t-slope but coupling is scaled by 1/4, multiplied by 10 for comparision
  pomeron_2S.set_params({10. * 0.379 / 4., 0.12});

  // ---------------------------------------------------------------------------
  // You shouldnt need to change anything below this line
  // ---------------------------------------------------------------------------

  // Amplitudes to plot
  std::vector<amplitude*> amps = {&pomeron_1S, &pomeron_2S};

  // Initialize objects to plot
  jpacGraph1D* plotter = new jpacGraph1D();

  // ---------------------------------------------------------------------------
  // Print the desired observable for each amplitude
  if (RATIO == false)
  {
    for (int n = 0; n < amps.size(); n++)
    {
      auto F = [&](double W)
      {
        double t = amps[n]->kinematics->t_man(W*W, theta);
        double result = amps[n]->differential_xsection(W*W, t);

        return result;
      };

      std::array<std::vector<double>, 2> x_fx = vec_fill(N, F, sqrt(amps[n]->kinematics->sth) + 0.01, max);
      plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->identifier);
    }
  }
  else
  {
    auto F = [&](double W)
    {
      double result  = pomeron_2S.differential_xsection(W*W, pomeron_2S.kinematics->t_man(W*W, theta));
             result /= pomeron_1S.differential_xsection(W*W, pomeron_1S.kinematics->t_man(W*W, theta));

      return result;
    };

    std::array<std::vector<double>, 2> x_fx = vec_fill(N, F, sqrt(pomeron_2S.kinematics->sth) + 0.01, max);
    plotter->SetLegend(false);
    plotter->AddEntry(x_fx[0], x_fx[1], "ratio");
  }

// ---------------------------------------------------------------------------
// Plotting Settings
// ---------------------------------------------------------------------------

// Tweak the axes
(custom_y == true) ? (plotter->SetYaxis(ylabel, y[0], y[1])) : (plotter->SetYaxis(ylabel));

plotter->SetLegend(0.2, 0.7);
plotter->SetXaxis("W   (GeV)", sqrt(ptr2S->sth), max);
plotter->Plot(filename);

  return 1.;
};