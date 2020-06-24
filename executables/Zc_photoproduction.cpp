// ---------------------------------------------------------------------------
// Photoproduction of Zc by a  charged pion exchange
// Reproduces the results from [1]
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:1503.02125 [hep-ph]
// ---------------------------------------------------------------------------
// COMMAND LINE OPTIONS:
// -c double          # Change CM angle in degree (default: 0)
// -n int             # Number of points to plot (default: 25)
// -m double          # Maximum CM angle to plot (default: 10 GeV)
// -diff              # Plot differential xsection (default: false)
// -y "[y1:y2]"       # Custom y bounds in output plot
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

  // ---------------------------------------------------------------------------
  // COMMAND LINE OPTIONS
  // ---------------------------------------------------------------------------

  double theta = 0.;
  double max = 25;
  double y[2]; bool custom_y = false;
  int N = 50;
  std::string ylabel = "#sigma  (nb)";
  std::string filename = "Zc_photoproduction.pdf";
  bool INTEG = true;

  // ---------------------------------------------------------------------------
  // Parse command line arguments
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-f")==0) filename = argv[i+1];
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-m")==0) max = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-n")==0) N = atoi(argv[i+1]);
    if (std::strcmp(argv[i],"-y")==0)
    {
      custom_y = true;
      y_range(argv[i+1], y);
    }
    if (std::strcmp(argv[i],"-diff")==0)
    {
      INTEG = false;
      ylabel = "d#sigma/dt  (#mub GeV^{-2})";
    }
  }

  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------

  // For regge amps need pion trajectory
  linear_trajectory alpha(1, -0.7*mPi*mPi, 0.7, "pionic trajectory");

  // ---------------------------------------------------------------------------
  // ZC(3900)

  // Kinematics for 3900 MeV vector
  reaction_kinematics * ptr3900 = new reaction_kinematics(3.9, "Z_{c}^{+}(3900)");

  // Amplitudes
  pseudoscalar_exchange Z3900(ptr3900, mPi, "Z_{c}^{+}(3900), #pi exchange");
  pseudoscalar_exchange Z3900R(ptr3900, &alpha, "Z_{c}^{+}(3900), Reggeon exchange");

  // Couplings for 4 MeV width
  Z3900.set_params({0.67 * 3.90, sqrt(4.*M_PI*14.4)});
  Z3900R.set_params({0.67 * 3.90, sqrt(4.*M_PI*14.4)});


  // ---------------------------------------------------------------------------
  // ZC(4200)

  // Kinematics for 4200
  reaction_kinematics * ptr4200 = new reaction_kinematics(4.20, "Z_{c}^{+}(4200)");

  // Amplitudes
  pseudoscalar_exchange Z4200(ptr4200, mPi, "Z_{c}^{+}(4200), #pi exchange");
  pseudoscalar_exchange Z4200R(ptr4200, &alpha, "Z_{c}^{+}(4200), Reggeon exchange");

  // Couplings
  Z4200.set_params({1.731 * 4.20, sqrt(4.*M_PI*14.4)});
  Z4200R.set_params({1.731 * 4.20, sqrt(4.*M_PI*14.4)});
  // ---------------------------------------------------------------------------

  // Add to a vector to plot them both
  std::vector<amplitude*> amps;
  amps.push_back(&Z4200);
  amps.push_back(&Z4200R);
  amps.push_back(&Z3900);
  amps.push_back(&Z3900R);

  // ---------------------------------------------------------------------------
  // You shouldnt need to change anything below this line
  // ---------------------------------------------------------------------------

  // Plotter object
  jpacGraph1D* plotter = new jpacGraph1D();

  // ---------------------------------------------------------------------------
  // Print the desired observable for each amplitude
  for (int n = 0; n < amps.size(); n++)
  {
    std::cout << std::endl << "Printing amplitude: " << amps[n]->identifier << "\n";
    auto F = [&](double W)
    {
      if (INTEG == false)
      {
        double t = amps[n]->kinematics->t_man(W*W, theta * deg2rad);
        return amps[n]->differential_xsection(W*W, t);
      }
      else
      {
        return amps[n]->integrated_xsection(W*W);
      }
    };


    std::array<std::vector<double>, 2> x_fx = vec_fill(N, F, sqrt(amps[n]->kinematics->sth) + EPS, max, true);
    plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->identifier);
  }

  // Set up X-axis
  plotter->SetXaxis("W  (GeV)", sqrt(ptr3900->sth), max);

  // To change the range of the Y-axis or the position of the Legend change the arguments here
  (custom_y == true) ? (plotter->SetYaxis(ylabel, y[0], y[1])) : (plotter->SetYaxis(ylabel));

  // Position of the legend
  plotter->SetLegend(0.5, 0.65);

  // Output to file
  plotter->Plot(filename);

  delete ptr4200, ptr3900, plotter;
  return 1.;
}
