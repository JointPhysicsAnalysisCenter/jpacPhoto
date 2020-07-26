// ---------------------------------------------------------------------------
// Photoproduction of Z states by a  charged pion exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:1503.02125 [hep-ph]
// [2] arXiv:1806.08414 [hep-ph]
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
  // Preliminaries
  // ---------------------------------------------------------------------------

  double g_NN = sqrt(4.* M_PI * 14.4); // Nucleon coupling same for all
  // double bPi = 4.557; // GeV^-2
  double bPi = 1. / (.9 * .9); 

  // Zc(3900)
  double mZc = 3.8884; // GeV
  reaction_kinematics * kZc = new reaction_kinematics(mZc, "Z_{c}(3900)");

  double gc_Psi = 1.91; // psi coupling before VMD scaling
  double gc_Gamma = e * fJpsi * gc_Psi / mJpsi;
  std::vector<double> Zc_couplings = {gc_Gamma, g_NN};

  // Zb(10610)
  double mZb = 10.6072;
  reaction_kinematics * kZb = new reaction_kinematics(mZb, "Z_{b}(10610)");

  double gb_Ups1 = 0.49, gb_Ups2 = 3.30, gb_Ups3 = 9.22;
  double gb_Gamma = e * (fUpsilon1S * gb_Ups1 / mUpsilon1S 
                       + fUpsilon2S * gb_Ups2 / mUpsilon2S
                       + fUpsilon3S * gb_Ups3 / mUpsilon3S);  
  std::vector<double> Zb_couplings = {gb_Gamma, g_NN};

  
  // Zb(10650)
  double mZbp = 10.6522;
  reaction_kinematics * kZbp = new reaction_kinematics(mZbp, "Z_{b}(10650)");

  double gbp_Ups1 = 0.21, gbp_Ups2 = 1.47, gbp_Ups3 = 4.8;
  double gbp_Gamma = e * (fUpsilon1S * gbp_Ups1 / mUpsilon1S 
                        + fUpsilon2S * gbp_Ups2 / mUpsilon2S
                        + fUpsilon3S * gbp_Ups3 / mUpsilon3S);  
  std::vector<double> Zbp_couplings = {gbp_Gamma, g_NN};
  
  // ---------------------------------------------------------------------------
  // Fixed-spin amplitudes
  // ---------------------------------------------------------------------------

  pseudoscalar_exchange Zc_fixedspin(kZc, mPi, "#it{Z_{c}}(3900)^{+}");
  Zc_fixedspin.set_params(Zc_couplings);
  Zc_fixedspin.set_formfactor(true, bPi);

  pseudoscalar_exchange Zb_fixedspin(kZb, mPi,  "#it{Z_{b}}(10610)^{+}");
  Zb_fixedspin.set_params(Zb_couplings);
  Zb_fixedspin.set_formfactor(true, bPi);

  pseudoscalar_exchange Zbp_fixedspin(kZbp, mPi, "#it{Z_{b}}(10650)^{+}");
  Zbp_fixedspin.set_params(Zbp_couplings);
  Zbp_fixedspin.set_formfactor(true, bPi);

  // ---------------------------------------------------------------------------
  // Regge-ized amplitudes
  // ---------------------------------------------------------------------------

  // Need pion trajectory
  int signature = +1;
  double alpha_prime = 0.7;
  double alpha_0 = - alpha_prime * mPi2;
  linear_trajectory * alpha = new linear_trajectory(signature, alpha_0, alpha_prime, "pion");
    
  pseudoscalar_exchange Zc_regge(kZc, alpha, "#it{Z_{c}}(3900)^{+}");
  Zc_regge.set_params(Zc_couplings);
  Zc_regge.set_formfactor(true, bPi);

  pseudoscalar_exchange Zb_regge(kZb, alpha,  "#it{Z_{b}}(10610)^{+}");
  Zb_regge.set_params(Zb_couplings);
  Zb_regge.set_formfactor(true, bPi);

  pseudoscalar_exchange Zbp_regge(kZbp, alpha,  "#it{Z_{b}}(10650)^{+}");
  Zbp_regge.set_params(Zbp_couplings);
  Zbp_regge.set_formfactor(true, bPi);

  // ---------------------------------------------------------------------------
  // Plotting options
  // ---------------------------------------------------------------------------
  
  // which amps to plot
  std::vector<amplitude*> amps;
  amps.push_back(&Zc_fixedspin);
  amps.push_back(&Zb_fixedspin);
  amps.push_back(&Zbp_fixedspin);
  // amps.push_back(&Zc_regge);
  // amps.push_back(&Zb_regge);
  // amps.push_back(&Zbp_regge);

  int N = 10;
 
  // ---------------------------------------------------------------------------
  double  xmin = 4;
  double  xmax = 20.;

  double  ymin = 2.E-2;
  double  ymax = 100;
  // ---------------------------------------------------------------------------
  // double  xmin = 20.;
  // double  xmax = 70.;

  // double  ymax = 1.;
  // double  ymin = 2.E-6;
  // ---------------------------------------------------------------------------

  std::string ylabel  = ROOT_italics("#sigma(#gamma p #rightarrow Z n)") + "  [nb]";
 
  std::string filename = "Z_FS.pdf";
  // std::string filename = "Z_regge.pdf";

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

    auto F = [&](double x)
    {
      return amps[n]->integrated_xsection(x*x);
    };

    std::array<std::vector<double>, 2> x_fx, x_fx1;
    if (xmin < amps[n]->kinematics->Wth)
    {
        double match =  amps[n]->kinematics->Wth + 3.;
        x_fx = vec_fill(25, F, amps[n]->kinematics->Wth + EPS, match, true);
        x_fx1 = vec_fill(N, F, match, xmax, true);

        for (int i = 0; i < x_fx1[0].size(); i++)
        {
          x_fx[0].push_back(x_fx1[0][i]);
          x_fx[1].push_back(x_fx1[1][i]);
        }

        x_fx[0].insert(x_fx[0].begin(), amps[n]->kinematics->Wth);
        x_fx[1].insert(x_fx[1].begin(), 0.);
    }
    else
    {
        x_fx = vec_fill(N, F, xmin, xmax, true);
    }

    plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->identifier);
  }

  plotter->SetXaxis(ROOT_italics("W_{#gammap}") +"   [GeV]", xmin, xmax);

  plotter->SetYaxis(ylabel, ymin, ymax);
  plotter->SetYlogscale(1);

  // plotter->SetLegend(0.2, 0.2);
  plotter->SetLegend(0.7, 0.67);


  // Output to file
  plotter->Plot(filename);

  delete kZc, kZb, kZbp, plotter, alpha;

  return 1.;
}
