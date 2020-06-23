// Parameterization of a resonant amplitude in the s-channel
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _RESONANCE_
#define _RESONANCE_

#include "amplitude.hpp"

// ---------------------------------------------------------------------------
// baryon_resonance class describes the amplitude corresponding to a narrow
// (Breit-Wigner) resonance in the channel. It is parameterized in terms of
// 3 functions:
//
// 1. Breit-Wigner pole for given mass and width
// 2. Hadronic decay coupling to J/psi p final state
// 3. Photo-excitation coupling to the gamma p initial state
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class baryon_resonance : public amplitude
  {
  private:
    int J, P, naturality; // (2xSpin) and parity of the resonance
    double mRes, gamRes; // Resonant mass and width

    int l_min; // lowest allowed relative angular momentum
    double P_t; // Combinatorial factor due to only transverse polarized J/psi contribute

    // Couplings
    double xBR; // Hadronic banching fraction to j/psi p
    double R_photo; // Photocoupling ratio

    // Initial and final CoM momenta evaluated at resonance energy.
    double pi_bar, pf_bar;

  public:
    // Constructor
    baryon_resonance(reaction_kinematics * xkinem, int j, int p, double mass, double width, std::string name = "")
    : amplitude(xkinem, name, 2), mRes(mass), gamRes(width), J(j), P(p),
      naturality(p * pow(-1, (j-1)/2))
    {
      pi_bar = - real(kinematics->initial.momentum("target", mass * mass));
      pf_bar = - real(kinematics->final.momentum("recoil", mass * mass));

      switch (J)
      {
        case 3:
        {
          if (P == -1)
            {l_min = 0; P_t = 2./3.;}
          else if (P == 1)
            {l_min = 1; P_t = 3./5.;}
          break;
        }
        case 5:
        {
          if (P == 1)
            {l_min = 1; P_t = 3./5.;}
          else if (P == -1)
            {l_min = 2; P_t = 1./3.;}
          break;
        }
        default:
        {
          std::cout << "\nbaryon_resonance: spin-parity combination for J = " << J << "/2 and P = " << P << " not available. ";
          std::cout << "Quiting... \n";
          exit(0);
        }
      }
    };

    // Copy Constructor
    baryon_resonance(const baryon_resonance & old)
    : amplitude(old),
      mRes(old.mRes), gamRes(old.gamRes), J(old.J), P(old.P), naturality(old.naturality),
      xBR(old.xBR), R_photo(old.R_photo), pi_bar(old.pi_bar), pf_bar(old.pf_bar)
    {};

    // Setting utility
    void set_params(std::vector<double> params)
    {
      check_Nparams(params);
      xBR = params[0];
      R_photo = params[1];
    };

    // Photoexcitation helicity amplitude for the process gamma p -> R
    std::complex<double> photo_coupling(int lam_i, double s);

    // Hadronic decay helicity amplitude for the R -> J/psi p process
    std::complex<double> hadronic_coupling(int lam_f, double s);

    // Ad-hoc threshold factor to kill the resonance at threshold
    double threshold_factor(double s, double beta);

    // Combined total amplitude including Breit Wigner pole
    std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double t);
  };
};
#endif
