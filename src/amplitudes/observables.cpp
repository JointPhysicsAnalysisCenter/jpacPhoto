// Abstract class for an amplitude. Used so we can easily build observables
// as the incoherent sum of amplitudes in s, t, and u channels.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/amplitude.hpp"

// ---------------------------------------------------------------------------
// Square of the spin averaged amplitude squared
double jpacPhoto::amplitude::probability_distribution(double s, double t)
{
  double sum = 0.;
  for (int i = 0; i < 24; i++)
  {
    std::complex<double> amp_i = helicity_amplitude(kinematics->helicities[i], s, t);
    sum += std::real(amp_i * conj(amp_i));
  }

  sum /= 4.; // Average over initial state helicites
  return sum;
};

// ---------------------------------------------------------------------------
// Differential cross section dsigma / dt
// in NANOBARN
double jpacPhoto::amplitude::differential_xsection(double s, double t)
{
  double sum = probability_distribution(s, t);

  double norm = 1.;
  norm /= 64. * M_PI * s;
  norm /= real(pow(kinematics->initial->momentum(s), 2.));
  norm /= (2.56819E-6); // Convert from GeV^-2 -> nb

  return norm * sum;
};

// ---------------------------------------------------------------------------
// Inegrated total cross-section
// IN NANOBARN
double jpacPhoto::amplitude::integrated_xsection(double s)
{
  auto F = [&](double t)
  {
    return differential_xsection(s, t);
  };

  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
  ROOT::Math::Functor1D wF(F);
  ig.SetFunction(wF);

  double t_min = kinematics->t_man(s, 0.);
  double t_max = kinematics->t_man(s, M_PI);

  return ig.Integral(t_max, t_min);
};

// ---------------------------------------------------------------------------
// Polarizatiopn asymmetry between beam and recoil proton
double jpacPhoto::amplitude::K_LL(double s, double t)
{
  double sigmapp = 0., sigmapm = 0.;
  for (int i = 0; i < 6; i++)
  {
    std::complex<double> squarepp, squarepm;

    // Amplitudes with lam_gam = + and lam_recoil = +
    squarepp = helicity_amplitude(kinematics->helicities[2*i+1], s, t);
    squarepp *= conj(squarepp);
    sigmapp += real(squarepp);

    // Amplitudes with lam_gam = + and lam_recoil = -
    squarepm = helicity_amplitude(kinematics->helicities[2*i], s, t);
    squarepm *= conj(squarepm);
    sigmapm += real(squarepm);
  }

  return (sigmapp - sigmapm) / (sigmapp + sigmapm);
}

// ---------------------------------------------------------------------------
// Polarization asymmetry between beam and target proton
double jpacPhoto::amplitude::A_LL(double s, double t)
{
  double sigmapp = 0., sigmapm = 0.;
  for (int i = 0; i < 6; i++)
  {
    std::complex<double> squarepp, squarepm;

    // Amplitudes with lam_gam = + and lam_targ = +
    squarepp = helicity_amplitude(kinematics->helicities[i+6], s, t);
    squarepp *= conj(squarepp);
    sigmapp += real(squarepp);

    // Amplitudes with lam_gam = + and lam_targ = -
    squarepm = helicity_amplitude(kinematics->helicities[i], s, t);
    squarepm *= conj(squarepm);
    sigmapm += real(squarepm);
  }

  return (sigmapp - sigmapm) / (sigmapp + sigmapm);
}

// ---------------------------------------------------------------------------
// Photon spin-density matrix elements
std::complex<double> jpacPhoto::amplitude::SDME(int alpha, int lam, int lamp, double s, double t)
{
  if (alpha < 0 || alpha > 2 || std::abs(lam) > 1 || std::abs(lamp) > 1)
  {
    std::cout << "\nError! Invalid parameter passed to SDME. Quitting...\n";
    exit(1);
  };

  // Phase and whether to conjugate at the end
  bool CONJ = false;
  double phase = 1.;

  // if first index smaller, switch them
  if (std::abs(lam) < std::abs(lamp))
  {
    int temp = lam;
    lam = lamp;
    lamp = temp;

    CONJ = true;
  }

  // if first index is negative, flip to positive
  if (lam < 0)
  {
    lam *= -1;
    lamp *= -1;

    phase *= pow(-1., double(lam - lamp));
  }

  // Normalization (sum over all amplitudes squared)
  double norm = probability_distribution(s, t);
  norm *= norm;

  // These are the indexes of the amplitudes in reaction_kinematics that have
  // lambda_V = +1
  std::vector<int> pos_iters = {0, 1, 6, 7, 12, 13, 18, 19};
  std::vector<int> neg_iters = {12, 13, 18, 19, 0, 1, 6, 7};

  // j and k filter the right helicity combinations for 00, 1-1, 10, 11
  int j, k;
  (lam == 0) ? (k = 2) : (k = 0);

  switch (lamp)
  {
    case -1: { j = 4; break; }
    case 0:  { j = 2; break; }
    case 1:  { j = 0; break; }
    default:
    {
     std::cout << "\nSDME: Invalid parameter. J/Psi helicity projection lamp = 0 or 1.";
     std::cout << " Quitting... \n";
     exit(1);
    }
  }

  // Sum over the appropriate amplitude combinations
  std::complex<double> result = 0.;
  for (int i = 0; i < 8; i++)
  {
    int index;
    (alpha == 0) ? (index = pos_iters[i]) : (index = neg_iters[i]);

    std::complex<double> amp_i, amp_j;
    amp_i = helicity_amplitude(kinematics->helicities[index + k], s, t);
    amp_j = helicity_amplitude(kinematics->helicities[pos_iters[i] + j], s, t);

    (alpha == 2) ? (amp_j *= xi * double(kinematics->helicities[pos_iters[i] + j][0])) : (amp_j *= xr);

    result += real(amp_i * conj(amp_j));
  }

  if (CONJ == true)
  {
    result = conj(result);
  }

  result /= norm;
  result *= phase;

  return result;
};

// ---------------------------------------------------------------------------
// Integrated beam asymmetry Sigma
double jpacPhoto::amplitude::beam_asymmetry(double s, double t)
{
  double rho100 = real(SDME(1, 0, 0, s, t));
  double rho111 = real(SDME(1, 1, 1, s, t));

  return - rho100 - 2. * rho111;
};

// ---------------------------------------------------------------------------
// Parity asymmetry P_sigma
double jpacPhoto::amplitude::parity_asymmetry(double s, double t)
{
  double rho100 = real(SDME(1, 0, 0, s, t));
  double rho11m1 = real(SDME(1, 1, -1, s, t));

  return 2. * rho11m1 - rho100;
};
