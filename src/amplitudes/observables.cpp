// Abstract class for an amplitude. Used so we can easily build observables
// as the incoherent sum of amplitudes in s, t, and u channels.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/amplitude.hpp"

// ---------------------------------------------------------------------------
// Differential cross section dsigma / dt
// in NANOBARN
double amplitude::differential_xsection(double s, double zs)
{
  double sum = 0.;
  for (int i = 0; i < 24; i++)
  {
    std::complex<double> amp_i = helicity_amplitude(kinematics->helicities[i], s, zs);
    sum += std::real(amp_i * conj(amp_i));
  }

  double norm = 1.;
  norm /= 64. * M_PI * s;
  norm /= real(pow(kinematics->initial.momentum("beam", s), 2.));
  norm /= (2.56819E-6); // Convert from GeV^-2 -> nb

  norm /= 4.; // Average over final state helicites

  return norm * sum;
};

// ---------------------------------------------------------------------------
// Inegrated total cross-section
// IN NANOBARN
double amplitude::integrated_xsection(double s)
{
  auto F = [&](double zs)
  {
    double jacobian = 2.; // 2. * k * q
    jacobian *= real(kinematics->initial.momentum("beam", s));
    jacobian *= real(kinematics->final.momentum(kinematics->vector_particle, s));

    return differential_xsection(s, zs) * jacobian;
  };

  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
  ROOT::Math::Functor1D wF(F);
  ig.SetFunction(wF);

  return ig.Integral(-1,1);
};

// ---------------------------------------------------------------------------
// Polarizatiopn asymmetry between beam and recoil proton
double amplitude::K_LL(double s, double zs)
{
  double sigmapp = 0., sigmapm = 0.;
  for (int i = 0; i < 6; i++)
  {
    std::complex<double> squarepp, squarepm;

    // Amplitudes with lam_gam = + and lam_recoil = +
    squarepp = helicity_amplitude(kinematics->helicities[2*i+1], s, zs);
    squarepp *= conj(squarepp);
    sigmapp += real(squarepp);

    // Amplitudes with lam_gam = + and lam_recoil = -
    squarepm = helicity_amplitude(kinematics->helicities[2*i], s, zs);
    squarepm *= conj(squarepm);
    sigmapm += real(squarepm);
  }

  return (sigmapp - sigmapm) / (sigmapp + sigmapm);
}

// ---------------------------------------------------------------------------
// Polarization asymmetry between beam and target proton
double amplitude::A_LL(double s, double zs)
{
  double sigmapp = 0., sigmapm = 0.;
  for (int i = 0; i < 6; i++)
  {
    std::complex<double> squarepp, squarepm;

    // Amplitudes with lam_gam = + and lam_targ = +
    squarepp = helicity_amplitude(kinematics->helicities[i+6], s, zs);
    squarepp *= conj(squarepp);
    sigmapp += real(squarepp);

    // Amplitudes with lam_gam = + and lam_targ = -
    squarepm = helicity_amplitude(kinematics->helicities[i], s, zs);
    squarepm *= conj(squarepm);
    sigmapm += real(squarepm);
  }

  return (sigmapp - sigmapm) / (sigmapp + sigmapm);
}

// ---------------------------------------------------------------------------
// Photon spin-density matrix elements
std::complex<double> amplitude::SDME(int alpha, int lam, int lamp, double s, double zs)
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

  // Normalization (sum over all amplitudes squared
  double norm = 0.;
  for (int i = 0; i < 24; i++)
  {
    std::complex<double> amp_i = helicity_amplitude(kinematics->helicities[i], s, zs);
    norm += std::real(amp_i * conj(amp_i));
  }

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
    amp_i = helicity_amplitude(kinematics->helicities[index + k], s, zs);
    amp_j = helicity_amplitude(kinematics->helicities[pos_iters[i] + j], s, zs);

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
// Integrated beam asymmetry Sigma_4pi
double amplitude::beam_asymmetry(double s, double zs)
{
  double rho100 = real(SDME(1, 0, 0, s, zs));
  double rho111 = real(SDME(1, 1, 1, s, zs));

  return - rho100 - 2. * rho111;
};

// ---------------------------------------------------------------------------
// Parity asymmetry P_sigma
double amplitude::parity_asymmetry(double s, double zs)
{
  double rho100 = real(SDME(1, 0, 0, s, zs));
  double rho11m1 = real(SDME(1, 1, -1, s, zs));

  return 2. * rho11m1 - rho100;
};
