// Axial-vector meson photoproduction proceeding through a vector meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/vector_exchange.hpp"

// ---------------------------------------------------------------------------
// Assemble the helicity amplitude by contracting the lorentz indices
std::complex<double> jpacPhoto::vector_exchange::helicity_amplitude(std::vector<int> helicities, double xs, double xt)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  // Update the saved energies and angles
  s = xs; t = xt;
  theta = kinematics->theta_s(xs, xt);
;
  std::complex<double> result = 0.;

  if (REGGE == false) // Use the covariant expression
  {
    // Need to contract the Lorentz indices
    for (int mu = 0; mu < 4; mu++)
    {
      for(int nu = 0; nu < 4; nu++)
      {
        std::complex<double> temp;
        temp  = top_vertex(mu, lam_gam, lam_vec);
        temp *= metric[mu];
        temp *= vector_propagator(mu, nu);
        temp *= metric[nu];
        temp *= bottom_vertex(nu, lam_targ, lam_rec);

        result += temp;
      }
    }
  }
  else // Use the helicity amplitude form in the t-channel
  {
    // NOTE THIS ONLY WORKS FOR UNPOLARIZED CROSS-SECTIONS
    // NEED TO WIGNER-ROTATE HELICITES TO S CHANNEL FOR SDMES

    // TODO: ADD CROSSING-MATRICES

    int lam  = lam_gam - lam_vec;
    int lamp = (lam_targ - lam_rec) / 2.;

    // Product of residues
    result  = top_residue(lam_gam, lam_vec);
    result *= bottom_residue(lam_targ, lam_rec);
    result *= regge_propagator(1, lam, lamp);
  }

  if (IF_FF == true)
  {
    double tprime = t - kinematics->t_man(s, 0.);
    result *= exp(b * tprime);
  }

  return result;
};

// ---------------------------------------------------------------------------
// REGGE EVALUATION
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Analytic residues for Regge form

// Photon - Axial - Vector
std::complex<double> jpacPhoto::vector_exchange::top_residue(int lam_gam, int lam_vec)
{
  // TODO: Explicit phases in terms of lam_gam and lam_vec instead of difference
  int lam = lam_gam - lam_vec;

  std::complex<double> result;
  switch (std::abs(lam))
  {
    case 0:
    {
      result = 1.;
      break;
    }
    case 1:
    {
      result = sqrt(xr * t) / kinematics->mVec;
      break;
    }
    case 2:
    {
      return gpGam * (t / kinematics->mVec2 - 1.);
    }
    default:
    {
      std::cout << "\nvector_exchange: invalid helicity flip lambda = " << lam << ". Quitting... \n";
      exit(0);
    }
  }

  std::complex<double> q = (t - kinematics->mVec2) / sqrt(4. * t * xr);
  return  result * q * gGam;
};

// Nucleon - Nucleon - Vector
std::complex<double> jpacPhoto::vector_exchange::bottom_residue(int lam_targ, int lam_rec)
{
  // TODO: Explicit phases in terms of lam_targ and lam_rec instead of difference
  int lamp = (lam_targ - lam_rec) / 2.;

  std::complex<double> vector, tensor;
  switch (std::abs(lamp))
  {
    case 0:
    {
      vector =  1.;
      tensor = sqrt(xr * t) / (2. * mPro);
      break;
    }
    case 1:
    {
      vector = sqrt(2.) * sqrt(xr * t) / (2. * mPro);
      tensor = sqrt(2.);
      break;
    }
    case 2:
    {
      return 0.;
    }
    default:
    {
      std::cout << "\nreggeon_exchange: invalid helicity flip lambda^prime = " << lamp << ". Quitting... \n";
      exit(0);
    }
  }

  std::complex<double> result;
  result = gV * vector + gT * tensor * sqrt(xr * t) / (2. * mPro);
  result *= 2. * mPro;

  return result;
};

// ---------------------------------------------------------------------------
// Reggeon Propagator
std::complex<double> jpacPhoto::vector_exchange::regge_propagator(int j, int lam, int lamp)
{
  int M = std::max(std::abs(lam), std::abs(lamp));

  if (M > j)
  {
    return 0.;
  }

  std::complex<double> alpha_t = alpha->eval(t);

  // the gamma function causes problesm for large t so
  if (std::abs(alpha_t) > 30.)
  {
    return 0.;
  }
  else
  {
    std::complex<double> result;
    result  = wigner_leading_coeff(j, lam, lamp);
    result /= barrier_factor(j, M);
    result *= half_angle_factor(lam, lamp);

    result *= - alpha->slope();
    result *= 0.5 * (double(alpha->signature) + exp(-xi * M_PI * alpha_t));
    result *= cgamma(1. - alpha_t);
    result *= pow(s, alpha_t - double(M));

    return result;
  }
};

//------------------------------------------------------------------------------
// Half angle factors
std::complex<double> jpacPhoto::vector_exchange::half_angle_factor(int lam, int lamp)
{
  std::complex<double> zt = kinematics->z_t(s, t);

  std::complex<double> sinhalf = sqrt((xr - zt) / 2.);
  std::complex<double> coshalf = sqrt((xr + zt) / 2.);

  std::complex<double> result;
  result  = pow(sinhalf, double(std::abs(lam - lamp)));
  result *= pow(coshalf, double(std::abs(lam + lamp)));

  return result;
};

//------------------------------------------------------------------------------
// Angular momentum barrier factor
std::complex<double> jpacPhoto::vector_exchange::barrier_factor(int j, int M)
{
  std::complex<double> q = (t - kinematics->mVec2) / sqrt(4. * t * xr);
  std::complex<double> p = sqrt(xr * t - 4.*mPro2) / 2.;

  std::complex<double> result = pow(xr * 2. * p * q, double(j - M));

  return result;
};

// ---------------------------------------------------------------------------
// FEYNMAN EVALUATION
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Photon - Axial Vector - Vector vertex
std::complex<double> jpacPhoto::vector_exchange::top_vertex(int mu, int lam_gam, int lam_vec)
{
  // Contract with LeviCivita
  std::complex<double> result = 0.;
  for (int alpha = 0; alpha < 4; alpha++)
  {
    for (int beta = 0; beta < 4; beta++)
    {
      for (int gamma = 0; gamma < 4; gamma++)
      {
        std::complex<double> temp;
        temp = levi_civita(mu, alpha, beta, gamma);
        temp *= metric[mu];
        temp *= kinematics->initial.component(alpha, "beam", s, 0.);
        temp *= kinematics->eps_gamma.component(beta, lam_gam, s, 0.);
        temp *= kinematics->eps_vec.component(gamma, lam_vec, s, theta);

        result += temp;
      }
    }
  }

  // Multiply by coupling
  return result * gGam;
};

// ---------------------------------------------------------------------------
// Nucleon - Nucleon - Vector vertex
std::complex<double> jpacPhoto::vector_exchange::bottom_vertex(int mu, int lam_targ, int lam_rec)
{
  // Vector coupling piece
  std::complex<double> vector = 0.;
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      std::complex<double> temp;
      temp = kinematics->recoil.adjoint_component(i, lam_rec, s, theta + M_PI); // theta_rec = theta + pi
      temp *= gamma_matrices[mu][i][j];
      temp *= kinematics->target.component(j, lam_targ, s, M_PI); // theta_targ = pi

      vector += temp;
    }
  }

  // Tensor coupling piece
  std::complex<double> tensor = 0.;
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      std::complex<double> sigma_q_ij = 0.;
      for (int nu = 0; nu < 4; nu++)
      {
        sigma_q_ij += sigma(mu, nu, i, j) * metric[nu] * exchange_momenta(nu) / (2. * mPro);
      }

      std::complex<double> temp;
      temp = kinematics->recoil.adjoint_component(i, lam_rec, s, theta + M_PI); // theta_rec = theta + pi
      temp *= sigma_q_ij;
      temp *= kinematics->target.component(j, lam_targ, s, M_PI); // theta_targ = pi

      tensor += temp;
    }
  }

  return gV * vector - gT * tensor;
};

// ---------------------------------------------------------------------------
// Four-momentum of the exchanged meson.
// Simply the difference of the photon and axial 4-momenta
std::complex<double> jpacPhoto::vector_exchange::exchange_momenta(int mu)
{
  std::complex<double> qGamma_mu, qA_mu;
  qGamma_mu = kinematics->initial.component(mu, "beam", s, 0.);
  qA_mu = kinematics->final.component(mu, kinematics->vector_particle, s, theta);

  return (qGamma_mu - qA_mu);
};

// ---------------------------------------------------------------------------
// Propagator of a massive spin-one particle
std::complex<double> jpacPhoto::vector_exchange::vector_propagator(int mu, int nu)
{
  // q_mu q_nu / mEx2 - g_mu nu
  std::complex<double> result;
  result = exchange_momenta(mu) * exchange_momenta(nu) / mEx2;

  if (mu == nu)
  {
    result -= metric[mu];
  }

  result /= t - mEx2;

  return result;
};