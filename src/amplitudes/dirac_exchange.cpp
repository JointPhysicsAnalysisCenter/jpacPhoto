// Spin-1/2 exchange ampltiude from perturbation theory
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/dirac_exchange.hpp"

//------------------------------------------------------------------------------
// Combine everything and contract indices
std::complex<double> jpacPhoto::dirac_exchange::helicity_amplitude(std::vector<int> helicities, double xs, double xt)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  // Store the invariant energies to avoid having to pass them around 
  s = xs; t = xt, theta = kinematics->theta_s(xs, xt);

  std::complex<double> result = 0.;
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      std::complex<double> temp;
      temp  = top_vertex(i, lam_gam, lam_rec);
      temp *= dirac_propagator(i, j);
      temp *= bottom_vertex(j, lam_vec, lam_targ);

      result += temp;
    }
  }

  return result;
};


//------------------------------------------------------------------------------
// Photon fermion fermion vertex
// (ubar epsilon-slashed)
std::complex<double> jpacPhoto::dirac_exchange::top_vertex(int i, int lam_gam, int lam_rec)
{
  if (ScTOP == true)
  {
    // Scalar for testing purposes
    return gGam * kinematics->recoil->adjoint_component(i, lam_rec, s, theta + M_PI);
  }

  std::complex<double> result = 0.;
  for (int k = 0; k < 4; k++)
  {
    std::complex<double> temp;
    temp  = kinematics->recoil->adjoint_component(k, lam_rec, s, theta + M_PI); // theta_recoil = theta + pi
    temp *= slashed_eps(k, i, lam_gam, kinematics->eps_gamma, false, s, 1.); // theta_gamma = 0

    result += temp;
  }

  return gGam * result;
};

//------------------------------------------------------------------------------
// Vector fermion fermion vertex
// (epsilon*-slashed u)
std::complex<double> jpacPhoto::dirac_exchange::bottom_vertex(int j, int lam_vec, int lam_targ)
{
  if (ScBOT == true)
  {
    // Scalar for testing purposes
    return gVec * kinematics->target->component(j, lam_targ, s , M_PI); // theta_targer = pi
  }

  std::complex<double> result = 0.;
  for (int k = 0; k < 4; k++)
  {
    std::complex<double> temp;
    temp  = slashed_eps(j, k, lam_vec, kinematics->eps_vec, true, s, theta); //theta_vec = theta
    temp *= kinematics->target->component(k, lam_targ, s , M_PI); // theta_target = pi

    result += temp;
  }

  return gVec * result;
};

//------------------------------------------------------------------------------
double jpacPhoto::dirac_exchange::exchange_mass()
{
    double result = 0.;
    for (int mu = 0; mu < 4; mu++)
    {
      std::complex<double> temp;
      temp  = exchange_momentum(mu);
      temp *= metric[mu];
      temp *= exchange_momentum(mu);

      result += real(temp);
    }
    return result;
}

std::complex<double> jpacPhoto::dirac_exchange::exchange_momentum(int mu)
{
  std::complex<double> qGamma_mu, qRec_mu;
  qGamma_mu   = kinematics->initial->q(mu, s, M_PI);
  qRec_mu     = kinematics->final->p(mu, s, theta + M_PI);

  return qRec_mu - qGamma_mu;
};

std::complex<double> jpacPhoto::dirac_exchange::slashed_exchange_momentum(int i, int j)
{
  std::complex<double> result = 0.;
  for (int mu = 0; mu < 4; mu++)
  {
    std::complex<double> temp;
    temp  = gamma_matrices[mu][i][j];
    temp *= metric[mu];
    temp *= exchange_momentum(mu);

    result += temp;
  }

  return result;
};

//------------------------------------------------------------------------------
// Slashed polarization vectors
std::complex<double> jpacPhoto::dirac_exchange::slashed_eps(int i, int j, double lam, polarization_vector * eps, bool STAR, double s, double theta)
{
  std::complex<double> result = 0.;
  for (int mu = 0; mu < 4; mu++)
  {
    std::complex<double> temp;
    if (STAR == false)
    {
      temp  = eps->component(mu, lam, s, theta);
    }
    else
    {
      temp = eps->conjugate_component(mu, lam, s, theta);
    }
    temp *= metric[mu];
    temp *= gamma_matrices[mu][i][j];

    result += temp;
  }

  return result;
};


//------------------------------------------------------------------------------
std::complex<double> jpacPhoto::dirac_exchange::dirac_propagator(int i, int j)
{
  std::complex<double> result;
  result = slashed_exchange_momentum(i, j);

  if (i == j)
  {
    result += mEx;
  }

  result /= exchange_mass() - mEx2;

  return result;
};
