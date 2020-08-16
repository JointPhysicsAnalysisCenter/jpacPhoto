// Spin-3/2 exchange amplitude from perturbation theory
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/rarita_exchange.hpp"

//------------------------------------------------------------------------------
// Combine everything and contract indices
std::complex<double> jpacPhoto::rarita_exchange::helicity_amplitude(std::vector<int> helicities, double xs, double xt)
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
      temp *= rarita_propagator(i, j);
      temp *= bottom_vertex(j, lam_vec, lam_targ);

      result += temp;
    }
  }

  return result;
};

//------------------------------------------------------------------------------
// rank-2 traceless tensor
std::complex<double> jpacPhoto::rarita_exchange::g_bar(int mu, int nu)
{
  std::complex<double> result;
  result = exchange_momentum(mu) * exchange_momentum(nu) / mEx2;

  if (mu == nu)
  {
    result -= metric[mu];
  }

  return result;
};

// g_bar contracted with gamma^nu
std::complex<double> jpacPhoto::rarita_exchange::slashed_g_bar(int mu, int i, int j)
{
  std::complex<double> result = 0.;

  for (int nu = 0; nu < 4; nu++)
  {
    std::complex<double> temp;
    temp  = g_bar(mu, nu);
    temp *= metric[nu];
    temp *= gamma_matrices[nu][i][j];

    result += temp;
  }

  return result;
};

//------------------------------------------------------------------------------
// Relative momentum either entering (top vertex) or exiting (bottom vertex) the propagator
std::complex<double> jpacPhoto::rarita_exchange::relative_momentum(int mu, std::string in_out)
{
  std::complex<double> q1_mu, q2_mu;

  if ((in_out == "in") || (in_out == "top") || (in_out == "initial") )
  {
    q1_mu = kinematics->initial->q(mu, s, 0.);
    q2_mu = kinematics->initial->p(mu, s, M_PI);
  }
  else if ((in_out == "out") || (in_out == "bot") || (in_out == "final"))
  {
    q1_mu = kinematics->final->q(mu, s, theta);
    q2_mu = kinematics->final->p(mu, s, theta + M_PI);
  }
  else
  {
    std::cout << "Error! Unkown parameter: " << in_out << "passed to relative_momentum. ";
    std::cout << "Quitting...";
    exit(1);
  }

  return q1_mu - q2_mu;
}

//------------------------------------------------------------------------------
// Rarita-Schwinger Propagator
std::complex<double> jpacPhoto::rarita_exchange::rarita_propagator(int i, int j)
{
  std::complex<double> result = 0.;

  for (int mu = 0; mu < 4; mu++)
  {
    for(int nu = 0; nu < 4; nu++)
    {
      std::complex<double> term_1;
      term_1  = relative_momentum(mu, "in");
      term_1 *= metric[mu];
      term_1 *= g_bar(mu, nu);
      term_1 *= metric[nu];
      term_1 *= relative_momentum(nu, "out");

      std::complex<double> term_2;
      term_2  = relative_momentum(mu, "in");
      term_2 *= metric[mu];
      term_2 *= slashed_g_bar(mu, i, j);
      term_2 *= slashed_g_bar(nu, i, j);
      term_2 *= metric[nu];
      term_2 *= relative_momentum(nu, "out");

      result += -term_1 + term_2 / 3.;
    }
  }

  result *= dirac_propagator(i, j);

  return result;
}
