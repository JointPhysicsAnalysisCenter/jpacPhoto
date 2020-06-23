// Vector meson photoproduction dynamics proceeding through a pomeron exchange
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/pomeron_exchange.hpp"

// ---------------------------------------------------------------------------
// Given a set of helicities for each particle, assemble the helicity amplitude by contracting Lorentz indicies
std::complex<double> jpacPhoto::pomeron_exchange::helicity_amplitude(std::vector<int> helicities, double s, double t)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  double theta = kinematics->theta_s(s,t);

  std::complex<double> result = 0.;
  for (int mu = 0; mu < 4; mu++)
  {
    std::complex<double> temp = regge_factor(s, t);
    temp *= top_vertex(mu, lam_gam, lam_vec, s, theta);
    temp *= metric[mu];
    temp *= bottom_vertex(mu, lam_targ, lam_rec, s, theta);

    result += temp;
  }

  return result;
};

// ---------------------------------------------------------------------------
// Bottom vertex coupling the target and recoil proton spinors to the vector pomeron
std::complex<double> jpacPhoto::pomeron_exchange::bottom_vertex(int mu, int lam_targ, int lam_rec, double s, double theta)
{
  std::complex<double> result = 0.;
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
    std::complex<double> temp;
    // Recoil oriented an angle theta + pi
    temp = kinematics->recoil.adjoint_component(i, lam_rec, s, theta + M_PI);

    // vector coupling
    temp *= gamma_matrices[mu][i][j];

    // target oriented in negative z direction
    temp *= kinematics->target.component(j, lam_targ, s, M_PI);

    result += temp;
    }
  }

  return result;
};

// ---------------------------------------------------------------------------
// Top vertex coupling the photon, pomeron, and vector meson.
std::complex<double> jpacPhoto::pomeron_exchange::top_vertex(int mu, int lam_gam, int lam_vec, double s, double theta)
{
  std::complex<double> sum1 = 0., sum2 = 0.;
  for (int nu = 0; nu < 4; nu++)
  {
    std::complex<double> temp1, temp2;

    temp1 = kinematics->initial.component(nu, "beam", s, 0.);
    temp1 *= metric[nu];
    temp1 *= kinematics->eps_vec.conjugate_component(nu, lam_vec, s, theta);
    sum1 += kinematics->eps_gamma.component(mu, lam_gam, s, 0.) * temp1;

    temp2 = kinematics->eps_gamma.component(nu, lam_gam, s, 0.);
    temp2 *= metric[nu];
    temp2 *= kinematics->eps_vec.conjugate_component(nu, lam_vec, s, theta);
    sum2 += kinematics->initial.component(mu, "beam", s, 0.) * temp2;
  }

  return -sum1 + sum2;
};

// ---------------------------------------------------------------------------
// Usual Regge power law behavior, s^alpha(t) with an exponential fall from the forward direction
std::complex<double> jpacPhoto::pomeron_exchange::regge_factor(double s, double t)
{
  if (s < kinematics->sth)
  {
    std::cout << " \n pomeron_exchange: Trying to evaluate below threshold (sqrt(s) = " << sqrt(s) << ")! Quitting... \n";
    exit(0);
  }

  if (s - kinematics->sth < 0.001)
  {
    return 0.;
  }

  double t_min = kinematics->t_man(s, 0.); // t_min = t(theta = 0)

  std::complex<double> result = exp(b0 * (t - t_min)) / s;
  result *= pow(s - kinematics->sth, pomeron_traj->eval(t));
  result *= xi * norm * sqrt(4. * M_PI * M_ALPHA);

  return result;
};
