// Vector meson photoproduction dynamics proceeding through a pomeron exchange
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/pomeron_exchange.hpp"

// ---------------------------------------------------------------------------
// Given a set of helicities for each particle, assemble the helicity amplitude by contracting Lorentz indicies
std::complex<double> jpacPhoto::pomeron_exchange::helicity_amplitude(std::vector<int> helicities, double xs, double xt)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  // Save energies 
  s = xs; t = xt; theta = kinematics->theta_s(xs, xt);

  std::complex<double> result = 0.;

  // IF using helicity conserving delta fuction model
  if (DELTA == true)
  {
    (lam_gam == lam_vec && lam_rec == lam_targ) ? (result = regge_factor()) : (result = 0.);
    return result;
  }

  // else use Lesniak-Szczepaniak Model
  for (int mu = 0; mu < 4; mu++)
  {
    std::complex<double> temp = regge_factor();
    temp *= top_vertex(mu, lam_gam, lam_vec);
    temp *= metric[mu];
    temp *= bottom_vertex(mu, lam_targ, lam_rec);

    result += temp;
  }

  return result;
};

// ---------------------------------------------------------------------------
// Bottom vertex coupling the target and recoil proton spinors to the vector pomeron
std::complex<double> jpacPhoto::pomeron_exchange::bottom_vertex(int mu, int lam_targ, int lam_rec)
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

  // Divide by s to remove the energy dependence left over by the spinors
  result /= s;

  return result;
};

// ---------------------------------------------------------------------------
// Top vertex coupling the photon, pomeron, and vector meson.
std::complex<double> jpacPhoto::pomeron_exchange::top_vertex(int mu, int lam_gam, int lam_vec)
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
std::complex<double> jpacPhoto::pomeron_exchange::regge_factor()
{
  if (s < kinematics->sth)
  {
    std::cout << " \n pomeron_exchange: Trying to evaluate below threshold (sqrt(s) = " << sqrt(s) << ")! Quitting... \n";
    exit(0);
  }

  double t_min = kinematics->t_man(s, 0.); // t_min = t(theta = 0)

  std::complex<double> result = exp(b0 * (t - t_min));
  
  // Use physical threshold
  if (DELTA == false)
  {
    result *= pow(s - kinematics->sth, pomeron_traj->eval(t));
  }
  // Unless using the old, helicity conserving model, in which case use fitted value.
  else
  {
    result *= pow( xr * (s - 16.8), pomeron_traj->eval(t));
  }
  
  result *= xi * norm * e;

  return result;
};
