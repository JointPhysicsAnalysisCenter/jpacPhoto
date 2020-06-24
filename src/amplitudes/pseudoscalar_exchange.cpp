// Charged axial-vector meson photoproduction proceeding through a pseudoscalar (pion) exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:1503.02125 [hep-ph]
// ---------------------------------------------------------------------------

#include "amplitudes/pseudoscalar_exchange.hpp"

//------------------------------------------------------------------------------
// Combine everything and contract indices
std::complex<double> jpacPhoto::pseudoscalar_exchange::helicity_amplitude(std::vector<int> helicities, double xs, double xt)
{
  int lam_gam = helicities[0];
  int lam_targ = helicities[1];
  int lam_vec = helicities[2];
  int lam_rec = helicities[3];

  s = xs; theta = kinematics->theta_s(xs, xt);

  // Because its a scalar exchange we dont have any loose indices to contract
  std::complex<double> result;
  result  = top_vertex(lam_gam, lam_vec);
  result *= scalar_propagator();
  result *= bottom_vertex(lam_rec, lam_targ);

  return result;
};

//------------------------------------------------------------------------------
// Pion form factors
double jpacPhoto::pseudoscalar_exchange::form_factor(double m)
{
  double result;
  result  = (m*m - mPi*mPi);
  result /= (m*m - kinematics->t_man(s, theta));

  return result;
};

//------------------------------------------------------------------------------
// Nucleon vertex
std::complex<double> jpacPhoto::pseudoscalar_exchange::bottom_vertex(double lam_rec, double lam_targ)
{
  std::complex<double> result = 0.;
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      std::complex<double> temp;
      temp  = kinematics->recoil.adjoint_component(i, lam_rec, s, theta + M_PI); // theta_recoil = theta + pi
      temp *= gamma_5[i][j];
      temp *= kinematics->target.component(j, lam_targ, s, M_PI); // theta_target = pi

      result += temp;
    }
  }

  result *= sqrt(2.) * gNN;
  result *= form_factor(LamPi);

  return result;
};

//------------------------------------------------------------------------------
// Photon vertex
std::complex<double> jpacPhoto::pseudoscalar_exchange::top_vertex(double lam_gam, double lam_vec)
{
  std::complex<double> term1 = 0., term2 = 0.;
  for (int mu = 0; mu < 4; mu++)
  {
    for (int nu = 0; nu < 4; nu++)
    {
      // (eps*_lam . eps_gam)(p . q)
      std::complex<double> temp1;
      temp1  = kinematics->eps_vec.conjugate_component(mu, lam_vec, s, theta);
      temp1 *= metric[mu];
      temp1 *= kinematics->eps_gamma.component(mu, lam_gam, s, 0.);
      temp1 *= kinematics->initial.component(nu, "beam", s, 0.);
      temp1 *= metric[nu];
      temp1 *= kinematics->final.component(nu, kinematics->vector_particle, s, theta);

      term1 += temp1;

      // (eps*_lam . p)(eps_gam . q)
      std::complex<double> temp2;
      temp2  = kinematics->eps_vec.conjugate_component(mu, lam_vec, s, theta);
      temp2 *= metric[mu];
      temp2 *= kinematics->initial.component(mu, "beam", s, 0.);
      temp2 *= kinematics->eps_gamma.component(nu, lam_gam, s, 0.);
      temp2 *= metric[nu];
      temp2 *= kinematics->final.component(nu, kinematics->vector_particle, s, theta);

      term2 += temp2;
    }
  }

  std::complex<double> result;
  result = term1 - term2;
  result *= gPsi / kinematics->mVec;
  result *= M_ALPHA / fJpsi;
  result *= form_factor(mJpsi);

  return result;
};

//------------------------------------------------------------------------------
// Simple pole propagator
std::complex<double> jpacPhoto::pseudoscalar_exchange::scalar_propagator()
{
  if (REGGE == false)
  {
    return 1. / (kinematics->t_man(s, theta) - mEx2);
  }

  // Else use the regge propagator
  std::complex<double> alpha_t = alpha->eval(kinematics->t_man(s, theta));

  // the gamma function causes problesm for large t so
  if (std::abs(alpha_t) > 30.)
  {
    return 0.;
  }
  else
  {
    std::complex<double> result;
    result  = alpha->slope();
    result *= 0.5 * (double(alpha->signature) + exp(-xi * M_PI * alpha_t));
    result *= cgamma(-alpha_t);
    result *= pow(xr * s, alpha_t);
    return result;
  }
};
