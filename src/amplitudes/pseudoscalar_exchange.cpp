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

  // Store the invariant energies to avoid having to pass them around 
  s = xs; t = xt, theta = kinematics->theta_s(xs, xt);

  // Because its a scalar exchange we dont have any loose indices to contract
  std::complex<double> result;
  // result  = top_vertex(lam_gam, lam_vec);
  // result *= scalar_propagator();
  // result *= bottom_vertex(lam_targ, lam_rec);
    

  if (lam_vec != lam_gam || lam_targ != lam_rec) 
  {
    return 0.; 
  }
  else
  {
    result  = sqrt(2.) * gNN;
    result *= gGamma / kinematics->mVec;
    result *= sqrt(xr * t) / 2.;
    result *= (kinematics->mVec2 - t);
    result *= scalar_propagator();
  }

  // Multiply by the optional expontial form factor
  if (IF_FF == true)
  {
    double tprime = t - kinematics->t_man(s, 0.);
    result *= exp(b * tprime);
  }

  return result;
};

//------------------------------------------------------------------------------
// Nucleon vertex
std::complex<double> jpacPhoto::pseudoscalar_exchange::bottom_vertex(double lam_targ, double lam_rec)
{
  std::complex<double> result = 0.;
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      // ubar(recoil) * gamma_5 * u(target)
      std::complex<double> temp;
      temp  = kinematics->recoil->adjoint_component(i, lam_rec, s, theta + M_PI); // theta_recoil = theta + pi
      temp *= gamma_5[i][j];
      temp *= kinematics->target->component(j, lam_targ, s, M_PI); // theta_target = pi

      result += temp;
    }
  }

  // Sqrt(2) from isospin considering a charged pion field
  // remove the Sqrt(2) if considering a neutral pion exchange
  result *= sqrt(2.) * gNN;

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
      // (eps*_lam . eps_gam)(q_vec . q_gam)
      std::complex<double> temp1;
      temp1  = kinematics->eps_vec->conjugate_component(mu, lam_vec, s, theta);
      temp1 *= metric[mu];
      temp1 *= kinematics->eps_gamma->component(mu, lam_gam, s, 0.);
      temp1 *= kinematics->initial->q(nu, s, 0.);
      temp1 *= metric[nu];
      temp1 *= kinematics->final->q(nu, s, theta);

      term1 += temp1;

      // (eps*_lam . q_gam)(eps_gam . q_vec)
      std::complex<double> temp2;
      temp2  = kinematics->eps_vec->conjugate_component(mu, lam_vec, s, theta);
      temp2 *= metric[mu];
      temp2 *= kinematics->initial->q(mu, s, 0.);
      temp2 *= kinematics->eps_gamma->component(nu, lam_gam, s, 0.);
      temp2 *= metric[nu];
      temp2 *= kinematics->final->q(nu, s, theta);

      term2 += temp2;
    }
  }

  std::complex<double> result;
  result = term1 - term2;

  // Coupling is normalized to the mass of the Axial vector particle
  result *= gGamma / kinematics->mVec;

  return result;
};

//------------------------------------------------------------------------------
// Simple pole propagator
std::complex<double> jpacPhoto::pseudoscalar_exchange::scalar_propagator()
{
  if (REGGE == false)
  {
    return 1. / (t - mEx2);
  }
  else
  {
    std::complex<double> alpha_t = alpha->eval(t);

    if (std::abs(alpha_t) > 20.) return 0.;

    // Else use the regge propagator
    std::complex<double> result = 1.;
    result  = - alpha->slope();
    result *= 0.5 * (double(alpha->signature) +  exp(-xi * M_PI * alpha_t));
    result *= cgamma(0. - alpha_t);
    result *= pow(s, alpha_t);
    return result;
  }
};
