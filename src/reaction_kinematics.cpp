// Class to contain all relevant kinematic quantities. The kinematics of the reaction
// gamma p -> V p' is entirely determined by specifying the mass of the vector particle
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "reaction_kinematics.hpp"

// ---------------------------------------------------------------------------
// s-channel variables from 4-vectors
double jpacPhoto::reaction_kinematics::s_man(event fvecs)
{
  return (fvecs.pGam + fvecs.pTarg).M2();
};

double jpacPhoto::reaction_kinematics::z_s(event fvecs)
{
  return z_s(fvecs.s_man(), fvecs.t_man());
};

double jpacPhoto::reaction_kinematics::z_s(double s, double t)
{
  std::complex<double> kq = initial.momentum("beam", s) * final.momentum(vector_particle, s);
  std::complex<double> E1E3 = initial.energy("beam", s) * final.energy(vector_particle, s);

  double result = t - mVec2 + 2.*abs(E1E3);
  result /= 2. * abs(kq);

  return result;
};

// ---------------------------------------------------------------------------
// Invariant variables
double jpacPhoto::reaction_kinematics::t_man(double s, double theta)
{
  std::complex<double> kq = initial.momentum("beam", s) * final.momentum(vector_particle, s);
  std::complex<double> E1E3 = initial.energy("beam", s) * final.energy(vector_particle, s);
  return mVec*mVec - 2. * abs(E1E3) + 2. * abs(kq) * cos(theta);
};

double jpacPhoto::reaction_kinematics::u_man(double s, double theta)
{
  double t = t_man(s, theta);

  return mVec2 + 2. * mPro2 - s - t;
};

// ---------------------------------------------------------------------------
// Scattering angles
double jpacPhoto::reaction_kinematics::theta_s(double s, double t)
{
  double zs = z_s(s, t);

  // Watch for rounding errors near the bounds
  if (std::abs(zs - 1.) < 0.0001)
  {
    return 0.;
  }
  else if (std::abs(zs + 1.) < 0.0001)
  {
    return 1.;
  }
  else
  {
    return acos(z_s(s,t));
  }
};

// Scattering angle in the t-channel
std::complex<double> jpacPhoto::reaction_kinematics::z_t(double s, double t)
{
  std::complex<double> p_t = sqrt(xr * t - 4. * mPro2) / 2.;
  std::complex<double> q_t = sqrt(xr * jpacPhoto::Kallen(t, mVec2, 0.)) / sqrt(xr * 4. * t);

  std::complex<double> result;
  result = 2. * s + t - 2. * mPro2 - mVec2; // s - u
  result /= 4. * p_t * q_t;

  return result;
};

std::complex<double> jpacPhoto::reaction_kinematics::z_u(double s, double theta)
{
  // TODO: fix this
  std::cout << "z_u not implimented yet i fucked up here :p\n";
  std::cout << "if youre seeing this message email me and yell at me because i forgot\n";
  return xr;
};

// ---------------------------------------------------------------------------
// Cosine of crossing an
std::complex<double> jpacPhoto::reaction_kinematics::crossing_angle(std::string particle, double s, double theta)
{
  double t = t_man(s, theta);

  std::complex<double> S_gp, S_vp, T_gv, T_pp;
  S_gp = sqrt(xr * jpacPhoto::Kallen(s, 0., mPro2));
  S_vp = sqrt(xr * jpacPhoto::Kallen(s, mVec2, mPro2));
  T_gv = sqrt(xr * jpacPhoto::Kallen(t, 0., mVec2));
  T_pp = sqrt(xr * jpacPhoto::Kallen(t, mPro2, mPro2));

  std::complex<double> result;
  if (particle == "beam")
  {
    result = - (s - mPro2) * (t - mVec2);
    result /= S_gp * T_gv;
  }
  else if (particle == vector_particle)
  {
    result = -(s + mVec2 - mPro2) * ( t + mVec2) + 2. * pow(mVec2, 2.);
    result /= S_vp * T_gv;
  }
  else if (particle == "target")
  {
    result = - (s + mPro2) * t + 2.*mPro2*mVec2;
    result /= S_gp * T_pp;
  }
  else if (particle == "recoil")
  {
    result = (s + mPro2 - mVec2) * t + 2.*mPro2*mVec2;
    result /= S_vp * T_pp;
  }
  else
  {
    std::cout << "\nreaction_kinematics: particle = " << particle << " not recognized in crossing_angle. Quitting... \n";
    exit(0);
  }

  return result;
};
