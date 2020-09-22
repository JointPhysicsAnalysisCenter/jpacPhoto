// Axial-vector meson photoproduction proceeding through the primakoff effect
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/primakoff_effect.hpp"

// ---------------------------------------------------------------------------
// Differential cross-sections with all the flux factors
double jpacPhoto::primakoff_effect::differential_xsection(double xs, double xt)
{
    // update saved energies
    s = xs; t = xt; theta = kinematics->theta_s(xs, xt);
    update_kinematics(); // calculate the other kinematics

    // Form factor
    F_0 = form_factor(t);
    
    // Amplitude depends on LT
    // double result = 1.;
    double result = amplitude_squared();

    // additional factors from nuclear tensor
    result *= 16. * M_PI;
    result *= Z * Z * 16. * mPro2 * mPro2 / (4. * M_PI);
    result *= F_0 * F_0;
    result /= (t - 4. * mPro2) * (t - 4. * mPro2);

    // from photonic tensor
    result *= (g / mX2) * (g / mX2);

    // photon propagator (squared)
    result *= e * e / (t * t);

    // Normalization for dxs
    // result /= 64. * M_PI * s;
    // result /= real(pow(kinematics->initial->momentum(s), 2.));
    // result /= (2.56819E-6); // Convert from GeV^-2 -> nb

    return result;
};

// ---------------------------------------------------------------------------
// Inegrated total cross-section
// IN NANOBARN
double jpacPhoto::primakoff_effect::integrated_xsection(double s)
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
// Normalization
void jpacPhoto::primakoff_effect::calculate_norm()
{
    auto F = [&](double r)
    {
        return r * r * charge_distribution(r);
    };

    ROOT::Math::GSLIntegrator ig(   ROOT::Math::IntegrationOneDim::kADAPTIVE,
                                    ROOT::Math::Integration::kGAUSS61);

    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    // Integrate [0:inf]
    rho_0 = 1. / ig.IntegralUp(0.);
};

// ---------------------------------------------------------------------------
// Fourier transform the charge_distribution
double jpacPhoto::primakoff_effect::form_factor(double x)
{
    // momentum in the t channel
    double q = sqrt(x * (x - 4. * mPro2)) / (2. * mPro);

    auto dF = [&] (double r)
    {
        return r * r * sin(q * r) * charge_distribution(r);
    };

    ROOT::Math::GSLIntegrator ig(   ROOT::Math::IntegrationOneDim::kADAPTIVE,
                                    ROOT::Math::Integration::kGAUSS61);

    ROOT::Math::Functor1D wF(dF);
    ig.SetFunction(wF);
    
    return rho_0 * ig.IntegralUp(0.) / q;

};

// ---------------------------------------------------------------------------
// Amplitude
double jpacPhoto::primakoff_effect::amplitude_squared()
{
    double result;

    switch (LT)
    {
        // Longitudinal photon
        case 0: 
        {
            result = p*p * Q2 *nu*nu * sqrt(1. - cX*cX);
            break;
        };
        // Transverse photon
        case 1:
        {
            result  = mX2 * (3.*p*p + 4.*(Q2 + nu*nu) + p*p*(cX*cX - sX*sX));
            result += 2. * p*p * (Q2 + nu*nu) * sX*sX;
            result -= 8. * mX2 *p * sqrt(Q2 + nu*nu) * cX;
            result *= Q2*Q2 / (4. *mX2);
            break;
        };
        default: result = 0.;
    }

    return result;
};