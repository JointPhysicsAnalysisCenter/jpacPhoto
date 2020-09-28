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
    s = xs; t = xt;
    update_kinematics(); // calculate the other kinematics

    // Form factor
    F_0 = form_factor(t);
    
    // output
    long double result = 1.;
    result  = M_ALPHA * g*g;
    result /= 8. * sqrt(mA2) * mX2 * mX2 * pGam * t*t;
    result /= (2. * sqrt(mA2) * nu - Q2);
    result *= W_00();

    // Amplitude depends on LT
    result *= amplitude_squared();
    
    // Convert from GeV^-2 -> nb
    result /= (2.56819E-6); 

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
  double t_max = kinematics->t_man(s, 1. * deg2rad); // Fall off is extremely fast in t so only integrate over that little bit

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
    double q = sqrt(x * (x - 4. * mA2)) / (2. * sqrt(mA2));

    auto dF = [&] (double r)
    {
        return r * sin(q * r) * charge_distribution(r);
    };

    ROOT::Math::GSLIntegrator ig(   ROOT::Math::IntegrationOneDim::kADAPTIVE,
                                    ROOT::Math::Integration::kGAUSS61);

    ROOT::Math::Functor1D wF(dF);
    ig.SetFunction(wF);
    
    return rho_0 * ig.IntegralUp(0.) / q;

};

// ---------------------------------------------------------------------------
// Amplitude
long double jpacPhoto::primakoff_effect::amplitude_squared()
{
    long double result;

    switch (LT)
    {
        // Longitudinal photon
        case 0: 
        {
            result = pX*pX * Q2 * EX*EX *  sX2;
            break;
        };
        // Transverse photon
        case 1:
        {
            double coshalf2 = (1. + cX) / 2.;
            double sinhalf2 = (1. - cX) / 2.;

            double symC = pX*pGam*(pX + pGam) + EX*nu*(pGam-pX) - 2.*pX*pGam*pGam*cX;
            double symS = pX*pGam*(pX - pGam) + EX*nu*(pGam+pX) - 2.*pX*pGam*pGam*cX;

            double temp = pow(pGam*(nu*(mX2+2.*pX*pX) - 2.*EX*pX*pGam*cX), 2.) / (2. * mX2);

            result  = pow(coshalf2 * symC, 2.);
            result += pow(sinhalf2 * symS, 2.);
            result += temp * sX2;

            break;
        };
        default: result = 0.;
    }

    return result;
};