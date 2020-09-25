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
    
    // output
    double result;

    // Additional factors
    result  = e / t; // photon progagator
    result *= Z * F_0; // bottom vertex
    result *= 8. * mA2 / (t - 4. * mA2); 
    result *= g / mX2; // from top vertex
    // result *= e;
    result *= result; //squared

    // Amplitude depends on LT
    result *= amplitude_squared();

    // Flux factor
    result /= 4.; 
    result /= 32. * M_PI * (nu*nu + Q2);
    result /= (2.56819E-6); // Convert from GeV^-2 -> nb

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
  double t_max = kinematics->t_man(s, 2. * deg2rad); // Fall off is extremely fast in t so only integrate over that little bit

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
double jpacPhoto::primakoff_effect::amplitude_squared()
{
    double result;

    switch (LT)
    {
        // Longitudinal photon
        case 0: 
        {
            result = pX*pX * Q2 *nu*nu * sX ;
            break;
        };
        // Transverse photon
        case 1:
        {
            result  = mX2 * (3.*pX*pX + 4.*(Q2 + nu*nu) + pX*pX*(cX*cX - sX*sX));
            result += 2. * pX*pX * (Q2 + nu*nu) * sX*sX;
            result -= 8. * mX2 *pX * sqrt(Q2 + nu*nu) * cX;
            result *= Q2*Q2 / (4. *mX2);
            break;
        };
        default: result = 0.;
    }

    return result;
};