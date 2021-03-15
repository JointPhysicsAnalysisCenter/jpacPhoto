// Vector production via a one-loop box diagram
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "one_loop/box_amplitude.hpp"

// ---------------------------------------------------------------------------
// Evaluate the entire loop diagram in terms of external helicites and invariant mass / momentum transfer of the gamma p -> jpsi p process
std::complex<double> jpacPhoto::box_amplitude::helicity_amplitude(std::array<int,4> helicities, double s, double t)
{
    // Store the invariant energies to avoid having to pass them around 
    _s = s; _t = t, _theta = _kinematics->theta_s(s, t);

    // Pass external values to the discontinuity
    _disc->set_externals(helicities, _theta);

    // Compute both parts of the integral
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kDEFAULT);
    
    auto F = [&](double sp)
    {
        double result = _disc->eval(sp);
        return result;
    };

    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double real = ig.IntegralCauchy(_s_thr + EPS, _s_cut, _s);
    double imag = - M_PI * _disc->eval(_s);

    std::complex<double> result =  (real + XI * imag) / M_PI;

    return result;
};