// Header file with global phyiscal constants.
// Everything is in GeV unless explicitly stated otherwise.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DEBUG_
#define _DEBUG_

#include <iostream>
#include <iomanip>

namespace jpacPhoto
{
    // little function for printing to screen instead of having to copy this line all the time
    template<typename T>
    inline void debug(T x)
    {
        std::cout << x << std::endl;
    };

    template<typename T, typename F>
    inline void debug(T x, F y)
    {
        std::cout << std::left << std::setw(15) << x;
        std::cout << std::left << std::setw(15) << y << std::endl;
    };

    template<typename T, typename F, typename G>
    inline void debug(T x, F y, G z)
    {
        std::cout << std::left << std::setw(15) << x;
        std::cout << std::left << std::setw(15) << y;
        std::cout << std::left << std::setw(15) << z << std::endl;
    };
};

#endif

#ifndef CONSTANT
#define CONSTANT

#include <cmath>
#include <complex>

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    const double deg2rad  = (M_PI / 180.);
    const double EPS      = 1.e-6;
    const double M_ALPHA  = 1. / 137.;
    const double e        = sqrt(4. * M_PI * M_ALPHA);

    const std::complex<double> xr(1., 0.);
    const std::complex<double> xi(0., 1.);
    const std::complex<double> ieps(0., EPS);

    // Meson masses in GeV
    const double mPi        = 0.13957000;
    const double mK         = 0.49367700;
    const double mEta       = 0.54753;
    const double mRho       = 0.77526;
    const double mOmega     = 0.78265;
    const double mPhi       = 1.01956;
    const double mJpsi      = 3.0969160;
    const double mPsi2S     = 3.686;
    const double mD         = 1.86965;
    const double mDstar     = 2.01026;
    const double mUpsilon1S = 9.4603;
    const double mUpsilon2S = 10.02336;
    const double mUpsilon3S = 10.3552;

    // Meson masses squared
    const double mPi2       = mPi * mPi;
    const double mJpsi2     = mJpsi * mJpsi;
    const double mD2        = mD * mD;
    const double mDstar2    = mDstar * mDstar; 

    // Baryon masses
    const double mPro       = 0.938272;
    const double mLambdaC   = 2.28646;

    // Baryon masses squared
    const double mPro2      = mPro * mPro;
    const double mLambdaC2   = mLambdaC * mLambdaC;

    // Decay constants in GeV
    const double fJpsi      = 0.278;
    const double fUpsilon1S = 0.23345;
    const double fUpsilon2S = 0.16563;
    const double fUpsilon3S = 0.1431;

    // Photon lab energy
    inline double E_beam(double W)
    {
        return (W*W / mPro - mPro) / 2.;
    };

    // Center of mass energy given beam energy
    inline double W_cm(double egam)
    {
        return sqrt(mPro * (2. * egam + mPro));
    };

};
// ---------------------------------------------------------------------------

#endif
