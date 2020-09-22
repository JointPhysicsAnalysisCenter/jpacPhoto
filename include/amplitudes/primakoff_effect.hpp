// Axial-vector meson photoproduction proceeding through the primakoff effect
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PRIMAKOFF_
#define _PRIMAKOFF_

#include "amplitude.hpp"
#include "gamma_technology.hpp"

namespace jpacPhoto
{
    class primakoff_effect : public amplitude
    {
        public:
        // Constructor 
        primakoff_effect(reaction_kinematics * xkinem, std::string amp_id = "")
        : amplitude(xkinem, amp_id, 4)
        {};

        void set_params(std::vector<double> params)
        {
            check_Nparams(params); 
            Z = params[0];
            R = params[1];
            a = params[2];
            g = params[3];

            calculate_norm();
        };

        inline void set_Q2(int _LT, double _Q2)
        {
            if (LT > 1 || LT < 0)
            {
                std::cout << "error! invalid parameter in set_Q2(). \n";
                std::cout << "LT = 0 for longitudinal and 1 for transverse photon.\n";
            };

            LT = _LT;
            Q2 = _Q2;
        };

        // individual helicity amplitudes not supported but need to provide definition for virtual class.
        inline std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double t)
        {
            std::cout << "Warning! Individual helicity amplitudes not supported by primakoff_effect!\n";
            return 0.;
        }

        // instead we override the definition of differential_xsection in amplitude.hpp
        double differential_xsection(double s, double t);
        double integrated_xsection(double s);

        // Fermi model nuclear charge distribution
        inline double charge_distribution(double r)
        {
            return 1. / ( 1. + exp((r - R) / a) );
        };

        // Normalized fourier transform of the above charge_distributions 
        double form_factor(double x);
        double F_0;

        private:

        // Parameters
        int    LT    = 0 ;  // longitudinal (0) or transverse (1) photon
        int    Z     = 0 ;  // atomic number
        double R     = 0.;  // radius parameter
        double a     = 0.;  // skin thickness parameter
        double g     = 0.;  // X -> gamma gamma* coupling
        double rho_0 = 0.;  // normalizaton
        double Q2    = 0.;  // Q^2 > 0 of the external photon

        // Calculate the normalization with above parameters
        void calculate_norm();     

        // Kinematic quantities   
        double mX2 = kinematics->mVec2;
        double nu, cX, sX, p;
        inline void update_kinematics()
        {
            // lab frame momentum transfer
            nu = (s - mPro2 + Q2) / (2. * mPro);

            // Momentum of the X
            p  = sqrt(t*t + 4.*mPro*t*nu + 4.*mPro2*(nu*nu - mX2));
            p /= 2. * mPro;

            // Cosine of scattering angle of the X in the lab frame
            cX  = t + Q2 - mX2 + 2.*sqrt(p*p + mX2);
            cX /= 2. * p * sqrt(nu*nu + Q2);

            // Sine of the above 
            sX  = sqrt(1. - cX * cX);
        };

        // Spin averaged amplitude squared
        double amplitude_squared();

        // Top tensor (spin averaged tensor of X -> gamma gamma* interaction)
        std::complex<double> top_tensor(int mu, int nu);

        // Hadronic bottom tensor (gamma* -> N Nbar)
        std::complex<double> bottom_tensor(int mu, int nu);       
    };
};

#endif