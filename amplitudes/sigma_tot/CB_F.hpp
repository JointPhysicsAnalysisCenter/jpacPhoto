// Implementation of unpolarized structure functions in the resonance region
// Uses the BW + background parameteriztion from [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] - https://arxiv.org/abs/0712.3731
// ------------------------------------------------------------------------------

#ifndef CB_F_HPP
#define CB_F_HPP

#include "constants.hpp"
#include "data_set.hpp"
#include "total_xsection.hpp"
#include <Math/Interpolator.h>

namespace jpacPhoto
{
    class CB_F : public raw_total_xsection
    {
        public:

        CB_F(unsigned x)
        {
            if (x != 1 && x != 2)
            {
                error("CB_F", "Integer argument must be 1 or 2 for F_1 and F_2 respectively! Defaulting to F_1...");
                return;
            }

            _mode = x;
        };

        inline double evaluate(double s, double q2)
        {
            // save external variables
            _s = s; _w = sqrt(s); _Q2 = -q2;
            _x = x();

            double sigmaT = sigma_T();
            
            if (_mode == 1) return prefactors() * sigmaT;

            double sigmaL = sigma_L();
            
            return prefactors() * (2*_x)/(1 + 4*M2_PROTON*_x*_x/_Q2) * (sigmaT + sigmaL);
        };

        // ----------------------------------------------------------------------
        // POLARIZED CROSS SECTIONS
        inline double sigma_T(double W, double Q2){ return sigmaT_NR(W, Q2) + sigmaT_R(W, Q2); };
        inline double sigma_T(){ return sigma_T(_w, _Q2); };

        inline double sigma_L(double W, double Q2){ return sigmaL_NR(W, Q2) + sigmaL_R(W, Q2); };
        inline double sigma_L(){ return sigma_L(_w, _Q2); };
    
        // ----------------------------------------------------------------------
        // NONRESONANT BACKGROUNDS 

        // TRANSVERSE contribution to nonresonant background
        inline double sigmaT_NR(double W, double Q2)
        {
            // Check threshold
            if (W <= M_PROTON + M_PION) return 0;

            // Parameters
            std::array<double,2> sigma0 = {246.1,   -89.4};
            std::array<double,2> a      = {0.0675,  0.2098};
            std::array<double,2> b      = {1.3501,  1.5715};
            std::array<double,2> c      = {0.1205,  0.0907};
            std::array<double,2> d      = {-0.0038, 0.0104};
            double Q20 = 0.05;

            // xp depends on Q20
            double xp = 1/(1+(W*W-(M_PROTON+M_PION)*(M_PROTON+M_PION))/(Q2 + Q20));

            double sigma = 0;
            for (int i = 0; i < 2; i++)
            {
                sigma += sigma0[i]*xp*pow(W-M_PROTON-M_PION,i+1.5)/pow(Q2 + a[i], b[i]+c[i]*Q2+d[i]*Q2*Q2);
            };

            return sigma;
        };
        inline double sigmaT_NR(){return sigmaT_NR(_w, _Q2);};


        // LONGITUDINAL PHOTONS

        inline double sigmaL_NR(double W, double Q2)
        {
            // Check threshold
            if (W <= M_PROTON + M_PION) return 0;

            // Parameters
            double sigma0 = 86.7;
            double a = 0.0, b = 4.0294, c = 3.1285, d = 0.3340, e = 4.9623;
            double Q20 = 0.125, m0 = 4.2802, mu = 0.33;

            // Intermediate variables
            double xp   = 1/(1 + (W*W - (M_PROTON+M_PION)*(M_PROTON+M_PION))/(Q2 + Q20));
            double tau  = log(log((Q2 + m0)/mu/mu ) / log(m0/mu/mu));

            double sigma;
            sigma  = sigma0 * pow(xp, d+e*tau);
            sigma *= pow(1-xp, a*tau+b)/(1-x());
            sigma *= pow(Q2, c) / pow(Q2 + Q20, 1+c);
            
            return sigma;
        };
        inline double sigmaL_NR(){return sigmaT_NR(_w, _Q2);};

        // ----------------------------------------------------------------------
        // RESONANCES

        inline double sigmaT_R(double W, double Q2)
        {
            double result = 0;
            for (auto res : _resonances)
            {
                result += W * res.BW(W)*res.A_T(Q2)*res.A_T(Q2);
            }
            return result;
        };

        inline double sigmaL_R(double W, double Q2)
        {
            double result = 0;
            for (auto res : _resonances)
            {
                result += W * res.BW(W)*res.A_L(Q2)*res.A_L(Q2);
            }
            return result;
        };

        // ----------------------------------------------------------------------
    
        private:

        // Whether outputing F_1 or F_2
        int _mode = 1;

        // saved external variables
        double _s, _w, _Q2;
        // Saved internal variables
        double _nu, _x;
        
        // Energy variables
        inline double nu(){ return (_s + _Q2 - M2_PROTON) / (2*M_PROTON); };
        inline double  x(){ return _Q2 / (_s + _Q2 - M2_PROTON); };
        
        // Kinematic prefactors connecting structure functions to sigma's
        inline double prefactors()
        {
            double K = (_s - M2_PROTON)/(2*M_PROTON);
            return M_PROTON * K / (4*PI*PI*ALPHA) / 389.39;
        };
        
       // Each resonance pre-calculates and stores its own quantities
        class resonance
        {
            public:
            
            resonance(int l, std::array<double,6> bw_pars, std::array<double,7> a_pars)
            : _ell(l), _m(bw_pars[0]), _gamma(bw_pars[1]),
              _beta({bw_pars[2], bw_pars[3], bw_pars[4]}), _x(bw_pars[5]),
              _AT0(a_pars[0]), _a(a_pars[1]), _b(a_pars[2]), _c(a_pars[3]),
              _AL0(a_pars[4]), _d(a_pars[5]), _e(a_pars[6])
            {
                // Evaluate and store all quantities that do not depend on W or Q2

                // Photon momenta at resonance mass
                _KR    = (_m*_m - M2_PROTON) / (2*M_PROTON);
                _KhatR = (_m*_m - M2_PROTON) / (2*_m);

                // Decay meson CM momenta at resonance mass
                for (int i = 0; i < 3; i++)
                {
                    if (_m <= M_PROTON + sqrt(_decay_m2[i])) continue;
                    _pR[i] = sqrt(Kallen(_m*_m, M2_PROTON, _decay_m2[i]))/(2*_m);
                };
            };

            // Return the amplitude squared from a BW amplitude
            inline double BW(double W)
            {
                check_W_cache(W);
                double bw = _gamma_tot * _gamma_photon / _gamma / (pow(W*W-_m*_m, 2) + _m*_m*_gamma_tot*_gamma_tot);   
                return _KR*_KhatR/(_K*_Khat) * bw;
            };

            // Transverse photocoupling
            inline double A_T(double Q2)
            {
                check_Q2_cache(Q2); 
                return _AT;
            };

            // Longitudinal photocoupling
            inline double A_L(double Q2)
            {
                check_Q2_cache(Q2); 
                return _AL;
            };
            
            private:

            // Use a cache system to only calculate things once at each W and Q2 step
            inline void check_W_cache(double W)
            {
                // If below tolerance do nothing
                if ( std::abs(_cachedW - W) < 1E-5 ) return;
                 
                // Else recalculate quantites
                _K    = (W*W - M2_PROTON) / (2*M_PROTON); // lab frame
                _Khat = (W*W - M2_PROTON) / (2*W);        // CM  frame

                // Virtual photon width
                _gamma_photon = _gamma*pow(_Khat/_KhatR,2)*(_KhatR*_KhatR+_x*_x)/(_Khat*_Khat+_x*_x);

                // Calculate Gamma total
                // 0 - pi 
                // 1 - pi pi p
                // 2 - eta p
                _gamma_tot = 0;
                for (int i = 0; i < 3; i++)
                {
                    // Check we are above threshold
                    if (W < sqrt(_decay_m2[i]) + M_PROTON) continue;
                    if ( is_zero(_beta[i]) ) continue;

                    // Meson momentum in CM frame
                    double p = sqrt(Kallen(W*W, _decay_m2[i], M2_PROTON))/(2*W);

                    // Energy dependent decay fraction
                    double betahat = _beta[i]*pow(p/_pR[i], 2*_ell+1+3*(i==1))*pow((_pR[i]*_pR[i]+_x*_x)/(p*p+_x*_x), _ell+2*(i==1));
                    if (i == 1) betahat *= (W/_m);
                    
                    _gamma_tot += betahat*_gamma;
                };
            };

            inline void check_Q2_cache(double Q2)
            {
                // If below tolerance do nothing
                if ( std::abs(_cachedQ2 - Q2) < 1E-5 ) return;

                // Else recalcualte A_T and A_L
                _AT = _AT0/pow(1+Q2/0.91, _c)*(1+ _a*Q2/(1+_b*Q2));
                _AL = _AL0*Q2/(1+_d*Q2)*exp(-_e*Q2);
            };

            // Cached values
            double _cachedW, _cachedQ2;

            // Energy dependent widths
            double _gamma_tot, _gamma_photon;

            // Photon momenta 
            double _K, _Khat;   // Evaluated at W
            double _KR, _KhatR; // Evaluated at W = m

            // Photocouplings
            double _AT, _AL;

            // BELOW ARE FIXED PARAMETERS

            // Intrinsic properties
            int _ell;
            double _m, _gamma;

            // Related to decay mesons 
            std::array<double,3> _beta, _pR, _decay_m2 = {M2_PION, 4*M2_PION, M_ETA*M_ETA};

            // Dampening factor 
            double _x;

            // A_T parameters
            double _AT0, _a, _b, _c;

            // A_L parameters
            double _AL0, _d, _e;
        };

        // BW parameters for all resonances
        std::array<double,6> P33_BW = {1.230, 0.136, 1.00, 0.00, 0.00, 0.1446};
        std::array<double,6> S11_BW = {1.530, 0.220, 0.45, 0.10, 0.45, 0.215};
        std::array<double,6> D13_BW = {1.506, 0.083, 0.65, 0.35, 0.00, 0.215};
        std::array<double,6> F15_BW = {1.698, 0.096, 0.65, 0.35, 0.00, 0.215};
        std::array<double,6> S15_BW = {1.665, 0.109, 0.40, 0.50, 0.10, 0.215};
        std::array<double,6> P11_BW = {1.433, 0.379, 0.65, 0.35, 0.00, 0.215};
        std::array<double,6> FXX_BW = {1.934, 0.380, 0.50, 0.50, 0.00, 0.215};

        // Photocoupling parameters for all resonances
        std::array<double,7> P33_A  = {7.780,  4.229,  1.260,   2.124, 29.4140, 19.910, 0.226};
        std::array<double,7> S11_A  = {6.335,  6823.2, 33521.0, 2.669, 0.0,     0.0,    0.0};
        std::array<double,7> D13_A  = {0.603,  21.240, 0.056,   2.489, 157.92,  97.046, 0.310};
        std::array<double,7> F15_A  = {2.330, -0.288,  0.186,   0.064, 4.216,   0.038,  1.218};
        std::array<double,7> S15_A  = {1.979, -0.562,  0.390,   0.549, 13.764,  0.314,  3.0};
        std::array<double,7> P11_A  = {0.0225, 462.13, 0.192,   1.914, 5.5124,  0.054,  1.309};
        std::array<double,7> FXX_A  = {3.419,  0.0,    0.0,     1.0,   11.0,    1.895,  0.514};

        // Container of all individual resonance contributions
        std::vector<resonance> _resonances = {
            resonance(1, P33_BW, P33_A ),
            resonance(0, S11_BW, S11_A ),
            resonance(2, D13_BW, D13_A ),
            resonance(3, F15_BW, F15_A ),
            resonance(0, S15_BW, S15_A ),
            resonance(1, P11_BW, P11_A ),
            resonance(3, FXX_BW, FXX_A )
        };
    };
};

#endif 