// Implementation of unpolarized structure functions at high energies
// Uses the Regge pole parameteriztion from [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] - https://arxiv.org/abs/hep-ph/0402081
// ------------------------------------------------------------------------------

#ifndef DL_F_HPP
#define DL_F_HPP

#include "constants.hpp"
#include "data_set.hpp"
#include "total_xsection.hpp"
#include <Math/Interpolator.h>

namespace jpacPhoto
{
    class DL_F : public raw_total_xsection
    {
        public: 

        DL_F(unsigned x)
        : raw_total_xsection({0, M_PROTON}) // Massless photon
        {
            if (x != 1 && x != 2)
            {
                error("DL_F", "Integer argument must be 1 or 2 for F_1 and F_2 respectively! Defaulting to F_2...");
                return;
            }

            // If F_1 requested, import R values
            if (x == 1)
            {
                auto Rs = import_data<2>("/data/other/R_vals.dat");
                _RBS.SetData(Rs[0], Rs[1]);
            };

            _mode = x;
        };

        inline double evaluate(double s, double q2)
        {
            // Save external variables
            _s = s; _w = sqrt(s); _Q2 = -q2;
            _nu = nu(); _x = x();

            // Calculate F2 by summing regge terms
            double F2 = 0;
            for (int i = 0; i < 3; i++)
            {
                double f_i = _A[i] * pow(_Q2 / (1.+ _Q2/_Q20[i]), 1.+_eps[i]);    
                if (i == 0)   f_i *= pow(1.+_Q2/_Q20[i], _eps[i]/2.);        

                F2 += f_i * pow(_x, -_eps[i]);    
            };

            if (_mode == 2) return F2;

            // Calculate F1 by rescaling
            return (M_PROTON*_nu/_Q2) * F2 / (1+R());
        };


        private:

        int _mode = 2;

        // Saved external variables
        double _s, _w, _Q2;
        // Saved internal variables
        double _nu, _x;

        // R = sigma_L / sigma_T parameters
        double _R0 = 0.23;
        // Instead of calcualting R from Christy & Bolsted,
        // Digitize and interpolate Fig. 15 of [1]
        ROOT::Math::Interpolator _RBS;
        
        inline double R()
        {
            return (_Q2 < 1.5) ? _RBS.Eval(_Q2) : _R0;
        };

        // parameters
        std::array<double,3> _A   = {0.00151, 0.658,   1.01};  // normalizations
        std::array<double,3> _eps = {0.452,   0.0667, -0.476}; // intercepts
        std::array<double,3> _Q20 = {7.85,    0.6,     0.214}; // scale parameters

        // Energy variables
        inline double nu(){ return (_s + _Q2 - M2_PROTON) / (2*M_PROTON); };
        inline double  x(){ return _Q2 / (_s + _Q2 - M2_PROTON); };
    };
};

#endif