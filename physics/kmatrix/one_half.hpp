// Implementation of a HPWA with spin-1/2 intermediate state
//
// This amplitude explicitly uses Eigen C++ to manipulate complex matrices,
// make sure an environment variable EIGEN points to the top level directory
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Instituto de Ciencias Nuclears (ICN),
//               Universidad Nacional Autonoma de México (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#ifndef ONE_HALF_HPP
#define ONE_HALF_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "utilities.hpp"
#include "partial_wave.hpp"

#include <Eigen/Dense>

namespace jpacPhoto
{
    namespace kmatrix
    {
        class one_half : public raw_partial_wave
        {
            public: 

            // Single channel K-matrix (no arguments needed at initialization)
            one_half(key k, kinematics xkinem)
            : raw_partial_wave(k, xkinem, 1, "kmatrix::one_half")
            {
                // 2 Production parameters and 3 elastic rescattering
                initialize(5);
                _mV = xkinem->get_meson_mass();
                _mN = xkinem->get_recoil_mass();
            };

            // -----------------------------------------------------------------------
            // Virtuals 

            // Specify vector quantum numbers
            inline std::vector<quantum_numbers> allowed_mesons() { return {VECTOR};   };
            // and proton target
            inline std::vector<quantum_numbers> allowed_baryons(){ return {HALFPLUS}; };
            // We have explicit helicity dependence in s-channel
            inline helicity_frame native_helicity_frame(){ return S_CHANNEL; };

            // Partial wave comes from K-matrix unitarized form
            inline complex partial_wave(std::array<int,4> helicities, complex sqrt_s)
            {
                bool need_to_recalculate =  !are_equal(sqrt_s, _cached_sqrt_s, _cache_tolerance) || _pars_changed;
                if ( need_to_recalculate ) recalculate(sqrt_s);

                auto hels = helicities;
                bool barred = (2*helicities[2] - helicities[3] < 0);
                if (barred) { hels[2] *= -1; hels[3] *= -1; };

                // Many years ago, I didnt think very hard about the order the helicities should be arranged
                // Because of that we need to remap them to match our matrices....
                // is there a more elegant way to do this? probably but idc
                int h = _kinematics->helicity_index(hels), i;
                switch (h)
                {   
                    // {+1, +1, +1, +1}
                    case 0: i = 0; break;
                    // {+1, +1,  0, -1}
                    case 3: i = 1; break;
                    default: return 0;
                };

                return (!barred) ? (_Fp(i)+_Fm(i))/2 : (_Fp(i)-_Fm(i))/2; 
            };

            // Evaluate on the real line (needed for raw_partial_wave::partial_wave)
            inline complex partial_wave(std::array<int,4> helicities, double s)
            {
                return partial_wave(helicities, sqrt(s)+_ieps);
            };
            
            protected:

            // Save the different parameters
            Eigen::Vector2cd _Q = Eigen::Vector2cd::Zero();
            std::array<double,3> _elastic_pars;
            inline void allocate_parameters(std::vector<double> pars)
            {
                for (int i = 0; i < pars.size(); i++)
                {
                    if (i < 2) _Q(i)               = pars[i]; 
                    else       _elastic_pars[i-2]  = pars[i];
                };
                _pars_changed = true;
            };
            inline std::vector<std::string> parameter_labels()
            {
                return {"f_1/2[0]", "f_1/2[1]", 
                        "a_1/2[0]", "a_1/2[1]", "a_1/2[2]"};
            };
            
            // Masses 
            double _mV, _mN; 
            // Small imaginary part to help control cut structures
            complex _ieps = I*1E-10;
            // Two body phase-space
            complex _rho  = 0;

            // Parameters to know when to recalculate amplitudes
            bool    _pars_changed    = false;
            complex _cached_sqrt_s   = 0.;
            double  _cache_tolerance = 1E-4;

            // Partial wave vectors
            Eigen::Vector2cd _Fp   = Eigen::Vector2cd::Zero(), _Fm = Eigen::Vector2cd::Zero();
            // Identity matrix
            Eigen::Matrix2cd _one  = Eigen::Matrix2cd::Identity();

            // Recalculate cached partial-waves if needed
            void recalculate(complex sqrt_s)
            {
                _cached_sqrt_s = sqrt_s;
                _rho          = rho(sqrt_s*sqrt_s);
                // Unitarized with the Q-vector of production parameters
                _Fp = (_one+I*_rho*T_matrix(+1))*M_matrix(+1)*_Q;
                // Now we do the exact same but with \eta = -1
                _Fm = (_one+I*_rho*T_matrix(-1))*M_matrix(-1)*_Q;
                // Once recalculated reset the pars changed flag
                _pars_changed = false;
            };  

            // Chew-Mandelstam phase-space 
            inline complex rho(complex s)
            {
                double m1 = _mV, m2 = _mN;
                complex rho, xi;
                complex result;

                rho    = csqrt(kallen(s, m1*m1, m2*m2)) / s;
                xi     = 1. - (m1+m2)*(m1+m2)/s;
                result = (rho*log((xi + rho) / (xi - rho)) - xi*(m2-m1)/(m2+m1)*log(m2/m1)) / PI;
                return I*result / (16.*PI);
            };

            // Initial fermion factor as function of sqrt_s
            inline complex omega_i(int pm, complex sqrt_s)
            { 
                complex Ei = (sqrt_s*sqrt_s + _mN*_mN)/(2*sqrt_s);
                return csqrt(Ei+pm*_mN);
            };
            // Recoil fermion factor as function of sqrt_s
            inline complex omega_f(int pm, complex sqrt_s)
            { 
                complex Ef = (sqrt_s*sqrt_s + _mN*_mN - _mV*_mV)/(2*sqrt_s);
                return csqrt(Ef+pm*_mN);
            };

            // Calculate 2x2 N-matrices
            inline Eigen::Matrix2cd N_matrix(unsigned int n, int eta)
            {              
                complex sqrt_s = (eta == +1) ? + _cached_sqrt_s : - _cached_sqrt_s;
                Eigen::Matrix2cd N = Eigen::Matrix2cd::Zero();

                // Factors of (2 p_s \sqrt{s}) 
                complex ps = csqrt(kallen(sqrt_s*sqrt_s,_mV*_mV, _mN*_mN))/2/sqrt_s;

                switch (n)
                {
                    case 0: 
                    {
                        N(1,1) = 1;
                        N *= pow(ps*sqrt_s, 2);
                        N *= omega_f(+1, sqrt_s) * omega_f(+1, sqrt_s);
                        return N;
                    };
                    case 1:
                    {
                        N(0,1) = 1; N(1,0) = 1;
                        N *= ps*sqrt_s;
                        N *= omega_f(+1, sqrt_s) * omega_f(-1, sqrt_s);
                        return N;
                    };
                    case 2:
                    {
                        complex ep = _mN + sqrt_s;
                        N(0,0) = 2*_mV*_mV;       N(0,1) = sqrt(2)*_mV*ep;
                        N(1,0) = sqrt(2)*_mV*ep;  N(1,1) = ep*ep;   
                        N *= omega_f(-1, sqrt_s) * omega_f(-1, sqrt_s);
                        return N;
                    };
                    default: return N;
                    // No eta dependence because macdowell phase is trivial
                };
            };

            // 2x2 K-matrix is sum of N's multiplied by the scattering lengths
            inline Eigen::Matrix2cd K_matrix(int eta)
            {
                Eigen::Matrix2cd K = Eigen::Matrix2cd::Zero();
                for (int n = 0; n < 3; n++) K += N_matrix(n, eta) * _elastic_pars[n];
                return K;
            };

            // 
            inline Eigen::Matrix2cd T_matrix(int eta)
            {
                Eigen::Matrix2cd K = K_matrix(eta);
                return (_one-I*_rho*K).inverse()*K;
            };
            
            // Mixing matrix between helicities and the 2 possible helicity structures
            inline Eigen::Matrix2cd M_matrix(int eta)
            {
                complex sqrt_s = (eta == +1) ? + _cached_sqrt_s : - _cached_sqrt_s;

                Eigen::Matrix2cd M = Eigen::Matrix2cd::Zero();
                complex ep = _mN + sqrt_s, em = _mN - sqrt_s;

                M(0,0) = 0;                M(0,1) = 2*sqrt(2)*_mV;
                M(1,0) = ep*ep - _mV*_mV;  M(1,1) = 2*ep;
                M *= omega_i(+1, sqrt_s) * omega_f(-1, sqrt_s) * em;

                // No eta dependence because macdowell phase is trivial
                return M;
            };
        };
    };
};

#endif