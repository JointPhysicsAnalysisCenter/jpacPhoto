// Implementation of a HPWA with spin-3/2 intermediate state. 
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

#ifndef THREE_HALFS_HPP
#define THREE_HALFS_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "utilities.hpp"
#include "partial_wave.hpp"

#include <Eigen/Dense>

namespace jpacPhoto
{
    namespace kmatrix
    {
        class three_halfs : public raw_partial_wave
        {
            public: 

            // Single channel K-matrix (no arguments needed at initialization)
            three_halfs(key k, kinematics xkinem)
            : raw_partial_wave(k, xkinem, 3, "kmatrix::three_halfs")
            {
                // 6 Production parameters and 6 elastic rescattering
                initialize(12);
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

                auto hels   = helicities;
                bool barred = (2*hels[2] - hels[3] < 0);
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
                    // {+1, +1, +1, -1}
                    case 1: i = 2; break;
                    // {+1, -1, +1, +1}
                    case 6: i = 3; break;
                    // {+1, -1,  0, -1}
                    case 9: i = 4; break;
                    // {+1, -1, +1, -1}
                    case 7: i = 5; break;
                    default: return 0;
                };

                return (!barred) ? (_FN(i)+_FU(i))/2 : (_FN(i)-_FU(i))/2; 
            };
            inline complex partial_wave(std::array<int,4> helicities, double s)
            {
                return partial_wave(helicities, sqrt(s)+_ieps);
            };

            protected:

            // Production parameters
            // In the J=3/2 case we decompose the 6x6 problem into two 3x3's
            Eigen::Vector3cd _fp = Eigen::Vector3cd::Zero(), _fm = Eigen::Vector3cd::Zero();
            // Elastic parameters
            std::array<double,6> _a;
            inline void allocate_parameters(std::vector<double> pars)
            {

                for (int i = 0; i < 3; i++) _fp(i) = pars[i];
                for (int i = 0; i < 3; i++) _fm(i) = pars[3+i];
                for (int i = 0; i < 6; i++) _a[i]  = pars[6+i];
                _pars_changed = true;
            };
            inline std::vector<std::string> parameter_labels()
            {
                std::vector<std::string> labels;
                for (int i = 0; i < 6; i++) labels.push_back("f_3/2[" + to_string(i) +"]");
                for (int i = 0; i < 6; i++) labels.push_back("a_3/2[" + to_string(i) +"]");
                return labels;
            };
            
            // Masses 
            double _mV, _mN; 
            // Small imaginary part to help control cut structures
            complex _ieps = I*1E-10;
            // Two body phase-space
            complex _rho = 0;

            // Parameters to know when to recalculate amplitudes
            bool    _pars_changed    = false;
            complex _cached_sqrt_s   = 0.;
            double  _cache_tolerance = 1E-4;

            // Partial wave vectors
            Eigen::VectorXcd _FN  = Eigen::VectorXcd::Zero(6); // Natural 
            Eigen::VectorXcd _FU  = Eigen::VectorXcd::Zero(6); // Unnatural
            // Identity matrix
            Eigen::Matrix3cd _one = Eigen::Matrix3cd::Identity();

            // Recalculate cached partial-waves if needed
            void recalculate(complex sqrt_s)
            {
                // Work in terms of \pm \sqrt{s}
                _cached_sqrt_s = sqrt_s;
                _rho           = rho(sqrt_s*sqrt_s);
                
                Eigen::Matrix3cd one_plus_irhot;
                std::array<Eigen::Matrix3cd,4> m;
                
                // First we calculate the Natural case (\eta_s = +1 or \sqrt{s}>1)
                one_plus_irhot = (_one+I*_rho*t_matrix(+1)); // Unitarized elastic amplitudes
                m = m_matrices(+1);                           // Production kinematic matrices
                _FN.head<3>()  = one_plus_irhot*(m[0]*_fp+m[1]*_fm);
                _FN.tail<3>()  = one_plus_irhot*(/* m[2]*_fp+ */m[3]*_fm);

                // Now we do the exact same but with \eta = -1
                one_plus_irhot = (_one+I*_rho*t_matrix(-1));                   
                m = m_matrices(-1); 
                _FU.head<3>()  = one_plus_irhot*(m[0]*_fp+m[1]*_fm);
                _FU.tail<3>()  = one_plus_irhot*(/* m[2]*_fp+ */m[3]*_fm);

                // Reset flag
                _pars_changed  = false;
            };  

            // Chew-Mandelstam phase-space 
            inline complex rho(complex s)
            {
                double m1 = _mV, m2 = _mN;
                complex rho, xi;
                complex result;

                rho    = csqrt(kallen(s, m1*m1, m2*m2)) / s;
                xi     = 1. - (m1+m2)*(m1+m2)/s;
                result = (rho*log((xi + rho)/(xi - rho))-xi*(m2-m1)/(m2+m1)*log(m2/m1))/PI;
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

            // Calculate 3x3 (little) n-matrices
            Eigen::Matrix3cd n_matrix(unsigned int i, int eta)
            {              
                complex sqrt_s  = (eta == +1) ? + _cached_sqrt_s : - _cached_sqrt_s;
                complex ps      = csqrt(kallen(sqrt_s*sqrt_s, _mV*_mV, _mN*_mN))/2/sqrt_s;
                complex omegafp = omega_f(+1, sqrt_s), omegafm = omega_f(-1, sqrt_s);

                Eigen::Matrix3cd n = Eigen::Matrix3cd::Zero();
                switch (i)
                {
                    case 0: 
                    {
                        n(1,1) = 1;
                        n *= pow(ps, 4) * sqrt_s*sqrt_s;
                        n *= omegafp * omegafp;
                        break;
                    };
                    case 1:
                    {
                        n(0,1) = 1; n(1,0) = 1;
                        n *= pow(ps, 3) * sqrt_s;
                        n *= omegafp * omegafm;
                        break;
                    };
                    case 2:
                    {
                        n(0,1) = 1;       n(1,0) = n(0,1);
                        n(1,2) = sqrt(3); n(2,1) = n(1,2);
                        n *= pow(ps, 2) * sqrt_s; 
                        n *= omegafp * omegafp;
                        break;
                    };
                    case 3:
                    {
                        n(0,0) = 2*_mV*_mV; 
                        n(1,1) = pow(_mN+sqrt_s, 2); 
                        n(0,1) = sqrt(2)*_mV*(_mN+sqrt_s);  n(1,0) = n(0,1);  
                        n *= pow(ps, 2); 
                        n *= omegafm * omegafm;
                        break;
                    };
                    case 4:
                    {
                        complex EV = (sqrt_s*sqrt_s+_mV*_mV-_mN*_mN)/2/sqrt_s;
                        n(0,1) = 2*EV-_mN-sqrt_s;       n(1,0) = n(0,1);
                        n(0,2) = sqrt(6)*_mV;           n(2,0) = n(0,2);
                        n(1,2) = -sqrt(3)*(_mN+sqrt_s); n(2,1) = n(1,2);
                        n *= ps; 
                        n *= omegafp * omegafm;
                        break;
                    };
                    case 5:
                    {
                        complex EV = (sqrt_s*sqrt_s+_mV*_mV-_mN*_mN)/2/sqrt_s;
                        n(0,0) = _mV*_mV;
                        n(1,1) = 2*EV*EV;
                        n(2,2) = 3*_mV*_mV;

                        n(0,1) = sqrt(2)*EV*_mV;   n(1,0) = n(0,1);
                        n(0,2) = sqrt(3)*_mV*_mV;  n(2,0) = n(0,2);
                        n(1,2) = sqrt(6)*EV*_mV;   n(2,1) = n(1,2);
                        n *= omegafp * omegafp;
                        break;
                    };
                    default: return n;
                };
                
                // If unnatural, multiply by (-1)^{\lam-lamp}
                if (eta == -1)
                {
                    n(0,2) *= -1; n(2,0) *= -1;
                    n(1,2) *= -1; n(2,1) *= -1;
                };

                return n;
            };

            // Assemble the 3x3 (little) k-matrix
            Eigen::Matrix3cd k_matrix(int eta)
            {
                Eigen::Matrix3cd k = Eigen::Matrix3cd::Zero();
                for (int n = 0; n < 6; n++) k += n_matrix(n, eta) * _a[n];
                return k;
            };

            // Assemble the 3x3 t-matrix
            // i.e. the unitarized elastic rescattering amplitudes
            Eigen::Matrix3cd t_matrix(int eta)
            {
                Eigen::Matrix3cd k   = k_matrix(eta);
                return (_one-I*_rho*k).inverse()*k;
            };
            
            // Diagonal 3x3 blocks of the 6x6 (big) M-matrix
            std::array<Eigen::Matrix3cd,4> m_matrices(int eta)
            {
                complex sqrt_s = (eta == +1) ? + _cached_sqrt_s : - _cached_sqrt_s;
                complex qs = csqrt(kallen(sqrt_s*sqrt_s, 0.,      _mN*_mN))/2/sqrt_s;
                complex ps = csqrt(kallen(sqrt_s*sqrt_s, _mV*_mV, _mN*_mN))/2/sqrt_s;
                complex EV = (sqrt_s*sqrt_s+_mV*_mV-_mN*_mN)/2/sqrt_s;
                complex x  = pow(_mN-sqrt_s, 2) - _mV*_mV;
                complex y  = (_mN+sqrt_s);

                // m+ (i.e. the upper 3x3 diagonal block)
                Eigen::Matrix3cd mp = Eigen::Matrix3cd::Zero();
                mp(0,1) = sqrt(2)*_mV*qs*x/sqrt_s;
                mp(0,2) = _mV*qs;
                mp(1,0) = ps*ps*qs*sqrt_s;
                mp(1,1) = qs*x*y/sqrt_s;
                mp(1,2) = sqrt(2)*EV*qs;
                mp(2,2) = sqrt(3)*_mV*qs;

                // m- (i.e. the lower 3x3 diagonal block)
                Eigen::Matrix3cd mm = Eigen::Matrix3cd::Zero();
                mm(0,1) = sqrt(6)*_mV*x/sqrt_s;
                mm(0,2) = sqrt(3)*_mV;
                mm(1,0) = sqrt(3)*ps*ps*sqrt_s;
                mm(1,1) = sqrt(3)*x*y/sqrt_s;
                mm(1,2) = sqrt(6)*EV;
                mm(2,2) = 3*_mV;

                // m'+ (i.e. the upper 3x3 off-diagonal block)
                Eigen::Matrix3cd mprimep = Eigen::Matrix3cd::Zero();
                mprimep(0,1) = sqrt(2)*_mN*_mV*x/sqrt_s/sqrt_s;
                mprimep(0,2) = _mN*_mV/sqrt_s;
                mprimep(1,0) = _mN*ps*ps;
                mprimep(1,1) = _mN*x*y/sqrt_s/sqrt_s;
                mprimep(1,2) = sqrt(2)*EV*_mN/sqrt_s;
                mprimep(2,2) = sqrt(3)*_mN*_mV/sqrt_s;

                // m'- (i.e. the lower 3x3 off-diagonal block)
                Eigen::Matrix3cd mprimem = Eigen::Matrix3cd::Zero();

                complex prefactor = omega_i(+1, sqrt_s)*omega_f(+1, sqrt_s)*(_mN-sqrt_s);
                mp *= prefactor; mprimep *= prefactor;
                mm *= prefactor; mprimem *= prefactor;
                
                if (eta == -1) 
                {
                    // MacDowell phase
                    for (int i = 0; i < 3; i++)
                    { 
                        // Multiply rows 2-4 of the 6x6 by -1
                        mp(2,i) *= -1;  mprimep(2,i) *= -1; 
                        mm(0,i) *= -1;  mprimem(0,i) *= -1;
                        mm(1,i) *= -1;  mprimem(1,i) *= -1; 
                    };
                };
                
                return {mp, mprimep, mprimem, mm};
            };
        };
    };
};

#endif