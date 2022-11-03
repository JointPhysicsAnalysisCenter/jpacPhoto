// Amplitude assuming a simple effective-range expansion near-threshold
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef SCATTERING_LENGTH
#define SCATTERING_LENGTH

#include "amplitude.hpp"

namespace jpacPhoto
{
    class scattering_length : public amplitude
    {
        // ---------------------------------------------------------------------------
        
        public:

        scattering_length(reaction_kinematics * xkinem, int lmax, std::string id = "scattering_length")
        : amplitude(xkinem, "scattering_length", id), _lmax(lmax)
        {
            set_nParams(3);
            check_JP(xkinem);
        };

        // Set the free parameters
        inline void set_params(std::vector<double> params)
        {
            check_nParams(params);

            _N  = params[0];
            _a  = params[1];
            _s0 = params[2];

            _extra_couplings.clear();
            for (int n_th = 0; n_th < _extra_thresholds.size(); n_th++)
            {
                _extra_couplings.push_back(params[3+n_th]);
            };
        };

        inline void add_threshold(double m1, double m2)
        {
            // Save the masses of the new threshold
            _extra_thresholds.push_back({m1, m2});

            // and add one to the number of free parameters
            set_nParams( get_nParams() + 1 );
        };  

        // Assemble the helicity amplitude by summing with appropriate angular funxtions
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double xs, double xt);

        // Expansion is directly in the S channel so helicities should always be defined in this frame
        helicity_channel helicity_CM_frame(){ return helicity_channel::S; };

        // only axial-vector, vector, and pseudo-scalar available
        inline std::vector<std::array<int,2>> allowed_meson_JP(){  return { {1, -1} }; };
        inline std::vector<std::array<int,2>> allowed_baryon_JP(){ return { {1, +1} }; };

        // ---------------------------------------------------------------------------

        protected:

        // Number of partial-wave terms to consider 
        int _lmax = 1;

        // Free parameters
        double _N;         // Normalization
        double _a;         // Scattering length
        double _s0;        // Scale parameter  

        // A vector containing all open thresholds
        std::vector<std::array<double,2>> _extra_thresholds; 
        std::vector<double>               _extra_couplings;

        // Legendre functions
        double P_l(int l, double z);
        inline std::complex<double> q2() { return Kallen(_s * XR, _mX*_mX * XR, _mR*_mR * XR) / (4.*_s); };

        // Partial wave amplitude
        std::complex<double> f_l(int l)
        {
            std::complex<double> K = -1./_a;
            double norm = _N ;

            return norm * K * barrier_factor(l) / (1. - XI * K * (rho(l) + rho_inelastic(l)) );
        };
        
        // Redefine momenta here instead of using the reaction_kinematics versions
        // so that they can be appropriately analytically continued below threshold
        inline std::complex<double> barrier_factor(int l, double m1, double m2, double szero)
        {         
            std::complex<double> pq = sqrt( Kallen(_s * XR, _mB*_mB * XR, _mT*_mT * XR)*Kallen(_s * XR, m1*m1 * XR, m2*m2 * XR) ) / (4.*_s);           
            return pow( pq / szero * XR, l );
        };
        inline std::complex<double> barrier_factor(int l){ return barrier_factor(l, _mX, _mR, _s0); };

        // Phase-space factor, this is defined with respect to a threshold to allow multiple thresholds to be considered
        inline std::complex<double> rho(int l)
        {
            return barrier_factor(l) * sqrt( Kallen(_s * XR, _mX*_mX * XR, _mR*_mR * XR) ) / _s;
        };

        inline std::complex<double> rho_inelastic(int l)
        {
            std::complex<double> result = 0.;
            for (int i = 0; i < _extra_thresholds.size(); i++)
            {
                double m1 = _extra_thresholds[i][0], m2 = _extra_thresholds[i][1];
                result   += _extra_couplings[i] * barrier_factor(l) * sqrt( Kallen(_s * XR, m1*m1 * XR, m2*m2 * XR) ) / _s;
            }

            return result;
        };
    };
};

#endif