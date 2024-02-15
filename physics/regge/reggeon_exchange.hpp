// Generic reggeon exchange model with simple vertices
// Adapted from the model considered in [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------
// References:
// [1] - https://arxiv.org/abs/1710.09394
// ------------------------------------------------------------------------------

#ifndef REGGEON_EXCHANGE_HPP
#define REGGEON_EXCHANGE_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto
{
    namespace regge
    {
        class reggeon_exchange : public raw_amplitude
        {
            public:

            // Constructor
            reggeon_exchange(key k, kinematics xkinem, int naturality, std::string id)
            : raw_amplitude(k, xkinem, id), _naturality(naturality)
            {  
                if (abs(naturality) > 1) warning("reggeon_exchange", "Invalid naturality passed to constructor!");
                initialize(8);
            };

            // ---------------------------------------------------------------------------
            // Defining the virtual functions required of an amplitude
            
            inline complex helicity_amplitude(std::array<int, 4> helicities, double s, double t)
            {
                // Save inputs
                store(helicities, s, t);
                auto result = exp(_b*t)*top()*propagator()*bottom();
                return result;
            };

            // Explicitly require t-channel helicities
            inline helicity_frame native_helicity_frame()        { return   S_CHANNEL;  };
            inline std::vector<quantum_numbers> allowed_mesons() { return {  VECTOR };  };
            inline std::vector<quantum_numbers> allowed_baryons(){ return { HALFPLUS }; };

            // -----------------------------------------------------------------------
            // Internal data members 

            protected:

            // Set parameters
            inline void allocate_parameters(std::vector<double> x)
            {
                _alpha0 = x[0];
                _alphaP = x[1];
                _b      = x[2];
                _gT[0]  = x[3]; _gT[1] = x[4]; _gT[2] = x[5];        
                _gB[0]  = x[6]; _gB[1] = x[7];         
            };

            //----------------------------------------------------
            int _signature  = +1; // Paper explicitly fixes all trajectories to + signature
            int _naturality = +1;
            double _s0 = 1;

            // Couplings
            double _b = 0.;           // Form factor cutoff
            std::array<double,3> _gT; // non-flip, single flip, double flip
            std::array<double,2> _gB; // non-flip and flip

            // Linear trajectory
            double _alpha0 = 0, _alphaP = 0;

            //----------------------------------------------------

            inline double trajectory(){ return _alpha0 + _alphaP * _t; };
            
            // Regge propagator 
            inline complex propagator()
            {
                double alpha = _alpha0 + _alphaP * _t;

                // Ghost killing factor for natural exchanges
                double gf = (_naturality == +1) ? alpha / _alpha0 : 1.; 
                if (is_zero(gf)) gf = 1/_alpha0/PI;
                else gf /= sin(PI*alpha);

                // signature factor
                complex sigf = (_signature + exp(-I*PI*alpha))/2.;

                return -_alphaP*sigf*gf*PI*pow(_s/_s0, alpha);
            };

            double top()
            { 
                int li = abs(_lamB - _lamX);
                int phase = pow(-_lamB, li + (_naturality < 0));
                return phase*_gT[li]*pow( sqrt(-_t)/_mX, li);
            };

            double bottom()
            {
                int lf = abs(_lamT - _lamR)/2;
                int phase = pow(-_lamT, lf + (_naturality < 0));
                return phase*_gB[lf] * pow( sqrt(-_t) / (_mT + _mR), lf);
            };
        };
    };
};

#endif