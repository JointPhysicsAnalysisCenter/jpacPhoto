// Electric component of the nucleon exchange for (charged) pion photoproduction
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz-Instituts für Strahlen- und Kernphysik (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef ELECTRIC_NUCLEON_HPP
#define ELECTRIC_NUCLEON_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto { namespace piN {

    class electric_nucleon : public raw_amplitude 
    {
        public:

        electric_nucleon(key k, kinematics xkinem, std::string id = "#it{N} (electric)")
        : raw_amplitude(k, xkinem, id)
        { initialize(0); };       

        inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            store(helicities, s, t);

            _kt   = _kinematics->initial_momentum_tframe(t);
            _pt   = _kinematics->final_momentum_tframe(t);
            _sint = csqrt(_kinematics->Kibble(s,t)/t) / (2*_kt*_pt);

            if (_lamT != _lamR) return 0.;
            return -2*(_ei*(_t - M2_PION)/(_s - M2_PROTON) - _ef*(_t - M2_PION)/(_u - M2_PROTON))*_gNNpi*_lamB*_lamR*_pt*_sint*csqrt(_t)/(_s - _u);
        };

        // Everything is evaluated in the t-channel
        inline helicity_frame native_helicity_frame(){ return T_CHANNEL; };
        inline std::vector<quantum_numbers> allowed_mesons() { return { PSEUDOSCALAR }; };
        inline std::vector<quantum_numbers> allowed_baryons(){ return { HALFPLUS }; };

        static const int kPiPlus  = 0;
        static const int kPiMinus = 1;
        inline void set_option(int opt)
        {
            switch (opt)
            {
                case (kPiPlus)  : {_ei = +E, _ef =  0; return; };
                case (kPiMinus) : {_ei =  0, _ef = +E; return; };
                default: return;
            }
        };

        private:
        
        // Couplings
        complex _kt = 0., _pt = 0., _sint = 0.;
        double _ei   = + E, _ef = 0.;
        double _gNNpi = 13.48;
    };
}; };

#endif