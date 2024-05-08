// Minimally gauge invariant pion exchange amplitude for (charged) pion photoproduction
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz-Instituts für Strahlen- und Kernphysik (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef MGI_PION_HPP
#define MGI_PION_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto { namespace piN {

    class mgi_pion : public raw_amplitude 
    {
        public:

        mgi_pion(key k, kinematics xkinem, std::string id = "m.g.i. #pi")
        : raw_amplitude(k, xkinem, id)
        { initialize(0); };       

        inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            store(helicities, s, t);

            _kt   = _kinematics->initial_momentum_tframe(t);
            _pt   = _kinematics->final_momentum_tframe(t);
            _sint = csqrt(_kinematics->Kibble(s,t)/t) / (2*_kt*_pt);

            return bare();
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
                case (kPiPlus)  : {_ePi = + E; return; };
                case (kPiMinus) : {_ePi = - E; return; };
                default: return;
            }
        };

        private:

        inline complex bare()
        {
            if (_lamT != _lamR) return 0.;
            return -4*_ePi*_gNNpi*_lamB*_lamR*_pt*_sint*csqrt(_t)/(_s - _u);
        };
        
        // Couplings
        complex _kt = 0., _pt = 0., _sint = 0.;
        double _ePi   = + E;
        double _gNNpi = 13.48;
    };
}; };

#endif