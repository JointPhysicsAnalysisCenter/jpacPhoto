// Electric component of the nucleon exchange for (charged) pion photoproduction
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz-Instituts für Strahlen- und Kernphysik (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef MAGNETIC_NUCLEON_HPP
#define MAGNETIC_NUCLEON_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto { namespace piN {

    class magnetic_nucleon : public raw_amplitude 
    {
        public:

        magnetic_nucleon(key k, kinematics xkinem, std::string id = "#it{N} (magnetic)")
        : raw_amplitude(k, xkinem, id)
        { initialize(1); };       

        inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            store(helicities, s, t);

            _kt   = (_t - M2_PION)/(2*csqrt(t));
            _pt   = csqrt(t/4. - M2_PROTON);
            _cost = (_s - _u)/(4*_kt*_pt);
            _sint = csqrt(_kinematics->Kibble(s,t)/t) / (2*_kt*_pt);

            complex dirac_structure = (_lamT == _lamR) ? -2.*(csqrt(t)/2. - _lamT*_lamB*_pt)*_kt*_sint/sqrt(2)
                                                       : -2*_lamT*M_PROTON*_kt*(_cost - _lamT*_lamB)/sqrt(2);
            
            double ff = (_ff) ? _lam2/(_lam2 - _t) : 1.;
            return sqrt(2)*_gNNpi*(_ei/(_s - M2_PROTON) + _ef/(_u - M2_PROTON))*dirac_structure*ff;
        };

        // Everything is evaluated in the t-channel
        inline helicity_frame native_helicity_frame(){ return T_CHANNEL; };
        inline std::vector<quantum_numbers> allowed_mesons() { return { PSEUDOSCALAR }; };
        inline std::vector<quantum_numbers> allowed_baryons(){ return { HALFPLUS }; };

        static const int kPiPlus   = 0;
        static const int kPiMinus  = 1;
        static const int kAddFF    = 2;
        static const int kRemoveFF = 3;
        inline void set_option(int opt)
        {
            switch (opt)
            {
                case (kPiPlus)  : {_ei = +E, _ef =  0; return; };
                case (kPiMinus) : {_ei =  0, _ef = +E; return; };
                case (kAddFF)   : {_ff = true;  set_N_pars(1); return; };
                case (kRemoveFF): {_ff = false; set_N_pars(0); return; };
                default: return;
            }
        };

        inline void allocate_parameters(std::vector<double> pars){ if (_ff) _lam2 = pars[0]; };

        private:
        
        // Couplings
        complex _kt, _pt, _sint, _cost;
        double _ei   = + E, _ef = 0., _lam2 = 1.;
        double _gNNpi = 13.48;
        bool _ff = false;
    };
}; };

#endif