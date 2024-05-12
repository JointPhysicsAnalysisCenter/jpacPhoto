// Pion photoproduction amplitude from [1]
// Instead of spliting different exchange contributions we just include all gauge invariant diagrams in one
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz-Instituts für Strahlen- und Kernphysik (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// References:
//
// [1] - 10.1016/S0375-9474(97)00612-X
// ------------------------------------------------------------------------------

#ifndef VGL_HPP
#define VGL_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto { namespace piN {

    class VGL : public raw_amplitude
    {
        public:

        VGL(key k, kinematics xkinem, std::string id = "VGL")
        : raw_amplitude(k, xkinem, id)
        { initialize(0); };

        inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            store(helicities, s, t);
            _covariants->update(helicities, s, t);

            complex reggeize = (_regge) ? (_t - M2_PION)*regge_propagator() : 1;
            return electric()*reggeize + magnetic();
        };

        // Covariants are always natively s-channel helicity amplitudes
        inline helicity_frame native_helicity_frame()        { return S_CHANNEL; };
        inline std::vector<quantum_numbers> allowed_mesons() { return { PSEUDOSCALAR }; };
        inline std::vector<quantum_numbers> allowed_baryons(){ return { HALFPLUS    }; };

        static const int kPiPlus   = 0;
        static const int kPiMinus  = 1;
        static const int kSimple   = 2;
        static const int kRegge    = 3;
        inline void set_option(int opt)
        {
            switch (opt)
            {
                case (kPiPlus)  : {_ei = +E, _ef =  0, _ep = +E; return; };
                case (kPiMinus) : {_ei =  0, _ef = +E, _ep = -E; return; };
                case (kSimple)  : {_regge = false; return; };
                case (kRegge)   : {_regge = true; return; };
                default: return;
            }
        };

        private:
        
        // Particle charges
        double _ei   = +E, _ef = 0., _ep = +E;

        // Pion nucleon coupling
        double _gNNpi = 13.48;
        double _s0    = 1.;

        // If reggeizing the pion 
        bool _regge = true;

        inline complex magnetic()
        {
            auto u     = _covariants->u();    // Target spinor
            auto ubar  = _covariants->ubar(); // Recoil spinor

            auto k     = _covariants->q();    // Photon momentum
            auto eps   = _covariants->eps();  // Photon polarization

            complex propagators = _ei/(_s - M2_PROTON) + _ef/(_u - M2_PROTON);

            return sqrt(2)*_gNNpi*propagators*contract(ubar, gamma_5()*slash(k)*slash(eps)*u);
        };

        inline complex electric()
        {
            auto u     = _covariants->u();        // Target spinor
            auto ubar  = _covariants->ubar();     // Recoil spinor

            auto eps   = _covariants->eps();      // Photon polarization
            auto p_pi  = _covariants->q_prime();  // Photon momentum
            auto p_i   = _covariants->p();        // incoming nucleon
            auto p_f   = _covariants->p_prime();  // outgoing nucleon

            complex propagators = _ep*contract(eps, p_pi)/(_t - M2_PION) 
                                + _ei*contract(eps, p_i) /(_s - M2_PROTON) 
                                + _ef*contract(eps, p_f) /(_u - M2_PROTON);

            return 2.*sqrt(2)*_gNNpi*propagators*contract(ubar, gamma_5()*u);
        };

        inline complex regge_propagator()
        {
            double alphaP = 0.7;
            double alpha  = alphaP*(_t - M2_PION);
            complex signature = (+1 + exp(-I*PI*alpha))/2;
            return -alphaP*signature*cgamma(-alpha)*pow(_s/_s0, alpha);
        };
    };
}; };

#endif