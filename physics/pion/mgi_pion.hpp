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
#include "cgamma.hpp"

#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace jpacPhoto { namespace piN {

    class mgi_pion : public raw_amplitude 
    {
        public:

        mgi_pion(key k, kinematics xkinem, int i = kFixedSpin, std::string id = "m.g.i. #pi")
        : raw_amplitude(k, xkinem, id)
        { initialize(0); set_option(i); };       

        inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            store(helicities, s, t);

            if (_lamT != _lamR) return 0.;

            _kt   = (_t - M2_PION)/(2*csqrt(t));
            _pt   = csqrt(t/4. - M2_PROTON);
            _cost = (_s - _u)/(4*_kt*_pt);
            _sint = csqrt(_kinematics->Kibble(s,t)/t) / (2*_kt*_pt);

            return -2*_ePi*_gNNpi*_lamB*_lamT*_t*_cost/_sint*propagator();
        };

        // Everything is evaluated in the t-channel
        inline helicity_frame native_helicity_frame(){ return T_CHANNEL; };
        inline std::vector<quantum_numbers> allowed_mesons() { return { PSEUDOSCALAR }; };
        inline std::vector<quantum_numbers> allowed_baryons(){ return { HALFPLUS }; };

        static const int kPiPlus           = 0;
        static const int kPiMinus          = 1;
        static const int kFixedSpin        = 2;
        static const int kHighEnergyLimit  = 3;
        static const int kSingleJPole      = 4;
        static const int kResummed         = 5;
        inline void set_option(int opt)
        {
            switch (opt)
            {
                case (kPiPlus)  : {_ePi = + E; return; };
                case (kPiMinus) : {_ePi = - E; return; };
                case (kFixedSpin)       : { set_N_pars(0); _option = opt; return;};
                case (kSingleJPole)     : { set_N_pars(0); _option = opt; return;}; 
                case (kHighEnergyLimit) : { set_N_pars(1); _option = opt; return;};
                case (kResummed)        : { set_N_pars(3); _option = opt; return;};
                default: option_error();
            }
            return;
        };

        inline void allocate_parameters(std::vector<double> pars)
        {
            if (_option == kFixedSpin || _option == kSingleJPole) return;
            _r2 = pars[0];
            if (_option != kResummed) return;
            _jp = pars[1]; _jz = pars[2];
        };  

        private:
        
        inline complex propagator()
        {
            switch (_option)
            {
                case kFixedSpin:       return bare_propagator();      
                case kHighEnergyLimit: return asymptotic_propagator();
                case kSingleJPole:     return single_Jpole();         
                case kResummed:        return resummed_propagator();
                default: return NaN<complex>();
            }
        };

        inline complex bare_propagator()      { return -1./(_t - M2_PION); };
        inline complex single_Jpole()         { return -_aP/alpha(); };
        inline complex asymptotic_propagator(){ return 2*_aP/sqrt(PI)*cgamma(alpha()+3/2.)*signature()*cgamma(-alpha())*pow(2*_s*_r2, alpha()); };

        inline complex resummed_propagator()
        {
            complex J0 = single_Jpole();
            complex K  = -2.*_pt*_kt*_r2;
            double a   = alpha();

            auto integrand = [&](double y)
            {
                complex term1 = -1/a;
                complex term2 = _jp*(a+_jz)*(a+1)*(2*a+1)/a/_jz/(a+_jp)  *pow(y,  -a);
                complex term3 = (_jz - _jp)*(1-_jp)*(1-2*_jp)/_jz/(_jp+a)*pow(y, _jp);

                complex Gp = sqrt(1-2*_cost*K*y+K*K*y*y);
                complex Gm = sqrt(1+2*_cost*K*y+K*K*y*y);

                complex resum = 2/Gp/(pow(1+Gp, 2) - K*K*y*y) - 2/Gm/(pow(1+Gm, 2) - K*K*y*y);

                return (term1 + term2 + term3) * resum;
            };
            complex integral = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(integrand, 0, 1, 25, 1.E-9, NULL);

            return (J0 - _sint*_sint/_cost*_aP*K*integral/2.);
        };


        inline double alpha()     { return _aP*(_t - M2_PION); };
        inline complex signature(){ return (1 + exp(-I*PI*alpha()))/2.; };
        
        // Couplings
        complex _kt = 0., _pt = 0., _sint = 0., _cost = 0.;
        double _ePi   = + E;
        double _gNNpi = 13.48;
        double _aP = 0.7, _r2 = 1/2;
        double _jp, _jz;
    };
}; };

#endif