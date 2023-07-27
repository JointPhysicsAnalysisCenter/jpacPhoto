// Semi-inclusive production of axial vectors vector exchange using proton structure
// functions
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#ifndef DIS_LIKE_HPP       
#define DIS_LIKE_HPP

#include "constants.hpp"
#include "inclusive_process.hpp"
#include "sigma_tot/ChristyBosted_F.hpp"

namespace jpacPhoto
{
    namespace inclusive
    {
        class DIS_like : public raw_inclusive_process
        {
            public: 

            DIS_like(inclusive_key key, double mX, std::string id = "")
            : raw_inclusive_process(key, mX, id),
              F1(new_inclusive_function<ChristyBosted_F>(1, kProton)),
              F2(new_inclusive_function<ChristyBosted_F>(2, kProton))
            {
                set_N_pars(4);
            };

            DIS_like(inclusive_key key, double mX, int pn, std::string id = "")
            : raw_inclusive_process(key, mX, id),
              F1(new_inclusive_function<ChristyBosted_F>(1, pn)),
              F2(new_inclusive_function<ChristyBosted_F>(2, pn))
            {
                set_N_pars(4);
            };

            // Minimum mass is the proton 
            inline double minimum_M2(){ return pow(M_PROTON + M_PION, 2); };

            // Only free parameters are the top photocoupling and the form factor cutoff
            inline void allocate_parameters(std::vector<double> pars)
            {
                _g     = pars[0];
                _eta   = pars[1];
                _mEx2  = pars[2]*pars[2];
                _lam2  = pars[3]*pars[3];
            };

            // The invariant cross section used S T and M2 as independent vareiables
            inline double invariant_xsection(double s, double t, double M2)
            {
                update(s, t, M2);

                // Contraction of top and bottom tensors
                double TdotW; 
                TdotW  = 3*_F1*_T1;
                TdotW += (_kdotq*_F1*_T2 + _pdotq*_T1*_F2)      /t;
                TdotW += (_kdotq*_pdotq  - 2*t*_pdotk)*_T2*_F2/t/t;
                TdotW += (_pdotk*_pdotk*_T2*_F2 - M2_PROTON*_kdotq*_T1*_F2)/_pdotq/_kdotq;

                // Flux factor
                double flux = 1/(2*sqrt(s)*qGamma());

                // Exchange propagator
                double P_ex    = 1/(_mEx2 - t);

                // Form factor in the case of massive vectors
                double tprime  = t - TMINfromM2( M2_PROTON );
                double beta_ex = (is_zero(_mEx2)) ? 1. : exp(tprime/_lam2)/pow(1-tprime/0.71,-2);

                // in nanobarn!!!!!
                return flux * pow(P_ex * beta_ex * _eta*_eta * E, 2) * TdotW / (8*PI*PI) / (2.56819E-6); 
            };

            // Options select proton or neutron target
            static const int kProton = 0, kNeutron = 1;
            inline void set_option (int opt)
            {
                F1 = new_inclusive_function<ChristyBosted_F>(1, opt);
                F2 = new_inclusive_function<ChristyBosted_F>(2, opt);
                return;
            };

            protected:

            // Calculate dot products and form factors
            inline void update(double s, double t, double M2)
            {
                // Dot products of all the relevant momenta
                _pdotk = (s  - M2_PROTON)     / 2;
                _pdotq = (M2 - M2_PROTON - t) / 2;
                _kdotq = (t - _mX2)           / 2;

                // Production form factors
                _prefactors = _g*_g/2*t*t/_mX2/_mX2/_mX2;
                _T1 = _prefactors * _kdotq*_kdotq;
                _T2 = _prefactors * _kdotq*(_mX2 - 2*_kdotq);

                // Hadronic structure function
                _F1 = F1->evaluate(M2, t);
                _F2 = F2->evaluate(M2, t);
            };  

            // If we are reggeized we use the high-energy approximation and use t & x
            inline bool use_TX(){ return false; };

            private:

            // Free parameters
            double _g    = 0; // Top couplings
            double _eta  = 1; // VMD coupling for exchange
            double _mEx2 = 0; // Mass (squared) of exchange
            double _lam2 = 1; // Form factor cutoff in GeV

            // Internal variables
            double _pdotk, _pdotq, _kdotq;
            double _prefactors, _T1, _T2;
            double _F1, _F2;

            // Proton form factors
            inclusive_function F1, F2;
        };
    };
};

#endif