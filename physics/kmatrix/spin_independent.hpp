// Implementation of a PWA in the scattering-length approximation with up to three
// coupled channels
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef SPIN_INDEPENDENT_HPP
#define SPIN_INDEPENDENT_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "utilities.hpp"
#include "partial_wave.hpp"

namespace jpacPhoto
{
    namespace kmatrix
    {
        struct arguments 
        {
            arguments( int j) : _spin(j) {};

            // Spin of the partial wave
            int  _spin;
            // How many terms to consider in the production vector
            int  _production_expansion  = 1; 
            // How many terms to include in the K-matrix for the diagonal and off-diagonal entries
            int  _diagonal_elastic_expansion = 1, _off_diagonal_elastic_expansion = 1;
            // If this is coupled channel or not
            void add_coupled_channel(double x, double y){ _coupled_channels.push_back({x,y}); };
            std::vector<std::array<double,2>> _coupled_channels;
        };

        class spin_independent : public raw_partial_wave
        {
            public: 

            // Single channel K-matrix
            spin_independent(key k, kinematics xkinem, arguments args)
            : raw_partial_wave(k, xkinem, args._spin, "kmatrix::spin_independent")
            {
                // Populate the thresholds
                _thresholds.push_back({xkinem->get_meson_mass(), xkinem->get_recoil_mass()});
                for (auto extra_threshold : args._coupled_channels)
                {
                    _thresholds.push_back({extra_threshold[0], extra_threshold[1]});
                };

                double nchan = _thresholds.size();
                int offdiags = (nchan+1)*nchan/2;
                _n_prod    = args._production_expansion;
                _n_diag    = args._diagonal_elastic_expansion;
                _n_offdiag = args._off_diagonal_elastic_expansion;
                
                int npars = (_n_prod + _n_diag)*nchan + _n_offdiag*offdiags;
                initialize(npars);
            };

            // -----------------------------------------------------------------------
            // Virtuals 

            // We can have any quantum numbers
            inline std::vector<quantum_numbers> allowed_mesons() { return {ANY}; };
            inline std::vector<quantum_numbers> allowed_baryons(){ return {ANY}; };

            // And helicity independent
            inline helicity_frame native_helicity_frame(){ return HELICITY_INDEPENDENT; };

            // These are projections onto the orbital angular momentum and therefore
            // helicity independent
            inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
            {
                // Save inputes
                store(helicities, s, t);
                
                return (_debug == 1) ? (2*_J+1) * legendre(_J, _kinematics->z_s(s,t)) * partial_wave(_s)
                                     : (2*_J+1) * legendre(_J, cos(_theta))           * partial_wave(_s);
            };

            inline complex partial_wave(double s)
            {
                return 1.;
            };

            inline void set_parameters(std::vector<double> x)
            {
                _production_pars.clear(); _elastic_pars.clear();
                for (int i = 0; i < x.size(); i++)
                {
                    if (x < _thresholds.size()) _production_pars.push_back(x[i]);
                    else                        _elastic_pars.push_back(x[i]);
                };
            };

            protected:

            int _n_prod, _n_diag _n_offdiag;

            // Save the different parameters
            std::vector<double> _production_pars, _elastic_pars;

            // Mass of intermediate coupled channels
            // we can have up to two additional channels
            std::vector<std::array<double,2>> _thresholds;
            
        };
    };
};

#endif