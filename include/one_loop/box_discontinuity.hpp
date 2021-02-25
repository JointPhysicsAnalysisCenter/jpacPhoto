// Discontinuity across unitarity cut for vector production via a one-loop box diagram
// Used as a container class in integration processes
// 
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _BOX_DISC_
#define _BOX_DISC_

#include "constants.hpp"
#include "amplitudes/amplitude.hpp"
#include "amplitudes/reaction_kinematics.hpp"

#include "Math/IntegratorMultiDim.h"

namespace jpacPhoto
{
    class box_discontinuity
    {
        public: 
        box_discontinuity(amplitude * left, amplitude * right)
        : _initialAmp(left), _finalAmp(right)
        {
            // Make sure the left and right amplitudes match!
            // Check the spins of the intermediate state
            _jp_left  = left->_kinematics->_jp;
            _jp_right = right->_kinematics->_jp;

            if (_jp_left != _jp_right) 
            {
                std::cout << std::left;
                std::cout << std::setw(40) << "box_amplitude: Intermediate state between sub-amplitudes dont match!" << std::endl;
                std::cout << std::setw(20) << _initialAmp->_identifier << ": \t (" << _jp_left[0]  << ", " << _jp_left[1]  << ")\n";
                std::cout << std::setw(20) << _finalAmp->_identifier   << ": \t (" << _jp_right[0] << ", " << _jp_right[1] << ")\n";
                std::cout << "Returning 0!" << std::endl;

                _matchError = true;
            };

            // But also masses
            if ((std::abs(left->_kinematics->_mX - right->_kinematics->_mX) > 1.E-4) 
             || (std::abs(left->_kinematics->_mR - right->_kinematics->_mR) > 1.E-4))
            {
                _matchError = true;
            };
            
            // IF they match, get the spin and therefor helicities of the intermediate meson
            _intermediate_helicities = get_helicities(_jp_left[0]);
        };

        // Evaluate the discontinuity integrated over intermediate phase space
        double eval(double s);

        // Pass values that come from the external gamma p -> V p reaction
        inline void set_externals(std::array<int,4> helicities, double theta)
        {
            _external_theta = theta; 
            _external_helicities = helicities;
        };

        private:
        double _external_theta;
        std::array<int,4> _external_helicities;

        // Individual tree amplitudes
        amplitude * _initialAmp;
        amplitude * _finalAmp;

        std::array<int,2> _jp_left, _jp_right;
        std::vector< std::array<int,4> > _intermediate_helicities;
        bool _matchError = false;
    };
};

#endif