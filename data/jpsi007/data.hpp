// Methods for quickly parsing J/psi-007 2022 data files into the data_set format
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef JPSI007_HPP
#define JPSI007_HPP

#include "data_set.hpp"
#include "constants.hpp"

namespace jpacPhoto
{
    namespace jpsi007
    {
        // Filter the data set by the E_idx setting 
        // which gives slices in energy variable 
        inline data_set slice(int Eidx)
        {
            // Import data and check theres no missing entries
            auto raw  = import_data<15>("data/jpsi007/jpsi007_electron.tsv");
            int  Ntot = check<15>(raw, "J/psi-007 all"); 

            // Destination for the reduced data set after filtering
            int N = 0;
            
            // Calcualte arithmetic mean of the energy values in E-bin
            double Eavg = 0;

            // Target location for 
            std::array<std::vector<double>,6> reduced;
            for (int i = 0; i < Ntot; i++)
            {
                if (raw[2][i] != Eidx) continue;

                // Begin summin # of points and energies to calc average
                N++;
                Eavg += raw[8][i];

                // Populate vectors
                reduced[0].push_back(raw[8][i]);
                reduced[1].push_back(raw[9][i]);
                reduced[2].push_back(raw[9][i] - raw[6][i]);
                reduced[3].push_back(raw[7][i] - raw[9][i]);
                reduced[4].push_back(raw[11][i]);
                reduced[5].push_back(raw[14][i]);
            };

            // Calculate average energy of this slice
            Eavg /= double(N);

            // Wrap the vectors above in a data_set object
            data_set wrapped;

            wrapped._id        = "Jpsi-007 (" + var_def("E", Eavg, "GeV") + ")";
            wrapped._N         = N;
            wrapped._x         = reduced[1];
            wrapped._xerr[0]   = reduced[2];
            wrapped._xerr[1]   = reduced[3];
            wrapped._y         = reduced[0];
            wrapped._z         = reduced[4];
            wrapped._zerr[0]   = reduced[5];
            wrapped._zerr[1]   = reduced[5];
            wrapped._extras.push_back(Eavg);
            wrapped._type = 2;

            return wrapped;
        };

        // Produces a vector containing all slices for the J/psi-007 data set
        inline std::vector<data_set> all()
        {
            std::vector<data_set> result;
            for (int i = 1; i <= 12; i++) result.push_back( jpsi007::slice(i) );
            return result;
        };
    };
};

#endif