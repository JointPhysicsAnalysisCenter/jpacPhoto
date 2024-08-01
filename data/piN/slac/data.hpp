// Interface functions for the Boyarski and GlueX data of pi delta photoproduction
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Universitat Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------

#ifndef PIN_DATA_HPP
#define PIN_DATA_HPP

#include "data_set.hpp"
#include "constants.hpp"

namespace jpacPhoto { namespace piN 
{
    inline data_set differential(int i)
    {
        std::string id = "SLAC (1968)";

        std::string filename = "/data/piN/slac/Boyarski_PRL20_table" + std::to_string(i) + ".txt";
        auto raw = import_data<6>(filename);
        int  N   = check<6>(raw, id);

        data_set out;
        out._id   = id;
        out._N    = N;
        out._x    = - raw[0];
        out._z    = raw[3];
        out._zerr = {raw[4],  -raw[5]};

        return out;
    };
    
    inline std::vector<data_set> differential()
    {
        std::vector<data_set> out;
        for (int i = 1; i <= 4; i++) out.push_back(differential(i));
        return out;
    };

}; };

#endif