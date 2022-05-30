// Arrays containing helicity combinations for indexing.
// Moved to a seperate file to not clutter up reaction_kinematics.hpp
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef HELIC_COMBO
#define HELIC_COMBO

#include <algorithm>
#include <iostream>
#include <vector>
#include <array>

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Massless helicity combinations
    const std::vector< std::array<int, 4> > SPIN_ZERO_HELICITIES =
    {
    //  {  γ,  p,  S,  p'}
        {  1, -1,  0, -1}, // 0
        {  1, -1,  0,  1}, // 1
        {  1,  1,  0, -1}, // 2
        {  1,  1,  0,  1}, // 3
        { -1, -1,  0, -1}, // 4
        { -1, -1,  0,  1}, // 5
        { -1,  1,  0, -1}, // 6
        { -1,  1,  0,  1}  // 7
    };

    // POS / NEG refer to ordering whether hel[0] is +-1 everything else equal
    const std::vector<int> SPIN_ZERO_POS_ITERS = {0, 1, 2, 3};
    const std::vector<int> SPIN_ZERO_NEG_ITERS = {4, 5, 6, 7};

    const std::vector< std::array<int, 4> > SPIN_ONE_HELICITIES =
    {
    //  {  γ,  p,  V,  p'}
        {  1, -1,  1, -1}, // 0
        {  1, -1,  1,  1}, // 1
        {  1, -1,  0, -1}, // 2
        {  1, -1,  0,  1}, // 3
        {  1, -1, -1, -1}, // 4
        {  1, -1, -1,  1}, // 5
        {  1,  1,  1, -1}, // 6
        {  1,  1,  1,  1}, // 7
        {  1,  1,  0, -1}, // 8
        {  1,  1,  0,  1}, // 9
        {  1,  1, -1, -1}, // 10
        {  1,  1, -1,  1}, // 11
        { -1, -1,  1, -1}, // 12
        { -1, -1,  1,  1}, // 13
        { -1, -1,  0, -1}, // 14
        { -1, -1,  0,  1}, // 15
        { -1, -1, -1, -1}, // 16
        { -1, -1, -1,  1}, // 17
        { -1,  1,  1, -1}, // 18
        { -1,  1,  1,  1}, // 19
        { -1,  1,  0, -1}, // 20
        { -1,  1,  0,  1}, // 21
        { -1,  1, -1, -1}, // 22
        { -1,  1, -1,  1}  // 23
    };

    // Correspond to amplitudes with hel[2] = +1
    // POS / NEG refer to ordering whether hel[0] is +-1 everything else equal
    const std::vector<int> SPIN_ONE_POS_ITERS = {0, 1, 6, 7, 12, 13, 18, 19};
    const std::vector<int> SPIN_ONE_NEG_ITERS = {12, 13, 18, 19, 0, 1, 6, 7};

    const std::vector< std::array<int, 4> > SPIN_TWO_HELICITIES =
    {
    //  {  γ,  p,  V,  p'}
        {  1, -1,  2, -1}, //  0
        {  1, -1,  2,  1}, //  1
        {  1, -1,  1, -1}, //  2
        {  1, -1,  1,  1}, //  3
        {  1, -1,  0, -1}, //  4
        {  1, -1,  0,  1}, //  5
        {  1, -1, -1, -1}, //  6
        {  1, -1, -1,  1}, //  7
        {  1, -1, -2, -1}, //  8
        {  1, -1, -2,  1}, //  9
        {  1,  1,  2, -1}, // 10
        {  1,  1,  2,  1}, // 11
        {  1,  1,  1, -1}, // 12
        {  1,  1,  1,  1}, // 13
        {  1,  1,  0, -1}, // 14
        {  1,  1,  0,  1}, // 15
        {  1,  1, -1, -1}, // 16
        {  1,  1, -1,  1}, // 17
        {  1,  1, -2, -1}, // 18
        {  1,  1, -2,  1}, // 19
        { -1, -1,  2, -1}, // 20
        { -1, -1,  2,  1}, // 21
        { -1, -1,  1, -1}, // 22
        { -1, -1,  1,  1}, // 23
        { -1, -1,  0, -1}, // 24
        { -1, -1,  0,  1}, // 25
        { -1, -1, -1, -1}, // 26
        { -1, -1, -1,  1}, // 27
        { -1, -1, -2, -1}, // 28
        { -1, -1, -2,  1}, // 29
        { -1,  1,  2, -1}, // 30
        { -1,  1,  2,  1}, // 31
        { -1,  1,  1, -1}, // 32
        { -1,  1,  1,  1}, // 33
        { -1,  1,  0, -1}, // 34
        { -1,  1,  0,  1}, // 35
        { -1,  1, -1, -1}, // 36
        { -1,  1, -1,  1}, // 37
        { -1,  1, -2, -1}, // 38
        { -1,  1, -2,  1}  // 39
    };

    // Correspond to amplitudes with hel[2] = +1
    // POS / NEG refer to ordering whether hel[0] is +-1 everything else equal
    const std::vector<int> SPIN_TWO_POS_ITERS = {2, 3, 12, 13, 22, 23, 32, 33};
    const std::vector<int> SPIN_TWO_NEG_ITERS = {22, 23, 32, 33, 2, 3, 12, 13};

    // ---------------------------------------------------------------------------
    // If we allow massive "beam"
    const std::vector< std::array<int, 4> > MASSIVE_SPIN_ZERO_HELICITIES =
    {
    //  {  γ,  p,  S,  p'}
        {  1, -1,  0, -1},
        {  1, -1,  0,  1},
        {  1,  1,  0, -1},
        {  1,  1,  0,  1},
        {  0, -1,  0, -1},
        {  0, -1,  0,  1},
        {  0,  1,  0, -1},
        {  0,  1,  0,  1},
        { -1, -1,  0, -1},
        { -1, -1,  0,  1},
        { -1,  1,  0, -1},
        { -1,  1,  0,  1}
    };
  
    const std::vector< std::array<int, 4> > MASSIVE_SPIN_ONE_HELICITIES =
    {
    //  {  γ,  p,  V,  p'}
        {  1, -1,  1, -1},
        {  1, -1,  1,  1},
        {  1, -1,  0, -1},
        {  1, -1,  0,  1},
        {  1, -1, -1, -1},
        {  1, -1, -1,  1},
        {  1,  1,  1, -1},
        {  1,  1,  1,  1},
        {  1,  1,  0, -1},
        {  1,  1,  0,  1},
        {  1,  1, -1, -1},
        {  1,  1, -1,  1},
        {  0, -1,  1, -1},
        {  0, -1,  1,  1},
        {  0, -1,  0, -1},
        {  0, -1,  0,  1},
        {  0, -1, -1, -1},
        {  0, -1, -1,  1},
        {  0,  1,  1, -1},
        {  0,  1,  1,  1},
        {  0,  1,  0, -1},
        {  0,  1,  0,  1},
        {  0,  1, -1, -1},
        {  0,  1, -1,  1},
        { -1, -1,  1, -1},
        { -1, -1,  1,  1},
        { -1, -1,  0, -1},
        { -1, -1,  0,  1},
        { -1, -1, -1, -1},
        { -1, -1, -1,  1},
        { -1,  1,  1, -1},
        { -1,  1,  1,  1},
        { -1,  1,  0, -1},
        { -1,  1,  0,  1}, 
        { -1,  1, -1, -1},
        { -1,  1, -1,  1}
    };

    inline std::vector<std::array<int, 4>> get_helicities(int mJ, int bJ, bool is_massless = true)
    {
        int mjbj = 100 * is_massless + 10 * mJ + bJ;

        switch (mjbj)
        {
            case (121): return SPIN_TWO_HELICITIES;
            case (111): return SPIN_ONE_HELICITIES;
            case (101): return SPIN_ZERO_HELICITIES;
            case ( 11): return MASSIVE_SPIN_ONE_HELICITIES;
            case (  1): return MASSIVE_SPIN_ZERO_HELICITIES;
            default: 
            {
                std::cout << "Fatal error! Can't find helicities for meson spin J = " << mJ << " and baryon spin " << bJ << "/2 not yet implemented.\n";
                std::cout << "Quitting...\n";
                exit(1);
            }
        }
    };

    inline std::array<std::vector<int>, 2> get_iters(int mJ, int bJ, bool massless = true)
    {
        int mjbj = 100 * !massless + 10 * mJ + bJ;
        switch (mjbj)
        {
            case (  1): return {SPIN_ZERO_POS_ITERS, SPIN_ZERO_NEG_ITERS};
            case ( 11): return {SPIN_ONE_POS_ITERS,  SPIN_ONE_NEG_ITERS};
            case ( 21): return {SPIN_TWO_POS_ITERS,  SPIN_TWO_NEG_ITERS};
            default:
            {
                std::cout << "Fatal error! Can't find iters for meson spin J = " << mJ << " and baryon spin " << bJ << "/2 not yet implemented.\n";
                std::cout << "Quitting...\n";
                exit(1);
            }
        };
    };

    inline int find_helicity(std::array<int, 4> helicities, int mj, int bj, bool is_massless = true)
    {
        std::vector<std::array<int,4>> hels = get_helicities(mj, bj, is_massless);

        auto iterator = std::find(hels.begin(), hels.end(), helicities);
        
        if (iterator != hels.end())
        {
            return iterator - hels.begin();
        }
        else
        {
            std::cout << "Error cannot find helicities! Returning -1... \n";  
            return -1;
        }
    };
};

#endif