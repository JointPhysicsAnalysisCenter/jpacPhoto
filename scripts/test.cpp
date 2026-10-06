
#include "constants.hpp"
#include "plotter.hpp"
#include "kmatrix/spin_independent.hpp"
#include "TMath.h"
#include "TMatrixD.h"

#include <cstring>
#include <iostream>
#include <iomanip>

void test()
{
    using namespace jpacPhoto;

    kmatrix::arguments swave_args(0);
    swave_args.add_coupled_channel(M_D, M_LAMBDAC);

    kinematics kJpsi = new_kinematics(M_JPSI, M_PROTON);
    kJpsi->set_meson_JP( {1, -1} );

    amplitude swave = new_amplitude<kmatrix::spin_independent>(kJpsi, swave_args);

    TMatrixDSym y(2);
    y[0][0] = 1;
    y[1][1] = 2;
    y[1][0] = 3;

    print("00", y[0][0]);
    print("11", y[1][1]);
    print("10", y[1][0]);
    print("01", y[0][1]);
};