
#include "constants.hpp"
#include "plotter.hpp"
#include "kmatrix/spin_independent.hpp"
#include "jpsip/gluex/plots.hpp"
#include <Eigen/Dense>

#include <cstring>
#include <iostream>
#include <iomanip>

void test()
{
    using namespace jpacPhoto;

    kinematics kJpsi = new_kinematics(M_JPSI, M_PROTON);
    kJpsi->set_meson_JP( {1, -1} );

    // // Single channel
    // amplitude s_1C = new_amplitude<kmatrix::spin_independent>(kJpsi, kmatrix::arguments(0, {1,2,1}));
    // s_1C->set_parameters({0.063170013, -418.239953512256, 320.763003700887});
    // amplitude p_1C = new_amplitude<kmatrix::spin_independent>(kJpsi, kmatrix::arguments(1));
    // p_1C->set_parameters({0.018325241, -133.771684286972});
    // amplitude d_1C = new_amplitude<kmatrix::spin_independent>(kJpsi, kmatrix::arguments(2));
    // d_1C->set_parameters({0.0030819158, -36.3240467916978});
    // amplitude f_1C = new_amplitude<kmatrix::spin_independent>(kJpsi, kmatrix::arguments(3));
    // f_1C->set_parameters({0.0008141322, -25.9098777574191});
    // amplitude sum_1C = s_1C + p_1C + d_1C + f_1C;

    // J/psi p & D* LambdaC
    kmatrix::arguments args_2C = kmatrix::arguments(0, {1,2,1});
    args_2C.add_coupled_channel(M_DSTAR, M_LAMBDAC);

    amplitude s_2C = new_amplitude<kmatrix::spin_independent>(kJpsi, args_2C);
    s_2C->set_parameters({0.10122243, 3.2136095, -219.68553, -181.30566, 47.09921, -145.67629, 4.9978015});
    // amplitude p_2C = new_amplitude<kmatrix::spin_independent>(kJpsi, kmatrix::arguments(1));
    // p_2C->set_parameters({0.014623196,   -43.997429});
    // amplitude d_2C = new_amplitude<kmatrix::spin_independent>(kJpsi, kmatrix::arguments(2));
    // d_2C->set_parameters({0.0030291014,  -2.3363761});
    // amplitude f_2C = new_amplitude<kmatrix::spin_independent>(kJpsi, kmatrix::arguments(3));
    // f_2C->set_parameters({0.00068965326, -6.0145674});
    // amplitude sum_2C = s_2C + p_2C + d_2C + f_2C;

    plotter plotter;
    plot p1 = gluex::plot_integrated(plotter);
    // p1.add_curve({8, 11.8}, [&](double Eg){ return sum_1C->integrated_xsection(s_cm(Eg)); }, solid(jpacColor::Blue, "Single channel (1C)"));
    p1.add_curve({8, 11.8}, [&](double Eg){ return s_2C->integrated_xsection(s_cm(Eg)); }, solid(jpacColor::Red,  "Two channels (2C)"));
    p1.save("swave.pdf");
};