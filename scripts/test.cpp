
#include "constants.hpp"
#include "plotter.hpp"
#include "kmatrix/spin_independent.hpp"
#include "kmatrix/K_matrix.hpp"
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

    amplitude s = new_amplitude<kmatrix::spin_independent>(kJpsi, kmatrix::arguments(0, {1,2,1}));
    s->set_parameters({0.063170013, -418.239953512256, 320.763003700887});
    amplitude p = new_amplitude<kmatrix::spin_independent>(kJpsi, kmatrix::arguments(1));
    p->set_parameters({0.018325241, -133.771684286972});
    amplitude d = new_amplitude<kmatrix::spin_independent>(kJpsi, kmatrix::arguments(2));
    d->set_parameters({0.0030819158, -36.3240467916978});
    amplitude f = new_amplitude<kmatrix::spin_independent>(kJpsi, kmatrix::arguments(3));
    f->set_parameters({0.0008141322, -25.9098777574191});

    amplitude total = s + p + d + f;

    plotter plotter;
    plot p1 = gluex::plot_integrated(plotter);
    // plot p1 = plotter.new_plot();
    p1.add_curve({8, 11.8}, [&](double Eg){ return total->integrated_xsection(s_cm(Eg)); });
    p1.save("swave.pdf");
};