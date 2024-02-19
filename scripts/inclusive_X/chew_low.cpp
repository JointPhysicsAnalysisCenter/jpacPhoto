// Draws the Chew-Low plot (physical region) for the X(3872) inclusive kinematics
//
// OUTPUT: chew_low.pdf 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] 	arXiv:2209.05882 [hep-ph]
// ------------------------------------------------------------------------------

#include "semi_inclusive.hpp"
#include "kinematics.cpp"
#include "semi_inclusive/phase_space.hpp"
#include "plotter.hpp"
#include "print.hpp"

using namespace jpacPhoto;

void chew_low()
{
    using namespace jpacPhoto;
    using complex = std::complex<double>;

    kinematics kX    = new_kinematics(M_X3872);
    semi_inclusive X = new_semi_inclusive<inclusive::phase_space>(kX, M_PROTON + M_PION, "X(3872) Phase space");

    double Wth = M_PROTON + M_PION;
    std::array<double,2> bounds;

    plotter plotter;
    plot p = plotter.new_plot();
    p.set_ranges({-0.1, 35}, {0.8, 3.5});

    auto add_curve = [&](double W)
    {
        auto M2max = [&](double mt)
        {
            return sqrt(X->M2MAXfromT(W*W, -mt));
        };
        auto M2min = [&](double mt)
        {
            return sqrt(X->M2MINfromT(W*W, -mt));
        };
        
        bounds = {-X->TMINfromM2(W*W, Wth*Wth) , -X->TMAXfromM2(W*W, Wth*Wth)};

        p.add_curve(bounds, M2max, var_def("#sqrt{#it{s}}", W, "GeV"));
        p.color_offset(-1);
        p.add_curve(bounds, M2min);
    };

    add_curve(7);
    add_curve(6);
    add_curve(5);

    p.set_legend(0.73, 0.65);
    p.set_labels( "#it{Q}^{2} [GeV^{2}]", "#it{W}  [GeV]");
    p.save("chew_low.pdf");
};