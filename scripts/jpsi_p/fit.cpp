// Script conducting a chi2 minimization fit of K-matrix partial waves
// to j/psi photoproduction data from both GlueX and J/psi-007.
//
// The S-wave model can be selected by uncommenting the appropriate definition.
// Produces fit results and parameters to screen and plots results in region
// of the data.
//
// OUTPUT: gluex_results.pdf, jpsi007_results.pdf
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------


#include "constants.hpp"
#include "kinematics.hpp"
#include "fitter.hpp"
#include "plotter.hpp"

#include "analytic/K_matrix.hpp"
#include "jpsip/gluex/data.hpp"
#include "jpsip/gluex/plots.hpp"
#include "jpsip/jpsi007/data.hpp"
#include "jpsip/jpsi007/plots.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

struct unpolarized
{
    static std::string data_type(int i)
    {
        switch (i)
        {
            case 0: return "integrated (Eg)";
            case 1: return "differential (Eg vs -t)";
            case 2: return "differential (Eg vs -t')";
            default: return "ERROR";
        };
        return "ERROR";
    };  

    static double fcn(const std::vector<data_set> &data_sets, amplitude amp)
    {
        double chi2 = 0;
        for (auto data : data_sets)
        {
            switch (data._type)
            {
                case 0:  { chi2 += chi2_int( data, amp); break; }
                case 1:  { chi2 += chi2_dif( data, amp); break; }
                case 2:  { chi2 += chi2_difp(data, amp); break; }
                default: { chi2 += NaN<double>(); };
            };
        };
        return chi2;
    };

    static double chi2_int(data_set &data, amplitude amp)
    {
        double chi2 = 0;
        for (int i = 0; i < data._N; i++)
        {
            double s = s_cm(data._x[i]);
            double sig_th = amp->integrated_xsection(s);
            double sig_ex = data._z[i];
            double error  = data._zerr[0][i];
            chi2 += pow((sig_th - sig_ex) / error, 2);
        };
        return chi2;
    };

    static double chi2_dif(data_set &data, amplitude amp)
    {
        double chi2 = 0;
        for (int i = 0; i < data._N; i++)
        {
            double s  = s_cm(data._y[i]);
            double t  = - data._x[i];

            double sig_th = amp->differential_xsection(s, t);
            double sig_ex = data._z[i];
            double error  = data._zerr[0][i];
            chi2 += pow((sig_th - sig_ex) / error, 2);
        };
        return chi2;
    };

    static double chi2_difp(data_set &data, amplitude amp)
    {
        double chi2 = 0.;
        for (int i = 0; i < data._N; i++)
        {
            double s  =   s_cm(data._y[i]);
            double t  = - data._x[i] + amp->get_kinematics()->t_min(s);

            double sig_th = amp->differential_xsection(s, t);
            double sig_ex = data._z[i];
            double error  = data._zerr[0][i];
            chi2 += pow((sig_th - sig_ex) / error, 2);
        };
        return chi2;
    };
};

void fit()
{
    using namespace jpacPhoto;
    using K_matrix = analytic::K_matrix;

    // ---------------------------------------------------------------------------
    // Amplitude setup
    // ---------------------------------------------------------------------------

    // J/psi proton final
    kinematics kJpsi = new_kinematics(M_JPSI, M_PROTON);
    kJpsi->set_meson_JP( {1, -1} );

    // Or with two-channels coupled
    std::array<double,2> lower  = {M_D,     M_LAMBDAC};
    std::array<double,2> higher = {M_DSTAR, M_LAMBDAC};
    std::array<std::array<double,2>,2> open_channels = {lower, higher};

    // // CHOOSE AN S-WAVE MODEL

    // Single-channel S-wave
    amplitude s = new_amplitude<K_matrix>(kJpsi, 0, "1-channel S-wave");

    // // Two-channel S-wave 
    // amplitude s = new_amplitude<K_matrix>(kJpsi, 0, higher, "2-channel S-wave");
    
    // // Three-channel S-wave
    // amplitude s = new_amplitude<K_matrix>(kJpsi, 0, open_channels, "3-channel S-wave");

    // EFFECTIVE RANGE?
    s->set_option(K_matrix::kEffectiveRange);

    // The rest of the waves are single channel and 
    amplitude p = new_amplitude<K_matrix>(kJpsi, 1, "P-wave");
    amplitude d = new_amplitude<K_matrix>(kJpsi, 2, "D-wave");
    amplitude f = new_amplitude<K_matrix>(kJpsi, 3, "F-wave");

    amplitude to_fit = s + p + d + f;
    to_fit->set_id("Sum");

    // // ---------------------------------------------------------------------------
    // //  Fitter setup
    // // ---------------------------------------------------------------------------

    fitter<unpolarized> fitter(to_fit, "Migrad", 1.E-4);

    // -----------------------------------------
    // Choose which data we want to fit against

    // All GlueX 2022 data
    std::vector<data_set> gluex   = gluex::all();
    std::vector<data_set> jpsi007 = jpsi007::all(); // And all Jpsi-007 data

    // Load everything into the fitter
    fitter.add_data(gluex);
    fitter.add_data(jpsi007);
    
    // Fit N times with randomly sampled inital parameters
    fitter.do_fit(1);

    // ---------------------------------------------------------------------------
    // Plot the results
    // ---------------------------------------------------------------------------

    plotter plotter;

    // -----------------------------------------
    // GlueX 

    std::vector<plot> gluex_plots;

     // Grab each pre-set plot but add the theory curve with add_curve
    plot pint = gluex::plot_integrated(plotter);

    auto sigma = [&] (double E)
    {
        double W = W_cm(E);
        return to_fit->integrated_xsection(W*W);
    };
    pint.add_curve({8., 11.8}, sigma, to_fit->id());
    gluex_plots.push_back(pint);

    // Do the same with the differential sets
    for (int i = 0; i <= 2; i++)
    {
        double Eavg = gluex[i]._extras[0];
        double Wavg = W_cm(Eavg);
        double tmin = -kJpsi->t_min(Wavg*Wavg);
        double tmax = -kJpsi->t_max(Wavg*Wavg); 
        
        auto dsig_dt = [&](double mt)
        {
            return to_fit->differential_xsection(Wavg*Wavg, -mt);
        };

        plot dif = gluex::plot_slice(plotter, i);
        dif.add_curve({tmin, tmax}, dsig_dt, to_fit->id());

        gluex_plots.push_back(dif);
    };
    
    // Print to file as a 2x2 grid
    plotter.combine({2,2}, gluex_plots, "gluex_results.pdf");
    
    // -----------------------------------------
    // J/psi-007
    
    std::vector<plot> jpsi007_plots;

    for (int i = 1; i <= 12; i++)
    {
        double Eavg = jpsi007[i-1]._extras[0];
        double Wavg = W_cm(Eavg);
        double tmin = -kJpsi->t_min(Wavg*Wavg);
        double tmax = -kJpsi->t_max(Wavg*Wavg);

        auto dsig_dtp = [&](double tp)
        {
            double t = tp + tmin;
            return to_fit->differential_xsection(Wavg*Wavg, -t);
        };

        plot dif =  jpsi007::plot_slice(plotter, i);
        dif.add_curve({0, tmax-tmin}, dsig_dtp, to_fit->id());
        jpsi007_plots.push_back(dif);
    };

    plotter.combine({4,3}, jpsi007_plots, "jpsi007_results.pdf");
};
