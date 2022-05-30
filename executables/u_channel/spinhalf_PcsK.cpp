#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/dirac_exchange.hpp"
#include "amplitudes/vector_exchange.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // Set up kinematics of the Pcs K final state
    double M_PCS = 4.5880;
    reaction_kinematics * kPcsK = new reaction_kinematics(M_KAON, M_PCS);
    kPcsK->set_meson_JP(  PSEUDO_SCALAR );
    kPcsK->set_baryon_JP( HALF_MINUS    );

    // Couplings 
    double gPhoton  =  0.0412501;
    double gNucleon = -6.4469; 
    double cutoff   = 0.9;

    // Create the amplitude
    dirac_exchange LamEx(kPcsK, M_LAMBDA, "#Lambda exchange");
    LamEx.set_params({gPhoton, gNucleon});
    LamEx.set_formfactor(3, cutoff);

    // // ---------------------------------------------------------------------------

    int   N = 50;
    bool PRINT_TO_COMMANDLINE  = true;
    std::string filename  = "PcsK.pdf";

    // ---------------------------------------------------------------------------


    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    double  xmin = 5.;
    double  xmax = 12.;

    double  ymin = 0.;
    double  ymax = 0.8;

    std::string ylabel    = "#sigma(#gamma p #rightarrow P K)  [nb]";
    std::string xlabel    = "W_{#gammap}  [GeV]";

    auto F = [&](double x)
    {
        return LamEx.integrated_xsection(x*x);
    };

    plotter->AddEntry(N, F, {xmin,xmax}, LamEx.get_id(), PRINT_TO_COMMANDLINE);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel,  ymin, ymax);
    plotter->SetLegend(0.7, 0.65);

    // Output to file
    plotter->Plot(filename);

};