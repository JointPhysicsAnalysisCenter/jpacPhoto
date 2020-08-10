// Sample event-generator for the gamma p -> R -> Jpsi p -> l+ l- p 
// reaction. Events are generated on a flat phase space and weighted by the probabilty distribution
// of a supplied amplitude.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _TOY_MC_
#define _TOY_MC_

#include "amplitudes/amplitude.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

namespace jpacPhoto
{
    class toy_monte_carlo
    {
        public:
        // constructor with only filename 
        toy_monte_carlo(std::string filename = "mc.root")
        {
            rnd = new TRandom3(0);
            outfile = new TFile(filename.c_str(), "RECREATE");
            set_up_kin();
        };

        // constructor if you want to specify seed
        toy_monte_carlo(int seed = 0, std::string filename = "mc.root")
        {
            rnd = new TRandom3(seed);
            outfile = new TFile(filename.c_str(), "RECREATE");
            set_up_kin();
        };

        // destructor
        ~toy_monte_carlo()
        {
            delete outfile;
            delete rnd;
        }
        
        inline void set_amplitude(amplitude * _amp)
        {
            amp = _amp;
            set_up_dyn();
        };
        
        // Produce root file with N events
        void generate(double _beam_energy, int N);

        private:

        // Set up the branch structure for the output root file
        void set_up_kin();
        void set_up_dyn();

        // Shorthand for generating a random number in a range 
        inline double random(double min, double max)
        {
            return min + (max - min) * rnd->Rndm();
        }
        
        // Generate event
        void generate_event();  // kinematics
        void generate_weight(); // dynamics

        // ROOT structures
        TFile * outfile = NULL;
        TTree * kin = NULL;
        TTree * dyn = NULL;
        TRandom3 * rnd = NULL;
        
        // Pointer to preset amplitude that will be the weighting function for events
        amplitude * amp = NULL;
        bool error_already_triggered = false;
        double beam_energy, W, s, t;
        double weight;

        // Helicities of the photon and 2 x lambda of the target and recoil protons
        int lam_gamma, lam_ptarg, lam_erel, lam_prec;

        // Lepton 4-momenta components 
        double ep_px, ep_py, ep_pz, ep_E;
		double em_px, em_py, em_pz, em_E;

        // Recoil and target proton components
		double prec_px, prec_py, prec_pz, prec_E;
		double ptarg_px, ptarg_py, ptarg_pz, ptarg_E;

        // Beam photon 4-momenta 
		double pgamma_px, pgamma_py, pgamma_pz, pgamma_E;

        // Angles
        double theta_ep, phi_ep; // of the lepton pair
        double theta_psi, phi_psi; // of the produced jpsi
    };
};

#endif 