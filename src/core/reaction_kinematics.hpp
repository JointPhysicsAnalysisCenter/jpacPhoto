// Class to contain all relevant kinematic quantities. The kinematics of the reaction
// gamma p -> X p' is entirely determined by specifying the mass of the vector particle.
//
// Additional options to include virtual photon and different baryons (e.g. gamma p -> X Lambda_c) 
// also available
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------


#ifndef KINEMATICS
#define KINEMATICS

#include "constants.hpp"
#include "helicities.hpp"

#include "TMath.h"

#include <array>
#include <vector>
#include <string>
#include <cmath>

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // The reaction kinematics object is intended to have all relevant kinematic quantities
    // forthe reaction. Here you'll find the momenta and energies of all particles, angles,
    // invariant, etc.
    //
    // These are the basis for more complicated structures such as covariants or the
    // amplitudes themselves
    // ---------------------------------------------------------------------------
    
    class reaction_kinematics
    {
        public: 
        // Empty constructor
        // defaults to compton scattering: gamma p -> gamma p
        reaction_kinematics()
        {};

        // Constructor to fully specify the final state
        reaction_kinematics(double mX, double mR = M_PROTON)
        : _mX(mX), _mX2(mX*mX), _mR(mR), _mR2(mR*mR),
          _mB(0.), _mB2(0.), _mT(M_PROTON), _mT2(M2_PROTON),
          _photon(true)
        {};

        // Constructor to fully specify both final and initial states
        reaction_kinematics(double mB, double mT, double mX, double mR)
        : _mX(mX), _mX2(mX*mX), _mR(mR), _mR2(mR*mR),
          _mB(mB), _mB2(mB*mB), _mT(mT), _mT2(mT*mT),
          _photon(!(mB>0.))
        {};

        // ---------------------------------------------------------------------------
        // Masses 

        inline double Wth(){ return (_mX + _mR); }; // square root of the threshold
        inline double sth(){ return Wth() * Wth(); }; // final state threshold

        // Assessor functions for masses
        inline double get_meson_mass() { return _mX; };
        inline double get_recoil_mass(){ return _mR; };
        inline double get_target_mass(){ return _mT; };
        inline double get_beam_mass()  { return _mB; };
        
        // Setters for masses
        inline void set_meson_mass(double x) { _mX = x; _mX2 = x*x; };
        inline void set_target_mass(double x){ _mT = x; _mT2 = x*x; };
        inline void set_recoil_mass(double x){ _mR = x; _mR2 = x*x; };
        inline void set_beam_mass(double x)  { _mB = x; _mB2 = x*x; };

        // Get whether current kinematics has photon beam
        inline bool is_photon() { return _photon;  };
        inline bool is_virtual(){ return _virtual; };
        inline void set_Q2(double x)
        {
            if (!is_photon()) 
            {
                std::cout << "Trying to set Q2 without initializing as a photon first! \n";
                std::cout << "Initialize reaction_kinematics with massless beam to indicate photon then use set_Q2()!" << std::endl;
                exit(1);
            }
            _virtual = true; 
            set_beam_mass(x);
        }

        // ---------------------------------------------------------------------------
        // Quantum numbers of produced meson. 
        
        inline std::array<int,2> get_meson_JP(){ return _mjp; };

        inline void set_meson_JP(int J, int P)
        { 
            _mjp = {J, P};
            _helicities = get_helicities(J, get_baryon_JP()[0], _photon);
            _nAmps = _helicities.size();
        };
        inline void set_meson_JP(std::array<int,2> jp){ set_meson_JP(jp[0], jp[1]); };

        
        // ---------------------------------------------------------------------------
        // Quantum numbers of produced baryon.
        // The baryon spin (and only this quantity) is multiplied by 2 to be saves as an int

        inline std::array<int,2> get_baryon_JP(){ return _bjp; };

        inline void set_baryon_JP(int J, int P)
        { 
            _bjp = {J, P};
            _helicities = get_helicities(get_meson_JP()[0], J, _photon);
            _nAmps = _helicities.size();
        };
        inline void set_baryon_JP(std::array<int,2> jp){ set_baryon_JP(jp[0], jp[1]); };

        // ---------------------------------------------------------------------------
        // Accessing the helicity combinations 
        inline int num_amps(){ return _nAmps; };
        inline std::array<int, 4> helicities(int i)
        {
            if (i < 0 || i >= _nAmps) 
            {
                std::cout << "Error! Can't find helicities with index " << i << "! Returning zeros..." << std::endl;
                return {0, 0, 0, 0};
            }

            return _helicities[i];
        };

        //--------------------------------------------------------------------------
        // Other quantities

        // Moduli of the initial and final state 3-momenta 

        inline double initial_momentum(double s)
        {
            return sqrt( Kallen(s, _mB2, _mT2)) / (2. * sqrt(s));
        };

        inline double final_momentum(double s)
        {
            return sqrt( Kallen(s, _mR2, _mX2)) / (2. * sqrt(s));
        };

        // Energies of all the particles

        inline double beam_energy(double s)
        {
            return (s - _mT2 + _mB2) / sqrt(4. * s);
        };

        inline double target_energy(double s)
        {
            return (s + _mT2 - _mB2) / sqrt(4. * s);
        };

        inline double meson_energy(double s)
        {
            return (s - _mR2 + _mX2) / sqrt(4. * s);
        };

        inline double recoil_energy(double s)
        {
            return (s + _mR2 - _mX2) / sqrt(4. * s);
        };

        // Get s-channel scattering angle from invariants
        inline double z_s(double s, double t)
        {
            std::complex<double> qdotqp = initial_momentum(s) * final_momentum(s);
            std::complex<double> E1E3   = beam_energy(s) * meson_energy(s);

            double result = t - _mX2 - _mB2 + 2.*real(E1E3);
            result /= 2. * real(qdotqp);

            return result;
        };

        // Scattering angle in the s-channel
        // Use TMath::ACos instead of std::acos because its safer at the end points
        inline double theta_s(double s, double t)
        {
            double zs = z_s(s, t);
            if (std::abs(zs - 1) < 1.E-5){zs =  1.;}
            if (std::abs(zs + 1) < 1.E-5){zs = -1.;}
            return TMath::ACos( zs );
        };

        // Invariant variables
        inline double t_man(double s, double theta)
        {
            std::complex<double> qdotqp = initial_momentum(s) * final_momentum(s);
            std::complex<double> E1E3   = beam_energy(s) * meson_energy(s);

            return _mX2 + _mB2 - 2. * real(E1E3) + 2. * real(qdotqp) * cos(theta);
        };

        inline double u_man(double s, double theta)
        { 
            return _mX2 + _mB2 + _mT2 + _mR2 - s - t_man(s, theta);
        };

        // Scattering angles in t and u channel frames
        inline std::complex<double> z_t(double s, double theta)
        {
            double t = t_man(s, theta);
            double u = u_man(s, theta);

            std::complex<double> result;
            result  = t * (s - u) + (_mB2 - _mX2) * (_mT2 - _mR2);
            result /=  sqrt(XR * Kallen(t, _mX2, _mB2)) * sqrt(XR * Kallen(t, _mT2, _mR2));

            return result;
        };

        inline std::complex<double> z_u(double s, double theta)
        {
            double t = t_man(s, theta);
            double u = u_man(s, theta);

            std::complex<double> result;
            result  = u * (t - s) + (_mB2 - _mR2) * (_mT2 - _mX2);
            result /=  sqrt(XR * Kallen(u, _mR2, _mB2)) * sqrt(XR * Kallen(u, _mT2, _mX2));

            return result;
        };

        inline int helicity_index(std::array<int,4> hel){ return find_helicity(hel, _mjp[0], _bjp[0], _photon); };

        // Phase relating lambda_gamma = +1 and lambda_gamma = -1 amplitudes 
        // Depends on the channel with respect to which the helicities are defined
        inline double intrinsic_parity(helicity_channel channel)
        {
            int s_a, s_b, s_c, s_d;
            int eta_a, eta_b, eta_c, eta_d;

            // a is always the photon
            s_a = 2; eta_a = -1; // spin multiplied by two because of spin 1/2 baryons

            switch (channel)
            {
                case helicity_channel::S :
                {
                    s_b =  1;            eta_b = +1;         // proton
                    s_c =  2*_mjp[0];    eta_c = _mjp[1];   // produced meson
                    s_d =  _bjp[0];      eta_d = _bjp[1];   // recoil baryon
                    break;
                }
                case helicity_channel::T :
                {
                    s_b =  2*_mjp[0];   eta_b = _mjp[1];    // produced meson
                    s_c =  1;           eta_c = +1;          // proton
                    s_d =  _bjp[0];     eta_d = _bjp[1];    // recoil baryon
                    break;
                }
                case helicity_channel::U :
                {
                    s_b =  _bjp[0];      eta_b = _bjp[1];    // recoil baryon
                    s_c =  1;            eta_c = +1;          // proton
                    s_d =  2*_mjp[0];    eta_d = _mjp[1];    // produced meson
                    break;
                }

                default: { return 0.; }
            };

            int eta = eta_a * eta_b * eta_c * eta_d * pow(-1., double( (s_c + s_d - s_a - s_b)/2 ));
            
            return double(eta);
        };

        inline double parity_phase(std::array<int, 4> helicities, helicity_channel channel)
        {
            int lam, lamp;
            switch (channel)
            {
                case helicity_channel::S :
                {
                    lam =  (2 * helicities[0] - helicities[1]);
                    lamp = (2 * helicities[2] - helicities[3]);
                    break;
                }
                case helicity_channel::T :
                {
                    lam =  (2 * (helicities[0] - helicities[2]));
                    lamp = (helicities[1] - helicities[3]);
                    break;
                }
                case helicity_channel::U :
                {
                    lam =  (2 * helicities[0] - helicities[3]);
                    lamp = (2 * helicities[2] - helicities[1]);
                    break;
                }

                default: { return 0.; }
            };

            double eta = intrinsic_parity(channel) *  pow(-1., double( (lam - lamp)/2 ));
            return eta;
        };

        inline double parity_phase(int i, helicity_channel channel)
        {
            return parity_phase(_helicities[i], channel);
        };

        inline double H_to_GJ_angle(double s, double t)
        {
            double beta = final_momentum(s) / meson_energy(s);
            double zs = z_s(s, t);

            double cosAlpha = (beta - zs) / (beta * zs - 1.);

            return TMath::ACos( cosAlpha );
        };

        private:

        // Masses are private to prevent them from being changed mid calculation 
        // Instead should be manipulated with public accessor and settor above
        bool _photon = true, _virtual = false;    // whether we have a photon and if its virtual
        double _mB = 0.,       _mB2 = 0.;         // mass and mass squared of the "beam" 
        double _mX = 0.,       _mX2 = 0.;         // mass and mass squared of the produced particle
        double _mT = M_PROTON, _mT2 = M2_PROTON;  // mass of the target, assumed to be proton unless overriden
        double _mR = M_PROTON, _mR2 = M2_PROTON;  // mass of the recoil baryon, assumed to be proton unless overriden

        // Quantum numbers of final state meson and baryon 
        std::array<int,2> _mjp{{1,1}};
        std::array<int,2> _bjp{{1,1}}; // Baryon is multiplied by two (J,P) = (1,1) -> 1/2+

        // Helicity configurations
        // Defaults to spin-1 meson, spin-1/2 baryon
        // Beam [0], Target [1], Produced Meson [2], Recoil Baryon [3]
        int _nAmps = SPIN_ONE_HELICITIES.size(); 
        std::vector< std::array<int, 4> > _helicities = SPIN_ONE_HELICITIES;
    };
};

#endif
