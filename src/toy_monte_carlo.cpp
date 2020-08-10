// Sample event-generator for the gamma p -> R -> Jpsi p -> l+ l- p 
// reaction. Events are generated on a flat phase space and weighted by the probabilty distribution
// of a supplied amplitude.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "toy_monte_carlo.hpp"

// ---------------------------------------------------------------------------
// main driver loop to generate N events 
void jpacPhoto::toy_monte_carlo::generate(double _beam_energy, int N)
{
    // TODO: can add non uniform beam profile here
    beam_energy = _beam_energy;
    W = W_cm(_beam_energy);
    s = W * W;

    for (int n = 0; n < N; n++)
    {
        generate_event();
        generate_weight();
    }

    kin->Write();

    if (dyn != NULL) dyn->Write();

    outfile->Close();
};

// ---------------------------------------------------------------------------
// Set up the tree structure of the output file
void jpacPhoto::toy_monte_carlo::set_up_kin()
{
    // tree for kinematics quantities
    kin = new TTree("kinematics", "kinematics");
    
    // Lepton pair 4-momentum components
    kin->Branch("ep_px",    &ep_px);
    kin->Branch("ep_py",    &ep_py);
    kin->Branch("ep_pz",    &ep_pz);
    kin->Branch("ep_E",     &ep_E);
    kin->Branch("em_px",    &em_px);
    kin->Branch("em_py",    &em_py);
    kin->Branch("em_pz",    &em_pz);
    kin->Branch("em_E",     &em_E);

    // Recoil proton 4-momenta
    kin->Branch("prec_px",  &prec_px);
    kin->Branch("prec_py",  &prec_py);
    kin->Branch("prec_pz",  &prec_pz);
    kin->Branch("prec_E",   &prec_E);

    // Angles
    kin->Branch("phi_psi",  &phi_psi);
    kin->Branch("theta_psi",&theta_psi);
    kin->Branch("phi_ep",   &phi_ep);
    kin->Branch("theta_ep", &theta_ep);
};

void jpacPhoto::toy_monte_carlo::set_up_dyn()
{
    // A second tree for dynamical quantities
    // Only need to have this if we have an amplitude supplied
    dyn = new TTree("dynamics", "dynamics");
    
    // Invariants
    dyn->Branch("s", &s);
    dyn->Branch("t", &t);

    // Helicities
    dyn->Branch("lambda_gamma", &lam_gamma);
    dyn->Branch("lambda_ptarg", &lam_ptarg);
    dyn->Branch("lambda_prec",  &lam_ptarg);
    dyn->Branch("lambda_erel",  &lam_erel);

    // Amplitude
    dyn->Branch("weight", &weight);
};

// ---------------------------------------------------------------------------
// Generate an event
void jpacPhoto::toy_monte_carlo::generate_event()
{
    // Generate a random set of angles
    // Make a flat phasespace by sampling costheta not theta directly
    phi_psi     = random(0., 2. * M_PI);
    theta_psi   = acos(random(-1., 1.));
    phi_ep      = random(0., 2. * M_PI);
    theta_ep    = acos(random(-1., 1.));

    // Initial leption momenta in the jpsi decay frame
    TLorentzVector p_psi(0., 0., 0., mJpsi);
    TLorentzVector p_ep(0., 0.,  mJpsi/2., mJpsi/2.);
    TLorentzVector p_em(0., 0., -mJpsi/2., mJpsi/2.);

    // Orient in the right direction before boosting
    p_ep.RotateY(theta_ep);
    p_ep.RotateZ(phi_ep);
    p_em.RotateY(theta_ep);
    p_em.RotateZ(phi_ep);

    // Boost to the psi-proton rest frame
    TVector3 boost_psi(0., 0., sqrt(Kallen(s, mJpsi2, mPro2)) / (s + mJpsi2 - mPro2));
    p_psi.Boost(boost_psi);
    p_ep.Boost(boost_psi);
    p_em.Boost(boost_psi);

    // Rotate by the psi angle
    p_psi.RotateY(theta_psi);
    p_psi.RotateZ(phi_psi);
    p_ep.RotateY(theta_psi);
    p_ep.RotateZ(phi_psi);
    p_em.RotateY(theta_psi);
    p_em.RotateZ(phi_psi);

    TLorentzVector p_prec(-p_psi.X(), -p_psi.Y(), -p_psi.Z(), W - p_psi.E());
    TLorentzVector p_ptarg(0., 0., (mPro2 - s)/ (2. * W), (s + mPro2)/ (2. * W));
    TLorentzVector p_gamma(0., 0., (s - mPro2)/ (2. * W), (s - mPro2)/ (2. * W));

    // Boost everything in the lab frame
    TVector3 boost_lab(0., 0., (s - mPro2) / (s + mPro2));
    p_ep.Boost(boost_lab);
    p_em.Boost(boost_lab);
    p_psi.Boost(boost_lab);
    p_ptarg.Boost(boost_lab);
    p_gamma.Boost(boost_lab);
    p_prec.Boost(boost_lab);
    
    // Save all the values
    ep_px = p_ep.X();  ep_py = p_ep.Y();  ep_pz = p_ep.Z();  ep_E = p_ep.E();
    em_px = p_em.X();  em_py = p_em.Y();  em_pz = p_em.Z();  em_E = p_em.E();
    prec_px = p_prec.X();  prec_py = p_prec.Y();  prec_pz = p_prec.Z();  prec_E = p_prec.E();
    pgamma_px = p_gamma.X();  pgamma_py = p_gamma.Y();  pgamma_pz = p_gamma.Z();  pgamma_E = p_gamma.E();
    ptarg_px = p_ptarg.X();  ptarg_py = p_ptarg.Y();  ptarg_pz = p_ptarg.Z();  ptarg_E = p_ptarg.E();

    kin->Fill();
};

// ---------------------------------------------------------------------------
// Generate the weights from the amplitude supplied
void jpacPhoto::toy_monte_carlo::generate_weight()
{
    if (dyn == NULL) 
    {
        if (error_already_triggered == false)
        {
            std::cout << "mc: No weighting amplitude specified. Skipping weight step... \n";
            error_already_triggered = true;
        }
        return;
    }

    // Generate random set of helicites
    lam_gamma   = 1 - 2 * rnd->Integer(2);
    lam_ptarg   = 1 - 2 * rnd->Integer(2);
    lam_erel    = 1 - 2 * rnd->Integer(2);
    lam_prec    = 1 - 2 * rnd->Integer(2);

    // Momentum transfer squared
    t = amp->kinematics->t_man(s, theta_psi);

    // We sum over the intermediate jpsi helicities
    int lam_psi[3] = {-1, 0, 1};
    std::complex<double> amp_sum = 0.;
    for (int i = 0; i < 3; i++)
    {
        amp_sum  = amp->helicity_amplitude({lam_gamma, lam_ptarg, lam_psi[i], lam_prec}, s, t);
        amp_sum *= wigner_d_int(1, lam_psi[i], lam_erel, theta_psi);
        amp_sum *= exp(- xi * double(lam_psi[i]) * phi_psi);
    }

    weight  = real( amp_sum * conj(amp_sum) );

    dyn->Fill();
};