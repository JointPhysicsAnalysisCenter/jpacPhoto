// Class which allows an amplitude to be fit to data based on chi2 minimization
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef FITTER_HPP
#define FITTER_HPP

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "data_set.hpp"
#include "print.hpp"

#include <chrono>
#include <string>
#include <vector>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom.h"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Structs for storing relevant info inside the fitter

    // Each free parameter of a model has associated with is a bunch of options
    class parameter
    {
        public:

        // Default initialized
        // This initializes a fixed parameter for normalizations
        parameter()
        : _fixed(true), _value(1.), 
          _label(default_label(0)), _i(-1),
          _lower(0), _upper(5)
        {};
        
        // Constructor for amplitude variables
        parameter(int i)
        : _i(i), _label(default_label(i))
        {};

        int         _i;
        std::string _label;
        std::string _message;
        bool   _fixed         = false;
        double _value         = 0;
        bool   _custom_limits = false;
        double _upper         = 0;
        double _lower         = 0;
        double _step          = 0.1;
        bool   _positive      = false;

        // If this parameter is synced to be equal to another
        bool   _synced        = false;
        int    _sync_to       = -1; // Which param its synced to

        static inline std::string default_label(int i)
        {
            return "par[" + std::to_string(i) + "]";
        };
    };

    // ---------------------------------------------------------------------------
    // Actual fitter object
    // This is templated because it requires an implementation of the chi2 function to be fit
    // The template F should contain the details of the fit and the following static functions
    // fcn(amplitude, std::vector<data_set>&) [function to be minimized, e.g. chi2]
    // 
    template<class F>
    class fitter
    {
        public: 

        // Basic constructor, only requires amplitude to be fit 
        // uses default settings for minuit
        fitter(amplitude amp_to_fit)
        : _amplitude(amp_to_fit),
          _minuit(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined"))
        {
            reset_parameters();
        };

        // Parameterized constructor 
        // with explicit choice of minimization strategy and tolerance of minuit routines
        fitter(amplitude amp_to_fit, std::string strategy, double tolerance = 1.E-6)
        : _amplitude(amp_to_fit), _tolerance(tolerance),
          _minuit(ROOT::Math::Factory::CreateMinimizer("Minuit2", strategy))
        {
            reset_parameters();
        };

        // -----------------------------------------------------------------------
        // Methods to add data to be fit against

        inline void add_data(data_set data){ _N += data._N; _data.push_back(data); };
        inline void add_data(std::vector<data_set> data){ for (auto datum : data) add_data(datum); };
        inline void clear_data(){ _data.clear(); _N = 0; };

        // -----------------------------------------------------------------------
        // Set limits, labels, and fix parameters

        // Reset labels, limits and options on all parameters
        inline void reset_parameters()
        {
            _pars.clear(); _Nfree = _amplitude->N_pars();
            for (int i = 0; i < _amplitude->N_pars(); i++) _pars.push_back(i);
            set_parameter_labels(_amplitude->parameter_labels());
        };

        // Give each parameter a label beyond their default par[i] name
        inline void set_parameter_labels(std::vector<std::string> labels)
        {
            if (labels.size() != _pars.size())
            {
                warning("fitter::set_parameter_labels", "Labels vector does not match number of parameters!");
                return;
            }
            for (int i = 0; i < _pars.size(); i++) _pars[i]._label = labels[i];
        };

        // Set limits and/or a custom stepsize
        inline void set_parameter_limits(parameter& par, std::array<double,2> bounds, double step = 0.1)
        {
            par._custom_limits = true;
            par._lower         = bounds[0];
            par._upper         = bounds[1];
            par._step          = step;
        };

        inline void set_parameter_limits(std::string label, std::array<double,2> bounds, double step = 0.1)
        {
            int index = find_parameter(label);
            if (index < 0) return;
            return set_parameter_limits(_pars[index], bounds, step);
        };

        inline void make_positive_definite(parameter& par, bool x = true){ par._positive = x; };
        inline void make_positive_definite(std::string label, bool x = true)
        {
            int index = find_parameter(label);
            if (index < 0) return;
            make_positive_definite(_pars[index], x);
        };

        inline void fix_parameter(parameter& par, double val)
        {
            // If parameter is already fixed, just update the fixed val
            // otherwise flip the fixed flag and update the number of free pars
            if (!par._fixed) _Nfree--;
            par._fixed = true;
            par._value = val;
            _fit = false;
        };

        inline void fix_parameter(std::string label, double val)
        {
            int index = find_parameter(label);
            if (index < 0) return;
            return fix_parameter(_pars[index], val);
        };

        inline void free_parameter(parameter& par)
        {
            // if not fixed, this does nothing
            if (!par._fixed) return;
            par._fixed = false;
            _fit = false;
            _Nfree++;
        };

        inline void free_parameter(std::string label)
        {
            int index = find_parameter(label);
            if (index < 0) return;
            free_parameter(_pars[index]);
        };

        inline void sync_parameter(std::string par, std::string synced_to)
        {
            int i            = find_parameter(par);
            int i_sync_to    = find_parameter(synced_to);
            _pars[i]._synced  = true;
            _pars[i]._sync_to = i_sync_to;

            // Also save present value in case there is one
            _pars[i]._value = _pars[i_sync_to]._value;
            if (!_pars[i_sync_to]._fixed)_Nfree--;
        };

        inline void unsync_parameter(std::string par){ _pars[ find_parameter(par) ]._synced = false; };

        // -----------------------------------------------------------------------
        // Methods related to fit options

        // Set the maximum number of calls minuit will do
        inline void set_max_calls(int n){ print(n); _max_calls = n; };
        
        // Message level for minuit (0-4)
        inline void set_print_level(int n){ _print_level = n; };

        // Change tolerance
        inline void set_tolerance(double tol){ _tolerance = tol; };

        // Change the guess range for initializing parameters
        inline void set_guess_range(std::array<double,2> bounds){ _guess_range = bounds; };

        // Actually do the fit given a vector of size amp->N_pars() as starting values
        // Prints results to command line but also returns the best-fit chi2 value
        inline void do_fit(std::vector<double> starting_guess, bool show_data = true)
        {
            if (starting_guess.size() != _Nfree) 
            {
                warning("fitter::do_fit", "Starting guess not the correct size! Expected " + std::to_string(_Nfree) + " parameters!");
                return;
            };

            set_up(starting_guess);

            if (show_data) { line(); data_info(); };
            parameter_info(starting_guess);

            auto start = std::chrono::high_resolution_clock::now();
            std::cout << "Beginning fit..." << std::flush; 

            if (_print_level != 0) line();   
            _minuit->Minimize();
            if (_print_level != 0) line();   

            std::cout << "Done! \n";

            // Timing info
            auto stop     = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast< std::chrono::seconds>(stop - start);
            std::cout << std::left << "Finished in " << duration.count() << " s" << std::endl;

            line();
            print_results();
        };

        // Same as above (do single fit) but initialize random parameters
        inline void do_fit()
        {            
            // Initial guess
            std::vector<double> guess;  

            // Initialize the guess for each parameter
            std::vector<int> synced_pars;
            for (auto par : _pars)
            {
                if (par._fixed) continue;
                if (par._synced){ synced_pars.push_back(par._i); continue; } // Flag any synced parameters, we'll come back to them

                if (par._custom_limits) guess.push_back(_guesser->Uniform(par._lower, par._upper)); 
                else if (par._positive) guess.push_back(_guesser->Uniform(0., _guess_range[1]));
                else                    guess.push_back(_guesser->Uniform(_guess_range[0], _guess_range[1]));
                par._value = guess.back();
            };  

            // After all randomized parameters have been set, rego through the synced ones
            for (int i : synced_pars)
            {
                int j = _pars[i]._sync_to;
                _pars[i]._value = _pars[j]._value;
            };

            do_fit(guess);
        };

        // Do N fits with random parameters and return the best fit found
        inline void do_fit(int N)
        {
            divider();
            std::cout << std::left << "Commencing N = " + std::to_string(N) + " fit iterations." << std::endl;

            // Initial guess
            std::vector<double> guess;

            for (int i = 1; i <= N; i++)
            {
                guess.clear();
                // Whether this is the first iteration
                bool first_fit =  !(i-1);
                // Initialize the guess for each parameter
                
                std::vector<int> synced_pars;
                for (auto par : _pars)
                {
                    if (par._fixed) continue;
                    if (par._synced){ synced_pars.push_back(par._i); continue; } // Flag any synced parameters, we'll come back to them

                    if (par._custom_limits) guess.push_back(_guesser->Uniform(par._lower, par._upper)); 
                    else if (par._positive) guess.push_back(_guesser->Uniform(0., _guess_range[1]));
                    else                    guess.push_back(_guesser->Uniform(_guess_range[0], _guess_range[1]));
                    par._value = guess.back();
                };

                // After all randomized parameters have been set, rego through the synced ones
                for (int i : synced_pars)
                {
                    int j = _pars[i]._sync_to;
                    _pars[i]._value = _pars[j]._value;
                };

                // Do out fit with this random guess
                if (!first_fit) std::cout << std::left << "Fit (" + std::to_string(i) + "/" + std::to_string(N) + ")" << std::endl;
                do_fit(guess, first_fit);


                // Compare with previous best and update
                if ( first_fit || fcn() < _best_fcn)
                {
                    _best_fcn_dof = _minuit->MinValue() / (_N - _minuit->NFree());
                    _best_fcn     = _minuit->MinValue();
                    _best_pars    = convert(_minuit->X());
                    _best_errs    = convert(_minuit->Errors());
                };
            };
            
            // After looping, set the best_pars to the amplitude
            _amplitude->set_parameters(_best_pars);

            // And set the global saved pars to the best_fit
            std::cout << std::left << "Best fit found after N = " + std::to_string(N) + " iterations" << std::endl;
            line();
            print_results(false);
        };

        // // Repeat do_fit N times and find the best fit
        // // Parameters are randomly initialized each time on the interval [-5, 5] unless custom limits are set
        // void do_fit(int N);

        // Return a vector of best-fit parameters from last fit
        inline std::vector<double> pars(){ return (_fit) ? _fit_pars : std::vector<double>(); };

        // Return value of fit function from last fit
        inline double fcn()     { return (_fit) ? _fcn     : NaN<double>(); };
        inline double fcn_dof() { return (_fit) ? _fcn_dof : NaN<double>(); };

        // Return the pointer to the amp being fit at its current state
        inline amplitude get_amplitude(){ return _amplitude; };

        private:

        // This ptr should point to the amplitude to be fit
        amplitude _amplitude = nullptr;

        // -----------------------------------------------------------------------
        // Data handling

        int _N = 0;  // Total number of data points
        std::vector<data_set> _data;  // Contain all data

        // -----------------------------------------------------------------------
        // MINUIT handling 

        int _print_level   = 0;     // Error code for MINUIT
        int _max_calls     = 1E6;   // Max calls allowed for minimization fcn
        double _tolerance  = 1.E-6; // Minimization tolerance

        ROOT::Math::Minimizer * _minuit;
        ROOT::Math::Functor _wfcn;

        // Random number generator for creating initial guesses;
        TRandom *_guesser = new TRandom(0);

        // Initialize minuit with all our parameter options etc
        inline void set_up(std::vector<double> starting_guess)
        {
            _minuit->Clear();
            _minuit->SetTolerance(_tolerance);
            _minuit->SetPrintLevel(_print_level);
            _minuit->SetMaxFunctionCalls(_max_calls);

            // Iterate over each _par but also keep track of the index in starting_guess 
            // because parameters might be fixed, these indexes dont necessarily line up
            int i = 0;
            for (auto par : _pars)
            {   
                if (par._fixed || par._synced) continue;
                _minuit->SetVariable(i, par._label, starting_guess[i], par._step);
                if (par._custom_limits) _minuit->SetVariableLimits(i, par._lower, par._upper);
                if (par._positive)      _minuit->SetVariableLowerLimit(i, 0.);
                i++; // move index up
            };
        
            _wfcn = ROOT::Math::Functor(this, &fitter::fit_fcn, _Nfree);
            _minuit->SetFunction(_wfcn);
        };

        // -----------------------------------------------------------------------
        // Calcualtions of chi-squared 
        
        // This is the actual function that gets called by minuit
        inline double fit_fcn(const double *cpars)
        { 
            // First convert the C string to a C++ vector
            std::vector<double> pars = convert(cpars);

            // Pass parameters to the amplitude
            _amplitude->set_parameters(pars);

            // Pass both this and data to fit function
            return F::fcn(_data, _amplitude); 
        };

        // Save of the last fit run
        bool _fit = false;      // Whether a fit has already been done or not yet
        double _fcn, _fcn_dof;  // Last saved value of fcn function
        std::vector<double> _fit_pars, _errors;

        // Save of the best fits found if running multiple times
        double _best_fcn, _best_fcn_dof;
        std::vector<double> _best_pars, _best_errs;

        // -----------------------------------------------------------------------
        // Parameter handling

        // Store of parameter info
        std::vector<parameter> _pars;

        // Number of free parameters
        int _Nfree = 0;

        // Default guess_range to initalize parameters
        std::array<double,2> _guess_range = {-5, 5};

        // Given a parameter label, find the corresponding index
        inline int find_parameter(std::string label)
        {
            for (auto par : _pars) if (par._label == label) return par._i;
            return error("fitter::find_parameter", "Cannot find parameter labeled " + label + "!", -1);
        };

        // Given a C-style array of size _Nfree
        // Convert to a C++ style std::vector and populate
        // fixed value parameters in the expected order
        inline std::vector<double> convert(const double * cpars)
        {
            std::vector<double> result;

            // Move along the pars index when a parameter is not fixed
            int i = 0;

            std::vector<int> synced_pars;
            for (auto par : _pars)
            {
                if (par._synced) synced_pars.push_back(par._i); 
                if (par._fixed || par._synced) result.push_back(par._value);
                else { result.push_back(cpars[i]); i++; };
            };

            for (int j : synced_pars)
            {
                _pars[j]._value = _pars[_pars[j]._sync_to]._value;
                result[j] = result[_pars[j]._sync_to];
            };

            if (i != _Nfree) warning("fitter::convert", "Something went wrong in converting parameter vector.");
            return result;
        };

        // -----------------------------------------------------------------------
        // Methods to print out status to command line

        // Summary of data sets that have been recieved
        inline void data_info()
        {
            using std::cout; using std::left; using std::endl; using std::setw;
            
            if (_data.size() == 0)
            {
                warning("fitter::data_info", "No data found!"); 
                return;
            };

            cout << left;
            divider();
            cout << "Fitting amplitude (\"" << _amplitude->id() << "\") to " << _N << " data points:" << endl;
            line();
            cout << setw(25) << "DATA SET"         << setw(30) << "TYPE     "      << setw(10) << "POINTS" << endl;
            cout << setw(25) << "----------------" << setw(30) << "--------------" << setw(10) << "-------" << endl;
            for (auto data : _data)
            {
                cout << setw(25) << data._id  << setw(30)  << F::data_type(data._type)  << setw(10) << data._N << endl;  
            };
        };

        // Similar summary for parameters
        // Display alongside a vector of current parameter values
        // bool start is whether this is the starting guess vector or the 
        // best fit results
        inline void parameter_info(std::vector<double> starting_guess)
        {
            using std::cout; using std::left; using std::endl; using std::setw;

            cout << std::setprecision(8);
            cout << left;

            line(); divider();
            // Print message at the beginning of the fit
            cout << "Fitting " + std::to_string(_Nfree) << " (of " << std::to_string(_pars.size()) << ") parameters" << endl;
            line();

            cout << left << setw(10) << "i"     << setw(17) << "PARAMETER"  << setw(20) << "START VALUE"  << endl;
            cout << left << setw(10) << "-----" << setw(17) << "----------" << setw(20) << "------------" << endl;

            // Moving index from the guess vector
            int i = 0;
            
            std::vector<int> synced_pars;
            std::vector<double> vals;
            for (auto &par : _pars)
            {
                if (par._synced)
                {
                    vals.push_back(0);
                    synced_pars.push_back(par._i);
                    par._message = "[= " + std::to_string(par._sync_to) + "]";
                    continue;
                };

                // Or is fixed
                if (par._fixed)
                {
                    vals.push_back(par._value);
                    par._message  =  "[FIXED]";
                    continue;
                }

                // if posdef
                if (par._positive) par._message  =  "[>= 0]";

                if (par._custom_limits)
                {   
                    std::stringstream ss;
                    ss << std::setprecision(5) << "[" << par._lower << ", " << par._upper << "]";
                    par._message = ss.str();
                };

                par._value = starting_guess[i];
                vals.push_back(par._value);
                i++;
            };

            for (int j : synced_pars)
            {
                int k = _pars[j]._sync_to;
                _pars[j]._value = _pars[k]._value;
                vals[j] = _pars[k]._value;
            };

            for (auto par : _pars)
            {
                cout << left << setw(10) << par._i << setw(17) << par._label << setw(20) << vals[par._i] << setw(20) << par._message << endl;
            };

            line(); divider(); line();
        };

        // After a fit return a summary of fit results
        // At the end of a fit, print out a table sumarizing the fit results
        // if last_fit == true, we grab the results from the most recent fit in _minuit
        // else we print out the ones saved in _best_fit
        inline void print_results(bool last_fit = true)
        {
            using std::cout; using std::left; using std::endl; using std::setw;

            cout << std::setprecision(8);
            cout << left;

            int dof                  = _N - _minuit->NFree();
            double fcn               = (last_fit) ? _minuit->MinValue()               : _best_fcn;
            double fcn_dof           = (last_fit) ? _minuit->MinValue() / double(dof) : _best_fcn_dof;
            std::vector<double> pars = (last_fit) ? convert(_minuit->X())             : _best_pars;
            std::vector<double> errs = (last_fit) ? convert(_minuit->Errors())        : _best_errs;


            divider();
            std::cout << std::left << std::setw(5) << "fcn = "       << std::setw(15) << fcn     << std::setw(5) << "";
            std::cout << std::left << std::setw(10) << "fcn/dof = "  << std::setw(15) << fcn_dof << "\n";

            line();

            cout << left << setw(10) << "i"     << setw(16) << "PARAMETER"  << setw(18) << "FIT VALUE"    << setw(18) << "ERROR"        << endl;
            cout << left << setw(10) << "-----" << setw(16) << "----------" << setw(18) << "------------" << setw(18) << "------------" << endl;

            for (int i = 0; i < pars.size(); i++) _pars[i]._value = pars[i];

            for (auto par : _pars)
            {
                double val;

                std::stringstream ss;
                ss << std::setprecision(8) << errs[par._i];
                std::string err = !(par._fixed || par._synced) ? ss.str() : par._message;

                cout << left << setw(10) << par._i << setw(16) << par._label << setw(18) << par._value << setw(18) << err << endl;
            };
            line(); divider(); line();
            
            // At the end update the amplitude parameters to include the fit results
            _amplitude->set_parameters(pars);

            _fcn      = fcn;
            _fcn_dof  = fcn_dof;
            _fit_pars = pars;

            // Let rest of the fitter that a fit result has been saved
            _fit = true;
        };
    };
};

#endif