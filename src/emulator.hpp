#pragma once

using namespace std;

#include "model.hpp"
#include "mpi.hpp"

class Emulator {
public:
	Emulator(const Model &model, int seed);
	void inference();
	
private:	
	double prior_dist_min;                         // Controls separation of points 

	vector <double> peak;                          // The current position of the peak
	
	unsigned int nem_var;
	vector <unsigned int> em_var;                  // Stores variables to be emulated
	
	unsigned int nvar;                             // Number of variables
	vector <Variable> var;                         // Stores all the variables in the system

	unsigned int nvar_param_list;
	vector <unsigned int> var_param_list;          // List of variables which are parameters

	vector <double> param_value;                   // The model parameter values
	vector <IndValue> ind_value;                   // Stores individual properties on the chain

	vector <GroupEvent> gr_ev;                     // Stores information about observed transition events 

	vector <unsigned int> param_var_ref;           // References variable from parameter
	vector < vector <unsigned int> > ind_effect_var_ref; // References variable from individual effect

	unsigned int nfe;
	vector <FE> fe;                                // Fixed effect in the model
	
	// These are variable use for fast evalulation of the emulator
	vector < vector <double> > Z, ZT;  
	vector < vector <double> > Dinv;	
	vector < vector <double> > Ainv;
	vector < vector <double> > Ainv_Z;
	vector <double> betahat;
	vector <double> q;                             // This stores Ainv*(g-Zbetahat)
	vector <double> g;
	double gmax; 
	double sigma_sq_fac;
	
	bool nugget;                                   // Set to true if nugget is on
	
	unsigned int nhyper;
	vector <double> hyper;                         // Correlation lengths (plus nugget at end)
		
	vector <DataSample> data_sample;               // Stores "data" samples regarding inference

	double laplace_approximation();
	double maximise_posterior(vector < vector <double> > &H);
	double maximise_posterior_BFGS(vector < vector <double> > &H);
	vector <double> get_value_from_param();
	void update_param_from_value(const vector <double> &value);
	void initialise_variable();
	vector < vector <double> > latin_hypercube_sample(unsigned int nsamp) const;
	vector < vector <double> > prior_sample(unsigned int nsamp, bool spacing) const;
	void add_data_samp(vector < vector <double> > &samp);
	void tune_emulator();
	double marginal();
	void set_Z_g();
	vector <double> finite_difference_grid_size() const;
	vector <double> calculate_gradient_hyper(vector <double> &hyper);
	vector < vector <double> > calculate_hessian_hyper(vector <double> &hyper);
	vector < vector <double> > calculate_A() const;
	EmulatorEstimate evalulate_emulator(const vector <double> &em_var_value, bool get_sd) const;
	vector < vector <double> > emulate_new_points(unsigned int g);
	double prior_dist(const vector <double> &value1, const vector <double> &value2) const;
	vector < vector <double> > mcmc(unsigned int nsamp, unsigned int nsamp_var, unsigned int thin, bool use_var, double phi);
	vector < vector <double> > MVN_param_sample(unsigned int nsamp);
	vector < vector <double> > sample_no_emulator(unsigned int nsamp);
	void trace_plot(const vector < vector <double> > &samp);
	void find_peak_cull_points();
	
	/// Used for checking
	double calculate_posterior(const vector <double> &value);
	void check_grad_H();
	vector < vector <double> > calculate_hessian_check(const vector <double> &mean, unsigned int var_max);
	vector <double> calculate_gradient_check(vector <double> value, const vector <unsigned int> &var_list);
	vector <double> finite_different_grid_size(const vector <double> &value);
	void simple();
	void show_fit();
	void minimum_prior_dist() const;
	void check_2d();
	void check_prediction_accuracy();
	void check_one_ind_effect();
	void check_scan();
	void phase_transition();
	
	const Model &model;
	Mpi mpi;
};
