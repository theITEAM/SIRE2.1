#pragma once

using namespace std;

#include "model.hpp"
#include "mpi.hpp"

class MCMC {
public:
	bool burnin;                                   // Set to true if in the burnin phase

	Anneal anneal;                                 // Stores the inverse temperatures
	
	vector <double> param_value;                   // The model parameter values

	vector <Jump> param_jump;                      // Information about univariate proposals

	vector <Jump> ind_effect_jump;                 // Information about individual effect proposals

	vector <Jump> var_ie_joint_jump;               // Information about ijoint ndividual effect variance proposals

	vector <Jump> event_jump;                      // Information about event time proposals

	vector <Jump> add_rem_jump;                    // Information about adding / removing infected individual proposals

	vector <Jump> mean_time_jump;                  // Proposals which change the mean time of transitions
	
	vector <IndValue> ind_value;                   // Stores individual properties on the chain

	vector <double> L_ind_effect;                  // The likelihood for the individual effects (for each covariance matrix)

	vector <double> L_inf_events;                  // The likelihood for the infection events (for each group)

	double L_trans_events;                         // The likelihood for the transition events

	vector <double> L_diag_test;                   // The likelihood for disease diagnostic test

	vector <double> L_cloglog;                     // The likelihood using a cloglog model

	double prior;                                  // The prior probability

	vector <IndPM> ind_effect_posterior_mean;      // Calculate posterior mean for individual effects
	int nind_effect_posterior_mean;

	vector <Sample> sample;                        // Stores samples from the chain (to generate diagnostics)

	double phi, phi_run;                           // The annealing temperature
	
	vector <double> L_samp;
	vector <double> burnin_Li;                     // Stores likelihood samples from the chain
	vector <double> burnin_phi;                    // Stores likelihood samples from the chain

	ofstream trace;                                // Used to output the trace plots

	MCMC(const Model &model);
	void set_likelihoods();
	void diff(unsigned int num);//zz
	void update();
	void trace_initialise(const string file);
	void trace_output(const int s);
	void output_combined_trace(const vector < vector <Sample> > &sample) const;
	void output_statistics();
	void store_sample();
	void store_burnin_Li();
	void ind_effect_posterior_mean_update();
	void initialise_inf_sampler(int s);
	void set_anneal(unsigned int s);
	void temp();                                  // In check.cc
	void check_annealing();
	//void optimum_quench_schedule();
	void check_chain(unsigned int num);           // In check.cc
	Statistics get_statistic(const vector <double> &vec) const;
	void update_cloglog();                        // In model_comp_log_log
	void set_likelihoods_cloglog();               // In model_comp_log_log
	void cloglog_check(unsigned int num); 
	void cloglog_trace_output(const int s);

private:
	
	void initialise_proposals();
	//double get_phi_from_schedule(unsigned int s) const;
	string stat(double num);
	
	Mpi mpi;
	const Model &model;
};

