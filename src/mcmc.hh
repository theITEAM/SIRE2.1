#ifndef SIRE__MCMC_HH
#define SIRE__MCMC_HH

using namespace std;

#include "model.hh"

class MCMC
{
public:
  bool burnin;                                   // Set to true if in the burnin phase 
  
  vector <double> param_value;                   // The model parameter values
  
  vector <Jump> param_jump;                      // Information about univariate proposals
  
  vector <Jump> ind_effect_jump;                 // Information about individual effect proposals
  
  vector <Jump> event_jump;                      // Information about event time proposals
  
  vector <Jump> add_rem_jump;                    // Information about adding / removing infected individual proposals
  
  vector <IndValue> ind_value;                   // Stores individual properties on the chain
  
  vector <double> L_ind_effect;                  // The likelihood for the individual effects (for each covariance matrix)

  vector <double> L_inf_events;                  // The likelihood for the infection events (for each group)
  
  double L_trans_events;                         // The likelihood for the transition events
 
	vector <double> L_diag_test;                   // The likelihood for disease diagnostic test
  
  double prior;                                  // The prior probability
  
  vector <IndPM> ind_effect_posterior_mean;      // Calculate posterior mean for individual effects
  int nind_effect_posterior_mean;
  
  vector <Sample> sample;                        // Stores samples from the chain (to generate diagnostics)
  
  ofstream trace;                                // Used to output the trace plots
  
  MCMC(const Model &model);
  void set_likelihoods();
  void update();
  void trace_initialise(const string file);
  void trace_output(const int s);
  void output_statistics(const string file) const;
  void store_sample();
  void ind_effect_posterior_mean_update();
	void initialise_inf_sampler(int s);
  void temp();                                  // In check.cc

private:
  Statistics get_statastic(const vector <double> &vec) const;
  void initialise_proposals();
  void check_chain();                           // In check.cc
  
  const Model &model;
};

#endif