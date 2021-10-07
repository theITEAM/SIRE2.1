/// This provides all the data structures used in the code

#ifndef SIRE__STRUCT_HH
#define SIRE__STRUCT_HH

struct Param{
  string name;                    // The name of the parameter
  ParamType type;                 // The type of parameter 
    
  PriorType prior_type;           // The type of prior

  int prior_sd_param;             // References a standard deviation parameter
  
  double prior_val1;              // Values used to specify the prior
  double prior_val2;
};

struct Comp{
  string name;                    // The name of the compartment
  
  double infectivity;             // Sets the relative infectivity of the compartment
  
  vector <int> ind_effect;        // References any individual effects acting on infectivity
  vector <int> fixed_effect;      // References any fixed effects acting on infectivity
	vector <int> snp_effect;        // References any snp effects acting on infectivity
	
	bool all_ind_effect_proposal;   // Set to true if using  proposals on all individual effects on compartment
};

struct Trans{
  string name;                    // The name of the transition

  int from;                       // From and to which compartment
  int to;

  TransType type;                 // The type of the transition 
  
  int mean_param;                 // Reference for the mean parameter
  int shape_param;                // Reference for the shape parameter
  
  vector <int> ind_effect;        // References any individual effects acting on transition rate
  vector <int> fixed_effect;      // References any fixed effects acting on transition rate
  vector <int> snp_effect;        // References any snp effects acting on transition rate
	 
  int data_column;                // The data column giving the time of the transition (if specified)
	
	bool all_ind_effect_proposal;   // Set to true if using  proposals on all individual effects on transition
};

struct TimeRange{
  double tmin;                    // Specifies a time range
  double tmax;
};

struct DiagnosticTest{
	string name;                    // Name of diagnostic test 
	vector <bool> comp;             // Set to true if the compartment 
	int Se_param;                   // The sensitivity parameter
	int Sp_param;                   // The specificity parameter
	
	int data_column;                // The data column giving the diagnostic test information
};

struct EvChange{
  int ind;                        // Specifies which individual must have transition times adaptively changed 
  int trans;                      // The transition of that individual
  EvSampType ev_samp_type;        // The type of event sampling to be used (either gamma/normal or uniform)
};

struct Group{
  string name;                    // The name of the group
  int index;                      // Sets the group infex number
  vector <int> ind_ref;           // The individuals contained in the group  
  int nind;                       // The number of individuals
  vector <EvChange> ev_change;    // Stores which events need to be change on which individuals 
  int nev_change;
  vector <int> unknown;           // A list of individuals with unknowb disease status
  int nunknown;
  TimeRange inference_range;      // The inference time range
  TimeRange observation_range;    // The observation time range
};

struct GroupEffect{
  bool on;                        // Set to true if group effect is on
  int sigma_param;                // Reference for the standard deviation parameter 
  vector <int> param;             // References the group effect parameters
};

struct Table{ 
  unsigned int ncol;              // The number of columns
  unsigned int nrow;              // The number of rows
  vector <vector <string> > ele;  // The elements of the table
};

struct DiagTestResult{
	bool positive;                  // Set to true if diagnistic test result is positive
	double time;                    // The time of the result
};

struct InfSampler{
	bool on;                        // Set to true of infetion sampler is used
	double tmin;                    // The minimum time for the sampler
	double tmax;                    // The maximum time for the sampler
	vector <int> bin;               // The number of time infection is in time range
	vector <double> log_prob;       // The log of the probability of infection occuring in a cetain time range 
	vector <double> prob_sum;       // The sum of the probability
	
	double sample_time();           // Samples a time from the sampler
	double sample_prob(double t);   // The log of the probability for a given time
};

struct Individual{
  string id;                      // The unique identified for an individual
  bool inside_group;              // This is set to true if the individual is inside a group
  int group;                      // The group to which the individual belongs
  int initial_comp;               // The initial compartment for the individual
  Status status;                  // Determines the infection status of the individual
  vector <TimeRange> trans_time_range;      // This stores the potential time range for a given transition 
  vector <double> ind_effect_value;         // Stores actual individual effect values (e.g breeding values) 
  vector <double> fixed_effect_X;           // The design matrix multiplying the fixed effects
  vector <double> fixed_effect_X_unshifted; // The unshifted design matrix 
	vector <SNP> SNP_genotype;                // The genotype for different SNPs
	vector <vector <DiagTestResult> > diag_test_result; // Series of diagnostic test results for an individual 
};

struct IndValue{
  bool infected;                  // Sets if the individual on the MCMC chain is infected
  int index;                      // Sets the individual number
  vector <double> trans_time;     // Sets the transition times for the individual
  vector <double> ind_effect;     // The values for the individual effects
                                  // These quantities are derived from other model parameters:
  double susceptibility;          // Individual susceptibiliy
  vector <double> infectivity;    // Individual infectivity (for each compartment)
  vector <double> trans_infectivity_change;    // Change in individual infectivity going down a transition
  vector <double> trans_mean;                  // The means of transitions
	
	InfSampler inf_sampler;         // An infection sampler (used if particle undergoes add / rem infection)        
};

struct Element{
  int i;
  double val;
};

struct Matrix{
  string name;                                 // The name of the matrix
  vector <double> Ainvdiag;                    // The diagonal members of the inverse matrix
  vector < vector <Element> > Ainvlist;        // Stores the inverse matrix in sparse format
  vector < vector <Element> > Ainvlist2;       // Stores only when i < j
  vector < vector <double> > A;                // Stores the full matrix
  vector < vector <double> > cholesky_matrix;  // Stores the Cholesky decomposition of A
};

struct Covariance{
  int matrix;                     // Reference to the relationship matrix
  int E;                          // The number of individual effects
  vector <int> ind_effect_ref;    // Stores which individual effects the matrix refers to 
  vector <int> var_param;         // References variance parameters
  vector < vector <int> > cor_param;  // References correlation parameters
  vector <int> param_list;        // Lists all parameters used in the covariance matrix
};

struct IndEffect{
  string name;                    // The name of individual effect
  int covar_ref;                  // The covariance the individual effect appears in
  int covar_num;                  // The position in the covariace matrix the individual effect appears   
  vector <int> susceptibility;    // Stores if changes susceptibility
  vector <int> infectivity_comp;  // Stores compartments it changes the infectivity of 
  vector <int> trans_mean;        // Stores which transitions it changes the mean of
	
	bool single_ind_effect_proposal;// Set to true if individual effect updated seperately
};

struct FixedEffect{
  string name;                    // The name of the fixed effect
  int param;                      // The parameter of the fixed effect
  bool L_inf_update;              // The likelihood for infection events must be updated
  bool L_trans_update;            // The likelihood for transition events must be updated
};

struct SNPEffect{
  string name;                    // The name of the SNP effect
  int param_mag;                  // The parameter of the SNP effect magnitude
	int param_dom;                  // The parameter of the SNP effect dominance
  bool L_inf_update;              // The likelihood for infection events must be updated
  bool L_trans_update;            // The likelihood for transition events must be updated
};

struct PredAcc{
  string name;                    // The name of the group for which we want to calculate predictrion accuracies
  vector <int> ind;               // A list of individuals
};

struct Jump{
  double size;                    // Stores the size of parameter jumps
  int nac;                        // The number of proposals accepted
  int ntr;                        // The number of proposals tried
  int nevent_tr;                  // The number of individual event proposals (event times only)
  int nevent_fa;                  // The number of individual event failures (event times only)
};

struct Event{
  double time;                    // The time of an event
  int trans;                      // The transition of an event
  int ind;                        // The individual the event happens to
};

struct PrecalcTransParam{         // This stores quantities for a fast calculation of the likelihood with transition params
  double shape_fac;               // Factor multiplying the shape 
  double shape_over_mean;         // Factor multiplying shape over mean
  double shape_log;               // Factor multiplying shape times log of b
  double lgamma_fac;              // Factor multiplying lgamma 
  double con_fac;                 // Constant 
};

struct PrecalcInfParam{           // Stores quantities for proposal in transmission rate
  double beta_fac;
  int log_beta_num;
  double con_fac;
};

struct PrecalcGroupEffect{        // Stores quantities for proposal in group effect
  double beta_fac;
  int log_beta_num;
  double con_fac;
};

struct PrecalcLiSus{              // Stores quantities for proposals to susceptibility individual effects
  double beta_fac;
  int log_beta_num;
};

struct InvCovMat{
  vector < vector <double> > M;   // Stores the inverse correlation matrix
};
  
struct PrecalcIndEff{             // Stores qauntities for fast calculation of Ind Eff likelihood
  double val_sq_fac;            
  double val_fac; 
  double mean;  
  double sd;
};

struct PrecalcLiInf{              // Stores quantities for proposals to infectivity individual effects
  double beta_fac;
  vector <int> inf_log_terms;
};

struct Statistics{                // Stores statistical information
  string mean;                    // The mean
  string CImin, CImax;            // The minimum and maximum of the 95% credible interval
  string ESS;                     // The estimated effective sample size
};

struct Sample{
  vector <double> param_value;    // Stores a sample of the parameter value
};

struct IndPM{
  vector <double> ind_effect_sum; // Keeps a sum of indiviual effects (to calculate average later)
};

struct TimeChange{                // Stores a change to a transition time
  int ind;
  int trans;
  double time;
};

struct AddRemChange{              // Stores a change to adding / removing an individual
  int ind;
  bool infected;
  vector <double> trans_time;
};


#endif
