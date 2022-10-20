/// This provides all the data structures used in the code

#pragma once
#include <vector>

struct Param {
	string name;                    // The name of the parameter
	ParamType type;                 // The type of parameter

	PriorType prior_type;           // The type of prior

	int prior_sd_param;             // References a standard deviation parameter

	double prior_val1;              // Values used to specify the prior
	double prior_val2;
};

struct Derived {
	string name;
	vector <unsigned int> numerator_param;
	vector <unsigned int> denominator_param;
};

struct Comp {
	string name;                    // The name of the compartment

	double infectivity;             // Sets the relative infectivity of the compartment
};

struct Infectivity {
	vector <int> ind_effect;        // References any individual effects acting on infectivity
	vector <int> fixed_effect;      // References any fixed effects acting on infectivity
	vector <int> snp_effect;        // References any snp effects acting on infectivity
};

struct Trans {
	string name;                    // The name of the transition

	int from;                       // From and to which compartment
	int to;

	TransType type;                 // The type of the transition

	int mean_param;                 // Reference for the mean parameter
	int shape_param;                // Reference for the shape parameter

	vector <int> ind_effect;        // References any individual effects acting on transition rate
	vector <int> fixed_effect;      // References any fixed effects acting on transition rate
	vector <int> snp_effect;        // References any snp effects acting on transition rate

	string data_column;             // The data column giving the time of the transition (if specified)
	string data_column_initial;     // The data column giving the time for initial MCMC chain (for testing)

	bool all_ind_effect_proposal;   // Set to true if using  proposals on all individual effects on transition
};

struct TimeRange {
	double tmin;                    // Specifies a time range
	double tmax;
};

struct DiagnosticTest {
	string name;                    // Name of diagnostic test
	vector <bool> comp;             // Set to true if the compartment
	int Se_param;                   // The sensitivity parameter
	int Sp_param;                   // The specificity parameter

	string data_column;             // The data column giving the diagnostic test information
};

struct EvChange {
	int ind;                        // Specifies which individual must have transition times adaptively changed
	int trans;                      // The transition of that individual
	EvSampType ev_samp_type;        // The type of event sampling to be used (either gamma/normal or uniform)
};

struct Group {
	string name;                    // The name of the group
	int index;                      // Sets the group index number
	vector <int> ind_ref;           // The individuals contained in the group
	int nind;                       // The number of individuals
	vector <EvChange> ev_change;    // Stores which events need to be change on which individuals
	int nev_change;
	vector <int> unknown;           // A list of individuals with unknowb disease status
	int nunknown;
	TimeRange inference_range;      // The inference time range
	TimeRange observation_range;    // The observation time range
};

struct GroupEffect {
	bool on;                        // Set to true if group effect is on
	int sigma_param;                // Reference for the standard deviation parameter
	vector <int> param;             // References the group effect parameters
};

struct Table {
	string file;                    // The file name
	unsigned int ncol;              // The number of columns
	unsigned int nrow;              // The number of rows
	vector <string> heading;        // The headings for the columns
	vector <vector <string> > ele;  // The elements of the table
};

struct DiagTestResult {
	bool positive;                  // Set to true if diagnistic test result is positive
	double time;                    // The time of the result
};

struct InfSampler {
	bool on;                        // Set to true of infetion sampler is used
	double tmin;                    // The minimum time for the sampler
	double tmax;                    // The maximum time for the sampler
	vector <int> bin;               // The number of time infection is in time range
	vector <double> log_prob;       // The log of the probability of infection occuring in a cetain time range
	vector <double> prob_sum;       // The sum of the probability

	double sample_time();           // Samples a time from the sampler
	double sample_prob(double t);   // The log of the probability for a given time
};

struct Individual {
	string id;                                // The unique identified for an individual
	string sire, dam;                         // Stores sire and dam information (if included)  
	bool inside_group;                        // This is set to true if the individual is inside a group
	int group;                                // The group to which the individual belongs
	int initial_comp;                         // The initial compartment for the individual
	Status status;                            // Determines the infection status of the individual
	vector <TimeRange> trans_time_range;      // This stores the potential time range for a given transition
	vector <double> trans_obs_time;           // The observed time for a transition
	vector <double> trans_time_initial;       // Set from 'data_column_initial'
	vector <double> ind_effect_value;         // Stores actual individual effect values (e.g breeding values)
	vector <double> fixed_effect_X;           // The design matrix multiplying the fixed effects
	vector <double> fixed_effect_X_unshifted; // The unshifted design matrix
	vector <SNP> SNP_genotype;                // The genotype for different SNPs
	vector <vector <DiagTestResult> > diag_test_result; // Series of diagnostic test results for an individual
};

struct IndValue {
	bool infected;                  // Sets if the individual on the MCMC chain is infected
	int index;                      // Sets the individual number
	vector <double> trans_time;     // Sets the transition times for the individual
	vector <double> ind_effect;     // The values for the individual effects
	// These quantities are derived from other model parameters:
	double susceptibility;          // Individual susceptibiliy
	double inf_single;              // Individual infectivity
	vector <double> trans_infectivity_change;    // Change in individual infectivity going down a transition
	vector <double> trans_mean;     // The means of transitions

	InfSampler inf_sampler;         // An infection sampler (used if particle undergoes add / rem infection)
};

struct Element {
	int i;
	double val;
};

struct Matrix {
	string name;                                 // The name of the matrix
	vector <double> Ainvdiag;                    // The diagonal members of the inverse matrix
	vector < vector <Element> > Ainvlist;        // Stores the inverse matrix in sparse format
	vector < vector <Element> > Ainvlist2;       // Stores only when i < j
	vector < vector <double> > A;                // Stores the full matrix
	vector < vector <double> > cholesky_matrix;  // Stores the Cholesky decomposition of A
};

struct SparseMatrixRow {
	double diag;
	vector <unsigned int> ele;
	vector <double> val;
};

struct Covariance {
	int matrix;                     // Reference to the relationship matrix
	int E;                          // The number of individual effects
	vector <int> ind_effect_ref;    // Stores which individual effects the matrix refers to
	vector <int> var_param;         // References variance parameters
	vector < vector <int> > cor_param;  // References correlation parameters
	vector <int> param_list;        // Lists all parameters used in the covariance matrix
};

struct IndEffect {
	string name;                    // The name of individual effect
	int covar_ref;                  // The covariance the individual effect appears in
	int covar_num;                  // The position in the covariace matrix the individual effect appears
	vector <int> susceptibility;    // Stores if changes susceptibility
	vector <int> infectivity_comp;  // Stores compartments it changes the infectivity of
	vector <int> trans_mean;        // Stores which transitions it changes the mean of

	bool single_ind_effect_proposal;// Set to true if individual effect updated seperately
};

struct FixedEffect {
	string name;                    // The name of the fixed effect
	int param;                      // The parameter of the fixed effect
	bool L_inf_update;              // The likelihood for infection events must be updated
	bool L_trans_update;            // The likelihood for transition events must be updated
};

struct SNPEffect {
	string name;                    // The name of the SNP effect
	int param_mag;                  // The parameter of the SNP effect magnitude
	int param_dom;                  // The parameter of the SNP effect dominance
	bool L_inf_update;              // The likelihood for infection events must be updated
	bool L_trans_update;            // The likelihood for transition events must be updated
};

struct PredAcc {
	string name;                    // The name of the group for which we want to calculate predictrion accuracies
	vector <int> ind;               // A list of individuals
};

struct Jump {
	double size;                    // Stores the size of parameter jumps
	int nac;                        // The number of proposals accepted
	int ntr;                        // The number of proposals tried
	int nfa;                        // The number of proposals failed
	int nevent_tr;                  // The number of individual event proposals (event times only)
	int nevent_fa;                  // The number of individual event failures (event times only)
};

struct Event {
	double time;                    // The time of an event
	int trans;                      // The transition of an event
	int ind;                        // The individual the event happens to
};

struct PrecalcTransParam {        // This stores quantities for a fast calculation of the likelihood with transition params
	double shape_fac;               // Factor multiplying the shape
	double shape_over_mean;         // Factor multiplying shape over mean
	double shape_log;               // Factor multiplying shape times log of b
	double lgamma_fac;              // Factor multiplying lgamma
	double con_fac;                 // Constant
};

struct PrecalcInfParam {          // Stores quantities for proposal in transmission rate
	double beta_fac;
	int log_beta_num;
	double con_fac;
};

struct PrecalcGroupEffect {       // Stores quantities for proposal in group effect
	double beta_fac;
	int log_beta_num;
	double con_fac;
};

struct PrecalcLiSus {             // Stores quantities for proposals to susceptibility individual effects
	double beta_fac;
	int log_beta_num;
};

struct InvCovMat {
	vector < vector <double> > M;   // Stores the inverse correlation matrix
};

struct PrecalcIndEff {            // Stores qauntities for fast calculation of Ind Eff likelihood
	double val_sq_fac;
	double val_fac;
	double mean;
	double sd;
};

struct LogTerm {                  // Stores information about log terms in likelihood
	int list;
	double relinf;
};

struct PrecalcLiInf {             // Stores quantities for proposals to infectivity individual effects
	double beta_fac;
	vector <LogTerm> inf_log_terms;
};

struct Statistics {               // Stores statistical information
	string mean;                    // The mean
	string CImin, CImax;            // The minimum and maximum of the 95% credible interval
	string ESS;                     // The estimated effective sample size
};

struct Sample {
	vector <double> param_value;    // Stores a sample of the parameter value
};

struct IndPM {
	vector <double> ind_effect_sum; // Keeps a sum of indiviual effects (to calculate average later)
	vector <double> ind_effect_sum2; // Keeps a square sum of indiviual effects 
};

struct TimeChange {               // Stores a change to a transition time
	int ind;
	int trans;
	double time;
};

struct AddRemChange {             // Stores a change to adding / removing an individual
	int ind;
	bool infected;
	vector <double> trans_time;
};

struct Anneal {                   // Parameters used for annealing
	double phi_L;                   // Inverse temperature on the event likelihood
	double phi_IE;                  // Inverse temperature on individual effects
	double phi_DT;                  // Inverse temperature on the disease diagnostic test likelihood
	double phi_Pr;                  // Inverse temperature on the prior 
};

struct Variable {                 // Variables (used in CMA-ES)
	VariableType type;              // Type of variable
	unsigned ind;                   // Individual (if applicable)
	unsigned ref;                   // Reference index
};

struct VariableSample{            // Stores information about a parameter sample from the posterior
	double L;                       // Stores the error function 
	vector <double> value;          // A parameter sample
};

struct Generation {               // Stores information about generation (in CMA-ES)
	double Lav;                     // The average likelihood
};

struct MatrixElement {            // Identifies an element of a matrix
	unsigned int j;                 // y position in matrix
	unsigned int i;                 // x position in matrix
};

struct TransSection {
	double dt;                     // The time for the section
	vector <unsigned int> trlist;  // List of possible transitions
	bool trans_end;                // Determines if section ends in a transition
};

struct IndEvent {
	vector <TransSection> trans_section;
};

struct GroupEvent {              // Stores group events (used in CMA-ES)
	unsigned int T;
	vector < vector <IndEvent> > ind_ev;  
};

struct GradCalcList {            // Sets which variable on which group to calculate gradient
	unsigned int g;
	unsigned int v;
	unsigned int i;
};

struct HessianCalcList {            // Sets which variable on which group to calculate gradient
	unsigned int g;
	unsigned int v1;
	unsigned int v2;
};

struct GradIECalcList {          // Sets which variable on which group to calculate gradient in IE
	unsigned int i;
	unsigned int ie;
};

struct GroupVariable {           // List variables associated with a group
	vector <unsigned int> var_list;// Variables 
	vector <bool> mask;            // Set to true for the variables in the list 
};

struct VarBlock {                // Varialbes in a block
	vector <unsigned int> var_list;// Lists varialbe to be updated in a block
};

struct DataSample {
	vector <double> em_var_value;
	vector <double> em_var_value_scale;
	double post;
};

struct FE {
	vector <unsigned int> em_var_list;
};

struct EmulatorEstimate {
	double val;
	double sd;
};

struct InfList {
	unsigned int k;
	double val;
};

struct DepList {
	unsigned int v;
	double val;
};

struct MeanTimeProposal {
	string name;
	unsigned int trans_move;
};

struct IeVarJointProposal{
	string name;
	unsigned int covar;
	unsigned int pos;
};

struct InfInd{
	unsigned int i;       // The individual which is infected
	double infectivity;       // The compartment that individual is in
};

struct CLogLogRef {
	vector <unsigned int> direct;
	vector <unsigned int> indirect;
};

struct L_CLogLog{
	unsigned int i;
	unsigned int tr;
	bool result;               // Set to true if the individual makes the transition
	double t;
	
	vector <InfInd> inf_ind;
};

struct CLogLog{
	bool on;                             // Sets if the complimentary log-log link function is used for the likelihood
	double DeltaT;                       // The timestep used for cloglog
	bool geometric_approx;               // Determines if the geometric mean approx used
	
	vector <L_CLogLog> L_list; 
	
	vector < vector <CLogLogRef> > L_ref;// Information about individuals
};


