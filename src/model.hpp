#pragma once

#include <vector>

#include "tinyxml2.h"

using namespace tinyxml2;

using namespace std;

#include "const.hpp"

class Model {
public:
	AlgorithmType algorithm;                     // Sets which algorithm is used 

	unsigned int beta_param;               // Reference for the transmission rate parameter
	InfModel inf_model;                    // Determines if density dependent or frequency dependent

	vector <Param> param;                  // The parameters in the model
	int nparam;                            // The number of parameters

	vector <Derived> derived;              // Quantities derived from paramaters
	int nderived;

	vector <Comp> comp;                    // The compartments in the model
	int ncomp;                             // The number of compartments

	Infectivity infectivity;               // Stores if infectivity has fixed or individual effects

	vector <Trans> trans;                  // The transitions in the model
	int ntrans;                            // The number of transitions

	vector <Group> group;                  // Stores information about each of the groups
	int ngroup;                            // The number of groups

	vector <Covariance> covariance;        // Stores a matrix giving the covariance between individual effects
	int ncovariance;

	vector <Matrix> matrix;                // Stores any relationship matrices

	vector <PredAcc> pred_acc;             // Prediction accuracy for different individual groups of individuals

	vector <Individual> individual;        // Store data about individuals
	int nindividual;

	vector <IndEffect> ind_effect;         // Makes a list of indivdidual effects
	int nind_effect;

	vector <FixedEffect> fixed_effect;     // Makes a list of fixed effects
	int nfixed_effect;

	vector <SNPEffect> snp_effect;         // Makes a list of snp effects
	int nsnp_effect;

	GroupEffect group_effect;              // Group effect

	vector <DiagnosticTest> diag_test;     // Diagnostic test performed on population
	int ndiag_test;

	string output_dir;                     // The output directory

	int nsample;                           // The number of mcmc samples
	int nburnin;                           // The number of burnin samples
	int nthin;                             // The number of samples per trace output

	int nquench;                           // The number of steps over which a quench is carried out
	int nprequench;                        // The number of steps before quench
	double quench_power;

	int nsample_per_gen;                   // The number of samples per generation
	double phi_final;                      // The final inverse temperature (pas only)
	
	bool set_ind_effect_initial;           // Determines if individual effects are set initially

	CLogLog cloglog;                       // Stores information about cloglog
	
	int N;                                 // The number of individuals

	double dt;                             // Timestep (used in emulation)
	unsigned int output_samples;           // The number of samples output (used in emulation)
	unsigned int ngeneration;              // Number of generations (used in emulation)
	unsigned int num_per_gen;              // Number of data samples per generation (used in emulation)
	string params;                         // The parameters which are used in emulation

	bool outp;                             // Set to true if on core 0 

	vector <MeanTimeProposal> mean_time_prop;// Potential mean time proposals

	vector <IeVarJointProposal> ie_var_joint_prop;// Potential mean time proposals

	// In model_initialise.cc
public:
	Model(string file);
	double calculate_derived(const Derived &der, const vector <double> &param_value) const;

private:
	void load_input_file(string file);
	void add_comp(XMLNode *child);
	void set_infectivity(XMLNode *child);
	void add_trans(XMLNode *child);
	void add_parameter(XMLNode *child);
	vector <unsigned int> get_param_sum(string st) const;
	void add_derived(XMLNode *child);
	void add_time_range(XMLNode *child);
	void add_datatable(XMLNode *child);
	void shift_fixed_effects();
	int add_ind_effect(string st);
	int add_fixed_effect(string st, string type);
	int add_snp_effect(string mag, string dom, string type);
	void add_pred_acc(XMLNode *child);
	void add_simulated_column(Table &tab) const;
	void add_diag_test(XMLNode *child);
	void ev_change_initialise();
	void unknown_initialise();
	void ind_effect_proposal_initialise();

	// In model_prior.cc
public:
	double calculate_prior(const vector <double> &param_value) const;
	double calculate_prior(const int th, const double value, const vector <double> &param_value) const;
	double calculate_prior_gradient(const int th, const vector <double> &param_value) const;
	vector <double> prior_sample() const;
	double calculate_prior_change(const int th, const double val_before, const vector <double> &param_value) const;
	void sample_initial_param(vector <double> &param_value, vector <IndValue> &ind_value) const;
	bool inbounds(const vector <double> &paramv) const;
	
	// In model_matrix.cc
public:
	void ind_effect_sample(vector <IndValue> &ind_value, const vector <double> &param_value) const;
	void set_ind_effect(vector <IndValue> &ind_value) const;
	bool check_valid_cov_matrix(const vector <double> &param_value) const;
	vector < vector <double> > set_covariance_matrix(const Covariance &cov, const vector <double> &param_value) const;
	double determinant(const vector < vector <double> > &M) const;
	vector < vector <double> > tensor_product(const vector < vector <double> > &mat1, const vector < vector <double> > &mat2) const;

private:
	void add_identity_matrix();
	void add_matrix(XMLNode *child);
	void add_covariance(XMLNode *child);
	vector < vector <double> > calculate_cholesky_matrix(const vector < vector <double> > &M, bool &illegal) const;
	void print_matrix(string name, const vector < vector <double> > &mat) const;

	// In model_ind_effects.cc
public:
	void set_individual_quantities(vector <IndValue> &ind_value, const vector <double> &param_value) const;
	void set_individual_quantities_group(unsigned int g, vector <IndValue> &ind_value, const vector <double> &param_value) const;
	vector <double> calculate_L_ind_effect(const vector <IndValue> &ind_value, const vector <double> &param_value) const;
	double likelihood_ind_effects(const Covariance &cov, const vector <IndValue> &ind_value, const vector <double> &param_value) const;
	void propose_covariance_matrices(const Covariance &cov, const vector <IndValue> &ind_value, vector <double> &param_value, double &L_ind_effect, double &prior, vector <Jump> &param_jump, const bool burnin, const Quench &quench) const;
	vector <InvCovMat> calculate_inv_cov_matrix(const vector <double> &param_value) const;
	double ind_effect_gradient(unsigned int i, unsigned int ie, const vector <IndValue> &ind_value, const vector <InvCovMat> &inv_cov_matrix) const;
	double covar_param_gradient(unsigned int th, const vector <IndValue> &ind_value, vector <double> param_value) const;
private:
	vector < vector <double> > set_precalc_cov(const Covariance &cov, const vector <IndValue> &ind_value) const;
	double Li_cov_fast(const Covariance &cov, const vector < vector <double> > &precalc, const vector <double> &param_value) const;
	
	// In model_ind_effects_proposals.cc
public:
	void propose_susceptibility_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_inf_events, vector <Jump> &ind_effect_jump, const Quench &quench) const;
	void propose_trans_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, double &L_trans_events, vector <Jump> &ind_effect_jump, const Quench &quench) const;
	void propose_infectivity_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_inf_events, vector <Jump> &ind_effect_jump, const Quench &quench) const;
	void propose_joint_ie_var_initialise();
	void propose_joint_ie_var(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_inf_events, double &L_trans_events, double &prior, vector <Jump> &pjie_jump, const bool burnin, const Quench &quench) const;
	
private:
	vector <PrecalcLiSus> set_precalc_likelihood_sus(const vector <IndValue> &ind_value, const vector <double> &param_value) const;
	double Li_change_sus(const double sus_before, const double sus_after, PrecalcLiSus &prec) const;
	void set_precalc_ind_eff(const int i, PrecalcIndEff &prec, int ie, const vector <IndValue> &ind_value, const vector <InvCovMat> &inv_cov_matrix) const;
	double Li_change_ind_eff(const double before, const double after, const PrecalcIndEff &prec) const;
	vector <PrecalcLiInf> set_precalc_likelihood_inf(vector <double> &group_const, vector < vector <double> > &I_profile, const vector <IndValue> &ind_value, const vector <double> &param_value) const;
	vector <double> likelihood_from_precalc(const vector <double> &group_const, const vector <PrecalcLiInf> &prec, const vector < vector <double> > &I_profilet, const vector <IndValue> &ind_value) const;
	double Li_change_inf(const double infectivity_change, const PrecalcLiInf &prec, const vector <double> &I_profile) const;
	void update_I_profile(const double infectivity_change, const PrecalcLiInf &prec, vector <double> &I_profile) const;

	// model_inf_events.cc
public:
	vector <double> calculate_L_inf_events(const vector <IndValue> &ind_value, const vector <double> &param_value) const;
	void set_initial_events(vector <IndValue> &ind_value, const vector <double> &param_value) const;
	double initial_event_sample(IndValue &ind, const vector <double> &param_value, EvInitSampType ev_samp) const;
	void propose_transmission_rate(const vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_inf_events, double &prior, vector <Jump> &param_jump, const bool burnin, const Quench &quench) const;
	void propose_event_times(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_inf_events, double &L_trans_events, vector <double> &L_diag_test, vector <Jump> &event_jump, const bool burnin, const Quench &quench) const;
	void propose_add_rem(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_inf_events, double &L_trans_events, vector <double> &L_diag_test, vector <Jump> &add_rem_jump, const bool burnin, const Quench &quench) const;
private:
	void get_event_sequence(vector <Event> &event, double &S, const Group &gr, const vector <IndValue> &ind_value) const;
	double likelihood_inf_events(const Group &gr, const vector <IndValue> &ind_value, const vector <double> &param_value) const;

	// model_mean_event_proposals.cc
public:
	void propose_mean_event_times(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_inf_events, double &L_trans_events, vector <double> &L_diag_test, double &prior, vector <Jump> &mean_time_jump, const bool burnin, const Quench &quench) const;
	void propose_mean_event_times_initialise();


	// model_trans_events.cc
public:
	double calculate_L_trans_events(const vector <IndValue> &ind_value, const vector <double> &param_value) const;
	void propose_trans_params(vector <IndValue> &ind_value, vector <double> &param_value, double &L_trans_events, double &prior, vector <Jump> &param_jump, const bool burnin, const Quench &quench) const;
private:
	vector <PrecalcTransParam> set_precalc_trans_param(const vector <IndValue> &ind_value, vector <double> &param_value) const;
	double calculate_L_trans_events_fast(const vector <PrecalcTransParam> &precalc, const vector <double> &param_value) const;

	// model_fixed_effects.cc
public:
	void propose_fixed_effects(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_inf_events, double &L_trans_events, double &prior, vector <Jump> &param_jump, const bool burnin, const Quench &quench) const;

	// model_snp_effects.cc
public:
	void propose_snp_effects(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_inf_events, double &L_trans_events, double &prior, vector <Jump> &param_jump, const bool burnin, const Quench &quench) const;

	// model_group_effect.cc
public:
	void switch_on_group_effect(string sigma);
	void propose_group_effect_sigma(vector <double> &param_value, double &prior, vector <Jump> &param_jump, const bool burnin, const Quench &quench) const;
	void propose_group_effect(const vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_inf_events, double &prior, vector <Jump> &param_jump, const bool burnin, const Quench &quench) const;
private:
	PrecalcGroupEffect set_precalc_group_effect(const Group &gr, const vector <IndValue> &ind_value, const vector <double> &param_value) const;
	double likelihood_inf_events_fast(const PrecalcGroupEffect &prec, const double value) const;

	// In model_diagnostic_test.cc
public:
	vector <double> calculate_L_diag_test(const vector <IndValue> &ind_value, const vector <double> &param_value) const;
	void propose_Se_Sp(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_diag_test, double &prior, vector <Jump> &param_jump, const bool burnin,const Quench &quench) const;
private:
	double likelihood_diag_test(const Group &gr, const vector <IndValue> &ind_value, const vector <double> &param_value) const;

	// In model_funcs.cc
private:
	string tab_name(XMLNode *node) const ;
	string get(XMLNode *node, string attr) const;
	double get_num(XMLNode *node, string attr) const;
	double get_int(XMLNode *node, string attr) const;
	bool exist(XMLNode *node, std::string attr) const;
	int find_comp(string st) const;
	int find_param(string st, ParamType type);
	int find_ind(string id) const;
	Table create_table(const string st) const;
	void remove_comments(string &str) const;
	string brackets(string st) const;
	void output_model() const;
	void check_model() const;
	void change_status(Individual &ind, const Status stat);
	
	// In model_likelihood.cpp
public:
	double likelihood(const vector <double> &param_value, vector <IndValue> &ind_value, const vector <GroupEvent> &gr_ev) const;
	double likelihood_group(unsigned int g, const vector <double> &param_value, vector <IndValue> &ind_value, const vector <GroupEvent> &gr_ev) const;
	void initialise_group_events(vector <GroupEvent> &gr_ev) const;
	vector <unsigned int> get_trlist(unsigned int c) const;
	double posterior_grad_H(vector <double> &grad, vector < vector <double> > &H, const vector <double> &param_value, vector <IndValue> &ind_value, const vector <GroupEvent> &gr_ev, const vector <unsigned int> &param_var_ref, const vector < vector <unsigned int> > &ind_effect_var_ref, unsigned int nvar, bool calc_grad, bool calc_H) const;
	vector < vector <double> > construct_initial_Binv(const vector <double> &param_value, const vector <unsigned int> &param_var_ref, const vector < vector <unsigned int> > &ind_effect_var_ref, unsigned int nvar) const;

	// In model_comp_log_log.cpp
public:
	void initialise_cloglog();
	vector <double> calculate_L_cloglog(const vector <IndValue> &ind_value, const vector <double> &param) const;
	double calculate_L_cloglog_ind(const vector <IndValue> &ind_value, const vector <double> &param, unsigned int m) const;
	void cloglog_propose_transmission_rate(const vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_cloglog, double &prior, vector <Jump> &param_jump, const bool burnin, const Quench &quench) const;
	void cloglog_propose_trans_params(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_cloglog, double &prior, vector <Jump> &param_jump, const bool burnin, const Quench &quench) const;
	void cloglog_propose_susceptibility_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_cloglog, vector <Jump> &ind_effect_jump, const Quench &quench) const;
	void cloglog_propose_infectivity_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_cloglog, vector <Jump> &ind_effect_jump, const Quench &quench) const;
	void cloglog_propose_joint_ie_var(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_cloglog, double &prior, vector <Jump> &pjie_jump, const bool burnin, const Quench &quench) const;
	void set_time_ranges();
};

