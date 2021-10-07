#ifndef SIRE__MODEL_HH
#define SIRE__MODEL_HH

#include <vector>

#include "tinyxml2.h"

using namespace tinyxml2;

using namespace std;

#include "const.hh"

class Model
{
  public:
    int beta_param;                        // Reference for the transmission rate parameter
    InfModel inf_model;                    // Determines if density dependent or frequency dependent
  
    vector <Param> param;                  // The parameters in the model
    int nparam;                            // The number of parameters
    
    vector <Comp> comp;                    // The compartments in the model
    int ncomp;                             // The number of compartments
    
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
    
    int N;                                 // The number of individuals
  
  // In model_initialise.cc
  public:
    Model(string file);
    
  private: 
    void load_input_file(string file);
    void add_comp(XMLNode *child);
    void add_trans(XMLNode *child);
    void add_parameter(XMLNode *child);
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
    vector <double> prior_sample() const;
    double calculate_prior_change(const int th, const double val_before, const vector <double> &param_value) const;
  
  // In model_matrix.cc
  public:
    void ind_effect_sample(vector <IndValue> &ind_value, const vector <double> &param_value) const;
    vector < vector <double> > set_covariance_matrix(const Covariance &cov, const vector <double> &param_value) const;
    double determinant(const vector < vector <double> > &M) const;
    vector < vector <double> > tensor_product(const vector < vector <double> > &mat1, const vector < vector <double> > &mat2) const;
	  vector <vector <double> > invert_matrix(const vector <vector <double> > &mat) const;  
  private:
    void add_identity_matrix();        
    void add_matrix(XMLNode *child);
    void add_covariance(XMLNode *child);  
    vector < vector <double> > calculate_cholesky_matrix(const vector < vector <double> > &M, bool &illegal) const;
    void print_matrix(string name, const vector < vector <double> > &mat) const;
    
  // In model_ind_effects.cc
  public:
    void set_individual_quantities(vector <IndValue> &ind_value, const vector <double> &param_value) const;
    vector <double> calculate_L_ind_effect(const vector <IndValue> &ind_value, const vector <double> &param_value) const;
    double likelihood_ind_effects(const Covariance &cov, const vector <IndValue> &ind_value, const vector <double> &param_value) const;
    void propose_covariance_matrices(const Covariance &cov, const vector <IndValue> &ind_value, vector <double> &param_value, double &L_ind_effect, double &prior, vector <Jump> &param_jump, const bool burnin) const;
  private:
    vector < vector <double> > set_precalc_cov(const Covariance &cov, const vector <IndValue> &ind_value) const;
    double Li_cov_fast(const Covariance &cov, const vector < vector <double> > &precalc, const vector <double> &param_value) const;
    vector <InvCovMat> calculate_inv_cov_matrix(const vector <double> &param_value) const;
    
  // In model_ind_effects_proposals.cc
  public:
    void propose_susceptibility_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_inf_events, vector <Jump> &ind_effect_jump) const;
    void propose_trans_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, double &L_trans_events, vector <Jump> &ind_effect_jump) const;
    void propose_infectivity_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_inf_events, vector <Jump> &ind_effect_jump) const;
  private:
    vector <PrecalcLiSus> set_precalc_likelihood_sus(const vector <IndValue> &ind_value, const vector <double> &param_value) const;
    double Li_change_sus(const double sus_before, const double sus_after, PrecalcLiSus &prec) const;
    void set_precalc_ind_eff(const int i, PrecalcIndEff &prec, int ie, const vector <IndValue> &ind_value, const vector <InvCovMat> &inv_cov_matrix) const;
    double Li_change_ind_eff(const double before, const double after, const PrecalcIndEff &prec) const;
    vector < vector <PrecalcLiInf> > set_precalc_likelihood_inf(vector < vector <double> > &I_profile, const vector <IndValue> &ind_value, const vector <double> &param_value) const;
    double Li_change_inf(const double infectivity_change, const PrecalcLiInf &prec, const vector <double> &I_profile) const;
    void update_I_profile(const double infectivity_change, const PrecalcLiInf &prec, vector <double> &I_profile) const;
		
  // model_inf_events.cc
  public:
    vector <double> calculate_L_inf_events(const vector <IndValue> &ind_value, const vector <double> &param_value) const;
    void set_initial_events(vector <IndValue> &ind_value, const vector <double> &param_value) const;
    double event_sample(IndValue &ind, const vector <double> &param_value, EvInitSampType ev_samp) const;
    void propose_transmission_rate(const vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_inf_events, double &prior, vector <Jump> &param_jump) const;
    void propose_event_times(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_inf_events,double &L_trans_events, vector <double> &L_diag_test, vector <Jump> &event_jump, const bool burnin) const;
    void propose_add_rem(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_inf_events, double &L_trans_events, vector <double> &L_diag_test, vector <Jump> &add_rem_jump, const bool burnin) const;
  private:
    void get_event_sequence(vector <Event> &event, double &S, const Group &gr, const vector <IndValue> &ind_value) const;
    double likelihood_inf_events(const Group &gr, const vector <IndValue> &ind_value, const vector <double> &param_value) const;
  
  // model_trans_events.cc
  public:
    double calculate_L_trans_events(const vector <IndValue> &ind_value, const vector <double> &param_value) const;
    void propose_trans_params(vector <IndValue> &ind_value, vector <double> &param_value, double &L_trans_events, double &prior, vector <Jump> &param_jump, const bool burnin) const;
  private:
    vector <PrecalcTransParam> set_precalc_trans_param(const vector <IndValue> &ind_value, vector <double> &param_value) const;
    double calculate_L_trans_events_fast(const vector <PrecalcTransParam> &precalc, const vector <double> &param_value) const;
    
  // model_fixed_effects.cc
  public:
    void propose_fixed_effects(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_inf_events, double &L_trans_events, double &prior, vector <Jump> &param_jump, const bool burnin) const;
    
	// model_snp_effects.cc
	public:
		void propose_snp_effects(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_inf_events, double &L_trans_events, double &prior, vector <Jump> &param_jump, const bool burnin) const;
		
  // model_group_effect.cc
  public:
    void switch_on_group_effect(string sigma);
    void propose_group_effect_sigma( vector <double> &param_value, double &prior, vector <Jump> &param_jump, const bool burnin) const;
    void propose_group_effect(const vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_inf_events, double &prior, vector <Jump> &param_jump, const bool burnin) const;
  private:
    PrecalcGroupEffect set_precalc_group_effect(const Group &gr, const vector <IndValue> &ind_value, const vector <double> &param_value) const;
    double likelihood_inf_events_fast(const PrecalcGroupEffect &prec, const double value) const;
    
	// In model_diagnostic_test.cc
	public:
		vector <double> calculate_L_diag_test(const vector <IndValue> &ind_value, const vector <double> &param_value) const;
		void propose_Se_Sp(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_diag_test, double &prior, vector <Jump> &param_jump, const bool burnin) const;
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
};

#endif
