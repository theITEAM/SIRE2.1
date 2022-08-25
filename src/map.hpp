#pragma once

using namespace std;

#include "model.hpp"
#include "mpi.hpp"

class MAP {
public:
	unsigned int nvar;                             // Number of variables
	vector <Variable> var;                         // Stores all the variables in the system
	int var_per_core;                              // Used in gradient decent
	
	vector <double> param_value;                   // The model parameter values
	vector <IndValue> ind_value;                   // Stores individual properties on the chain

	vector <GroupEvent> gr_ev;                    // Stores information about observed transition events 

	vector <GroupVariable> group_variable;        // List all variables which affect likelihood in a group

	vector <unsigned int> param_var_ref;          // References variable from parameter
	vector < vector <unsigned int> > ind_effect_var_ref; // References variable from individual effect

	MAP(const Model &model);
	void cmaes(int seed);
	void gradient_descent();
	void REML();
	
private:
	void initialise_variable();
	vector < vector <double> > sample_param(const vector <double> &mean, const vector < vector <double> > &C, const double sigma, const unsigned int num) const;
	void generate_samples(vector <VariableSample> &ps_per_core, const vector <double> &mean, const vector < vector <double> > &C, const double sigma);
	vector < vector <double> > calculate_covariance_martrix(const vector <double> &mean);
	bool terminate_generation(const vector <double> &EFbest_store);
	double scale_covariance_martrix(const vector <double> &mean, const vector < vector <double> > &C);
	vector < vector <double> > calculate_hessian(const vector <double> &mean, MatrixType mattype);
	vector <double> get_value_from_param();
	vector < vector <double> > get_C_from_param();
	void update_param_from_value(const vector <double> &value);
	void update_param_from_value(const vector <double> &value, unsigned int v);
	double calculate_posterior(const vector <double> &value);
	void output_parameter_samples(const vector <double> &mean, const vector < vector <double> > &C, const double sigma);
	vector <double> finite_different_grid_size(const vector <double> &value);
	vector <double> calculate_gradient(vector <double> value, const vector <unsigned int> &var_list);
	vector <double> calculate_gradient_slow(vector <double> value, const vector <unsigned int> &var_list);
	void var_add(vector <unsigned int> &vec, unsigned int th) const;
	void indeff_add(vector <unsigned int> &vec, unsigned int i, unsigned int ie) const;
	void check_gradient();
	void check_gradient_ie();
	vector <VarBlock> set_var_block() const;
	
	// Used in REML
	double likelihood_ie_marginal();
	vector < vector <double> > calculate_gradient_ie();
	vector < vector <double> > calculate_hessian_fast(vector <double> value);
	void check_grad_H();

	const Model &model;
	Mpi mpi;
};
