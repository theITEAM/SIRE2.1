/// This gives model functions related to the prior

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "model.hpp"
#include "utils.hpp"


/// Samples a set of values from the prior
vector <double> Model::prior_sample() const {
	// cout << "Model::prior_sample()" << endl; // DEBUG
	vector <double> param_value(nparam);

	for (auto th = 0; th < nparam; th++) {
		const auto &par = param[th];
		switch (par.prior_type) {
			case FIXED_PRIOR:
				param_value[th] = par.prior_val1;
				break;
				
			case FLAT_PRIOR:
				param_value[th] = par.prior_val1 + ran() * (par.prior_val2 - par.prior_val1);
				break;

			case NORMAL_FROM_SD_PRIOR:
				param_value[th] = normal_sample(0, param_value[par.prior_sd_param]);
				break;

			default:
				emsg("Prior not specified");
				break;
		}
	}

	return param_value;
}


/// Calculates the log of the prior probability
double Model::calculate_prior(const vector <double> &param_value) const {
	// cout << "Model::calculate_prior()" << endl; // DEBUG
	auto prior = 0.0;
	for (auto th = 0; th < nparam; th++)
		prior += calculate_prior(th, param_value[th], param_value);
	return prior;
}


/// Determines if the parameters is within the bounds set by the prior
bool Model::inbounds(const vector <double> &paramv) const 
{
	for (auto th = 0; th < nparam; th++){
		const auto &par = param[th];
		auto value = paramv[th];
		
		switch (par.prior_type) {
			case FIXED_PRIOR: 
				if(value !=  par.prior_val1) return false;
				break;
			
			case FLAT_PRIOR: 
				if (value < par.prior_val1 || value > par.prior_val2) return false;
				break;

			case NORMAL_FROM_SD_PRIOR:
				break;
		}
	}
	return true;	
}


/// Calculates the prior for a particular parameter th with a value
double Model::calculate_prior(const int th, const double value, const vector <double> &param_value) const {
	// cout << "Model::calculate_prior()" << endl; // DEBUG
	const auto &par = param[th];
	switch (par.prior_type) {
		case FIXED_PRIOR: 
			if(value !=  par.prior_val1) return ZERO_PRIOR;
			else return 0;
		
		case FLAT_PRIOR: {
			auto val1 = par.prior_val1;
			auto val2 = par.prior_val2;

			if (value < val1 || value > val2)
				return ZERO_PRIOR;
			else
				return log(1.0 / (val2 - val1));
		}

		case NORMAL_FROM_SD_PRIOR: {
			auto sd = param_value[par.prior_sd_param];
			if (sd <= 0)
				return ZERO_PRIOR;
			else
				return normal_probability(value, 0, sd);
		}

		default:
			emsg("Prior not specified");
			break;
	}
	return 0.0;
}


/// Returns the gradient in the log of the prior probability
double Model::calculate_prior_gradient(const int th, const vector <double> &param_value) const 
{
	const auto &par = param[th];
	switch (par.prior_type) {
		case NORMAL_FROM_SD_PRIOR:
			{
				auto sd = param_value[par.prior_sd_param];
				return -param_value[th]/(sd*sd);
			}
	}
	
	if(group_effect.on == true && th == group_effect.sigma_param){
		auto sd = param_value[th];
		
		auto grad = -ngroup/sd;
		for(auto gr_th : group_effect.param){
			grad += param_value[gr_th]*param_value[gr_th]/(sd*sd*sd);
		}
		
		return grad;
	}
	
	return 0;
}


/// Calculates the change in prior when a particular parameter th changes
double Model::calculate_prior_change(const int th, const double val_before, const vector <double> &param_value) const {
	// cout << "Model::calculate_prior_change()" << endl; // DEBUG
	auto prior_after = calculate_prior(th, param_value[th], param_value);
	if (prior_after == ZERO_PRIOR)
		return ZERO_PRIOR;
	auto prior_before = calculate_prior(th, val_before, param_value);

	return prior_after - prior_before;
}


/// Samples the initial parameter set
void Model::sample_initial_param(vector <double> &param_value, vector <IndValue> &ind_value) const 
{
	auto loop = 0;
	do {                                                               // Repeat until a valid parameter set generated
		param_value = prior_sample();                             // Parameters sampled from prior
		loop++;
	} while (loop < INITIAL_PARAM_SAMPLE && check_valid_cov_matrix(param_value) == false);
	if (loop == INITIAL_PARAM_SAMPLE)
		emsg("Valid covariance matrices could not be specified. Please check the prior.");
	
	ind_value.resize(N);                                        // Sets a variable for individuals information
	for (auto i = 0; i < N; i++)
		ind_value[i].index = i;        // Sets index for each individual

	ind_effect_sample(ind_value, param_value);                 // Individual effects sampled from prior
	
	if(set_ind_effect_initial == true){                       	// Individual effects specfied 
		set_ind_effect(ind_value);
	}	
	
	// Arbitarily sets susceptibility and infectivity for individuals not in groups 
	for(auto i = 0; i < N; i++){  
		if(individual[i].group == UNSET){
			ind_value[i].susceptibility = 1;
			ind_value[i].inf_single = 1;
			ind_value[i].trans_infectivity_change.resize(ntrans);
		}
	}
	
	set_individual_quantities(ind_value, param_value);          // Sets individual transition parameters and infectivity
}

