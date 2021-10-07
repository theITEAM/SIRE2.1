/// This provides model functions related to diagnostic test results

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "model.hh"
#include "timers.hh"
#include "utils.hh"

/// Calculates the log likelihood for the diagnistic test results and returns a vector (for each group)
vector <double> Model::calculate_L_diag_test(const vector <IndValue> &ind_value, const vector <double> &param_value) const
{
  vector <double> L_diag_test(ngroup);
  for(auto g = 0u; g < ngroup; g++){
    L_diag_test[g] = likelihood_diag_test(group[g],ind_value,param_value);
  }

  return L_diag_test;
}


/// Calculates the log likelihood for the diagnostic test results in a given group
double Model::likelihood_diag_test(const Group &gr, const vector <IndValue> &ind_value, const vector <double> &param_value) const
{
	auto trange = gr.inference_range;
	
	auto L = 0.0;
	for(auto dt = 0; dt < ndiag_test; dt++){ 
		const auto &dia_tes = diag_test[dt];
	
		auto N_pos_inf = 0u, N_neg_inf = 0u, N_pos_notinf = 0u, N_neg_notinf = 0u; 
		for(auto i : gr.ind_ref){
			const auto &indiv = individual[i];
			const auto &ind = ind_value[i];
			const auto &idr = indiv.diag_test_result;
			
			auto c = indiv.initial_comp;
			auto tr = 0;
			for(auto dt = 0; dt < ndiag_test; dt++){
				for(auto &idr : indiv.diag_test_result[dt]){
					auto tt = idr.time;
					while(tr < ntrans && ind.trans_time[tr] < tt){
						c = trans[tr].to; 
						tr++;
					}
					if(dia_tes.comp[c] == true){  // The test is sensitive to the comparment
						if(idr.positive == true) N_pos_inf++;
						else N_neg_inf++;
					}
					else{                         // The test is insensitive to the compartment
						if(idr.positive == true) N_pos_notinf++;
						else N_neg_notinf++;
					}
				}
			}
		}
		
		auto Se = param_value[dia_tes.Se_param];
		auto Sp = param_value[dia_tes.Sp_param];
		
		L += N_pos_inf*log(Se) + N_neg_inf*log(1-Se) + N_pos_notinf*log(1-Sp) + N_neg_notinf*log(Sp);
	}
	
	return L;
}


/// This makes proposals to diagnostic test parameters for sensitivity and specificity
void Model::propose_Se_Sp(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_diag_test, double &prior, vector <Jump> &param_jump, const bool burnin) const
{
  timer[TIME_SE_SP].start();
  for(const auto &dt : diag_test){
		for(auto p = 0u; p < 2; p++){
			auto th = dt.Se_param;
			if(p == 1) th = dt.Sp_param;
			auto &jump = param_jump[th];
    
			for(auto loop = 0; loop < se_sp_prop; loop++){
				auto param_store = param_value[th];

				// Makes a change to one of the covariance matrix parameters
				param_value[th] += normal_sample(0,jump.size);
			
				// Calculates the Metropolis-Hastings probability
				auto L_diag_test_prop = calculate_L_diag_test(ind_value,param_value);
	
				auto prior_change = calculate_prior_change(th,param_store,param_value);
				
				auto sum = prior_change;	
				for(auto g = 0u; g < ngroup; g++) sum += L_diag_test_prop[g] - L_diag_test[g];
			
				auto al = exp(sum);
		
				jump.ntr++;
				if(ran() < al){
					jump.nac++;
					
					L_diag_test = L_diag_test_prop;
					prior += prior_change;
					if(burnin == true) jump.size *= prop_up; 
				}
				else{
					param_value[th] = param_store;
					if(burnin == true) jump.size *= prop_down; 
				}
			}
		}
  }
  timer[TIME_SE_SP].stop();
}


