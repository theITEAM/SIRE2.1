/// This provides model functions related to diagnostic test results

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "model.hpp"
#include "timers.hpp"
#include "utils.hpp"

/// Calculates the log likelihood for the diagnistic test results and returns a vector (for each group)
vector <double> Model::calculate_L_diag_test(const vector <IndValue> &ind_value, const vector <double> &param_value) const {
	// cout << "Model::calculate_L_diag_test()" << endl; // DEBUG
	vector <double> L_diag_test(ngroup);
	for (auto g = 0; g < ngroup; g++)
		L_diag_test[g] = likelihood_diag_test(group[g], ind_value, param_value);

	return L_diag_test;
}


/// Calculates the log likelihood for the diagnostic test results in a given group
double Model::likelihood_diag_test(const Group &gr, const vector <IndValue> &ind_value, const vector <double> &param_value) const {
	// cout << "Model::likelihood_diag_test()" << endl; // DEBUG
	auto trange = gr.inference_range;

	auto L = 0.0;
	for (auto dt = 0; dt < ndiag_test; dt++) {
		const auto &dia_tes = diag_test[dt];

		auto N_pos_inf = 0, N_neg_inf = 0, N_pos_notinf = 0, N_neg_notinf = 0;
		for (auto i : gr.ind_ref) {
			const auto &indiv = individual[i];
			const auto &ind = ind_value[i];
			const auto &idr = indiv.diag_test_result;

			auto c = indiv.initial_comp;
			auto tr = 0;
			for (auto dt = 0; dt < ndiag_test; dt++) {
				for (auto &idr : indiv.diag_test_result[dt]) {
					auto tt = idr.time;
					while (tr < ntrans && ind.trans_time[tr] < tt) {
						c = trans[tr].to;
						tr++;
					}
					if (dia_tes.comp[c] == true) { // The test is sensitive to the comparment
						if (idr.positive == true)
							N_pos_inf++;
						else
							N_neg_inf++;
					} else {                      // The test is insensitive to the compartment
						if (idr.positive == true)
							N_pos_notinf++;
						else
							N_neg_notinf++;
					}
				}
			}
		}

		auto Se = param_value[dia_tes.Se_param];
		auto Sp = param_value[dia_tes.Sp_param];

		L += N_pos_inf * log(Se) + N_neg_inf * log(1 - Se) + N_pos_notinf * log(1 - Sp) + N_neg_notinf * log(Sp);
	}

	return L;
}


/// This makes proposals to diagnostic test parameters for sensitivity and specificity
void Model::propose_Se_Sp(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_diag_test, double &prior, vector <Jump> &param_jump, const bool burnin, const Anneal &anneal) const {
	// cout << "Model::propose_Se_Sp()" << endl; // DEBUG
	timer[TIME_SE_SP].start();
	for (const auto &dt : diag_test) {
		for (auto p = 0; p < 2; p++) {
			auto th = dt.Se_param;
			if (p == 1)
				th = dt.Sp_param;
			auto &jump = param_jump[th];

			for (auto loop = 0; loop < se_sp_prop; loop++) {
				auto param_store = param_value[th];

				// Makes a change to one of the covariance matrix parameters
				param_value[th] += normal_sample(0, jump.size);

				if(inbounds(param_value) == false) 	param_value[th] = param_store;
				else{
					// Calculates the Metropolis-Hastings probability
					auto L_diag_test_prop = calculate_L_diag_test(ind_value, param_value);

					auto sum = 0.0;
					for (auto g = 0; g < ngroup; g++)
						sum += L_diag_test_prop[g] - L_diag_test[g];

					auto prior_change = calculate_prior_change(th, param_store, param_value);

					auto al = exp(anneal.phi_DT*sum + anneal.phi_Pr*prior_change);

					jump.ntr++;
					if (MH_proposal(al,11)) {
						jump.nac++;

						L_diag_test = L_diag_test_prop;
						prior += prior_change;
						if (burnin == true)
							jump.size *= prop_up;
					} 
					else {
						param_value[th] = param_store;
						if (burnin == true)
							jump.size *= prop_down;
					}
				}
			}
		}
	}
	timer[TIME_SE_SP].stop();
}


