/// This provides model functions related to any snp effects in the model

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "model.hpp"
#include "timers.hpp"
#include "utils.hpp"


/// This makes proposals to snp effect parameters
void Model::propose_snp_effects(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_inf_events, double &L_trans_events, double &prior, vector <Jump> &param_jump, const bool burnin, const Quench &quench) const {
	// cout << "Model::propose_snp_effects()" << endl; // DEBUG
	timer[TIME_SNP_EFFECTS].start();
	for (auto se = 0; se < nsnp_effect; se++) {
		auto &snp_ef = snp_effect[se];
		for (auto p = 0; p < 2; p++) {
			auto th = snp_ef.param_mag;
			if (p == 1)
				th = snp_ef.param_dom;

			auto &jump = param_jump[th];

			for (auto loop = 0; loop < snp_effect_prop; loop++) {
				auto param_store = param_value[th];

				// Makes a change to one of the covariance matrix parameters
				param_value[th] += normal_sample(0, jump.size);

				set_individual_quantities(ind_value, param_value);

				// Calculates the Metropolis-Hastings probability
				vector <double> L_inf_events_prop;
				double L_trans_events_prop;

				auto sum = 0.0;

				if (snp_ef.L_inf_update == true) {
					L_inf_events_prop = calculate_L_inf_events(ind_value, param_value);
					for (auto g = 0; g < ngroup; g++)
						sum += L_inf_events_prop[g] - L_inf_events[g];
				}

				if (snp_ef.L_trans_update == true) {
					L_trans_events_prop = calculate_L_trans_events(ind_value, param_value);
					sum += L_trans_events_prop - L_trans_events;
				}

				auto prior_change = calculate_prior_change(th, param_store, param_value);
				auto al = exp(quench.phi_L*sum + quench.phi_Pr*prior_change);

				jump.ntr++;
				if (MH_proposal(al,14)) {
					jump.nac++;

					if (snp_ef.L_inf_update == true)
						L_inf_events = L_inf_events_prop;
					if (snp_ef.L_trans_update == true)
						L_trans_events = L_trans_events_prop;
					prior += prior_change;
					if (burnin == true)
						jump.size *= prop_up;
				} else {
					param_value[th] = param_store;
					set_individual_quantities(ind_value, param_value);
					if (burnin == true)
						jump.size *= prop_down;
				}
			}
		}
	}
	timer[TIME_SNP_EFFECTS].stop();
}
