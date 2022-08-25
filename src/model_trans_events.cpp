/// This includes model functions related to transition events (not infection transitions)

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "model.hpp"
#include "timers.hpp"
#include "utils.hpp"


/// Calculates the log likelihood for the transition events (not infection)
double Model::calculate_L_trans_events(const vector <IndValue> &ind_value, const vector <double> &param_value) const {
	// cout << "Model::calculate_L_trans_events()" << endl; // DEBUG
	auto L = 0.0;
	for (auto tr = 1u; tr < ntrans; tr++) {
		switch (trans[tr].type) {
			case GAMMA: case EXP: {
				auto shape = 1.0; if(trans[tr].type == GAMMA) shape = param_value[trans[tr].shape_param];
				for (const auto &ind : ind_value) {
					if (ind.infected == true) {
						auto dt = ind.trans_time[tr] - ind.trans_time[tr - 1];
						auto mean = ind.trans_mean[tr];
						L += gamma_probability(dt, mean, shape);
					}
				}
			}
			break;

			default:
				emsg("Not set");
				break;
		}
	}

	return L;
}


/// This makes proposals to transition parameters
void Model::propose_trans_params(vector <IndValue> &ind_value, vector <double> &param_value, double &L_trans_events, double &prior, vector <Jump> &param_jump, const bool burnin, const Quench &quench) const {
	// cout << "Model::propose_trans_params()" << endl; // DEBUG
	timer[TIME_TRANS_PARAM_INIT].start();
	auto precalc = set_precalc_trans_param(ind_value, param_value);
	timer[TIME_TRANS_PARAM_INIT].stop();

	timer[TIME_TRANS_PARAM].start();
	for (auto loop = 0; loop < trans_param_prop; loop++) {
		for (auto th = 0; th < nparam; th++) {
			if (param[th].type == TRANS_MEAN || param[th].type == TRANS_SHAPE) {
				auto &jump = param_jump[th];
				auto param_store = param_value[th];

				// Makes a change to one of the transition parameters
				param_value[th] += normal_sample(0, jump.size);

				// Calculates the Metropolis-Hastings acceptance probability
				auto prior_change = calculate_prior_change(th, param_store, param_value);
				auto L_propose = calculate_L_trans_events_fast(precalc, param_value);

				auto al = exp(quench.phi_L*(L_propose - L_trans_events) + quench.phi_Pr*prior_change);
				if (prior_change == ZERO_PRIOR) al = 0;

				jump.ntr++;
				if (MH_proposal(al,16)) {
					jump.nac++;
					L_trans_events = L_propose;
					prior += prior_change;
					if (burnin == true)
						jump.size *= prop_up;
				} else {
					param_value[th] = param_store;
					if (burnin == true)
						jump.size *= prop_down;
				}
			}
		}
	}

	set_individual_quantities(ind_value, param_value);          // Sets individual transition parameters and infectivity

	timer[TIME_TRANS_PARAM].stop();
}


/// Pre-calculates quantities for a fast likelihood calculation
vector <PrecalcTransParam> Model::set_precalc_trans_param(const vector <IndValue> &ind_value, vector <double> &param_value) const {
	// cout << "Model::set_precalc_trans_param()" << endl; // DEBUG
	vector <PrecalcTransParam> precalc(ntrans);

	auto L = 0.0;
	for (auto tr = 1u; tr < ntrans; tr++) {
		auto &prec = precalc[tr];
		const auto &tra = trans[tr];
		prec.shape_fac = 0;
		prec.shape_over_mean = 0;
		prec.shape_log = 0;
		prec.lgamma_fac = 0;
		prec.con_fac = 0;

		switch (trans[tr].type) {
			case GAMMA: case EXP: {
				auto mean_av = param_value[tra.mean_param];
				auto shape = 1.0; if(trans[tr].type == GAMMA) shape = param_value[tra.shape_param];
				for (const auto &ind : ind_value) {
					if (ind.infected == true) {
						auto dt = ind.trans_time[tr] - ind.trans_time[tr - 1];
						auto mean = ind.trans_mean[tr];
						auto lg = log(dt);
						prec.shape_fac += lg + log(mean_av / mean);
						prec.shape_over_mean -= dt * mean_av / mean;
						prec.shape_log++;
						prec.lgamma_fac--;
						prec.con_fac -= lg;
					}
				}
			}
			break;

			default:
				emsg("Not set");
				break;
		}
	}

	return precalc;
}


/// Calculates likelihood fast (using pre-calculated quantities)
double Model::calculate_L_trans_events_fast(const vector <PrecalcTransParam> &precalc, const vector <double> &param_value) const {
	// cout << "Model::calculate_L_trans_events_fast()" << endl; // DEBUG
	auto L = 0.0;
	for (auto tr = 1u; tr < ntrans; tr++) {
		const auto &prec = precalc[tr];

		switch (trans[tr].type) {
			case GAMMA: case EXP: {
				auto mean_av = param_value[trans[tr].mean_param];
				auto shape = 1.0; if(trans[tr].type == GAMMA) shape = param_value[trans[tr].shape_param];
				L += prec.con_fac + shape * prec.shape_fac + (shape / mean_av) * prec.shape_over_mean +
				     shape * log(shape / mean_av) * prec.shape_log + prec.lgamma_fac * lgamma(shape);
			}
			break;

			default:
				emsg("Not set");
				break;
		}
	}

	return L;
}
