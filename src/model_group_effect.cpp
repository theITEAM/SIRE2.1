/// This gives model functions related to the group effect (if it exists)

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "model.hpp"
#include "timers.hpp"
#include "utils.hpp"

/// This switches on a group effect and adds extra model parameters for each of the group effects
void Model::switch_on_group_effect(string sigma) {
	// cout << "Model::switch_on_group_effect()" << endl; // DEBUG
	group_effect.on = true;
	group_effect.sigma_param = find_param(sigma, GROUP_SD);

	for (auto gr = 0; gr < ngroup; gr++) {
		group_effect.param.push_back(nparam);

		Param par;
		par.name = "Group effect " + to_string(gr);
		par.type = GROUP_EFFECT;

		par.prior_type = NORMAL_FROM_SD_PRIOR;
		par.prior_sd_param = group_effect.sigma_param;
		par.prior_val1 = UNSET;
		par.prior_val2 = UNSET;

		param.push_back(par);
		nparam++;
	}
}


/// <akes proposals to the parameter sigma (standard deviation in the group effects)
void Model::propose_group_effect_sigma(vector <double> &param_value, double &prior, vector <Jump> &param_jump, const bool burnin, const Anneal &anneal) const {
	// cout << "Model::propose_group_effect_sigma()" << endl; // DEBUG
	timer[TIME_GROUP_SD].start();

	auto th = group_effect.sigma_param;

	auto &jump = param_jump[th];

	for (auto loop = 0; loop < group_effect_sigma_prop; loop++) {
		auto param_store = param_value[th];

		param_value[th] += normal_sample(0, jump.size);

		auto prior_propose = calculate_prior(param_value);
		auto al = exp(anneal.phi_Pr*(prior_propose - prior));

		jump.ntr++;
		if (MH_proposal(al,13)) {
			jump.nac++;
			prior = prior_propose;
			if (burnin == true)
				jump.size *= prop_up;
		} else {
			param_value[th] = param_store;
			if (burnin == true)
				jump.size *= prop_down;
		}
	}
	timer[TIME_GROUP_SD].stop();
}


/// Makes proposals to the group effects
void Model::propose_group_effect(const vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_inf_events, double &prior, vector <Jump> &param_jump, const bool burnin, const Anneal &anneal) const {
	// cout << "Model::propose_group_effect()" << endl; // DEBUG
	for (auto g = 0; g < ngroup; g++) {
		const auto &gr = group[g];

		timer[TIME_GROUP_EFFECT_INIT].start();
		auto precalc = set_precalc_group_effect(gr, ind_value, param_value);
		timer[TIME_GROUP_EFFECT_INIT].stop();

		timer[TIME_GROUP_EFFECT].start();

		auto th = group_effect.param[g];
		auto &jump = param_jump[th];

		for (auto loop = 0; loop < group_effect_prop; loop++) {
			auto param_store = param_value[th];

			// Makes a change to one of the covariance matrix parameters
			param_value[th] += normal_sample(0, jump.size);

			// Calculates the Metropolis-Hastings acceptance probability
			auto L_propose = likelihood_inf_events_fast(precalc, param_value[th]);
			auto prior_change = calculate_prior_change(th, param_store, param_value);
			auto al = exp(anneal.phi_L*(L_propose - L_inf_events[g]) + anneal.phi_Pr*prior_change);

			jump.ntr++;
			if (MH_proposal(al,20)) {
				jump.nac++;
				L_inf_events[g] = L_propose;
				prior += prior_change;
				if (burnin == true)
					jump.size *= prop_up;
			} else {
				param_value[th] = param_store;
				if (burnin == true)
					jump.size *= prop_down;
			}
		}
		timer[TIME_GROUP_EFFECT].stop();
	}
}


/// Pre-calculates quantities for fast calculation of the likelihood
PrecalcGroupEffect Model::set_precalc_group_effect(const Group &gr, const vector <IndValue> &ind_value, const vector <double> &param_value) const {
	// cout << "Model::set_precalc_group_effect()" << endl; // DEBUG
	PrecalcGroupEffect prec;

	prec.beta_fac = 0;
	prec.log_beta_num = 0;
	prec.con_fac = 0;

	auto I = 0.0;             // The total infectivity
	auto S = 0.0;             // The total susceptibility

	// Gets a sorted list of events for the entire group
	vector <Event> event;
	get_event_sequence(event, S, gr, ind_value);

	auto tmin = gr.inference_range.tmin;
	auto tmax = gr.inference_range.tmax;

	auto beta = param_value[beta_param];
	if (inf_model == FREQ_DEP)
		beta /= gr.nind;     // If frequency dependent divided by group size

	auto t = tmin;                                 // Sets a time variable
	for (const auto &ev : event) {
		const auto &ind = ind_value[ev.ind];
		auto tev = ev.time;
		auto tr = ev.trans;

		prec.beta_fac -= beta * S * I * (tev - t);
		t = tev;

		if (tr == 0) {
			if (tev > tmin) {
				prec.log_beta_num++;
				prec.con_fac += log(beta * ind.susceptibility * (I + EXT_FOI));
			}
			S -= ind.susceptibility;
		}

		I += ind.trans_infectivity_change[tr];
	}

	if (t > tmax)
		emsg("Likelihood error");
	prec.beta_fac -= beta * S * I * (tmax - t);

	return prec;
}


/// Fast calculation of the log likelihood for the transition events (using pre-calculated qautnties)
double Model::likelihood_inf_events_fast(const PrecalcGroupEffect &prec, const double value) const {
	// cout << "Model::likelihood_inf_events_fast()" << endl; // DEBUG
	auto group_factor = exp(value);
	return prec.con_fac + prec.log_beta_num * log(group_factor) + prec.beta_fac * group_factor;
}

