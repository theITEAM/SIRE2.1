/// This includes model functions related to infection events

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


/// This sets transition means, susceptibility and compartmental infectivity from fixed, SNP, and individual effects
void Model::set_individual_quantities(vector <IndValue> &ind_value, const vector <double> &param_value) const 
{
	for(auto g = 0u; g < ngroup; g++){
		set_individual_quantities_group(g,ind_value,param_value);
	}
}

void Model::set_individual_quantities_group(unsigned int g, vector <IndValue> &ind_value, const vector <double> &param_value) const 
{
	// cout << "Model::set_individual_quantities()" << endl; // DEBUG
	vector <double> fixed_effect_value(nfixed_effect);
	for (auto fe = 0; fe < nfixed_effect; fe++)
		fixed_effect_value[fe] = param_value[fixed_effect[fe].param];

	vector <double> snp_effect_value(nsnp_effect), snp_effect_dom(nsnp_effect);
	for (auto se = 0; se < nsnp_effect; se++) {
		const auto &snp_eff = snp_effect[se];
		snp_effect_value[se] = param_value[snp_eff.param_mag];
		snp_effect_dom[se] = snp_effect_value[se] * param_value[snp_eff.param_dom];
	}

	for(auto i : group[g].ind_ref){
		auto &ind = ind_value[i];
		const auto &indiv = individual[i];

		ind.trans_mean.resize(ntrans);
		for (auto tr = 0; tr < ntrans; tr++) {
			const auto &tra = trans[tr];

			auto sum = 0.0;
			for (auto ie : tra.ind_effect)
				sum += ind.ind_effect[ie];
			for (auto fe : tra.fixed_effect)
				sum += indiv.fixed_effect_X[fe] * fixed_effect_value[fe];
			
			for (auto se : tra.snp_effect) {
				switch (indiv.SNP_genotype[se]) {
					case AA:
						sum += snp_effect_value[se];
						break;
					case AB:
						sum += snp_effect_dom[se];
						break;
					case BB:
						sum -= snp_effect_value[se];
						break;
				}
			}

			switch (tra.type) {
				case INFECTION:
					ind.susceptibility = exp(sum);                        // This sets individual-based susceptibility
					ind.trans_mean[tr] = UNSET;
					break;

				case GAMMA: case EXP:
					ind.trans_mean[tr] = param_value[tra.mean_param] * exp(sum); // This sets individual-based transition mean
					break;
			}
		}

		auto sum = 0.0;
		for (auto ie : infectivity.ind_effect) sum += ind.ind_effect[ie];
		
		for (auto fe : infectivity.fixed_effect) sum += indiv.fixed_effect_X[fe]*fixed_effect_value[fe];
		
		for (auto se : infectivity.snp_effect) {
			switch (indiv.SNP_genotype[se]) {
				case AA:
					sum += snp_effect_value[se];
					break;
				case AB:
					sum += snp_effect_dom[se];
					break;
				case BB:
					sum -= snp_effect_value[se];
					break;
			}
		}
		ind.inf_single = exp(sum);
		
		/*
		ind.infectivity.resize(ncomp);
		for (auto c = 0; c < ncomp; c++) {                          // This sets individual-based infectivity
			auto inf = comp[c].infectivity;
			if (inf != 0) {
				auto sum = 0.0;   // TO BE REMOVED
				for (auto ie : comp[c].ind_effect)
					sum += ind.ind_effect[ie];
				for (auto fe : comp[c].fixed_effect)
					sum += indiv.fixed_effect_X[fe] * fixed_effect_value[fe];
				for (auto se : comp[c].snp_effect) {
					switch (indiv.SNP_genotype[se]) {
						case AA:
							sum += snp_effect_value[se];
							break;
						case AB:
							sum += snp_effect_dom[se];
							break;
						case BB:
							sum -= snp_effect_value[se];
							break;
					}
				}

				if (sum != 0) inf *= exp(sum);
				
				inf *= ind.inf_single;
			}
			ind.infectivity[c] = inf;
		}
		*/

		ind.trans_infectivity_change.resize(ntrans);                 // This sets change in infectivity from a transition
		for (auto tr = 0; tr < ntrans; tr++){
			ind.trans_infectivity_change[tr] = ind.inf_single*(comp[trans[tr].to].infectivity - comp[trans[tr].from].infectivity);
		}
		
		if (false) {
			cout << ind.susceptibility << " ";
			cout << "   Transition mean:";
			for (auto tr = 1u; tr < ntrans; tr++)
				cout << ind.trans_mean[tr] << " ";
			cout << "   Transition infectivity change: ";
			for (auto tr = 0; tr < ntrans; tr++)
				cout << ind.trans_infectivity_change[tr] << " ";
			cout << endl;
		}
	}
}


/// Calculates the log likelihood for the events and returns a vector
vector <double> Model::calculate_L_inf_events(const vector <IndValue> &ind_value, const vector <double> &param_value) const {
	// cout << "Model::calculate_L_inf_events()" << endl; // DEBUG
	vector <double> L_inf_events(ngroup);
	for (auto g = 0; g < ngroup; g++)
		L_inf_events[g] = likelihood_inf_events(group[g], ind_value, param_value);

	return L_inf_events;
}


/// Used for time-ordering events
bool EV_ord(Event ev1, Event ev2) {
	return (ev1.time < ev2.time);
};


/// Constructs an event sequence for a given group
void Model::get_event_sequence(vector <Event> &event, double &S, const Group &gr, const vector <IndValue> &ind_value) const {
	auto tmax = gr.inference_range.tmax;

	for (auto i : gr.ind_ref) {
		const auto &ind = ind_value[i];
		if (ind.infected == true) {
			for (auto tr = 0; tr < ntrans; tr++) {
				auto t = ind.trans_time[tr];
				if (t < tmax) {
					Event ev;
					ev.ind = ind.index;
					ev.trans = tr;
					ev.time = t;
					event.push_back(ev);
				}
			}
		}
		S += ind.susceptibility;
	}

	sort(event.begin(), event.end(), EV_ord);
}


/// Calculates the log likelihood for the transition events
double Model::likelihood_inf_events(const Group &gr, const vector <IndValue> &ind_value, const vector <double> &param_value) const {
	// cout << "Model::likelihood_inf_events()" << endl; // DEBUG
	auto I = 0.0;                                                    // The total infectivity
	auto S = 0.0; 	                                                   // The total susceptibility

	// Gererates a sorted list of events for the entire group
	vector <Event> event;
	get_event_sequence(event, S, gr, ind_value);

	auto tmin = gr.inference_range.tmin;
	auto tmax = gr.inference_range.tmax;

	auto beta = param_value[beta_param];
	if (inf_model == FREQ_DEP)
		beta /= gr.nind;                        // If frequency dependent divided by group size
	if (group_effect.on == true)                                      // Adds in the group effect
		beta *= exp(param_value[group_effect.param[gr.index]]);

	auto L = 0.0;                                                     // The log of the likelihood

	auto t = tmin;                                                    // Sets a time variable
	for (const auto &ev : event) {                                    // Goes through each event and caluculates likelhood
		const auto &ind = ind_value[ev.ind];
		auto tev = ev.time;
		auto tr = ev.trans;

		L -= beta * S * I * (tev - t);
		t = tev;

		if (tr == 0) {
			if (tev > tmin)
				L += log(beta * ind.susceptibility * (I + EXT_FOI));
			S -= ind.susceptibility;
		}

		I += ind.trans_infectivity_change[tr];
	}

	if (t > tmax)
		emsg("Likelihood error");
	L -= beta * S * I * (tmax - t);

	return L;
}


/// This resamples the transmission rate parameter
void Model::propose_transmission_rate(const vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_inf_events, double &prior, vector <Jump> &param_jump, const bool burnin, const Quench &quench) const {
	if(param[beta_param].prior_type == FIXED_PRIOR) return;

	// cout << "Model::propose_transmission_rate()" << endl; // DEBUG
	timer[TIME_TRANS_RATE].start();

	if (trans[0].type != INFECTION)
		emsg("Transition should be infection");

	auto beta = param_value[beta_param];

	vector <PrecalcInfParam> precalc(ngroup);

	auto beta_fac_total = 0.0;
	auto log_beta_num_total = 0.0;
	auto con_fac_total = 0.0;
	for (auto g = 0; g < ngroup; g++) {
		auto &prec = precalc[g];
		const auto &gr = group[g];

		auto fac = 1.0;
		if (inf_model == FREQ_DEP)
			fac /= gr.nind;         // Add frequency dependence
		if (group_effect.on == true)                                      // Adds in the group effect
			fac *= exp(param_value[group_effect.param[gr.index]]);

		prec.beta_fac = 0;
		prec.log_beta_num = 0;
		prec.con_fac = 0;

		auto I = 0.0;                                                     // The total infectivity
		auto S = 0.0;                                                     // The total susceptibility

		vector <Event> event;
		get_event_sequence(event, S, gr, ind_value);                      // Gererates a sorted list of events for the entire group

		auto tmin = gr.inference_range.tmin;
		auto tmax = gr.inference_range.tmax;

		auto t = tmin;                                                    // Sets a time variable
		for (const auto &ev : event) {
			const auto &ind = ind_value[ev.ind];
			auto tev = ev.time;
			auto tr = ev.trans;

			prec.beta_fac -= fac * S * I * (tev - t);
			t = tev;

			if (tr == 0) {
				if (tev > tmin) {
					prec.log_beta_num++;
					prec.con_fac += log(fac * ind.susceptibility * (I + EXT_FOI));
				}
				S -= ind.susceptibility;
			}

			I += ind.trans_infectivity_change[tr];
		}

		prec.beta_fac -= fac * S * I * (tmax - t);

		beta_fac_total += prec.beta_fac;
		log_beta_num_total += prec.log_beta_num;
		con_fac_total += prec.con_fac;
	}
	
	if(false){
		auto beta = param_value[beta_param];
		auto Li = 0.0;
		for (auto g = 0; g < ngroup; g++) Li += L_inf_events[g];
		auto Lp = con_fac_total + log_beta_num_total * log(beta) + beta_fac_total * beta;
		cout << Li << " " << Lp << " " << Li-Lp << " check L\n";
		emsg("done");
	}
	
	auto phi = quench.phi_L; 
	
	if(phi < 0.001 || ran() < 0.1){                     // This performs a random walk proposal
		auto th = beta_param;
		auto &jump = param_jump[th];
	
		// Makes a change to beta
		auto beta_prop = beta + normal_sample(0, jump.size);
		param_value[th] = beta_prop;

		if(inbounds(param_value) == false){
			param_value[th] = beta;
		}
		else{			
			// Calculates the Metropolis-Hastings acceptance probability
			auto prior_change = calculate_prior_change(th, beta, param_value);
			
			auto dL = 0.0;
			vector <double> L_inf_events_prop(ngroup,0);
			for (auto g = 0; g < ngroup; g++) {
				auto &prec = precalc[g];
				L_inf_events_prop[g] = prec.con_fac + prec.log_beta_num * log(beta_prop) + prec.beta_fac * beta_prop;
				dL += L_inf_events_prop[g] - L_inf_events[g];
			}
		
			auto inside = quench.phi_L*dL + quench.phi_Pr*prior_change;
			if(inside > 100.0) inside = 100;
				
			auto al = exp(inside);
			if (prior_change == ZERO_PRIOR) al = 0;

			jump.ntr++;
			if (MH_proposal(al,46)) {
				jump.nac++;
				L_inf_events = L_inf_events_prop;
				
				prior += prior_change;
				if (burnin == true) jump.size *= prop_up;
			}
			else {
				param_value[th] = beta;
				if (burnin == true) jump.size *= prop_down;
			}
		}
	}
	else{                                      // This proposal samples from a gamma distribution
		/// Beta is gamma distributed, so we sample from this distribution
		auto shape = 1 + log_beta_num_total;
		auto mean = -shape / beta_fac_total;
	
		if(phi != 1){
			auto shape_new = 1+(shape-1)*phi;
			auto mean_new = mean*shape_new/(shape*phi);
			mean = mean_new; shape = shape_new;
		}
		
		param_value[beta_param] = gamma_sample(mean, shape);

		auto prior_change = calculate_prior_change(beta_param, beta, param_value);

		auto al = exp(quench.phi_Pr*prior_change);
	
		auto &jump = param_jump[beta_param];

		if (false) { // This checks that Gibbs sampling is being performed correctly
			auto beta_prop = param_value[beta_param];
			auto Li = 0.0;
			for (auto g = 0; g < ngroup; g++)
				Li += L_inf_events[g];
			auto Lp = con_fac_total + log_beta_num_total * log(beta_prop) + beta_fac_total * beta_prop;
			auto probfi = gamma_probability(beta, mean, shape);
			auto probif = gamma_probability(param_value[beta_param], mean, shape);
			cout << Lp-Li << " " 
					 << quench.phi_L*(Lp-Li) + probfi-probif + quench.phi_Pr*prior_change  
					 << " check" << endl;
			//emsg("do");
		}

		jump.ntr++;
		if (MH_proposal(al,1)) {
			jump.nac++;

			// This updates the likelihoods in each of the groups
			auto beta_prop = param_value[beta_param];
			for (auto g = 0; g < ngroup; g++) {
				auto &prec = precalc[g];
				L_inf_events[g] = prec.con_fac + prec.log_beta_num * log(beta_prop) + prec.beta_fac * beta_prop;
			}
			prior += prior_change;
		} else
			param_value[beta_param] = beta;
	}
	
	timer[TIME_TRANS_RATE].stop();
}


/// Sets the events for the initial state on the chain
void Model::set_initial_events(vector<IndValue> &ind_value, const vector <double> &param_value) const {
	// cout << "Model::set_initial_events()" << endl; // DEBUG
	for (auto i = 0; i < N; i++) {
		ind_value[i].trans_time.resize(ntrans);

		auto loop = 0;
		do {
			auto ev_samp = SAMPLE_ONLY;
			if (loop > 0.5 * INITIAL_TRANS_SAMPLE)
				ev_samp = SAMPLE_MIX;
			if (loop > 0.8 * INITIAL_TRANS_SAMPLE)
				ev_samp = SAMPLE_RANGE;

			if (initial_event_sample(ind_value[i], param_value, ev_samp) != -LARGE)
				break;
			loop++;
		} while (loop < INITIAL_TRANS_SAMPLE);

		if (loop == INITIAL_TRANS_SAMPLE)
			emsg("Cannot find initial sequence");
	}
}


/// Samples an initial possible set of events for each individual
double Model::initial_event_sample(IndValue &ind, const vector <double> &param_value, EvInitSampType ev_samp) const {
	// cout << "Model::event_sample()" << endl; // DEBUG
	auto probif = 0.0;

	const auto &indiv = individual[ind.index];

	switch (indiv.status) {
		case INFECTED:
			ind.infected = true;
			for (auto tr = 0; tr < ntrans; tr++) {
				auto &tr_time = ind.trans_time[tr];
				const auto &trange = indiv.trans_time_range[tr];
				auto tmin = trange.tmin, tmax = trange.tmax;
				
				 // Sets the initial transition time if set in 'data_column_initial'
				if(indiv.trans_time_initial[tr] != UNSET){           
					tmin = indiv.trans_time_initial[tr]; tmax = tmin;
				}
				
				for(auto tr2 = tr+1; tr2 < ntrans; tr2++){
					auto tr_tinit = indiv.trans_time_initial[tr2];
					if(tr_tinit != UNSET && tr_tinit < tmax) tmax = tr_tinit;
				}
				
				if (tmin != tmax) {
					switch (trans[tr].type) {
						case INFECTION:
							tr_time = tmin + ran() * (tmax - tmin);
							break;

						case GAMMA: case EXP: {
							auto mean = ind.trans_mean[tr];
							auto shape = 1.0; if(trans[tr].type == GAMMA) shape = param_value[trans[tr].shape_param];
							auto dt = gamma_sample(mean, shape);

							switch (ev_samp) {
								case SAMPLE_ONLY:
									tr_time = ind.trans_time[tr - 1] + dt;
									break;

								case SAMPLE_MIX:
									if (ran() < 0.5)
										tr_time = tmin + dt;
									else
										tr_time = ind.trans_time[tr - 1] + dt;
									break;

								case SAMPLE_RANGE:
									tr_time = tmin + ran() * (tmax - tmin);
									break;
							}

							if (tr_time < ind.trans_time[tr - 1])
								return -LARGE;
							if (tr_time <= tmin || tr_time >= tmax)
								return -LARGE;
							probif += gamma_probability(dt, mean, shape);
						}
						break;

						default:
							emsg("Not implemented");
							break;
					}
				} else
					tr_time = tmin;
			}

			break;

		case NOT_INFECTED:
			ind.infected = false;
			for (auto tr = 0; tr < ntrans; tr++)
				ind.trans_time[tr] = UNSET;
			break;

		case UNKNOWN:
			ind.infected = false;
			for (auto tr = 0; tr < ntrans; tr++)
				ind.trans_time[tr] = UNSET;
			break;
	}

	return probif;
}


/// This makes changes to event times
void Model::propose_event_times(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_inf_events, double &L_trans_events, vector <double> &L_diag_test, vector <Jump> &event_jump, const bool burnin, const Quench &quench) const 
{
	// cout << "Model::propose_event_times()" << endl; // DEBUG
	timer[TIME_EVENTS].start();
	for (const auto &gr : group) {
		if (gr.nev_change > 0) {
			auto g = gr.index;
			auto &jump = event_jump[g];

			vector <TimeChange> time_change;                      // Keeps track of any changes made to event times

			auto probif = 0.0, probfi = 0.0;
			auto L_trans_events_prop = L_trans_events;

			auto loop_max = int(jump.size);                       // The number of event changes per proposal is dynamically changed
			if (loop_max == 0)
				loop_max = 1;
	
			for (auto loop = 0; loop < loop_max; loop++) {
				const auto &ev_samp = gr.ev_change[int(ran() * gr.nev_change)];

				auto i = ev_samp.ind;
				auto tr = ev_samp.trans;
				auto type = ev_samp.ev_samp_type;

				auto &ind = ind_value[i];
				if (ind.infected == true) {
					const auto &indiv = individual[i];
					const auto &tran = indiv.trans_time_range[tr];

					auto t = ind.trans_time[tr];

					double dprobif, dprobfi;
					double t_prop;

					switch (type) {
						case PEAKED_SAMPLE:
							if (tr == 0) {                        // This uses the next transition time and samples backwards in time
								auto tr_next = tr + 1;
								switch (trans[tr_next].type) {
									case GAMMA: case EXP: {
										auto shape = 1.0; if(trans[tr_next].type == GAMMA) shape = param_value[trans[tr_next].shape_param];
										auto mean = ind.trans_mean[tr_next];

										auto t_next = ind.trans_time[tr_next];
										auto dt = t_next - t;
										dprobfi = gamma_probability(dt, mean, shape);

										auto dt_prop = gamma_sample(mean, shape);
										dprobif = gamma_probability(dt_prop, mean, shape);

										t_prop = t_next - dt_prop;
									}
									break;

									default:
										emsg("Not supported");
										break;
								}
							} else {                                // This uses the last transition time and samples forward in time
								if (tr == ntrans - 1) {
									switch (trans[tr].type) {
										case GAMMA: case EXP: {
											auto shape = 1.0; if(trans[tr].type == GAMMA) shape = param_value[trans[tr].shape_param];
											auto mean = ind.trans_mean[tr];

											auto t_last = ind.trans_time[tr - 1];
											auto dt = t - t_last;
											dprobfi = gamma_probability(dt, mean, shape);

											auto dt_prop = gamma_sample(mean, shape);
											dprobif = gamma_probability(dt_prop, mean, shape);
											t_prop = t_last + dt_prop;
										}
										break;

										default:
											emsg("Not supported");
											break;
									}
								} else {                                    // This is for a trasition time between two transitions
									auto tr_next = tr + 1;
									if (!(trans[tr].type == GAMMA || trans[tr].type == EXP)
										|| !(trans[tr + 1].type == GAMMA || trans[tr + 1].type == EXP)){
										emsg("Problem");
									}
										
									auto mean1 = ind.trans_time[tr - 1] + ind.trans_mean[tr];
									auto sh1 = 1.0; if(trans[tr].type == GAMMA) sh1 = param_value[trans[tr].shape_param];
									auto var1 = mean1 * mean1 / sh1;

									auto mean2 = ind.trans_time[tr + 1] - ind.trans_mean[tr + 1];
									auto sh2 = 1.0; if(trans[tr+1].type == GAMMA) sh2 = param_value[trans[tr + 1].shape_param];
									auto var2 = mean2 * mean2 / sh2;

									auto sd = sqrt(var1 * var2 / (var1 + var2)); // A normal approximation is made for the gamma distributions
									auto mean = (mean1 * var2 + mean2 * var1) / (var1 + var2);

									dprobfi = normal_probability(t, mean, sd);

									t_prop = normal_sample(mean, sd);
									dprobif = normal_probability(t_prop, mean, sd);
								}
							}
							break;

						case UNIFORM_SAMPLE:                            // A transtion time is randomly selected from the possible range
							auto tmin = tran.tmin;
							auto tmax = tran.tmax;
							if (tr > 0) {
								auto t_last = ind.trans_time[tr - 1];
								if (t_last > tmin)
									tmin = t_last;
							}

							if (tr < ntrans - 1) {
								auto t_next = ind.trans_time[tr + 1];
								if (t_next < tmax)
									tmax = t_next;
							}

							t_prop = tmin + ran() * (tmax - tmin);
							dprobif = 0;
							dprobfi = 0;
							break;
					}

					auto fail = false;
					if (t_prop <= tran.tmin || t_prop >= tran.tmax)
						fail = true;  // The resample fails because out of time range
					if (tr > 0) {
						if (t_prop <= ind.trans_time[tr - 1])
							fail = true;   // Fails because time ordering incorrect
					}
					if (tr < ntrans - 1) {
						if (t_prop >= ind.trans_time[tr + 1])
							fail = true;
					}

					jump.nevent_tr++;
					if (fail == true)
						jump.nevent_fa++;
					else {
						if (tr > 0) {                                               // Calculates the change in transition likelihood
							const auto &tra = trans[tr];
							switch (tra.type) {
								case GAMMA: case EXP: {
									auto shape = 1.0; if(tra.type == GAMMA) shape = param_value[tra.shape_param];
									auto mean = ind.trans_mean[tr];
									auto t_last = ind.trans_time[tr - 1];
									L_trans_events_prop += gamma_probability(t_prop - t_last, mean, shape) - gamma_probability(t - t_last, mean, shape);
								}
								break;

								default:
									emsg("Not supported");
									break;
							}
						}

						if (tr < ntrans - 1) {
							const auto &tra = trans[tr + 1];
							switch (tra.type) {
								case GAMMA: case EXP: {
									auto shape = 1.0; if(tra.type == GAMMA) shape = param_value[tra.shape_param];
									auto mean = ind.trans_mean[tr + 1];
									auto t_next = ind.trans_time[tr + 1];
									L_trans_events_prop += gamma_probability(t_next - t_prop, mean, shape) - gamma_probability(t_next - t, mean, shape);
								}
								break;

								default:
									emsg("Not supported");
									break;
							}
						}

						TimeChange tch;
						tch.ind = i;
						tch.trans = tr;
						tch.time = t;
						time_change.push_back(tch);

						ind.trans_time[tr] = t_prop;

						probif += dprobif;
						probfi += dprobfi;
					}
				}
			}

			timer[TIME_EVENTS_LIKE].start();
			auto L_inf_events_prop = likelihood_inf_events(gr, ind_value, param_value);
			auto L_diag_test_prop = likelihood_diag_test(gr, ind_value, param_value);
			timer[TIME_EVENTS_LIKE].stop();

			if (false)
				cout << L_trans_events_prop - L_trans_events + probfi - probif << " shoudl be zero" << endl;

			auto al = exp(quench.phi_L*(L_inf_events_prop - L_inf_events[g]) +  quench.phi_DT*(L_diag_test_prop - L_diag_test[g]) + quench.phi_L*(L_trans_events_prop - L_trans_events) + probfi - probif);
			
if(std::isnan(al)) cout << L_inf_events_prop << " " << L_inf_events[g] << " " << L_diag_test_prop  << " " << L_diag_test[g] << " " << L_trans_events_prop << " " << L_trans_events << " hh\n";
			jump.ntr++;
			if (MH_proposal(al,2)) {
				jump.nac++;
				L_inf_events[g] = L_inf_events_prop;
				L_trans_events = L_trans_events_prop;
				L_diag_test[g] = L_diag_test_prop;
				if (burnin == true) {
					jump.size *= prop_up;
					if (jump.size > gr.nev_change)
						jump.size = gr.nev_change;
				}
			} else {
				for (int j = time_change.size() - 1; j >= 0; j--) {        // Resets the event times
					const auto &tch = time_change[j];

					ind_value[tch.ind].trans_time[tch.trans] = tch.time;
				}
				if (burnin == true)
					jump.size *= prop_down;
			}
		}
	}
	timer[TIME_EVENTS].stop();
}


/// Samples a time from an infection sampler
double InfSampler::sample_time() {
	auto z = ran();
	auto b = 0;
	while (b < nbin && z > prob_sum[b])
		b++;

	return tmin + (b + ran()) * (tmax - tmin) / nbin;
}


/// Gives the probability of sampling a certain infection time from the sampler
double InfSampler::sample_prob(double t) {
	auto b = int(nbin * (t - tmin) / (tmax - tmin));
	if (b < 0 || b >= nbin)
		emsg("Sampler out of range");
	return log_prob[b];
}


/// Makes proposals to add and remove infected individuals (whose infection status is unknown)
/// NOTE: This could be improved by doing a better job samping the initial infection time
void Model::propose_add_rem(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_inf_events, double &L_trans_events, vector <double> &L_diag_test, vector <Jump> &add_rem_jump, const bool burnin, const Quench &quench) const {
	// cout << "Model::propose_add_rem()" << endl; // DEBUG
	timer[TIME_ADD_REM].start();

	for (const auto &gr : group) {
		if (gr.nunknown > 0) {
			auto g = gr.index;
			auto &jump = add_rem_jump[g];

			auto loop_max = int(jump.size);
			if (loop_max == 0)
				loop_max = 1;

			vector <double> trans_time(ntrans);
			for (auto loop = 0; loop < loop_max; loop++) {
				vector <int> infected, uninfected;
				for (auto i : gr.unknown) {
					if (ind_value[i].infected == true)
						infected.push_back(i);
					else
						uninfected.push_back(i);
				}
				auto ninfected = infected.size(), nuninfected = uninfected.size();

				auto probif = 0.0, probfi = 0.0;
				auto L_trans_events_prop = L_trans_events;
				vector <AddRemChange> addrem_change;

				jump.nevent_tr++;
				if (ran() < 0.5) {                                           // Add a new infected
					if (nuninfected > 0) {
						auto j = int(ran() * nuninfected);                       // Randomly select from the list

						auto i = uninfected[j];
						auto &ind = ind_value[i];
						const auto &indiv = individual[i];
						auto &timer = indiv.trans_time_range;

						auto t_inf = ind.inf_sampler.sample_time();              // Samples infection time from the infection sampler

						bool fail = false;

						auto dL = 0.0, dprobif = 0.0;

						auto t = t_inf;
						trans_time[0] = t;
						for (auto tr = 1u; tr < ntrans; tr++) {                  // Generates other event times by samping from model
							const auto &tra = trans[tr];
							if (tra.type == GAMMA || tra.type == EXP) {
								auto shape = 1.0; if(tra.type == GAMMA) shape = param_value[tra.shape_param];
								auto mean = ind.trans_mean[tr];
								auto t_next = t + gamma_sample(mean, shape);
								if (t_next < timer[tr].tmin || t_next > timer[tr].tmax) {
									fail = true;
									break;
								}

								trans_time[tr] = t_next;

								auto prob = gamma_probability(t_next - t, mean, shape);
								dL += prob;
								dprobif += prob;
								t = t_next;
							}
						}

						if (fail == false) {
							L_trans_events_prop += dL;
							probif += dprobif;

							probif += log(1.0 / nuninfected);
							probif += ind.inf_sampler.sample_prob(t_inf);

							uninfected.erase(uninfected.begin() + j);
							nuninfected--;
							infected.push_back(i);
							ninfected++;

							AddRemChange arch;
							arch.ind = i;
							arch.infected = ind.infected;
							arch.trans_time = ind.trans_time;
							addrem_change.push_back(arch);

							ind.infected = true;
							ind.trans_time = trans_time;
							probfi += log(1.0 / ninfected);
						}
					} else
						jump.nevent_fa++;
				} else {                                              // Remove an infected individual
					if (ninfected > 0) {
						auto j = int(ran() * ninfected);                  // Randomly select infected individual

						auto i = infected[j];
						auto &ind = ind_value[i];
						const auto &indiv = individual[i];
						auto &timer = indiv.trans_time_range;

						probif += log(1.0 / ninfected);

						infected.erase(infected.begin() + j);
						ninfected--;
						uninfected.push_back(i);
						nuninfected++;

						probfi += log(1.0 / nuninfected);                    // Calculates reverse proposal probability
						probfi += ind.inf_sampler.sample_prob(ind.trans_time[0]);

						for (auto tr = 1u; tr < ntrans; tr++) {
							const auto &tra = trans[tr];
							if (tra.type == GAMMA || tra.type == EXP) {
								auto shape = 1.0; if(tra.type == GAMMA) shape = param_value[tra.shape_param];
								auto mean = ind.trans_mean[tr];
								auto prob = gamma_probability(ind.trans_time[tr] - ind.trans_time[tr - 1], mean, shape);
								L_trans_events_prop -= prob;
								probfi += prob;
							}
						}

						AddRemChange arch;
						arch.ind = i;
						arch.infected = ind.infected;
						arch.trans_time = ind.trans_time;
						addrem_change.push_back(arch);

						for (auto tr = 0; tr < ntrans; tr++)
							trans_time[tr] = UNSET;

						ind.infected = false;
						ind.trans_time = trans_time;
					}
				}

				if (false) {                                          // Checks lists are updated correctly
					if (ninfected != infected.size())
						emsg("error1");
					if (nuninfected != uninfected.size())
						emsg("error2");
					if (ninfected + nuninfected != gr.nunknown)
						emsg("error3");

					for (auto i : gr.unknown) {
						if (ind_value[i].infected == true) {
							if (find_in(infected, i) == UNSET)
								emsg("error4");
						} else {
							if (find_in(uninfected, i) == UNSET)
								emsg("error5");
						}
					}
				}

				if (addrem_change.size() > 0) {
					timer[TIME_ADD_REM_LIKE].start();
					auto L_inf_events_prop = likelihood_inf_events(gr, ind_value, param_value);
					auto L_diag_test_prop = likelihood_diag_test(gr, ind_value, param_value);
					timer[TIME_ADD_REM_LIKE].stop();

					auto al = exp(quench.phi_L*(L_inf_events_prop - L_inf_events[g] + L_trans_events_prop - L_trans_events) +  quench.phi_DT*(L_diag_test_prop - L_diag_test[g])
					               + probfi - probif);

					jump.ntr++;
					if (MH_proposal(al,3)) {
						jump.nac++;
						L_inf_events[g] = L_inf_events_prop;
						L_trans_events = L_trans_events_prop;
						L_diag_test[g] = L_diag_test_prop;
						if (burnin == true) {
							jump.size *= prop_up;
							if (jump.size > gr.nunknown * 2)
								jump.size = gr.nunknown * 2;
						}
					} else {
						for (int j = addrem_change.size() - 1; j >= 0; j--) {
							const auto &arch = addrem_change[j];
							ind_value[arch.ind].infected = arch.infected;
							ind_value[arch.ind].trans_time = arch.trans_time;
						}
						if (burnin == true)
							jump.size *= prop_down;
					}
				}
			}
		}
	}

	timer[TIME_ADD_REM].stop();
}
