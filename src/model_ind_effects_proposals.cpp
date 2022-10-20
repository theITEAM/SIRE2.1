/// This stores functions giving proposals on individual effects

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "model.hpp"
#include "timers.hpp"
#include "utils.hpp"


/// This makes proposals to susceptibility individual effects
void Model::propose_susceptibility_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_inf_events, vector <Jump> &ind_effect_jump, const Anneal &anneal) const {
	// cout << "Model::propose_susceptibility_ind_effects()" << endl; // DEBUG
	const auto &ind_eff = trans[0].ind_effect;
	auto nie = ind_eff.size();
	if (nie == 0)
		return;

	timer[TIME_SUS_IND_EFFECT_INIT].start();
	auto inv_cov_matrix = calculate_inv_cov_matrix(param_value);
	auto precalc = set_precalc_likelihood_sus(ind_value, param_value);
	timer[TIME_SUS_IND_EFFECT_INIT].stop();

	timer[TIME_SUS_IND_EFFECT].start();
	vector <PrecalcIndEff> prec_ind_eff(nie);

	for (auto i = 0; i < N; i++) {
		auto &prec = precalc[i];
		auto &ind = ind_value[i];
		auto &indiv = individual[i];
		auto g = indiv.group;

		for (auto j = 0; j < nie; j++)
			set_precalc_ind_eff(i, prec_ind_eff[j], ind_eff[j], ind_value, inv_cov_matrix);

		auto sus_fac = log(ind.susceptibility);
		for (auto loop = 0; loop < ind_effect_sus_prop; loop++) {
			for (auto j = 0; j < nie; j++) {
				const auto &prec_ie = prec_ind_eff[j];
				auto ie = ind_eff[j];

				auto sd = prec_ie.sd, mean = prec_ie.mean;
				if(anneal.phi_IE != 1) sd /= sqrt(anneal.phi_IE);
				
				auto ind_eff_store = ind.ind_effect[ie];
				auto ind_eff_prop = normal_sample(mean, sd);

				auto sus_fac_store = sus_fac;
				sus_fac += ind_eff_prop - ind_eff_store;

				auto dL_inf_ev = Li_change_sus(sus_fac_store, sus_fac, prec);
				
				auto inside = anneal.phi_L*dL_inf_ev; if(inside > 100) inside = 100;
				
				auto al = exp(inside);
	
				if (false) { // Checks for correct sampling from individual effect distribution
					auto probfi = normal_probability(ind_eff_store,mean,sd);
					auto probif = normal_probability(ind_eff_prop, mean,sd);
					cout << probfi - probif + anneal.phi_IE*Li_change_ind_eff(ind_eff_store, ind_eff_prop, prec_ie) << "zero" << " " << ind_eff_store << " " <<  ind_eff_prop << endl;
				}
					
				if(indiv.group == UNSET) check_one(al,1);

				ind_effect_jump[ie].ntr++;
				if (MH_proposal(al,5)) {
					ind_effect_jump[ie].nac++;
					ind.ind_effect[ie] = ind_eff_prop;

					auto cv = ind_effect[ie].covar_ref;
					L_ind_effect[cv] += Li_change_ind_eff(ind_eff_store, ind_eff_prop, prec_ie);
					if (g != UNSET)
						L_inf_events[g] += dL_inf_ev;
				} else
					sus_fac = sus_fac_store;
			}
		}

		ind.susceptibility = exp(sus_fac);
	}
	timer[TIME_SUS_IND_EFFECT].stop();
}


/// This pre-calculates quantities so likelihood can be calculated from susceptibility
vector <PrecalcLiSus> Model::set_precalc_likelihood_sus(const vector <IndValue> &ind_value, const vector <double> &param_value) const {
	// cout << "Model::set_precalc_likelihood_sus()" << endl; // DEBUG
	vector <PrecalcLiSus> prec(N);
	for (auto &pre : prec) {
		pre.beta_fac = 0;
		pre.log_beta_num = 0;
	}

	vector <int> infected(N);
	for (auto &inf : infected)
		inf = 0;

	for (const auto &gr : group) {
		auto I = 0.0;                                           // The total infectivity
		auto S = 0.0;                                           // The total susceptibility

		auto sus_list = gr.ind_ref;                             // Makes a list of susceptible individuals

		vector <Event> event;
		get_event_sequence(event, S, gr, ind_value);            // Gets a sorted list of events for the entire group

		auto tmin = gr.inference_range.tmin;
		auto tmax = gr.inference_range.tmax;

		auto beta = param_value[beta_param];
		if (inf_model == FREQ_DEP)
			beta /= gr.nind;              // If frequency dependent divided by group size
		if (group_effect.on == true)                            // Adds in the group effect
			beta *= exp(param_value[group_effect.param[gr.index]]);

		auto t = tmin;                                          // Sets a time variable
		for (const auto &ev : event) {
			const auto &ind = ind_value[ev.ind];
			auto tev = ev.time;
			auto tr = ev.trans;

			auto val = -beta * I * (tev - t);
			for (auto i : sus_list)
				if (infected[i] == 0)
					prec[i].beta_fac += val;
			t = tev;

			if (tr == 0) {
				auto i = ev.ind;
				if (tev > tmin)
					prec[i].log_beta_num++;
				infected[i] = 1;
			}

			I += ind.trans_infectivity_change[tr];
		}

		auto val = -beta * I * (tmax - t);
		for (auto i : sus_list)
			if (infected[i] == false)
				prec[i].beta_fac += val;
	}

	return prec;
}


/// Calculates the change in event likelihood as a result of susceptibility changing
double Model::Li_change_sus(const double sus_before, const double sus_after, PrecalcLiSus &prec) const {
	// cout << "Model::Li_change_sus()" << endl; // DEBUG
	return prec.log_beta_num * (sus_after - sus_before) + prec.beta_fac * (exp(sus_after) - exp(sus_before));
}


/// This pre-calculates quantities so changes in individual effect likelihood is fast
void Model::set_precalc_ind_eff(const int i, PrecalcIndEff &prec, int ie, const vector <IndValue> &ind_value, const vector <InvCovMat> &inv_cov_matrix) const {
	// cout << "Model::set_precalc_ind_eff()" << endl; // DEBUG
	auto val_sq_fac = 0.0, val_fac = 0.0;

	const auto cv = ind_effect[ie].covar_ref;
	const auto &cov = covariance[cv];
	const auto e_val = ind_effect[ie].covar_num;
	const auto &M = inv_cov_matrix[cv].M;
	const auto &mat = matrix[cov.matrix];
	auto E = cov.E;

	const auto &ie_ref = cov.ind_effect_ref;

	auto val = mat.Ainvdiag[i];
	val_sq_fac -= 0.5 * val * M[e_val][e_val];

	const auto &ind_eff = ind_value[i].ind_effect;
	auto sum = 0.0;
	for (auto ei = 0; ei < E; ei++) {
		if (ei != e_val)
			sum += M[e_val][ei] * ind_eff[ie_ref[ei]];
	}
	val_fac -= val * sum;

	for (const auto &ele : mat.Ainvlist[i]) {
		const auto &ind_eff = ind_value[ele.i].ind_effect;
		auto sum = 0.0;
		for (auto ei = 0; ei < E; ei++)
			sum += M[e_val][ei] * ind_eff[ie_ref[ei]];
		val_fac -= ele.val * sum;
	}

	prec.val_sq_fac = val_sq_fac;
	prec.val_fac = val_fac;

	// Calculates the mean and standard deviation in the distribution
	auto var = -1.0 / (2 * val_sq_fac);
	prec.sd = sqrt(var);
	prec.mean = val_fac * var;
}


/// Calculates the change in individual effect likelihood
double Model::Li_change_ind_eff(const double before, const double after, const PrecalcIndEff &prec) const {
	// cout << "Model::Li_change_ind_eff()" << endl; // DEBUG
	return prec.val_sq_fac * (after * after - before * before) +  prec.val_fac * (after - before);
}


/// This makes proposals to transition individual effects
void Model::propose_trans_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, double &L_trans_events, vector <Jump> &ind_effect_jump, const Anneal &anneal) const {
	// cout << "Model::propose_trans_ind_effects()" << endl; // DEBUG
	// Checks if individual effects exist
	auto exist = false;
	for (auto tr = 1u; tr < ntrans; tr++) {
		if (trans[tr].ind_effect.size() > 0)
			exist = true;
	}
	if (exist == false)
		return;

	timer[TIME_TRANS_IND_EFFECT].start();
	auto inv_cov_matrix = calculate_inv_cov_matrix(param_value);

	// This changes individual effect for a given transition

	for (auto tr = 1u; tr < ntrans; tr++) {
		const auto &tra = trans[tr];
		const auto &ind_eff = tra.ind_effect;

		if (tra.all_ind_effect_proposal == true) {
			auto nie = ind_eff.size();
			vector <PrecalcIndEff> prec_ind_eff(nie);

			auto mean_av = param_value[tra.mean_param];
			auto shape = 1; if(tra.type == GAMMA) shape = param_value[tra.shape_param];  // Loads

			double L, L_prop, dt;

			for (auto i = 0; i < N; i++) {
				auto &ind = ind_value[i];
				auto &indiv = individual[i];

				for (auto j = 0; j < nie; j++)
					set_precalc_ind_eff(i, prec_ind_eff[j], ind_eff[j], ind_value, inv_cov_matrix);

				auto mean = ind.trans_mean[tr];
				auto mean_fac = log(mean / mean_av);
				if (ind.infected == true) {
					dt = ind.trans_time[tr] - ind.trans_time[tr - 1];
					L = gamma_probability(dt, mean, shape);
				}

				for (auto loop = 0; loop < ind_effect_trans_prop; loop++) {
					for (auto j = 0; j < nie; j++) {
						const auto &prec_ie = prec_ind_eff[j];
						auto ie = ind_eff[j];

						auto sd = prec_ie.sd;
						if(anneal.phi_IE != 1) sd /= sqrt(anneal.phi_IE);
			
						auto ind_eff_store = ind.ind_effect[ie];
						auto ind_eff_prop = normal_sample(prec_ie.mean,sd);

						auto mean_fac_store = mean_fac;
						mean_fac += ind_eff_prop - ind_eff_store;

						auto al = 1.0;
						if (ind.infected == true) {
							L_prop = gamma_probability(dt, mean_av * exp(mean_fac), shape);
							al = exp(anneal.phi_L*(L_prop - L));
						}

						if(indiv.group == UNSET) check_one(al,20);
	
						ind_effect_jump[ie].ntr++;
						if (MH_proposal(al,6)) {
							ind_effect_jump[ie].nac++;

							ind.ind_effect[ie] = ind_eff_prop;

							auto cv = ind_effect[ie].covar_ref;
							L_ind_effect[cv] += Li_change_ind_eff(ind_eff_store, ind_eff_prop, prec_ie);
							if (ind.infected == true) {
								L_trans_events += L_prop - L;
								L = L_prop;
							}
						} else
							mean_fac = mean_fac_store;
					}
				}

				ind.trans_mean[tr] = mean_av * exp(mean_fac);
			}
		}
	}

	// This changes individual effects for multiple transitions
	for (auto ie = 0; ie < nind_effect; ie++) {
		const auto &ind_eff = ind_effect[ie];
		const auto &tr_mean = ind_eff.trans_mean;
		auto ntr_mean = tr_mean.size();
		if (ntr_mean > 0 && ind_eff.single_ind_effect_proposal == true) {
			PrecalcIndEff prec_ind_eff;

			vector <double> mean_av(ntrans), shape(ntrans), mean_fac(ntrans), mean_fac_store(ntrans), dt(ntrans);

			for (auto tr : tr_mean) {
				mean_av[tr] = param_value[trans[tr].mean_param];
				shape[tr] = 1; if(trans[tr].type == GAMMA) shape[tr] = param_value[trans[tr].shape_param];
			}

			for (auto i = 0; i < N; i++) {
				auto &ind = ind_value[i];
				auto &indiv = individual[i];

				set_precalc_ind_eff(i, prec_ind_eff, ie, ind_value, inv_cov_matrix);

				auto L = 0.0;

				for (auto tr : tr_mean) {
					if(indiv.inside_group == false) mean_fac[tr] = 0;
					else{
						auto mean = ind.trans_mean[tr];
						mean_fac[tr] = log(mean / mean_av[tr]);
						if (ind.infected == true) {
							dt[tr] = ind.trans_time[tr] - ind.trans_time[tr - 1];
							L += gamma_probability(dt[tr], mean, shape[tr]);
						}
					}
				}

				for (auto loop = 0; loop < ind_effect_trans_prop2; loop++) {
					auto sd = prec_ind_eff.sd;
					if(anneal.phi_IE != 1) sd /= sqrt(anneal.phi_IE);
			
					auto ind_eff_store = ind.ind_effect[ie];
					auto ind_eff_prop = normal_sample(prec_ind_eff.mean,sd);
					auto change_ind_eff = ind_eff_prop - ind_eff_store;

					if(indiv.inside_group == true){
						for (auto tr : tr_mean) {
							mean_fac_store[tr] = mean_fac[tr];
							mean_fac[tr] += change_ind_eff;
						}
					}
					
					if (false) { // Checks for correct sampling from individual effect distribution
						auto probfi = normal_probability(ind_eff_store, prec_ind_eff.mean, prec_ind_eff.sd);
						auto probif = normal_probability(ind_eff_prop, prec_ind_eff.mean, prec_ind_eff.sd);
						cout << probfi - probif + Li_change_ind_eff(ind_eff_store, ind_eff_prop, prec_ind_eff) << "zero" << endl;
					}

					auto al = 1.0;

					auto L_prop = 0.0;
					if (ind.infected == true) {
						for (auto tr : tr_mean)
							L_prop += gamma_probability(dt[tr], mean_av[tr] * exp(mean_fac[tr]), shape[tr]);
						al = exp(anneal.phi_L*(L_prop - L));
					}
					
					if(indiv.group == UNSET) check_one(al,3);
	
					ind_effect_jump[ie].ntr++;
					if (MH_proposal(al,7)) {
						ind_effect_jump[ie].nac++;

						ind.ind_effect[ie] = ind_eff_prop;

						auto cv = ind_effect[ie].covar_ref;
						L_ind_effect[cv] += Li_change_ind_eff(ind_eff_store, ind_eff_prop, prec_ind_eff);
						if (ind.infected == true) {
							L_trans_events += L_prop - L;
							L = L_prop;
						}
					} else {
						if(indiv.inside_group == true){
							for (auto tr : tr_mean) mean_fac[tr] = mean_fac_store[tr];
						}
					}
				}

				if(indiv.inside_group == true){
					for (auto tr : tr_mean) ind.trans_mean[tr] = mean_av[tr] * exp(mean_fac[tr]);
				}
			}
		}
	}
	timer[TIME_TRANS_IND_EFFECT].stop();
}


/// This makes proposals to infectivity individual effects on a given compartment
void Model::propose_infectivity_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_inf_events, vector <Jump> &ind_effect_jump, const Anneal &anneal) const {
	// cout << "Model::propose_infectivity_ind_effects()" << endl; // DEBUG
	
	if(infectivity.ind_effect.size() == 0) return;

	timer[TIME_INF_IND_EFFECT_INIT].start();
	auto inv_cov_matrix = calculate_inv_cov_matrix(param_value);
	vector <vector <double> > I_profile;
	vector <double> group_const;
	auto precalc = set_precalc_likelihood_inf(group_const,I_profile, ind_value, param_value);
	timer[TIME_INF_IND_EFFECT_INIT].stop();

	if(false){
		set_individual_quantities(ind_value,param_value);
		
		for (auto &ind : ind_value) {
			for (auto tr = 0; tr < ntrans; tr++){
				ind.trans_infectivity_change[tr] = ind.inf_single*(comp[trans[tr].to].infectivity - comp[trans[tr].from].infectivity);
			}
		}
	}
	
	if(false){
		auto L_inf_events_check = likelihood_from_precalc(group_const,precalc,I_profile,ind_value);
		
		auto dsum = 0.0;
		for (auto g = 0; g < ngroup; g++) {
			auto d = L_inf_events[g] -  L_inf_events_check[g];
			if(d < 0) d = -d;
			dsum += d;
			cout << g << " " << L_inf_events[g] << " "<< L_inf_events_check[g] << " " << d << " co\n";
		}
		cout << dsum << " dsum\n";
		if(dsum > 0.0001) emsg("done");//zz
	}

	timer[TIME_INF_IND_EFFECT].start();


	// This changes infectivity individual effects
	const auto &ind_eff = infectivity.ind_effect;
	auto nie = ind_eff.size();
	vector <PrecalcIndEff> prec_ind_eff(nie);

	for (auto i = 0; i < N; i++) {
		auto &prec = precalc[i];
		auto &ind = ind_value[i];
		auto &indiv = individual[i];
		auto g = indiv.group;

		for (auto j = 0; j < nie; j++)
			set_precalc_ind_eff(i, prec_ind_eff[j], ind_eff[j], ind_value, inv_cov_matrix);

		double infectivity_change, dL_inf_ev;

		auto inf_fac = log(ind.inf_single);
		for (auto loop = 0; loop < ind_effect_inf_prop; loop++) {
			for (auto j = 0; j < nie; j++) {
				const auto &prec_ie = prec_ind_eff[j];
				auto ie = ind_eff[j];

				auto sd = prec_ie.sd;
				if(anneal.phi_IE != 1) sd /= sqrt(anneal.phi_IE); 
			
				auto ind_eff_store = ind.ind_effect[ie];
				auto ind_eff_prop = normal_sample(prec_ie.mean,sd);

				auto inf_fac_store = inf_fac;
				inf_fac += ind_eff_prop - ind_eff_store;

				if (g != UNSET) {
					//infectivity_change = inf_av * (exp(inf_fac) - exp(inf_fac_store));
					infectivity_change = exp(inf_fac) - exp(inf_fac_store);
					dL_inf_ev = Li_change_inf(infectivity_change, prec, I_profile[g]);
				} 
				else dL_inf_ev = 0;

				auto inside = anneal.phi_L*dL_inf_ev; if(inside > 100) inside = 100;
				
				auto al = exp(inside);


				if(indiv.group == UNSET) check_one(al,1);

				ind_effect_jump[ie].ntr++;
				if (MH_proposal(al,8)) {
					ind_effect_jump[ie].nac++;

					ind.ind_effect[ie] = ind_eff_prop;

					auto cv = ind_effect[ie].covar_ref;
					L_ind_effect[cv] += Li_change_ind_eff(ind_eff_store, ind_eff_prop, prec_ie);
					if (g != UNSET) {
						L_inf_events[g] += dL_inf_ev;
						update_I_profile(infectivity_change, prec, I_profile[g]);
					}
				} else
					inf_fac = inf_fac_store;
			}
		}
		ind.inf_single = exp(inf_fac);
	}

	/// Sets the change in infectivity down a transition
	for (auto &ind : ind_value) {
		for (auto tr = 0; tr < ntrans; tr++){
			ind.trans_infectivity_change[tr] = ind.inf_single*(comp[trans[tr].to].infectivity - comp[trans[tr].from].infectivity);
		}
	}
	
	if(false){ // Checks that I_profile has been correctly updated
		vector <vector <double> > I_profile_ch;
		vector <double> group_const_ch;
		auto precalc_ch = set_precalc_likelihood_inf(group_const_ch,I_profile_ch, ind_value, param_value);
		for(auto g = 0; g < ngroup; g++){
			if(different(group_const[g],group_const_ch[g]) == true) emsg("group_const error");
			for(auto i = 0; i < I_profile[g].size(); i++){
				auto d = I_profile[g][i] - I_profile_ch[g][i];
				if(different(I_profile[g][i],I_profile_ch[g][i]) == true){
					emsg("I_profile error");
				}
			}
		}
		cout << "hh\n";
		//emsg("check");
	}
	
	timer[TIME_INF_IND_EFFECT].stop();
}
 

/// This pre-calculates quantities so likelihood can be calculated from infectivity
vector <PrecalcLiInf> Model::set_precalc_likelihood_inf(vector <double> &group_const, vector < vector <double> > &I_profile, const vector <IndValue> &ind_value, const vector <double> &param_value) const {
	// cout << "Model::set_precalc_likelihood_inf()" << endl; // DEBUG
	vector <PrecalcLiInf> prec;
	prec.resize(N);
	for (auto i = 0; i < N; i++) {
		prec[i].beta_fac = 0;
	}

	vector <int> state(N);
	for (auto &sta : state)
		sta = 0;

	group_const.resize(ngroup);
	for(auto &val : group_const) val = 0;

	I_profile.resize(ngroup);
	for (const auto &gr : group) {
		auto &Iprof = I_profile[gr.index];

		auto I = 0.0;                                            // The total infectivity
		auto S = 0.0;                                            // The total susceptibility

		auto ind_list = gr.ind_ref;                              // Makes a list of group individuals

		vector <Event> event;
		get_event_sequence(event, S, gr, ind_value);             // Gererates a sorted list of events for the entire group

		auto tmin = gr.inference_range.tmin;
		auto tmax = gr.inference_range.tmax;

		auto beta = param_value[beta_param];
		if (inf_model == FREQ_DEP)
			beta /= gr.nind;               // If frequency dependent divided by group size
		if (group_effect.on == true)                             // Adds in the group effect
			beta *= exp(param_value[group_effect.param[gr.index]]);

		auto t = tmin;                                           // Sets a time variable
		for (const auto &ev : event) {
			const auto &ind = ind_value[ev.ind];
			auto tev = ev.time;
			auto tr = ev.trans;

			auto val = -beta * S * (tev - t);
			
			if(val != 0){
				for (auto i : ind_list)
					prec[i].beta_fac += val*comp[state[i]].infectivity;
			}
			t = tev;

			auto i = ev.ind;
			if (tr == 0) {
				if (tev > tmin) {
					group_const[gr.index] += log(beta * ind.susceptibility);
					
					auto num = Iprof.size();
					for (auto i : ind_list) {
						auto c = state[i];
						if (comp[c].infectivity != 0){
							LogTerm lt; lt.list = num; lt.relinf = comp[c].infectivity;
							prec[i].inf_log_terms.push_back(lt);
						}
					}
					Iprof.push_back(I + EXT_FOI);
				}

				S -= ind.susceptibility;
			}

			I += ind.trans_infectivity_change[tr];
			state[i] = trans[tr].to;
		}

		auto val = -beta * S * (tmax - t);
		for (auto i : ind_list)
			prec[i].beta_fac += val*comp[state[i]].infectivity;
	}

	return prec;
}


/// Calculates L_inf_events from precalc (used for testing)
vector <double> Model::likelihood_from_precalc(const vector <double> &group_const, const vector <PrecalcLiInf> &precalc, const vector < vector <double> > &I_profile, const vector <IndValue> &ind_value) const
{
	vector <double> L_inf_events(ngroup);
	
	for (const auto &gr : group) {
		auto g = gr.index;
		
		auto L = group_const[g];
		auto &ind_list = gr.ind_ref;
		
		for(auto fac : I_profile[g]) L += log(fac);
		
		for (auto i : ind_list){	
			const auto &ind = ind_value[i];
			L += precalc[i].beta_fac*ind.inf_single;
		}
	
		L_inf_events[g] = L;
	}

	return L_inf_events;
}


/// Calculates the change in event likelihood as a result of infectivity changing
double Model::Li_change_inf(const double infectivity_change, const PrecalcLiInf &prec, const vector <double> &I_profile) const {
	// cout << "Model::Li_change_inf()" << endl; // DEBUG
	auto dL = prec.beta_fac*infectivity_change;
	for (auto lt : prec.inf_log_terms)
		dL += log(1 + (lt.relinf*infectivity_change/I_profile[lt.list]));
	
	//dL += log(I_profile[e] + infectivity_change) - log(I_profile[e]);

	return dL;
}


/// Updates I_profile (this stores how infectivity changes for all the infection events in the system)
void Model::update_I_profile(const double infectivity_change, const PrecalcLiInf &prec, vector <double> &I_profile) const {
	// cout << "Model::update_I_profile()" << endl; // DEBUG
	for (auto lt : prec.inf_log_terms)
		I_profile[lt.list] += lt.relinf*infectivity_change;
}


/// This works how what type of mean/event time proposals can be made
void Model::propose_joint_ie_var_initialise()
{
	for(auto cv = 0; cv < covariance.size(); cv++){
		for(auto i = 0; i <  covariance[cv].var_param.size(); i++){
			IeVarJointProposal ivjp;
			ivjp.name = "Covariance matrix: "+matrix[covariance[cv].matrix].name+"  Param:"+param[covariance[cv].var_param[i]].name;
			ivjp.covar = cv;
			ivjp.pos = i;
			
			ie_var_joint_prop.push_back(ivjp);
		}
	}
}


/// CHecks how ie/var joint proposals can be made
void Model::propose_joint_ie_var(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_inf_events, double &L_trans_events, double &prior, vector <Jump> &pjie_jump, const bool burnin, const Anneal &anneal) const
{
	timer[TIME_JOINT_IEVAR].start();
		
	for(auto i = 0; i < ie_var_joint_prop.size(); i++){
		const auto &ievj = ie_var_joint_prop[i];
		auto &jump = pjie_jump[i];
		
		const auto &cov = covariance[ievj.covar];
		const auto var_th = cov.var_param[ievj.pos];
		auto ie = cov.ind_effect_ref[ievj.pos];
		
		auto fac = exp(normal_sample(0,jump.size));
		
		auto param_st = param_value[var_th];
		
		param_value[var_th] *= fac*fac;
		
		if(inbounds(param_value) == false){
			param_value[var_th] = param_st;
		}			
		else{
			for(auto i = 0; i < N; i++) ind_value[i].ind_effect[ie] *= fac;
		
			set_individual_quantities(ind_value, param_value);
			
			auto L_inf_events_prop = calculate_L_inf_events(ind_value, param_value);

			auto L_trans_events_prop = calculate_L_trans_events(ind_value, param_value); 

			auto dL = 0.0;
			for(auto i = 0; i < L_inf_events.size(); i++) dL += L_inf_events_prop[i] - L_inf_events[i];
			
			dL += L_trans_events_prop - L_trans_events;

			auto prior_prop = calculate_prior(param_value);
			
			auto al = exp(anneal.phi_L*dL + anneal.phi_Pr*(prior_prop-prior) + 2*log(fac));
	
			jump.ntr++;
			if (MH_proposal(al,120)) {
				jump.nac++;

				L_inf_events = L_inf_events_prop;
				L_trans_events = L_trans_events_prop;
				prior = prior_prop;
				
				L_ind_effect = calculate_L_ind_effect(ind_value, param_value);

				if (burnin == true) jump.size *= prop_up;
			} 
			else{
				param_value[var_th] = param_st;
				for(auto i = 0; i < N; i++) ind_value[i].ind_effect[ie] /= fac;
				set_individual_quantities(ind_value, param_value);
		
				if (burnin == true) jump.size *= prop_down;
			}
		}
	}
	timer[TIME_JOINT_IEVAR].stop();
}

