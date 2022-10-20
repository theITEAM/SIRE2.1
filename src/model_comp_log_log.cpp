// This uses the complimentary log-log likelihood function for observations

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "model.hpp"
#include "mcmc.hpp"
#include "utils.hpp"
#include "timers.hpp"

void MCMC::update_cloglog()
{
	auto check_all = false;
	// cout << "MCMC::update()" << endl; // DEBUG

	timer[TIME_UPDATE].start();
	
	model.cloglog_propose_transmission_rate(ind_value, param_value, L_cloglog, prior, param_jump, burnin, anneal);
	
	if(check_all == true) check_chain(1);
	
	model.cloglog_propose_trans_params(ind_value, param_value, L_cloglog, prior, param_jump, burnin, anneal);

	if(check_all == true) check_chain(2);

	for (auto c = 0; c < model.ncovariance; c++){
		model.propose_covariance_matrices(model.covariance[c], ind_value, param_value, L_ind_effect[c], prior, param_jump, burnin, anneal);
	}
	
	if(check_all == true) check_chain(3);
	
	model.cloglog_propose_susceptibility_ind_effects(ind_value, param_value, L_ind_effect, L_cloglog, ind_effect_jump, anneal);

	if(check_all == true) check_chain(3);
	
	model.cloglog_propose_infectivity_ind_effects(ind_value, param_value, L_ind_effect, L_cloglog, ind_effect_jump, anneal);

	if(check_all == true) check_chain(4);

	model.cloglog_propose_joint_ie_var(ind_value, param_value, L_ind_effect, L_cloglog, prior, var_ie_joint_jump, burnin, anneal);
	
	if(check_all == true) check_chain(5);

	check_chain(100);
}

void MCMC::set_likelihoods_cloglog()
{
	L_ind_effect = model.calculate_L_ind_effect(ind_value, param_value); // Sets likelihoods for individual effects

	L_cloglog = model.calculate_L_cloglog(ind_value,param_value);

	prior = model.calculate_prior(param_value);                          // Sets the prior

	check_chain(0);
}


/// Calculates the likelihood in the case of a cloglog link function
vector <double> Model::calculate_L_cloglog(const vector <IndValue> &ind_value, const vector <double> &param) const
{
	vector <double>  L;

	auto M = cloglog.L_list.size();

	L.resize(M);
	for(auto m = 0; m < M; m++){
		L[m] = calculate_L_cloglog_ind(ind_value,param,m);
	}
	
	return L;	
}

double Model::calculate_L_cloglog_ind(const vector <IndValue> &ind_value, const vector <double> &param_value, unsigned int m) const
{
	auto beta = param_value[beta_param];
	const auto &ll = cloglog.L_list[m];
	
	auto i = ll.i;
	auto gr = group[individual[i].group];
	
	if(inf_model == FREQ_DEP) beta /= gr.nind;        // If frequency dependent divided by group size
	if (group_effect.on == true) beta *= exp(param_value[group_effect.param[gr.index]]);

	const auto &indv = ind_value[i];
	auto tr = ll.tr;
	
	double r;
	
	switch(trans[tr].type){
		case INFECTION:
			{
				r = beta*indv.susceptibility;
		
				if(cloglog.geometric_approx == false){
					auto sum = 0.0;
					for(const auto &infi : ll.inf_ind){
						sum += ind_value[infi.i].inf_single*infi.infectivity;
					}
					r *= sum; 
				}
				else{
					auto II = ll.inf_ind.size();
					
					auto sum = 0.0;
					for(const auto &infi : ll.inf_ind){
						sum += log(ind_value[infi.i].inf_single*infi.infectivity);
					}
				
					sum /= II;
					r *= exp(sum)*II;
				}
			}
			break;
		
		case EXP:
			r = 1.0/ind_value[i].trans_mean[tr];
			break;
			
		case GAMMA:
			emsg("Gamma distribution not supported");
			break;
	}
	
	
	
	if(ll.result == true){
		auto P = 1-exp(-r*cloglog.DeltaT);
		if(P < 0.000000001) P = 0.000000001;
		return log(P);
	}
	else return -r*cloglog.DeltaT;
}


/// Initialises the likelihood
void Model::initialise_cloglog()
{
	for(auto g = 0; g < ngroup; g++){
		const auto &gr = group[g];
		
		auto n = gr.ind_ref.size();
		
		const auto &tra = gr.observation_range;
				
		auto tmin = tra.tmin, tmax = tra.tmax;
		
		if(tmax > 10000){
			tmax = tmin;
			for(auto j = 0; j < n; j++){
				auto i = gr.ind_ref[j];
				const auto &in = individual[i];
				if(in.status == INFECTED){
					for(auto trai = 0; trai < ntrans; trai++){
						const auto &tran = individual[i].trans_time_range[trai];
						auto tma = tran.tmax;
						if(tma != UNSET && tma > tmax) tmax = tma;
					}
				}					
			}
		}

		auto T = int(2 + (tmax-tmin)/cloglog.DeltaT);
		
		vector < vector < vector <int> > > state;
		state.resize(T);
		
		for(auto ti = 0; ti < T; ti++){
			auto t = tmin + ti*cloglog.DeltaT;
				
			state[ti].resize(ncomp);
			
			for(auto j = 0; j < n; j++){
				auto i = gr.ind_ref[j];
				const auto &in = individual[i];
				
				auto c = in.initial_comp;
					
				switch(in.status){
					case INFECTED:
						{
							for(auto trai = c; trai < ntrans; trai++){
								if(trans[trai].from != c) emsg("Problem");
												
								const auto &tran = individual[i].trans_time_range[trai];
								if(tran.tmin > t) break;
												
								c = trans[trai].to;
							}

							state[ti][c].push_back(i);
						}
						break;
						
					case NOT_INFECTED: 
						if(c != 0) emsg("Prob");
						state[ti][0].push_back(i);
						break;
										
					case UNKNOWN: emsg("Problem"); break;
				}
			}
		}
			
		for(auto ti = 0; ti < T-1; ti++){
			auto t = tmin + ti*cloglog.DeltaT;
			
			vector <InfInd> inf_list;       // Gets a list of infected individuals
			for(auto c = 1; c < ncomp; c++){
				if(comp[c].infectivity != 0){
					for(auto i : state[ti][c]){
						InfInd infind; infind.i = i; infind.infectivity = comp[c].infectivity;
						inf_list.push_back(infind);
					}
				}
			}		

			for(auto tr = 0; tr < ntrans; tr++){
				auto c = trans[tr].from;
				
				for(auto i : state[ti][c]){
					L_CLogLog L;
					L.t = t; L.i = i; L.tr = tr; L.result = false;
					if(find_in(state[ti+1][c],i) == UNSET) L.result = true;
					if(trans[tr].type == INFECTION) L.inf_ind = inf_list;
					
					cloglog.L_list.push_back(L);
				}
			}
		}
		
		cloglog.L_ref.resize(N);
		for(auto i = 0; i < N; i++){
			cloglog.L_ref[i].resize(ntrans);		
		}
	}
	
	for(auto k = 0; k < cloglog.L_list.size(); k++){
		auto &ll = cloglog.L_list[k];
	
		cloglog.L_ref[ll.i][ll.tr].direct.push_back(k);
		
		for(auto infi : ll.inf_ind){
			if(ll.tr != 0) emsg("prob");
			cloglog.L_ref[infi.i][ll.tr].indirect.push_back(k);
		}
	}
}


/// Proposes changes to beta
void Model::cloglog_propose_transmission_rate(const vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_cloglog, double &prior, vector <Jump> &param_jump, const bool burnin, const Anneal &anneal) const 
{
	auto th = beta_param;
	auto &jump = param_jump[th];
	
	auto beta = param_value[th];
	
	// Makes a change to beta
	auto beta_prop = beta + normal_sample(0, jump.size);
	param_value[th] = beta_prop;

	// Calculates the Metropolis-Hastings acceptance probability
	auto prior_change = calculate_prior_change(th, beta, param_value);
	
	if(inbounds(param_value) == false){
		param_value[th] = beta;
		return;
	}
		
	auto L_cloglog_prop = calculate_L_cloglog(ind_value,param_value);
		
	auto dL = 0.0;
	for(auto m = 0; m < L_cloglog.size(); m++){
		dL += L_cloglog_prop[m] - L_cloglog[m]; 	
	}
	if(dL > 100) dL = 100;
	
	auto al = exp(anneal.phi_L*dL + anneal.phi_Pr*prior_change);
		
	jump.ntr++;
	if (MH_proposal(al,16)) {
		jump.nac++;
		L_cloglog = L_cloglog_prop;
			
		prior += prior_change;
		if (burnin == true) jump.size *= prop_up;
	}
	else {
		param_value[th] = beta;
		if (burnin == true) jump.size *= prop_down;
	}
}


/// This makes proposals to transition parameters
void Model::cloglog_propose_trans_params(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_cloglog, double &prior, vector <Jump> &param_jump, const bool burnin, const Anneal &anneal) const 
{	
	timer[TIME_TRANS_PARAM].start();

	for (auto th = 0; th < nparam; th++) {
		if (param[th].type == TRANS_MEAN) {
			auto &jump = param_jump[th];
			auto param_store = param_value[th];

			// Makes a change to one of the transition parameters
			param_value[th] += normal_sample(0, jump.size);

			// Calculates the Metropolis-Hastings acceptance probability
			auto prior_change = calculate_prior_change(th, param_store, param_value);
			
			if(inbounds(param_value) == false){
				param_value[th] = param_store;
			}
			else{
				set_individual_quantities(ind_value, param_value); 
				
				auto L_cloglog_prop = calculate_L_cloglog(ind_value,param_value);
		
				auto dL = 0.0;
				for(auto m = 0; m < L_cloglog.size(); m++){
					dL += L_cloglog_prop[m] - L_cloglog[m]; 	
				}
				if(dL > 100) dL = 100;

				auto al = exp(anneal.phi_L*dL + anneal.phi_Pr*prior_change);
		
				jump.ntr++;
				if (MH_proposal(al,17)) {
					jump.nac++;
					L_cloglog = L_cloglog_prop;
					prior += prior_change;
					if (burnin == true) jump.size *= prop_up;
				} else {
					param_value[th] = param_store;
					if (burnin == true) jump.size *= prop_down;
				}
			}
		}
	}

	set_individual_quantities(ind_value, param_value);          // Sets individual transition parameters and infectivity

	timer[TIME_TRANS_PARAM].stop();
}


/// This makes proposals to susceptibility individual effects
void Model::cloglog_propose_susceptibility_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_cloglog, vector <Jump> &ind_effect_jump, const Anneal &anneal) const {
	// cout << "Model::propose_susceptibility_ind_effects()" << endl; // DEBUG
	const auto &ind_eff = trans[0].ind_effect;
	auto nie = ind_eff.size();
	if (nie == 0) return;

	timer[TIME_SUS_IND_EFFECT_INIT].start();
	auto inv_cov_matrix = calculate_inv_cov_matrix(param_value);
	timer[TIME_SUS_IND_EFFECT_INIT].stop();

	timer[TIME_SUS_IND_EFFECT].start();
	vector <PrecalcIndEff> prec_ind_eff(nie);

	for (auto i = 0; i < N; i++) {
		auto &ind = ind_value[i];
		auto &indiv = individual[i];
		auto g = indiv.group;
	
		for (auto j = 0; j < nie; j++){
			set_precalc_ind_eff(i, prec_ind_eff[j], ind_eff[j], ind_value, inv_cov_matrix);
		}
		
		for (auto j = 0; j < nie; j++) {
			const auto &prec_ie = prec_ind_eff[j];
			auto ie = ind_eff[j];

			auto sd = prec_ie.sd, mean = prec_ie.mean;
			if(anneal.phi_IE != 1) sd /= sqrt(anneal.phi_IE);
			
			auto ind_eff_store = ind.ind_effect[ie];
			auto ind_eff_prop = normal_sample(mean, sd);

			auto sus_store = ind.susceptibility;
			
			ind.susceptibility *= exp(ind_eff_prop - ind_eff_store);
			
			const auto &direct = cloglog.L_ref[i][0].direct;
			
			vector <double> store(direct.size());
			
			auto dL = 0.0;
			
			for(auto j = 0; j < direct.size(); j++){
				auto m = direct[j];

				/*
				if(individual[i].status == NOT_INFECTED){
					cout << individual[cloglog.L_list[m].i].id << " " << cloglog.L_list[m].t << " " << m << " " << cloglog.L_list[m].tr << " " <<   cloglog.L_list[m].inf_ind.size() << "\n";
					for(auto infi : cloglog.L_list[m].inf_ind) cout << individual[infi.i].id << "," << infi.infectivity << "  ";
					cout << "inf\n";
				}				
				*/
				
				store[j] = calculate_L_cloglog_ind(ind_value,param_value,m);
				dL += store[j] - L_cloglog[m];
			}
			
			auto inside = anneal.phi_L*dL; if(inside > 100) inside = 100;
			
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
				
				for(auto j = 0; j < direct.size(); j++){
					auto m = direct[j];
					L_cloglog[m] = store[j];
				}
			} 
			else{
				ind.susceptibility = sus_store;
			}
		}
	}
	timer[TIME_SUS_IND_EFFECT].stop();
}


/// This makes proposals to susceptibility individual effects
void Model::cloglog_propose_infectivity_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_cloglog, vector <Jump> &ind_effect_jump, const Anneal &anneal) const {
	// cout << "Model::propose_susceptibility_ind_effects()" << endl; // DEBUG
	const auto &ind_eff = infectivity.ind_effect;
	auto nie = ind_eff.size();
	if (nie == 0) return;

	timer[TIME_INF_IND_EFFECT_INIT].start();
	auto inv_cov_matrix = calculate_inv_cov_matrix(param_value);
	timer[TIME_INF_IND_EFFECT_INIT].stop();

	timer[TIME_INF_IND_EFFECT].start();
	vector <PrecalcIndEff> prec_ind_eff(nie);

	for (auto i = 0; i < N; i++) {
		auto &ind = ind_value[i];
		auto &indiv = individual[i];
		auto g = indiv.group;
	
		for (auto j = 0; j < nie; j++){
			set_precalc_ind_eff(i, prec_ind_eff[j], ind_eff[j], ind_value, inv_cov_matrix);
		}
		
		for (auto j = 0; j < nie; j++) {
			const auto &prec_ie = prec_ind_eff[j];
			auto ie = ind_eff[j];

			auto sd = prec_ie.sd, mean = prec_ie.mean;
			if(anneal.phi_IE != 1) sd /= sqrt(anneal.phi_IE);
			
			auto ind_eff_store = ind.ind_effect[ie];
			auto ind_eff_prop = normal_sample(mean, sd);

			auto inf_store = ind.inf_single;
			
			ind.inf_single *= exp(ind_eff_prop - ind_eff_store);
			
			const auto &indirect = cloglog.L_ref[i][0].indirect;
			
			vector <double> store(indirect.size());
			
			auto dL = 0.0;
			
			/*
			cout << indiv.id << ":|n";
			for(auto j = 0; j < indirect.size(); j++){
					auto m = indirect[j];
				cout << individual[cloglog.L_list[m].i].id << " " << cloglog.L_list[m].t << endl;
			}
			cout << "\n";
			*/
			
			for(auto j = 0; j < indirect.size(); j++){
				auto m = indirect[j];

				/*
				if(individual[i].status == NOT_INFECTED){
					cout << individual[cloglog.L_list[m].i].id << " " << cloglog.L_list[m].t << " " << m << " " << cloglog.L_list[m].tr << " " <<   cloglog.L_list[m].inf_ind.size() << "\n";
					for(auto infi : cloglog.L_list[m].inf_ind) cout << individual[infi.i].id << "," << infi.infectivity << "  ";
					cout << "inf\n";
				}				
				*/
				
				store[j] = calculate_L_cloglog_ind(ind_value,param_value,m);
				dL += store[j] - L_cloglog[m];
			}
			
			auto inside = anneal.phi_L*dL; if(inside > 100) inside = 100;
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
				
				for(auto j = 0; j < indirect.size(); j++){
					auto m = indirect[j];
					L_cloglog[m] = store[j];
				}
			} 
			else{
				ind.inf_single = inf_store;
			}
		}
	}
	
	/// Sets the change in infectivity down a transition
	for (auto &ind : ind_value) {
		for (auto tr = 0; tr < ntrans; tr++){
			ind.trans_infectivity_change[tr] = ind.inf_single*(comp[trans[tr].to].infectivity - comp[trans[tr].from].infectivity);
		}
	}
	
	timer[TIME_INF_IND_EFFECT].stop();
}


/// CHecks how ie/var joint proposals can be made
void Model::cloglog_propose_joint_ie_var(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_cloglog, double &prior, vector <Jump> &pjie_jump, const bool burnin, const Anneal &anneal) const
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
			
			auto L_cloglog_prop = calculate_L_cloglog(ind_value,param_value);
		
			auto dL = 0.0;
			for(auto m = 0; m < L_cloglog.size(); m++){
				dL += L_cloglog_prop[m] - L_cloglog[m]; 	
			}
			if(dL > 100) dL = 100;
		
			auto prior_prop = calculate_prior(param_value);
			
			auto al = exp(anneal.phi_L*dL + anneal.phi_Pr*(prior_prop-prior) + 2*log(fac));
	
			jump.ntr++;
			if (MH_proposal(al,120)) {
				jump.nac++;

				L_cloglog = L_cloglog_prop;
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



/// Outputs the trace plot
void MCMC::cloglog_trace_output(const int s) 
{	
	trace << s;
	for (auto val : param_value) trace << '\t' << val;
	
	for (auto der : model.derived) trace << '\t' << model.calculate_derived(der,param_value);
	
	auto PP = 0.0;
	for (auto c = 0; c < model.ncovariance; c++){
		trace << '\t' << L_ind_effect[c];
		PP += L_ind_effect[c];
	}
	
	auto sum = 0.0;
	
	trace << '\t' << sum;
	trace << '\t' << 0;
	
	sum = 0.0;
	
	trace << '\t' << sum;
	trace << '\t' << prior;
	PP += prior;
	
	trace << '\t' << PP;
	auto num_infected = 0;
	
	trace << '\t' << num_infected;
	trace << '\t' << log(anneal.phi_L); 
	trace << endl;
}


/// Checks working correctly
void MCMC::cloglog_check(unsigned int num) 
{
	timer[TIME_CHECK].start();

	// Checks that individual quantities are correctly set
	if (model.nind_effect > 0) {
		auto ind_value_check = ind_value;
		model.set_individual_quantities(ind_value, param_value);
		for (auto i = 0; i < model.N; i++) {
			if(model.individual[i].group != UNSET){
				auto &ind = ind_value[i];
				const auto &ind_check = ind_value_check[i];

				if (different(ind.susceptibility, ind_check.susceptibility) == true)
					emsg("Susceptibility problem");
				ind.susceptibility = ind_check.susceptibility;

				if (different(ind.inf_single, ind_check.inf_single) == true)
					emsg("Susceptibility problem");
				ind.inf_single = ind_check.inf_single;
			
				for (auto tr = 0; tr < model.ntrans; tr++) {
					if (different(ind.trans_mean[tr], ind_check.trans_mean[tr]) == true){
						cout << model.individual[i].id << " " << ind.trans_mean[tr] << " " << ind_check.trans_mean[tr] << "\n";
						emsg("Trans_mean problem");
					}
					
					ind.trans_mean[tr] = ind_check.trans_mean[tr];

					if (different(ind.trans_infectivity_change[tr], ind_check.trans_infectivity_change[tr]) == true) {
						cout << tr << " " << i << " " << ind.trans_infectivity_change[tr] << " " << ind_check.trans_infectivity_change[tr] << " ch\n";
						emsg("trans_infectivity_change problem");
					}
					ind.trans_infectivity_change[tr] = ind_check.trans_infectivity_change[tr];
				}
			}
		}
	}

	auto L_cloglog_check = model.calculate_L_cloglog(ind_value,param_value);
	for(auto m = 0; m < L_cloglog.size(); m++){	
		if (different(L_cloglog[m], L_cloglog_check[m]) == true) {
			emsg("L_cloglog problem");
		}
	}
	L_cloglog = L_cloglog_check;

	// Checks the prior is correctly set
	auto prior_check = model.calculate_prior(param_value);
	if (different(prior, prior_check) == true) emsg("prior problem");
	prior = prior_check;

	// Checks the individual effect likelihoods are correctly set
	auto L_ind_effect_check = model.calculate_L_ind_effect(ind_value, param_value);
	
	for (auto c = 0; c < model.ncovariance; c++) {
		if (different(L_ind_effect[c], L_ind_effect_check[c]) == true)
			emsg("L_ind_effect problem");
	}
	L_ind_effect = L_ind_effect_check;

	// Check for not a number
	if (std::isnan(prior) || std::isinf(prior))
		emsg("Prior is not a number");

	for (auto c = 0; c < model.ncovariance; c++) {
		if (std::isnan(L_ind_effect[c]) || std::isinf(L_ind_effect[c]))
			emsg("L_ind_effect is not a number");
	}

	timer[TIME_CHECK].stop();
}

void Model::check_relationship(string ind1, string ind2, const vector < vector <double> > A, double expected)
{
	auto i1 = 0; while(i1 < N && individual[i1].id != ind1) i1++;
	auto i2 = 0; while(i2 < N && individual[i2].id != ind2) i2++;
	
	cout << A[i1][i2] << " " << A[i2][i1] << " should be " << expected << "\n";
}


/// Atificially sets time ranges to match with 
void Model::set_time_ranges()
{
	auto DeltaT = 2.0;
	
	for(auto &ind : individual){
		if(ind.status == INFECTED){
			for(auto tr = 0; tr < ntrans; tr++){
				auto &tra = ind.trans_time_range[tr];
				if(tra.tmin != 0){
					tra.tmin = int(tra.tmin/DeltaT)*DeltaT;
					tra.tmax = tra.tmin + DeltaT;
				}
			}
		}
	}
}
