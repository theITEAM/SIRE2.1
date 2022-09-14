/// This provides functions for the MCMC chain

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "mcmc.hpp"
#include "model.hpp"
#include "check.hpp"
#include "utils.hpp"
#include "timers.hpp"


/// This initialises the mcmc chain
MCMC::MCMC(const Model &model) : model(model) 
{
	if(mpi.core == 0){
		cout << "Initialising MCMC chain..." << endl;
		if(mpi.ncore > 1 && model.algorithm == ALG_MCMC) emsg("Only one core is needed to run MCMC");
	}
	mpi.barrier();
				
	string file = "trace.txt";
	if(MCMC_PARA == true){
		auto seed = 10+mpi.core*1000;
		set_seed(seed);                             // Sets the random seed
		srand(seed);                      

		stringstream ss; ss << "trace_" << mpi.core << ".txt";
		file = ss.str();
	}
	
	trace_initialise(file);                                           // Initialises the trace plot

	model.sample_initial_param(param_value,ind_value);
	
	ind_effect_posterior_mean.resize(model.N);                        // Initialises individual effect means
	for (auto &ie : ind_effect_posterior_mean) {
		ie.ind_effect_sum.resize(model.nind_effect);
		for (auto &val : ie.ind_effect_sum) val = 0;
		ie.ind_effect_sum2.resize(model.nind_effect);
		for (auto &val : ie.ind_effect_sum2) val = 0;
	}
	nind_effect_posterior_mean = 0;

	model.set_initial_events(ind_value, param_value);                 // Sample initial individual events

	initialise_proposals();                                           // Initialises MCMC proposals

	phi = 0; phi_run = 1;                                             // Initial quench temperature set to zero
	
	//load_phi_schedule();
}


/// Sets all the likelihoods at the start of the chain
void MCMC::set_likelihoods() {
	if(model.cloglog.on == true){ set_likelihoods_cloglog(); return;}
	
	// cout << "MCMC::set_likelihoods()" << endl; // DEBUG
	L_ind_effect = model.calculate_L_ind_effect(ind_value, param_value); // Sets likelihoods for individual effects

	L_inf_events = model.calculate_L_inf_events(ind_value, param_value); // Sets likelihoods for infection events

	L_trans_events = model.calculate_L_trans_events(ind_value, param_value); // Sets likelihoods for infection events

	L_diag_test = model.calculate_L_diag_test(ind_value, param_value); // The likelihood for diagnostic test results

	prior = model.calculate_prior(param_value);                        // Sets the prior

	check_chain(0);
}


/// Used for checking when L_inf_events is incorrect
void MCMC::diff(unsigned int num)
{
	auto L_inf_events_check = model.calculate_L_inf_events(ind_value, param_value);
	cout << num << " " << L_inf_events[14] - L_inf_events_check[14] << " dif\n";
}


/// Performs a series of MCMC proposal which consist of an "update"
void MCMC::update() 
{
	if(model.cloglog.on == true){ update_cloglog(); return;}
	
	auto check_all = false;
	// cout << "MCMC::update()" << endl; // DEBUG

	timer[TIME_UPDATE].start();
	
	model.propose_transmission_rate(ind_value, param_value, L_inf_events, prior, param_jump, burnin, quench);
	
	if(check_all == true) check_chain(1);

	for (auto c = 0; c < model.ncovariance; c++)
		model.propose_covariance_matrices(model.covariance[c], ind_value, param_value, L_ind_effect[c], prior, param_jump, burnin, quench);

	if(check_all == true) check_chain(2);
	
	model.propose_trans_params(ind_value, param_value, L_trans_events, prior, param_jump, burnin, quench);

	if(check_all == true) check_chain(3);

	if (model.group_effect.on == true) {
		model.propose_group_effect_sigma(param_value, prior, param_jump, burnin, quench);
		model.propose_group_effect(ind_value, param_value, L_inf_events, prior, param_jump, burnin, quench);
	}

	if(check_all == true) check_chain(4);

	model.propose_fixed_effects(ind_value, param_value, L_inf_events, L_trans_events, prior, param_jump, burnin,quench);

	if(check_all == true) check_chain(5);

	model.propose_snp_effects(ind_value, param_value, L_inf_events, L_trans_events, prior, param_jump, burnin, quench);
	
	if(check_all == true) check_chain(6);

	model.propose_event_times(ind_value, param_value, L_inf_events, L_trans_events, L_diag_test, event_jump, burnin, quench);

	if(check_all == true) check_chain(7);
	
	model.propose_mean_event_times(ind_value, param_value, L_inf_events, L_trans_events, L_diag_test, prior, mean_time_jump, burnin, quench);
	
	if(check_all == true) check_chain(8);
	
	model.propose_add_rem(ind_value, param_value, L_inf_events, L_trans_events, L_diag_test, add_rem_jump, burnin, quench);
	
	if(check_all == true) check_chain(9);
	
	model.propose_susceptibility_ind_effects(ind_value, param_value, L_ind_effect, L_inf_events, ind_effect_jump, quench);

	if(check_all == true) check_chain(10);

	model.propose_infectivity_ind_effects(ind_value, param_value, L_ind_effect, L_inf_events, ind_effect_jump, quench);

	if(check_all == true) check_chain(11);

	model.propose_trans_ind_effects(ind_value, param_value, L_ind_effect, L_trans_events, ind_effect_jump, quench);

	if(check_all == true) check_chain(12);

	model.propose_joint_ie_var(ind_value, param_value, L_ind_effect, L_inf_events, L_trans_events, prior, var_ie_joint_jump, burnin, quench);

	if(check_all == true) check_chain(120);

	model.propose_Se_Sp(ind_value, param_value, L_diag_test, prior, param_jump, burnin, quench);

	if(check_all == true) check_chain(13);

	timer[TIME_UPDATE].stop();

	check_chain(0);
}


/// Sets inverse temperatures for quench
void MCMC::set_quench(unsigned int s)
{
	if(model.nquench != UNSET){
		if(s < model.nprequench) phi = 0;
		else{
			if(phi_schedule.size() > 0){
				phi = get_phi_from_schedule(s-model.nprequench);
				//if(s%10 == 0) cout << s << " " << phi << " qq\n";
			}
			else{
				double frac;
				if(s < model.nprequench + model.nquench) frac = double(s-model.nprequench)/model.nquench;
				else frac = 1;
		
				phi = pow(frac,model.quench_power);
			}
		}
		
		/*
		auto frac = double(s+1)/model.nquench;
		if(frac > 1) frac = 1;
		phi = pow(frac,model.quench_power);
		
		phi *= phi_ch[mpi.core];
		*/
		
		//phi *= 0.8;
		//cout << phi << "\n";
		/*
		if(s < 100) phi = 0;
		else{
			if(s >= model.nburnin) phi = phi_run;
			else{
				int smin = s - 500; if(smin < 0) smin = 0;
			
				auto av = 0.0, av2 = 0.0; 
				for(auto ss = smin; ss < s; ss++){
					av += burnin_Li[ss];
					av2 += burnin_Li[ss]*burnin_Li[ss];
				}
				auto var = (av2/(s-smin)) - (av/(s-smin))*(av/(s-smin));

				cout << s << " " << var << " " << phi << "var\n";
				phi +=	0.25/var;
				if(phi > phi_run) phi = phi_run;
			}
		}
		*/
	}
	else phi = 1;

	quench.phi_L = phi*phi_run;
	quench.phi_DT = phi*phi_run;
	quench.phi_IE = 1;
	quench.phi_Pr = 1;
}


/// Works out the optimum schedule for quenching
void MCMC::optimum_quench_schedule()
{
	auto ra = 200;

	vector <double> var_st, phi_st;
	for(auto s = model.nprequench; s < model.nprequench+model.nquench; s++){		
		int smin = s - ra; if(smin < 0) smin = 0;
		int smax = s + ra; if(smax > model.nburnin) smax = model.nburnin;
		
		/*
		auto av = 0.0, av2 = 0.0; 
		for(auto ss = smin; ss < smax; ss++){
			av += burnin_Li[ss];
			av2 += burnin_Li[ss]*burnin_Li[ss];
		}
		auto var = (av2/(smax-smin)) - (av/(smax-smin))*(av/(smax-smin));
	
		*/
		
		auto av = 0.0; 
		for(auto ss = smin; ss < smax; ss++) av += burnin_Li[ss];
		av /= smax-smin;
		
		auto av2 = 0.0, nav = 0.0; 
		for(auto ss = smin; ss < smax; ss++){
			if(burnin_Li[ss] > av){ av2 += (burnin_Li[ss]-av)*(burnin_Li[ss]-av); nav++;}
		}
		av2 /= nav;
		auto var = av2;
		
		var_st.push_back(var);
		phi_st.push_back(burnin_phi[s]);
	}
	
	ofstream phi_schedule_out(model.output_dir + "/phi_schedule.txt");
	auto A = 0.0;
	for(auto i = 0; i < var_st.size()-1; i++){
		A += 0.5*(var_st[i] + var_st[i+1])*(phi_st[i+1]-phi_st[i]);
	}
	
	phi_schedule_out << 0 << endl;
	
	auto ndiv = 1000;
	auto ph = 0.0;
	auto dA = A/ndiv;
	auto i = 0;
	auto Adone = 0.0;
	for(auto d = 0; d < ndiv; d++){
		do{
			auto frac = (ph-phi_st[i])/(phi_st[i+1]-phi_st[i]);
			auto var = var_st[i]*(1-frac) + var_st[i+1]*frac;
		
			auto Aleft = 0.5*(var + var_st[i+1])*(phi_st[i+1]-ph);
			if(Adone + Aleft > dA){
				auto AA = dA-Adone;
				
				auto a = (var_st[i+1]-var)/(2*(phi_st[i+1]-ph));
				auto b = var;
				auto c = -AA;
				
				auto dph = (-b+sqrt(b*b-4*a*c))/(2*a);
				
				auto fracnew = ((ph+dph)-phi_st[i])/(phi_st[i+1]-phi_st[i]);
				auto varnew = var_st[i]*(1-fracnew) + var_st[i+1]*fracnew;
				Adone = 0; ph += dph; 
				phi_schedule_out << ph << endl;
				break;
			}
			
			if(Aleft < 0) emsg("PP");
			Adone += Aleft; ph = phi_st[i+1]; i++;
		}while(i < phi_st.size());
	}
	
	phi_schedule_out << 1 << endl;
	phi_schedule_out.close();
	
	
	load_phi_schedule();
	ofstream p1("p1.txt");
	for(auto s = 0; s < var_st.size()-1; s++){
		p1 << s << " " << (phi_st[s+1]-phi_st[s])*var_st[s]	<< endl;
	}
	
	ofstream p2("p2.txt");
	for(auto s = 0; s < phi_schedule.size()-1; s++){
		auto ph = phi_schedule[s];
		
		auto ss = 0; while(ss < phi_st.size() && ph > phi_st[ss]) ss++;
		
		p2 << s << " " << (phi_schedule[s+1]-phi_schedule[s])*var_st[ss] << " " << phi_schedule[s] << " " << phi_st[ss] << endl;
	}
	
	ofstream phiout(model.output_dir + "/phi_variation.txt");
	auto P = model.quench_power;
	for(auto s = 0; s < var_st.size()-1; s++){
		auto f = double(s)/var_st.size();
		phiout << s << " " << phi_st[s] << " " << get_phi_from_schedule(s) << " " << var_st[s] << " " << phi_st[s+1]-phi_st[s]  << endl;
	}
}


/// If available loads schedule in phi
void MCMC::load_phi_schedule()
{
	ifstream phi_schedule_in(model.output_dir + "/phi_schedule.txt");
	if(!phi_schedule_in) return;
		
	//cout << "Loading phi schedule..." << endl;
	
	phi_schedule.clear();
	while(!phi_schedule_in.eof()){
		double phi;
		phi_schedule_in >> phi;
		phi_schedule.push_back(phi);
	}
}


/// Gets the value of pho from the loaded scedule
double MCMC::get_phi_from_schedule(unsigned int s) const
{
	if(s >= model.nquench) return 1;
	
	auto f = (phi_schedule.size()-1)*(double(s)/model.nquench);
	auto i = int(f);
	auto frac = f-i;
	return phi_schedule[i]*(1-frac) + phi_schedule[i+1]*frac;
}

 
/// Intialises the trace plot
void MCMC::trace_initialise(const string file) {
	// cout << "MCMC::trace_initialise()" << endl; // DEBUG
	trace.open(model.output_dir + "/" + file);
	trace << "state";
	for (auto par : model.param) trace << '\t' << par.name;
	for (auto der : model.derived) trace << '\t' << der.name;
	for (auto c = 0; c < model.ncovariance; c++)
		trace << "\tL_ind_effect " << c;
	trace << "\tL_inf_events"
		<< "\tL_trans_events"
		<< "\tL_diag_test"
		<< "\tPrior"
		<< "\tPosterior"
		<< "\tNumber infected"
		<< "\tlog(phi)"
		<< endl;
}


/// Outputs the trace plot
void MCMC::trace_output(const int s) 
{
	// cout << "MCMC::trace_output()" << endl; // DEBUG
	
	if(model.cloglog.on == true){ cloglog_trace_output(s); return;}
	
	trace << s;
	for (auto val : param_value) trace << '\t' << val;
	
	for (auto der : model.derived) trace << '\t' << model.calculate_derived(der,param_value);
	
	auto PP = 0.0;
	for (auto c = 0; c < model.ncovariance; c++){
		trace << '\t' << L_ind_effect[c];
		PP += L_ind_effect[c];
	}
	
	auto sum = 0.0;
	for (auto g = 0; g < model.ngroup; g++){
		sum += L_inf_events[g];
		PP += L_inf_events[g];
	}
	
	trace << '\t' << sum;
	trace << '\t' << L_trans_events;
	PP += L_trans_events;
	
	sum = 0.0;
	for (auto g = 0; g < model.ngroup; g++){
		sum += L_diag_test[g];
		PP += L_diag_test[g];
	}
	
	trace << '\t' << sum;
	trace << '\t' << prior;
	PP += prior;
	
	trace << '\t' << PP;
	auto num_infected = 0;
	for (const auto &ind : ind_value) {
		if (ind.infected == true)
			num_infected++;
	}
	trace << '\t' << num_infected;
	trace << '\t' << log(quench.phi_L); 
	trace << endl;
}


/// Combines trace for each of the chaings
void MCMC::output_combined_trace(const vector < vector <Sample> > &sample) const
{
	ofstream trace(model.output_dir + "/trace_combine.txt");
	trace << "state";
	for (auto par : model.param) trace << '\t' << par.name;
	trace	<< endl;
		
	auto s = 0;
	for(auto ch = 0; ch < mpi.ncore; ch++){
		for(auto i = 0; i < sample[ch].size(); i++){
			trace << s;
			for (auto val : sample[ch][i].param_value) trace << '\t' << val;
			trace	<< endl;
			s++;
		}
	}
}


/// Given a vector of values this returns statistical information
Statistics MCMC::get_statistic(const vector <double> &vec) const {
	// cout << "MCMC::get_statistic()" << endl; // DEBUG
	Statistics stat;

	auto n = vec.size();
	if (n == 0) {
		stat.mean = "---";
		stat.CImin = "---";
		stat.CImax = "---";
		stat.ESS = "---";
	} else {
		auto sum = 0.0, sum2 = 0.0;
		for (auto i = 0; i < vec.size(); i++) {
			sum += vec[i];
			sum2 += vec[i] * vec[i];
		}
		sum /= n;
		sum2 /= n;

		stat.mean = to_string(sum);

		vector <double> vec2 = vec;
		sort(vec2.begin(), vec2.end());

		if (n >= 2) {
			auto i = (unsigned int)((n - 1) * 0.025);
			auto f = (n - 1) * 0.025 - i;
			stat.CImin = to_string(vec2[i] * (1 - f) + vec2[i + 1] * f);

			i = (unsigned int)((n - 1) * 0.975);
			f = (n - 1) * 0.975 - i;
			stat.CImax = to_string(vec2[i] * (1 - f) + vec2[i + 1] * f);
		} else {
			stat.CImin = to_string(vec2[0]);
			stat.CImax = to_string(vec2[0]);
		}

		vec2 = vec;
		auto var = sum2 - sum * sum;
		if (var <= 0.0000000001 || n <= 2)
			stat.ESS = "---";
		else {
			auto sd = sqrt(var);
			for (auto i = 0; i < vec2.size(); i++)
				vec2[i] = (vec2[i] - sum) / sd;

			auto sum = 1.0;
			for (auto d = 1u; d < n / 2; d++) {         // Calculates the effective sample size
				auto a = 0.0;
				for (auto i = 0; i < n - d; i++)
					a += vec2[i] * vec2[i + d];
				auto cor = a / (n - d);
				if (cor < 0)
					break;
				sum += 2 * cor;
			}
			stat.ESS = to_string(int(n / sum));
		}
	}

	return stat;
}


/// Outputs statistics at the end of the run
void MCMC::output_statistics()
{
	// cout << "MCMC::output_statistics()" << endl; // DEBUG
	
	auto chain_sample = mpi.gather_sample(sample);
	
	auto chain_ind_effect_posterior_mean = mpi.gather_indPM(ind_effect_posterior_mean);
	
	ofstream fout;
	if(mpi.core == 0) fout.open(model.output_dir + "/diagnostics.txt");
	
	if(mpi.core == 0){
		output_combined_trace(chain_sample);
		
		// Calculates prediction accuracies

		if (model.pred_acc.size() > 0) {    // save prediction accuracies to a CSV file
			ofstream pa_csv(model.output_dir + "/" + "pred_accs.csv");
			pa_csv << "group, trait, value\n";

			fout << "Prediction accuracies:\n" << endl;
			for (const auto &pa : model.pred_acc) {
				fout << pa.name << '\t';

				for (auto ie = 0; ie < model.nind_effect; ie++) {
					if (model.individual[0].ind_effect_value[ie] != UNSET) {
						auto av_PM = 0.0, av_PM2 = 0.0;
						auto av_actual = 0.0, av_actual2 = 0.0;
						auto av_PM_actual = 0.0;
						auto nav = 0.0;

						for (auto i : pa.ind) {
							// Calculate the posterior mean and SD for the individual effect
							auto PM = 0.0;
							for(auto ch = 0; ch < mpi.ncore; ch++){
								const auto &iepm = chain_ind_effect_posterior_mean[ch][i];
								PM += iepm.ind_effect_sum[ie] / nind_effect_posterior_mean;
							}
							PM /= mpi.ncore;
							
							// The actual value for the individual effect (from the data table)
							auto actual = model.individual[i].ind_effect_value[ie];

							av_PM += PM;
							av_PM2 += PM * PM;
							av_actual += actual;
							av_actual2 += actual * actual;
							av_PM_actual += PM * actual;
							nav++;
						}
						av_PM /= nav;
						av_PM2 /= nav;
						av_actual /= nav;
						av_actual2 /= nav;
						av_PM_actual /= nav;

						auto var_PM = av_PM2 - av_PM * av_PM;
						auto var_actual = av_actual2 - av_actual * av_actual;

						// Calculates the Pearson correlation coefficient
						auto cor = (av_PM_actual - av_PM * av_actual) / sqrt(var_PM * var_actual);

						fout << model.ind_effect[ie].name << ": " << cor << '\t';
						pa_csv << pa.name << ", " << model.ind_effect[ie].name << ", " << cor << '\n';
					}
				}
				fout << '\n';
			}
			fout << '\n' << endl;
			pa_csv.close();
		}

		// Save EBVs to a CSV file
		if(model.nind_effect > 0){
			ofstream ebvs_csv(model.output_dir + "/" + "ebvs.csv");
			ebvs_csv << "id";
			for (auto ie = 0; ie < model.nind_effect; ++ie){
				if (model.individual[0].ind_effect_value[ie] != UNSET){
					ebvs_csv << ", " << model.ind_effect[ie].name;
				}
			}
			ebvs_csv << '\n';

			for (auto i = 0; i < model.individual.size(); ++ i) {
				ebvs_csv << i + 1;
				for (auto ie = 0; ie < model.nind_effect; ++ie){
					if (model.individual[0].ind_effect_value[ie] != UNSET){
						auto av = 0.0;
						for(auto ch = 0; ch < mpi.ncore; ch++){
							const auto &iepm = chain_ind_effect_posterior_mean[ch][i];
							av += iepm.ind_effect_sum[ie] / nind_effect_posterior_mean;
						}
						av /= mpi.ncore;
						
						ebvs_csv << ", " << av;
				
						//ebvs_csv << ", " << ind_effect_posterior_mean[i].ind_effect_sum[ie] / nind_effect_posterior_mean;
					}
				}
				ebvs_csv << '\n';
			}
			ebvs_csv.close();
		}

		// Calculates posterior parameter statistics

		ofstream post_csv(model.output_dir + "/" + "posterior.csv");

		post_csv << "parameter, mean, CI 95 min, CI 95 max, ESS, GR\n";
		for (auto th = 0; th < model.nparam; th++) {
			vector <double> vec;
			for(auto ch = 0; ch < chain_sample.size(); ch++){
				for (const auto &samp : chain_sample[ch]) vec.push_back(samp.param_value[th]);
			}
			
			auto stat = get_statistic(vec);

			post_csv << model.param[th].name << ", " << stat.mean << ", " << stat.CImin << ", "  << stat.CImax << ", ";
			post_csv << stat.ESS << ", ";
			
			if(mpi.ncore == 1) post_csv << "---";
			else{
				auto nrun = mpi.ncore;
				auto N = chain_sample[0].size();
				double mu[nrun], vari[nrun];
					
				auto muav = 0.0;
				for(auto ru = 0u; ru < nrun; ru++){ 
					auto valav = 0.0; 
					for(auto i = 0u; i < N; i++){
						valav += chain_sample[ru][i].param_value[th]/N;
					}
					
					auto varr = 0.0; 
					for(auto i = 0u; i < N; i++){
						auto val = chain_sample[ru][i].param_value[th];
						varr += (val-valav)*(val-valav)/(N-1);
					}
					
					mu[ru] = valav;
					vari[ru] = varr;
					muav += mu[ru]/nrun;
				}
				auto W = 0.0; for(auto ru = 0u; ru < nrun; ru++) W += vari[ru]/nrun;
				auto B = 0.0; for(auto ru = 0u; ru < nrun; ru++) B += (mu[ru]-muav)*(mu[ru]-muav)*N/(nrun-1);
				post_csv << sqrt(((1-1.0/N)*W + B/N)/W);
			}	
			post_csv << endl;
		}
		
		for (auto k = 0; k < model.nderived; k++) {
			const auto &der = model.derived[k];
			
			vector <double> vec;
			for(auto ch = 0; ch < chain_sample.size(); ch++){
				for (const auto &samp : chain_sample[ch]) vec.push_back(model.calculate_derived(der,samp.param_value));
			}
			
			auto stat = get_statistic(vec);

			post_csv << der.name << ", " << stat.mean << ", " << stat.CImin << ", "  << stat.CImax << ", ";
			post_csv << stat.ESS << ", ";
			
			post_csv << "---";
			post_csv << endl;
		}
			
		post_csv.close();
	}
	
	// Displays MCMC acceptance probabilities

	if(mpi.core == 0) fout << "Model parameter acceptance probabilities:" << endl << endl;
	for (auto th = 0; th < model.nparam; th++) {
		
		const auto &jump = param_jump[th];
		
		auto ac = stat(int(100.0 * jump.nac / jump.ntr));
		auto si = stat(jump.size);
		
		if(mpi.core == 0) fout << model.param[th].name << ":\t" << ac << "%\tSize: " << si << '\n';
	}

	if(mpi.core == 0) fout << '\n';

	if (model.nind_effect > 0) {
		if(mpi.core == 0) fout << "Individual effect acceptance probabilities:\n\n";
		for (auto ie = 0; ie < model.nind_effect; ie++) {
			const auto &jump = ind_effect_jump[ie];
			auto ac = stat(int(100.0 * jump.nac / jump.ntr));
			if(mpi.core == 0) fout << model.ind_effect[ie].name << ": " << ac << "%\n";
		}
		if(mpi.core == 0) fout << '\n' << endl;
	}

	auto flag = false;
	for (auto &jump : event_jump)
	{
		if (jump.ntr > 0) flag = true;
	}
	
	if (flag == true) {
		if(mpi.core == 0) fout << "Event time acceptance probabilities:\n\n";
		for (auto g = 0; g < model.ngroup; g++) {
			const auto &jump = event_jump[g];
			if(mpi.core == 0) fout << model.group[g].name << ":\t";
			
			if (jump.ntr == 0){
				if(mpi.core == 0) fout << "No proposals\n";
			}
			else {
				auto ac = stat(int(100.0 * jump.nac / jump.ntr));
				auto fa = stat(int(100.0 * jump.nevent_fa / jump.nevent_tr));
				auto si = stat(jump.size);
				if(mpi.core == 0){
					fout << ac << "%\t";
					fout << "Event failure: " << fa << "%\t";
					fout << "Changes per proposal: " << si << '\n';
				}
			}
		}
		if(mpi.core == 0) fout << '\n' << endl;
	}

	flag = false;
	for (auto &jump : add_rem_jump){
		if (jump.ntr > 0) flag = true;
	}
	
	if (flag == true) {
		if(mpi.core == 0) fout << "Add/remove acceptance probabilities:\n\n";
		for (auto g = 0; g < model.ngroup; g++) {
			const auto &jump = event_jump[g];
			if(mpi.core == 0) fout << model.group[g].name << ":\t";
			if (jump.ntr == 0){
				if(mpi.core == 0) fout << "No proposals\n";
			}
			else {
				auto ac = stat(int(100.0 * jump.nac / jump.ntr));
				auto fa = stat(int(100.0 * jump.nevent_fa / jump.nevent_tr));
				auto si = stat(jump.size);
				
				if(mpi.core == 0){
					fout << ac << "%   ";
					fout << "Failure: " << fa << "%\t";
					fout << "Changes per proposal: " << si << '\n';
				}
			}
		}
	}

	flag = false;
	for (auto &jump : mean_time_jump){
		if (jump.ntr > 0) flag = true;
	}
	
	if (flag == true) {
		if(mpi.core == 0) fout << "Mean Time acceptance probabilities:\n\n";
		
		for (auto mtp = 0; mtp < model.mean_time_prop.size(); mtp++) {
			const auto &jump = mean_time_jump[mtp];
			
			auto ac = stat(int(100.0 * jump.nac / jump.ntr));
			auto fa = stat(int(100.0 * jump.nfa / jump.ntr));
			auto si = stat(jump.size);
			
			if(mpi.core == 0){
				fout << model.mean_time_prop[mtp].name << ": " 
							<< ac << "%    fail: "
							<< fa << "%    size: "
							<< si << "\n";
			}
		}
		if(mpi.core == 0) fout << '\n' << endl;
	}
	
	flag = false;
	for (auto &jump : var_ie_joint_jump){
		if (jump.ntr > 0) flag = true;
	}
	
	if (flag == true) {
		if(mpi.core == 0) fout << "Joint inidividual effect / variance acceptance probabilities:\n\n";
		
		for (auto j = 0; j < var_ie_joint_jump.size(); j++) {
			const auto &jump = var_ie_joint_jump[j];
			
			auto ac = stat(int(100.0 * jump.nac / jump.ntr));
			auto si = stat(jump.size);
			
			if(mpi.core == 0){
				fout << model.ie_var_joint_prop[j].name << ": " 
							<< ac << "%    size: "
							<< si << "\n";
			}
		}
		if(mpi.core == 0) fout << '\n' << endl;
	}
	
	if(mpi.core == 0) fout << endl;
	
	if(false){
		ofstream ie_plot(model.output_dir + "/ie_plot.txt"); 
		for (auto i = 0; i < model.N; i++) {
			ie_plot << model.individual[i].id;
			for(auto ie = 0; ie < model.nind_effect; ie++){
				// Calculate the posterior mean and SD for the individual effect
				const auto &iepm = ind_effect_posterior_mean[i];
				auto PM = iepm.ind_effect_sum[ie] / nind_effect_posterior_mean;
				auto SD = sqrt(iepm.ind_effect_sum2[ie] / nind_effect_posterior_mean - PM*PM);
					
				// The actual value for the individual effect (from the data table)
				auto actual = model.individual[i].ind_effect_value[ie];
				
				ie_plot << " " << actual << " " << PM << " " << SD;
			}
			ie_plot << endl;
		}
	}	
}


// Calcualtes the range of a a statistic
string MCMC::stat(double num)
{
	auto num_tot = mpi.gather(num);
	
	sort(num_tot.begin(),num_tot.end());
	stringstream ss;
	ss << num_tot[num_tot.size()/2] << "(" << num_tot[0] << " - " << num_tot[num_tot.size()-1] << ")";
	
	return ss.str();
}


/// Generates a sample of the current state of the chain (used for generating diganositc information later)
void MCMC::store_sample() {
	// cout << "MCMC::store_sample()" << endl; // DEBUG
	Sample samp;
	samp.param_value = param_value;
	sample.push_back(samp);
}


/// Generates a sample of the likelihood (to aid calculation of quenching
void MCMC::store_burnin_Li() 
{
	auto sum = 0.0;
	if(model.cloglog.on == false){
		for (auto g = 0; g < model.ngroup; g++) sum += L_inf_events[g];
		sum += L_trans_events;
	}
	
	burnin_Li.push_back(sum);
	burnin_phi.push_back(phi);
}


/// Sums up individual effects (such that posterior averages can be caluculated later)
void MCMC::ind_effect_posterior_mean_update() {
	// cout << "MCMC::ind_effect_posterior_mean_update()" << endl; // DEBUG
	for (auto i = 0; i < model.N; i++) {
		for (auto ie = 0; ie < model.nind_effect; ie++){
			auto val = ind_value[i].ind_effect[ie];
			ind_effect_posterior_mean[i].ind_effect_sum[ie] += val;
			ind_effect_posterior_mean[i].ind_effect_sum2[ie] += val*val;
		}
	}
	nind_effect_posterior_mean++;
}


/// Initialises MCMC proposals
void MCMC::initialise_proposals() {
	// cout << "MCMC::initialise_proposals()" << endl; // DEBUG
	param_jump.resize(model.nparam);                                  // Initialises univariate proposals
	for (auto &jump : param_jump) {
		jump.size = 0.1;
		jump.nac = 0;
		jump.ntr = 0;
	}

	ind_effect_jump.resize(model.nparam);                             // Initialises individual effects proposals
	for (auto &jump : ind_effect_jump) {
		jump.nac = 0;
		jump.ntr = 0;
	}

	event_jump.resize(model.ngroup);                                  // Initialises event jumping proposals
	for (auto &jump : event_jump) {
		jump.size = 2;
		jump.nac = 0;
		jump.ntr = 0;
		jump.nevent_fa = 0;
		jump.nevent_tr = 0;
	}

	add_rem_jump.resize(model.ngroup);                                // Initialises add/remove jumping proposals
	for (auto &jump : add_rem_jump) {
		jump.size = 2;
		jump.nac = 0;
		jump.ntr = 0;
		jump.nevent_fa = 0;
		jump.nevent_tr = 0;
	}

	mean_time_jump.resize(model.ngroup);                  // Initialises mean time jumping proposals
	for (auto &jump : mean_time_jump) {
		jump.size = 0.01;
		jump.nac = 0;
		jump.ntr = 0;
		jump.nfa = 0;
	}
	
	var_ie_joint_jump.resize(model.ie_var_joint_prop.size());                  // Initialises mean time jumping proposals
	for (auto &jump : var_ie_joint_jump) {
		jump.size = 0.01;
		jump.nac = 0;
		jump.ntr = 0;
		jump.nfa = 0;
	}
	
	for (auto i = 0; i < model.N; i++) {                             // Initialises infection sampler for add/rem
		auto &inf_samp = ind_value[i].inf_sampler;
		auto &indiv = model.individual[i];
		if (indiv.status == UNKNOWN) {
			inf_samp.on = true;
			auto &timer = indiv.trans_time_range;
			inf_samp.tmin = timer[0].tmin;
			inf_samp.tmax = timer[0].tmax;
			if (inf_samp.tmax == LARGE)
				emsg("The inference time range is too large for infection sampler");

			inf_samp.bin.resize(nbin);
			inf_samp.log_prob.resize(nbin);
			inf_samp.prob_sum.resize(nbin);
			for (auto b = 0; b < nbin; b++)
				inf_samp.bin[b] = 10;
		} else
			inf_samp.on = false;
	}

	initialise_inf_sampler(0);
}


/// Initialises samplers used to sample infection time (this is used to add and removed infected individuals)
void MCMC::initialise_inf_sampler(int s) {
	// cout << "MCMC::initialise_inf_sampler()" << endl; // DEBUG
	
	if(model.cloglog.on == true) return;
	
	timer[TIME_INF_SAMP_UPDATE].start();
	for (auto &ind : ind_value) {
		auto &inf_samp = ind.inf_sampler;
		if (inf_samp.on == true) {
			if (ind.infected == true) {
				auto b = int(nbin * (ind.trans_time[0] - inf_samp.tmin) / (inf_samp.tmax - inf_samp.tmin));
				if (b < 0 || b >= nbin)
					emsg("Infection outside of sampler range");

				inf_samp.bin[b]++;
			}

			// Sets up the sampler
			if (s % 100 == 0) {
				vector <double> bin(nbin);
				auto total = 0.0;
				for (auto b = 0; b < nbin; b++) {
					bin[b] = inf_samp.bin[b];
					total += bin[b];
				}

				// Adds some over dispersion
				auto add = 0.3 * total / nbin;
				for (auto b = 0; b < nbin; b++)
					bin[b] += add;
				total += add * nbin;

				auto factor = nbin / (inf_samp.tmax - inf_samp.tmin);
				auto sum = 0.0;
				for (auto b = 0; b < nbin; b++) {
					auto prob = bin[b] / total;
					inf_samp.log_prob[b] = log(prob * factor);
					sum += prob;
					inf_samp.prob_sum[b] = sum;
				}
			}
		}
	}

	timer[TIME_INF_SAMP_UPDATE].stop();
}


/// Checks the effect of quenching
void MCMC::check_quenching()
{
	auto chain_sample = mpi.gather_sample(sample);
	
	if(mpi.core == 0){
	
		for (auto th = 0; th < model.nparam; th++) {
			ofstream opt(model.param[th].name+".txt");
			
			for(auto ch = 0; ch < chain_sample.size(); ch++){
				opt << ch << " ";
				opt << phi_ch[ch] << " ";
				
				vector <double> vec;
				for (const auto &samp : chain_sample[ch])
				vec.push_back(samp.param_value[th]);
				auto stat = get_statistic(vec);
		
				opt << stat.mean <<" " << stat.CImin << " " << stat.CImax << " " << stat.ESS << "\n";
			}
		}
	}
}

