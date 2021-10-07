/// This provides functions for the MCMC chain

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "mcmc.hh"
#include "model.hh"
#include "check.hh"
#include "utils.hh"
#include "timers.hh"

/// This initialises the mcmc chain
MCMC::MCMC(const Model &model) : model(model)
{
	cout << "Initialising MCMC chain..." << endl;
	
  trace_initialise("trace.txt");                                    // Initialises the trace plot
  
  param_value = model.prior_sample();                               // Parameters sampled from prior
  
  ind_effect_posterior_mean.resize(model.N);                        // Initialises individual effect means
  for(auto &ie : ind_effect_posterior_mean){
    ie.ind_effect_sum.resize(model.nind_effect);
    for(auto &val : ie.ind_effect_sum) val = 0;
  }
  nind_effect_posterior_mean = 0;
  
  ind_value.resize(model.N);                                        // Sets a variable for individuals information
  for(auto i = 0u; i < model.N; i++) ind_value[i].index = i;        // Sets index for each individual
  
  model.ind_effect_sample(ind_value,param_value);                   // Individual effects sampled from prior

  model.set_individual_quantities(ind_value,param_value);           // Sets individual transition parameters and infectivity
  
  model.set_initial_events(ind_value,param_value);                  // Sample initial individual events
	
	initialise_proposals();                                           // Initialises MCMC proposals
}


/// Sets all the likelihoods at the start of the chain
void MCMC::set_likelihoods()
{
  L_ind_effect = model.calculate_L_ind_effect(ind_value,param_value);// Sets likelihoods for individual effects 

  L_inf_events = model.calculate_L_inf_events(ind_value,param_value);// Sets likelihoods for infection events
  
  L_trans_events = model.calculate_L_trans_events(ind_value,param_value);// Sets likelihoods for infection events
    
	L_diag_test = model.calculate_L_diag_test(ind_value,param_value);  // The likelihood for diagnostic test results
	
  prior = model.calculate_prior(param_value);                        // Sets the prior
  
  check_chain();
}


/// Performs a series of MCMC proposal which consist of an "update" 
void MCMC::update()
{
  model.propose_transmission_rate(ind_value,param_value,L_inf_events,prior,param_jump);

  for(auto c = 0u; c < model.ncovariance; c++){
    model.propose_covariance_matrices(model.covariance[c],ind_value,param_value,L_ind_effect[c],prior,param_jump,burnin);
  }

  model.propose_trans_params(ind_value,param_value,L_trans_events,prior,param_jump,burnin);

  if(model.group_effect.on == true){
    model.propose_group_effect_sigma(param_value,prior,param_jump,burnin);
    model.propose_group_effect(ind_value,param_value,L_inf_events,prior,param_jump,burnin);
  }

  model.propose_fixed_effects(ind_value,param_value,L_inf_events,L_trans_events,prior,param_jump,burnin);

	model.propose_snp_effects(ind_value,param_value,L_inf_events,L_trans_events,prior,param_jump,burnin);
 
  model.propose_event_times(ind_value,param_value,L_inf_events,L_trans_events,L_diag_test,event_jump,burnin);

  model.propose_add_rem(ind_value,param_value,L_inf_events,L_trans_events,L_diag_test,add_rem_jump,burnin);

  model.propose_susceptibility_ind_effects(ind_value,param_value,L_ind_effect,L_inf_events,ind_effect_jump);

	model.propose_infectivity_ind_effects(ind_value,param_value,L_ind_effect,L_inf_events,ind_effect_jump);

  model.propose_trans_ind_effects(ind_value,param_value,L_ind_effect,L_trans_events,ind_effect_jump);

	model.propose_Se_Sp(ind_value,param_value,L_diag_test,prior,param_jump,burnin); 
	
  check_chain();
}


/// Intialises the trace plot
void MCMC::trace_initialise(const string file)
{
  trace.open(model.output_dir+"/"+file);
  trace << "state";
  for(auto par : model.param) trace << "\t" << par.name;
  for(auto c = 0u; c < model.ncovariance; c++) trace << "\tL_ind_effect " << c;
  trace << "\tL_inf_events";
  trace << "\tL_trans_events";
	trace << "\tL_diag_test";
  trace << "\tPrior";
  trace << "\tNumber infected";
  trace << endl;
}


/// Outputs the trace plot
void MCMC::trace_output(const int s)
{
  trace << s;
  for(auto val : param_value) trace << "\t" << val;
  for(auto c = 0u; c < model.ncovariance; c++) trace << "\t" << L_ind_effect[c];
  auto sum = 0.0; for(auto g = 0u; g < model.ngroup; g++) sum += L_inf_events[g];
  trace << "\t" << sum;
  trace << "\t" << L_trans_events;
	sum = 0.0; for(auto g = 0u; g < model.ngroup; g++) sum += L_diag_test[g];
	trace << "\t" << sum;
  trace << "\t" << prior;
  auto num_infected = 0u; for(const auto &ind : ind_value){ if(ind.infected == true) num_infected++;}
  trace << "\t" << num_infected;
  trace << endl;
}


/// Given a vector of values this returns statistical information
Statistics MCMC::get_statastic(const vector <double> &vec) const              
{
  Statistics stat;
  
  auto n = vec.size();
  if(n == 0){
    stat.mean = "---"; stat.CImin = "---"; stat.CImax = "---"; stat.ESS = "---"; 
  }
  else{
    auto sum = 0.0, sum2 = 0.0; 
    for(auto i = 0u; i < vec.size(); i++){ sum += vec[i]; sum2 += vec[i]*vec[i];}
    sum /= n; sum2 /= n;
    
    stat.mean = to_string(sum); 
    
    vector <double> vec2 = vec;
    sort(vec2.begin(),vec2.end());
  
    if(n >= 2){
      auto i = (unsigned int)((n-1)*0.025); auto f = (n-1)*0.025 - i;
      stat.CImin = to_string(vec2[i]*(1-f) + vec2[i+1]*f);
        
      i = (unsigned int)((n-1)*0.975); f = (n-1)*0.975 - i;
      stat.CImax = to_string(vec2[i]*(1-f) + vec2[i+1]*f);
    }
    else{
      stat.CImin = to_string(vec2[0]);
      stat.CImax = to_string(vec2[0]);
    }

    vec2 = vec;
    auto var = sum2 - sum*sum;
    if(var <= 0.0000000001 || n <= 2) stat.ESS = "---";
    else{ 
      auto sd = sqrt(var);
      for(auto i = 0u; i < vec2.size(); i++) vec2[i] = (vec2[i]-sum)/sd;
        
      auto sum = 1.0;
      for(auto d = 1u; d < n/2; d++){             // Calculates the effective sample size
        auto a = 0.0; for(auto i = 0u; i < n-d; i++) a += vec2[i]*vec2[i+d]; 
        auto cor = a/(n-d);
        if(cor < 0) break;
        sum += 2*cor;     
      }
      stat.ESS = to_string(int(n/sum));
    }
  }
  
  return stat;
}
  
  
/// Outputs statisitcs at the end of the run
void MCMC::output_statistics(const string file) const
{
  ofstream fout(model.output_dir+"/"+file);
  
  // Calculates prediction accuracies
  
  if(model.pred_acc.size() > 0) fout << "Prediction accuracies:" << endl << endl;
  for(const auto &pa : model.pred_acc){
    fout << pa.name << "   ";
    
    for(auto ie = 0u; ie < model.nind_effect; ie++){
      if(model.individual[0].ind_effect_value[ie] != UNSET){
        auto av_PM = 0.0, av_PM2 = 0.0;
        auto av_actual = 0.0, av_actual2 = 0.0;
        auto av_PM_actual = 0.0;
        auto nav = 0.0;
        
        for(auto i : pa.ind){
          // Calculate the posterior mean for the individual effect
          auto PM = ind_effect_posterior_mean[i].ind_effect_sum[ie]/nind_effect_posterior_mean;
    
          // The actual value for the individual effect (from the data table)
          auto actual = model.individual[i].ind_effect_value[ie];
          
          av_PM += PM; av_PM2 += PM*PM;
          av_actual += actual; av_actual2 += actual*actual;
          av_PM_actual += PM*actual;
          nav++;
        }
        av_PM /= nav; av_PM2 /= nav;
        av_actual /= nav; av_actual2 /= nav;
        av_PM_actual /= nav;
        
        auto var_PM = av_PM2 - av_PM*av_PM;
        auto var_actual = av_actual2 - av_actual*av_actual;
        
        // Calculates the Pearson correlation coefficient
        auto cor = (av_PM_actual - av_PM*av_actual)/sqrt(var_PM*var_actual);
        
        fout << model.ind_effect[ie].name << ": " << cor << "  ";
      }
    }
    fout << endl;
  }
  fout << endl << endl;
  
  // Calculates parameter statistics
  
  fout << "Posterior parameter posterior estimates:" << endl << endl;
  fout << "Name\tMean\t95% CI (min - max)\tESS" << endl;
  for(auto th = 0u; th < model.nparam; th++){
    vector <double> vec; for(const auto &samp : sample) vec.push_back(samp.param_value[th]);
    auto stat = get_statastic(vec);
    
    fout << model.param[th].name << "\t" << stat.mean << "\t"
          << stat.CImin << " - "  << stat.CImax << "\t" << stat.ESS << endl; 
  }
  fout << endl;
  fout << "ESS above 200 for every parameter indicates convergence" << endl;
  fout << endl << endl;
  
  // Displays MCMC acceptance probabilities
  
  fout << "Model parameter acceptance probabilities:" << endl << endl;
  for(auto th = 0u; th < model.nparam; th++){
    const auto &jump = param_jump[th];
    fout << model.param[th].name << ": " << int(100.0*jump.nac/jump.ntr) << "%   Size: " << jump.size << endl;
  } 
  
  fout << endl;
  
  if(model.nind_effect > 0){
    fout << "Individual effect acceptance probabilities:" << endl << endl;
    for(auto ie = 0u; ie < model.nind_effect; ie++){
      const auto &jump = ind_effect_jump[ie];
      fout << model.ind_effect[ie].name << ": " << int(100.0*jump.nac/jump.ntr) << "%" << endl;
    } 
  }
  
  auto flag = false;
  for(auto &jump : event_jump) if(jump.ntr > 0) flag = true; 
  if(flag == true){
    fout << endl << endl;
    fout << "Event time acceptance probabilities:" << endl << endl;
    for(auto g = 0u; g < model.ngroup; g++){
      const auto &jump = event_jump[g];
      fout << model.group[g].name << ":    ";
      if(jump.ntr == 0) fout << "No proposals" << endl;
      else{
        fout << int(100.0*jump.nac/jump.ntr) << "%   ";
        fout << "Event failure: " << int(100.0*jump.nevent_fa/jump.nevent_tr) << "%   ";
        fout << "Changes per proposal: " << jump.size << endl;
      }
    } 
  }
  
  flag = false;
  for(auto &jump : add_rem_jump) if(jump.ntr > 0) flag = true; 
  if(flag == true){
    fout << endl << endl;
    fout << "Add/remove acceptance probabilities:" << endl << endl;
    for(auto g = 0u; g < model.ngroup; g++){
      const auto &jump = event_jump[g];
      fout << model.group[g].name << ":    ";
      if(jump.ntr == 0) fout << "No proposals" << endl;
      else{
        fout << int(100.0*jump.nac/jump.ntr) << "%   ";
        fout << "Failure: " << int(100.0*jump.nevent_fa/jump.nevent_tr) << "%   ";
        fout << "Changes per proposal: " << jump.size << endl;
      }
    } 
  }
}


/// Generates a sample of the current state of the chain (used for generating diganositc information later)
void MCMC::store_sample()
{
  Sample samp;
  samp.param_value = param_value;
  sample.push_back(samp);
}


/// Sums up individual effects (such that posterior averages can be caluculated later)
void MCMC::ind_effect_posterior_mean_update()
{
  for(auto i = 0u; i < model.N; i++){
    for(auto ie = 0u; ie < model.nind_effect; ie++){
      ind_effect_posterior_mean[i].ind_effect_sum[ie] += ind_value[i].ind_effect[ie];
    }   
  }
  nind_effect_posterior_mean++;
}


/// Initialises MCMC proposals
void MCMC::initialise_proposals()
{
  param_jump.resize(model.nparam);                                  // Initialises univariate proposals
  for(auto &jump : param_jump){
    jump.size = 0.1; jump.nac = 0; jump.ntr = 0;
  }
  
  ind_effect_jump.resize(model.nparam);                             // Initialises individual effects proposals
  for(auto &jump : ind_effect_jump){
    jump.nac = 0; jump.ntr = 0;
  }
  
  event_jump.resize(model.ngroup);                                  // Initialises event jumping proposals
  for(auto &jump : event_jump){
    jump.size = 2; jump.nac = 0; jump.ntr = 0; jump.nevent_fa = 0; jump.nevent_tr = 0;
  }
  
  add_rem_jump.resize(model.ngroup);                                // Initialises add/remove jumping proposals
  for(auto &jump : add_rem_jump){
    jump.size = 2; jump.nac = 0; jump.ntr = 0; jump.nevent_fa = 0; jump.nevent_tr = 0;
  }
	
	for(auto i = 0u; i < model.N; i++){                               // Initialises infection sampler for add/rem
		auto &inf_samp = ind_value[i].inf_sampler;
		auto &indiv = model.individual[i];
		if(indiv.status == UNKNOWN){	
			inf_samp.on = true;
			auto &timer = indiv.trans_time_range;
      inf_samp.tmin = timer[0].tmin;
			inf_samp.tmax = timer[0].tmax;  
			if(inf_samp.tmax == LARGE) emsg("The inference time range is too large for infection sampler");
		
			inf_samp.bin.resize(nbin);
			inf_samp.log_prob.resize(nbin);
			inf_samp.prob_sum.resize(nbin);
			for(auto b = 0u; b < nbin; b++) inf_samp.bin[b] = 10;
		}
		else inf_samp.on = false;
	}
	
	initialise_inf_sampler(0);
}


/// Initialises samplers used to sample infection time (this is used to add and removed infected individuals)
void MCMC::initialise_inf_sampler(int s)
{
	timer[TIME_INF_SAMP_UPDATE].start();
	for(auto &ind : ind_value){
		auto &inf_samp = ind.inf_sampler;
		if(inf_samp.on == true){
			if(ind.infected == true){
				auto b = int(nbin*(ind.trans_time[0] - inf_samp.tmin)/(inf_samp.tmax - inf_samp.tmin));
				if(b < 0 || b >= nbin) emsg("Infection outside of sampler range");
				
				inf_samp.bin[b]++;
			}
			
			// Sets up the sampler
			if(s%100 == 0){
				vector <double> bin(nbin);
				auto total = 0.0;
				for(auto b = 0u; b < nbin; b++){
					bin[b] = inf_samp.bin[b];
					total += bin[b];
				}
				
				// Adds some over dispersion
				auto add = 0.3*total/nbin;  
				for(auto b = 0u; b < nbin; b++) bin[b] += add; total += add*nbin;
				
				auto factor = nbin/(inf_samp.tmax - inf_samp.tmin);
				auto sum = 0.0;
				for(auto b = 0u; b < nbin; b++){
					auto prob = bin[b]/total;
					inf_samp.log_prob[b] = log(prob*factor);					
					sum += prob;
					inf_samp.prob_sum[b] = sum; 
				}
			}
		}
	}
	
	timer[TIME_INF_SAMP_UPDATE].stop();
}
