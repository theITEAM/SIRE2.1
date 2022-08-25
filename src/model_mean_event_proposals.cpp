/// Proposals for joint mean and event updates

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "model.hpp"
#include "timers.hpp"
#include "utils.hpp"

/// This makes simultaneous changes to event times and distribution means
void Model::propose_mean_event_times(vector <IndValue> &ind_value, vector <double> &param_value, vector <double> &L_inf_events, double &L_trans_events, vector <double> &L_diag_test, double &prior, vector <Jump> &mean_time_jump, const bool burnin, const Quench &quench) const 
{
	timer[TIME_MEAN_EVENTS].start();
	
	for(auto mtp = 0; mtp < mean_time_prop.size(); mtp++){
		auto tr = mean_time_prop[mtp].trans_move;
		auto &jump = mean_time_jump[mtp];
		
		if(tr > 0 && tr < ntrans-1){
			auto th1 = trans[tr].mean_param;
			auto th2 = trans[tr+1].mean_param;
			
			auto lam1 = param_value[th1];
			auto lam2 = param_value[th2];
			auto d = normal_sample(0,jump.size);
			
			auto lam1_prop = lam1 + d;
			auto lam2_prop = lam2 - d;
			
			auto param_prop = param_value;
			
			param_prop[th1] = lam1_prop;
			param_prop[th2] = lam2_prop;
			
			jump.ntr++;
			if(inbounds(param_prop) == false) jump.nfa++;
			else{
				auto fac = lam1_prop/lam1;
				
				auto ind_value_prop = ind_value;
				
				auto flag = false;
				auto num = 0;
		
				for(auto &ind : ind_value_prop){
					auto t = ind.trans_time[tr];
					if(t != UNSET){
						num++;
						
						auto t_b = ind.trans_time[tr-1];
						auto t_a = ind.trans_time[tr+1];
						auto t_prop = t_b + fac*(t - t_b);
						if(t_prop > t_a) flag = true;
						
						ind.trans_time[tr] = t_prop;
					}
				}
			
				if(flag == true) jump.nfa++;
				else{
					set_individual_quantities(ind_value_prop, param_prop);
					
					auto prior_prop = calculate_prior(param_prop);
		
					auto sum = 0.0;
					auto L_inf_events_prop = calculate_L_inf_events(ind_value_prop, param_prop);
					for (auto g = 0; g < ngroup; g++) sum += L_inf_events_prop[g] - L_inf_events[g];
			
					auto L_trans_events_prop = calculate_L_trans_events(ind_value_prop, param_prop);
					sum += L_trans_events_prop - L_trans_events;
			
					//auto val = quench.phi_L*sum + quench.phi_Pr*(prior_prop-prior);
					//if(val > 10) val = 10;
					//auto al = exp(val);
				
					auto al = exp(quench.phi_L*sum + quench.phi_Pr*(prior_prop-prior) + num*log(fac));

					if (MH_proposal(al,12)) {
						jump.nac++;

						param_value = param_prop;
						ind_value = ind_value_prop;

						L_inf_events = L_inf_events_prop;
						L_trans_events = L_trans_events_prop;
						prior = prior_prop;
						
						if (burnin == true) jump.size *= prop_up;
					} 
					else{
						if (burnin == true) jump.size *= prop_down;
					}
				}
			}
		}
	}

	timer[TIME_MEAN_EVENTS].stop();
}



/// This works how what type of mean/event time proposals can be made
void Model::propose_mean_event_times_initialise()
{
	vector <bool> cannot_move(ntrans,false);
	
	for(auto tr = 0; tr < ntrans; tr++){
		auto cannot_move = false;
		
		for(const auto &ind : individual){
			auto &tra = ind.trans_time_range[tr];
			if(tra.tmin == tra.tmax && tra.tmin != UNSET) cannot_move = true;
		}
		
		if(cannot_move == false){
			MeanTimeProposal m_prop;
			m_prop.name = "Mean time proposal for "+trans[tr].name;
			m_prop.trans_move = tr;
			
			mean_time_prop.push_back(m_prop);
		}
	}
}
