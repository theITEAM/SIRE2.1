/// This stores functions giving proposals on individual effects

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "model.hh"
#include "timers.hh"
#include "utils.hh"


/// This makes proposals to susceptibility individual effects
void Model::propose_susceptibility_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_inf_events, vector <Jump> &ind_effect_jump) const
{ 
  const auto &ind_eff = trans[0].ind_effect;
  auto nie = ind_eff.size();
  if(nie == 0) return;
    
  timer[TIME_SUS_IND_EFFECT_INIT].start();    
  auto inv_cov_matrix = calculate_inv_cov_matrix(param_value);
  auto precalc = set_precalc_likelihood_sus(ind_value,param_value);
  timer[TIME_SUS_IND_EFFECT_INIT].stop();
  
  timer[TIME_SUS_IND_EFFECT].start();
  vector <PrecalcIndEff> prec_ind_eff(nie);
    
  for(auto i = 0u; i < N; i++){
    auto &prec = precalc[i];
    auto &ind = ind_value[i];
    auto &indiv = individual[i];
    auto g = indiv.group;
    
    for(auto j = 0u; j < nie; j++){
      set_precalc_ind_eff(i,prec_ind_eff[j],ind_eff[j],ind_value,inv_cov_matrix);
    }
     
    auto sus_fac = log(ind.susceptibility);
    for(auto loop = 0; loop < ind_effect_sus_prop; loop++){
      for(auto j = 0u; j < nie; j++){
        const auto &prec_ie = prec_ind_eff[j];
        auto ie = ind_eff[j];
        
        auto ind_eff_store = ind.ind_effect[ie];
        auto ind_eff_prop = normal_sample(prec_ie.mean,prec_ie.sd);
        
        auto sus_fac_store = sus_fac;
        sus_fac += ind_eff_prop - ind_eff_store;
          
        auto dL_inf_ev = Li_change_sus(sus_fac_store,sus_fac,prec);
        
        auto al = exp(dL_inf_ev);
        
        ind_effect_jump[ie].ntr++;
        if(ran() < al){
          ind_effect_jump[ie].nac++;
          
          ind.ind_effect[ie] = ind_eff_prop;
          
          auto cv = ind_effect[ie].covar_ref;
          L_ind_effect[cv] += Li_change_ind_eff(ind_eff_store,ind_eff_prop,prec_ie);
					
          if(g != UNSET) L_inf_events[g] += dL_inf_ev;
        }
        else{
          sus_fac = sus_fac_store;
        }
      }
    }
    
    ind.susceptibility = exp(sus_fac);
  }
  timer[TIME_SUS_IND_EFFECT].stop();  
}


/// This pre-calculates quantities so likelihood can be calculated from susceptibility 
vector <PrecalcLiSus> Model::set_precalc_likelihood_sus(const vector <IndValue> &ind_value, const vector <double> &param_value) const
{
  vector <PrecalcLiSus> prec(N);
  for(auto &pre : prec){ pre.beta_fac = 0; pre.log_beta_num = 0;}
  
  vector <int> infected(N);
  for(auto &inf : infected) inf = 0;
  
  for(const auto &gr : group){
    auto I = 0.0;                                           // The total infectivity
    auto S = 0.0;                                           // The total susceptibility
    
    auto sus_list = gr.ind_ref;                             // Makes a list of susceptible individuals
    
    vector <Event> event;
    get_event_sequence(event,S,gr,ind_value);               // Gets a sorted list of events for the entire group
      
    auto tmin = gr.inference_range.tmin;
    auto tmax = gr.inference_range.tmax;

    auto beta = param_value[beta_param];
    if(inf_model == FREQ_DEP) beta /= gr.nind;              // If frequency dependent divided by group size
    if(group_effect.on == true){                            // Adds in the group effect
      beta *= exp(param_value[group_effect.param[gr.index]]);
    }
    
    auto t = tmin;                                          // Sets a time variable
    for(const auto &ev : event){
      const auto &ind = ind_value[ev.ind];
      auto tev = ev.time;
      auto tr = ev.trans;
      
      auto val = -beta*I*(tev-t);
      for(auto i : sus_list) if(infected[i] == 0) prec[i].beta_fac += val;    
      t = tev;  
      
      if(tr == 0){
        auto i = ev.ind;
        if(tev > tmin) prec[i].log_beta_num++;      
        infected[i] = 1; 
      }
      
      I += ind.trans_infectivity_change[tr];
    }
    
    auto val = -beta*I*(tmax-t);
    for(auto i : sus_list) if(infected[i] == false) prec[i].beta_fac += val;
  }
  
  return prec;
}


/// Calculates the change in event likelihood as a result of susceptibility changing
double Model::Li_change_sus(const double sus_before, const double sus_after, PrecalcLiSus &prec) const
{
  return prec.log_beta_num*(sus_after-sus_before) + prec.beta_fac*(exp(sus_after)-exp(sus_before));
}


/// This pre-calculates quantities so changes in individual effect likelihood is fast
void Model::set_precalc_ind_eff(const int i, PrecalcIndEff &prec, int ie, const vector <IndValue> &ind_value, const vector <InvCovMat> &inv_cov_matrix) const
{
  auto val_sq_fac = 0.0, val_fac = 0.0;
  
  const auto cv = ind_effect[ie].covar_ref;
  const auto &cov = covariance[cv];
  const auto e_val = ind_effect[ie].covar_num;
  const auto &M = inv_cov_matrix[cv].M; 
  const auto &mat = matrix[cov.matrix];
  auto E = cov.E;

  const auto &ie_ref = cov.ind_effect_ref;

  auto val = mat.Ainvdiag[i];
  val_sq_fac -= 0.5*val*M[e_val][e_val];
  
  const auto &ind_eff = ind_value[i].ind_effect;
  auto sum = 0.0;
  for(auto ei = 0u; ei < E; ei++){
    if(ei != e_val) sum += M[e_val][ei]*ind_eff[ie_ref[ei]];
  }
  val_fac -= val*sum;
  
  for(const auto &ele : mat.Ainvlist[i]){
    const auto &ind_eff = ind_value[ele.i].ind_effect;
    auto sum = 0.0;
    for(auto ei = 0u; ei < E; ei++){
      sum += M[e_val][ei]*ind_eff[ie_ref[ei]];
    }
    val_fac -= ele.val*sum;
  }
    
  prec.val_sq_fac = val_sq_fac; prec.val_fac = val_fac;
  
  // Calculates the mean and standard deviation in the distribution
  auto var = -1.0/(2*val_sq_fac);           
  prec.sd = sqrt(var);
  prec.mean = val_fac*var;
}


/// Calculates the change in individual effect likelihood
double Model::Li_change_ind_eff(const double before, const double after, const PrecalcIndEff &prec) const
{
  return prec.val_sq_fac*(after*after-before*before) +  prec.val_fac*(after-before);
}


/// This makes proposals to transition individual effects
void Model::propose_trans_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, double &L_trans_events, vector <Jump> &ind_effect_jump) const
{
  // Checks if individual effects exist 
  auto exist = false; for(auto tr = 1u; tr < ntrans; tr++){ if(trans[tr].ind_effect.size() > 0) exist = true;}
  if(exist == false) return;
  
  timer[TIME_TRANS_IND_EFFECT].start();
  auto inv_cov_matrix = calculate_inv_cov_matrix(param_value);
  
	// This changes individual effect for a given transition				
					
  for(auto tr = 1u; tr < ntrans; tr++){
    const auto &tra = trans[tr];
    const auto &ind_eff = tra.ind_effect;
   
    if(tra.all_ind_effect_proposal == true){
			auto nie = ind_eff.size();
      vector <PrecalcIndEff> prec_ind_eff(nie);
      
      auto mean_av = param_value[tra.mean_param];
      auto shape = param_value[tra.shape_param];
              
      double L, L_prop, dt;
      
      for(auto i = 0u; i < N; i++){
        auto &ind = ind_value[i];
        auto &indiv = individual[i];
        
        for(auto j = 0u; j < nie; j++){
          set_precalc_ind_eff(i,prec_ind_eff[j],ind_eff[j],ind_value,inv_cov_matrix);
        }
          
        auto mean = ind.trans_mean[tr];
        auto mean_fac = log(mean/mean_av);
        if(ind.infected == true){
          dt = ind.trans_time[tr]-ind.trans_time[tr-1];
          L = gamma_probability(dt,mean,shape);
        }
        
        for(auto loop = 0; loop < ind_effect_trans_prop; loop++){
          for(auto j = 0u; j < nie; j++){
            const auto &prec_ie = prec_ind_eff[j];
            auto ie = ind_eff[j];
            
            auto ind_eff_store = ind.ind_effect[ie];
            auto ind_eff_prop = normal_sample(prec_ie.mean,prec_ie.sd);
            
            auto mean_fac_store = mean_fac;
            mean_fac += ind_eff_prop - ind_eff_store;
            
            auto al = 1.0;
            if(ind.infected == true){
              L_prop = gamma_probability(dt,mean_av*exp(mean_fac),shape);
              al = exp(L_prop-L);
            }
          
            ind_effect_jump[ie].ntr++;
            if(ran() < al){
              ind_effect_jump[ie].nac++;
              
              ind.ind_effect[ie] = ind_eff_prop;
              
              auto cv = ind_effect[ie].covar_ref;
              L_ind_effect[cv] += Li_change_ind_eff(ind_eff_store,ind_eff_prop,prec_ie);
						
              if(ind.infected == true){
                L_trans_events += L_prop-L;
                L = L_prop;
              }
            }
            else{
              mean_fac = mean_fac_store;
            }
          }
        }
        
        ind.trans_mean[tr] = mean_av*exp(mean_fac);
      }
    }
  }

	// This changes individual effects for multiple transitions
	for(auto ie = 0u; ie < nind_effect; ie++){
		const auto &ind_eff = ind_effect[ie];
		const auto &tr_mean = ind_eff.trans_mean;
		auto ntr_mean = tr_mean.size();
		if(ntr_mean > 0 && ind_eff.single_ind_effect_proposal == true){
			PrecalcIndEff prec_ind_eff;
		
			vector <double> mean_av(ntrans), shape(ntrans), mean_fac(ntrans), mean_fac_store(ntrans), dt(ntrans);
			
			for(auto tr : tr_mean){
			  mean_av[tr] = param_value[trans[tr].mean_param];
        shape[tr] = param_value[trans[tr].shape_param];
			}
			
			for(auto i = 0u; i < N; i++){
        auto &ind = ind_value[i];
        auto &indiv = individual[i];
        
        set_precalc_ind_eff(i,prec_ind_eff,ie,ind_value,inv_cov_matrix);
        
				auto L = 0.0;
			
				for(auto tr : tr_mean){
					auto mean = ind.trans_mean[tr];
					mean_fac[tr] = log(mean/mean_av[tr]);
					if(ind.infected == true){
						dt[tr] = ind.trans_time[tr]-ind.trans_time[tr-1];
						L += gamma_probability(dt[tr],mean,shape[tr]);
					}
				}
        
        for(auto loop = 0; loop < ind_effect_trans_prop2; loop++){  
          auto ind_eff_store = ind.ind_effect[ie];
          auto ind_eff_prop = normal_sample(prec_ind_eff.mean,prec_ind_eff.sd);
          auto change_ind_eff = ind_eff_prop - ind_eff_store;
					
					for(auto tr : tr_mean){
            mean_fac_store[tr] = mean_fac[tr];
            mean_fac[tr] += change_ind_eff;
					}
            
					if(false){ // Checks for correct sampling from individual effect distribution
						auto probfi = normal_probability(ind_eff_store,prec_ind_eff.mean,prec_ind_eff.sd);
						auto probif = normal_probability(ind_eff_prop,prec_ind_eff.mean,prec_ind_eff.sd);
						cout << probfi - probif + Li_change_ind_eff(ind_eff_store,ind_eff_prop,prec_ind_eff) << "zero" << endl;
					}
					
					auto al = 1.0;
					
					auto L_prop = 0.0;
					if(ind.infected == true){
						for(auto tr : tr_mean) L_prop += gamma_probability(dt[tr],mean_av[tr]*exp(mean_fac[tr]),shape[tr]);
						al = exp(L_prop-L);
					}
				
					ind_effect_jump[ie].ntr++;
					if(ran() < al){
						ind_effect_jump[ie].nac++;
						
						ind.ind_effect[ie] = ind_eff_prop;
						
						auto cv = ind_effect[ie].covar_ref;
						L_ind_effect[cv] += Li_change_ind_eff(ind_eff_store,ind_eff_prop,prec_ind_eff);
						if(ind.infected == true){
							L_trans_events += L_prop-L;
							L = L_prop;
						}
					}
					else{
						for(auto tr : tr_mean) mean_fac[tr] = mean_fac_store[tr];
					}
				}
        
        for(auto tr : tr_mean) ind.trans_mean[tr] = mean_av[tr]*exp(mean_fac[tr]);
      }		
		}
	}
  timer[TIME_TRANS_IND_EFFECT].stop();  
}


/// This makes proposals to infectivity individual effects on a given compartment
void Model::propose_infectivity_ind_effects(vector <IndValue> &ind_value, const vector <double> &param_value, vector <double> &L_ind_effect, vector <double> &L_inf_events, vector <Jump> &ind_effect_jump) const
{ 
  // Checks if individual effects exist 
  auto exist = false; for(const auto &co : comp){ if(co.ind_effect.size() > 0) exist = true;}
  if(exist == false) return;
  
  timer[TIME_INF_IND_EFFECT_INIT].start();    
  auto inv_cov_matrix = calculate_inv_cov_matrix(param_value);
  vector <vector <double> > I_profile;
  auto precalc = set_precalc_likelihood_inf(I_profile,ind_value,param_value);
  timer[TIME_INF_IND_EFFECT_INIT].stop();
  
  timer[TIME_INF_IND_EFFECT].start();
	
	// This changes infectivity individual effect for a given compartment
	for(auto c = 1u; c < ncomp; c++){
    const auto &co = comp[c];
    const auto &ind_eff = co.ind_effect;
		
    if(co.all_ind_effect_proposal == true){  
		  auto nie = ind_eff.size();
      vector <PrecalcIndEff> prec_ind_eff(nie);
  
      auto inf_av = co.infectivity;
  
      for(auto i = 0u; i < N; i++){
        auto &prec = precalc[i][c];
        auto &ind = ind_value[i];
        auto &indiv = individual[i];
        auto g = indiv.group;
        
        for(auto j = 0u; j < nie; j++){
          set_precalc_ind_eff(i,prec_ind_eff[j],ind_eff[j],ind_value,inv_cov_matrix);
        }
        
        double infectivity_change, dL_inf_ev;
        
        auto inf_fac = log(ind.infectivity[c]/inf_av);
        for(auto loop = 0; loop < ind_effect_inf_prop; loop++){
          for(auto j = 0u; j < nie; j++){
            const auto &prec_ie = prec_ind_eff[j];
            auto ie = ind_eff[j];
            
            auto ind_eff_store = ind.ind_effect[ie];
            auto ind_eff_prop = normal_sample(prec_ie.mean,prec_ie.sd);
              
            auto inf_fac_store = inf_fac;
            inf_fac += ind_eff_prop - ind_eff_store;
              
            if(g != UNSET){
              infectivity_change = inf_av*(exp(inf_fac) - exp(inf_fac_store));
              dL_inf_ev = Li_change_inf(infectivity_change,prec,I_profile[g]);
            }
            else dL_inf_ev = 0;
            
            auto al = exp(dL_inf_ev);
            
            ind_effect_jump[ie].ntr++;
            if(ran() < al){
              ind_effect_jump[ie].nac++;
              
              ind.ind_effect[ie] = ind_eff_prop;
              
              auto cv = ind_effect[ie].covar_ref;
              L_ind_effect[cv] += Li_change_ind_eff(ind_eff_store,ind_eff_prop,prec_ie);
              if(g != UNSET){
                L_inf_events[g] += dL_inf_ev;
                update_I_profile(infectivity_change,prec,I_profile[g]);
              }
            }
            else{
              inf_fac = inf_fac_store;
            }
          }
        }
        ind.infectivity[c] = inf_av*exp(inf_fac);
      }
    }
  }
	
	// This changes infectivity individual effects from multiple compartments
	for(auto ie = 0u; ie < nind_effect; ie++){
		const auto &ind_eff = ind_effect[ie];
		const auto &inf_comp = ind_eff.infectivity_comp;
		auto nc = inf_comp.size();
		if(nc > 0 && ind_eff.single_ind_effect_proposal == true){
			PrecalcIndEff prec_ind_eff;
			
			vector <double> inf_fac(ncomp), inf_fac_store(ncomp), infectivity_change(ncomp);
			
			for(auto i = 0u; i < N; i++){
        auto &prec = precalc[i];
        auto &ind = ind_value[i];
        const auto &indiv = individual[i];
        auto g = indiv.group;
				
        set_precalc_ind_eff(i,prec_ind_eff,ie,ind_value,inv_cov_matrix);
				
				for(auto c : inf_comp) inf_fac[c] = log(ind.infectivity[c]/comp[c].infectivity);
					
				for(auto loop = 0; loop < ind_effect_inf_prop2; loop++){
					auto ind_eff_store = ind.ind_effect[ie];
					auto ind_eff_prop = normal_sample(prec_ind_eff.mean,prec_ind_eff.sd);
					 
					auto change = ind_eff_prop - ind_eff_store;
					for(auto c : inf_comp){
						inf_fac_store[c] = inf_fac[c];
						inf_fac[c] += change;
					}				
						
					auto dL_inf_ev = 0.0;
					if(g != UNSET){
						for(auto c : inf_comp){
							infectivity_change[c] = comp[c].infectivity*(exp(inf_fac[c]) - exp(inf_fac_store[c]));
							dL_inf_ev += Li_change_inf(infectivity_change[c],prec[c],I_profile[g]);
						}
					}
            
					if(false){ // Checks for correct sampling from individual effect distribution
						auto probfi = normal_probability(ind_eff_store,prec_ind_eff.mean,prec_ind_eff.sd);
						auto probif = normal_probability(ind_eff_prop,prec_ind_eff.mean,prec_ind_eff.sd);
						cout << probfi - probif + Li_change_ind_eff(ind_eff_store,ind_eff_prop,prec_ind_eff) << "zero" << endl;
					}
					
					auto al = exp(dL_inf_ev);
					
					ind_effect_jump[ie].ntr++;
					if(ran() < al){
						ind_effect_jump[ie].nac++;
						
						ind.ind_effect[ie] = ind_eff_prop;
						
						auto cv = ind_effect[ie].covar_ref;
						L_ind_effect[cv] += Li_change_ind_eff(ind_eff_store,ind_eff_prop,prec_ind_eff);
						if(g != UNSET){
							L_inf_events[g] += dL_inf_ev;
							for(auto c : inf_comp) update_I_profile(infectivity_change[c],prec[c],I_profile[g]);
						}
					}
					else{
						for(auto c : inf_comp) inf_fac[c] = inf_fac_store[c];
					}
				}
				
				for(auto c : inf_comp) ind.infectivity[c] = comp[c].infectivity*exp(inf_fac[c]);
			}
		}
	}
	
	/// Sets the change in infectivity down a transition
  for(auto &ind : ind_value){
    for(auto tr = 0u; tr < ntrans; tr++){
      ind.trans_infectivity_change[tr] = ind.infectivity[trans[tr].to] - ind.infectivity[trans[tr].from]; 
    }
  } 
  timer[TIME_INF_IND_EFFECT].stop();  
}


/// This pre-calculates quantities so likelihood can be calculated from infectivity
vector < vector <PrecalcLiInf> > Model::set_precalc_likelihood_inf(vector < vector <double> > &I_profile, const vector <IndValue> &ind_value, const vector <double> &param_value) const
{
  vector < vector <PrecalcLiInf> > prec;
  prec.resize(N);
  for(auto i = 0u; i < N; i++){
    prec[i].resize(ncomp);
    for(auto c = 0u; c < ncomp; c++) prec[i][c].beta_fac = 0;
  }
  
  vector <int> state(N);
  for(auto &sta : state) sta = 0;
  
  I_profile.resize(ngroup);
  for(const auto &gr : group){
    auto &Iprof = I_profile[gr.index];
    
    auto I = 0.0;                                            // The total infectivity
    auto S = 0.0;                                            // The total susceptibility
    
    auto ind_list = gr.ind_ref;                              // Makes a list of group individuals
    
    vector <Event> event;
    get_event_sequence(event,S,gr,ind_value);                // Gererates a sorted list of events for the entire group
      
    auto tmin = gr.inference_range.tmin;
    auto tmax = gr.inference_range.tmax;

    auto beta = param_value[beta_param];
    if(inf_model == FREQ_DEP) beta /= gr.nind;               // If frequency dependent divided by group size
    if(group_effect.on == true){                             // Adds in the group effect
      beta *= exp(param_value[group_effect.param[gr.index]]);
    }
    
    auto t = tmin;                                           // Sets a time variable
    for(const auto &ev : event){
      const auto &ind = ind_value[ev.ind];
      auto tev = ev.time;
      auto tr = ev.trans;
      
      auto val = -beta*S*(tev-t);
      for(auto i : ind_list) prec[i][state[i]].beta_fac += val;   
      t = tev;  
      
      auto i = ev.ind;
      if(tr == 0){
        if(tev > tmin){
          auto num = Iprof.size();
          for(auto i : ind_list){
            auto c = state[i];
            if(comp[c].infectivity != 0){
              prec[i][c].inf_log_terms.push_back(num);
            }
          }     
          Iprof.push_back(I+EXT_FOI);         
        }
        
        S -= ind.susceptibility;
      }
      
      I += ind.trans_infectivity_change[tr];
      state[i] = trans[tr].to;
    }
  
    auto val = -beta*S*(tmax-t);
    for(auto i : ind_list) prec[i][state[i]].beta_fac += val;   
  }
  
  return prec;
}


/// Calculates the change in event likelihood as a result of infectivity changing
double Model::Li_change_inf(const double infectivity_change, const PrecalcLiInf &prec, const vector <double> &I_profile) const
{
  auto dL = prec.beta_fac*infectivity_change;
  for(auto e : prec.inf_log_terms){
    dL += log(I_profile[e]+infectivity_change) - log(I_profile[e]);
  }
  
  return dL;
}


/// Updates I_profile (this stores how infectivity changes for all the infection events in the system)
void Model::update_I_profile(const double infectivity_change, const PrecalcLiInf &prec, vector <double> &I_profile) const
{
  for(auto e : prec.inf_log_terms) I_profile[e] += infectivity_change;
}
