/// This provides routines that check the algorithm is performing correctly

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "mcmc.hh"
#include "utils.hh"
#include "check.hh"
#include "timers.hh"

/// Checks the quantities on the chain are correctly set
void MCMC::check_chain()
{
  timer[TIME_CHECK].start();
  
  // Checks that individual quantities are correctly set
  if(model.nind_effect > 0){
	  auto ind_value_check = ind_value;
    model.set_individual_quantities(ind_value,param_value);
    for(auto i = 0u; i < model.N; i++){
      auto &ind = ind_value[i];
      const auto &ind_check = ind_value_check[i];
      
      if(different(ind.susceptibility,ind_check.susceptibility) == true){
        emsg("Susceptibility problem");
      }
      ind.susceptibility = ind_check.susceptibility;
      
      for(auto c = 0u; c < model.ncomp; c++){
        if(different(ind.infectivity[c],ind_check.infectivity[c]) == true){
          emsg("Infectivity problem");
        }
        ind.infectivity[c] = ind_check.infectivity[c];
      }
      
      for(auto tr = 0u; tr < model.ntrans; tr++){  
        if(different(ind.trans_mean[tr],ind_check.trans_mean[tr]) == true){
          emsg("Trans_mean problem");
        }
        ind.trans_mean[tr] = ind_check.trans_mean[tr];
        
        if(different(ind.trans_infectivity_change[tr],ind_check.trans_infectivity_change[tr]) == true){
          cout << tr << " " << i << " " << ind.trans_infectivity_change[tr] << " " << ind_check.trans_infectivity_change[tr]<< " ch\n";
          emsg("trans_infectivity_change problem");
        }
        ind.trans_infectivity_change[tr] = ind_check.trans_infectivity_change[tr];
      }
    }
  }
  
  // Checks the prior is correctly set
  auto prior_check = model.calculate_prior(param_value);
  if(different(prior,prior_check) == true) emsg("prior problem");
  prior = prior_check;

  // Checks the infection event likelihoods are correctly set
  auto L_inf_events_check = model.calculate_L_inf_events(ind_value,param_value);
  for(auto g = 0; g < model.ngroup; g++){
    if(different(L_inf_events[g],L_inf_events_check[g]) == true) emsg("L_inf_events problem");
  }
  L_inf_events = L_inf_events_check;
  
	 /// Checks the individual effect likelihoods are correctly set
  auto L_ind_effect_check = model.calculate_L_ind_effect(ind_value,param_value);
  for(auto c = 0; c < model.ncovariance; c++){
    if(different(L_ind_effect[c],L_ind_effect_check[c]) == true) emsg("L_ind_effect problem");
  }
  L_ind_effect = L_ind_effect_check;
  
  // Checks the transtion event likelihood is correctly set
  auto L_trans_events_check = model.calculate_L_trans_events(ind_value,param_value);
  if(different(L_trans_events,L_trans_events_check) == true) emsg("L_trans_events problem");
  L_trans_events = L_trans_events_check;

	 // Checks the diagnostic test likelihoods are correctly set
  auto L_diag_test_check = model.calculate_L_diag_test(ind_value,param_value);
  for(auto g = 0; g < model.ngroup; g++){
    if(different(L_diag_test[g],L_diag_test_check[g]) == true) emsg("L_diag_test problem");
  }
  L_diag_test = L_diag_test_check;
  

  // Checks that events sequences for individuals are consistent
  for(auto i = 0u; i < model.N; i++){
    const auto &ind = ind_value[i];
    const auto &indiv = model.individual[i];
    if(ind.infected == false){
      for(auto tr = 0u; tr < model.ntrans; tr++){
        if(ind.trans_time[tr] != UNSET) emsg("Trans time should be unset");
      }
    }
    else{
      for(auto tr = 1u; tr < model.ntrans; tr++){
        if(ind.trans_time[tr] < ind.trans_time[tr-1]) emsg("Trans time in wrong order");
      }
      
      for(auto tr = 1u; tr < model.ntrans; tr++){
        auto t = ind.trans_time[tr];
        const auto &trange = indiv.trans_time_range[tr];
        if(t < trange.tmin || t > trange.tmax){
          emsg("Transition not in time range");
        }
      }     
    }   
  }
  
  // Check for not a number
  if(std::isnan(prior) || std::isinf(prior)) emsg("Prior is not a number");

  for(auto c = 0; c < model.ncovariance; c++){
    if(std::isnan(L_ind_effect[c]) || std::isinf(L_ind_effect[c])) emsg("L_ind_effect is not a number");
  }
  
  for(auto g = 0; g < model.ngroup; g++){
    if(std::isnan(L_inf_events[g]) || std::isinf(L_inf_events[g])) emsg("L_inf_events is not a number");
  }
  
  if(std::isnan(L_trans_events) || std::isinf(L_trans_events)) emsg("L_trans_events is not a number");

	for(auto g = 0; g < model.ngroup; g++){
    if(std::isnan(L_diag_test[g]) || std::isinf(L_diag_test[g])) emsg("L_diag_test is not a number");
  }
	
  timer[TIME_CHECK].stop();
}


/// This performs general checks
void check(const Model &model)
{
  vector< vector <double> > A { { -2, -1,2 },
                                { 2, 1,4 }, {-3,3,-1}};
  
	vector< vector <double> > B { { -2, -1},
                                { 2, 3 }};
  auto C = model.tensor_product(A,B);
    
  cout << model.determinant(A) << " " << model.determinant(B) << " " << model.determinant(C) <<  " " << pow(model.determinant(A),2)*pow(model.determinant(B),3) << " det" << endl;
}


/// This checks that the individual effect sampler is working
void check_ind_effect_sample(const Model &model)
{
  cout << "Checking ind_effect_sampler" << endl;
  
  auto param_value = model.prior_sample();
  
  vector <IndValue> ind_value(model.N);
  
  const auto &cov = model.covariance[0];
  
  /*
  for(auto i =0u; i < 300; i++){
    for(auto j =0u; j < 300; j++){
      auto val = model.matrix[cov.matrix].A[i][j];
      if(val != 0){
        cout << val << "val" << endl;
        if(val > 0.24 && val < 0.26) cout << i << " " << j << "  comb\n";
      }
    }
  }
  */
  
  auto ind1 = 245u, ind2 = 245u;
  //auto ind1 = 298u, ind2 = 106u;
  //auto ind1 = 145u, ind2 = 321u;
  
  auto E = cov.E;
  
  auto cov_mat = model.set_covariance_matrix(cov,param_value);
  
  auto nav = 0.0;
  vector <double> av1(E), av2(E);
  vector < vector <double> > av12;
  av12.resize(E);
  for(auto j = 0u; j < E; j++){
    av1[j] = 0; av2[j] = 0;
    av12[j].resize(E); for(auto i = 0u; i < E; i++) av12[j][i] = 0;
  }
  
  const auto loopmax = 100;
  for(auto loop = 0u; loop < loopmax; loop++){
    cout << loop << " Sample" << endl;
    model.ind_effect_sample(ind_value,param_value);
    
    for(auto j = 0u; j < E; j++){
      av1[j] += ind_value[ind1].ind_effect[cov.ind_effect_ref[j]];
      av2[j] += ind_value[ind2].ind_effect[cov.ind_effect_ref[j]];
      for(auto i = 0u; i < E; i++){
        av12[j][i] += ind_value[ind1].ind_effect[cov.ind_effect_ref[j]]*
                      ind_value[ind2].ind_effect[cov.ind_effect_ref[i]];
      }
    }
    nav++;
  }
  
  cout << "Expected:" << endl;
  for(auto j = 0u; j < E; j++){
    for(auto i = 0u; i < E; i++){
      cout << model.matrix[cov.matrix].A[ind1][ind2]*cov_mat[j][i] << " ";
    }
    cout << endl;
  } 
  
  cout << endl;
  cout << "Observed:" << endl;
  for(auto j = 0u; j < E; j++){
    for(auto i = 0u; i < E; i++){
      cout << av12[j][i]/nav - (av1[j]/nav)*(av2[i]/nav) << " ";
    }
    cout << endl;
  }
  
  cout << endl;
  cout << "Observed mean:" << endl;
  for(auto j = 0u; j < E; j++){
    cout << av1[j]/nav << " " << av2[j]/nav << endl;
  } 
}


/// A temporary function for testing ideas
void MCMC::temp()
{
	// Checks reconstruction from pedigree is correct
	int N = 6;
	vector < vector <double> > Ainv;
	
	Ainv.resize(N);
	for(auto j = 0; j < N; j++){
		Ainv[j].resize(N);
		for(auto i = 0u; i < N; i++){
			if(i == j) Ainv[j][i] = 1;
			else Ainv[j][i] = 0;
		}
	}
	
	//vector <int> ind_list = {0,1,2,3};
	//vector <int> par1_list = {UNSET,UNSET,0,1};
	//vector <int> par2_list = {UNSET,UNSET,UNSET,UNSET};
	
	/*
	vector <int> ind_list = {0,1,2,3};
	vector <int> par1_list = {UNSET,UNSET,0,0};
	vector <int> par2_list = {UNSET,UNSET,1,1};
	*/
	/*
	vector <int> ind_list = {0,1,2,3,4};
	vector <int> par1_list = {UNSET,UNSET,UNSET,0,1};
	vector <int> par2_list = {UNSET,UNSET,UNSET,1,2};
	*/
	vector <int> ind_list = {0,1,2,3,4,5};
	vector <int> par1_list = {UNSET,UNSET,UNSET,0,1,3};
	vector <int> par2_list = {UNSET,UNSET,UNSET,1,2,4};
	
	
	for(auto k = 0u; k < ind_list.size(); k++){
		auto i = ind_list[k];
		auto par1 = par1_list[k];
		auto par2 = par2_list[k];
		
		if(par1 == UNSET && par2 == UNSET){  // Both parents unknown
		}
		else{
			if(par1 == UNSET || par2 == UNSET){ // One parent known
				auto p = par1; if(par1 == UNSET) p = par2;
				Ainv[p][p] += 1.0/3;
				Ainv[i][i] += 1.0/3;
				Ainv[p][i] -= 2.0/3;
				Ainv[i][p] -= 2.0/3;
			}
			else{                              // Both parents known
				Ainv[par1][par1] += 1.0/2;
				Ainv[par2][par2] += 1.0/2;		
				Ainv[i][i] += 1.0;
				
				Ainv[par1][par2] += 1.0/2;
				Ainv[par2][par1] += 1.0/2;
				Ainv[par1][i] -= 1.0;
				Ainv[i][par1] -= 1.0;
				Ainv[par2][i] -= 1.0;
				Ainv[i][par2] -= 1.0;
			}
		}	
	}
	
	auto A = model.invert_matrix(Ainv);
	
	for(auto j = 0; j < N; j++){
		for(auto i = 0u; i < N; i++){
			cout << Ainv[j][i] << " ";
		}
		cout << "   Ainv" << endl;
	}
	
	for(auto j = 0; j < N; j++){
		for(auto i = 0u; i < N; i++){
			cout << A[j][i] << " ";
		}
		cout << endl;
	}
	
	emsg("relationship");
	/*
  ofstream Lplot("Lplot.txt");
  auto sum = 0.0; for(auto val : L_inf_events) sum += val;
  
  for(param_value[0] = 0.0;  param_value[0] < 0.1; param_value[0]+= 0.001){
    auto L_inf_events_new = model.calculate_L_inf_events(ind_value,param_value);
    auto sum_new = 0.0; for(auto val : L_inf_events_new) sum_new += val;
    Lplot << param_value[0] << " " << exp(sum_new-sum) << endl;
    cout << param_value[0] << " " << exp(sum_new-sum) << endl;
  } 
	*/
}
