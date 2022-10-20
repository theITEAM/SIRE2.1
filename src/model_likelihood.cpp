/// This gives model functions related to the prior

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "model.hpp"
#include "utils.hpp"
#include "timers.hpp"
#include "matrix.hpp"

/// Calculates the likelihood of the data given the model
double Model::likelihood(const vector <double> &param_value, vector <IndValue> &ind_value, const vector <GroupEvent> &gr_ev) const
{	
	auto L = 0.0;
	for(auto g = 0; g < ngroup; g++){
		L += likelihood_group(g,param_value,ind_value,gr_ev);
	}
	
	return L;
}


/// Calculate the likelihood for a given group zz1
double Model::likelihood_group(unsigned int g, const vector <double> &param_value, vector <IndValue> &ind_value, const vector <GroupEvent> &gr_ev) const
{
	set_individual_quantities_group(g,ind_value, param_value); 
	
	auto L = 0.0;
	
	const auto &gr = group[g];
	const auto &ge = gr_ev[g];
	auto beta = param_value[beta_param];
	if (inf_model == FREQ_DEP) beta /= gr.nind;        // If frequency dependent divided by group size
	if (group_effect.on == true)                       // Adds in the group effect
		beta *= exp(param_value[group_effect.param[gr.index]]);
	
	vector < vector <double> > comp_prob;
	comp_prob.resize(gr.nind);
	for(auto k = 0; k < gr.nind; k++){
		const auto &ind = individual[gr.ind_ref[k]];
		auto c_init = ind.initial_comp; 
		
		comp_prob[k].resize(ncomp);
		for(auto c = 0; c < ncomp; c++){
			if(c == c_init) comp_prob[k][c] = 1;
			else comp_prob[k][c] = 0;
		}
	}

	auto T = ge.T;

	auto dtsum = 0.0;

	for(auto ti = 0; ti < T; ti++){
		auto tmi = ti*dt, tma = (ti+1)*dt;
		
		// Works out total infectivity
		auto I = 0.0;
		for(auto k = 0; k < gr.nind; k++){
			const auto inf = ind_value[gr.ind_ref[k]].inf_single;
			for(auto c = 0; c < ncomp; c++) I += comp_prob[k][c]*comp[c].infectivity*inf;	
		}
		if(ti == 0) cout << I << "I\n";
		//cout << g << " " << ti << " " << I << " " << gr.nind << " " << tmi << " " << tma << " ti\n";
		
		for(auto k = 0; k < gr.nind; k++){
			auto &cp = comp_prob[k];
		
			// Accounts for compartmental model
			auto i = gr.ind_ref[k];
			const auto &ind = individual[i];
			const auto &indv = ind_value[i];
			
			for(const auto &sec : ge.ind_ev[k][ti].trans_section){
				auto dtt = sec.dt;
				dtsum += dtt;
				
				auto mmax = sec.trlist.size();
				for(auto m = 0; m < mmax; m++){
					const auto j = sec.trlist[m];
					const auto &tr = trans[j];
				
					auto v = cp[tr.from];
					if(v != 0){
						double prob;
						if(tr.type == INFECTION) prob = beta*I*indv.susceptibility*dtt;
						else prob = (1.0/indv.trans_mean[j])*dtt;
							
						//if(prob > 1) prob = 1;	qq
					
						v *= prob;
					}

					if(m == mmax-1){
						L -= v;
					}
					else{
						cp[tr.from] -= v;
						cp[tr.to] += v;
					}
				}
			
				if(sec.trans_end == true){
					const auto j = sec.trlist[mmax-1];
					const auto &tr = trans[j];
					
					auto v = cp[tr.from];
					if(v != 0){
						double prob;
						if(tr.type == INFECTION) prob = beta*I*indv.susceptibility;
						else prob = 1.0/indv.trans_mean[j];
	
						v *= prob;
					}
					
					if(v == 0) emsg("Zero probability");
					
					L += log(v + TINY);
					
					auto c_to = tr.to;
					for(auto c = 0; c < ncomp; c++){
						if(c == c_to) cp[c] = 1; else cp[c] = 0;
					}
				}
			}
		}
		//cout << ti << " " << L << "L\n";
	}
	
	/// Checks no individuals in an infectious state
	if(false){
		auto still_sus = false;
		for(auto k = 0; k < gr.nind; k++){
			if(comp_prob[k][0] != 0) still_sus = true;
		}
		
		if(still_sus == true){
			for(auto k = 0; k < gr.nind; k++){
				for(auto c = 0; c < ncomp; c++){
					if(comp_prob[k][c]*comp[c].infectivity != 0) emsg("Still infecious individuals");
				}
			}
		}
		
		auto dtcompare = gr.nind*dt*T;
		if(dtsum < dtcompare-SMALL || dtsum > dtcompare+SMALL){
			cout << g << " " << dtsum << " " << dtcompare << " Times\n";
			emsg("dt not agree");
		}	
	}
	
	return L;
}


/// Initialises group events (for fast evalulation of the likelihood)
void Model::initialise_group_events(vector <GroupEvent> &gr_ev) const
{
	for(auto g = 0; g < ngroup; g++){
		const auto &gr = group[g];
		
		auto tmin = 0.0, tmax = 0.0; 
		for(auto k = 0; k < gr.nind; k++){
			const auto &ind = individual[gr.ind_ref[k]];
			for(auto j = 0; j < ntrans; j++){
				auto t = ind.trans_obs_time[j];
				if(t != UNSET && t != LARGE && t > tmax) tmax = t;
			}
		}

		auto T = (unsigned int)(1+(tmax-tmin)/dt);
				
		GroupEvent ge;
		ge.T = T;
		
		ge.ind_ev.resize(gr.nind);
		for(auto k = 0; k < gr.nind; k++){
			ge.ind_ev[k].resize(T);
			
			const auto &ind = individual[gr.ind_ref[k]];		
			
			auto stat = ind.initial_comp;
			for(auto ti = 0; ti < T; ti++){
				auto tmi = ti*dt, tma = (ti+1)*dt;
			
				TransSection sec;
				sec.dt = dt;
				sec.trlist = get_trlist(stat);
				sec.trans_end = false;
				ge.ind_ev[k][ti].trans_section.push_back(sec);
				
				for(auto j = 0; j < ntrans; j++){
					auto t = ind.trans_obs_time[j];
					if(t != UNSET && t != LARGE && t > 0 && t >= tmi && t < tma){  
						auto &sec_last = ge.ind_ev[k][ti].trans_section[ge.ind_ev[k][ti].trans_section.size()-1];
					
						if(sec_last.trlist[sec_last.trlist.size()-1] != j) emsg("Problem");
						
						stat = trans[j].to;
					
						auto dt_change = tma-t;
					
						sec_last.trans_end = true;
						sec_last.dt -= dt_change;
						
						TransSection sec;
						sec.dt = dt_change;
						sec.trlist = get_trlist(stat);
						sec.trans_end = false;
						ge.ind_ev[k][ti].trans_section.push_back(sec);
					}
				}
			}
		}
		
		if(false){
			for(auto k = 0; k < gr.nind; k++){
				const auto &ind = individual[gr.ind_ref[k]];		
				
				for(auto ti = 0; ti < T; ti++){
					cout << ind.id << " " << k << " " << ti*dt << ": \n";
					for(auto sec : ge.ind_ev[k][ti].trans_section){
						cout << sec.dt << " ";
						for(auto j : sec.trlist) cout << j << ",";
						if(sec.trans_end == false) cout << " no event\n";
						else cout << "event\n";
					}
				}
			}
			emsg("Check");
		}
		
		gr_ev.push_back(ge);
	}
}


/// Based on a given compartment gives the transitions which are not known
vector <unsigned int> Model::get_trlist(unsigned int c) const
{
	vector <unsigned int> trlist;   // Makes a list of all transitions which could be passed down consistent with the data
	auto j = c;
	while(j < ntrans){
		trlist.push_back(j);
		if(trans[j].data_column != "") break;
		j++;
	}
	
	return trlist;
}


/// Calculates the posterior, gradient and hessian matrix zz2
double Model::posterior_grad_H(vector <double> &grad, vector < vector <double> > &H, const vector <double> &param_value, vector <IndValue> &ind_value, const vector <GroupEvent> &gr_ev, const vector <unsigned int> &param_var_ref, const vector < vector <unsigned int> > &ind_effect_var_ref, unsigned int nvar, bool calc_grad, bool calc_H) const
{
	set_individual_quantities(ind_value, param_value); 
	
	if(calc_grad == true){
		grad.resize(nvar);
		for(auto v = 0; v < nvar; v++) grad[v] = 0;
	}
	
	if(calc_H == true){
		H.resize(nvar);
		for(auto v = 0; v < nvar; v++){
			H[v].resize(nvar);
			for(auto vv = 0; vv < nvar; vv++) H[v][vv] = 0;
		}
	}
	
	auto L = 0.0;

	for(auto g = 0; g < ngroup; g++){
		timer[TIME_GRAD_H1].start();
		
		const auto &gr = group[g];
		const auto &ge = gr_ev[g];
		auto beta = param_value[beta_param];
		auto beta_fac = 1.0;
		if (inf_model == FREQ_DEP) beta_fac /= gr.nind;        // If frequency dependent divided by group size
		if (group_effect.on == true)                       // Adds in the group effect
			beta_fac *= exp(param_value[group_effect.param[gr.index]]);
		
		vector < vector <double> > comp_prob;
		comp_prob.resize(gr.nind);
		for(auto k = 0; k < gr.nind; k++){
			const auto &ind = individual[gr.ind_ref[k]];
			auto c_init = ind.initial_comp; 
			
			comp_prob[k].resize(ncomp);
			for(auto c = 0; c < ncomp; c++){
				if(c == c_init) comp_prob[k][c] = 1;
				else comp_prob[k][c] = 0;
			}
		}

		auto N = gr.nind;

		vector < vector <double> > dL_dr, d2L_dr2;
		dL_dr.resize(ntrans); d2L_dr2.resize(ntrans);
		for(auto tr = 1; tr < ntrans; tr++){
			dL_dr[tr].resize(N); d2L_dr2[tr].resize(N);
			for(auto k = 0; k < N; k++){ dL_dr[tr][k] = 0; d2L_dr2[tr][k] = 0;} 
		}
	
		vector <double> dL_dg(N,0);
		vector <double> dL_df(N,0);
		vector <double> d2L_dg2(N,0);
		vector < vector <double> > d2L_df2;
		vector < vector <double> > d2L_dg_df;
		
		if(calc_H == true){
			d2L_dg_df.resize(N); d2L_df2.resize(N);
	
			for(auto i = 0; i < N; i++){
				d2L_dg_df[i].resize(N); d2L_df2[i].resize(N); 
				for(auto j = 0; j < N; j++){ d2L_dg_df[i][j] = 0; d2L_df2[i][j] = 0;}
			}
		}
		
		auto T = ge.T;

		auto dtsum = 0.0;

		for(auto ti = 0; ti < T; ti++){
			auto tmi = ti*dt, tma = (ti+1)*dt;
			
			// Works out total infectivity
			
			vector <InfList> inf_list;  // Lists infected individuals
			
			auto I = 0.0;
			for(auto k = 0; k < N; k++){
				auto sum = 0.0;
				for(auto c = 0; c < ncomp; c++) sum += comp_prob[k][c]*comp[c].infectivity;
				
				if(sum > 0){
					sum *= ind_value[gr.ind_ref[k]].inf_single;
					I += sum; 
					InfList il; il.k = k; il.val = sum;
					inf_list.push_back(il);
				}
			}
		
			for(auto k = 0; k < N; k++){
				auto &cp = comp_prob[k];
			
				// Accounts for compartmental model
				auto i = gr.ind_ref[k];
				const auto &ind = individual[i];
				const auto &indv = ind_value[i];
				
				for(const auto &sec : ge.ind_ev[k][ti].trans_section){
					auto dtt = sec.dt;
					dtsum += dtt;
					
					auto mmax = sec.trlist.size();
					for(auto m = 0; m < mmax; m++){
						const auto j = sec.trlist[m];
						const auto &tr = trans[j];
				
						auto v = cp[tr.from];
						double prob;
						if(v != 0){
							if(tr.type == INFECTION) prob = beta*beta_fac*I*indv.susceptibility*dtt;
							else prob = (1.0/indv.trans_mean[j])*dtt;
								
							//if(prob > 1) prob = 1; qq
						
							v *= prob;
						}

						if(m == mmax-1){
							//if(calc_grad == true && prob < 1){ qq
							if(calc_grad == true){
								if(tr.type == INFECTION){
									dL_dg[k] -= v;		
									for(const auto &il : inf_list) dL_df[il.k] -= (v/I)*il.val;
		
									if(calc_H == true){
										d2L_dg2[k] -= v;
										for(const auto &il : inf_list){											
											auto kk = il.k;
											d2L_df2[kk][kk] -= (v/I)*il.val;
											d2L_dg_df[k][kk] -= (v/I)*il.val;
										}
									}
								}
								else{
									dL_dr[j][k] += v;
									if(calc_H == true) d2L_dr2[j][k] -= v;
								}
							}
							
							L -= v;
							if(std::isinf(L)) emsg("p1");
	
						}
						else{
							cp[tr.from] -= v;
							cp[tr.to] += v;
						}
					}
				
					if(sec.trans_end == true){
						const auto j = sec.trlist[mmax-1];
						const auto &tr = trans[j];
						
						auto v = cp[tr.from];
						if(v != 0){
							double prob;
							if(tr.type == INFECTION) prob = beta*beta_fac*I*indv.susceptibility;
							else prob = 1.0/indv.trans_mean[j];
						
							v *= prob;
						}
						
						L += log(v + TINY);
						if(std::isinf(L)){ 
							if(tr.type == INFECTION) cout << "inf\n";
							else cout << "notinf\n";
							
							for(auto i = 0; i < nindividual; i++){
								cout << individual[i].id << " " <<  ind_value[i].inf_single << " " << ind_value[i].ind_effect[0] <<  " sing\n";
								if(std::isinf(ind_value[i].inf_single )) emsg("do");
							}
							cout << v << " " << I << "v\n"; emsg("p2");
						}
	
						/*
						if(v == 0){
							cout << ind.id << " " << tmi << endl;
							if(tr.type == INFECTION) cout << "inf\n";
							else cout << "notinf\n";
							cout << cp[tr.from] << " " << indv.trans_mean[j] << "fr\n";
							emsg("Zero probability for transition");
						}
						*/
						
						if(calc_grad == true){
							if(tr.type == INFECTION){
								dL_dg[k]++;
								for(const auto &il : inf_list) dL_df[il.k] += il.val/I;
								
								if(calc_H == true){
									for(const auto &il : inf_list){	
										auto kk = il.k;
										d2L_df2[kk][kk] += il.val/I;
										for(const auto &il2 : inf_list){	
											d2L_df2[kk][il2.k] -= (il.val*il2.val)/(I*I);
										}
									}
								}
							}
							else{
								dL_dr[j][k]--;						
							}
						}
						
						auto c_to = tr.to;
						for(auto c = 0; c < ncomp; c++){
							if(c == c_to) cp[c] = 1; else cp[c] = 0;
						}
					}
				}
			}
		}
		timer[TIME_GRAD_H1].stop();
		
		timer[TIME_GRAD_H2].start();
		
		if(false){
			for(auto i = 0; i < N; i++){
				for(auto j = 0; j < N; j++){
					cout << d2L_df2[i][j] << " ";
				}
				cout << "d2L_df2" << endl;
			}
			
			for(auto i = 0; i < N; i++){
				for(auto j = 0; j < N; j++){
					cout << d2L_dg_df[i][j] << " ";
				}
				cout << "d2L_dg_df" << endl;
			}
			emsg("P");
		}
		
		if(calc_grad == true){	
			auto v_beta = param_var_ref[beta_param];
			
			auto v_group = UNSET;
			if (group_effect.on == true) v_group = param_var_ref[group_effect.param[g]];
					
			/// susceptibility terms
			for(auto k = 0; k < N; k++){
				auto gradient = dL_dg[k];
				
				if(v_beta != UNSET){
					grad[v_beta] += gradient/beta;
					if(calc_H == true) H[v_beta][v_beta] -= gradient/(beta*beta);
				}
				
				if(v_group != UNSET) grad[v_group] += gradient;
				
				auto i = gr.ind_ref[k];
				for(auto ie : trans[0].ind_effect){		
					auto v_ie = ind_effect_var_ref[i][ie];
					if(v_ie != UNSET) grad[v_ie] += gradient;
				}
				
				for(auto fe : trans[0].fixed_effect){
					auto v_fix = param_var_ref[fixed_effect[fe].param];
					if(v_fix != UNSET) grad[v_fix] += individual[i].fixed_effect_X[fe]*gradient;
				}	
			}
			
			if(calc_H == true){
				for(auto k = 0; k < N; k++){
					auto i = gr.ind_ref[k];
				
					auto hess = d2L_dg2[k];
					
					vector <DepList> dep_list;
					
					if(v_beta != UNSET){ DepList dl; dl.v = v_beta; dl.val = 1.0/beta; dep_list.push_back(dl);}
					
					if(v_group != UNSET){ DepList dl; dl.v = v_group; dl.val = 1.0; dep_list.push_back(dl);}
						
					for(auto ie : trans[0].ind_effect){		
						auto v_ie = ind_effect_var_ref[i][ie];
						if(v_ie != UNSET){ DepList dl; dl.v = v_ie; dl.val = 1.0; dep_list.push_back(dl);}
					}
					
					for(auto fe : trans[0].fixed_effect){
						auto v_fix = param_var_ref[fixed_effect[fe].param];
						if(v_fix != UNSET){ DepList dl; dl.v = v_fix; dl.val = individual[i].fixed_effect_X[fe]; dep_list.push_back(dl);}
					}	
				
					for(const auto &dl : dep_list){
						for(const auto &dl2 : dep_list){
							H[dl.v][dl2.v] += dl.val*dl2.val*hess;
						}
					}
				}
			}
			
			// Infectivity effects
			for(auto k = 0; k < N; k++){
				auto gradient = dL_df[k];
				
				auto i = gr.ind_ref[k];
				for(auto ie : infectivity.ind_effect){		
					auto v_ie = ind_effect_var_ref[i][ie];
					if(v_ie != UNSET) grad[v_ie] += gradient;
				}
				
				for(auto fe : infectivity.fixed_effect){
					auto v_fix = param_var_ref[fixed_effect[fe].param];
					if(v_fix != UNSET) grad[v_fix] += individual[i].fixed_effect_X[fe]*gradient;
				}	
			}
			
			if(calc_H == true){
				for(auto k = 0; k < N; k++){
					auto i = gr.ind_ref[k];
					for(auto kk = 0; kk < N; kk++){
						auto ii = gr.ind_ref[kk];
					
						auto hess = d2L_df2[k][kk];
						if(hess != 0){
							vector <DepList> dep_list, dep_list2;
								
							for(auto ie : infectivity.ind_effect){		
								auto v_ie = ind_effect_var_ref[i][ie];
								if(v_ie != UNSET){ DepList dl; dl.v = v_ie; dl.val = 1.0; dep_list.push_back(dl);}
								
								v_ie = ind_effect_var_ref[ii][ie];
								if(v_ie != UNSET){ DepList dl; dl.v = v_ie; dl.val = 1.0; dep_list2.push_back(dl);}
							}
							
							for(auto fe : infectivity.fixed_effect){
								auto v_fix = param_var_ref[fixed_effect[fe].param];
								if(v_fix != UNSET){ 
									DepList dl; dl.v = v_fix;
									dl.val = individual[i].fixed_effect_X[fe]; 
									dep_list.push_back(dl);
									dl.val = individual[ii].fixed_effect_X[fe]; 
									dep_list2.push_back(dl);
								}
							}	
						
							for(const auto &dl : dep_list){
								for(const auto &dl2 : dep_list2){
									H[dl.v][dl2.v] += dl.val*dl2.val*hess;
								}
							}
						}
					}
				}
			}
			
			// Infectivity / susceptibility effects
			
			if(calc_H == true){
				for(auto k = 0; k < N; k++){
					auto i = gr.ind_ref[k];
					for(auto kk = 0; kk < N; kk++){
						auto ii = gr.ind_ref[kk];
					
						auto hess = d2L_dg_df[k][kk];
						if(hess != 0){
							vector <DepList> dep_list, dep_list2;
							
							if(v_beta != UNSET){ DepList dl; dl.v = v_beta; dl.val = 1.0/beta; dep_list.push_back(dl);}
					
							if(v_group != UNSET){ DepList dl; dl.v = v_group; dl.val = 1.0; dep_list.push_back(dl);}
					
							for(auto ie : trans[0].ind_effect){		
								auto v_ie = ind_effect_var_ref[i][ie];
								if(v_ie != UNSET){ DepList dl; dl.v = v_ie; dl.val = 1.0; dep_list.push_back(dl);}
							}
							
							for(auto fe : trans[0].fixed_effect){
								auto v_fix = param_var_ref[fixed_effect[fe].param];
								if(v_fix != UNSET){ DepList dl; dl.v = v_fix; dl.val = individual[i].fixed_effect_X[fe]; dep_list.push_back(dl);}
							}	
							
							for(auto ie : infectivity.ind_effect){		
								auto v_ie = ind_effect_var_ref[ii][ie];
								if(v_ie != UNSET){ DepList dl; dl.v = v_ie; dl.val = 1.0; dep_list2.push_back(dl);}
							}
							
							for(auto fe : infectivity.fixed_effect){
								auto v_fix = param_var_ref[fixed_effect[fe].param];
								if(v_fix != UNSET){ 
									DepList dl; dl.v = v_fix; dl.val = individual[ii].fixed_effect_X[fe]; 
									dep_list2.push_back(dl);
								}
							}	
						
							for(const auto &dl : dep_list){
								auto val = dl.val*hess;
								for(const auto &dl2 : dep_list2){
									auto val2 = val*dl2.val;
									H[dl.v][dl2.v] += val2;
									H[dl2.v][dl.v] += val2;
								}
							}
						}
					}
				}
			}
			
			// Transition terms
			for(auto tr = 1; tr < ntrans; tr++){
				const auto &tra = trans[tr];
				for(auto k = 0; k < N; k++){
					auto gradient = dL_dr[tr][k];
					
					auto th_mean = tra.mean_param;
					auto v_mean = param_var_ref[th_mean];
					if(v_mean != UNSET){
						auto mean = param_value[th_mean];
						grad[v_mean] += gradient/mean;
						if(calc_H == true) H[v_mean][v_mean] -= gradient/(mean*mean);
					}
				
					auto i = gr.ind_ref[k];
					for(auto ie : tra.ind_effect){		
						auto v_ie = ind_effect_var_ref[i][ie];
						if(v_ie != UNSET) grad[v_ie] += gradient;
					}
					
					for(auto fe : tra.fixed_effect){
						auto v_fix = param_var_ref[fixed_effect[fe].param];
						if(v_fix != UNSET) grad[v_fix] += individual[i].fixed_effect_X[fe]*gradient;
					}	
				}
				
				if(calc_H == true){
					for(auto k = 0; k < N; k++){
						auto i = gr.ind_ref[k];
					
						auto hess = d2L_dr2[tr][k];
						
						vector <DepList> dep_list;
						
						auto th_mean = trans[tr].mean_param;
						auto v_mean = param_var_ref[th_mean];
						if(v_mean != UNSET){
							auto mean = param_value[th_mean];
							DepList dl; dl.v = v_mean; dl.val = 1.0/mean; dep_list.push_back(dl);
						}
			
						for(auto ie : tra.ind_effect){		
							auto v_ie = ind_effect_var_ref[i][ie];
							if(v_ie != UNSET){ DepList dl; dl.v = v_ie; dl.val = 1.0; dep_list.push_back(dl);}
						}
						
						for(auto fe : tra.fixed_effect){
							auto v_fix = param_var_ref[fixed_effect[fe].param];
							if(v_fix != UNSET){ DepList dl; dl.v = v_fix; dl.val = individual[i].fixed_effect_X[fe]; dep_list.push_back(dl);}
						}	
					
						for(const auto &dl : dep_list){
							for(const auto &dl2 : dep_list){
								H[dl.v][dl2.v] += dl.val*dl2.val*hess;
							}
						}
					}
				}
			}
		}
		timer[TIME_GRAD_H2].stop();
	}
	
	timer[TIME_GRAD_H3].start();
	
	if(calc_grad == true){	
		auto inv_cov_matrix = calculate_inv_cov_matrix(param_value);
		
		for(auto i = 0; i < N; i++){
			for(auto ie = 0; ie < nind_effect; ie++){
				auto v = ind_effect_var_ref[i][ie];
				grad[v] += ind_effect_gradient(i,ie,ind_value,inv_cov_matrix);
			}
		}
	
		if(calc_H == true){
			for (auto cv = 0; cv < ncovariance; cv++) {
				const auto &cov = covariance[cv];
				
				const auto &M = matrix[cov.matrix];
				for(auto i1 = 0; i1 < cov.E; i1++){
					auto ie1 = cov.ind_effect_ref[i1];
					for(auto i2 = 0; i2 < cov.E; i2++){
						auto ie2 = cov.ind_effect_ref[i2];

						auto variance = inv_cov_matrix[cv].M[i1][i2];
						for(auto i = 0; i < N; i++){
							auto v1 = ind_effect_var_ref[i][ie1];
							auto v2 = ind_effect_var_ref[i][ie2];
							
							H[v1][v2] -= M.Ainvdiag[i]*variance;
							for(const auto el : M.Ainvlist[i]){
								auto j = el.i;
								auto v2 = ind_effect_var_ref[j][ie2];
								H[v1][v2] -= el.val*variance;
							}
						}
					}
				}
			}
		}
	}
	
	auto L_ind_eff = calculate_L_ind_effect(ind_value, param_value);
	for(auto val : L_ind_eff) L += val;
	if(std::isinf(L)) emsg("p2");
	
	L += calculate_prior(param_value);
	if(std::isinf(L)) emsg("p3");
	
	for(auto th = 0; th < nparam; th++){   // Adds any gradients in the prior
		auto v = param_var_ref[th];
		if(v != UNSET){ 
			const auto &par = param[th];
			if(par.prior_type == NORMAL_FROM_SD_PRIOR){
				auto sd = param_value[par.prior_sd_param];
				
				if(calc_grad == true) grad[v] -= param_value[th]/(sd*sd);
				if(calc_H == true) H[v][v] -= 1.0/(sd*sd);
			}
		}
	}
	
	timer[TIME_GRAD_H3].stop();
	
	return L;
}


/// This constructs an estimate for Binv used in the 
vector < vector <double> > Model::construct_initial_Binv(const vector <double> &param_value, const vector <unsigned int> &param_var_ref, const vector < vector <unsigned int> > &ind_effect_var_ref, unsigned int nvar) const
{
	vector < vector <double> > Binv;
	
	Binv.resize(nvar);
	for(auto v = 0; v < nvar; v++){
		Binv[v].resize(nvar);
		for(auto vv = 0; vv < nvar; vv++){
			Binv[v][vv] = 0;
		}
	}
	
		for (auto cv = 0; cv < ncovariance; cv++) {
		const auto &cov = covariance[cv];
		auto cov_matrix = set_covariance_matrix(cov, param_value);
		
		const auto &M = matrix[cov.matrix];
		
		for(auto i = 0; i < N; i++){
			for(auto j = 0; j < N; j++){	
				auto A = M.A[i][j];
				if(A != 0){
					for(auto i1 = 0; i1 < cov.E; i1++){
						auto ie1 = cov.ind_effect_ref[i1];
						auto v1 = ind_effect_var_ref[i][ie1];
			
						for(auto i2 = 0; i2 < cov.E; i2++){
							auto ie2 = cov.ind_effect_ref[i2];
							auto v2 = ind_effect_var_ref[j][ie2];
				
							auto val = A*cov_matrix[i1][i2];
							Binv[v1][v2] += val;
						}
					}
				}
			}
		}
	}
	
	if(false){  // This code is used to check that Binv has been created correctly
		vector < vector <double> > H;
		H.resize(nvar);
		for(auto v = 0; v < nvar; v++){
			H[v].resize(nvar);
			for(auto vv = 0; vv < nvar; vv++){
				H[v][vv] = 0;
			}
		}
		
		auto inv_cov_matrix = calculate_inv_cov_matrix(param_value);
			
		for (auto cv = 0; cv < ncovariance; cv++) {
			const auto &cov = covariance[cv];
				
			const auto &M = matrix[cov.matrix];
			for(auto i1 = 0; i1 < cov.E; i1++){
				auto ie1 = cov.ind_effect_ref[i1];
				for(auto i2 = 0; i2 < cov.E; i2++){
					auto ie2 = cov.ind_effect_ref[i2];

					auto variance = inv_cov_matrix[cv].M[i1][i2];
					for(auto i = 0; i < N; i++){
						auto v1 = ind_effect_var_ref[i][ie1];
						auto v2 = ind_effect_var_ref[i][ie2];
						
						H[v1][v2] += M.Ainvdiag[i]*variance;
						for(const auto el : M.Ainvlist[i]){
							auto j = el.i;
							auto v2 = ind_effect_var_ref[j][ie2];
							H[v1][v2] += el.val*variance;
						}
					}
				}
			}
		}
		
		auto Hinv = invert_matrix(H);
	
		for(auto v = 0; v < nvar; v++){
			for(auto vv = 0; vv < nvar; vv++){
				auto d = Binv[v][vv] - Hinv[v][vv];
				if(d*d > SMALL) emsg("Problem with Binv");
			}
		}
	}
	
	return Binv;
}
