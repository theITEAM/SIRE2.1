// Functions used to check the emulator is working correctly

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "emulator.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include "matrix.hpp"

/// Estimates posterior from a parameter set
double Emulator::calculate_posterior(const vector <double> &value)
{
	update_param_from_value(value);
	
	auto L = model.likelihood(param_value,ind_value,gr_ev);
	//L = 0.0;
	auto L_ind_eff = model.calculate_L_ind_effect(ind_value, param_value);
	for(auto val : L_ind_eff) L += val;
	
	L += model.calculate_prior(param_value);

	return L;	
}


/// Checks gradients and Hessian are correctly calculatued
void Emulator::check_grad_H()
{
	vector <double> grad;
	vector < vector <double> > H;
	
	for(auto i = 0; i < model.nindividual; i++){
		for(auto ie = 0; ie < model.nind_effect; ie++){
			ind_value[i].ind_effect[ie] = normal_sample(0,1);
		}
	}
				
	auto post = model.posterior_grad_H(grad,H,param_value,ind_value,gr_ev,param_var_ref,ind_effect_var_ref,var.size(),true,true);
	
	auto value = get_value_from_param();
	
	if(false){
		auto post_ch = calculate_posterior(value);
		cout << post_ch << " " << post << " Check posterior\n";
		auto d = post_ch - post; if(d*d > TINY) emsg("Posterior problem");
	}
	
	if(true){
		vector <unsigned int> var_list;
		auto vmax = nvar; if(vmax > 120) vmax = 120;
		for(auto v = 0; v < vmax; v++) var_list.push_back(v);
		
		auto grad_ch = calculate_gradient_check(value,var_list);
		
		for(auto v = 0; v < vmax; v++){
			cout << v << "   check:" << grad_ch[v] << "    actual:"<< grad[v] << " Check Gradient\n";
			
			auto thresh = 0.001*grad[v]; if(thresh < 0) thresh = -thresh;
			auto d = grad_ch[v] - grad[v];
			if(d > thresh || d < -thresh){
				cout << v << " " << grad_ch[v] << " " <<  grad[v] << "grad wrong\n";
				//cout << "IND EFFECT " << var[v].ind << " " << var[v].ref << " ind\n";
				emsg("wrong grad");	
			}			
		}
	}
	
	if(true){
		auto var_max = nvar; 
		var_max = 200; 
		if(var_max > nvar) var_max = nvar;
		
		auto H_ch = calculate_hessian_check(value,var_max);
	
		if(mpi.core == 0){
			ofstream Hcomp("Hcomp.txt");
			for(auto v = 0; v < var_max; v++){
				for(auto vv = 0; vv < var_max; vv++){
					auto thresh = 0.001*H[v][vv]; if(thresh < 0) thresh = -thresh;
					auto thresh2 = 0.001*H_ch[v][vv]; if(thresh2 < 0) thresh2 = -thresh2;
					if(thresh2 > thresh) thresh = thresh2;
					if(thresh < 0.1) thresh = 0.1;
					
					
					//cout << v << " " << vv << " CHeck:" <<  H_ch[v][vv] << "  Actual:" << H[v][vv] << " Check Hessian\n";
					auto thr = 0.001;
					if(H_ch[v][vv] < -thr || H_ch[v][vv] > thr || H[v][vv] < -thr || H[v][vv] > thr){
						Hcomp << v << " " << vv << " CHeck:" <<  H_ch[v][vv] << "  Actual:" << H[v][vv] << endl;
					}
					
					auto d = H_ch[v][vv] - H[v][vv];
					if(d > thresh || d < -thresh){
						//cout << "check:" << H_ch[v][vv] << " " << H[v][vv] << " compare\n";
						//emsg("Wrong Hessian");	
						cout << v << " " << vv << " CHeck:" <<  H_ch[v][vv] << "  Actual:" << H[v][vv] << " Wrong Hessian" << endl;
					}				
				}
			}
		}
	}
	
	/*
	auto vv = 103;
	for(auto v = 0; v < 200; v++){
		cout << v << " " <<  H[v][vv] << " " << H_ch[v][vv] << "check H\n";
	}
	*/
	
}


/// Calculates the Hessian and uses the fact the fact that likelihood independent between groups
vector < vector <double> > Emulator::calculate_hessian_check(const vector <double> &mean, unsigned int var_max)
{
	auto d = finite_different_grid_size(mean);          // Determines grid size to use
	for(auto &val: d) val *= 0.1;
	
	// Should check not overlapping boundaries...
	
	auto testfl = false;
	
	vector <MatrixElement> matrix_element;
	if(true){
		for(auto j = 0; j < var_max; j++){
			for(auto i = j; i < var_max; i++){
				MatrixElement el; el.j = j; el.i = i;
				matrix_element.push_back(el);
			}
		}
	}
	else{
		testfl = true;
		for(auto j = 103; j < 104; j++){ 
			for(auto i = 0; i < 200; i++){
				MatrixElement el; el.j = j; el.i = i;
				matrix_element.push_back(el);
			}
		}
	}
	
	auto element_per_core = (unsigned int)((matrix_element.size()+mpi.ncore-1)/mpi.ncore);
	
	vector <double> ur, dr, ul, dl;
	
	for(auto k = 0; k < element_per_core; k++){
		if(mpi.core == 0) cout << k << " / " << element_per_core << endl;
		
		auto l = mpi.core*element_per_core + k;
		if(l < matrix_element.size()){	
			auto &el = matrix_element[l];
	
			auto vari = el.i, varj = el.j;

			auto value = mean; value[vari] += d[vari]; value[varj] += d[varj];
			ur.push_back(calculate_posterior(value));
		
			value = mean; value[vari] += d[vari]; value[varj] -= d[varj];
			dr.push_back(calculate_posterior(value));
		
			value = mean; value[vari] -= d[vari]; value[varj] += d[varj];
			ul.push_back(calculate_posterior(value));
		
			value = mean; value[vari] -= d[vari]; value[varj] -= d[varj];
			dl.push_back(calculate_posterior(value));
		}
		else{
			ur.push_back(UNSET); dr.push_back(UNSET); ul.push_back(UNSET); dl.push_back(UNSET);
		}
	}

	auto ur_tot = mpi.gather(ur);
	auto dr_tot = mpi.gather(dr);
	auto ul_tot = mpi.gather(ul);
	auto dl_tot = mpi.gather(dl);
	
	vector < vector <double> > H;
	H.resize(nvar);
	for(auto j = 0; j < nvar; j++){
		H[j].resize(nvar); for(auto i = 0; i < nvar; i++) H[j][i] = UNSET;
	}
		
	if(mpi.core == 0){	
		for(auto k = 0; k < matrix_element.size(); k++){
			const auto &el = matrix_element[k];
			H[el.j][el.i] = (ur_tot[k] + dl_tot[k] - ul_tot[k] - dr_tot[k])/(4*d[el.i]*d[el.j]);
			H[el.i][el.j] = H[el.j][el.i];
		}
		
		if(testfl == false){
			for(auto j = 0; j < var_max; j++){
				for(auto i = 0; i < var_max; i++){
					if(H[j][i] == UNSET || std::isnan(H[j][i])) emsg("Hessian prob");
				}
			}
		}		
	}		

	mpi.bcast(H);

	return H;
}


/// Directly calculates the gradient in posterior probability (used for testing)
vector <double> Emulator::calculate_gradient_check(vector <double> value, const vector <unsigned int> &var_list) 
{	
	vector <double> gradient(var_list.size());
	
	auto d = finite_different_grid_size(value);
	for(auto i = 0u; i < var_list.size(); i++){
		auto v = var_list[i];
		cout << i << " " << d[v] << " i\n";
	
		auto value_st = value[v];
		value[v] += d[v]; 
		auto val_up = calculate_posterior(value);
		
		value[v] -= 2*d[v];
		auto val_down = calculate_posterior(value);
		
		gradient[i] = (val_up-val_down)/(2*d[v]);
			
		value[v] = value_st;
		update_param_from_value(value);
	}
	
	return gradient;
}


/// Works out a suitable size for doing finite difference estimates
vector <double> Emulator::finite_different_grid_size(const vector <double> &value)
{
	update_param_from_value(value);
	
	vector <double> d(nvar);
	for(auto i = 0; i < nvar; i++){  
		auto sd = 1.0;
		if(var[i].type == PARAM){
			const auto &par = model.param[var[i].ref];
			switch(par.prior_type){
				case FLAT_PRIOR: 	
					sd = (par.prior_val2-par.prior_val1)/2;
					break;
					
				case NORMAL_FROM_SD_PRIOR:
					sd = param_value[par.prior_sd_param];
					break;
					
				default:
					emsg("TO DO3");
					break;
			}		
		}
		
		d[i] = 0.0001*sd; 
		if(d[i] == 0) emsg("MAP");
	}

	return d;
}


/// A simple one dimentional emulator
void Emulator::simple()
{
	const auto T = 100;
	const auto step = 50;
	vector <double> data_t;
	vector <double> data;
		
	auto x = 0.0;
	for(auto t = 0; t <= T; t++){
		cout << t << " " << x << " x\n";
		if(t%step == 0){
			data_t.push_back(t);
			data.push_back(x);
		}
		
		x += normal_sample(0,1);
	}
	//emsg("d");
	
	auto ndata = data.size();
	vector < vector <double> > A;
	A.resize(ndata);
	for(auto d = 0; d < ndata; d++){
		A[d].resize(ndata);
	}
	
	auto delta = 0.0;

	for(delta = 1; delta < 100.0; delta += 0.1){
		for(auto d = 0; d < ndata; d++){
			for(auto dd = 0; dd < ndata; dd++){
				A[d][dd] = exp(-(data_t[d]-data_t[dd])*(data_t[d]-data_t[dd])/(delta*delta));
			}
		}
		auto A_inv = invert_matrix(A);

		//print_matrix("A",A);
		//emsg("L");
		auto sum = 0.0;
		for(auto d = 0; d < ndata; d++){
			for(auto dd = 0; dd < ndata; dd++){
				sum += data[d]*A_inv[d][dd]*data[dd];
			}
		}
		
		auto L = -(ndata/2.0)*log(2*M_PI) - 0.5*determinant_fast(A) - 0.5*sum;
		cout << delta << " " << L << " " << determinant_fast(A)  <<" del\n";
	}
	
	
	delta = 60;
	for(auto d = 0; d < ndata; d++){
		for(auto dd = 0; dd < ndata; dd++){
			A[d][dd] = exp(-(data_t[d]-data_t[dd])*(data_t[d]-data_t[dd])/(delta*delta));
		}
	}
	auto A_inv = invert_matrix(A);

	ofstream dat("data.txt");
	for(auto d = 0; d < ndata; d++){
		dat << data_t[d] << " " << data[d] << "\n";
	}
	
	ofstream em("em.txt");
	for(auto t = 0; t < T; t++){
		vector <double> c(ndata);
		for(auto d = 0; d < ndata; d++) c[d] = exp(-(t-data_t[d])*(t-data_t[d])/(delta*delta));
		
		auto mean = 0.0;
		for(auto d = 0; d < ndata; d++){
			for(auto dd = 0; dd < ndata; dd++){
				mean += c[d]*A_inv[d][dd]*data[dd];
			}
		}
		
		auto var = 1.0;
		for(auto d = 0; d < ndata; d++){
			for(auto dd = 0; dd < ndata; dd++){
				var -= c[d]*A_inv[d][dd]*c[dd];
			}
		}
		
		auto sd = sqrt(var+TINY);
		cout << t << " " << mean << " " << sd << " " << var << " \n";
		em << t << " " << mean << " " << mean-sd << " " <<  mean+sd << "\n";
	}
}


/// Shows the emulator fit for the simple 1d case
void Emulator::show_fit()
{
	tune_emulator();
	
	auto nds = data_sample.size();
	if(data_sample[0].em_var_value.size() != 1) emsg("Not one dimensional");
	
	if(mpi.core == 0){
		ofstream data("data.txt");
		for(const auto &ds : data_sample){
			data << ds.em_var_value[0] << " " << ds.post << endl;
		}
		
		ofstream em("em.txt");
		vector <double> em_var_value(nem_var);
		//for(auto var = 0.001; var < 1; var += 0.001){
		for(auto var = 0.001; var < 2; var += 0.001){
			em_var_value[0] = var;
			
			auto ev_em = evalulate_emulator(em_var_value,true);
		
			em << var << " " << ev_em.val << " " << ev_em.val - ev_em.sd << " " << ev_em.val + ev_em.sd << endl;
		}
	}
}


/// Calculates the mimimum 
void Emulator::minimum_prior_dist() const 
{
	auto min = LARGE;
	auto av = 0.0;
	for(auto d = 0; d < data_sample.size(); d++){
		auto min_p = LARGE;
		for(auto dd = 0; dd < data_sample.size(); dd++){
			if(d != dd){
				auto W = prior_dist(data_sample[d].em_var_value,data_sample[dd].em_var_value);
				if(W < min_p) min_p = W;
			}
		}
		av += min_p;
		if(min_p < min) min = min_p;
	}
	
	cout << "Minimum distance: " << min << endl;
	cout << "Average distance: " << av/data_sample.size() << endl;
}


/// Checks a two dimentional emulator
void Emulator::check_2d()
{
	auto LX = 20, LY = 20;

	for(auto v = 0; v < nhyper; v++) cout << hyper[v] << " hyper\n";

	auto beta_min = 0.04, beta_max = 0.06;
	auto m_min = 5.0, m_max = 9.0;
	
	cout << "beta,"<< beta_min << "," << beta_max << ",m,"<<m_min<< "," << m_max << endl;
	for(auto &ds : data_sample){
		cout << ds.em_var_value[0] << "," << ds.em_var_value[1] << "|"; 
	}
	cout << endl;
	
	for(auto j = 0; j < LY; j++){
		for(auto i = 0; i < LX; i++){
			param_value[0] = beta_min+ ((beta_max-beta_min)*i)/LX;
			param_value[1] = m_min+ ((m_max-m_min)*j)/LY;
			
			auto value = get_value_from_param();
			if(i != 0) cout << ",";
			
			vector <double> em_var_value(nem_var);
			 em_var_value[0]= param_value[0];
			 em_var_value[1]= param_value[1];
			 
			auto ev_em = evalulate_emulator(em_var_value,false);
			
			//cout << calculate_posterior(value) << " " << ev_em.val << "dif\n"; 
			cout << ev_em.val;
			//cout << calculate_posterior(value);
		}
		cout << "\n";
	}
	
}


/// Works out the prediction accuracy
void Emulator::check_prediction_accuracy()
{
	if(mpi.core == 0){
		ofstream scatter("scatter.txt");
		for(auto i = 0; i < model.nindividual; i++){
			scatter << model.individual[i].id << " ";
			for(auto ie = 0; ie < model.nind_effect; ie++){
				scatter << ind_value[i].ind_effect[ie] << " " << model.individual[i].ind_effect_value[ie] << " ";
			}
			scatter << endl;
		}
	}
	
	vector <double> grad;
	vector < vector <double> > H;
	
	auto post = model.posterior_grad_H(grad,H,param_value,ind_value,gr_ev,param_var_ref,ind_effect_var_ref,var.size(),true,true);
	
	
	cout << grad[100] << " grad\n";
}


/// Looks at how the likelihood changes for a given individual effect
void Emulator::check_one_ind_effect()
{
	auto i = 100;
	auto ie = 0;
	
	vector <double> grad;
	vector < vector <double> > H;
	cout << ind_value[i].ind_effect[ie] << "value\n";
	
	ofstream Lie("Lie.txt");
	for(auto val = -2.0; val < 2; val += 0.01){
		cout << val << "val\n";
		ind_value[i].ind_effect[ie] = val;
		
		auto post = model.posterior_grad_H(grad,H,param_value,ind_value,gr_ev,param_var_ref,ind_effect_var_ref,var.size(),false,false);
		
		Lie << val << " " << post << endl;
	}
}


/// Scans one of the variables
void Emulator::check_scan()
{
	/*
	param_value[em_var[0]] = 2;
	
	auto val = laplace_approximation();
	cout << val << "val\n";
	
	vector <double> indeff_store;
	
	
	for(auto loop = 0; loop < 2; loop++){
		auto val = laplace_approximation();
		cout << val << "val\n";
		
		ofstream ieout("ieout"+to_string(loop));
		for(auto i = 0; i < model.nindividual; i++){
			ieout << i << " " << model.individual[i].ind_effect_value[0] << " " << ind_value[i].ind_effect[0] << endl; 
		}
		
		if(loop == 0){
			indeff_store.resize(model.nindividual);
			for(auto i = 0; i < model.nindividual; i++){
				indeff_store[i] =  ind_value[i].ind_effect[0];
			}
		}
		
		if(loop == 1){
			for(auto i = 0; i < model.nindividual; i++){
				auto d = ind_value[i].ind_effect[0] - indeff_store[i];
				if(d > 0.1 || d < -0.1){
					cout << i << " " << d << " " << indeff_store[i] << " " << ind_value[i].ind_effect[0] << " d\n";
				}
			}
		}
	}
	*/
	
	
	ofstream scan("scan2.txt");
	for(auto var = 0.5; var < 2; var += 0.05){
		cout << var << " var\n";
		
		param_value[em_var[0]] = var; 
		auto val = laplace_approximation();
		
		scan << var << " " << val << "\n";
	}
}


/// Looks to understand how individual effects have different values
void Emulator::phase_transition()
{
	auto C = 2.3;
	auto sigma = 1.0;
	ofstream vout("vout.txt");
	for(auto v1 = -0.0; v1 < 1.0; v1+= 0.01){
		if(exp(v1) < C){
			auto v2 = log(C-exp(v1)); 
			auto L = -0.5*v1*v1/(sigma*sigma) -0.5*v2*v2/(sigma*sigma);
			vout << v1 << " " << v2 << " " << L << "\n"; 
		}
	}
	
	ofstream bifer("bifer.txt");
	for(auto C = 0.1; C < 10; C += 0.1){
		double Lmax = -LARGE, v1max, v2max;
		for(auto v1 = -6.0; v1 < 6.0; v1+= 0.001){
			if(exp(v1) < C){
				auto v2 = log(C-exp(v1)); 
				auto L = -0.5*v1*v1/(sigma*sigma) -0.5*v2*v2/(sigma*sigma);
				if(L > Lmax){ Lmax = L; v1max = v1; v2max = v2;}
			}
		}
		if(v1max > v2max){ auto temp = v1max; v1max = v2max; v2max = temp;}
		
		bifer << C << " " << exp(v1max) << " " << exp(v2max) << "\n";
		cout <<C << " " << exp(v1max) << " " <<  exp(v2max) << " fir\n";
	}
	
	ofstream bifer2("bifer2.txt");
	for(auto logC = -8.0; logC < 3; logC += 0.01){
	//for(auto logC = -0.0; logC < 0.0001; logC += 0.01){
		auto C = exp(logC);
		
		double Lmax = -LARGE, v1max, v2max;
		for(auto v1 = -16.0; v1 < 16.0; v1+= 0.001){
		
			if(exp(v1) < C){
				auto v2 = log(C-exp(v1)); 
				auto L = -0.5*v1*v1/(sigma*sigma) -0.5*v2*v2/(sigma*sigma);
				if(L >Lmax){ Lmax = L; v1max = v1; v2max = v2;}
			}
		}
		if(v1max > v2max){ auto temp = v1max; v1max = v2max; v2max = temp;}
		
		cout  << C << " " << exp(v1max) << " " << exp(v2max) << " " << Lmax << " logC\n";
		
		bifer2 << logC << " " << v1max << " " << v2max << "\n";
	}
}

	