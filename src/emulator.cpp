// This uses an emulator to construct the probability distribution for variance parameters

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
#include "timers.hpp"
#include "const.hpp"

Emulator::Emulator(const Model &model, int seed) : model(model)
{
	if(mpi.core == 0) cout << "Starting Emulator..." << endl;
	
	initialise_variable();

	prior_dist_min = 0.001;

	if(model.num_per_gen%mpi.ncore) emsg("'npoint' must be a multiple of the number of cores.");
	
	if(model.num_per_gen <= nfe+2) emsg("'npoint' must be higher");
		
	model.sample_initial_param(param_value,ind_value);
	
	model.initialise_group_events(gr_ev);
	
	seed += 10000*mpi.core;
	 
	set_seed(seed);                             // Sets the random seed
	srand(seed); 
}

void Emulator::inference()
{
	//phase_transition(); emsg("done");
	//check_grad_H(); emsg("done");
	//check_scan(); emsg("done");
	
	//param_value[em_var[0]] = 1; 
	//auto val = laplace_approximation();
	//check_prediction_accuracy(); emsg("done");
	
	//check_one_ind_effect(); emsg("done");

	//simple(); emsg("done");
	//check_2d(); emsg("done");
	
	//param_value[em_var[0]] = 1; 
	//auto val = laplace_approximation();
	//emsg("doo");
		
	if(nem_var == 0){
		if(mpi.core == 0) cout << "Emulator not needed for this model" << endl;
		auto samp = sample_no_emulator(model.output_samples);
		if(mpi.core == 0) trace_plot(samp);
		return;
	}
	
	for(auto g = 0; g < model.ngeneration; g++){
		if(mpi.core == 0) cout << endl << "Generation " << g << endl;
	
		vector < vector <double> > samp;
		if(mpi.core == 0){
			if(g == 0) samp = latin_hypercube_sample(model.num_per_gen);
			else samp = emulate_new_points(g);
		}
		mpi.bcast(samp);
	
		if(false){
			cout << mpi.core << ": ";
			for(auto s = 0; s < samp.size(); s++){
				cout << samp[s][0] << " ";
			}
			cout << "sa\n";
		}
		
		add_data_samp(samp);
		
		tune_emulator();
	}
	show_fit(); emsg("doe");
	
	//check_2d(); emsg("doe");
	
	//minimum_prior_dist(); 
	
	if(mpi.core == 0) cout << endl << "Run MCMC on emulator" << endl;
	
	if(mpi.core == 0){ 
		auto samp = mcmc(model.output_samples,1,10,false,1); 
		trace_plot(samp);
	}
}


/// Adds a set of data samples 
void Emulator::add_data_samp(vector < vector <double> > &samp)
{
	timer[TIME_ADD_DATA_SAMP].start();
	
	auto nsamp = samp.size();
	auto per_core = nsamp/mpi.ncore;
	
	vector <double> post(per_core);
	for(auto s = 0; s < per_core; s++){
		if(mpi.core == 0) cout << s << " / " << per_core << endl;
		auto ss = mpi.core*per_core+s;
		
		for(auto j = 0; j < nem_var; j++) param_value[em_var[j]] = samp[ss][j]; 
		
		post[s] = laplace_approximation();
	}
	
	auto post_tot = mpi.gather(post);
	
	if(mpi.core == 0){
		for(auto s = 0; s < nsamp; s++){
			DataSample ds;
			ds.em_var_value = samp[s];
			ds.post = post_tot[s];
			data_sample.push_back(ds);
		}
	}

	mpi.bcast(data_sample);

	timer[TIME_ADD_DATA_SAMP].stop();
}


/// Initialises a vector of model parameters
void Emulator::initialise_variable()
{
	auto params = split(model.params,',');
	
	param_var_ref.resize(model.nparam);
	for(auto th = 0; th < model.nparam; th++){
		param_var_ref[th] = UNSET;
			
		const auto &par = model.param[th];
		if(par.prior_type != FIXED_PRIOR && 
		  !(par.prior_type == FLAT_PRIOR && par.prior_val1 == par.prior_val2)){
				
			auto i = 0; while(i < params.size() && params[i] != par.name) i++;
			if(i < params.size()){
				params.erase(params.begin()+i);
				em_var.push_back(th);
			}
			else{
				param_var_ref[th] = var.size();
				Variable v;
				v.type = PARAM;
				v.ind = UNSET;
				v.ref = th;
				var.push_back(v);
			}
		}
	}
	nem_var = em_var.size();
	
	if(params.size() != 0) emsg("Not all variables in 'params' are recognised");
	
	ind_effect_var_ref.resize(model.N);
	for(auto i = 0; i < model.N; i++){
		ind_effect_var_ref[i].resize(model.nind_effect);
		for(auto ie = 0; ie < model.nind_effect; ie++){
			ind_effect_var_ref[i][ie] = var.size();
			
			Variable v;
			v.type = IND_EFFECT;
			v.ind = i;
			v.ref = ie;
			var.push_back(v);
		}
	}
	nvar = var.size();

	if(false){
		for(auto v = 0; v < 200; v++){
			if(var[v].type == IND_EFFECT){
				cout << v << " " << model.individual[var[v].ind].id << " " << var[v].ref << " h\n";
			}
		}
		emsg("do");
	}
	
	// constant, linear and quadratic fixed effects
	FE fix; fe.push_back(fix);
	
	for(auto v = 0; v < nem_var; v++){
		FE fix; fix.em_var_list.push_back(v);
		fe.push_back(fix);
	}
	
	for(auto v = 0; v < nem_var; v++){
		for(auto vv = v; vv < nem_var; vv++){
			FE fix; 
			fix.em_var_list.push_back(v);
			fix.em_var_list.push_back(vv);
			fe.push_back(fix);
		}
	}
	nfe = fe.size();
	
	// Stores which variables are parameters
	for(auto v = 0; v < nvar; v++){
		if(var[v].type == PARAM) var_param_list.push_back(v);
	}
	nvar_param_list = var_param_list.size();


	nugget = false;
	
	nhyper = nem_var; // If nugget included then nem_var+1;
	if(nugget == true) nhyper++;
	
	hyper.resize(nhyper);
	for(auto v = 0; v < nem_var; v++) hyper[v] = 1.0;
	if(nugget == true) hyper[nem_var] = 0;
	
	if(false){
		auto vmax = nvar; if(vmax > 100) vmax = 100;
		for(auto v = 0; v < vmax; v++){
			cout << v << " ";
			switch(var[v].type){
				case PARAM:
					cout << model.param[var[v].ref].name;
					break;
				
				case IND_EFFECT:
					cout << model.individual[var[v].ind].id << " " << model.ind_effect[var[v].ref].name;
					break;
			}
			cout << endl;
		}
		emsg("Done");
	}
}


/// Fits a MVN to the posterior (excluding variance components)
double Emulator::laplace_approximation()
{
	vector < vector <double> > H;
	
	auto init_type = 0;

	switch(init_type){
		case 0:  // Set to true initial values
			for(auto i = 0; i < model.nindividual; i++){
				for(auto ie = 0; ie < model.nind_effect; ie++){
					ind_value[i].ind_effect[ie] = model.individual[i].ind_effect_value[ie];
				}
			}
			break;
		
		case 1:
			model.ind_effect_sample(ind_value, param_value);      
			break;
	}
	
	timer[TIME_MAXIMISE].start();
	//auto post = maximise_posterior(H);
	auto post = maximise_posterior_BFGS(H);
	timer[TIME_MAXIMISE].stop();

	timer[TIME_DET].start();
	auto det = determinant_sparse(H);
	timer[TIME_DET].stop();
cout << post + 0.5*nvar*log(2*M_PI) - 0.5*det << "ans\n";
	return post + 0.5*nvar*log(2*M_PI) - 0.5*det;
}


/// Maximises the posterior probability
double Emulator::maximise_posterior(vector < vector <double> > &H)
{
	//-8062.72 
	auto value = get_value_from_param();
	
	if(nvar == 0) return calculate_posterior(value);
	
	vector <double> grad;
	double post, post_new;
	
	auto eta_ga = 1.0;
	
	auto loop = 0;
	do{
		post = model.posterior_grad_H(grad,H,param_value,ind_value,gr_ev,param_var_ref,ind_effect_var_ref,var.size(),true,true);
		if(std::isinf(post)){
			for(auto th = 0; th < model.nparam; th++) cout << model.param[th].name << " " << param_value[th] << "\n";
				post = model.posterior_grad_H(grad,H,param_value,ind_value,gr_ev,param_var_ref,ind_effect_var_ref,var.size(),true,true);
	
			emsg("inf");
		}
		cout << loop << " " << post << " " << " loop\n";

		//auto Hinv = invert_matrix_sparse(H);
		auto Hinv = invert_matrix(H);

		auto d = matrix_mult(Hinv,grad);
		for(auto &val : d) val *= -1;
		
		//cout << grad[0] << " " << d[0] << " " << Hinv[0][0] << " kk\n";
		auto eta = 1.0;
		auto type = NEWTONS_METHOD;
		if(dot_prod(d,grad) < 0){ type = GRADIENT_ASCENT; d = grad;}
		
		auto value_st = value;
		
		auto loop2 = 0;
		do{
			if(type == GRADIENT_ASCENT) eta = eta_ga;
			//cout << loop2 << " " << eta << "Loop2\n";

			for(auto v = 0; v < nvar; v++) value[v] += eta*d[v];
			
			update_param_from_value(value);
			//cout << d[0] << " " << value[0] << " "<< ind_value[0].ind_effect[0] << "indval\n"; 
	
			if(model.inbounds(param_value) == true){
				post_new = model.posterior_grad_H(grad,H,param_value,ind_value,gr_ev,param_var_ref,ind_effect_var_ref,var.size(),false,false);
			
				if(post_new > post-SMALL){
					if(type == GRADIENT_ASCENT) eta_ga *= 2;
					break;	
				}
			}
			
			value = value_st;
			eta *= 0.5;
			if(type == GRADIENT_ASCENT) eta_ga *= 0.25;
			loop2++; if(loop2 > 100) emsg("Problem with convergence1");
		}while(true);
	
		if(post_new - post < 0.0001 && eta == 1) break;
	
		loop++;
		if(loop > 100) break;
	}while(true);

	return post;
}


/// Maximises the posterior probability
double Emulator::maximise_posterior_BFGS(vector < vector <double> > &H)
{
	auto GD_mode = false;//true;
	
	//auto true_post = maximise_posterior(H);
	//cout << true_post << "true\n";
	//emsg("P");
	
	auto value = get_value_from_param();
	
	if(nvar == 0) return calculate_posterior(value);
	
	vector <double> grad, grad_new;
	double post, post_new;
	vector <double> s(nvar), y(nvar);
	
	//check_grad_H();
	
	auto Binv = model.construct_initial_Binv(param_value,param_var_ref,ind_effect_var_ref,nvar);
	
	if(GD_mode == true){
		for(auto v = 0; v < nvar; v++){
			for(auto vv = 0; vv < nvar; vv++){
				if(v == vv) Binv[v][vv] = 1;
				else Binv[v][vv] = 0;
			}
		}
	}
	
	auto eta = 1.0;
	
	post = model.posterior_grad_H(grad,H,param_value,ind_value,gr_ev,param_var_ref,ind_effect_var_ref,var.size(),true,false);
	for(auto &val : grad) val *= -1;
	
	//auto thresh = 0.0001;
	auto thresh = 0.000001;
	
	auto loop = 0;
	do{ 
		//if(loop < 2 && GD_mode == false) eta = 1;
		
		auto grmax = 0.0;
		for(auto val : grad){
			if(val > grmax) grmax = val;
			if(-val > grmax) grmax = -val;
		}
		//cout << endl << loop << " " << post << " " << eta << " " << grmax << " loop\n";
	
		auto d = matrix_mult(Binv,grad);
		for(auto &val : d) val *= -1;
		
		auto value_st = value;
		
		auto loop2 = 0;
		do{
			for(auto v = 0; v < nvar; v++){
				s[v] = eta*d[v];
				value[v] += s[v];
			}
			
			update_param_from_value(value);
			if(model.inbounds(param_value) == true){
				//cout << eta << "eta\n";
				//for(auto th = 0; th < model.nparam; th++) cout << model.param[th].name << " " << param_value[th] << "\n";
				
				post_new = model.posterior_grad_H(grad_new,H,param_value,ind_value,gr_ev,param_var_ref,ind_effect_var_ref,var.size(),true,false);
			
				//if(post_new > post-SMALL){
				if(post_new >= post){
					eta *= 2; if(eta > 1 && GD_mode == false) eta = 1;
					break;	
				}
			}
			
			value = value_st;
			eta *= 0.5; //cout << "fail\n";
			loop2++; if(loop2 > 100) emsg("Problem with convergence1");
		//}while(eta > SMALL || post_new - post > thresh);
		//}while(eta > SMALL || post_new - post >= 0);
		}while(true);
		
		//cout << post_new - post << " " << eta << "pos\n";
		for(auto &val : grad_new) val *= -1;
		
		if(false){
			for(auto v = 0; v < nvar; v++) y[v] = grad_new[v]-grad[v];
			
			auto sy = dot_prod(s,y);
			auto Binv_y = matrix_mult(Binv,y);
			auto y_Binv_y = dot_prod(y,Binv_y);
			
			auto fac = (sy + y_Binv_y)/(sy*sy);
			auto fac2 = 1.0/sy;
			for(auto v = 0; v < nvar; v++){
				for(auto vv = 0; vv < nvar; vv++){
					Binv[v][vv] += fac*s[v]*s[vv] - fac2*(Binv_y[v]*s[vv] + s[v]*Binv_y[vv]);
				}
			}
		}
		
		if(post_new - post < thresh) break;
		
		post = post_new;
		grad = grad_new;
		
		loop++;
		if(loop > 200) break;
	}while(true);
	cout << loop<< "loop\n";
	post = model.posterior_grad_H(grad,H,param_value,ind_value,gr_ev,param_var_ref,ind_effect_var_ref,var.size(),true,true);
	
	if(false){
		for(auto i = 200; i < 210; i++){
			for(auto j = 200; j < 210; j++){
				cout << Binv[j][i] << " ";
			}
			cout << "B\n";
		}

		for(auto v = 0; v < nvar; v++){
			for(auto vv = 0; vv < nvar; vv++){
				H[v][vv] *= -1;
			}
		}
			
		auto Hinv = invert_matrix(H);
		for(auto i = 200; i < 210; i++){
			for(auto j = 200; j < 210; j++){
				cout << Hinv[j][i] << " ";
			}
			cout << "H\n";
			
		}
		cout <<  " la\n";

		emsg("do");
	}
	
	for(auto i = 0; i < nvar; i++){
		//if(grad[i] < -0.00000001 || grad[i] > 0.00000001)  cout << i << " " << grad[i] << "grad\n";
	}
	
	return post;
}


/// Return a set of variable value based on param_value and ind_value
vector <double> Emulator::get_value_from_param()
{
	vector <double> value(nvar);
	
	for(auto i = 0; i < nvar; i++){
		switch(var[i].type){
			case PARAM:
				{
					auto th = var[i].ref;
					const auto &par = model.param[th];	
					value[i] = param_value[th];
				}
				break;
				
			case IND_EFFECT:
				auto ie = var[i].ref;
				value[i] = ind_value[var[i].ind].ind_effect[ie];
				break;
		}
	}
	
	return value;
}


// Uses the vector of varaible values to update param_value and ind_value
void Emulator::update_param_from_value(const vector <double> &value)
{
	for(auto v = 0; v < nvar; v++){
		switch(var[v].type){
			case PARAM:
				{
					auto th =  var[v].ref;
					param_value[th] = value[v];
				}
				break;
				
			case IND_EFFECT:
				{
					auto ie = var[v].ref;
					ind_value[var[v].ind].ind_effect[ie] = value[v];
				}
				break;
		}
	}
}


/// Samples from a latin hypercube using the prior
vector < vector <double> > Emulator::latin_hypercube_sample(unsigned int nsamp) const
{
	vector < vector <double> > result;
	
	vector < vector <double> > sub_samp;
	sub_samp.resize(nem_var);
	for (auto v = 0; v < nem_var; v++) {
		auto th = em_var[v];
		const auto &par = model.param[th];
		
		if(par.prior_type != FLAT_PRIOR) emsg("Should be flat prior");
		
		for(auto s = 0; s < nsamp; s++){
			auto d = (par.prior_val2 - par.prior_val1)/nsamp;
			auto min = par.prior_val1;
			sub_samp[v].push_back(min + (s+ran())*d);
		}
	}
		
	for(auto s = 0; s < nsamp; s++){
		vector <double> samp(nem_var);
		
		for (auto v = 0; v < nem_var; v++){
			auto i = int(ran()*sub_samp[v].size());
			samp[v] = sub_samp[v][i];
			sub_samp[v].erase(sub_samp[v].begin()+i);
		}
		
		result.push_back(samp);
	}
	
	return result;
}


/// Samples using the prior
vector < vector <double> > Emulator::prior_sample(unsigned int nsamp, bool spacing) const
{
	vector < vector <double> > result;
	
	for(auto s = 0; s < nsamp; s++){
		vector <double> samp(nem_var);
		
		bool too_close;
		auto loop = 0;
		do{
			for (auto v = 0; v < nem_var; v++) {
				auto th = em_var[v];
				const auto &par = model.param[th];
				
				if(par.prior_type != FLAT_PRIOR) emsg("Should be flat prior");
				samp[v] = par.prior_val1 + ran() * (par.prior_val2 - par.prior_val1);
			}
			
			too_close = false;
			if(spacing == true){
				for(auto ss = 0; ss < s; ss++){
					if(prior_dist(result[ss],samp) < prior_dist_min){ too_close = true; break;}
				}
			
				for(const auto &ds : data_sample){
					if(prior_dist(ds.em_var_value,samp) < prior_dist_min){ too_close = true; break;}
				}
			}
			
			loop++; if(loop == 100) emsg("Cannot find separated point");
		}while(too_close == true);
		
		result.push_back(samp);
	}
	
	return result;
}


/// Calculates the distance between two points based on prior ranges 
double Emulator::prior_dist(const vector <double> &value1, const vector <double> &value2) const
{
	auto W = 0.0;
	for(auto v = 0; v < nem_var; v++){
		auto th = em_var[v];
		const auto &par = model.param[th];		
		if(par.prior_type != FLAT_PRIOR) emsg("Should be flat prior");
		
		auto d = (value1[v] - value2[v])/(par.prior_val2 - par.prior_val1);
		W += d*d;
	}
	
	return sqrt(W);
}

// Used to order particles by EF
bool ds_ord(DataSample p1, DataSample p2)                      
{ return (p1.post > p2.post); };  


/// Culls points which are far from the peak
void Emulator::find_peak_cull_points()
{
	sort(data_sample.begin(),data_sample.end(),ds_ord);
	
	peak = data_sample[0].em_var_value;
	
	// This cut-off is based on the gamma distributed profile in the
	//negative log-probability for a MVN distribution
	
	auto dif = 5*(0.5*nem_var);
	if(dif < 10) dif = 10;
	
	auto post_min = data_sample[0].post - dif;
	auto d = 0; while(d < data_sample.size() && data_sample[d].post > post_min) d++;
	
	if(d < 2*nfe) d = 2*nfe;
	if(d < data_sample.size()) data_sample.resize(d);
	
	for(auto &ds : data_sample){
		ds.em_var_value_scale = ds.em_var_value;
		for(auto v = 0; v < nem_var; v++) ds.em_var_value_scale[v] /= peak[v];
	}
}


/// Based on data samples this tunes the emulator
void Emulator::tune_emulator()
{
	if(mpi.core == 0) cout << "Tune emulator" << endl;
	
	timer[TIME_TUNE].start();
	//for(auto v = 0; v < nem_var; v++) hyper[v] = 1;
	
	find_peak_cull_points();
	
	set_Z_g();
	
	auto L = marginal();
	double L_new;
	
	/*
	{
	ofstream hyp("hyp.txt");
	for(hyper[0] = 0.001; hyper[0] < 1; hyper[0] += 0.001){
		cout << hyper[0] << " " <<marginal() << "\n";
		hyp << hyper[0] << " " << marginal() << "\n";
	}
	}
	*/
	//emsg("P");
	
	auto eta_ga = 1.0;
	
	auto loop = 0;
	do{
		//cout << loop << " "<< hyper[0] << " " << L << " loop tune\n";
		
		auto grad = calculate_gradient_hyper(hyper);
		auto H = calculate_hessian_hyper(hyper);
		auto Hinv = invert_matrix(H);
		
		//print_vector("grad",grad);
		//print_matrix("H",H);
		
		if(false){
			auto d = 0.001;
			auto val = marginal();
			
			print_vector("grad",grad);
			print_matrix("H",H);
				
			if(true){
				hyper[0] += d; 
				auto valr = marginal();
				hyper[0] -= 2*d; 
				auto vall = marginal();
				hyper[0] += d; 
				
				cout << valr << " " << vall << " r l\n";
				cout << (valr-vall)/(2*d) << " gr\n";
				cout << (valr + vall - 2*val)/(d*d) << " h\n";
			}
			else{
				hyper[0] += d; 
				auto valr = marginal();
				hyper[0] -= 2*d; 
				auto vall = marginal();
				hyper[0] += d; 
				
				hyper[1] += d; 
				auto valu = marginal();
				hyper[1] -= 2*d; 
				auto vald = marginal();
				hyper[1] += d; 
				
				hyper[1] += d; 
				
				hyper[0] += d; 
				auto valur = marginal();
				hyper[0] -= 2*d; 
				auto valul = marginal();
				hyper[0] += d; 
					
				hyper[1] -= 2*d; 
						
				hyper[0] += d; 
				auto valdr = marginal();
				hyper[0] -= 2*d; 
				auto valdl = marginal();
				hyper[0] += d; 
				
				//cout << valup << " " << valdown << " " << val << " val\n";
				
				cout << (valr-vall)/(2*d) << " gr1\n";
				cout << (valu-vald)/(2*d) << " gr2\n";
				cout << (valr+vall - 2*val)/(d*d) << " h00\n";
				cout << (valu+vald - 2*val)/(d*d) << " h11\n";
				cout << 0.25*((valur+valdl) - (valul+valdr))/(d*d) << " h01\n";
			}
			//cout << grad[0] << " " << H[0][0] << " real\n";
			//cout << (valup-valdown)/(2*d) << " " << (valup+valdown-2*val)/(d*d) << "comp\n";
		
			//emsg("jj");
		}
		
		auto exit = false;
		
		if(mpi.core == 0){	
			auto d = matrix_mult(Hinv,grad);
			for(auto &val : d) val *= -1;
			
			auto eta = 1.0;
			auto type = NEWTONS_METHOD;
		
			if(dot_prod(d,grad) < 0){ type = GRADIENT_ASCENT; d = grad;}
			
			//cout << "TYPE " << type << endl;
			//cout << hyper[0] << " " << d[0] << " " << grad[0] << " do\n"; emsg("P");
		
			auto hyper_st = hyper;
		
			auto loop2 = 0;
			do{		
				if(type == GRADIENT_ASCENT) eta = eta_ga;
				//cout << eta << "eta\n";
			
				auto ill = false;
				for(auto v = 0; v < nhyper; v++){
					hyper[v] += eta*d[v];
					if(hyper[v] < 0) ill = true;
				}
				
				if(ill == false){
					L_new = marginal();
					if(L_new < L-SMALL) ill = true;
				}
				
				if(ill == false){
					if(type == GRADIENT_ASCENT) eta_ga *= 2;
					break;
				}
				
				hyper = hyper_st;
				eta *= 0.5;
				if(type == GRADIENT_ASCENT) eta_ga *= 0.25;
				
				loop2++; if(loop2 > 100) emsg("Problem with hyper convergence");
			}while(true);
			
			if(L_new - L < 0.0001 && type == NEWTONS_METHOD && (eta == 1 || loop > 3)) exit = true;
	
			L = L_new;
		}
		mpi.bcast(exit);
		
		if(exit == true) break;
		
		loop++; if(loop > 100) emsg("Problem with hyper convergence2");
	}while(true);
	
	if(mpi.core == 0) print_vector("hyper",hyper);
	
	/// Precalculates all quantities used to evaluate emulator
	auto nds = data_sample.size();
	
	auto A = calculate_A();
	Ainv = invert_matrix(A);
	Ainv_Z = matrix_mult(Ainv,Z);
	auto D = matrix_mult(ZT,Ainv_Z);
	Dinv = invert_matrix(D);
	auto Dinv_ZT = matrix_mult(Dinv,ZT);
	auto Dinv_ZT_Ainv = matrix_mult(Dinv_ZT,Ainv);
	betahat = matrix_mult(Dinv_ZT_Ainv,g);
	
	auto Z_betahat = matrix_mult(Z,betahat);
	auto dif = g;
	for(auto d = 0; d < nds; d++) dif[d] -= Z_betahat[d];
	q = matrix_mult(Ainv,dif);
	
	auto Ainv_g = matrix_mult(Ainv,g);
	
	auto Ainv_Z = matrix_mult(Ainv,Z);
	auto Ainv_Z_betahat = matrix_mult(Ainv_Z,betahat);
	
	auto sigma_sq = dot_prod(g,Ainv_g) - dot_prod(g,Ainv_Z_betahat);
	sigma_sq_fac = sigma_sq/(nds - nfe - 2);
	
	timer[TIME_TUNE].stop();
}


/// Sets fixed matrixes
void Emulator::set_Z_g()
{
	auto nds = data_sample.size();
	
	Z.resize(nds);
	for(auto d = 0; d < nds; d++){
		Z[d].resize(nfe);
		for(auto f = 0; f < nfe; f++){
			auto val = 1.0; 
			for(auto v : fe[f].em_var_list) val *= data_sample[d].em_var_value_scale[v];
			Z[d][f] = val;			
		}
	}

	ZT = transpose(Z);

	gmax = -LARGE;
	for(auto d = 0; d < nds; d++){
		if(data_sample[d].post > gmax) gmax = data_sample[d].post;
	}
	
	g.resize(nds);
	for(auto d = 0; d < nds; d++){
		g[d] = data_sample[d].post-gmax;
	}
}


/// Calculates the A matrix
vector < vector <double> > Emulator::calculate_A() const
{
	vector < vector <double> > A;
	
	auto nds = data_sample.size();
	
	A.resize(nds);
	for(auto d = 0; d < nds; d++){
		A[d].resize(nds);
		for(auto dd = 0; dd < nds; dd++){
			auto sum = 0.0;
			for(auto v = 0; v < nem_var; v++){
				auto dif = data_sample[d].em_var_value_scale[v] - data_sample[dd].em_var_value_scale[v];
				sum += (dif*dif)/(hyper[v]*hyper[v]);
			}
			
			A[d][dd] = exp(-sum);
			//if(nugget == true && d == dd) A[d][dd] += hyper[nem_var];
			//if(d == dd) A[d][dd] += 0.01;
			if(d == dd) A[d][dd] += 1;
		}
	}

	return A;
}


/// This calculates the GP likelihood marginilised over beta 
double Emulator::marginal()
{
	auto nds = data_sample.size();
	
	auto A = calculate_A();

	auto Ainv = invert_matrix(A);
	auto ZT_Ainv = matrix_mult(ZT,Ainv);
	auto D = matrix_mult(ZT_Ainv,Z);
	auto det_D = determinant_fast(D);
	auto Dinv = invert_matrix(D);
	auto det_A = determinant_fast(A);
	auto Dinv_ZT = matrix_mult(Dinv,ZT);
	auto Dinv_ZT_Ainv = matrix_mult(Dinv_ZT,Ainv);
	
	auto betahat = matrix_mult(Dinv_ZT_Ainv,g);
	
	auto Ainv_g = matrix_mult(Ainv,g);
	
	auto Ainv_Z = matrix_mult(Ainv,Z);
	auto Ainv_Z_betahat = matrix_mult(Ainv_Z,betahat);
	auto sigma_sq = dot_prod(g,Ainv_g) - dot_prod(g,Ainv_Z_betahat);
	if(sigma_sq < -SMALL) return -LARGE;
	if(sigma_sq < SMALL) sigma_sq = SMALL;

	return -0.5*det_A - 0.5*det_D - 0.5*(nds-nfe)*log(sigma_sq);
}

/// Determines the size for the grid on which gradients are calculated
vector <double> Emulator::finite_difference_grid_size() const
{
	vector <double> d(nhyper);
	
	for(auto v = 0; v < nem_var; v++){
		auto th = em_var[v];
		const auto &par = model.param[th];
			
		if(par.prior_type != FLAT_PRIOR) emsg("Should be flat prior");
		d[v] = 0.001*(par.prior_val2 - par.prior_val1)/data_sample.size();
	}
	
	if(nugget == true) d[nem_var] = 0.001;
	
	return d;
}

 
/// Calculates the gradient in the marginal likelihood
vector <double> Emulator::calculate_gradient_hyper(vector <double> &hyper)
{	
	mpi.bcast(hyper);

	auto d = finite_difference_grid_size();          // Determines grid size to use		
	
	auto per_core = (unsigned int)((nhyper+mpi.ncore-1)/mpi.ncore);
	
	vector <double> grad(per_core);
	
	for(auto k = 0u; k < per_core; k++){
		auto j = mpi.core*per_core + k;
	
		if(j < nhyper){
			hyper[j] += d[j]; 
			auto val_up = marginal();
			
			hyper[j] -= 2*d[j]; 
			auto val_down = marginal();
			
			grad[k] = (val_up-val_down)/(2*d[j]);
			
			hyper[j] += d[j]; 
		}
		else{
			grad[k] = UNSET;
		}
	}
	
	auto grad_tot = mpi.gather(grad);
	
	grad_tot.resize(nhyper);
	
	return grad_tot;
}

vector < vector <double> > Emulator::calculate_hessian_hyper(vector <double> &hyper)
{	
	mpi.bcast(hyper);

	auto d = 0.00001;      // Determines grid size to use	
	
	vector <HessianCalcList> hessian_calc_list; 
	for(auto i = 0; i < nhyper; i++){
		for(auto j = i; j < nhyper; j++){
			HessianCalcList gcl;
			gcl.v1 = i;
			gcl.v2 = j; 				
			hessian_calc_list.push_back(gcl);
		}
	}
	auto nhessian_calc_list = hessian_calc_list.size();
	
	auto per_core = (unsigned int)((nhessian_calc_list+mpi.ncore-1)/mpi.ncore);
	
	vector <double> hessian(per_core);
	
	for(auto k = 0u; k < per_core; k++){
		auto kk = mpi.core*per_core + k;
	
		if(kk < nhessian_calc_list){
			auto v1 = hessian_calc_list[kk].v1;
			auto v2 = hessian_calc_list[kk].v2;
			
			hyper[v1] += d; hyper[v2] += d;
			auto val_ur = marginal();
			hyper[v1] -= d; hyper[v2] -= d;
			
			hyper[v1] += d; hyper[v2] -= d;
			auto val_dr = marginal();
			hyper[v1] -= d; hyper[v2] += d;
			
			hyper[v1] -= d; hyper[v2] += d;
			auto val_ul = marginal();
			hyper[v1] += d; hyper[v2] -= d;
		
			hyper[v1] -= d; hyper[v2] -= d;
			auto val_dl = marginal();
			hyper[v1] += d; hyper[v2] += d;
			
			hessian[k] = (val_ur + val_dl - val_ul - val_dr)/(4*d*d);
		}
		else{
			hessian[k] = UNSET;
		}
	}
	
	auto hessian_tot = mpi.gather(hessian);
	
	vector < vector <double> > H;
	if(mpi.core == 0){
		H.resize(nhyper);
		for(auto v = 0; v < nhyper; v++) H[v].resize(nhyper);
	
		for(auto j = 0; j < hessian_calc_list.size(); j++){	
			auto v1 = hessian_calc_list[j].v1;
			auto v2 = hessian_calc_list[j].v2;
			H[v1][v2] = hessian_tot[j];
			if(v1 != v2) H[v2][v1] = hessian_tot[j];
		}
	}
		
	return H;
}


/// Evalulates the emulator at a new point
EmulatorEstimate Emulator::evalulate_emulator(const vector <double> &em_var_value, bool get_sd) const
{
	EmulatorEstimate em_est;
	if(nem_var == 0){
		em_est.val = 0;
		em_est.sd = 0;
		return em_est;
	}
	
	auto em_var_value_scale = em_var_value;
	for(auto v = 0; v < nem_var; v++) em_var_value_scale[v] /= peak[v];
	
	auto nds = data_sample.size();
	
	vector <double>  h(nfe);
	for(auto f = 0; f < nfe; f++){
		auto val = 1.0; for(auto v : fe[f].em_var_list) val *= em_var_value_scale[v];
		h[f] = val;
	}
	
	vector <double> c(nds);
	for(auto d = 0; d < nds; d++){
		auto sum = 0.0;
		for(auto v = 0; v < nem_var; v++){
			auto dif = data_sample[d].em_var_value_scale[v] - em_var_value_scale[v];
			sum += (dif*dif)/(hyper[v]*hyper[v]);
		}
		
		c[d] = exp(-sum);
	}
		
	em_est.val = gmax + dot_prod(h,betahat) + dot_prod(c,q);
	auto Ainv_c = matrix_mult(Ainv,c);
	
	if(get_sd == true){
		if(nfe != Ainv_Z[0].size()) emsg("Emulator problem 2");

		auto n = h;
		for(auto j = 0; j < nfe; j++){
			auto sum = 0.0;
			for(auto k = 0; k < nds; k++) sum += c[k]*Ainv_Z[k][j];
			n[j] -= sum; 
		}
		
		auto Dinv_n = matrix_mult(Dinv,n);
		
		auto var = sigma_sq_fac*(1- dot_prod(c,Ainv_c) + dot_prod(n,Dinv_n));
		if(var < 0) var = 0;
		
		em_est.sd = sqrt(var);
	}
	else em_est.sd = UNSET;
	
	return em_est;
}


/// This uses the emulator to add new points
vector < vector <double> > Emulator::emulate_new_points(unsigned int g)
{
	timer[TIME_EMULATE_NEW_POINT].start();
	
	auto d = 2*nfe; if(d >= data_sample.size()) d = data_sample.size()-1;
	auto difL = data_sample[0].post - data_sample[d].post;
	auto phi = 2*nem_var/(2.0*difL);
	if(phi > 1) phi = 1;
	
	//if(g == model.ngeneration-1) phi = 1;
	
	//cout << phi << " " << difL <<  " phi\n";
	
	//emsg("JJ");
	
	auto samp = mcmc(model.num_per_gen,0,10,true,phi);
	//print_matrix("new points added",samp);
	
	vector < vector <double> > samp_new;
	for(auto s = 0; s < model.num_per_gen; s++){
		vector <double> sa(nem_var);
		for(auto v = 0; v < nem_var; v++) sa[v] = samp[s][em_var[v]];
		
		samp_new.push_back(sa);
	}
	
	timer[TIME_EMULATE_NEW_POINT].stop();
	
	return samp_new;
}


/// Performs mcmc based on emulator results
vector < vector <double> > Emulator::mcmc(unsigned int nsamp, unsigned int nsamp_var, unsigned int thin, bool use_var, double phi)
{
	timer[TIME_EMULATOR_MCMC].start();
	    
	vector < vector <double> > cholesky_matrix;
	vector <double> em_var_value;
	
	if(nem_var > 0){
		vector < vector <double> > H;
		H.resize(nem_var);
		for(auto v = 0; v < nem_var; v++){
			H[v].resize(nem_var);
			for(auto vv = 0; vv < nem_var; vv++) H[v][vv] = 0;
		}
			
		for(auto f = 0; f < nfe; f++){
			if(fe[f].em_var_list.size() == 2){
				auto v = fe[f].em_var_list[0];
				auto vv = fe[f].em_var_list[1];
				auto val = phi*betahat[f]/(peak[v]*peak[vv]);
				H[v][vv] -= val;
				H[vv][v] -= val;
			}
		}
		
		auto C = invert_matrix(H);
		//print_matrix("C",C);
		//emsg("C");
		auto ill = false;
		for(auto v = 0; v < nem_var; v++){
			if(C[v][v] < 0) ill = true;
			
			auto th = em_var[v];
			const auto &par = model.param[th];
			auto d = par.prior_val2 - par.prior_val1;
			
			if(C[v][v] > d*d) ill = true;
			
			for(auto vv = v+1; vv < nem_var; vv++){
				if(C[v][vv]*C[v][vv] > C[v][v]*C[vv][vv]) ill = true;
			}
		}
		
		if(ill == true){
			for(auto v = 0; v < nem_var; v++){
				for(auto vv = 0; vv < nem_var; vv++){
					if(v == vv){ 
						auto th = em_var[v];
						const auto &par = model.param[th];
						auto d = par.prior_val2 - par.prior_val1;
						C[v][vv] = 0.01*d*d;
					}
					else C[v][vv] = 0;
				}
			}
		}

		cholesky_matrix = calculate_cholesky(C);
		
		auto max = -LARGE;
		for(const auto &ds : data_sample){
			if(ds.post > max){ max = ds.post;	em_var_value = ds.em_var_value;}
		}
		if(max == -LARGE) emsg("Emulator problem 9");
	}
	
	//print_matrix("C",C);

	auto j = 1.0;
	
	auto ntr = 0.0, nac = 0.0; 

	auto ev_em = evalulate_emulator(em_var_value,use_var);
	auto L = ev_em.val; if(use_var == true) L += ev_em.sd;
 
	vector < vector <double> > samp;
	
	auto burnin = int(0.1*nsamp*thin);
	
	for(auto s = 0; s < nsamp*thin; s++){
		//if(s%(nsamp*thin/10) == 0) cout << s << endl; 
		
		if(s%thin == 0){
			timer[TIME_EMULATOR_MCMC_SAMPLE].start();

			for(auto j = 0; j < nem_var; j++) param_value[em_var[j]] = em_var_value[j]; 
	
			if(nsamp_var > 0 && nvar_param_list > 0){
				auto psamp = MVN_param_sample(nsamp_var);
				for(const auto &ps : psamp) samp.push_back(ps);
			}
			else{
				samp.push_back(param_value);
			}
			
			timer[TIME_EMULATOR_MCMC_SAMPLE].stop();
		}
		
		auto em_var_value_prop = em_var_value;
		
		vector <double> norm(nem_var);	
		for(auto v = 0; v < nem_var; v++) norm[v] = normal_sample(0,1);
		
		auto shift = matrix_mult(cholesky_matrix,norm);
		
		for(auto v = 0; v < nem_var; v++) em_var_value_prop[v] += j*shift[v];
		//cout << em_var_value[0] << " " << em_var_value_prop[0] << "Prop\n";
		
		bool ill = false;
		for(auto v = 0; v < nem_var; v++){
			auto th = em_var[v];
			const auto &par = model.param[th];
			if(em_var_value_prop[v] > par.prior_val2 || em_var_value_prop[v] < par.prior_val1) ill = true;
		}			

		double al = 0.0, L_prop;
		if(ill == false){
			auto ev_em = evalulate_emulator(em_var_value_prop,use_var);
			L_prop = ev_em.val; if(use_var == true) L_prop += ev_em.sd;

			al = exp(phi*(L_prop-L));
		}
		
		ntr++;
		if(ran() < al){
			if(s < burnin) j *= 1.2;
			
			nac++;
			L = L_prop;
			em_var_value = em_var_value_prop;
		}
		else{
			if(s < burnin) j *= 0.9;
		}
	}
	
	cout << "MCMC Acceptance rate: " << nac/ntr << endl;

	timer[TIME_EMULATOR_MCMC].stop();

	return samp;
}


/// Given a fixed set of em_var this samples from the MVN distribution which approximates var
vector < vector <double> > Emulator::MVN_param_sample(unsigned int nsamp)
{	
	vector < vector <double> > H;
	auto post = maximise_posterior(H);

	for(auto v = 0; v < nvar; v++){
		for(auto vv = 0; vv < nvar; vv++){
			H[v][vv] *= -1;
		}
	}

	auto C = invert_matrix_sparse(H);

	// Makes a covariance matrix just for the parameters
	vector < vector <double> > C_param;
	C_param.resize(nvar_param_list);
	for(auto v = 0; v < nvar_param_list; v++){
		C_param[v].resize(nvar_param_list);
		for(auto vv = 0; vv < nvar_param_list; vv++){
			C_param[v][vv] = C[var_param_list[v]][var_param_list[vv]];
		}
	}
	
	//print_matrix("C",C); 
	//print_matrix("C_param",C_param); 
	
	auto Z = calculate_cholesky(C_param);
	
	vector < vector <double> > samp;
	
	for(auto k = 0; k < nsamp; k++){
		auto par = param_value;
		
		vector <double> norm(nvar_param_list);	
		for(auto v = 0; v < nvar_param_list; v++) norm[v] = normal_sample(0,1);

		auto shift = matrix_mult(Z,norm);
	
		for(auto v = 0; v < nvar_param_list; v++){
			auto vv = var_param_list[v];

			if(var[vv].type != PARAM) emsg("Should be parameter");
		
			par[var[vv].ref] += shift[v];
		}						
		samp.push_back(par);
	}
	
	return samp;
}


/// Samples in the case in which where are no emulator variables
vector < vector <double> > Emulator::sample_no_emulator(unsigned int nsamp)
{
	auto samp = MVN_param_sample(nsamp);
	return samp;
}


/// Generates a trace plot from a series of examples
void Emulator::trace_plot(const vector < vector <double> > &samp)
{	
	ofstream trace(model.output_dir + "/trace.txt");
		
	trace << "state"; 
	for (auto par : model.param) trace << '\t' << par.name; 
	trace << endl;
	
	for(auto s = 0; s < samp.size(); s++){
		trace << s;
		for(auto th = 0; th < model.nparam; th++) trace << "\t" << samp[s][th];
		trace << endl;
	}
}

