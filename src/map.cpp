// This uses a MAP approach using the CMA-ES algorithm

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "map.hpp"
#include "mpi.hpp"
#include "matrix.hpp"
#include "utils.hpp"

MAP::MAP(const Model &model) : model(model)
{
	if(mpi.core == 0) cout << "Starting MAP..." << endl;
	initialise_variable();

	model.sample_initial_param(param_value,ind_value);
	
	model.initialise_group_events(gr_ev);
	
	if(false){
		for(auto th = 0; th < model.nparam; th++){
			cout << model.param[th].name << " " << param_value[th] << " \n";
		}
		
		for(auto v = 0; v < nvar; v++){
			cout << v << " ";
			switch(var[v].type){
				case PARAM:
					cout << "PARAM " << model.param[var[v].ref].name;
					break;
				
				case IND_EFFECT:
					cout << "IND EFFECT " << var[v].ind << " " << var[v].ref;
					break;
			}
			cout << endl;
		}
		emsg("Done");
	}
	
	var_per_core = (unsigned int)((nvar+mpi.ncore)/mpi.ncore);
	//var_per_core2 = (unsigned int)((2*nvar+mpi.ncore-1)/mpi.ncore);
}


/// Initialises a vector of model parameters
void MAP::initialise_variable()
{
	param_var_ref.resize(model.nparam);
	for(auto th = 0; th < model.nparam; th++){
		const auto &par = model.param[th];
		if(par.prior_type != FIXED_PRIOR && 
		  !(par.prior_type == FLAT_PRIOR && par.prior_val1 == par.prior_val2) &&
			par.type != COVAR_MATRIX){	
			param_var_ref[th] = var.size();
			Variable v;
			v.type = PARAM;
			v.ind = UNSET;
			v.ref = th;
			var.push_back(v);
		}
		else{
			param_var_ref[th] = UNSET;
		}
	}
	
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
	
	// Sets up group_variable, which lists all variable which affect a group
	group_variable.resize(model.ngroup);
	for(auto g = 0; g < model.ngroup; g++){
		var_add(group_variable[g].var_list,model.beta_param);
		
		for(const auto &fe : model.fixed_effect){
			var_add(group_variable[g].var_list,fe.param);
		}
			
		for(const auto &tr : model.trans){
			if(tr.type != INFECTION) var_add(group_variable[g].var_list,tr.mean_param);
		}
		
		if(model.group_effect.on == true){
			var_add(group_variable[g].var_list,model.group_effect.param[g]);
		}
		
		for(auto i : model.group[g].ind_ref){
			for(auto ie = 0; ie < model.nind_effect; ie++){
				indeff_add(group_variable[g].var_list,i,ie);
			}
		}
	}
	
	for(auto &gv : group_variable){
		gv.mask.resize(nvar);
		for(auto v = 0; v < nvar; v++) gv.mask[v] = false;
		for(auto v : gv.var_list) gv.mask[v] = true;
	}
	
	if(false){
		for(auto g = 0; g < model.ngroup ; g++){
			cout << "Group " << g << ": ";
			const auto &gv = group_variable[g];
			for(auto j = 0; j < gv.var_list.size(); j++){
				auto v = gv.var_list[j];
				
				switch(var[v].type){
					case PARAM:
						cout << model.param[var[v].ref].name;
						break;
					
					case IND_EFFECT:
						cout << model.individual[var[v].ind].id << " " << model.ind_effect[var[v].ref].name;
						break;
				}
				cout << ", ";
			}
			cout << endl;
			/*
			cout << "Mask: ";
			for(auto v = 0; v < nvar; v++){
				if(gv.mask[v] == true) cout << "true, ";
				else cout << "false, ";
			}
			cout << endl;
			*/
		}
		emsg("Done");
	}
}


/// Adds a value onto a vector (if one of the parameters)
void MAP::var_add(vector <unsigned int> &vec, unsigned int th) const 
{
	auto v = 0; while(v < nvar && !(var[v].type == PARAM && var[v].ref == th)) v++;
	if(v < nvar){
		vec.push_back(v);
	}
}


/// Adds a value onto a vector (if one of the parameters)
void MAP::indeff_add(vector <unsigned int> &vec, unsigned int i, unsigned int ie) const 
{
	auto v = 0; 
	while(v < nvar && !(var[v].type == IND_EFFECT && var[v].ind == i && var[v].ref == ie)) v++;
	if(v < nvar){
		vec.push_back(v);
	}
}

// Used to order particles by EF
bool VariableSample_ord(VariableSample p1,VariableSample p2)                      
{ return (p1.L > p2.L); };  


/// Performs the CMA-ES algorithm
void MAP::cmaes(int seed)
{
	seed += 10000*mpi.core;
	set_seed(seed);                             // Sets the random seed
	srand(seed); 
	
	auto npart = (unsigned int)(4+int(5*log(nvar))); 
	npart = int((npart+mpi.ncore-1)/mpi.ncore)*mpi.ncore;
			
	const auto c_c_limit = 0.9;              // These values are used to limit the system in the case of few parameters
	const auto c_sigma_limit = 0.9;
	const auto c1_limit = 0.3;
	const auto c_mu_limit = 0.3;
	const auto c_s_limit = 0.3;

	auto sigma = 1.0;
	vector <double> p_sigma(nvar), p_c(nvar);
	for(auto i = 0; i < nvar; i++){
		p_sigma[i] = 0.0; p_c[i] = 0.0;
	}
	
	cout << model.likelihood(param_value,ind_value,gr_ev) << " Like before" << endl;
	 
	/*
	ofstream lp("lp2.txt");
	for(mean[0] = 0.01; mean[0] < 0.1; mean[0] += 0.001){
		cout << mean[0] << " " << calculate_posterior(mean) << "\n";
		lp << mean[0] << " " << calculate_posterior(mean) << "\n";
	}
	return;
	emsg("do");
	*/
	auto mean = get_value_from_param();
	auto C = get_C_from_param();
	
	update_param_from_value(mean);
	cout << model.likelihood(param_value,ind_value,gr_ev) << " Like after\n";
	
	if(false){
		for(auto i = 0; i < nvar; i++){
			switch(var[i].type){
				case PARAM:
					cout << model.param[var[i].ref].name << " " << mean[i] << " " << C[i][i] << "\n";
					break;
			}
		}
		emsg("do");
	}
	
	/* Calculates quantities used later */
	auto mu = (unsigned int)(npart/2);
	vector <double> w(mu);
	//for(auto i = 0; i < mu; i++) w[i] = mu-i;
	for(auto i = 0; i < mu; i++) w[i] = log(mu+0.5) - log(i+1);

	auto sum = 0.0; for(auto i = 0; i < mu; i++) sum += w[i];
	for(auto i = 0; i < mu; i++) w[i] /= sum;
	
	auto su = 0.0; for(auto i = 0; i < mu; i++) su += w[i]*w[i];
	auto mu_w = 1.0/su;

	auto npart_per_core = (unsigned int)(npart/mpi.ncore);
	vector <VariableSample> ps_per_core(npart_per_core);
		
	if(mpi.core == 0) cout << endl << "Start Inference..." << endl;
	
	vector <double> Lbest_store;
	
	vector <Generation> generation;
	
	auto g = 0;
	do{
		mpi.bcast(mean);
		mpi.bcast(C);
		mpi.bcast(sigma);

		Generation gen;
	
		/* Samples parameters */

		generate_samples(ps_per_core,mean,C,sigma);
		
		auto ps = mpi.gather_psamp(ps_per_core);

		if(mpi.core == 0){
			sort(ps.begin(),ps.end(),VariableSample_ord);   
			
			auto Lav = 0.0; for(auto i = 0; i < mu; i++) Lav += w[i]*ps[i].L;
	
			cout << "Generation " << g << " - Best log(Post. prob.): " << ps[0].L <<  "   Sigma: " << sigma << endl;
		
			Lbest_store.push_back(ps[0].L);
			
			gen.Lav = Lav;
		
			auto mean_store = mean;

			/* Updates the mean */
			for(auto i = 0; i < nvar; i++){
				auto sum = 0.0; for(auto k = 0; k < mu; k++) sum += ps[k].value[i]*w[k];
				mean[i] = sum; 
			}		
			
			/* Update p_sigma */
			vector <double> dif(nvar);
			for(auto i = 0; i < nvar; i++) dif[i] = mean[i] - mean_store[i]; 
			
			auto Cinvsqrt = invert_matrix_square_root(C);
			
			auto vec = matrix_mult(Cinvsqrt,dif);
			
			auto c_sigma = 3.0/nvar; if(c_sigma > c_sigma_limit) c_sigma = c_sigma_limit;
			
			vector <double> p_sigma_new(nvar);
			for(auto i = 0; i < nvar; i++){
				p_sigma_new[i] = (1-c_sigma)*p_sigma[i] + sqrt(1-(1-c_sigma)*(1-c_sigma))*sqrt(mu_w)*vec[i]/sigma;
			}
			p_sigma = p_sigma_new;
			
			/* Update p_c */
			auto alpha = 1.5;
			
			auto p_sigma_mag = 0.0; for(auto i = 0; i < nvar; i++) p_sigma_mag += p_sigma[i]*p_sigma[i];
			p_sigma_mag = sqrt(p_sigma_mag);
			
			auto ind = 1u; if(p_sigma_mag > alpha*sqrt(nvar)) ind = 0;
				
			auto c_c = 4.0/nvar; if(c_c > c_c_limit) c_c = c_c_limit;
				
			vector <double> p_c_new(nvar);
			for(auto i = 0; i < nvar; i++){
				p_c_new[i] = (1-c_c)*p_c[i] + ind*sqrt(1-(1-c_c)*(1-c_c))*sqrt(mu_w)*dif[i]/sigma;
			}
			p_c = p_c_new;
			
			/* Update C */
			auto c1 = 2.0/(nvar*nvar); if(c1 > c1_limit) c1 = c1_limit;
			auto c_mu = mu_w/(nvar*nvar); if(c_mu > c_mu_limit) c_mu = c_mu_limit; if(c_mu > 1-c1) c_mu = 1-c1;
			auto c_s = (1-ind)*c1*c_c*(2-c_c); if(c_s > c_s_limit) c_s = c_s_limit;
			
			vector < vector <double> > C_new;
			C_new.resize(nvar);
			for(auto j = 0; j < nvar; j++){
				C_new[j].resize(nvar);
				for(auto i = 0; i < nvar; i++){
					auto sum = 0.0;
					for(auto k = 0; k < mu; k++){
						sum += w[k]*(ps[k].value[i]-mean_store[i])*(ps[k].value[j]-mean_store[j])/(sigma*sigma);
					}
					
					C_new[j][i] = (1-c1-c_mu+c_s)*C[j][i] + c1*p_c[i]*p_c[j] + c_mu*sum;
				}
			}
			C = C_new;
			
			/* Update sigma */
			//auto d_sigma = 1;
			auto d_sigma = 1 + 2*max(0.0,sqrt((mu_w-1.0)/(nvar+1.0))-1) + c_sigma;
		
			auto EN = sqrt(nvar)*(1-(1.0/(4*nvar)) + 1.0/(21*nvar*nvar));
			
			//EN *= 2;
			auto fac = (c_sigma/d_sigma)*((p_sigma_mag/EN) - 1);
			//fac -= 0.01;
			if(fac > 0.2) fac = 0.2;
			sigma *= exp(fac);
		}
		
		generation.push_back(gen);

		g++; //if(g == 20) return;
	}while(terminate_generation(Lbest_store) == false);

	if(mpi.core == 0) cout << "Generating posterior samples..." << endl;

	mpi.bcast(mean);
	mpi.bcast(C);

	mpi.barrier();
	if(mpi.core == 0) cout << "Calculate covariance..." << endl;
	auto C_new = calculate_covariance_martrix(mean);	

	auto CC = calculate_cholesky(C_new);	
	if(CC[0][0] == UNSET){
		if(mpi.core == 0) cout << "Cannot calculate directly, scale covariance estimate from CMA-ES..." << endl;
		sigma = scale_covariance_martrix(mean,C);
	}
	else{
		C = C_new; sigma = 1;
	}
	
	if(mpi.core == 0){
		output_parameter_samples(mean,C,sigma);
	}
}

/// Generates samples based on a mean and covariance matrix
void MAP::generate_samples(vector <VariableSample> &ps_per_core, const vector <double> &mean, const vector < vector <double> > &C, const double sigma)
{
	auto psamp = sample_param(mean,C,sigma,ps_per_core.size());

	for(auto i = 0; i < ps_per_core.size(); i++){
		const auto &param = psamp[i];
		ps_per_core[i].value = psamp[i];
		ps_per_core[i].L = calculate_posterior(psamp[i]);
	}
}


/// Samples a set of parameters from a given set of priors
vector < vector <double> > MAP::sample_param(const vector <double> &mean, const vector < vector <double> > &C, const double sigma, const unsigned int num) const
{
	auto cholesky_matrix = calculate_cholesky(C);
	 
	auto nparam = model.param.size();
	
	vector <double> norm(nvar);	
		
	vector < vector <double> > psamp;
	
	for(auto i = 0; i < num; i++){
		vector <double> value(nvar);
		
		for(auto v = 0; v < nvar; v++) norm[v] = normal_sample(0,1);

		for(auto v = 0; v < nvar; v++){
			auto dva = 0.0; for(auto v2 = 0; v2 <= v; v2++) dva += cholesky_matrix[v][v2]*norm[v2];

			value[v] += mean[v] + sigma*dva;
		}
		
		for(auto v = 0; v < nvar; v++){   // This reflects values to ensure inside prior
			if(var[v].type == PARAM){
				const auto &pa = model.param[var[v].ref];
				if(pa.prior_type == FLAT_PRIOR){
					auto val = value[v];
					auto val1 = pa.prior_val1, val2 = pa.prior_val2;
			
					bool reflect;
					do{
						reflect = false;
						if(val > val2){
							val = val2 - (val-val2);
							reflect = true;
						}
						else{
							if(val < val1){
								val = val1 + (val1-val);
								reflect = true;
							}
						}
					}while(reflect == true);
					value[v] = val;
				}
			}
		}

		psamp.push_back(value);
	}
	
	return psamp;
}

/// Determines when generations are stopped
bool MAP::terminate_generation(const vector <double> &EFbest_store)
{
	auto term = false;
	if(mpi.core == 0){
		auto ng = EFbest_store.size();
		if(ng >= ML_GENERATION_TERM_COND){
			auto av = 0.0;
			for(auto i = ng-ML_GENERATION_TERM_COND; i < ng; i++) av += EFbest_store.size();
			av /= ML_GENERATION_TERM_COND;
			
			auto tol = 0.01;
			auto value = EFbest_store[ng-1];
			auto i = ng-ML_GENERATION_TERM_COND; 
			while(i < ng && EFbest_store[i] > value-tol && EFbest_store[i] < value+tol) i++;
			if(i == ng) term = true;
		}
	}
	mpi.bcast(term);
	
	return term;
}


/// Calculates the Hessian matrix and uses this to estimate the covariance matrix
vector < vector <double> > MAP::calculate_covariance_martrix(const vector <double> &mean)
{
	if(mpi.core == 0) cout << "Calculating Hessian" << endl;	
	//auto H = calculate_hessian(mean,FULL);
	auto H = calculate_hessian(mean,DIAG);
	for(auto &Hrow : H){
		for(auto &val : Hrow){
			if(val == UNSET) val = 0;
			//else cout << val << "\n";
		}
	}
	
	if(mpi.core == 0) cout << "Invert matrix" << endl;	
	auto C = invert_matrix(H);
	for(auto &Crow : C){
		for(auto &val : Crow) val *= -1;
	}
	
	return C;
}


/// Return a set of variable value based on param_value and ind_value
vector <double> MAP::get_value_from_param()
{
	vector <double> value(nvar);
	
	for(auto i = 0; i < nvar; i++){
		switch(var[i].type){
			case PARAM:
				{
					auto th = var[i].ref;
					const auto &par = model.param[th];
					
					value[i] = param_value[th];
					/*
					if(par.prior_type == NORMAL_FROM_SD_PRIOR){
						auto sd = param_value[par.prior_sd_param];
						value[i] /= sd;
					}
					*/
				}
				break;
				
			case IND_EFFECT:
				auto ie = var[i].ref;
				//const auto &ind_eff = model.ind_effect[ie];
				//const auto th = model.covariance[ind_eff.covar_ref].var_param[ind_eff.covar_num];
				//cout << model.param[th].name << "nan\n";
				
				value[i] = ind_value[var[i].ind].ind_effect[ie];
				break;
		}
	}
	
	return value;
}


/// Gets the initial covariance matrix based on the current parameter values
vector < vector <double> > MAP::get_C_from_param()
{
	vector < vector <double> > C;
	C.resize(nvar); 
	for(auto j = 0; j < nvar; j++){
		C[j].resize(nvar);
		for(auto i = 0; i < nvar; i++) C[j][i] = 0;
	}

	for(auto i = 0; i < nvar; i++){
		switch(var[i].type){
			case PARAM:
				{
					auto th = var[i].ref;
					const auto &par = model.param[th];
				
					switch(par.prior_type){
						case FLAT_PRIOR:
							C[i][i] = (par.prior_val2-par.prior_val1)*(par.prior_val2-par.prior_val1)/4;
							break;
					
						case NORMAL_FROM_SD_PRIOR:
							C[i][i] = 1;
							break;
						
						default:
							emsg("TO DO1");
							break;
					}
				}
				break;
				
			case IND_EFFECT:
				C[i][i] = 1;
				break;
		}
	}
	
	return C;
}

	
// Uses the vector of varaible values to update param_value and ind_value
void MAP::update_param_from_value(const vector <double> &value)
{
	for(auto v = 0; v < nvar; v++){
		update_param_from_value(value,v);
	}
}

void MAP::update_param_from_value(const vector <double> &value, unsigned int v)
{
	switch(var[v].type){
		case PARAM:
			{
				auto th =  var[v].ref;
				param_value[th] = value[v];
				
				//const auto &par = model.param[th];
				/*
				if(par.prior_type == NORMAL_FROM_SD_PRIOR){
					auto sd = param_value[par.prior_sd_param];
					param_value[th] *= sd;
				}
				*/
			}
			break;
			
		case IND_EFFECT:
			auto ie = var[v].ref;
			ind_value[var[v].ind].ind_effect[ie] = value[v];
			break;
	}
}


/// Estimates posterior from a parameter set
double MAP::calculate_posterior(const vector <double> &value)
{
	update_param_from_value(value);
	
	auto L = model.likelihood(param_value,ind_value,gr_ev);
//L = 0.0;
	auto L_ind_eff = model.calculate_L_ind_effect(ind_value, param_value);
	for(auto val : L_ind_eff) L += val;
	
	L += model.calculate_prior(param_value);

	return L;	
}


/// Works out the scaling factor for the covaiance matrix
double MAP::scale_covariance_martrix(const vector <double> &mean, const vector < vector <double> > &C)
{
	auto H_diag = calculate_hessian(mean,DIAG);

	double ratio;
	if(mpi.core == 0){
		auto Cinv = invert_matrix(C);
	
		if(false){
			print_matrix("Cinv",Cinv);
			print_matrix("H_diag",H_diag);
		
			for(auto i = 0; i < nvar; i++){
				cout << Cinv[i][i] << " " << H_diag[i][i] << " " << -H_diag[i][i]/Cinv[i][i] << " compare" << endl;
			}
		}
		
		auto sum_log_ratio = 0.0;
		for(auto i = 0; i < nvar; i++){
			auto ra = -Cinv[i][i]/H_diag[i][i];
			if(ra < 0.0001) ra = 0.0001;
			if(ra > 10000) ra = 10000;
			sum_log_ratio += log(ra);
		}
		sum_log_ratio /= nvar;
		
		ratio = exp(sum_log_ratio);
	}
	auto sigma = sqrt(ratio);
	
	mpi.bcast(sigma);

	return sigma;
}


/// Works out a suitable size for doing finite difference estimates
vector <double> MAP::finite_different_grid_size(const vector <double> &value)
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
		
		d[i] = 0.001*sd;
		if(d[i] == 0) emsg("MAP");
	}

	return d;
}

/// Calculates the Hessian matrix
vector < vector <double> > MAP::calculate_hessian(const vector <double> &mean, MatrixType mattype)
{
	auto d = finite_different_grid_size(mean);          // Determines grid size to use
	
	// Should check not overlapping boundaries...
	
	auto testfl = false;
	
	vector <MatrixElement> matrix_element;
	if(false){
		for(auto j = 0; j < nvar; j++){
			auto imax = nvar;
			if(mattype == DIAG) imax = j+1; 
			for(auto i = j; i < imax; i++){
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
			for(auto j = 0; j < nvar; j++){
				switch(mattype){
					case FULL:
						for(auto i = 0; i < nvar; i++){
							if(H[j][i] == UNSET || std::isnan(H[j][i])) emsg("MAP");
						}
						break;
					
					case DIAG:
						if(H[j][j] == UNSET || std::isnan(H[j][j])) emsg("MAP");
						break;
				}
			}
		}		
	}		

	mpi.bcast(H);

	return H;
}

	
/// Outputs randomly generated samples
void MAP::output_parameter_samples(const vector <double> &mean, const vector < vector <double> > &C, const double sigma)
{
	auto psamp = sample_param(mean,C,sigma,model.output_samples);
	
	ofstream samples(model.output_dir + "/param_samples.txt");

	samples << "state";
	for(auto th = 0; th < model.nparam; th++){
		samples << "\t" << model.param[th].name;
	}
	samples << endl;
	
	for(auto i = 0; i < model.output_samples; i++){
		update_param_from_value(psamp[i]);
		
		samples << i;
		for(auto th = 0; th < model.nparam; th++){
			samples << "\t" << param_value[th];
		}
		samples << endl;
	}
}

/// Performs the gradient descent algorithm
void MAP::gradient_descent()
{
	auto var_block = set_var_block();
	auto nvar_block = var_block.size();
	
	vector <double> eta(nvar_block);  // Sets the initial value for the step size
	for(auto &val : eta) val = 1;
		
	if(false){
		ifstream param_store("param_store.txt");
		for(auto th = 0; th < model.nparam; th++){
			param_store >> param_value[th];
		}
	}

	auto value = get_value_from_param();
	
	//check_gradient();
	
	auto L = calculate_posterior(value);
	
/*
for(auto sig = 0.01; sig < 0.5; sig += 0.01){
	value[0] = sig;
	update_param_from_value(value);
	cout << sig << " " << param_value[3] << " ch\n";
	auto L_ind_eff = model.calculate_L_ind_effect(ind_value, param_value);
	auto Lin = 0.0; for(auto val : L_ind_eff) Lin += val;
	
	cout << sig << "  post: " << calculate_posterior(value) << "   like:" 
	    <<  model.likelihood(param_value,ind_value,gr_ev) << " prio:"
			<<  model.calculate_prior(param_value) << " " << Lin << " g\n";
}

emsg("do");
*/

	vector <double> store; 

	auto accept_flag = false;

	auto loop = 0u;
	do{
		for(auto b = 0; b < nvar_block; b++){
		//for(auto b = 0; b < 1; b++){
			if(mpi.core == 0) cout << "Block: " << b << "  Iteration: " << loop << "  Postetior: " << L << "  Jump size: " << eta[b] << endl;
			
			auto L_store = L;
			auto value_store = value;
			const auto &var_list = var_block[b].var_list;
			
			auto grad =	calculate_gradient(value,var_list);
		
			if(false && loop%10 == 0){
				auto grad_slow = calculate_gradient_slow(value,var_list);
				check_same(grad,grad_slow);
			}
			
			auto fail = false;
			do{
				cout << eta[b] << " eta\n";
				fail = false;
				
				if(mpi.core == 0){
					for(auto i = 0; i < var_list.size(); i++){
						auto v = var_list[i];
						value[v] += eta[b]*grad[i];
					}
					
					update_param_from_value(value);
					if(model.inbounds(param_value) != true) fail = true;
				}
				mpi.bcast(fail);
			
				if(mpi.core == 0){
					if(fail == false){
						L = calculate_posterior(value);
						if(L < L_store) fail = true;	
					}
				}
				mpi.bcast(fail);
		
				if(mpi.core == 0){
					if(fail == true){
						eta[b] *= 0.5; 
						L = L_store; 
						value = value_store;
					}
					else{
						accept_flag = true;
						eta[b] *= 1.5;
					}
				}
			}while(fail == true);
		}
		
		auto finish = false;
		if(mpi.core == 0){
			store.push_back(L);
			if(store.size() > 10 && accept_flag == true){
				if(L < store[store.size()-10] +0.001) finish = true;
			}
		}
		mpi.bcast(finish);
		if(finish == true) break;
		
		loop++;
	}while(loop < 1000);
	
	/*
	mpi.bcast(value);

	mpi.barrier();
	if(mpi.core == 0) cout << "Calculate covariance..." << endl;
	vector < vector <double> > C;
	C.resize(nvar);
	for(auto v = 0; v < nvar; v++){
		C[v].resize(nvar);
		for(auto vv = 0; vv < nvar; vv++){
			if(v == vv) C[v][vv] = 0.01*value[v]*value[v];
			else C[v][vv] = 0;
		}
	}
	*/
	
	/*
	auto C = calculate_covariance_martrix(value);	
	auto CC = calculate_cholesky(C);	
	if(CC[0][0] == UNSET) cout << "cholesky problem" << endl;
	*/

/*
	if(mpi.core == 0){
		if(false){
			ofstream param_store("param_store.txt");
			param_store.precision(17);
			for(auto th = 0; th < model.nparam; th++){
				param_store << param_value[th] << endl;
			}
		}
		
		output_parameter_samples(value,C,1);
	}
	*/
}


/// Calculates the gradient in the likelihood for the non-fixed parameters
vector <double> MAP::calculate_gradient(vector <double> value, const vector <unsigned int> &var_list) 
{	
	mpi.bcast(value);
	auto d = finite_different_grid_size(value);          // Determines grid size to use
	
	vector <GradCalcList> grad_calc_list; 
	for(auto g = 0; g < model.ngroup; g++){
		for(auto i = 0; i < var_list.size(); i++){
			auto v = var_list[i];
			if(group_variable[g].mask[v] == true){
				GradCalcList gcl; gcl.g = g; gcl.v = v; gcl.i = i;
				grad_calc_list.push_back(gcl);
			}
		}
	}
		
	auto list_per_core = (unsigned int)((grad_calc_list.size()+mpi.ncore-1)/mpi.ncore);
	
	vector <double> grad(list_per_core);
	
	for(auto k = 0u; k < list_per_core; k++){
		auto j = mpi.core*list_per_core + k;
	
		if(j < grad_calc_list.size()){
			auto g = grad_calc_list[j].g;
			auto v = grad_calc_list[j].v;
			
			auto value_st = value[v];
			value[v] += d[v]; 
			update_param_from_value(value,v);
			auto val_up = model.likelihood_group(g,param_value,ind_value,gr_ev);
			
			value[v] -= 2*d[v]; 
			update_param_from_value(value,v);
			auto val_down = model.likelihood_group(g,param_value,ind_value,gr_ev);
			
			grad[k] = (val_up-val_down)/(2*d[v]);
			
			value[v] = value_st;
			update_param_from_value(value,v);
		}
		else{
			grad[k] = UNSET;
		}
	}
	
	auto grad_tot = mpi.gather(grad);
	
	vector <double> gradient(var_list.size());
	if(mpi.core == 0){
		for(auto &val : gradient) val = 0;

		for(auto k = 0; k < grad_calc_list.size(); k++){
			gradient[grad_calc_list[k].i] += grad_tot[k];
		}
		
		auto ind_effect_flag = false;
		
		for(auto i = 0; i < var_list.size(); i++){
			auto v = var_list[i];
			switch(var[v].type){
				case PARAM:
					gradient[i] += model.calculate_prior_gradient(var[v].ref,param_value);
					break;
					
				case IND_EFFECT:
					ind_effect_flag = true;
			}
		}
		
		if(ind_effect_flag == true){
			auto inv_cov_matrix = model.calculate_inv_cov_matrix(param_value);
		
			for(auto i = 0; i < var_list.size(); i++){
				auto v = var_list[i];
				switch(var[v].type){
					case IND_EFFECT:
						gradient[i] += model.ind_effect_gradient(var[v].ind,var[v].ref,ind_value,inv_cov_matrix);
						break;
						
					case PARAM:
						auto th = var[v].ref;
						if(model.param[th].type == COVAR_MATRIX){
							gradient[i] += model.covar_param_gradient(th,ind_value,param_value);
						}
						break;
				}
			}
		}
	}
		
	return gradient;
}


/// Calculates the Hessian and uses the fact the fact that likelihood independent between groups
vector < vector <double> > MAP::calculate_hessian_fast(vector <double> value) 
{	
	mpi.bcast(value);
	auto d = finite_different_grid_size(value);          // Determines grid size to use
	
	vector <HessianCalcList> hessian_calc_list; 
	for(auto g = 0; g < model.ngroup; g++){
		for(auto i = 0; i < group_variable[g].var_list.size(); i++){
			for(auto j = i; j < group_variable[g].var_list.size(); j++){
				HessianCalcList gcl;
				gcl.g = g;
				gcl.v1 = group_variable[g].var_list[i];
				gcl.v2 = group_variable[g].var_list[j]; 				
				hessian_calc_list.push_back(gcl);
			}
		}
	}
	
	auto list_per_core = (unsigned int)((hessian_calc_list.size()+mpi.ncore-1)/mpi.ncore);
	
	vector <double> hessian(list_per_core);
	
	for(auto k = 0u; k < list_per_core; k++){
		auto j = mpi.core*list_per_core + k;
	
		if(j < hessian_calc_list.size()){
			auto g = hessian_calc_list[j].g;
			auto v1 = hessian_calc_list[j].v1;
			auto v2 = hessian_calc_list[j].v2;
				
			value[v1] += d[v1]; value[v2] += d[v2];
			update_param_from_value(value,v1); update_param_from_value(value,v2);
			auto val_ur = model.likelihood_group(g,param_value,ind_value,gr_ev);
			value[v1] -= d[v1]; value[v2] -= d[v2];
			
			value[v1] += d[v1]; value[v2] -= d[v2];
			update_param_from_value(value,v1); update_param_from_value(value,v2);
			auto val_dr = model.likelihood_group(g,param_value,ind_value,gr_ev);
			value[v1] -= d[v1]; value[v2] += d[v2];
			
			value[v1] -= d[v1]; value[v2] += d[v2];
			update_param_from_value(value,v1); update_param_from_value(value,v2);
			auto val_ul = model.likelihood_group(g,param_value,ind_value,gr_ev);
			value[v1] += d[v1]; value[v2] -= d[v2];
		
			value[v1] -= d[v1]; value[v2] -= d[v2];
			update_param_from_value(value,v1); update_param_from_value(value,v2);
			auto val_dl = model.likelihood_group(g,param_value,ind_value,gr_ev);
			value[v1] += d[v1]; value[v2] += d[v2];
			hessian[k] = (val_ur + val_dl - val_ul - val_dr)/(4*d[v1]*d[v2]);
			update_param_from_value(value,v1); update_param_from_value(value,v2);
		}
		else{
			hessian[k] = UNSET;
		}
	}
	
	auto hessian_tot = mpi.gather(hessian);
	
	vector < vector <double> > H;
	if(mpi.core == 0){
		H.resize(nvar);
		for(auto v = 0; v < nvar; v++){
			H[v].resize(nvar);
			for(auto vv = 0; vv < nvar; vv++) H[v][vv] = 0;
		}
	
		for(auto j = 0; j < hessian_calc_list.size(); j++){	
			auto v1 = hessian_calc_list[j].v1;
			auto v2 = hessian_calc_list[j].v2;
			//cout << hessian_tot[j] << " tot\n";
			if(v1 == 103) cout << v2 << " " << hessian_tot[j] << "cont\n";
			H[v1][v2] += hessian_tot[j];
			//if(v1 != v2) H[v2][v1] += hessian_tot[j];
		}
		
		auto inv_cov_matrix = model.calculate_inv_cov_matrix(param_value);
		for (auto cv = 0; cv < model.ncovariance; cv++) {
			const auto &cov = model.covariance[cv];
			
			const auto &M = model.matrix[cov.matrix];
			for(auto i1 = 0; i1 < cov.E; i1++){
				auto ie1 = cov.ind_effect_ref[i1];
				for(auto i2 = 0; i2 < cov.E; i2++){
					auto ie2 = cov.ind_effect_ref[i2];

					auto variance = inv_cov_matrix[cv].M[i1][i2];
					cout << variance << "va\n";
					for(auto i = 0; i < model.N; i++){
						auto v1 = ind_effect_var_ref[i][ie1];
						auto v2 = ind_effect_var_ref[i][ie2];
						
						H[v1][v2] -= M.Ainvdiag[i]*variance;// zz
						//cout << v1 << " " << v2 << " " << M.Ainvdiag[i] << " " << variance << " diag\n";
						for(const auto el : M.Ainvlist[i]){
							auto j = el.i;
							auto v2 = ind_effect_var_ref[j][ie2];
							H[v1][v2] -= el.val*variance;
							//cout << v1 << " " << v2 << " " <<el.val << " " << variance << " nondiag\n";
						}
					}
				}
			}
		}
	}
		
	return H;
}


/// Directly calculates the gradient in posterior probability (used for testing)
vector <double> MAP::calculate_gradient_slow(vector <double> value, const vector <unsigned int> &var_list) 
{	
	vector <double> gradient(var_list.size());
	
	auto d = finite_different_grid_size(value);
	for(auto i = 0u; i < var_list.size(); i++){
		cout << i << " i\n";
		auto v = var_list[i];
		
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


/// Checks the gradient in the posterior is correctly being calculated
void MAP::check_gradient()
{
	cout << "Check gradient\n";
	auto value = get_value_from_param();
	
	vector <unsigned int> var_list;
	auto vmax = nvar; if(vmax > 20) vmax = 20;
	for(auto v = 0; v < vmax; v++) var_list.push_back(v);
	
	auto grad_slow = calculate_gradient_slow(value,var_list);
		
	auto grad = calculate_gradient(value,var_list);
	
	if(mpi.core == 0){
		for(auto i = 0; i < var_list.size(); i++){
			auto v = var_list[i];
			if(var[v].type == PARAM) cout << model.param[var[v].ref].name;
			else cout << "INF EFF";
			cout << " " << v << " " << grad[v] << " " << grad_slow[v] << " Compare gradients" << endl;
		}
	}
	mpi.barrier();
	emsg("done");
}


/// Checks the gradient in the ie is correctly being calculated 
// (note all non liklihood terms must be turned off in the gradient calculation)
void MAP::check_gradient_ie()
{
	auto grad =	calculate_gradient_ie();
	
	vector <unsigned int> var_list;
	auto vmax = nvar; if(vmax > 120) vmax = 120;
	for(auto v = 0; v < vmax; v++) var_list.push_back(v);
	
	auto value = get_value_from_param();
	
	auto grad_comp = calculate_gradient(value,var_list);
	if(mpi.core == 0){
		for(auto i = 0; i < var_list.size(); i++){
			auto v = var_list[i];
			if(var[v].type == IND_EFFECT){
				cout <<  grad_comp[v] << " " << grad[var[v].ind][var[v].ref] << " Compare\n";
			}
		}
	}
	
	mpi.barrier();
	emsg("done");
}


/// Set variable blocks
vector <VarBlock> MAP::set_var_block() const
{
	vector <VarBlock> var_block;
	
	auto GEflag = false;
	auto IEflag = false;
	for(auto v = 0; v < nvar; v++){
		switch(var[v].type){
			case PARAM:
				{
					const auto &par = model.param[var[v].ref];
					if(par.type == GROUP_EFFECT) GEflag = true;
					else{
						VarBlock vb; vb.var_list.push_back(v);
						var_block.push_back(vb);
					}
				}
				break;
				
			case IND_EFFECT:
				IEflag = true;
				break;
		}
	}
	
	if(GEflag == true){
		VarBlock vb;
		for(auto v = 0; v < nvar; v++){
			if(var[v].type == PARAM){
				const auto &par = model.param[var[v].ref];
				if(par.type == GROUP_EFFECT) vb.var_list.push_back(v);
			}
		}
		var_block.push_back(vb);
	}	
	
	if(IEflag == 1){
		VarBlock vb;
		for(auto v = 0; v < nvar; v++){
			if(var[v].type == IND_EFFECT){
				vb.var_list.push_back(v);
			}
		}
		var_block.push_back(vb);
	}
	
	if(false){
		for(const auto &bl : var_block){
			cout << "Block: ";
			for(auto v : bl.var_list) cout << v << ", "; cout << endl;  
		}
		emsg("done");
	}
	
	return var_block;
}


/// Performs a restricted likelihood (by integrating out latent variables)
void MAP::REML()
{

}


/// Calculates the likelihood with the individual effects marginailised
double MAP::likelihood_ie_marginal()
{
	for(auto i = 0; i < model.N; i++){
		for(auto ie = 0; ie < model. nind_effect; ie++){
			ind_value[i].ind_effect[ie] = 0;
		}
	}

	auto L0 = model.likelihood(param_value,ind_value,gr_ev);
	
	auto grad =	calculate_gradient_ie();
	
	cout << model.matrix.size() << " " << grad[0].size() << " mat\n";
	
	const auto &covar = model.covariance[0];
	
	auto var = param_value[covar.var_param[0]];
	cout <<  var << " " << model.matrix[0].name << " " <<  model.matrix[1].name << " na\n";
	const auto &M = model.matrix[covar.matrix].A;
	
	vector <double> delta_a(model.N);
	for(auto j = 0; j < model.N; j++){
		auto sum = 0.0;
		for(auto i = 0; i < model.N; i++)	sum += M[j][i]*grad[i][0];
		delta_a[j] = sum*var;		
	}
	
	auto i = 102;
	
	auto Lmax = -LARGE;
	auto Prmax = -LARGE;
	auto postmax = -LARGE;
	for(auto a = -2.0; a < 2; a += 0.01){
		cout << a << " a\n";
		ind_value[i].ind_effect[0] = a;
		auto L = model.likelihood(param_value,ind_value,gr_ev);
		if(L > Lmax) Lmax = L;
			
		auto L_ind_eff = model.calculate_L_ind_effect(ind_value, param_value);
		if(L_ind_eff[0] > Prmax) Prmax = L_ind_eff[0];
		
		if(L + L_ind_eff[0] > postmax) postmax = L + L_ind_eff[0];
	}
	
	ofstream ypl("ypl.txt");
	
	for(auto a= -2.0; a < 2.0; a += 0.01){
		ind_value[i].ind_effect[0] = a;
		
		auto L = model.likelihood(param_value,ind_value,gr_ev);
		auto L_ind_eff = model.calculate_L_ind_effect(ind_value, param_value);
		
		ypl << a << " " << L-Lmax << " " << L_ind_eff[0]-Prmax << " " << L+L_ind_eff[0]-postmax << " " << grad[i][0]*a << "\n";
		//cout << a << " "<< model.likelihood(param_value,ind_value,gr_ev) << "\n";
	}
}


/// Calculates the gradient in the likelihood for individual effects
vector < vector <double> > MAP::calculate_gradient_ie() 
{	
	auto d = 0.001;
	
	vector <GradIECalcList> grad_ie_calc_list; 
	
	for(auto i = 0; i < model.N; i++){
		for(auto ie = 0; ie < model.nind_effect; ie++){
			GradIECalcList gcl; gcl.i = i; gcl.ie = ie; 
			grad_ie_calc_list.push_back(gcl);
		}
	}
		
	auto list_per_core = (unsigned int)((grad_ie_calc_list.size()+mpi.ncore-1)/mpi.ncore);
	
	vector <double> grad(list_per_core);
	
	for(auto k = 0u; k < list_per_core; k++){
		auto j = mpi.core*list_per_core + k;
	
		if(j < grad_ie_calc_list.size()){
			auto i = grad_ie_calc_list[j].i;
			auto ie = grad_ie_calc_list[j].ie;
			
			auto g = model.individual[i].group;
			if(g == UNSET) grad[k] = 0;
			else{
				auto value_st = ind_value[i].ind_effect[ie];
				ind_value[i].ind_effect[ie] += d; 
				auto val_up = model.likelihood_group(g,param_value,ind_value,gr_ev);
			
				ind_value[i].ind_effect[ie] -= 2*d; 
				auto val_down = model.likelihood_group(g,param_value,ind_value,gr_ev);
			
				grad[k] = (val_up-val_down)/(2*d);
			
				ind_value[i].ind_effect[ie] = value_st;
			}
		}
		else{
			grad[k] = UNSET;
		}
	}
	
	auto grad_tot = mpi.gather(grad);
	
	vector < vector <double> > gradient;
	
	if(mpi.core == 0){
		gradient.resize(model.N);
		for(auto i = 0; i < model.N; i++){
			gradient[i].resize(model.nind_effect);
			for(auto &val : gradient[i]) val = 0;
		}
		
		for(auto k = 0; k < grad_ie_calc_list.size(); k++){
			gradient[grad_ie_calc_list[k].i][grad_ie_calc_list[k].ie] += grad_tot[k];
		}
	}
		
	return gradient;
}


/// Checks gradients and Hessian are correctly calculatued
void MAP::check_grad_H()
{
	vector <double> grad;
	vector < vector <double> > H;
	
	auto post = model.posterior_grad_H(grad,H,param_value,ind_value,gr_ev,param_var_ref,ind_effect_var_ref,var.size(),true,true);
	
	auto value = get_value_from_param();
	
	vector <unsigned int> var_list;
	auto vmax = nvar; //if(vmax > 120) vmax = 120;
	for(auto v = 0; v < vmax; v++) var_list.push_back(v);
	
	auto grad_ch2 = calculate_gradient_slow(value,var_list);
	auto grad_ch = calculate_gradient(value,var_list);
	
	for(auto v = 0; v < vmax; v++){
		cout << v << " "<< grad_ch[v] << " "<< grad[v] << " " << grad_ch2[v] << "\n";
		auto d = grad_ch[v] - grad[v];
		if(d > SMALLISH || d < -SMALLISH){
			cout << v << " " << grad_ch[v] << " " <<  grad[v] << "grad\n";
			cout << "IND EFFECT " << var[v].ind << " " << var[v].ref << " ind\n";
			emsg("wrong grad");	
		}			
	}
	
	cout << "Hessian\n";
	auto H_ch = calculate_hessian_fast(value);
	
	for(auto v = 0; v < nvar; v++){
		for(auto vv = 0; vv < nvar; vv++){
			auto d = H_ch[v][vv] - H[v][vv];
			if(d > SMALLISH || d < -SMALLISH) emsg("wrong H");		
		}
	}
	
	/*
	auto vv = 103;
	for(auto v = 0; v < 200; v++){
		cout << v << " " <<  H[v][vv] << " " << H_ch[v][vv] << "check H\n";
	}
	*/
	
}


