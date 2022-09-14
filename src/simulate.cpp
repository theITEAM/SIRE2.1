/// This includes functions for simulating from the model (for testing purposes)

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "simulate.hpp"
#include "model.hpp"
#include "timers.hpp"
#include "utils.hpp"
#include "matrix.hpp"


/// Constructor
Simulate::Simulate(Model &model) : model(model) {}


/// Runs a simulation
void Simulate::run() 
{
	param_value = model.prior_sample();

	set_value("beta", 0.05);
	set_value("m", 7);

	/*
	auto pheno_var = 1;
	auto h2 = 0.0001;
	set_value("omega_gg", h2*pheno_var);
	set_value("sigma_gg", (1-h2)*pheno_var);
	*/
	
	//set_value("omega_gg", 1.5);
	//set_value("omega_ff", 1.5);
	auto pheno_var_f = 3;
	auto h2 = 0.999;
	set_value("omega_ff", h2*pheno_var_f);
	set_value("sigma_ff", (1-h2)*pheno_var_f);
	
	/*
	auto pheno_var_g = 1;
	auto pheno_var_f = 3;
	auto h2 = 0.9999;
	
	set_value("omega_gg", h2*pheno_var_g);
	set_value("sigma_gg", (1-h2)*pheno_var_g);
		
	set_value("omega_ff", h2*pheno_var_f);
	set_value("sigma_ff", (1-h2)*pheno_var_f);
	
	set_value("omega_cor_gf", 0);
	set_value("sigma_cor_gf", 0);
		*/
		
	
	

	//set_value("beta", 0.05);
	//set_value("m", 7);
	//set_value("sigma",0.5);
	
	//set_value("omega_ff",1.0);
	//set_value("sigma_ff",1.0);
	   
	/*
	set_value("k", 3);
	*/
	/*
	set_value("beta", 0.2);
	set_value("trial_i", 0.4);
	set_value("trial_l", 0.2);
	set_value("latent_period", 4);
	set_value("detection_period", 6);
	set_value("recovery_period", 3);
	
	set_value("eta_shape", 3);
	set_value("rho_shape", 5);
	set_value("gamma_shape", 10);
	*/
	
	/*
	    set_value("a_g",0.3);
	    set_value("delta_g",0.5);

	    set_value("a_f",-0.4);
	    set_value("delta_f",0.1);

	    set_value("a_r",0.2);
	    set_value("delta_r",-0.9);
	*/
	/*

	    set_value("tau_EI",3);
	    set_value("k_EI",3);
	    set_value("tau_IV",4);
	    set_value("k_IV",3);
	    set_value("tau_VR",3);
	    set_value("k_VR",3);
	*/
	
	/*
	    set_value("omega_gg",1.0);
	    set_value("omega_ff",1.3);
	    set_value("omega_cor_gf",0.4);

	    set_value("sigma_gg",0.7);
	    set_value("sigma_ff",1.0);
	    set_value("sigma_cor_gf",0.2);
	*/

	/*
		set_value("omega_gg",1.0);
	    set_value("omega_rr",1.3);
	    set_value("omega_cor_gr",0.4);

	    set_value("sigma_gg",0.7);
	    set_value("sigma_rr",1.0);
	    set_value("sigma_cor_gr",0.2);
	*/


	/*
	    set_value("tau_EI",3);
	    set_value("k_EI",3);
	    set_value("tau_IR",7);
	    set_value("k_IR",3);
	*/

	//set_value("sigma",0.3);
	//set_value("nu",0.5);
	//set_value("nu_g",0.5);
	//set_value("nu_f",0.3);
	//set_value("nu_r",-0.2);

	//set_value("A_gg",0.7);
	//set_value("E_gg",1.0);



	// Used for E24
	/*
set_value("beta", 0.2);
	 set_value("latent_period",3);
	    set_value("detection_period",4);
			set_value("recovery_period",8);

		set_value("donor_inf",0.3);
		set_value("donor_rec",-0.2);

	set_value("sigma",0.1);
	
	 set_value("omega_gg",1.0);
	    set_value("omega_ff",1.3);
	    set_value("omega_rr",0.5);
	    set_value("omega_cor_gf",0.4);
	    set_value("omega_cor_gr",0.1);
	    set_value("omega_cor_fr",-0.2);

	    set_value("sigma_gg",0.7);
	    set_value("sigma_ff",1.0);
	    set_value("sigma_rr",0.9);
	    set_value("sigma_cor_gf",0.2);
	    set_value("sigma_cor_gr",-0.1);
	    set_value("sigma_cor_fr",-0.4);
*/

/*
	// Used for E23
set_value("beta", 0.2);
	 set_value("latent_period",3);
	    set_value("detection_period",4);
			set_value("recovery_period",8);

		set_value("donor_inf",0.3);
		set_value("donor_rec",-0.2);
		*/
		


/*
	// Used for E22
		set_value("beta", 0.2);
	 set_value("latent_period",3);
	    set_value("detection_period",4);
			set_value("recovery_period",8);
	*/		
/*  // Used for EX19-21
	set_value("beta", 0.05);
	 set_value("tau_EI",3);
	    set_value("tau_ID",4);
			set_value("tau_DR",3);
	    set_value("omega_gg",1.0);
	    set_value("omega_ff",1.3);
	    set_value("omega_rr",0.5);
	    set_value("omega_cor_gf",0.4);
	    set_value("omega_cor_gr",0.1);
	    set_value("omega_cor_fr",-0.2);

	    set_value("sigma_gg",0.7);
	    set_value("sigma_ff",1.0);
	    set_value("sigma_rr",0.9);
	    set_value("sigma_cor_gf",0.2);
	    set_value("sigma_cor_gr",-0.1);
	    set_value("sigma_cor_fr",-0.4);
	*/
	

	for (auto th = 0; th < model.nparam; th++) {
		const auto &par = model.param[th];
		if(par.prior_type == NORMAL_FROM_SD_PRIOR){
			param_value[th] = normal_sample(0, param_value[par.prior_sd_param]);
		}
	}	

	for (auto th = 0; th < model.nparam; th++)
		cout << model.param[th].name << " " << param_value[th] << " simulated value\n";

	ind_value.resize(model.N);
	for (auto i = 0; i < model.N; i++)
		ind_value[i].index = i;        // Sets index for each individual

	model.ind_effect_sample(ind_value, param_value);                  // Individual effects sampled from prior
	
	/*
	for(auto &ind : ind_value){
		for(auto val : ind.ind_effect) cout << val << "val\n";
	}
	emsg("P");
	*/

	model.set_individual_quantities(ind_value, param_value);          // Sets individual transition parameters and infectivity

	for (auto &gr : model.group) {                                    // Sets unlimited observation and inference rage
		gr.inference_range.tmax = LARGE;
		gr.observation_range.tmax = LARGE;
		//gr.inference_range.tmax = 10;
		//gr.observation_range.tmax = 10;
	}

	for (auto &ind : ind_value) {                                     // Removes all infections
		ind.infected = false;
		ind.trans_time.resize(model.ntrans);
		for (auto tr = 0; tr < model.ntrans; tr++)
			ind.trans_time[tr] = UNSET;
	}

	auto tmax = LARGE;                                                // Sets the end time

	for (auto &indiv : model.individual) {                            // Removes information from model
		indiv.status = NOT_INFECTED;
		//indiv.initial_comp = UNSET;
		for (auto tr = 0; tr < model.ntrans; tr++) {
			indiv.trans_time_range[tr].tmin = 0;
			indiv.trans_time_range[tr].tmax = tmax;
		}
	}

	for (const auto &gr : model.group) {
		auto N = gr.nind;                                               // This simulates from the model

		auto beta = param_value[model.beta_param];
		if (model.inf_model == FREQ_DEP)
			beta /= N;                      // If frequency dependent divided by group size
		if (model.group_effect.on == true)                              // Adds in the group effect
			beta *= exp(param_value[model.group_effect.param[gr.index]]);

		auto I = 0.0;                                                   // The number of infectious individuals
		for (auto i : gr.ind_ref) {
			if (model.individual[i].initial_comp != 0) {
				auto c = model.individual[i].initial_comp;
				if (c != 0) {
					add_infected(0, i, c);
					for (auto tr = 0; tr < c; tr++)
						I += ind_value[i].trans_infectivity_change[tr];
				}
			}
		}

		auto t = 0.0;
		do {
			auto next_time = LARGE;                                       // Store next non-infection event
			auto next_ind = 0;
			auto next_tr = 0;
			//cout << I << " I\n";
			auto S = 0.0;
			vector <double> S_st(N);
			for (auto j = 0; j < N; j++) {
				auto i = gr.ind_ref[j];
				auto &ind = ind_value[i];
				if (ind.infected == false)
					S += ind.susceptibility;
				else {
					auto tr = 0;
					while (tr < model.ntrans && ind.trans_time[tr] <= t)
						tr++;
					if (tr < model.ntrans) {
						if (ind.trans_time[tr] < next_time) {
							next_time = ind.trans_time[tr];
							next_ind = i;
							next_tr = tr;
						}
					}
				}
				S_st[j] = S;
			}

			double tt = LARGE;
			if (S > 0 && I > TINY)
				tt = t + exp_sample(beta * S * I);         // Samples time of next infection event

			if (next_time >= tmax && tt >= tmax)
				break;                  // Stops when final time reached

			if (tt < next_time) {
				auto z = ran() * S;
				auto j = 0;
				while (j < N && z > S_st[j])
					j++;
				if (j == N)
					emsg("Problem");

				auto i = gr.ind_ref[j];
				add_infected(tt, i, 1);
				I += ind_value[i].trans_infectivity_change[0];
				t = tt;
			} 
			else {                                                            // Performs an non-infection event
				//cout <<next_ind << " " <<  ind_value[next_ind].inf_single << " val\n";
			
				I += ind_value[next_ind].trans_infectivity_change[next_tr];
				t = next_time;
			}
		} while (true);
	}

/*
  for(auto i = 0; i < model.N; i++){			
		auto &ind = ind_value[i];	
		cout << i << ": ";
		for(auto j = 0; j < ind.trans_time.size(); j++) cout << " " << ind.trans_time[j];
		cout << "k\n";
	}
	*/
	
	output_datatable();
	//output_matrix();
	//output_pedigree();
}


/// Adds an infeced individual to the system
void Simulate::add_infected(const double t_inf, const int i, const int comp_init) {
	auto &ind = ind_value[i];
	auto &indiv = model.individual[i];

	ind.infected = true;

	indiv.status = INFECTED;

	for (auto tr = 0; tr < comp_init; tr++) {
		ind.trans_time[tr] = t_inf;
		indiv.trans_time_range[tr].tmin = t_inf;
		indiv.trans_time_range[tr].tmax = t_inf;
	}

	auto t = t_inf;
	for (auto tr = comp_init; tr < model.ntrans; tr++) {
		switch (model.trans[tr].type) {
			case GAMMA: {
				auto shape = param_value[model.trans[tr].shape_param];
				auto mean = ind.trans_mean[tr];
				t += gamma_sample(mean, shape);
			}
			break;
			
			case EXP: {
				auto mean = ind.trans_mean[tr];
				t += exp_sample(1.0/mean);
			}
			break;

			default:
				emsg("Not supported");
				break;
		}
		ind.trans_time[tr] = t;
		indiv.trans_time_range[tr].tmin = t;
		indiv.trans_time_range[tr].tmax = t;
	}
}


/// Sets the value of a parameter
void Simulate::set_value(string name, double val) {
	auto th = 0;
	while (th < model.nparam && model.param[th].name != name)
		th++;
	if (th == model.nparam)
		emsg("Cannot find '" + name + "'");

	param_value[th] = val;
}


/// Sets the intial value of the chain to that used in the simulation
void Simulate::set_initial_state(MCMC &mcmc) const {
	mcmc.param_value = param_value;
	mcmc.ind_value = ind_value;
}


/// This generates a file containing a datatable
void Simulate::output_datatable() {
	ofstream fout(model.output_dir + "/datatable.txt");

	fout << "<datatable id='1' group='2' initial_comp='3'";
	auto col = 4;

	if (model.nfixed_effect > 0) {
		for (auto fi = 0; fi < model.nfixed_effect; fi++){
			fout << " " << model.fixed_effect[fi].name << "='" << col << "'";
			col++;
		}
	}

	if (model.nsnp_effect > 0) {
		for (auto se = 0; se < model.nsnp_effect; se++){
			fout << " " << model.snp_effect[se].name << "='" << col << "'";
			col++;
		}
	}

	//fout << " comp_status='" << col << "'";
	//col++;
	col += model.ntrans;

	vector <int> ie_plot;
	for (auto ie = 0; ie < model.nind_effect; ie++) {
		//if (model.ind_effect[ie].covar_ref == 0) {
		ie_plot.push_back(ie);
		fout << " " << model.ind_effect[ie].name << "='" << col << "'";
		col++;
		//}
	}
	fout << ">" << endl;

	auto num_infected = 0;
	for (auto i = 0; i < model.N; i++) {
		auto &indiv = model.individual[i];
		auto &ind = ind_value[i];
		if (ind.infected == true)
			num_infected++;
	}
	cout << num_infected << " num_infected\n";

	for (auto i = 0; i < model.N; i++) {
		auto &indiv = model.individual[i];
		auto &ind = ind_value[i];

		fout << indiv.id;

		auto g = indiv.group;
		if (g == UNSET)
			fout << "\t" << "NA";
		else
			fout << "\t" << g + 1;

		auto c = indiv.initial_comp;

		if (c == UNSET)
			fout << "\tNA";
		else
			fout << "\t" << model.comp[c].name;

		if (model.nfixed_effect > 0)
			fout << "\t" << indiv.fixed_effect_X_unshifted[0];

		if(false){
			if(c == UNSET) fout << "\tNA";
			else{
				if(model.comp[c].name == "E") fout << "\t1"; else fout << "\t0";
			}
		}
		
		if (model.nsnp_effect > 0) {
			fout << "\t";
			switch (indiv.SNP_genotype[0]) {
				case AA:
					fout << "AA";
					break;
				case AB:
					fout << "AB";
					break;
				case BB:
					fout << "BB";
					break;
			}
		}

		/*
		    fout << "\t";
		    if(indiv.initial_comp ==  UNSET) fout << "NA";
		    else{
		    for(auto t = 0; t < 20; t += 2){
		    auto c = indiv.initial_comp;
		    auto tr = 0; while(tr < model.ntrans &&  ind.trans_time[tr] < t){ c = model.trans[tr].to; tr++;}
		    if(t != 0) fout << ",";
		    fout << "[";
		    if( model.comp[c].name == "S" ||  model.comp[c].name == "I") fout << "S|I";
		    else fout << "R";

		    // fout << model.comp[c].name << "
		    fout << ":" << t << "]";
		    }
		    }
		*/

		/*
		    auto Se = 0.8, Sp = 0.95;
		    fout << "\t";
		    if(indiv.initial_comp ==  UNSET) fout << "NA";
		    else{
		    for(auto t = 0; t < 20; t += 2){
		    auto c = indiv.initial_comp;
		    auto tr = 0; while(tr < model.ntrans &&  ind.trans_time[tr] < t){ c = model.trans[tr].to; tr++;}
		    if(t != 0) fout << ",";
		    fout << "[";
				if(c == 1){
					if(ran() < Se) fout << "+"; else fout << "-";
				}
				else{
					if(ran() < Sp) fout << "-"; else fout << "+";
				}

		    // fout << model.comp[c].name << "
		    fout << ":" << t << "]";
		    }
		    }
		*/

		for (auto tr = 0; tr < model.ntrans; tr++) {
			if (g == UNSET)
				fout << "\tNA";
			else {
				auto &gr = model.group[g];

				auto t = ind.trans_time[tr];
				if (t != UNSET && t > gr.observation_range.tmin && t < gr.observation_range.tmax) {
					if (ran() < 0.0)
						fout << "\t.";
					else
						fout << "\t" << t;
				} else
					fout << "\tno";
			}
		}

		for (auto ie : ie_plot)
			fout << "\t" << ind.ind_effect[ie];

		/*
		    auto q = 0.3; // A allele frequency
		    auto num = 0;
		    if(ran() < q) num++;
		    if(ran() < q) num++;

		    switch(num){
			case 0: fout << "\tAA"; break;
			case 1: fout << "\tAB"; break;
			case 2: fout << "\tBB"; break;
		    }
		*/

		fout << endl;
	}
	fout << "</datatable>" << endl;
}


/// Generates the pedigree from the the relationship matrix
void Simulate::output_pedigree() {
	ofstream fout(model.output_dir + "/ped.txt");
	cout << model.matrix[1].name << "na\n";

	auto A = model.matrix[1].A;
	auto Ainv = invert_matrix(A);

	fout << "<pedigree>" << endl;
	for (auto j = 0; j < model.N; j++) {
		vector <int> par;
		for (auto i = 0; i < model.N; i++) {
			if (Ainv[j][i] != 0) {
				if (j > i)
					par.push_back(i);


				cout << j << " " << i << " " << Ainv[j][i] << "\n";
			}
		}

		fout << model.individual[j].id << "\t";
		if (par.size() == 0)
			fout << ".\t.";
		if (par.size() == 1)
			fout << model.individual[par[0]].id << "\t.";
		if (par.size() == 2)
			fout << model.individual[par[0]].id << "\t" << model.individual[par[1]].id;
		if (par.size() > 2)
			emsg("Parent problem");
		fout << endl;
	}
	fout << "</pedigree>" << endl;
}


/// This generates a file containing the matrix is different formats
void Simulate::output_matrix() const {
	ofstream fout(model.output_dir + "/A_nonzero.txt");

	auto &M = model.matrix[1].A;

	for (auto j = 0; j < model.N; j++) {
		for (auto i = 0; i < model.N; i++) {
			if (M[j][i] < 1.001 && M[j][i] > 0.9999)
				M[j][i] = 1;
			if (M[j][i] < 0.501 && M[j][i] > 0.4999)
				M[j][i] = 0.5;
			if (M[j][i] < 0.2501 && M[j][i] > 0.24999)
				M[j][i] = 0.25;
		}
	}

	fout << "<A_nonzero name='A'>" << endl;
	for (auto j = 0; j < model.N; j++) {
		for (auto i = 0; i < model.N; i++) {
			if (M[j][i] != 0)
				fout << j << "\t" << i << "\t" << M[j][i] << endl;
		}
	}
	fout << "</A_nonzero>" << endl;

	ofstream fout2(model.output_dir + "/A.txt");

	fout2 << "<A name='A'>" << endl;
	for (auto j = 0; j < model.N; j++) {
		for (auto i = 0; i < model.N; i++) {
			if (i != 0)
				fout2 << "\t";
			fout2 << M[j][i];
		}
		fout2 << endl;
	}
	fout2 << "</A>" << endl;

	ofstream fout3(model.output_dir + "/Alist.txt");
	for (auto j = 0; j < model.N; j++) {
		fout3 << j << ": " << endl;
		fout3 << model.matrix[1].Ainvdiag[j] << endl;
		for (auto val : model.matrix[1].Ainvlist[j])
			fout3 << val.i << "," << val.val << "  ";
		fout3 << endl;
	}
}
