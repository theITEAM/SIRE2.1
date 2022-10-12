/// This includes functions for simulating from the model (for testing purposes)

// To generate files from the template use:
// ./sire cloglog/sus_template.xml 0 0

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
void Simulate::run(string file_in) 
{
	param_value = model.prior_sample();

	set_value("beta", 0.05);
	//set_value("m", 7);

	/*
	auto pheno_var = 1;
	auto h2 = 0.0001;
	set_value("omega_gg", h2*pheno_var);
	set_value("sigma_gg", (1-h2)*pheno_var);
	*/
	
	//set_value("omega_gg", 1.5);
	//set_value("omega_ff", 1.5);
	/*
	auto pheno_var_f = 3;
	auto h2 = 0.999;
	set_value("omega_ff", h2*pheno_var_f);
	set_value("sigma_ff", (1-h2)*pheno_var_f);
*/

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
	
	gillespie();

	auto options = "output_dir='cloglog/Output_files' nsample='100000' burnin='25000' thin='10' cloglog='on' deltaT='2'";
	
	output_datatable(file_in,model.output_dir + "/datatable.txt",options,UNSET);
	
	//output_matrix();
	//output_pedigree();
}
	
	

/// This scans through different h2  
void Simulate::scan_h2(string file_in)
{
	//population_structure(file_in); return;
	
	cout << file_in << " file\n";
	
	auto pheno_var = 1;

	auto nsamp = 200000;
	auto burnin = 50000;

	enum TempType { SUS, INF, SUSINF_FIXINF, SUSINF_FIXSUS, SUSINF_SIZE};
	TempType temptype;
	
	if(file_in == "cloglog/sus_template.xml") temptype = SUS;
	else{
		if(file_in == "cloglog/inf_template.xml") temptype = INF;
		else{
			if(file_in == "cloglog/susinf_fixinf_template.xml") temptype = SUSINF_FIXINF;
			else{
				if(file_in == "cloglog/susinf_fixsus_template.xml") temptype = SUSINF_FIXSUS;
				else{
					if(file_in.substr(0,19) == "cloglog/susinf_size") temptype = SUSINF_SIZE;
					else emsg("Template not recognised");
				}
			}
		}
	}
	
	for(auto h2loop = 0.0; h2loop <= 1.001; h2loop += 0.05){
		cout << h2loop << " h2" << endl;
		auto h2 = h2loop; 
		if(h2 == 0) h2 = 0.001;
		if(h2 > 0.999 && h2 < 1.001) h2 = 0.999;
		
		param_value = model.prior_sample();

		set_value("beta", 0.1);
		set_value("m", 7);

		string type;
		switch(temptype){
			case SUS:
				type = "sus";
				set_value("omega_gg", h2*pheno_var);
				set_value("sigma_gg", (1-h2)*pheno_var);
				break;
				
			case INF:
				type = "inf";
				set_value("omega_ff", h2*pheno_var);
				set_value("sigma_ff", (1-h2)*pheno_var);
				break;
				
			case SUSINF_FIXINF:
				type = "susinf_fixinf";
				set_value("omega_gg", h2*pheno_var);
				set_value("sigma_gg", (1-h2)*pheno_var);
				set_value("omega_ff", 0.5*pheno_var);
				set_value("sigma_ff", 0.5*pheno_var);
				break;
				
			case SUSINF_FIXSUS:
				type = "susinf_fixsus";
				set_value("omega_gg", 0.5*pheno_var);
				set_value("sigma_gg", 0.5*pheno_var);
				set_value("omega_ff", h2*pheno_var);
				set_value("sigma_ff", (1-h2)*pheno_var);
				break;
				
			case SUSINF_SIZE:
				type = "susinf_size";
				set_value("omega_gg", 0.5*pheno_var);
				set_value("sigma_gg", 0.5*pheno_var);
				set_value("omega_ff", 0.5*pheno_var);
				set_value("sigma_ff", 0.5*pheno_var);
				break;
		}
		
		gillespie();
	
		
		stringstream label; label << type << "_h2_" << h2; 
		
		for(auto alg = 1; alg <= 4; alg++){
			stringstream ssroot; ssroot << "ALG_" << alg << "_" << label.str(); 
		
			auto outdir = "cloglog/Data_files/";
			stringstream ss; 
	
			ss << "output_dir='" << outdir << ssroot.str() << "' ";
			ss << "nsample= '" << nsamp << "' burnin='" <<burnin << "' ";
			ss << "thin='100' ";
			auto deltaT_check = UNSET;
			switch(alg){
				case 1: break;
				case 2: deltaT_check = 2; break;
				case 3: ss << "cloglog='on' deltaT='2'"; break;
				case 4: ss << "cloglog='on' deltaT='2' geometric_approx='on'"; break;
			}
			output_datatable(file_in,outdir+ssroot.str()+".xml",ss.str(),deltaT_check);
		}
			
		output_simple_datatable("cloglog/Ricardo_data_files/"+label.str()+".csv");
	}
}


/// Finds a given number within a document
double Simulate::find_file_text_num(string file, string st, string st2, unsigned int col, char delimiter)
{
	//cout << file << " file\n";
	ifstream fin(file);
	
	while(!fin.eof()){
		string li;
		getline(fin,li);
		auto spl = split(li,delimiter);
		if(spl.size() > 0){
			for(auto k = 0; k < spl.size(); k++){
				if(spl[k].substr(0,1) == " " ) spl[k] = spl[k].substr(1);
			}
			
			if(spl[0] == st && (spl[1] == st2 || st2 == "")){
				return number(spl[col]);
			}
		}
	}
	cout << "In file '" << file << " could not find " << st << endl;
}
				

/// Gathers together results from all the inferences done
void Simulate::gather_results()
{
	auto mod_type = {"sus","inf","susinf_fixinf","susinf_fixsus"};
	
	for(string type : mod_type){
		cout << type << "type" << endl;
		for(auto alg = 1; alg <=4; alg++){
			stringstream ss; ss << "cloglog/graph_data_ALG_" << alg << "_" << type << ".txt";
			ofstream fout(ss.str());
		
			stringstream ss2; ss2 << "cloglog/graph_data_ALG_" << alg << "_" << type << "_PA.txt";
			ofstream fout2(ss2.str());
		
			stringstream ss3; ss3 << "cloglog/graph_data_ALG_" << alg << "_" << type << "_CPU.txt";
			ofstream fout3(ss3.str());
		
			for(auto h2loop = 0.0; h2loop <= 1.001; h2loop += 0.05){
				auto h2 = h2loop; 
				if(h2 == 0) h2 = 0.001;
				if(h2 > 0.999 && h2 < 1.001) h2 = 0.999;
				
				stringstream sss; sss << "ALG_" << alg << "_" << type << "_h2_" << h2 << ".txt";
				cout << sss.str() << "file\n";
				ifstream fin(sss.str());
				string last;
				while(!fin.eof()){
					string li;
					getline(fin,li);
					if(fin.eof()) break;
					last = li;
				}
	
				stringstream ss; ss << "cloglog/Data_files/ALG_" << alg << "_" << type << "_h2_" << h2 << "/";
				auto dir = ss.str();
				
				fout << h2;
					
				// Loads trace file
				double mean, CImin, CImax;	
				if(type == "sus" || type == "susinf_fixinf" || type == "susinf_fixsus"){
					load_trace_file(dir+"trace_0.txt","omega_gg",mean,CImin,CImax);
					fout << "\t" << mean << "\t" << CImin << "\t" << CImax;
				}
				if(type == "inf" || type == "susinf_fixinf" || type == "susinf_fixsus"){
					load_trace_file(dir+"trace_0.txt","omega_ff",mean,CImin,CImax);
					fout << "\t" << mean << "\t" << CImin << "\t" << CImax;
				}
				fout << endl;
					
				auto i = 0; while(i < last.length()-4 && last.substr(i,4) != "Time") i++;
				if(i == last.length()-4){
				//if(true){
					cout << "Problem with run " << sss.str() << endl;
				}
				else{
					auto algtime = number(last.substr(0,i-1));
					
					fout2 << h2;
					fout3 << h2;
				
					auto ESSmin = LARGE;
					if(type == "sus" || type == "susinf_fixinf" || type == "susinf_fixsus"){
						auto mean = find_file_text_num(dir+"posterior.csv","omega_gg","",1,',');
						auto CImin = find_file_text_num(dir+"posterior.csv","omega_gg","",2,',');
						auto CImax = find_file_text_num(dir+"posterior.csv","omega_gg","",3,',');
					
						auto ESS = find_file_text_num(dir+"posterior.csv","omega_gg","",4,',');
						if(ESS < ESSmin) ESSmin = ESS;
					
						//fout << "\t" << mean << "\t" << CImin << "\t" << CImax;
					}
					
					if(type == "inf" || type == "susinf_fixinf" || type == "susinf_fixsus"){
						auto mean = find_file_text_num(dir+"posterior.csv","omega_ff","",1,',');
						auto CImin = find_file_text_num(dir+"posterior.csv","omega_ff","",2,',');
						auto CImax = find_file_text_num(dir+"posterior.csv","omega_ff","",3,',');
					
						auto ESS = find_file_text_num(dir+"posterior.csv","omega_ff","",4,',');
				
						if(ESS < ESSmin) ESSmin = ESS;
						
						//fout << "\t" << mean << "\t" << CImin << "\t" << CImax;
					}

					if(type == "sus" || type == "susinf_fixinf" || type == "susinf_fixsus"){
						auto pa_sire_g = find_file_text_num(dir+"pred_accs.csv","sire","g_a",2,',');
						auto pa_dam_g = find_file_text_num(dir+"pred_accs.csv","dam","g_a",2,',');
						auto pa_prog_g = find_file_text_num(dir+"pred_accs.csv","prog","g_a",2,',');
						
						fout2 << "\t" << pa_sire_g << "\t" << pa_dam_g << "\t" << pa_prog_g;
					}
										
					if(type == "inf" || type == "susinf_fixinf" || type == "susinf_fixsus"){
						auto pa_sire_f = find_file_text_num(dir+"pred_accs.csv","sire","f_a",2,',');
						//cout << h2 << " " << pa_sire_f << "j\n";
						auto pa_dam_f = find_file_text_num(dir+"pred_accs.csv","dam","f_a",2,',');
						auto pa_prog_f= find_file_text_num(dir+"pred_accs.csv","prog","f_a",2,',');
						
						fout2 << "\t" << pa_sire_f << "\t" << pa_dam_f << "\t" << pa_prog_f;
					}
					
				
					fout3 << "\t" << algtime*200/ESSmin << "\t" << ESSmin;
					// Need to calculate time taken
					
					fout << endl;
					fout2 << endl;
					fout3 << endl;
				}
			}
			//if(alg == 2) break;
		}
	}
}

/// Extracts information from trace file
void Simulate::load_trace_file(string file, string col, double &mean, double &CImin, double &CImax)
{
	ifstream fin(file);
	
	string li;
	getline(fin,li);
	auto spl = split(li,'\t');
	auto c = 0; while(c < spl.size() && spl[c] != col) c++;
	if(c == spl.size()){ cout << file << "   "<< col << "col\n"; emsg("problem");}
		
	vector <double> result;
	while(!fin.eof()){
		string li;
		getline(fin,li);
		auto spl = split(li,'\t');
		if(spl.size() > c) result.push_back(number(spl[c]));
	}
	if(result.size() == 0) return;
	vector <double> vec2;
	for(auto i = result.size()/4; i < result.size(); i++){
		vec2.push_back(result[i]);
	}
	
	mean = 0.0;
	for(auto val :vec2) mean += val;
	mean /= vec2.size();
	cout << mean << " " << vec2.size() << "mea\n";
	auto n = vec2.size();
	sort(vec2.begin(), vec2.end());

	auto i = (unsigned int)((n - 1) * 0.025);
	auto f = (n - 1) * 0.025 - i;
	CImin = vec2[i] * (1 - f) + vec2[i + 1] * f;

	i = (unsigned int)((n - 1) * 0.975);
	f = (n - 1) * 0.975 - i;
	CImax = vec2[i] * (1 - f) + vec2[i + 1] * f;
}


/// Given a set of model parameters this simulates the system using a modified Gillespie algorithm
void Simulate::gillespie()
{
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
				while (j < N && z > S_st[j]) j++;
				
				if (j == N) emsg("Problem here");

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
void Simulate::output_datatable(string file_in, string file_out, string options, double deltaT_check) 
{
	vector <string> line;
	ifstream fin(file_in);
	while(!fin.eof()){
		string li;
		getline(fin,li);
		
		if(li.length() > 0){
			if(li.substr(li.length()-1,1) == "\r") li = li.substr(0,li.length()-1);
		}
		
		if(li.substr(0,6) == "<mcmc "){
			li = "<mcmc "+options+"/>";
		}
		
		if(deltaT_check != UNSET){
			auto len = 13;
			if(li.length() > len){
				int i = 0; while(i < li.length()-len && li.substr(i,len) != "data_column='") i++;
			
				if(i >= 0 && i < li.length()-len){
					//cout << l
					li = li.substr(0,i)+li.substr(i+len+2);
				}
			}
		}
		
		line.push_back(li);
	}
	
	ofstream fout(file_out);

	auto lin = 0; 
	while(lin < line.size()){
		if(line[lin].substr(0,10) == "<datatable") break;
		fout << line[lin] << endl;
		lin++;
	}
		
	fout << "<datatable id='1' sire='2' dam='3' group='4' initial_comp='5'";
	auto col = 6;

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

	if(deltaT_check != UNSET){
		fout << " comp_status='" << col << "'"; col++;
	}
	else{
		col += model.ntrans;
	}
	
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

		fout << "\t" << indiv.sire;
		fout << "\t" << indiv.dam;
		
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

		if(deltaT_check != UNSET){
			fout << "\t";
			if(indiv.initial_comp ==  UNSET) fout << "NA";
			else{
				for(auto t = 0; t < 20; t += deltaT_check){
					auto c = indiv.initial_comp;
					auto tr = 0; while(tr < model.ntrans &&  ind.trans_time[tr] < t){ c = model.trans[tr].to; tr++;}
					if(t != 0) fout << ",";
					fout << "[";
					//if( model.comp[c].name == "S" ||  model.comp[c].name == "I") fout << "S|I";
					//else fout << "R";

					fout << model.comp[c].name;
					fout << ":" << t << "]";
				}
			}
		
	
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
		}
		else{
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

	fout << "</datatable>";

	vector <string> sire_id;
	vector <string> dam_id;
	vector <string> prog_id;
	
	for (auto i = 0; i < model.N; i++) {
		auto id = model.individual[i].id;
		if(i < 1000) sire_id.push_back(id);
		else{
			if(i < 3000) dam_id.push_back(id);
			else prog_id.push_back(id);
		}
	}
	
	for(auto loop = 0; loop < 3; loop++){
		vector <string> list;
		string str;
		switch(loop){
			case 0: str = "sire"; list = sire_id; break;
			case 1: str = "dam"; list = dam_id; break;
			case 2: str = "prog"; list = prog_id; break;
		}
		
		fout << endl;
	
		fout << "<prediction_accuracy name='" << str << "' ind='";
		for(auto i = 0; i < list.size(); i++){
			if(i != 0) fout << ",";
			fout << list[i];
		}
		fout << "'/>" << endl;
	}
	
	fout << endl << "</SIRE>" << endl << endl;
	
	if(false){
		while(lin < line.size()){
			if(line[lin].substr(0,11) == "</datatable") break;
			lin++;
		}
		
		while(lin < line.size()){
			fout << line[lin] << endl;
			lin++;
		}
	}
}

void Simulate::output_simple_datatable(string file)
{
	ofstream fout(file);

	fout << "id,sire,dam,group,donor,It,Rt";
	
	for (auto ie = 0; ie < model.nind_effect; ie++) {
		fout << "," << model.ind_effect[ie].name;
	}
	fout << endl;
	
	for(auto i = 0; i < model.individual.size(); i++){
		auto &indiv = model.individual[i];
		auto &ind = ind_value[i];

		fout << indiv.id << "," << indiv.sire << "," << indiv.dam;
		if(indiv.group == UNSET) fout << ",NA";
		else fout << "," << indiv.group+1;
		
		if(indiv.group == UNSET) fout << ",NA";
		else{
			if(indiv.initial_comp == 0) fout << ",0";
			else fout << ",1";
		}
		
		for (auto tr = 0; tr < model.ntrans; tr++) {
			auto g = indiv.group;
			if (g == UNSET)	fout << ",NA";
			else {
				auto &gr = model.group[g];

				auto t = ind.trans_time[tr];
				if (t != UNSET && t > gr.observation_range.tmin && t < gr.observation_range.tmax){	
					fout << "," << t;
				}
				else fout << ",no";
			}
		}
		
		for (auto ie = 0; ie < model.nind_effect; ie++) {
			fout << "," << ind.ind_effect[ie];
		}

		fout << endl;
	}
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

/// This simulates a new population structure
void Simulate::population_structure(string file_in)
{
	cout << "Population structure simulated!" << endl;
	
	enum TempType { SUS, INF, SUSINF_FIXINF, SUSINF_FIXSUS, SUSINF_SIZE};
	TempType temptype;
	
	if(file_in == "cloglog/sus_template.xml") temptype = SUS;
	else{
		if(file_in == "cloglog/inf_template.xml") temptype = INF;
		else{
			if(file_in == "cloglog/susinf_fixinf_template.xml") temptype = SUSINF_FIXINF;
			else{
				if(file_in == "cloglog/susinf_fixsus_template.xml") temptype = SUSINF_FIXSUS;
				else{
					if(file_in == "cloglog/susinf_size_template.xml") temptype = SUSINF_SIZE;
					else emsg("Template not recognised");
				}
			}
		}
	}
	
	auto nsire = 1000;
	
	auto	ndam = 2;
	auto nprog = 5;
	
	auto num_per_group = 5;
	auto nseed = 1;
	
	struct Ind{
		string id; 
		string sire;
		string dam;
		unsigned int group;
	};
	
	vector <Ind> ind;
	
	auto IDnum = 1;
	
	for(auto i = 0; i < nsire; i++){
		stringstream ss; ss << IDnum;cout << ss.str() << "sire\n"; IDnum++;
		Ind in; in.id = ss.str(); in.sire = "."; in.dam = "."; in.group = -1;
		ind.push_back(in);
	}
	
	for(auto i = 0; i < ndam*nsire; i++){
		stringstream ss; ss << IDnum; IDnum++;cout << ss.str() << "dam\n"; 
		Ind in; in.id = ss.str(); in.sire = "."; in.dam = "."; in.group = -1;
		ind.push_back(in);
	}
	
	auto num_pr = 0;
	for(auto i = 0; i < nsire; i++){
		for(auto d = 0; d < ndam; d++){
			for(auto k = 0; k < nprog; k++){
				stringstream ss; ss << IDnum; IDnum++;
				Ind in; in.id = ss.str(); in.sire = ind[i].id; in.dam = ind[nsire+i*2+d].id;
				cout << ss.str() << " " << ind[i].id << " " << ind[nsire+i*2+d].id << "kk\n";
				in.group = -1;
				ind.push_back(in);
				num_pr++;
			}
		}
	}
	
	// Creates datatable
	ofstream fout("datatab.txt");
	fout << "<datatable id='1' sire='2' dam='3' group='4' initial_comp='5'";
	
	switch(temptype){
		case SUS: fout << " g_a='8' g_e='9'"; break;
		case INF: fout << " f_a='8' f_e='9'"; break;
		case SUSINF_FIXINF: fout << " f_a='8' f_e='9' g_a='10' g_e='11'"; break;
		case SUSINF_FIXSUS: fout << " f_a='8' f_e='9' g_a='10' g_e='11'"; break;
		case SUSINF_SIZE: fout << " f_a='8' f_e='9' g_a='10' g_e='11'"; break;
	}
	cout << temptype << "temptype\n";
	fout << ">" << endl;
	
	for(auto i = 0; i < nsire+ndam*nsire; i++){
		fout << ind[i].id << "\t" << ind[i].sire << "\t" << ind[i].dam <<	"\tNA\tNA\tNA\tNA\t";
		switch(temptype){
			case SUS: case INF:	fout << "0\t0"<< endl; break;
			case SUSINF_FIXINF: case SUSINF_FIXSUS: fout << "0\t0\t0\t0" << endl; break;
			case SUSINF_SIZE: fout << "0\t0\t0\t0" << endl; break;
		}
	}
	
	vector <unsigned int> list;
	for(auto i = 0; i < num_pr; i++){
		list.push_back(nsire+ndam*nsire+i);
	}
	
	auto Z = num_pr/num_per_group;
	for(auto z = 0; z < Z; z++){
		stringstream ss; ss << z+1;
		for(auto j = 0; j < num_per_group; j++){
			auto init = "S";
			if(j == 0) init = "I";
			
			auto k = int(ran()*list.size());
			auto i = list[k];
			
			list.erase(list.begin()+k);
			fout << ind[i].id << "\t" << ind[i].sire << "\t" << ind[i].dam << "\t" 
			      << ss.str() << "\t" << init << "\t";
			if(init == "I") fout << "no\t" << ran() << "\t";
			else fout << ran() << "\t" << 1+ ran() << "\t";
				
			switch(temptype){
			case SUS: case INF:	fout << "0\t0"<< endl; break;
			case SUSINF_FIXINF: case SUSINF_FIXSUS: fout << "0\t0\t0\t0" << endl; break;
			case SUSINF_SIZE: fout << "0\t0\t0\t0" << endl; break;
			}
		}
	}
	fout << "</datatable>" << endl;
	
	if(list.size() != 0) emsg("List not zero!");
	
	cout << Z << " number of groups" << endl;
	cout << ind.size() << " number of individuals" << endl;
}

