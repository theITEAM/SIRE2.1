/// Particle annealed sampling

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
const auto ie_fac = 1.0;

using namespace std;

#include "pas.hpp"
#include "show_progress.hpp"
#include "timers.hpp"

PAS::PAS(const Model &model) : model(model), mcmc(model)
{
	mcmc.set_likelihoods();  
	//simulate.set_initial_state(mcmc);         // Sets initial chain state to simulation
	
	phi = 0;

	mcmc.burnin = true;

	mcmc.quench.phi_L = phi;
	mcmc.quench.phi_DT = phi;
	mcmc.quench.phi_IE = ie_fac + (1-ie_fac)*phi;
	mcmc.quench.phi_Pr = 1;
}

void PAS::run()
{
	if(mpi.core == 0) initialise_gen_plot();
	
	timer[TIME_QUENCH].start();
		
	auto g = 0;
	do{
		if(mpi.core == 0) cout << "Generation: " << g << "  phi " << phi << endl;
		
		mcmc.L_samp.clear();
		gen_sample.clear();
		
		for(auto s = 0; s < model.nsample_per_gen; s++){   // Iterates over MCMC samples
			mcmc.update();
			
			auto L = mcmc.L_trans_events;
			for(auto val : mcmc.L_inf_events) L += val;
			mcmc.L_samp.push_back(L);
			
			Sample samp; samp.param_value = mcmc.param_value;
			gen_sample.push_back(samp);
		}
		gen_plot(g);

		if(phi == model.phi_final) break;
		
		bootstrap();
		g++;
	}while(true);

	timer[TIME_QUENCH].stop();
	
	for(auto s = 0; s < model.nsample; s++) {
		if(mpi.core == 0) show_progress(s, model.nsample, 100);
		if (s < model.nburnin) mcmc.burnin = true; // Determines if burnin or not
		else mcmc.burnin = false;

		mcmc.update();                            // Performs MCMC updates

		if (mcmc.burnin == true)                  // Update the infection sampler (used to add and remove infected individuals)
			mcmc.initialise_inf_sampler(s);

		if (s % model.nthin == 0) {
			mcmc.trace_output(s);                   // Outputs to the trace plot (every nthin steps)

			if(mcmc.burnin == false) mcmc.store_sample(); // Stores a parameter sample (for statistical analysis later)
		}
		
	
		if (mcmc.burnin == false)                 // Updates posterior average of individual effects (to calcualte PA later)
			mcmc.ind_effect_posterior_mean_update();
	}

	auto dir = model.output_dir;

	if(mpi.core == 0) cout << endl << "Outputs are placed in the '" << dir << "' directory" << endl;

	mcmc.output_statistics();  // Outputs diagnostic information
}


/// Initialises a file which 
void PAS::initialise_gen_plot()
{
	genout.open(model.output_dir +"/generations.csv");
	genout << "Generation, phi";
	for(auto th = 0; th < model.nparam; th++){
		genout << ", " << model.param[th].name << ", 95% CI min, 95% CI max";
	}
	genout << endl;
}


/// Initialises a file which 
void PAS::gen_plot(unsigned int g)
{
	auto chain_sample = mpi.gather_sample(gen_sample);

	if(mpi.core == 0){
		genout << g << ", " << phi;
		for(auto th = 0; th < model.nparam; th++){
			vector <double> vec;
			for(auto c = 0; c < mpi.ncore; c++){
				for(auto i = 0; i < gen_sample.size(); i++){
					vec.push_back(chain_sample[c][i].param_value[th]);
				}
			}
			auto stat = mcmc.get_statistic(vec);

			genout << ", " << stat.mean << ", " << stat.CImin << ", "  << stat.CImax;
		}
		genout << endl;
	}
}

		
/// Performs a bootstrap step
void PAS::bootstrap()
{
	auto av = 0.0;
	for(auto val : mcmc.L_samp) av += val;
	av /= mcmc.L_samp.size();
	
	auto var = 0.0, nav = 0.0;
	for(auto val : mcmc.L_samp){
		if(val > av){ var += (val-av)*(val-av); nav++;}
	}
	var /= nav;
	
	auto dph = 1.0/sqrt(var);
	//cout << mpi.core << " " << dph << "dph\n";
	
	auto dph_tot = mpi.gather(dph);

	auto dph_av = 0.0; 
	if(mpi.core == 0){
		for(auto val : dph_tot) dph_av += val;
		dph_av /= dph_tot.size();
		if(phi + dph_av > model.phi_final) dph_av = model.phi_final-phi;
	}
	mpi.bcast(dph_av);
	
	auto L = mcmc.L_samp[mcmc.L_samp.size()-1];

	auto L_tot = mpi.gather(L);
	
	vector <long> map(mpi.ncore);
	
	if(mpi.core == 0){
		auto Lav = 0.0; for(auto val : L_tot) Lav += val;
		Lav /= L_tot.size();
		
		//for(auto val : L_tot) cout << val << " L\n";
		
		vector <double> prob_sum(mpi.ncore);
		vector <int> num(mpi.ncore,0);
		
		auto sum = 0.0;
		for(auto c = 0; c < mpi.ncore; c++){
			auto pr = exp((L_tot[c]-Lav)*dph_av);			
			sum += pr;
			prob_sum[c] = sum;
		}
		
		for(auto c = 0; c < mpi.ncore; c++){
			auto z = ran()*sum;
			auto cc = 0; while(cc < mpi.ncore && z > prob_sum[cc]) cc++;
			if(cc == mpi.ncore) emsg("Problem with bootstrap");
		
			num[cc]++;
		}
		
		vector <unsigned int> list;
		for(auto c = 0; c < mpi.ncore; c++){
			if(num[c] > 1){
				for(auto i = 0; i < num[c]-1; i++) list.push_back(c);
			}
		}		
		
		for(auto c = 0; c < mpi.ncore; c++){
			if(num[c] == 0){
				map[c] = list[list.size()-1];
				list.pop_back();
			}
			else map[c] = c;
		}
	}
	mpi.bcast(map);
		
	//if(mpi.core == 0){
		//cout << "Bootstrap " <<  mpi.core << ":  "; for(auto val : map) cout << val << "  "; cout << endl;
//	}
	
	mpi.mcmc_boostrap(mcmc.param_value, mcmc.ind_value, mcmc.L_ind_effect, mcmc.L_inf_events, mcmc.L_trans_events, mcmc.L_diag_test, mcmc.prior, map);
	
	mcmc.check_chain(0);
	
	phi += dph_av;
	mcmc.quench.phi_L = phi;
	mcmc.quench.phi_DT = phi;
	mcmc.quench.phi_IE = ie_fac + (1-ie_fac)*phi;;
}

