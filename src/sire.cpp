// SIRE 2.0 stands for “Susceptibility, Infectivity and Recoverability Estimation”.


// This is a software tool which allows simultaneous estimation of the genetic
// effect of a single nucleotide polymorphism (SNP), as well as additive
// genetic, environmental and non-genetic influence on host susceptibility,
// infectivity and recoverability.

// This work is described in:

// "Estimating individuals’ genetic and non-genetic effects underlying
// infectious disease transmission from temporal epidemic data"

// Christopher M. Pooley 1,2*, Glenn Marion 2&, Stephen C. Bishop&, Richard I.
// Bailey, and Andrea B. Doeschl-Wilson 1& 1 The Roslin Institute, The
// University of Edinburgh, Midlothian, EH25 9RG, UK.  2 Biomathematics and
// Statistics Scotland, James Clerk Maxwell Building, The King's Buildings,
// Peter Guthrie Tait Road, Edinburgh, EH9 3FD, UK & Deceased


// This code takes an init.xml input file and performs MCMC analysis.
// This produces posterior samples for both paramters and event
// sequences which are stored in the output directory.


// Load mpi: module load mpi/openmpi-x86_64

// Compile using:    make

// Run using: mpirun -n 20 ./sire examples/example1.xml 0
// Note, -n denotes the number of cores (and chains) and 0 is the random seed used to generate MCMC. This can be changed to any other integer.

// mpirun -n 20 ./sire Fishboost/scen-3-1.xml 0

// Simulate using: ./sire examples/example1.xml 0 0
// Simulate using: ./sire cloglog/sus_template.xml 0 0
// Simulate using: ./sire cloglog/susinf_size_template.xml 0 0
// Simulate using: ./sire cloglog/susinf_size10_template.xml 0 0

// Inference using: nohup ./sire cloglog/Data_files/ALG_1_sus_h2_0.001.xml 0 > alg1
// ./sire cloglog/Data_files/test.xml 0

// Inference using: nohup ./sire cloglog/Data_files/ALG_2_susinf_size10.xml 0 > size10&


#include <iostream>
#include <sstream>
#include <math.h>
#include <map>
#include <algorithm>
#include <vector>
#include <iterator>
#include <fstream>
#include <iomanip>
#include <string>
#include <mpi.h>

using namespace std;

#include "model.hpp"
#include "mcmc.hpp"
#include "map.hpp"
#include "emulator.hpp"
#include "pas.hpp"
#include "simulate.hpp"
#include "utils.hpp"
#include "show_progress.hpp"
#include "timers.hpp"

int main(int argc, char *argv[]) 
{
	if (argc < 3 || argc > 4) {
		int ncore;
		MPI_Comm_size(MPI_COMM_WORLD,&ncore); 
		if(ncore == 0) cout << "You must include an input xml file name and a seed number." << endl;
		return 0;
	}

  MPI_Init(&argc,&argv);
	
	timers_init();                              // Initialises timers

	auto seed = atoi(argv[2]);
	set_seed(seed);                             // Sets the random seed
	srand(seed);                     

	auto sim = false;
	if (argc == 4) sim = true;

	string file(argv[1]);
	Model model(file,sim);                          // This loads up the model from the input XML file

	Simulate simulate(model);
	//simulate.gather_results(); return 0;
	
	if (sim == true) {
		simulate.run(argv[1]);    // This simulates from the model (used only for testing)
		//simulate.scan_h2(argv[1]);
		return 0;
	}

	timer[TIME_TOTAL].start();
	
	switch(model.algorithm){
		case ALG_MAP:
			{
				MAP map(model);
				//map.cmaes(seed);
				//map.gradient_descent(); 
				map.REML();
			}
			break;
			
		case ALG_EMULATOR:
			{
				Emulator em(model,seed);
				em.inference();
			}
			break;
			
		case ALG_MCMC:
			{
				MCMC mcmc(model);                           // Initialises an mcmc chain based on the model

				//simulate.set_initial_state(mcmc);         // Sets initial chain state to simulation

				mcmc.set_likelihoods();                     // Sets the likelihoods at the start of the chain

				for(auto s = 0; s < model.nsample; s++) {   // Iterates over MCMC samples
					show_progress(s, model.nsample, 100);
					//if (s % 1000 == 0) cout << "Sample: " << s << " / " << model.nsample << " " << mcmc.quench.phi_L <<  endl;

					if (s < model.nburnin)
						mcmc.burnin = true;                     // Determines if burnin or not
					else
						mcmc.burnin = false;

					mcmc.set_quench(s);                       // Sets inverse temperatures for quenching 
					
					mcmc.update();                            // Performs MCMC updates

					if (mcmc.burnin == true)                  // Update the infection sampler (used to add and remove infected individuals)
						mcmc.initialise_inf_sampler(s);

					if (s % model.nthin == 0) {
						mcmc.trace_output(s);                   // Outputs to the trace plot (every nthin steps)

						if(mcmc.burnin == false) mcmc.store_sample(); // Stores a parameter sample (for statistical analysis later)
					}
					
					if(mcmc.burnin == true) mcmc.store_burnin_Li();

					if (mcmc.burnin == false)                 // Updates posterior average of individual effects (to calcualte PA later)
						mcmc.ind_effect_posterior_mean_update();
					
					//if(s == model.nburnin) mcmc.optimum_quench_schedule();
				}
			
				auto dir = model.output_dir;

				cout << endl << "Outputs are placed in the '" << dir << "' directory" << endl;

				mcmc.output_statistics();  // Outputs diagnostic information
				
				mcmc.check_quenching();
			}
			break;
			
		case ALG_PAS:
			PAS pas(model);
			pas.run();
			break;
	}
	
	timer[TIME_TOTAL].stop();

	output_timers(model.output_dir + "/CPU_timings.txt");    // Outputs information about CPU timings
		
	MPI_Finalize();
}
