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

// Compile using:    make
// Run using: ./sire examples/example1.xml 0 

// Note, 0 is the random seed used to generate MCMC. This can be changed to any other integer.

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

using namespace std;

#include "model.hh"
#include "mcmc.hh"
#include "simulate.hh"
#include "utils.hh"
#include "timers.hh"

int main(int argc, char *argv[]) 
{
  if(argc != 3 && argc != 4) {
    cout << "You must include an input xml file name and a seed number." << endl;
    return 0;
  }
 
  timers_init();                              // Initialises timers
  
  set_seed(atoi(argv[2]));                    // Sets the random seed
  srand(atoi(argv[2]));                       // Sets the random seed
  
  string file(argv[1]);
  Model model(file);                          // This loads up the model from the input XML file
  
  auto sim = false; if(argc == 4) sim = true;
  
  Simulate simulate(model);
  if(sim == true){ simulate.run(); return 0;} // This simulates from the model (used only for testing)
  
  MCMC mcmc(model);                           // Initialises an mcmc chain based on the model
   
  //simulate.set_initial_state(mcmc);         // Sets initial chain state to simulation
  
  mcmc.set_likelihoods();                     // Sets the likelihoods at the start of the chain
  
  timer[TIME_TOTAL].start();
  for(auto s = 0u; s < model.nsample; s++){   // Iterates over MCMC samples
    if(s%1000 == 0) cout << "Sample: " << s << " / " << model.nsample << endl;
    
    if(s < model.nburnin) mcmc.burnin = true; // Determines if burnin or not 
    else mcmc.burnin = false;
    
    mcmc.update();                            // Performs MCMC updates
    
		if(mcmc.burnin == true){                  // Update the infection sampler (used to add and remove infected individuals)
			mcmc.initialise_inf_sampler(s);
		}
		
    if(s%model.nthin == 0){                  
      mcmc.trace_output(s);                   // Outputs to the trace plot (every nthin steps)
    
      if(mcmc.burnin == true){               
        mcmc.store_sample();                  // Stores a parameter sample (for statistical analysis later)
      }       
    }
    
    if(mcmc.burnin == true){                  // Updates posterior average of individual effects (to calcualte PA later)
      mcmc.ind_effect_posterior_mean_update();
    }   
  }
  timer[TIME_TOTAL].stop();
	
	auto dir = model.output_dir;
	
	cout << endl << "Outputs are placed in the '" << dir << "' directory" << endl;
		
  mcmc.output_statistics("diagnostics.txt");  // Outputs diagnostic information
  
  output_timers(dir+"/CPU_timings.txt");      // Outputs information about CPU timings
}
