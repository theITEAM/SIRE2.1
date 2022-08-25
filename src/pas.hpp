#pragma once

using namespace std;

#include "model.hpp"
#include "mcmc.hpp"
#include "mpi.hpp"
#include "utils.hpp"

class PAS {
	public:
		PAS(const Model &model);
		void run();
	
	private:
		double phi;
	
		vector <Sample> gen_sample; 
		
		ofstream genout;
	
		void initialise_gen_plot();
		void gen_plot(unsigned int g);
		void bootstrap();
	
		MCMC mcmc;		 
		const Model &model;
		Mpi mpi;
};
