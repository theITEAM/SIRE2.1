/// Stores the CPU clock times for different parts of the algorithm (for diagnostic purposes)

#include <fstream>
#include <sstream>

using namespace std;

#include "timers.hpp"
#include "const.hpp"
#include "mpi.hpp"

vector <Timer> timer;

void Timer::start() {
	val -= clock();
}

void Timer::stop() {
	val += clock();
}

void timers_init() {
	timer.resize(TIMERMAX);

	for (auto &ti : timer)
		ti.val = 0;
}


/// Outputs CPU timing information to a file
void output_timers(string file)
{
	Mpi mpi;

	vector <double> timer_tot(TIMERMAX);
	
	for(auto sel = 0; sel < TIMERMAX; sel++){
		auto tot = mpi.gather(timer[sel].val);
		if(mpi.core == 0){
			timer_tot[sel] = 0;
			for(auto co = 0; co < mpi.ncore; co++) timer_tot[sel] += tot[co];
		}
	}
	
	if(mpi.core == 0){	
		ofstream dia(file);
		if (!dia)
			emsg("Cannot open the file '" + file + "'");

		dia.precision(3);

		dia << "Timings for different parts of the algorithm:" << endl << endl;
		double total = timer_tot[TIME_TOTAL];

		for(auto sel = 0; sel < TIMERMAX; sel++){
			string text;
			switch(sel){
				case TIME_TOTAL: break;
				case TIME_UPDATE: text = "Overall MCMC update"; break; 
				case TIME_anneal: text = "anneal"; break; 
				case TIME_COVAR: text = "Covariate"; break; 
				case TIME_COVAR_INIT: text = "Covariate initialise"; break; 
				case TIME_TEMP: text = "Temp"; break; 
				case TIME_TRANS_PARAM_INIT: text = "Transition parameters initialise"; break; 
				case TIME_TRANS_PARAM: text = "Transition parameters"; break; 
				case TIME_GROUP_EFFECT_INIT: text = "Group effect initialise"; break; 
				case TIME_GROUP_EFFECT: text = "Group effect"; break; 
				case TIME_GROUP_SD: text = "Group standard deviation"; break; 
				case TIME_FIXED_EFFECTS: text = "Fixed effects"; break; 
				case TIME_SNP_EFFECTS: text = "SNP effects"; break; 
				case TIME_SUS_IND_EFFECT_INIT: text = "Susceptibility individual effect initialsie"; break; 
				case TIME_SUS_IND_EFFECT: text = "Susceptibility individual effect"; break; 
				case TIME_TRANS_IND_EFFECT: text = "Transition individual effect"; break; 
				case TIME_INF_IND_EFFECT_INIT: text = "Infectivity individual effect initialise"; break; 
				case TIME_INF_IND_EFFECT: text = "Infectivity individual effect"; break; 
				case TIME_JOINT_IEVAR: text = "Joint individual effect and variance updates"; break; 
				case TIME_TRANS_RATE: text = "Transition rate"; break; 
				case TIME_EVENTS: text = "Events"; break; 
				case TIME_EVENTS_LIKE: text = "Events likelihood"; break; 
				case TIME_MEAN_EVENTS: text = "Mean Events"; break; 
				case TIME_CHECK: text = "Checking"; break; 
				case TIME_ADD_REM: text = "Add and removal"; break; 
				case TIME_ADD_REM_LIKE: text = "Add and removal likelihood"; break; 
				case TIME_SE_SP: text = "Se and Sp"; break; 
				case TIME_INF_SAMP_UPDATE: text = "Infection sample update"; break; 
				case TIME_EMULATE_NEW_POINT: text = "Emulate new points"; break; 
				case TIME_ADD_DATA_SAMP: text = "Add data samples"; break; 
				case TIME_TUNE: text = "Tune emulator"; break; 
				case TIME_EMULATOR_MCMC: text = "MCMC on emulator"; break; 
				case TIME_EMULATOR_MCMC_SAMPLE: text = "MCMC on emulator sample"; break; 
				case TIME_GRAD_H1: text = "Calcuate gradient and H part 1"; break;  
				case TIME_GRAD_H2: text = "Calculate gradient and H part 2"; break;  
				case TIME_GRAD_H3: text = "Calculate gradient and H part 3"; break;
				case TIME_MAXIMISE: text = "Maximise likelihood"; break;
				case TIME_DET: text = "Calculate determinant"; break;		
			}
			
			if(text != ""){
				auto per = 100 * timer_tot[sel] / total;
				if(per > 0) dia << int(per) << "% " << text << endl;
			}
		}
	
		cout << timer_tot[TIME_TOTAL] / (60.0 * CLOCKS_PER_SEC) << " Time in minutes" << endl;
	}
}
