/// Stores the CPU clock times for different parts of the algorithm (for diagnostic purposes)

#include <fstream>
#include <sstream>

using namespace std;

#include "timers.hh"
#include "const.hh"

vector <Timer> timer;

void Timer::start()
{
  val -= clock();
}

void Timer::stop()
{
  val += clock();
}

void timers_init()
{
  timer.resize(TIMERMAX);
  
  for(auto &ti : timer) ti.val = 0;
}


/// Outputs CPU timing information to a file
void output_timers(string file)
{
  ofstream dia(file); if(!dia) emsg("Cannot open the file '"+file+"'");

  dia.precision(3);
    
  dia << "Timings for different parts of the algorithm:" << endl << endl;
  double total = timer[TIME_TOTAL].val;
  
  dia << int(100*timer[TIME_COVAR_INIT].val/total) << "% Covariance proposals initialise" << endl;
  dia << int(100*timer[TIME_COVAR].val/total) << "% Covariance proposals" << endl;

  dia << int(100*timer[TIME_TRANS_RATE].val/total) << "% Transition rate proposals" << endl;

  dia << int(100*timer[TIME_TRANS_PARAM_INIT].val/total) << "% Trans parameter initialise" << endl;
  dia << int(100*timer[TIME_TRANS_PARAM].val/total) << "% Trans param proposals" << endl;
  
  dia << int(100*timer[TIME_GROUP_EFFECT_INIT].val/total) << "% Group effect initialise" << endl;
  dia << int(100*timer[TIME_GROUP_EFFECT].val/total) << "% Group effect proposals" << endl;
  dia << int(100*timer[TIME_GROUP_SD].val/total) << "% Group effect sd proposals" << endl;
    
  dia << int(100*timer[TIME_FIXED_EFFECTS].val/total) << "% Fixed effects proposals" << endl;
  
	dia << int(100*timer[TIME_SNP_EFFECTS].val/total) << "% SNP effects proposals" << endl;
 
	dia << int(100*timer[TIME_SE_SP].val/total) << "% sensitivty and specificity proposals" << endl;
  
  dia << int(100*timer[TIME_SUS_IND_EFFECT_INIT].val/total) << "% Susceptible ind effect initialise" << endl;
  dia << int(100*timer[TIME_SUS_IND_EFFECT].val/total) << "% Susceptible ind effect proposals" << endl;

  dia << int(100*timer[TIME_INF_IND_EFFECT_INIT].val/total) << "% Infectivity ind effect initialise" << endl;
  dia << int(100*timer[TIME_INF_IND_EFFECT].val/total) << "% Infectivity ind effect proposals" << endl;

  dia << int(100*timer[TIME_TRANS_IND_EFFECT].val/total) << "% Trans ind effect proposals" << endl;

  dia << int(100*timer[TIME_EVENTS].val/total) << "% Events proposals" << endl;
  dia << int(100*timer[TIME_EVENTS_LIKE].val/total) << "% Events proposals likelihood" << endl;

  dia << int(100*timer[TIME_ADD_REM].val/total) << "% Add/remove proposals" << endl;
  dia << int(100*timer[TIME_ADD_REM_LIKE].val/total) << "% Add/remove likelihood" << endl;

	dia << int(100*timer[TIME_INF_SAMP_UPDATE].val/total) << "% Infection sampler update" << endl;
	
  dia << int(100*timer[TIME_CHECK].val/total) << "% Checking" << endl;
  
  //dia << int(100*timer[TIME_TEMP].val/total) << "% Temp" << endl;
  
  dia << endl;
  
  cout << timer[TIME_TOTAL].val/(60.0*CLOCKS_PER_SEC) << " Time in minutes" << endl;
}
