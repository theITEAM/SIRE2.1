#ifndef SIRE__CONST_HH
#define SIRE__CONST_HH

const int UNSET = 999999;                            // Value used to represent an unset quantity
const double LARGE = 99999999;                       // A token large amount
const double TINY = 0.00000000001;                   // Used to represent a tiny number
const double SMALL = 0.000001;                       // Used to represent a small number
const double ZERO_PRIOR = -LARGE;                    // Sets an effective zero prior probability
const double EXT_FOI = 0.0000001;                    // A token external force of infection

const int INITIAL_TRANS_SAMPLE = 1000;               // The number of times to initial try to sample a transitiom

const auto covar_prop = 3u;                          // Number of covariance proposals
const auto trans_param_prop = 10u;                   // Number of transition parameter proposals
const auto group_effect_sigma_prop = 10u;            // Number of group effect sigma proposals
const auto group_effect_prop = 10u;                  // Number of group effect proposals
const auto fixed_effect_prop = 1u;                   // Number of fixed effect proposals
const auto snp_effect_prop = 1u;                     // Number of snp effect proposals
const auto ind_effect_sus_prop = 1u;                 // Number of individual effect proposals (susceptibility)
const auto ind_effect_inf_prop = 1u;                 // Number of individual effect proposals (infectivity)
const auto ind_effect_inf_prop2 = 1u;                // Number of individual effect proposals (infectivity multiple)
const auto ind_effect_trans_prop = 1u;               // Number of individual effect proposals (transitions)
const auto ind_effect_trans_prop2 = 1u;              // Number of individual effect proposals (transitions multiple)
const auto se_sp_prop = 1u;                          // Number of Se Sp proposals
const auto nbin = 50;                                // The number of bins used for the infection samplier 
                                          
const auto prop_up = 1.01;                           // Dynamically alters proposal jump size
const auto prop_down = 0.995;
	
enum PriorType { FLAT_PRIOR, NORMAL_FROM_SD_PRIOR};  // Different types of prior

enum ParamType { TRANS_INFRATE, TRANS_MEAN,          // Different types of parameter
                 TRANS_SHAPE, COVAR_MATRIX, FIXED_EFFECT,
								 SNP_EFFECT_MAG, SNP_EFFECT_DOM, GROUP_SD, GROUP_EFFECT, DIAG_TEST, NO_PARAMTYPE};

enum TransType { INFECTION, GAMMA};                  // Different types of transition

enum SNP { AA, AB, BB};                              // Potential values a SNP can take

enum Status { INFECTED, NOT_INFECTED, UNKNOWN};      // Defines individual disease status

enum InfModel { FREQ_DEP, DENSITY_DEP};              // Type of model used for infection process

enum EvInitSampType { SAMPLE_ONLY, SAMPLE_MIX, SAMPLE_RANGE};// Different sampling for initial state

enum EvSampType { PEAKED_SAMPLE, UNIFORM_SAMPLE};    // Type of sampling performed when changing event times

enum Timers { TIME_TOTAL, TIME_COVAR, TIME_COVAR_INIT,
					TIME_TEMP, TIME_TRANS_PARAM_INIT, TIME_TRANS_PARAM, 
					TIME_GROUP_EFFECT_INIT, TIME_GROUP_EFFECT, TIME_GROUP_SD,
					TIME_FIXED_EFFECTS, TIME_SNP_EFFECTS, 
					TIME_SUS_IND_EFFECT_INIT, TIME_SUS_IND_EFFECT,
					TIME_TRANS_IND_EFFECT,
					TIME_INF_IND_EFFECT_INIT, TIME_INF_IND_EFFECT,
					TIME_TRANS_RATE, TIME_EVENTS, TIME_EVENTS_LIKE,
					TIME_CHECK,
					TIME_ADD_REM, TIME_ADD_REM_LIKE,
					TIME_SE_SP,
					TIME_INF_SAMP_UPDATE,
					TIMERMAX};
#include "struct.hh"

#endif