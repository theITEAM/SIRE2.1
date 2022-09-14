#pragma once

const vector <double> phi_ch = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.82, 0.85, 0.87, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99, 1.0};

const int UNSET = 1e6 - 1;                           // Value used to represent an unset quantity
const double LARGE = 1e8 - 1;                        // A token large amount
const double TINY = 1e-11;                           // Used to represent a tiny number
const double SMALL = 1e-8;                           // Used to represent a small number
const double SMALLISH = 1e-3;                        // Used to represent a smallish number
const double ZERO_PRIOR = -LARGE;                    // Sets an effective zero prior probability
const double EXT_FOI = 1e-7;                         // A token external force of infection

const auto MCMC_PARA = true;

const int INITIAL_TRANS_SAMPLE = 1000;               // The number of times to initially try to sample a transitiom
const int INITIAL_PARAM_SAMPLE = 1000;               // The number of tines to initially sample from prior

const auto covar_prop = 3u;                          // Number of covariance proposals
const auto trans_param_prop = 10;                   // Number of transition parameter proposals
const auto group_effect_sigma_prop = 10;            // Number of group effect sigma proposals
const auto group_effect_prop = 10;                  // Number of group effect proposals
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

const unsigned int ML_GENERATION_TERM_COND = 10;     // The number of generation used in termination (CMAES)

enum PriorType { FIXED_PRIOR, FLAT_PRIOR, NORMAL_FROM_SD_PRIOR};  // Different types of prior

enum ParamType { TRANS_INFRATE, TRANS_MEAN,          // Different types of parameter
                 TRANS_SHAPE, COVAR_MATRIX, FIXED_EFFECT,
                 SNP_EFFECT_MAG, SNP_EFFECT_DOM, GROUP_SD, GROUP_EFFECT, DIAG_TEST, NO_PARAMTYPE
               };

enum TransType { INFECTION, GAMMA, EXP};                  // Different types of transition

enum AlgorithmType { ALG_UNSET, ALG_MCMC, ALG_MAP, ALG_EMULATOR, ALG_PAS};  // Different inference algorithms
	
enum MaximisationType { NEWTONS_METHOD, GRADIENT_ASCENT};

enum VariableType { PARAM, IND_EFFECT};              // Variable type (CMA-ES)
	
enum SNP { AA, AB, BB};                              // Potential values a SNP can take

enum Status { INFECTED, NOT_INFECTED, UNKNOWN};      // Defines individual disease status

enum InfModel { FREQ_DEP, DENSITY_DEP};              // Type of model used for infection process

enum EvInitSampType { SAMPLE_ONLY, SAMPLE_MIX, SAMPLE_RANGE};// Different sampling for initial state

enum EvSampType { PEAKED_SAMPLE, UNIFORM_SAMPLE};    // Type of sampling performed when changing event times

enum MatrixType { DIAG, FULL};

enum Timers { TIME_TOTAL,
							TIME_UPDATE,
							TIME_QUENCH,
              TIME_TEMP, TIME_TRANS_PARAM_INIT, TIME_TRANS_PARAM,
              TIME_GROUP_EFFECT_INIT, TIME_GROUP_EFFECT, TIME_GROUP_SD,
              TIME_FIXED_EFFECTS, TIME_SNP_EFFECTS,
              TIME_SUS_IND_EFFECT_INIT, TIME_SUS_IND_EFFECT,
							TIME_COVAR, TIME_COVAR_INIT,
							TIME_JOINT_IEVAR,
              TIME_TRANS_IND_EFFECT,
              TIME_INF_IND_EFFECT_INIT, TIME_INF_IND_EFFECT,
              TIME_TRANS_RATE, TIME_EVENTS, TIME_EVENTS_LIKE, TIME_MEAN_EVENTS, 
              TIME_CHECK,
              TIME_ADD_REM, TIME_ADD_REM_LIKE,
              TIME_SE_SP,
              TIME_INF_SAMP_UPDATE,
							TIME_EMULATE_NEW_POINT,
							TIME_ADD_DATA_SAMP,
							TIME_TUNE,
							TIME_EMULATOR_MCMC,
							TIME_EMULATOR_MCMC_SAMPLE,
							TIME_GRAD_H1, TIME_GRAD_H2, TIME_GRAD_H3, 
							TIME_MAXIMISE, TIME_DET,
              TIMERMAX
            };

#include "struct.hpp"
#include <vector>
