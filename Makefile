CXX := g++
CXXFLAGS :=  -std=c++11 -fmax-errors=3 -O3
exe := sire

sire: src/sire.cc src/model_ind_effects_proposals.cc src/model_initialise.cc src/model_prior.cc src/model_funcs.cc src/model_matrix.cc src/mcmc.cc src/utils.cc src/check.cc src/model_ind_effects.cc src/model_inf_events.cc src/model_trans_events.cc src/model_fixed_effects.cc  src/model_snp_effects.cc src/model_group_effect.cc src/model_diagnostic_test.cc src/simulate.cc src/timers.cc src/tinyxml2.cpp
	$(CXX) -o sire src/sire.cc src/model_ind_effects_proposals.cc src/model_initialise.cc src/model_prior.cc src/model_funcs.cc src/model_matrix.cc src/mcmc.cc src/utils.cc src/check.cc src/model_ind_effects.cc src/model_inf_events.cc src/model_trans_events.cc src/model_fixed_effects.cc src/model_snp_effects.cc src/model_group_effect.cc src/model_diagnostic_test.cc src/simulate.cc src/timers.cc  src/tinyxml2.cpp -std=c++11 -fmax-errors=3 -O3 -g