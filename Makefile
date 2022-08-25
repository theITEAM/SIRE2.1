CXX = mpicxx
CXXFLAGS = -O3 -std=c++11 -march=native -MMD -MP   -fmax-errors=3 -g
#CXXLFAGS += -Wall -Wextra 
LIBRARIES = 
SOURCES = src/check.cpp src/emulator.cpp src/emulator_check.cpp src/matrix.cpp src/mcmc.cpp src/model_diagnostic_test.cpp src/map.cpp src/model_fixed_effects.cpp src/model_funcs.cpp src/model_group_effect.cpp src/model_ind_effects.cpp src/model_ind_effects_proposals.cpp src/model_inf_events.cpp src/model_mean_event_proposals.cpp src/model_initialise.cpp src/model_likelihood.cpp src/model_matrix.cpp src/model_prior.cpp src/model_snp_effects.cpp src/model_trans_events.cpp src/mpi.cpp src/show_progress.cpp src/simulate.cpp src/sire.cpp src/pas.cpp src/timers.cpp src/tinyxml2.cpp src/utils.cpp

OBJECTS = $(SOURCES:.cpp=.o)

.PHONY: clean all
.DEFAULT_GOAL := all

all: sire

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(LIBRARIES) $< -o $@ -c

sire: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LIBRARIES) $^ -o $@
	cp sire ..

clean:
	rm *.o sire *.d

-include $(OBJECTS:.o=.d)

