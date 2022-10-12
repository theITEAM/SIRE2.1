#pragma once

using namespace std;

#include "model.hpp"
#include "mcmc.hpp"

class Simulate {
public:
	vector <double> param_value;                   // The model parameter values

	vector <IndValue> ind_value;                   // Stores individual properties

	Simulate(Model &model);
	void run(string file);
	void scan_h2(string file);
	void set_initial_state(MCMC &mcmc) const;

private:
	void add_infected(const double t_inf, const int i, const int comp_init);
	void set_value(string name, double val);
	void output_datatable(string file_in, string file_out, string options, double deltaT_check);
	void output_simple_datatable(string file);
	void output_matrix() const;
	void output_pedigree();
	void gillespie();

public:  // Used for generating the results in the paper
	void population_structure(string file_in);                    
	void gather_results();   
	double find_file_text_num(string file, string st, string st2, unsigned int col, char delimiter);
	void load_trace_file(string file, string col, double &mean, double &CImin, double &CImax);
	
	Model &model;
};

