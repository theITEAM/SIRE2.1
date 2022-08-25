#pragma once

using namespace std;

#include "model.hpp"
#include "mcmc.hpp"

class Simulate {
public:
	vector <double> param_value;                   // The model parameter values

	vector <IndValue> ind_value;                   // Stores individual properties

	Simulate(Model &model);
	void run();
	void set_initial_state(MCMC &mcmc) const;

private:
	void add_infected(const double t_inf, const int i, const int comp_init);
	void set_value(string name, double val);
	void output_datatable();
	void output_matrix() const;
	void output_pedigree();

	Model &model;
};

