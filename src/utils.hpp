#pragma once

#include <string>
#include <vector>
#include "const.hpp"

using namespace std;

vector <string> split(const string &s, char delimiter);
double number(const string &s);
string num_out(const double x);
void add_to_list(vector <int> &vec, const int val);
void set_seed(const int i);
double ran();
double normal_sample(const double mu, const double sd);
double normal_probability(const double x, const double mean, const double sd);
double gamma_sample(const double mu, const double shape);
double gamma_probability(const double x, const double mu, const double shape);
double exp_sample(const double rate);

bool different(double val1, double val2, double tol = SMALL);
int find_in(const vector <int> &vec, const int val);
void check_not_repeated(vector <int> &vec, string name);
double min(double a, double b);
double max(double a, double b);
void check_same(const vector <double> &vec1, const vector <double> &vec2);
bool MH_proposal(double al, unsigned int num);
void check_one(double al, unsigned int num);

void emsg(const string &msg);

