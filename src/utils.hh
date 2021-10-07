#ifndef SIRE__UTILS_HH
#define SIRE__UTILS_HH

#include <string>
#include <vector>

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

bool different(double val1, double val2);
int find_in(const vector <int> &vec, const int val);
void check_not_repeated(vector <int> &vec, string name);

void emsg(const string& msg);

#endif