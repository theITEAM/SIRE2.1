/// This contains generic useful functions

#include <iostream>
#include <string>
#include <vector>
//#include <bits/stdc++.h>
#include <random>
#include <sstream>

using namespace std;

#include "utils.hpp"
#include "const.hpp"

//default_random_engine generator;
mt19937_64 generator;


/// Sets the random seed
void set_seed(const int i) {
	srand(i);
	generator.seed(i);
}


/// Split up a string using a specified delimiter
vector<string> split(string s, char delimiter) 
{
	if(s.length()>1 && s.substr(s.length()-1,1) == "\r") s = s.substr(0,s.length()-1);
	
	vector<string> splits;
	stringstream ss(s);
	string item;
	while (getline(ss, item, delimiter)) {
		if (item != "")
			splits.push_back(item);
	}

	return splits;
}


/// Checks a string has a certain set of characters
bool allow_string(const string st, const string ok_char) {
	for (auto i = 0; i < st.length(); i++) {
		auto j = 0;
		while (j == ok_char.length() && st.substr(i, 1) != ok_char.substr(j, 1))
			j++;
		if (j == ok_char.length())
			return false;
	}
	return true;
}


/// Converts a string to a number
double number(const string &s) {
	string st = s;
	if (st == "infinity")
		return LARGE;
	if (st == "-infinity")
		return -LARGE;

	if (allow_string(st, "-0123456789.") == false)
		emsg("The expression '" + st + "' is not a number.");
	else {
		char *endptr;
		double val = strtod(&st[0], &endptr);
		ptrdiff_t j = endptr - &st[0];
		if (j != (int)st.length())
			emsg("The expression '" + st + "' is not a number.");
		return val;
	}
	return 0.0;
}


/// Converts a number to a string
string num_out(const double x) {
	if (x == LARGE)
		return "infinity";
	if (x == -LARGE)
		return "-infinity";
	return to_string(x);
}


/// Adds a value to a list (if it doesn't exist already)
void add_to_list(vector <int> &vec, const int val) {
	auto i = 0;
	while (i < vec.size() && vec[i] != val)
		i++;
	if (i == vec.size())
		vec.push_back(val);
}


/// Draws a random number between 0 and 1
double ran() {
	return double(0.999999999 * rand()) / RAND_MAX;
}


/// Draws a normally distributed number with mean mu and standard deviation sd
double normal_sample(const double mu, const double sd) {
	normal_distribution<double> distribution(mu, sd);
	return distribution(generator);
}


/// The log of the probability from the normal distribution
double normal_probability(const double x, const double mean, const double sd) {
	if (sd <= 0)
		emsg("Normal probability error");

	return -0.5 * log(2 * M_PI * sd * sd) - (x - mean) * (x - mean) / (2 * sd * sd);
}


/// Draws a sample from the gamma distribution with mean and shape parameter
double gamma_sample(const double mu, const double shape) {
	gamma_distribution<double> distribution(shape, mu / shape);
	return distribution(generator);
}


/// The log of the probability from the gamma distribution
double gamma_probability(const double x, const double mu, const double shape) {
	auto b = shape / mu;
	return (shape - 1) * log(x) - b * x + shape * log(b) - lgamma(shape);
}


/// Samples from the exponential distribution using a rate
double exp_sample(const double rate) {
	exponential_distribution<double> distribution(rate);
	return distribution(generator);
}


/// Returns true if two values are significantly different
bool different(double val1, double val2, double tol) {
	auto av = (val1+val2)/2;
	if(av > 1 || av < -1){ val1 /= av; val2 /= av;}
	if (val1 < val2 - tol || val1 > val2 + tol)
		return true;
	return false;
}


/// Finds the index of a value in a vector
int find_in(const vector <int> &vec, const int val) {
	for (auto i = 0; i < vec.size(); i++) {
		if (vec[i] == val)
			return i;
	}
	return UNSET;
}


/// Checks if a value is repeated multiple times in the model
void check_not_repeated(vector <int> &vec, string name) {
	for (auto i = 0; i < vec.size(); i++) {
		for (auto j = i + 1; j < vec.size(); j++) {
			if (vec[i] == vec[j])
				emsg("Error '" + name + "' is repeated multiple times in the model");
		}
	}
}


/// Returns the mimimum of two numbers
double min(double a, double b)
{
	if(a < b) return a; 
	return b;
}


/// Returns the mimimum of two numbers
double max(double a, double b)
{
	if(a > b) return a;
	return b;
}


/// Checks if two vectors are the same
void check_same(const vector <double> &vec1, const vector <double> &vec2)
{
	auto thresh = 0.001; 
	if(vec1.size() != vec2.size()) emsg("Vectors not the same");
	
	for(auto i = 0; i < vec1.size(); i++){
		if(vec2[i] != 0){
			auto frac = vec1[i]/vec2[i];
			//if(vec1[i] < vec2[i] - thresh || vec1[i] > vec2[i] + thresh){
			if(frac < 1-thresh && frac  > 1+thresh){
				for(auto j = 0; j < vec1.size(); j++){
					cout << i << " " << vec1[j] << " " << vec2[j] << endl;
				}
				emsg("Vectors are different");
			}
		}			
	}
}

	
/// Performs a Metropolis-Hastings proposal
bool MH_proposal(double al, unsigned int num)
{
	if(std::isnan(al)) emsg("Metropolis-Hastings problem nan  "+to_string(num));
	//if(std::isinf(al)) emsg("Metropolis-Hastings problem inf "+to_string(num));
	//if(al == 0) emsg("Metropolis-Hastings problem zero "+to_string(num));

	if(ran() < al) return true;
	return false;
}


/// Checks that a proposal is exactly one (e.g. individual effects on non-group individuals_
void check_one(double al, unsigned int num)
{
	if(al < 0.999 || al > 1.001) emsg("Proposal not one "+num);
}


/// Displays an error message
void emsg(const string &msg) {
	cout << msg << endl;
	exit(EXIT_FAILURE);
}
