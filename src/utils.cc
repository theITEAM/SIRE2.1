/// This contains generic useful functions

#include <iostream>
#include <string>
#include <vector>
#include <bits/stdc++.h>

using namespace std;

#include "utils.hh"
#include "const.hh"

default_random_engine generator;


/// Sets the random seed
void set_seed(const int i)
{
  srand(i);
  generator.seed(i);
}

  
/// Split up a string using a specified delimiter
vector<string> split(const string &s, char delimiter)                                             
{                              
  vector<string> splits;                       
  stringstream ss(s);
  string item;
  while(getline(ss,item,delimiter)){
    if(item != "") splits.push_back(item);
  }
  
  return splits;                                           
}


/// Checks a string has a certain set of characters
bool allow_string(const string st, const string ok_char)
{
  for(auto i = 0u; i < st.length(); i++){
    auto j = 0u; while(j == ok_char.length() && st.substr(i,1) != ok_char.substr(j,1)) j++;
    if(j == ok_char.length()) return false;
  }
  return true;
}


/// Converts a string to a number
double number(const string &s)
{
  string st = s;
  if(st == "infinity") return LARGE;
  if(st == "-infinity") return -LARGE;
  
  if(allow_string(st,"-0123456789.") == false) emsg("The expression '"+st+"' is not a number.");
  else{
    char* endptr;
    double val = strtod(&st[0],&endptr);
    ptrdiff_t j = endptr-&st[0];
    if(j != (int)st.length()) emsg("The expression '"+st+"' is not a number.");
    return val;
  }
}


/// Converts a number to a string
string num_out(const double x)
{
  if(x == LARGE) return "infinity";
  if(x == -LARGE) return "-infinity";
  return to_string(x);
}


/// Adds a value to a list (if it doesn't exist already)
void add_to_list(vector <int> &vec, const int val)
{
  auto i = 0u; while(i < vec.size() && vec[i] != val) i++;
  if(i == vec.size()) vec.push_back(val);
}


/// Draws a random number between 0 and 1
double ran()
{
  return double(0.999999999*rand())/RAND_MAX;
}


/// Draws a normally distributed number with mean mu and standard deviation sd
double normal_sample(const double mu, const double sd)
{
  normal_distribution<double> distribution(mu,sd);
  return distribution(generator);
}


/// The log of the probability from the normal distribution
double normal_probability(const double x, const double mean, const double sd)
{
  if(sd <= 0) emsg("Normal probability error");
  return -0.5*log(2*M_PI*sd*sd) - (x-mean)*(x-mean)/(2*sd*sd);
}


/// Draws a sample from the gamma distribution with mean and shape parameter
double gamma_sample(const double mu, const double shape)
{
  gamma_distribution<double> distribution(shape,mu/shape);
  return distribution(generator);
}


/// The log of the probability from the gamma distribution
double gamma_probability(const double x, const double mu, const double shape)
{
  auto b = shape/mu;
  return (shape-1)*log(x) - b*x + shape*log(b) - lgamma(shape);
}


/// Samples from the exponential distribution using a rate
double exp_sample(const double rate)
{
  exponential_distribution<double> distribution(rate);
  return distribution(generator);
}


/// Returns true if two values are significantly different
bool different(double val1, double val2)
{
  if(val1 < val2-SMALL || val1 > val2+SMALL) return true;
  return false;
}


/// Finds the index of a value in a vector
int find_in(const vector <int> &vec, const int val)
{
  for(auto i = 0u; i < vec.size(); i++){
    if(vec[i] == val) return i;
  }
  return UNSET;
}


/// Checks if a value is repeated multiple times in the model
void check_not_repeated(vector <int> &vec, string name)
{
	for(auto i = 0u; i < vec.size(); i++){
		for(auto j = i+1; j < vec.size(); j++){
			if(vec[i] == vec[j]) emsg("Error '"+name+"' is repeated multiple times in the model");
		}
	}
}


/// Displays an error message
void emsg(const string& msg)
{
  cout << msg << endl;
  exit (EXIT_FAILURE);
}
