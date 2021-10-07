/// This gives model functions related to the prior

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "model.hh"
#include "utils.hh"


/// Samples a set of values from the prior
vector <double> Model::prior_sample() const
{
  vector <double> param_value(nparam);
  
  for(auto th = 0u; th < nparam; th++){
    const auto &par = param[th];
    switch(par.prior_type){
      case FLAT_PRIOR:
        param_value[th] = par.prior_val1 + ran()*(par.prior_val2 - par.prior_val1);
        break;
      
      case NORMAL_FROM_SD_PRIOR:
        param_value[th] = normal_sample(0,param_value[par.prior_sd_param]);
        break;

      default: 
        emsg("Prior not specified");
        break;
    }
  }
  
  return param_value;
}


/// Calculates the log of the prior probability
double Model::calculate_prior(const vector <double> &param_value) const
{
  auto prior = 0.0;
  for(auto th = 0u; th < nparam; th++){
    prior += calculate_prior(th,param_value[th],param_value);   
  }
  return prior;
}


/// Calculates the prior for a particular parameter th with a value 
double Model::calculate_prior(const int th, const double value, const vector <double> &param_value) const
{
  const auto &par = param[th];
  switch(par.prior_type){
    case FLAT_PRIOR:
      {
        auto val1 = par.prior_val1;
        auto val2 = par.prior_val2;
  
        if(value < val1 || value > val2) return ZERO_PRIOR;
        else return log(1.0/(val2-val1));
      }
      
    case NORMAL_FROM_SD_PRIOR:
      {
        auto sd = param_value[par.prior_sd_param];
        if(sd <= 0) return ZERO_PRIOR;
        else return normal_probability(value,0,param_value[par.prior_sd_param]);
      }
      
    default: 
      emsg("Prior not specified");
      break;
  }
}


/// Calculates the change in prior when a particular parameter th changes
double Model::calculate_prior_change(const int th, const double val_before, const vector <double> &param_value) const
{
  auto prior_after = calculate_prior(th,param_value[th],param_value);
  if(prior_after == ZERO_PRIOR) return ZERO_PRIOR;
  auto prior_before = calculate_prior(th,val_before,param_value);
  
  return prior_after - prior_before;
}
