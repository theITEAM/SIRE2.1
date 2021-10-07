/// This stores functions related to individual effects

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "model.hh"
#include "timers.hh"
#include "utils.hh"


/// Calculates the log likelihood for the individual effects and returns a vector
vector <double> Model::calculate_L_ind_effect(const vector <IndValue> &ind_value, const vector <double> &param_value) const
{
  vector <double> L_ind_effect(ncovariance);
  for(auto c = 0u; c < ncovariance; c++){
    L_ind_effect[c] = likelihood_ind_effects(covariance[c],ind_value,param_value);
  }
  return L_ind_effect;
}


/// Calculates the log likelihood for the individual effects for a given covariance matrix
double Model::likelihood_ind_effects(const Covariance &cov, const vector <IndValue> &ind_value, const vector <double> &param_value) const
{
  auto L = 0.0;
  const auto &mat = matrix[cov.matrix];
  
  auto cov_matrix = set_covariance_matrix(cov,param_value);
  auto inverse_cov_matrix = invert_matrix(cov_matrix);
  auto E = cov.E;
  vector <double> vecj(E), veci(E);
  
  for(auto j = 0u; j < N; j++){
    for(auto e = 0u; e < E; e++) vecj[e] = ind_value[j].ind_effect[cov.ind_effect_ref[e]];
    
    // Does the diagonal elements
    auto val = mat.Ainvdiag[j];
    
    auto sum = 0.0;
    for(auto ej = 0u; ej < E; ej++){
      for(auto ei = 0u; ei < E; ei++){
        sum += vecj[ej]*inverse_cov_matrix[ej][ei]*vecj[ei];
      }
    }
    L -= 0.5*sum*val;
  
    // Does the off-diagnoal elements
    for(const auto &ele : mat.Ainvlist2[j]){
      auto i = ele.i;
      auto val = ele.val;
      
      for(auto e = 0u; e < E; e++) veci[e] = ind_value[i].ind_effect[cov.ind_effect_ref[e]];
    
      auto sum = 0.0;
      for(auto ej = 0u; ej < E; ej++){
        for(auto ei = 0u; ei < E; ei++){
          sum += vecj[ej]*inverse_cov_matrix[ej][ei]*veci[ei];
        }
      }
      L -= sum*val;
    }
  }
  auto det = determinant(cov_matrix);
  L -= 0.5*N*log(det);    
  
  return L;
}


/// This stores the inverse covariance matrices
vector <InvCovMat> Model::calculate_inv_cov_matrix(const vector <double> &param_value) const
{
  vector <InvCovMat> inv_cov_mat;
  for(auto cv = 0u; cv < ncovariance; cv++){
    auto cov_matrix = set_covariance_matrix(covariance[cv],param_value);
    InvCovMat inv; inv.M = invert_matrix(cov_matrix);
    inv_cov_mat.push_back(inv);
  }
  
  return inv_cov_mat;
}


/// This makes proposals to parameters in the covariance matrices
void Model::propose_covariance_matrices(const Covariance &cov, const vector <IndValue> &ind_value, vector <double> &param_value, double &L_ind_effect, double &prior, vector <Jump> &param_jump, const bool burnin) const
{
  timer[TIME_COVAR_INIT].start();
  auto precalc = set_precalc_cov(cov,ind_value);
  timer[TIME_COVAR_INIT].stop();
  
  timer[TIME_COVAR].start();
  auto L_init = Li_cov_fast(cov,precalc,param_value);
  auto L = L_init;
  for(auto loop = 0; loop < covar_prop; loop++){
    for(auto th : cov.param_list){
      auto &jump = param_jump[th];
      auto param_store = param_value[th];

      // Makes a change to one of the covariance matrix parameters
      param_value[th] += normal_sample(0,jump.size);
      
      // Calculates the Metropolis-Hastings acceptance probability
      auto al = 0.0, L_propose = L; 
      auto prior_change = calculate_prior_change(th,param_store,param_value);
      if(prior_change != ZERO_PRIOR){
        L_propose = Li_cov_fast(cov,precalc,param_value);
        al = exp(L_propose - L + prior_change);
      }
    
      jump.ntr++;
      if(ran() < al){
        jump.nac++;
        L = L_propose;
        prior += prior_change;
        if(burnin == true) jump.size *= prop_up; 
      }
      else{
        param_value[th] = param_store;
        if(burnin == true) jump.size *= prop_down; 
      }
    }
  }
  
  L_ind_effect += L-L_init;
  
  timer[TIME_COVAR].stop();
}


/// This precalculates a matrix which can be used to calculate the likelihood faster
vector < vector <double> > Model::set_precalc_cov(const Covariance &cov, const vector <IndValue> &ind_value) const
{
  vector < vector <double> > precalc;
  auto E = cov.E;
  vector <double> vecj(E), veci(E);
  precalc.resize(E);
  for(auto ej = 0u; ej < E; ej++){
    precalc[ej].resize(E);
    for(auto ei = 0u; ei < E; ei++) precalc[ej][ei] = 0;
  }
  
  const auto &mat = matrix[cov.matrix];
  for(auto j = 0u; j < N; j++){
    for(auto e = 0u; e < E; e++) vecj[e] = ind_value[j].ind_effect[cov.ind_effect_ref[e]];
    
    // Does the diagonal elements
    auto val = mat.Ainvdiag[j];
    
    auto sum = 0.0;
    for(auto ej = 0u; ej < E; ej++){
      for(auto ei = 0u; ei < E; ei++){
        precalc[ej][ei] -= 0.5*val*vecj[ej]*vecj[ei];
      }
    }
    
    // Does the off-diagnoal elements
    for(const auto &ele : mat.Ainvlist2[j]){
      auto i = ele.i;
      auto val = ele.val;
      
      for(auto e = 0u; e < E; e++) veci[e] = ind_value[i].ind_effect[cov.ind_effect_ref[e]];
    
      auto sum = 0.0;
      for(auto ej = 0u; ej < E; ej++){
        for(auto ei = 0u; ei < E; ei++){
          precalc[ej][ei] -= val*vecj[ej]*veci[ei];
        }
      }
    }
  }
  
  return precalc;
}


/// This performs a fast likelihood calculation (based onpre-calculated quantities)
double Model::Li_cov_fast(const Covariance &cov, const vector < vector <double> > &precalc, const vector <double> &param_value) const
{
  auto cov_matrix = set_covariance_matrix(cov,param_value);
  auto inverse_cov_matrix = invert_matrix(cov_matrix);
  auto E = cov.E;

  auto L = 0.0;
  for(auto ej = 0u; ej < E; ej++){
    for(auto ei = 0u; ei < E; ei++){
      L += precalc[ej][ei]*inverse_cov_matrix[ej][ei];
    }
  }   
  auto det = determinant(cov_matrix);
  L -= 0.5*N*log(det);    
  
  return L;
}

