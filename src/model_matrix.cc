/// This provides code in the model class which is related to matrix manipulation

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include "tinyxml2.h"

using namespace std;

using namespace tinyxml2;

#include "model.hh"
#include "utils.hh"
#include "const.hh"


/// Adds the identity matrix to the model
void Model::add_identity_matrix()
{
  Matrix mat;
  mat.name = "I";
  mat.Ainvlist.resize(N);
  mat.Ainvlist2.resize(N);
  mat.Ainvdiag.resize(N);
  for(auto i = 0u; i < N; i++) mat.Ainvdiag[i] = 1;
  
  mat.A.resize(N);
  for(auto j = 0; j < N; j++){
    mat.A[j].resize(N);
    for(auto i = 0u; i < N; i++){
      if(i == j) mat.A[j][i] = 1;
      else mat.A[j][i] = 0;
    }
  }
  
  bool illegal;
  mat.cholesky_matrix = calculate_cholesky_matrix(mat.A,illegal);
  if(illegal == true) emsg("Cannot calculate the cholesky matrix for '"+mat.name+"'");
    
  matrix.push_back(mat);
}


/// Adds a relationship matrix to the model
void Model::add_matrix(XMLNode *child)
{
  auto text = child->FirstChild()->Value();
  auto tab = create_table(text);
  
  Matrix mat;
  mat.name = get(child,"name");
  
	vector < vector <double> > Ainv; 
	
  // Loads the matrix
    
	if(tab_name(child) == "pedigree"){
		Ainv.resize(N);
		for(auto j = 0; j < N; j++){
			Ainv[j].resize(N);
			for(auto i = 0u; i < N; i++){
				if(i == j) Ainv[j][i] = 1;
				else Ainv[j][i] = 0;
			}
		}
		
		if(tab.ncol != 3) emsg("For 'pedigree' there should be three columns: first the individuals and then the two parents");
		for(auto r = 0u; r < tab.nrow; r++){
			auto i = find_ind(tab.ele[r][0]);
			auto par1 = find_ind(tab.ele[r][1]);
			auto par2 = find_ind(tab.ele[r][2]);
			
			if(i == UNSET) emsg("In 'pedigree' the indiviudal id must be set");
			
			if(par1 == UNSET && par2 == UNSET){  // Both parents unknown
			}
			else{
				if(par1 == UNSET || par2 == UNSET){ // One parent known
					auto p = par1; if(par1 == UNSET) p = par2;
					Ainv[p][p] += 1.0/3;
					Ainv[i][i] += 1.0/3;
					Ainv[p][i] -= 2.0/3;
					Ainv[i][p] -= 2.0/3;
				}
				else{                              // Both parents known
					Ainv[par1][par1] += 1.0/2;
					Ainv[par2][par2] += 1.0/2;		
					Ainv[i][i] += 1.0;
					
					Ainv[par1][par2] += 1.0/2;
					Ainv[par2][par1] += 1.0/2;
					Ainv[par1][i] -= 1.0;
					Ainv[i][par1] -= 1.0;
					Ainv[par2][i] -= 1.0;
					Ainv[i][par2] -= 1.0;
				}
			}	
		}
		
		mat.A = invert_matrix(Ainv);
	}
	else{
		mat.A.resize(N);
		for(auto j = 0; j < N; j++) mat.A[j].resize(N);
 
		if(tab_name(child) == "A"){
			for(auto j = 0; j < N; j++){
				for(auto i = 0u; i < N; i++){
					mat.A[j][i] = number(tab.ele[j][i]);
				}
			}
		}
		
		if(tab_name(child) == "A_nonzero"){
			for(auto j = 0; j < N; j++){
				for(auto i = 0u; i < N; i++) mat.A[j][i] = 0;
			}
			
			if(tab.ncol != 3) emsg("For 'A_nonzero' there should be three columns");
			
			for(auto r = 0u; r < tab.nrow; r++){
				auto j = int(number(tab.ele[r][0]));
				auto i = int(number(tab.ele[r][1]));
				auto val = number(tab.ele[r][2]);
				if(j < 0 || j >= N || i < 0 || i >= N) emsg("Loading matrix problem");
				mat.A[j][i] = val;
			}
		}
		
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i < N; i++){
				if(mat.A[j][i] != mat.A[i][j]){
					emsg("The matrix '"+mat.name+"' specified by the '"+tab_name(child)+"' tag must be symmetric");
				}
			}
		}
		
		Ainv = invert_matrix(mat.A);
	}
	
  // Generates non-zero elements of the inverse matrix 
 
  mat.Ainvlist.resize(N);
  mat.Ainvlist2.resize(N);
  mat.Ainvdiag.resize(N);
  for(auto j = 0; j < N; j++){
    for(auto i = 0u; i < N; i++){
      auto val = Ainv[j][i];
      if(j != i){
        if(val != 0){
          Element ele; ele.i = i; ele.val = val;
          mat.Ainvlist[j].push_back(ele);
          if(i < j) mat.Ainvlist2[j].push_back(ele);
        }
      } 
      else mat.Ainvdiag[j] = val;
    }
  }
  
  bool illegal;
  mat.cholesky_matrix = calculate_cholesky_matrix(mat.A,illegal);
  if(illegal == true) emsg("Cannot calculate the cholesky matrix for '"+mat.name+"'");
  
  matrix.push_back(mat);
}


/// Adds a covariance matrix to the model
void Model::add_covariance(XMLNode *child)
{
  Covariance covar;
  auto vec = split(get(child,"individual_effect"),',');
  auto E = vec.size();
  covar.E = E;

  for(auto e = 0; e < E; e++){
    auto i = 0u; while(i < ind_effect.size() && vec[e] != ind_effect[i].name) i++;
    if(i == ind_effect.size()) emsg("Cannot find individual effect '"+vec[e]+"'.");
    
    if(ind_effect[i].covar_ref != UNSET){
      emsg("Individual effect '"+ind_effect[i].name+"' cannot appear in more than one covariance matrix");
    }
    
    ind_effect[i].covar_ref = covariance.size();
    ind_effect[i].covar_num = e;
    
    covar.ind_effect_ref.push_back(i);
  }
  
  auto Aname = get(child,"relationship_matrix");
  auto i = 0u; while(i < matrix.size() && Aname != matrix[i].name) i++;
  if(i == matrix.size()) emsg("Cannot find matrix '"+Aname+"'");
  covar.matrix = i;
  
  for(auto tag = child->FirstChild(); tag; tag = tag->NextSibling()){
    if(tab_name(tag) == "variance"){
      auto text = tag->FirstChild()->Value();
      auto tab = create_table(text);
      covar.var_param.resize(E);
      for(auto row = 0u; row < E; row++){
        auto th = find_param(tab.ele[row][0],COVAR_MATRIX);
        covar.var_param[row] = th;
        add_to_list(covar.param_list,th);
      }
    }
    
    if(tab_name(tag) == "correlation"){
      auto text = tag->FirstChild()->Value();
      auto tab = create_table(text);
      covar.cor_param.resize(E);
      for(auto row = 0u; row < E; row++){
        covar.cor_param[row].resize(E);
        for(auto col = 0u; col < E; col++){
          if(row == col){
            if(tab.ele[row][col] != "1"){
              emsg("Correlation matrices should have 1 on the diagonal");
            }
            covar.cor_param[row][col] = UNSET;
          }
          else{
            auto th = find_param(tab.ele[row][col],COVAR_MATRIX);
            covar.cor_param[row][col] = th;
            add_to_list(covar.param_list,th);
          }
        }
      }
    }
  }
  
  if(covar.var_param.size() == 0){
    emsg("The covariance matrix must have parameters specifying the variance");
  }
  
  if(covar.cor_param.size() == 0){
    if(E == 1){
      covar.cor_param.resize(1);
      covar.cor_param[0].resize(1);
      covar.cor_param[0][0] = UNSET;
    }
    else{
      emsg("A correlation matrix must be set.");
    }
  }
  
  covariance.push_back(covar);
}


/*
/// Checks that the covariance matrices are all valid 
bool Model::covariance_valid(const vector <double> param_value) const
{
  for(const auto &cov : covariance){
    auto cov_mat = set_covariance_matrix(cov,param_value);
  
    bool illegal;
    calculate_cholesky_matrix(cov_mat,illegal);
    if(illegal == true) return false;
  }
  return true;
}
*/


/// Inverts a matrix
vector <vector <double> > Model::invert_matrix(const vector <vector <double> > &mat) const
{
  unsigned int nvar = mat.size();

  vector <vector <double> > inv_M, A2;
  
  inv_M.resize(nvar); A2.resize(nvar);
  for(auto i = 0u; i < nvar; i++){
    inv_M[i].resize(nvar); A2[i].resize(nvar);
    for(auto j = 0u; j < nvar; j++){
      A2[i][j] = mat[i][j];
      if(i == j) inv_M[i][j] = 1; else inv_M[i][j] = 0;
    }
  }

  for(auto ii = 0u; ii < nvar; ii++){
    double r = A2[ii][ii];
    for(auto i = 0u; i < nvar; i++){
      A2[ii][i] /= r; inv_M[ii][i] /= r; 
    }

    for(auto jj = ii+1; jj < nvar; jj++){
      double r = A2[jj][ii];
      if(r != 0){
        for(auto i = 0u; i < nvar; i++){ 
          A2[jj][i] -= r*A2[ii][i];
          inv_M[jj][i] -= r*inv_M[ii][i];
          
          if(inv_M[jj][i] > -TINY && inv_M[jj][i] < TINY) inv_M[jj][i] = 0;
          if(A2[jj][i] > -TINY && A2[jj][i] < TINY) A2[jj][i] = 0;
        }
      }
    }
  }

  for(int ii = nvar-1; ii > 0; ii--){
    for(int jj = ii-1; jj >= 0; jj--){
      double r = A2[jj][ii];
      if(r != 0){
        for(auto i = 0u; i < nvar; i++){ 
          A2[jj][i] -= r*A2[ii][i];
          inv_M[jj][i] -= r*inv_M[ii][i];
          
          if(inv_M[jj][i] > -TINY && inv_M[jj][i] < TINY) inv_M[jj][i] = 0;
          if(A2[jj][i] > -TINY && A2[jj][i] < TINY) A2[jj][i] = 0;
        }
      }
    }
  }
  
  if(false){
    ofstream ch("ch.txt");
    for(auto j = 0u; j < nvar; j++){
      for(auto i = 0u; i < nvar; i++){
        ch << inv_M[j][i]<< " ";
      }
      ch << endl;
    }
  }
  
  if(false){ // checks inverse
    for(auto j = 0u; j < nvar; j++){
      for(auto i = 0u; i < nvar; i++){
        double sum = 0; for(auto ii = 0u; ii < nvar; ii++) sum += mat[j][ii]*inv_M[ii][i];
        if(i != j){ if(sum < -TINY || sum > TINY) emsg("Inverse Problem");}
        else{ if(sum < 1-TINY || sum > 1+TINY) emsg("Inverse Problem");}    
      }
    }
  }
  
  return inv_M;
}


/// Calculates a lower diagonal matrix used in Cholesky decomposition
vector < vector <double> > Model::calculate_cholesky_matrix(const vector < vector <double> > &M, bool &illegal) const 
{ 
  auto nvar = M.size();

  illegal = false;
   
  vector <vector <double> > L;
  L.resize(nvar);
  for(auto i = 0u; i < nvar; i++) L[i].resize(nvar); 
  
  for(auto j = 0u; j < nvar; j++){
    for(auto i = 0u; i <= j; i++){
      auto sum = 0.0;
      for(auto k = 0u; k < i; k++) sum += L[j][k]*L[i][k];

      if(j == i){
        auto val = M[j][j] - sum;
        if(val < 0){ illegal = true; return L;}
        L[j][i] = sqrt(val);
      }
      else L[j][i] = (1.0 / L[i][i] * (M[j][i] - sum));
    } 
  }

  if(false){   // Check that M = L*LT
    for(auto j = 0u; j < nvar; j++){
      for(auto i = 0u; i < nvar; i++){
        auto sum = 0.0; for(auto ii = 0u; ii < nvar; ii++) sum += L[j][ii]*L[i][ii];
    
        if(M[j][i] < sum-TINY || M[j][i] > sum+TINY) emsg("Inverse Problem");
      }
    }
  }
  
  return L;
}


/// Returns the tensor product of two matrices
vector < vector <double> > Model::tensor_product(const vector < vector <double> > &mat1, const vector < vector <double> > &mat2) const
{
  auto N1 = mat1.size();
  auto N2 = mat2.size();
  auto N = N1*N2;
    
  vector < vector <double> > mat;
  mat.resize(N);
  for(auto j = 0u; j < N; j++){
    auto jj = int(j/N2);
    auto jjj = j%N2;
    
    mat[j].resize(N);
    for(auto i = 0u; i < N; i++){
      auto ii = int(i/N2);
      auto iii = i%N2;
      
      mat[j][i] = mat1[jj][ii]*mat2[jjj][iii];
    }
  } 
  
  return mat;
}


/// Sets the covariance matrix using the parameter value vector
vector < vector <double> > Model::set_covariance_matrix(const Covariance &cov, const vector <double> &param_value) const
{
  auto E = cov.E;
  vector < vector <double> > cov_mat;
  cov_mat.resize(E);
  for(auto j = 0u; j < E; j++){
    cov_mat[j].resize(E);
    for(auto i = 0u; i < E; i++){
      if(i == j){
        cov_mat[i][i] = param_value[cov.var_param[i]];
      }
      else{
        auto var_i = param_value[cov.var_param[i]];
        auto var_j = param_value[cov.var_param[j]];
        auto cor = param_value[cov.cor_param[j][i]];
        cov_mat[j][i] = cor*sqrt(var_i*var_j);
      }
    }
  }
  
  return cov_mat;
}


/// Calculates the determinant of a matrix
double Model::determinant(const vector < vector <double> > &M) const
{
  auto n = M.size();
  if(n == 1) return M[0][0];
  if(n == 2) return M[0][0]*M[1][1]-M[0][1]*M[1][0];

  vector < vector <double> > temp;
  temp.resize(n-1); for(auto j = 0u; j < n-1; j++) temp[j].resize(n-1); 
  
  auto det = 0.0;
  for(auto p = 0u; p < n; p++) {
    auto h = 0u, k = 0u;
    for(auto i = 1u; i < n; i++){
      for(auto j = 0u; j < n; j++){
        if(j != p){
          temp[h][k] = M[i][j];
          k++;
          if(k == n-1) {
            h++;
            k = 0;
          }
        }
      }
    }
    det = det + M[0][p]*pow(-1,p)*determinant(temp);
  }
  
  return det;
}


/// Prints a matrix
void Model::print_matrix(string name, const vector < vector <double> > &mat) const 
{
  cout << "Matrix " << name << endl;
  for(auto j = 0u; j < mat.size(); j++){
    for(auto i = 0u; i < mat.size(); i++){
      cout << mat[j][i] <<" ";
    }
    cout << endl;
  }
}


/// Samples individual effects for a given individual i
void Model::ind_effect_sample(vector <IndValue> &ind_value, const vector <double> &param_value) const
{
  // Initialises the vector of individual effects for each individuals 
  for(auto &ind_val : ind_value){
    ind_val.ind_effect.resize(nind_effect);
  }
  
  // Sets the individual effects based on the specified covariance matrices
  for(const auto &cov : covariance){
    auto cov_mat = set_covariance_matrix(cov,param_value);
  
    bool illegal;
    auto cov_cholesky_matrix = calculate_cholesky_matrix(cov_mat,illegal);
    if(illegal == true) emsg("Problem with getting cholesky_matrix");
    
    auto M = tensor_product(matrix[cov.matrix].cholesky_matrix,cov_cholesky_matrix);
    
    auto E = cov.E;
    auto Z = N*E;
    if(Z != M.size()) emsg("Matrix not the right size");
    
    /// This generate a random sample from the MVN distribution
    vector <double> z(Z);
    for(auto i = 0u; i < Z; i++) z[i] = normal_sample(0,1);
    
    for(auto j = 0u; j < Z; j++){
      auto sum = 0.0;
      for(auto i = 0u; i < Z; i++) sum += M[j][i]*z[i];
      
      auto i = int(j/E);
      auto e = j%E;
      
      ind_value[i].ind_effect[cov.ind_effect_ref[e]] = sum;
    }   
  }
}
