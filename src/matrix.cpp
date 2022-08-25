/// Routines for matrix manipulation

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "const.hpp"
#include "matrix.hpp"
#include "utils.hpp"

/// Inverts a matrix
vector < vector <double> > invert_matrix(vector <vector <double> > M)   
{
	int nvar = M.size();
	
	vector <vector <double> > inv_M;
	
	inv_M.resize(nvar);
  for(auto i = 0; i < nvar; i++) inv_M[i].resize(nvar,0);
	for(auto i = 0; i < nvar; i++) inv_M[i][i] = 1;
   
  for(auto ii = 0; ii < nvar; ii++){
		auto &Mii = M[ii], &inv_Mii = inv_M[ii];
		
    auto r = Mii[ii];
    for(auto i = ii; i < nvar; i++) Mii[i] /= r; 
	  for(auto i = 0; i <= ii; i++) inv_Mii[i] /= r; 
	
    for(auto jj = ii+1; jj < nvar; jj++){
      auto r = M[jj][ii];
			if(r != 0){
				auto &Mjj = M[jj], &inv_Mjj = inv_M[jj];
			
				for(auto i = ii; i < nvar; i++) Mjj[i] -= r*Mii[i];
				for(auto i = 0; i <= ii; i++) inv_Mjj[i] -= r*inv_Mii[i];
			}
    }
  }

  for(int ii = nvar-1; ii > 0; ii--){
		//auto &Mii = M[ii];
		auto &inv_Mii = inv_M[ii];
    for(int jj = ii-1; jj >= 0; jj--){
      auto r = M[jj][ii];
			if(r != 0){
				//auto &Mjj = M[jj];
				auto &inv_Mjj = inv_M[jj];
		
				//for(auto i = 0; i < nvar; i++) Mjj[i] -= r*Mii[i];
				for(auto i = 0; i < nvar; i++) inv_Mjj[i] -= r*inv_Mii[i];
			}
    }
  }

	if(false){ // checks inverse
		cout << "Check matrix" << endl;
		for(auto j = 0; j < nvar; j++){
			for(auto i = 0; i < nvar; i++){
				double sum = 0; for(auto ii = 0; ii < nvar; ii++) sum += M[j][ii]*inv_M[ii][i];
				
				cout << sum << " ";
				if(i != j){ if(sum < -TINY || sum > TINY) emsg("Matrix1");}
				else{ if(sum < 1-TINY || sum > 1+TINY) emsg("Matrix2");}		
			}
			cout << endl;
		}
	}
	
	return inv_M;
}


/// Inverts a sparse matrix
vector < vector <double> > invert_matrix_sparse(vector <vector <double> > M)   
{
	int nvar = M.size();
	//auto M_st = M;
	vector <vector <double> > inv_M;
	
	inv_M.resize(nvar);
  for(auto i = 0; i < nvar; i++) inv_M[i].resize(nvar,0);
	for(auto i = 0; i < nvar; i++) inv_M[i][i] = 1;
  
	vector < vector <unsigned int> > M_ele, inv_M_ele; // Store non-diagonal zero elements
	vector < vector <unsigned int> > M_map, inv_M_map; // References position on list
	
	M_ele.resize(nvar); inv_M_ele.resize(nvar);
	M_map.resize(nvar); inv_M_map.resize(nvar);
	for(auto j = 0; j < nvar; j++){cout << j << " nvar\n";
		M_map[j].resize(nvar,UNSET); inv_M_map[j].resize(nvar,UNSET);
			
		add_mat_element(j,inv_M[j],inv_M_ele[j],inv_M_map[j]);
		for(auto i = 0; i < nvar; i++){
			if(M[j][i] != 0) add_mat_element(i,M[j],M_ele[j],M_map[j]);
		}
	}
	 
	//check_sparse("M",M,M_ele,M_map); check_sparse("M_inv",inv_M,inv_M_ele,inv_M_map);
	
  for(auto ii = 0; ii < nvar; ii++){cout <<ii << " nvar down\n";
		auto &Mii = M[ii], &inv_Mii = inv_M[ii];
		
    auto r = Mii[ii];
		for(auto i : M_ele[ii]) Mii[i] /= r; 
	  for(auto i : inv_M_ele[ii]) inv_Mii[i] /= r; 
	
    for(auto jj = ii+1; jj < nvar; jj++){
      auto r = M[jj][ii];
			if(r != 0){
				auto &Mjj = M[jj], &inv_Mjj = inv_M[jj];
				auto &M_elejj = M_ele[jj], &inv_M_elejj = inv_M_ele[jj];
				auto &M_mapjj = M_map[jj], &inv_M_mapjj = inv_M_map[jj];
			
				for(auto i : M_ele[ii]){
					Mjj[i] -= r*Mii[i]; add_mat_element(i,Mjj,M_elejj,M_mapjj);
				}	
				rem_mat_element(ii,Mjj,M_elejj,M_mapjj);
				
				for(auto i : inv_M_ele[ii]){
					inv_Mjj[i] -= r*inv_Mii[i];
					add_mat_element(i,inv_Mjj,inv_M_elejj,inv_M_mapjj);
				}
			}
    }
  }
	
		cout << sparsity(inv_M) << "invM\n";
	//check_sparse("M",M,M_ele,M_map); check_sparse("M_inv",inv_M,inv_M_ele,inv_M_map);
	
  for(int ii = nvar-1; ii > 0; ii--){ cout << ii << " back\n";
		auto &inv_Mii = inv_M[ii];
    for(int jj = ii-1; jj >= 0; jj--){
      auto r = M[jj][ii];
			if(r != 0){
				auto &inv_Mjj = inv_M[jj];
				auto &inv_M_elejj = inv_M_ele[jj];
				auto &inv_M_mapjj = inv_M_map[jj];
			
				//for(auto i = 0; i < nvar; i++) M[jj][i] -= r*M[ii][i];
			
				for(auto i : inv_M_ele[ii]){
					inv_Mjj[i] -= r*inv_Mii[i];
					add_mat_element(i,inv_Mjj,inv_M_elejj,inv_M_mapjj);
				}
			}
    }
  }
	cout << sparsity(inv_M) << "invM\n";
	emsg("L");
	//check_sparse("M_inv",inv_M,inv_M_ele,inv_M_map);
	
	if(false){ // checks inverse
		auto M_st = M;
		cout << "Check matrix" << endl;
		for(auto j = 0; j < nvar; j++){
			for(auto i = 0; i < nvar; i++){
				double sum = 0; for(auto ii = 0; ii < nvar; ii++) sum += M_st[j][ii]*inv_M[ii][i];
				
				if(i != j){ if(sum < -TINY || sum > TINY) emsg("Matrix3");}
				else{ if(sum < 1-TINY || sum > 1+TINY) emsg("Matrix4");}		
			}
		}
		emsg("do");
	}
	
	return inv_M;
}


/// For a sparse matrix adds a matrix element
void add_mat_element(unsigned int i, vector <double> &Mjj, vector <unsigned int> &M_elejj,  vector <unsigned int> &M_mapjj)
{	
	if(Mjj[i] != 0 && M_mapjj[i] == UNSET){
		M_mapjj[i] = M_elejj.size(); M_elejj.push_back(i);
	}
}


/// For a sparse matrix removes a matrix element
void rem_mat_element(unsigned int i, vector <double> &Mjj, vector <unsigned int> &M_elejj,  vector <unsigned int> &M_mapjj)
{
	if(Mjj[i] < -TINY || Mjj[i] > TINY) emsg("prb");
	Mjj[i] = 0;
	auto k = M_mapjj[i]; if(k == UNSET) emsg("problem");
	M_mapjj[i] = UNSET;
	auto kmax = M_elejj.size();
	if(k < kmax-1){
		M_elejj[k] = M_elejj[kmax-1];
		M_mapjj[M_elejj[k]] = k;
	}
	M_elejj.pop_back();
}


/// Checks properties of the sparse matrix
void check_sparse(string name, const vector < vector <double> > &M, const vector < vector <unsigned int> > &M_ele, const vector < vector <unsigned int> > &M_map)
{
	auto nvar = M.size();
	
	for(auto j = 0; j < nvar; j++){
		for(auto k = 0; k < M_ele[j].size(); k++){
			auto i = M_ele[j][k];
			if(M_map[j][i] != k) emsg(name+"prob1");
		}
		
		for(auto i = 0; i < nvar; i++){
			if(M[j][i] != 0 && M_map[j][i] == UNSET) emsg(name+"prob2");
			if(M_map[j][i] != UNSET && M_map[j][i] >= M_ele[j].size()) emsg(name+"prob3");
		}
	}		
}


/// Calculates the square root of the inverse matrix using Denman–Beavers iteration
vector < vector <double> > invert_matrix_square_root(const vector < vector <double> > &M)
{
	auto n = M.size();
	
	auto Y = M;
	
	vector < vector <double> > Z;  // Sets the identity matrix
	Z.resize(n);
	for(auto j = 0; j < n; j++){
		Z[j].resize(n);
		for(auto i = 0; i < n; i++){
			if(i == j) Z[j][i] = 1;
			else Z[j][i] = 0;
		}
	}

	auto limit = TINY;
	auto loop = 0, loopmax = 1000;
	do{
		auto Yinv = invert_matrix(Y);
		auto Zinv = invert_matrix(Z);
		
		vector < vector <double> > Ynew;
		Ynew.resize(n);
		for(auto j = 0; j < n; j++){
			Ynew[j].resize(n);
			for(auto i = 0; i < n; i++){
				Ynew[j][i] = 0.5*(Y[j][i] + Zinv[j][i]);
			}
		}
		
		vector < vector <double> > Znew;
		Znew.resize(n);
		for(auto j = 0; j < n; j++){
			Znew[j].resize(n);
			for(auto i = 0; i < n; i++){
				Znew[j][i] = 0.5*(Z[j][i] + Yinv[j][i]);
			}
		}
		
		auto dif = LARGE;
		for(auto j = 0; j < n; j++){
			for(auto i = 0; i < n; i++){
				auto d = Znew[j][i]-Z[j][i];
				if(d < 0) d = -d;
				if(d < dif) dif = d;
			}
		}
		
		Y = Ynew; Z = Znew;
		if(dif < limit) break;
		loop++; if(loop%100 == 0) limit *= 10;
	}while(loop < loopmax);
	if(loop == loopmax) emsg("Matrix5");
	
	if(false){
		auto Minv = invert_matrix(M);
		print_matrix("Minv",Minv);
		
		auto N = matrix_mult(Z,Z);
		print_matrix("N",N);
		
		for(auto j = 0; j < n; j++){
			for(auto i = 0; i < n; i++){
				auto d = Minv[j][i] - N[j][i];
				if(d < -TINY || d > TINY){
					emsg("Matrix6");
				}
			}
		}		
		emsg("matrix check");
	}
	
	return Z;
}


/// Prints a matrix
void print_matrix(string name, const vector < vector <double> > &M)
{
	cout << name << ":" << endl;	
	for(auto j = 0; j < M.size(); j++){
		cout << "   ";
		for(auto i = 0; i < M[j].size(); i++){
			cout << M[j][i] << " ";
		}
		cout << endl;
	}
}


/// Prints a vector
void print_vector(string name, const vector <double> &v)
{
	cout << name << ": ";	
	for(auto i = 0; i < v.size(); i++){
		cout << v[i] << " ";
	}
	cout << endl;
}


/// Multiplies a matrix and a vector
vector <double> matrix_mult(const vector < vector <double> > &M, const vector <double> &vec)
{
	const auto NY = M.size();
	if(NY == 0) emsg("Matrix has zero size");
	const auto NX = M[0].size();
	vector <double> result(NY);
	
	if(vec.size() != NX) emsg("Matrix7");
	
	for(auto j = 0; j < NY; j++){
		auto sum = 0.0;
		for(auto i = 0; i < NX; i++) sum += M[j][i]*vec[i];
		result[j] = sum;
	}
	
	return result;
}

/// Multiplies a vector and a matrix
vector <double> matrix_mult(const vector <double> &vec, const vector < vector <double> > &M)
{
	const auto NY = M.size();
	if(NY == 0) emsg("Matrix has zero size");
	const auto NX = M[0].size();
	vector <double> result(NY);
	
	if(vec.size() != NY) emsg("Matrix7");
	
	for(auto i = 0; i < NX; i++){
		auto sum = 0.0;
		for(auto j = 0; j < NY; j++) sum += vec[j]*M[j][i];
		result[i] = sum;
	}
	
	return result;
}


/// Multiplies a matrix and a matrix
vector < vector <double> > matrix_mult(const vector < vector <double> > &M1, const vector < vector <double> > &M2)
{	
	const auto NY = M1.size(), NX = M1[0].size();
	const auto NY2 = M2.size(), NX2 = M2[0].size();
	vector < vector <double> > result;
	
	if(NX != NY2) emsg("Matrix8");
	
	result.resize(NY);
	for(auto j = 0; j < NY; j++){
		result[j].resize(NX2);
		for(auto k = 0; k < NX2; k++){
			auto sum = 0.0;
			for(auto i = 0; i < NX; i++) sum += M1[j][i]*M2[i][k];
			result[j][k] = sum;
		}
	}	

	return result;
}

/// Calculates the cholesky decomposition of a matrix
vector < vector <double> > calculate_cholesky(const vector < vector <double> > &M)
{
	vector < vector <double> > CM;
	auto nvar = M.size();
	
	CM.resize(nvar);
	for(auto v1 = 0; v1 < nvar; v1++) CM[v1].resize(nvar);
			
	for(auto i = 0; i < nvar; i++) {
		for(auto j = 0; j <= i; j++) {
			auto sum = 0.0; for(auto k = 0; k < j; k++) sum += CM[i][k]*CM[j][k];

			if (i == j){
				auto val = M[i][i] - sum;
				if(val < -TINY){ CM[0][0] = UNSET; return CM;}
				else{
					if(val < 0) val = 0;
					CM[i][j] = sqrt(val);
				}
			}
			else{
				auto denom = CM[j][j];
				if(denom == 0){ CM[0][0] = UNSET; return CM;}
				CM[i][j] = ((1.0/denom)*(M[i][j] - sum));
			}
		}
	}
	
	if(false){
		for(auto j = 0; j < nvar; j++){
			for(auto i = 0; i < nvar; i++){
				auto sum = 0.0; for(auto k = 0; k < nvar; k++) sum += CM[j][k]*CM[i][k];
				cout << sum << " ";
			}
			cout << " Check "<< endl;
		}
	}
	
	return CM;
}

/// Performs LU decomposition to return the log of the determinant
double determinant_fast(const vector < vector <double> > &a) 
{
	auto n = a.size();  
	
	vector < vector <double> > l, u; // Performs LU decomposition (converts into a lower and an upper triagular matrix
	
	l.resize(n); u.resize(n);
	for(auto i = 0u; i < n; i++){ l[i].resize(n,0); u[i].resize(n,0);}
	
	auto i = 0u, j = 0u, k = 0u;
	for(i = 0; i < n; i++){
		auto &ui = u[i], &li = l[i];
		
		for(j = i; j < n; j++){
			auto &lj = l[j];
			auto val = a[j][i];
			for(k = 0; k < i; k++) val -= lj[k]*ui[k];
			lj[i] = val;
		}
		
		u[i][i] = 1;
		for(j = i+1; j < n; j++){
			auto &uj = u[j];
			
			if(li[i] == 0) emsg("Matrix9");
			auto fac = 1.0/li[i];
			auto val = fac*a[i][j];
			for(k = 0; k < i; k++) val -= fac*li[k]*uj[k];
			u[j][i] = val;
		}
	}
	
	auto det = 0.0;
	auto neg = 0u;
	for(i = 0; i < n; i++){
		auto valu = u[i][i]; if(valu < 0){ valu = -valu; neg++;}
		auto vall = l[i][i]; if(vall < 0){ vall = -vall; neg++;}
		det += log(valu*vall);
	}
	
	return det;
}

double sparsity(const vector < vector <double> > &a) 
{
	auto num_zero = 0.0, num = 0.0;
	for(auto v = 0; v < a.size(); v++){
		for(auto vv = 0; vv < a.size(); vv++){
			if(a[v][vv] == 0) num_zero++;
			num++;
		}
	}
	
	return num_zero/num;
}


/// Performs LU decomposition to return the log of the determinant
double determinant_sparse(const vector < vector <double> > &a) 
{
	auto n = a.size();  
	
	vector < vector <double> > l, u; // Performs LU decomposition (converts into a lower and an upper triagular matrix
	
	vector < vector <unsigned int> > l_ele, u_ele;
	l_ele.resize(n); u_ele.resize(n); 
	
	l.resize(n); u.resize(n);
	for(auto i = 0u; i < n; i++){ l[i].resize(n,0); u[i].resize(n,0);}
	
	auto i = 0u, j = 0u, k = 0u;
	for(i = 0; i < n; i++){
		auto &ui = u[i], &li = l[i];
		
		for(j = i; j < n; j++){
			auto &lj = l[j];
			auto &l_elej = l_ele[j];
			auto val = a[j][i];
			for(auto k : l_elej) val -= lj[k]*ui[k];
			lj[i] = val;
			if(val != 0) l_elej.push_back(i);
		}
		
		u[i][i] = 1;
		for(j = i+1; j < n; j++){
			auto &uj = u[j];
			auto &u_elej = u_ele[j];
			if(li[i] == 0) emsg("Matrix10");
			auto fac = 1.0/li[i];
			auto val = fac*a[i][j];
			for(auto k : u_elej) val -= fac*li[k]*uj[k];
			u[j][i] = val;
			if(val != 0) u_elej.push_back(i);
		}
	}
	
	//cout << sparsity(a) << " " << sparsity(u) << " " << sparsity(l) << " spar\n";

	auto det = 0.0;
	auto neg = 0u;
	for(i = 0; i < n; i++){
		auto valu = u[i][i]; if(valu < 0){ valu = -valu; neg++;}
		auto vall = l[i][i]; if(vall < 0){ vall = -vall; neg++;}
		det += log(valu*vall);
	}
	
	return det;
}


/// Transposes a matrix
vector < vector <double> > transpose(const vector < vector <double> > &M)
{
	vector < vector <double> > M_T;
	auto X = M.size(); if(X == 0) emsg("Matrix transpose problem"); 
	auto Y = M[0].size();
	
	M_T.resize(Y);
	for(auto j = 0; j < Y; j++){
		M_T[j].resize(X);
		for(auto i = 0; i < X; i++){
			M_T[j][i] = M[i][j];
		}
	}
	
	return M_T;
}


/// Takes the dot product of two vectors
double dot_prod(const vector <double> &vec1, const vector <double> &vec2)
{
	if(vec1.size() != vec2.size()) emsg("Dot product error");
	
	auto sum = 0.0;
	for(auto i = 0; i < vec1.size(); i++) sum += vec1[i]*vec2[i];
	
	return sum;
}
