#ifndef SIRE__MPI_HH
#define SIRE__MPI_HH

#include <vector>

using namespace std;

#include "const.hpp"
#include "struct.hpp"

struct Mpi {
public:
	unsigned int ncore;                                          // The number of cores that MPI is using
	unsigned int core;      

	Mpi();
	void bcast(double &val);
	void bcast(bool &val);
	void bcast(unsigned int &val);
	void bcast(vector <unsigned int> &vec);
	void bcast(vector <double> &vec);
	void bcast(vector <long> &vec);
	void bcast(vector < vector <double> > &M);
	void bcast(vector <DataSample> &ds);
	vector <VariableSample> gather_psamp(const vector <VariableSample> &psample);
	vector <double> gather(const vector <double> &vec);
	vector <double> gather(const double val);
	vector <long> gather(const long val);
	vector < vector <Sample> > gather_sample(const vector <Sample> &sample);
	vector < vector <IndPM> > gather_indPM(const vector <IndPM> &sample);
	void barrier();
	void mcmc_boostrap(vector <double> &param_value, vector <IndValue> &ind_value, vector <double> &L_ind_effect, vector <double> &L_inf_events, double &L_trans_events, vector <double> &L_diag_test, double &prior, const vector <long> &map);
	

private:
	vector<double> buffer;                                              // Stores packed up information to be sent between cores
	unsigned int k;                                                     // Indexes the buffer
	
	void pack_recv(const unsigned int co);
	void pack_initialise(const size_t size);
	void unpack_check();
	void pack(const vector <VariableSample> &vec);
	void unpack(vector<VariableSample> &vec);
	void pack(const vector<Sample> &vec);
	void unpack(vector<Sample> &vec);
	void pack(const vector<IndPM> &vec);
	void unpack(vector<IndPM> &vec);
	void pack_send(const unsigned int co);
	void pack(const unsigned int num);
	void unpack(unsigned int &num);
	void pack(const double num);
	void unpack(double &num);
	void pack(const bool num);
	void unpack(bool &num);
	void pack(const vector<double> &vec);
	void unpack(vector<double> &vec);

	double* packbuffer();
	size_t packsize();
};

#endif