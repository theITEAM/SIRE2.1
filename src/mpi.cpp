/// Information and routines for transferring data between cores using MPI

#include <sstream>
#include <mpi.h>

using namespace std;

#include "mpi.hpp"
#include "utils.hpp"

Mpi::Mpi()
{
	int num;
	MPI_Comm_size(MPI_COMM_WORLD,&num); ncore = (unsigned int) num;
  MPI_Comm_rank(MPI_COMM_WORLD,&num); core = (unsigned int) num;
}

/// Copies a variable in core 0 to all the other cores
void Mpi::bcast(double &val)
{
	MPI_Bcast(&val,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}


/// Copies a variable in core 0 to all the other cores
void Mpi::bcast(bool &val)
{
	MPI_Bcast(&val,1,MPI_CXX_BOOL,0,MPI_COMM_WORLD);
}


/// Copies a variable in core 0 to all the other cores
void Mpi::bcast(unsigned int &val)
{
	MPI_Bcast(&val,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
}


/// Copies a vector in core 0 to all the other cores
void Mpi::bcast(vector <unsigned int> &vec)
{
	auto n = (unsigned int)(vec.size());
	bcast(n);
	vec.resize(n);
	MPI_Bcast(&vec[0],vec.size(),MPI_UNSIGNED,0,MPI_COMM_WORLD);
}
	

/// Copies a vector in core 0 to all the other cores
void Mpi::bcast(vector <double> &vec)
{
	auto n = (unsigned int)(vec.size());
	bcast(n);
	vec.resize(n);
	MPI_Bcast(&vec[0],vec.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}


/// Copies a vector in core 0 to all the other cores
void Mpi::bcast(vector <long> &vec)
{
	auto n = (unsigned int)(vec.size());
	bcast(n);
	vec.resize(n);
	MPI_Bcast(&vec[0],vec.size(),MPI_LONG,0,MPI_COMM_WORLD);
}


/// Copies a matrix in core 0 to all the other cores
void Mpi::bcast(vector < vector <double> > &M)
{
	auto n = (unsigned int)(M.size());
	bcast(n);
	M.resize(n);
	for(auto j = 0; j < n; j++) bcast(M[j]);
}


/// Copies data_sample in core 0 to all the other cores
void Mpi::bcast(vector <DataSample> &ds)
{
	auto n = (unsigned int)(ds.size());
	bcast(n);
	ds.resize(n);
	
	auto nv = (unsigned int)(ds[0].em_var_value.size());
	bcast(nv);
	
	for(auto j = 0; j < n; j++){
		ds[j].em_var_value.resize(nv);
		bcast(ds[j].em_var_value);
		bcast(ds[j].post);
	}
}


/// Gathers together results parameter samples from particles all cores onto core 0
vector <VariableSample> Mpi::gather_psamp(const vector <VariableSample> &psample)
{
	vector <VariableSample> psamp;
	
	if(core == 0){
		psamp = psample;
		
		for(auto co = 1u; co < ncore; co++){
			vector <VariableSample> psa;
			pack_recv(co);
			unpack(psa);
			for(const auto &ps : psa) psamp.push_back(ps);
			unpack_check();
		}
	}
	else{
		pack_initialise(0);
		pack(psample);
		pack_send(0);
	}
	
	return psamp;
}


/// Gathers samples from different chains
vector < vector <Sample> > Mpi::gather_sample(const vector <Sample> &sample)
{
	vector < vector <Sample> > samp;
	
	if(core == 0){
		samp.push_back(sample);
		
		for(auto co = 1u; co < ncore; co++){
			vector <Sample> psa;
			pack_recv(co);
			unpack(psa);
			samp.push_back(psa);
			unpack_check();
		}
	}
	else{
		pack_initialise(0);
		pack(sample);
		pack_send(0);
	}
	
	return samp;
}


vector < vector <IndPM> > Mpi::gather_indPM(const vector <IndPM> &sample)
{
	vector < vector <IndPM> > samp;
	
	if(core == 0){
		samp.push_back(sample);
		
		for(auto co = 1u; co < ncore; co++){
			vector <IndPM> psa;
			pack_recv(co);
			unpack(psa);
			samp.push_back(psa);
			unpack_check();
		}
	}
	else{
		pack_initialise(0);
		pack(sample);
		pack_send(0);
	}
	
	return samp;
}


/// Gathers a double vector across all cores and returns the combined vector to core 0
vector <double> Mpi::gather(const vector <double> &vec)
{
	vector <double> vectot(vec.size()*ncore);
	
	MPI_Gather(&vec[0],vec.size(),MPI_DOUBLE,&vectot[0],vec.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	return vectot;
}


/// Gathers a long across all cores and returns the combined vector to core 0
vector <double> Mpi::gather(const double val)
{
	vector <double> valtot;
	valtot.resize(ncore);
	
	MPI_Gather(&val,1,MPI_DOUBLE,&valtot[0],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	return valtot;
}


/// Gathers a long across all cores and returns the combined vector to core 0
vector <long> Mpi::gather(const long val)
{
	vector <long> valtot;
	valtot.resize(ncore);
	
	MPI_Gather(&val,1,MPI_LONG,&valtot[0],1,MPI_LONG,0,MPI_COMM_WORLD);
	
	return valtot;
}


/// Recieves a message from core co and places it into the buffer
void Mpi::pack_recv(const unsigned int co)
{
	unsigned int si;
	MPI_Recv(&si,1,MPI_UNSIGNED,co,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	pack_initialise(si);
	MPI_Recv(packbuffer(),si,MPI_DOUBLE,co,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}	

/// Initialises the buffer
void Mpi::pack_initialise(const size_t size)
{
	k = 0;
	buffer.resize(size);
}

/// Checks all the values have been read from the buffer
void Mpi::unpack_check()
{
	if(k != buffer.size()) emsg("Mpi");
}

void Mpi::pack(const vector<VariableSample> &vec)
{
	unsigned int imax = vec.size();
	pack(imax);
	for(auto i = 0u; i < imax; i++){
		pack(vec[i].value);
		pack(vec[i].L);
	}
}

void Mpi::unpack(vector<VariableSample> &vec)
{
	unsigned int imax; 
	unpack(imax); vec.resize(imax);
	for(auto i = 0u; i < imax; i++){
		unpack(vec[i].value);
		unpack(vec[i].L);
	}
}

void Mpi::pack(const vector<Sample> &vec)
{
	unsigned int imax = vec.size();
	pack(imax);
	for(auto i = 0u; i < imax; i++){
		pack(vec[i].param_value);
	}
}

void Mpi::unpack(vector<Sample> &vec)
{
	unsigned int imax; 
	unpack(imax); vec.resize(imax);
	for(auto i = 0u; i < imax; i++){
		unpack(vec[i].param_value);
	}
}

void Mpi::pack(const vector<IndPM> &vec)
{
	unsigned int imax = vec.size();
	pack(imax);
	for(auto i = 0u; i < imax; i++){
		pack(vec[i].ind_effect_sum);
		pack(vec[i].ind_effect_sum2);
	}
}

void Mpi::unpack(vector<IndPM> &vec)
{
	unsigned int imax; 
	unpack(imax); vec.resize(imax);
	for(auto i = 0u; i < imax; i++){
		unpack(vec[i].ind_effect_sum);
		unpack(vec[i].ind_effect_sum2);
	}
}

/// Sends the buffer to core co
void Mpi::pack_send(const unsigned int co)
{
	unsigned int si = packsize();
	MPI_Send(&si,1,MPI_UNSIGNED,co,0,MPI_COMM_WORLD);
	MPI_Send(packbuffer(),si,MPI_DOUBLE,co,0,MPI_COMM_WORLD);
}


/// The pointer to the buffer
double* Mpi::packbuffer()
{
	return &buffer[0];
}

void Mpi::pack(const bool num)
{
	buffer.push_back(num); k++;
}

void Mpi::unpack(bool &num)
{
	num = buffer[k]; k++;
}

void Mpi::pack(const unsigned int num)
{
	buffer.push_back(num); k++;
}

void Mpi::unpack(unsigned int &num)
{
	num = buffer[k]; k++;
}

void Mpi::pack(const double num)
{
	buffer.push_back(num); k++;
}

void Mpi::unpack(double &num)
{
	num = buffer[k]; k++;
}


void Mpi::pack(const vector<double> &vec)
{
	pack((unsigned int)(vec.size()));
	for (auto &item : vec) pack(item);
}

void Mpi::unpack(vector<double> &vec)
{
	unsigned int size;
	unpack(size);
	vec.resize(size);
	for (auto &item : vec) unpack(item);
}

/// Returns the size of the buffer
size_t Mpi::packsize()
{
	return k;
}


/// Calculates the time taken for other cores to finish what they are doing
void Mpi::barrier()
{
	MPI_Barrier(MPI_COMM_WORLD);  
}


/// Swaps infomation between mcmc chains when doing a bootstrap step
void Mpi::mcmc_boostrap(vector <double> &param_value, vector <IndValue> &ind_value, vector <double> &L_ind_effect, vector <double> &L_inf_events, double &L_trans_events, vector <double> &L_diag_test, double &prior, const vector <long> &map)
{
	if(map[core] == core){
		vector <unsigned int> send;
		for(auto c = 0; c < ncore; c++){
			if(c != core && map[c] == core) send.push_back(c);
		}
		
		if(send.size() > 0){
			pack_initialise(0);
			pack(param_value);
			for(auto i = 0; i < ind_value.size(); i++){
				const auto &ind = ind_value[i];
				pack(ind.infected);
				pack(ind.trans_time);
				pack(ind.ind_effect);
				pack(ind.susceptibility);
				pack(ind.inf_single);
				pack(ind.trans_infectivity_change);
				pack(ind.trans_mean);
			}
			pack(L_ind_effect);
			pack(L_inf_events);
			pack(L_trans_events);
			pack(L_diag_test);
			pack(prior);
			//for(auto c : send) cout << core << "  send to "<< c << "\n"; 
			for(auto c : send) pack_send(c);
		}
	}
	else{
		//cout << core << " recieve from " << map[core] << endl;
		pack_recv(map[core]);
		
		unpack(param_value);
		for(auto i = 0; i < ind_value.size(); i++){
			auto &ind = ind_value[i];
			unpack(ind.infected);
			unpack(ind.trans_time);
			unpack(ind.ind_effect);
			unpack(ind.susceptibility);
			unpack(ind.inf_single);
			unpack(ind.trans_infectivity_change);
			unpack(ind.trans_mean);
		}
		unpack(L_ind_effect);
		unpack(L_inf_events);
		unpack(L_trans_events);
		unpack(L_diag_test);
		unpack(prior);
		
		unpack_check();
	}
}

