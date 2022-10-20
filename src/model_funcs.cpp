/// This useful function in the model class not directly related to the actual model definition

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include "tinyxml2.h"

using namespace std;

using namespace tinyxml2;

#include "model.hpp"
#include "utils.hpp"
#include "const.hpp"


/// Removes comments from a string
void Model::remove_comments(string &str) const {
	// cout << "Model::remove_comments()" << endl; // DEBUG
	auto i = 0;
	do {
		while (i < str.size() - 4 && str.substr(i, 4) != "<!--")
			++i;
		if (i == str.size() - 4)
			break;

		auto ist = i;
		i += 4;
		while (i < str.size() - 3 && str.substr(i, 3) != "-->")
			++i;
		if (i == str.size() - 3)
			emsg("Problem reading file");

		str = str.substr(0, ist) + str.substr(i + 3, str.size() - (i + 3));

		i = ist;
	} while (true);
}


/// Finds the compartment number from the name
int Model::find_comp(string st) const {
	// cout << "Model::find_comp()" << endl; // DEBUG
	for (auto c = 0; c < comp.size(); c++) {
		if (comp[c].name == st)
			return c;
	}
	emsg("Cannot find compartment '" + st + "'");
	return -1;
}


/// Finds the parameter number from the name
int Model::find_param(string st, ParamType type) {
	// cout << "Model::find_param()" << endl; // DEBUG
	for (auto th = 0; th < param.size(); th++) {
		if (param[th].name == st) {
			if (param[th].type != NO_PARAMTYPE && param[th].type != type)
				emsg("Parameter '" + st + "' type set more than once");
			param[th].type = type;
			return th;
		}
	}
	emsg("Cannot find parameter '" + st + "'");
	return -1;
}


/// Finds an individual from the id
int Model::find_ind(string id) const {
	// cout << "Model::find_ind()" << endl; // DEBUG
	for (auto i = 0; i < individual.size(); i++) {
		if (individual[i].id == id)
			return i;
	}
	if (id == ".")
		return UNSET;
	emsg("Cannot find individual '" + id + "'");
	return -1;
}


/// Returns the name of a node
string Model::tab_name(XMLNode *node) const {
	// cout << "Model::tab_name()" << endl; // DEBUG
	return node->Value();
}


/// Gets an attribute from an XML node
string Model::get(XMLNode *node, string attr) const {
	// cout << "Model::get()" << endl; // DEBUG
	if (node->ToElement()->Attribute(attr.c_str()))
		return node->ToElement()->Attribute(attr.c_str());
	else
		emsg("Cannot find '" + attr + "'");
	return "";
}


/// Gets a number attribute from an XML node
double Model::get_num(XMLNode *node, string attr) const {
	return number(get(node, attr));
}


/// Gets an integer attribute from an XML node
double Model::get_int(XMLNode *node, string attr) const {
	return atoi(get(node, attr).c_str());
}


/// Checks if an attibute exists on an XML node
bool Model::exist(XMLNode *node, std::string attr) const {
	if (node->ToElement()->Attribute(attr.c_str()))
		return true;
	else
		return false;
}


/// Creates a table from a string
Table Model::create_text_table(const string st) const {
	Table tab;

	tab.ncol = UNSET;

	stringstream ss(st);
	do {
		string line;
		getline(ss, line);
		if (ss.eof())
			break;

		auto vec = split(line, '\t');
		if (vec.size() != 0) {
			if (tab.ncol == UNSET)
				tab.ncol = vec.size();
			else {
				if (vec.size() != tab.ncol)
					emsg("Rows in datatable do not all share the same number of columns.");
			}

			tab.ele.push_back(vec);
		}
	} while (true);
	tab.nrow = tab.ele.size();

	return tab;
}


/// Loads a table from a file
Table Model::load_table(const string file, const bool head)
{
	if(file.length() < 4) emsg("File '"+file+"' must be '.csv', '.tsv' or '.txt'");
	char sep;
	auto end = file.substr(file.length()-4,4);
	if(end == ".csv") sep = ',';
	else{
		if(end == ".txt" || end == ".tsv" ) sep = '\t';
		else emsg("File '"+file+"' must be '.csv', '.tsv' or '.txt'");
	}
	
	ifstream in;

  auto used_file = file;
	in.open(used_file.c_str());
	if(!in) emsg("Cannot open the file '"+used_file+"'.");
	
	Table tab;
	tab.file = file;
	
	string line;
	do{
		getline(in,line);
	}while(line.substr(0,1) == "#");
		
	if(head == true){
		tab.heading = split(line,sep);
		tab.ncol = tab.heading.size();
	}
	else{
		tab.ncol = 0;
	}
	
	do{
		getline(in,line);
		if(in.eof()) break;

		auto vec = split(line,sep);
	
		if(tab.ncol == 0) tab.ncol = vec.size();
		else{
			if(vec.size() != tab.ncol) emsg("Rows in the file '"+file+"' do not all share the same number of columns.");
		}
		
		tab.ele.push_back(vec);
	}while(true);
	tab.nrow = tab.ele.size();
	
	return tab;
}


// Finds a column in a table (but doesn't return an error if can't find)
unsigned int Model::find_column(const Table &tab, const string name)
{
	for(auto c = 0u; c < tab.ncol; c++){
		if(tab.heading[c] == name) return c;
	}
		
	emsg("Cannot find column: "+name);
	return UNSET;	
}


/// Places a string in the format <----- string ---->
string Model::brackets(string st) const {
	const auto size = 80;
	int first = (size - 2 - st.size()) / 2;
	string s;
	s = "<";
	for (auto i = 0; i < first; i++)
		s += "-";
	s += " " + st + " ";
	for (auto i = 0; i < size - first - 2 - st.size(); i++)
		s += "-";
	s += ">";
	return s;
}


/// This outpus a summary of the model
void Model::output_model() const {
	struct stat st;
	if (stat(output_dir.c_str(), &st) == -1) {  // Directory not found
		int ret = mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (ret == -1)
			emsg("Error creating the directory '" + output_dir + "'");
	}

	// This generates a file that outputs a description of the model
	ofstream model_out(output_dir + "/model.txt");

	model_out << "# Individuals: " << N << endl;
	model_out << "# Groups: " << ngroup << endl;
	model_out << "# Samples: " << nsample << endl;
	model_out << "# Burnin: " << nburnin << endl;
	model_out << "Output Directory: " << output_dir << endl;

	model_out << endl;

	model_out << brackets("Compartments") << endl;
	for (auto co : comp) {
		model_out << co.name << "  ";
		if (co.infectivity != 0) {
			model_out << "Relative Infectivity: " << co.infectivity << "  ";
		}
		model_out << endl;
	}

	model_out << endl;

	if(infectivity.ind_effect.size() + infectivity.fixed_effect.size() 
		+ infectivity.snp_effect.size() > 0){		
		model_out << brackets("Infectivity")  << endl;
		if (infectivity.ind_effect.size() > 0) {
			model_out << "  Individual effect: ";
			for (auto i : infectivity.ind_effect)
				model_out << ind_effect[i].name + " ";
			model_out << endl;
		}

		if (infectivity.fixed_effect.size() > 0) {
			model_out << "  Fixed effect: ";
			for (auto i : infectivity.fixed_effect)
				model_out << fixed_effect[i].name + " ";
			model_out << endl;
		}

		if (infectivity.snp_effect.size() > 0) {
			model_out << "  SNP effect: ";
			for (auto i : infectivity.snp_effect)
				model_out << snp_effect[i].name + " ";
			model_out << endl;
		}
	}
		
	model_out << brackets("Transitions")  << endl;
	for (auto tr : trans) {
		model_out << tr.name << "   From: " << comp[tr.from].name << "  To: " << comp[tr.to].name << "  ";
		switch (tr.type) {
			case INFECTION:
				model_out << "Infection  beta: "  << param[beta_param].name;
				switch (inf_model) {
					case FREQ_DEP:
						model_out << "  Frequency dependent";
						break;
					case DENSITY_DEP:
						model_out << "  Density dependent";
						break;
				}
				break;

			case GAMMA:
				model_out << "Gamma  mean: " << param[tr.mean_param].name << "  shape: " << param[tr.shape_param].name;
				break;
			
			case EXP:
				model_out << "Exponential  mean: " << param[tr.mean_param].name;
				break;
		}
		model_out << "   ";

		if (tr.ind_effect.size() > 0) {
			model_out << "  Individual effect: ";
			for (auto ie : tr.ind_effect)
				model_out << ind_effect[ie].name << " ";
		}

		if (tr.fixed_effect.size() > 0) {
			model_out << "  Fixed effect: ";
			for (auto ie : tr.fixed_effect)
				model_out << fixed_effect[ie].name << " ";
		}

		if (tr.snp_effect.size() > 0) {
			model_out << "  SNP effect: ";
			for (auto ie : tr.snp_effect)
				model_out << snp_effect[ie].name << " ";
		}

		model_out << "   Data column: ";
		if (tr.data_column == "")
			model_out << "Unset";
		else
			model_out << tr.data_column;

		model_out << endl;
	}

	model_out << endl;

	model_out << brackets("Group effects") << endl;
	if (group_effect.on == false)
		model_out << "OFF";
	else
		model_out << "ON   Sigma: " <<  param[group_effect.sigma_param].name; ;
	model_out << endl;

	model_out << endl;

	if (nind_effect > 0) {
		model_out << brackets("Individual effects") << endl;
		for (auto ie : ind_effect) {
			model_out << ie.name << "   covariance matrix: " << ie.covar_ref << "     column/row:" << ie.covar_num << "     ";
			if (ie.infectivity_comp.size() > 0) {
				model_out << "Infectivity compartments: ";
				for (auto c : ie.infectivity_comp)
					model_out << comp[c].name << " ";
			}
			model_out << endl;
		}

		model_out << endl;
	}

	if (nfixed_effect > 0) {
		model_out << brackets("Fixed effects") << endl;
		for (auto fe : fixed_effect)
			model_out << fe.name << endl;

		model_out << endl;
	}

	if (nsnp_effect > 0) {
		model_out << brackets("SNP effects") << endl;
		for (auto se : snp_effect)
			model_out << se.name << endl;

		model_out << endl;
	}

	if (matrix.size() > 0) {
		model_out << brackets("Matrices") << endl;
		for (auto mat : matrix)
			model_out << mat.name << endl;

		model_out << endl;
	}

	if (ncovariance > 0) {
		model_out << brackets("Covariance matrices") << endl;
		for (auto cv : covariance) {
			model_out << "Matrix: " << matrix[cv.matrix].name << endl;
			model_out << "Individul effects: ";
			for (auto i : cv.ind_effect_ref)
				model_out << ind_effect[i].name << " ";
			model_out << endl;
			model_out << "Variance: ";
			for (auto i = 0; i < cv.E; i++)
				model_out << param[cv.var_param[i]].name << "  ";
			model_out << endl;
			model_out << "Correlation:" << endl;
			for (auto j = 0; j < cv.E; j++) {

				for (auto i = 0; i < cv.E; i++) {
					if (i == j)
						model_out << "1\t";
					else
						model_out << param[cv.cor_param[j][i]].name << "\t";
				}
				model_out << endl;
			}
			model_out << "All parameters used: ";
			for (auto th : cv.param_list)
				model_out << param[th].name << " ";
			model_out << endl;
			model_out << endl;
		}

		model_out << endl;
	}

	model_out << endl;

	if (ndiag_test > 0) {
		model_out << brackets("Diagnostic tests") << endl;

		for (auto dt : diag_test) {
			model_out << dt.name << "  ";
			model_out << "Test sensitive compartments: ";
			for (auto c = 0; c < ncomp; c++)
				if (dt.comp[c] == true)
					model_out << comp[c].name << " ";
			model_out << "   Se: " << param[dt.Se_param].name << "   Sp: " << param[dt.Sp_param].name << endl;
		}

		model_out << endl;
	}

	model_out << brackets("Parameters") << endl;

	for (auto par : param) {
		model_out << par.name << "  ";
		switch (par.type) {
			case TRANS_INFRATE:
				model_out << "Transmission rate";
				break;
			case TRANS_MEAN:
				model_out << "Transition mean";
				break;
			case TRANS_SHAPE:
				model_out << "Shape";
				break;
			case COVAR_MATRIX:
				model_out << "Covariance matrix";
				break;
			case GROUP_SD:
				model_out << "Group effect SD";
				break;
			case NO_PARAMTYPE:
				emsg("Parameter '" + par.name + "' is not used");
				break;
			// default:
				// cout << "Unknown parameter '" + par.name + "'" << endl;
				// break;
		}
		model_out << endl;
	}

	model_out << endl;

	if (pred_acc.size() > 0) {
		model_out << brackets("Prediction accuracy") << endl;
		for (auto pa : pred_acc) {
			model_out << pa.name << ":  ";
			for (auto i : pa.ind)
				model_out << individual[i].id << " ";
			model_out << endl << endl;
		}
		model_out << endl;
	}

	model_out << brackets("Groups") << endl;
	for (auto gr : group) {
		model_out << "< " << gr.name << " >" << endl;
		model_out << "Inference time range: ";
		model_out << gr.inference_range.tmin << " - " <<  num_out(gr.inference_range.tmax) << endl;
		model_out << "Observation time range: ";
		model_out << gr.observation_range.tmin << " - " <<  num_out(gr.observation_range.tmax) << endl;

		model_out << "Ind: ";
		for (auto i : gr.ind_ref)
			model_out << individual[i].id << " ";
		model_out << endl;
		model_out << endl;
	}

	model_out << brackets("Individuals") << endl;
	for (auto ind : individual) {
		model_out << ind.id << ":  ";
		if (ind.inside_group == false)
			model_out << "Not in a group";
		else {
			model_out << " Group: " << ind.group << "  ";

			model_out << " Initial compartment: " << comp[ind.initial_comp].name << "  ";

			switch (ind.status) {
				case NOT_INFECTED:
					model_out << "NOT INFECTED    ";
					break;

				case INFECTED:
					model_out << "INFECTED    ";
					break;
				case UNKNOWN:
					model_out << "UNKNOWN INFECTION STATUS    ";
					break;
			}

			if (ind.status != NOT_INFECTED) {
				for (auto tr = 0; tr < ntrans; tr++) {
					auto tmin = ind.trans_time_range[tr].tmin;
					auto tmax = ind.trans_time_range[tr].tmax;
					model_out << trans[tr].name << ": " << num_out(tmin);
					if (tmax != tmin)
						model_out << " - " << num_out(tmax);
					model_out << "     ";
				}
			}

			for (auto fe = 0; fe < nfixed_effect; fe++)
				model_out << fixed_effect[fe].name << ": " << ind.fixed_effect_X[fe] << "   ";

			for (auto se = 0; se < nsnp_effect; se++) {
				model_out << snp_effect[se].name << ": ";
				switch (ind.SNP_genotype[se]) {
					case AA:
						model_out << "AA ";
						break;
					case AB:
						model_out << "AB ";
						break;
					case BB:
						model_out << "BB ";
						break;
				}
			}
		}
		model_out << endl;
	}
}


/// This checks the model is correctly specified
void Model::check_model() const {
	if (ntrans != ncomp - 1)
		emsg("For " + to_string(ncomp) + " compartments " + to_string(ncomp - 1) + " transitions should be specified");

	if (trans[0].type != INFECTION)
		emsg("The first transition much be an infection");

	for (auto tr = 0; tr < ntrans; tr++) {
		if (trans[tr].from != tr)
			emsg("The transitions must go from one compartment to the next");
		if (trans[tr].to != tr + 1)
			emsg("The transitions must go from one compartment to the next");
	}

	for (const auto &gr : group) {
		if (gr.inference_range.tmin != gr.observation_range.tmin)
			emsg("In the current version the observations and inference must start at the same time.");
	}

	for (const auto &ie : ind_effect) {
		if (ie.covar_ref == UNSET)
			emsg("The individual effect '" + ie.name + "' must be specified in 'covariance'");
	}
}


/// Changes the infection status of an individual
void Model::change_status(Individual &ind, const Status stat) {
	if (ind.status != UNKNOWN) {
		if (ind.status != stat)
			emsg("Inconsistent disease status for individual '" + ind.id + "'");
	}

	ind.status = stat;
}
