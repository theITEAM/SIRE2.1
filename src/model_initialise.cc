/// This provides code which takes the input xml file and sets the properties of the model class

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include "tinyxml2.h"

using namespace std;

using namespace tinyxml2;

#include "mcmc.hh"
#include "utils.hh"
#include "const.hh"


/// The constructor for the model class
Model::Model(string file)
{
	cout << "Initialising model..." << endl;

  load_input_file(file);
  
  output_model();
}


/// Loads the input xml file
void Model::load_input_file(string file)
{
  ifstream inp(file.c_str());
  stringstream strStream;
  strStream << inp.rdbuf();
  string str = strStream.str();
  if(str.size() == 0) emsg("No input file");
  
  remove_comments(str);
  
  char *xml = new char[str.size() + 1];
  xml[str.size()] = 0;
  memcpy(xml, str.c_str(), str.size());
  XMLDocument doc;
  doc.Parse(xml);
  
  vector <XMLNode *> child_list;
  for(auto child = doc.FirstChildElement()->FirstChild(); child; child = child->NextSibling()){
    child_list.push_back(child);
  }
  
  group_effect.on = false;              // By default the group effect is off
  
  // Adds parameters
  for(auto child : child_list){
    if(tab_name(child) == "prior") add_parameter(child);
  }
  nparam = param.size();
  
	// Gets basic properties of the MCMC chain
  for(auto child : child_list){
    auto val = tab_name(child);
    if(val == "mcmc"){
      output_dir = get(child,"output_dir");
      nsample = get_int(child,"nsample");
      nburnin = get_int(child,"burnin");
      nthin = 1; if(exist(child,"thin")) nthin = get_int(child,"thin");
    }
    
    if(val == "data"){
      if(exist(child,"N")) N = get_int(child,"N");
      if(exist(child,"Z")) ngroup = get_int(child,"Z");
    }
  }
  
  // Sets up epidemiological groups
  group.resize(ngroup);
  for(auto g = 0; g < ngroup; g++){
    group[g].name = "Group "+to_string(g);
    group[g].nind = 0;
    group[g].index = g;
  }
  
	// Finds observation and inference time ranges for the groups
  for(auto child : child_list){
    auto val = tab_name(child);
    if(val == "inference" || val == "observation") add_time_range(child);
  }
  
  // Add groups effects to the model
  for(auto child : child_list){
    if(tab_name(child) == "group_effect") switch_on_group_effect(get(child,"sigma"));
  }
    
  // Adds comparments to the model
  for(auto child : child_list){
    if(tab_name(child) == "comp") add_comp(child);
  }
  ncomp = comp.size();

  // Adds transitions to the model
  for(auto child : child_list){
    if(tab_name(child) == "trans") add_trans(child);
  }
  ntrans = trans.size();

  nind_effect = ind_effect.size();
  nfixed_effect = fixed_effect.size();
  nsnp_effect = snp_effect.size();
	
	// Initialise individual effect proposals
	ind_effect_proposal_initialise();
	
	// Adds diagnostic test to the model
  for(auto child : child_list){
    if(tab_name(child) == "diagnostic_test") add_diag_test(child);
  }
  ndiag_test = diag_test.size();
	
	// Adds information from the data table
  for(auto child : child_list){
    if(tab_name(child) == "datatable") add_datatable(child);
  }

  // Shifts fixed effects so they have an average of zero
  shift_fixed_effects();

  // Loads matrices
  if(nind_effect > 0) add_identity_matrix();
  for(auto child : child_list){
    auto val = tab_name(child);
    if(val == "A" || val == "A_nonzero" || val == "pedigree") add_matrix(child);
  }

  // Loads covariance matrices
  for(auto child : child_list){
    if(tab_name(child) == "covariance") add_covariance(child);
  }
  ncovariance = covariance.size();
  
  // Loads groups for calculating prediction accuracy
  for(auto child : child_list){
    if(tab_name(child) == "prediction_accuracy") add_pred_acc(child);
  }

  // Initialise proposals on event times
  ev_change_initialise();
  
  // Initialise add/remove infected individual proposals
  unknown_initialise();
	 
  // Checks the model is correctly specified
  check_model();
}


/// Add compartment to the model
void Model::add_comp(XMLNode *child)
{
  Comp co;
  co.name = get(child,"name");
  
  co.infectivity = 0;
  
  if(exist(child,"infectivity")){
    co.infectivity = get_num(child,"infectivity");
  
    if(exist(child,"individual_effect")){
      auto vec = split(get(child,"individual_effect"),',');
      for(auto v : vec){
        auto ie = add_ind_effect(v);
        ind_effect[ie].infectivity_comp.push_back(comp.size());
        co.ind_effect.push_back(ie);
      }
    }
    
    if(exist(child,"fixed_effect")){
      auto vec = split(get(child,"fixed_effect"),',');
      for(auto v : vec){
        auto fe = add_fixed_effect(v,"infectivity");
        co.fixed_effect.push_back(fe);
      }
    }
		
		if(exist(child,"snp_effect")){
      auto vec_mag = split(get(child,"snp_effect"),',');
			if(!exist(child,"snp_dominance")) emsg("There must be 'snp_dominance' and 'snp_effect'");
			auto vec_dom = split(get(child,"snp_dominance"),',');
			if(vec_mag.size() != vec_dom.size()) emsg("'snp_effect' and 'snp_dominance' must have the same size");
			
			for(auto j = 0u; j < vec_mag.size(); j++){ 	
        auto se = add_snp_effect(vec_mag[j],vec_dom[j],"infectivity");
        co.snp_effect.push_back(se);
      }
    }
  } 
    
  comp.push_back(co);
}


/// Add transition to the model
void Model::add_trans(XMLNode *child)
{
  Trans tr;
  
  auto from = get(child,"from");
  tr.from = find_comp(from);
  
  auto to = get(child,"to");
  tr.to = find_comp(to);
    
  tr.name = from+"->"+to;
  
  auto type = get(child,"type");
  if(type == "infection"){
    beta_param = find_param(get(child,"beta"),TRANS_INFRATE);
    auto inf_model_str = get(child,"inf_model");
    if(inf_model_str == "density dependent") inf_model = DENSITY_DEP;
    else{
      if(inf_model_str == "frequency dependent") inf_model = FREQ_DEP;
      else{
        emsg("'inf_model' must be set to 'density dependent' or 'frequency dependent'");
      }
    }
    
    tr.mean_param = UNSET;
    tr.shape_param = UNSET;
    tr.type = INFECTION;
  }
  else{
    if(type == "gamma"){
      tr.type = GAMMA;
      tr.mean_param = find_param(get(child,"mean"),TRANS_MEAN);
      tr.shape_param = find_param(get(child,"shape"),TRANS_SHAPE);
    }
  }
  
  if(exist(child,"individual_effect")){
    auto vec = split(get(child,"individual_effect"),',');
    for(auto v : vec){
      auto ie = add_ind_effect(v);
			if(trans.size() == 0) ind_effect[ie].susceptibility.push_back(trans.size());
			else ind_effect[ie].trans_mean.push_back(trans.size());
      tr.ind_effect.push_back(ie);
    }
  }
  
  if(exist(child,"fixed_effect")){
    auto vec = split(get(child,"fixed_effect"),',');
    for(auto v : vec){
      auto fe = add_fixed_effect(v,type);
      tr.fixed_effect.push_back(fe);
    }
  }
	
	if(exist(child,"snp_effect")){
		auto vec_mag = split(get(child,"snp_effect"),',');
		if(!exist(child,"snp_dominance")) emsg("There must be 'snp_dominance' and 'snp_effect'");
		auto vec_dom = split(get(child,"snp_dominance"),',');
		if(vec_mag.size() != vec_dom.size()) emsg("'snp_effect' and 'snp_dominance' must have the same size");
		
		for(auto j = 0u; j < vec_mag.size(); j++){ 	
			auto se = add_snp_effect(vec_mag[j],vec_dom[j],type);
			tr.snp_effect.push_back(se);
		}
	}
  
  tr.data_column = UNSET;
  if(exist(child,"data_column")){
    tr.data_column = get_int(child,"data_column")-1;
  }
  
  trans.push_back(tr);
}


/// Add diagnostic test to the model
void Model::add_diag_test(XMLNode *child)
{
	DiagnosticTest dt;
	dt.name = get(child,"name");
	
	// Sets which compartments are test sensitive
	dt.comp.resize(ncomp); for(auto c = 0u; c < ncomp; c++) dt.comp[c] = false;
	auto vec = split(get(child,"comp"),',');
	for(auto val : vec) dt.comp[find_comp(val)] = true;
	
	if(!exist(child, "sensitivity")) emsg("'sensitivity' must be specified for a diagnostic test");
	dt.Se_param = find_param(get(child,"sensitivity"),DIAG_TEST);
	
	if(!exist(child, "specificity")) emsg("'specificity' must be specified for a diagnostic test");
	dt.Sp_param = find_param(get(child,"specificity"),DIAG_TEST);
	
	if(!exist(child, "data_column")) emsg("'data_column' must be specified for a diagnostic test");
	dt.data_column = get_int(child,"data_column")-1;
	
	diag_test.push_back(dt);
}


/// Add parameter to the model
void Model::add_parameter(XMLNode *child)
{
  Param par;
  par.name = get(child,"parameter");
  par.type = NO_PARAMTYPE;
  
  auto type = get(child,"type");
  if(type == "Flat"){
    par.prior_type = FLAT_PRIOR;
    par.prior_val1 = get_num(child,"val1");
    par.prior_val2 = get_num(child,"val2");
  }
  else{
    emsg("Cannot find prior type '"+type+"'");
  }
  
  param.push_back(par);
}


/// Add time range to the model
void Model::add_time_range(XMLNode *child)
{
  auto tmin = get_num(child,"tmin");
  auto tmax = get_num(child,"tmax");
  auto gmin = 0, gmax = ngroup; 
  if(exist(child, "group")){
    gmin = int(get_int(child, "group")) - 1;
    gmax = gmin + 1;
  }
  for(auto g = gmin; g < gmax; g++){
    if(tab_name(child) == "inference"){
      group[g].inference_range.tmin = tmin; 
      group[g].inference_range.tmax = tmax; 
    }
    else{
      group[g].observation_range.tmin = tmin; 
      group[g].observation_range.tmax = tmax; 
    }
  }
}


/// Adds a simulated column to the datatable (used to generate simulated examples)
void Model::add_simulated_column(Table &tab) const
{
  for(auto r = 0u; r < tab.nrow; r++){
    tab.ele[r].push_back(to_string(int(ran()*2)));
  }
  tab.ncol++;
}


/// Add datatable
void Model::add_datatable(XMLNode *child)
{
  auto text = child->FirstChild()->Value();
  auto tab = create_table(text);
  
  //add_simulated_column(tab);
  
  auto id_column = get_int(child,"id")-1; 
  if(id_column < 0 || id_column >= tab.ncol) emsg("'id' column out of range");
  
  auto group_column = UNSET;
  if(ngroup > 1){
    group_column = get_int(child,"group")-1;
    if(group_column < 0 || group_column >= tab.ncol) emsg("'group' column out of range");
  }
  
  auto init_comp_column = get_int(child,"initial_comp")-1;
  if(init_comp_column < 0 || init_comp_column >= tab.ncol) emsg("'initial_comp' column out of range");
      
  auto comp_status_column = UNSET;
  if(exist(child,"comp_status")) comp_status_column = get_int(child,"comp_status")-1;
  
  vector <int> ind_effect_column(nind_effect);
  for(auto e = 0u; e < nind_effect; e++){
    if(exist(child,ind_effect[e].name)){
      ind_effect_column[e] = get_int(child,ind_effect[e].name)-1;
      if(ind_effect_column[e] < 0 || ind_effect_column[e] >= tab.ncol) emsg("'"+ind_effect[e].name+"' column out of range");
    }
    else ind_effect_column[e] = UNSET;
  }
    
  vector <int> fixed_effect_column(nfixed_effect);
  for(auto fe = 0u; fe < nfixed_effect; fe++){
    auto name = fixed_effect[fe].name;
    if(exist(child,name)){
      fixed_effect_column[fe] = get_int(child,name)-1;
      if(fixed_effect_column[fe] < 0 || fixed_effect_column[fe] >= tab.ncol) emsg("'"+name+"' column out of range");
    }
    else emsg("There must be a column specified for '"+name+"'");
  }
	
	vector <int> snp_effect_column(nsnp_effect);
  for(auto se = 0u; se < nsnp_effect; se++){
    auto name = snp_effect[se].name;
    if(exist(child,name)){
      snp_effect_column[se] = get_int(child,name)-1;
      if(snp_effect_column[se] < 0 || snp_effect_column[se] >= tab.ncol) emsg("'"+name+"' column out of range");
    }
    else emsg("There must be a column specified for '"+name+"'");
  }
	
  if(tab.nrow != N) emsg("The number of individuals does not agree with the size of the datatable");

  for(auto r = 0u; r < tab.nrow; r++){            // We go through each row in the datatable
    Individual ind;
    
    ind.id = tab.ele[r][id_column];               // The id for the individual
    
    ind.inside_group = true;                      // This determines if an individual belongs to a group
    ind.group = UNSET;
    ind.status = UNKNOWN;

    auto g = 0;                                   // Defines the group to which the individual belongs
    if(group_column != UNSET){
      auto val = tab.ele[r][group_column];
      if(val == "NA"){
        ind.inside_group = false;
        change_status(ind,NOT_INFECTED);
      }
      else{
        g = number(val)-1;
        if(g < 0 || g >= group.size()) emsg("In datatable group value is out of range"); 
        ind.group = g;
        group[g].ind_ref.push_back(individual.size());
        group[g].nind++;
      }
    } 

    auto &timer = ind.trans_time_range;                // This finds the time range consistent with observations in the datatable
    timer.resize(ntrans); 
    for(auto tr = 0u; tr < ntrans; tr++){
      timer[tr].tmin = UNSET;
      timer[tr].tmax = UNSET;
    }
    
    if(init_comp_column != UNSET){
      if(ind.inside_group == true){
        ind.initial_comp = find_comp(tab.ele[r][init_comp_column]);
        
        if(ind.initial_comp > 0){                       // If initially infected then set times
          auto initial_time = group[g].inference_range.tmin;
          for(auto tr = 0u; tr < ind.initial_comp; tr++){   
            timer[tr].tmin = initial_time;
            timer[tr].tmax = initial_time;
          }
          change_status(ind,INFECTED);
        }
      }
      else{
        ind.initial_comp = UNSET;
        if(tab.ele[r][init_comp_column] != "NA"){
          emsg("Cannot specify initial state for individual not in a group");
        }
      }
    }
    else{
      emsg("The 'init_comp' column must be spcified.");
    }

    /// Sets the time ranges over which transition can occur
    for(auto tr = 0u; tr < ntrans; tr++){
      if(ind.inside_group == true){
        auto data_column = trans[tr].data_column; 
        if(data_column != UNSET){
          if(data_column < 0 || data_column >= tab.ncol) emsg("'"+trans[tr].name+"' data column out of range");
        }
        
        if(timer[tr].tmin != UNSET){
          if(data_column != UNSET){                 // The data column provides information about this transition
            auto val = tab.ele[r][data_column];
            if(val != "no"){
              emsg("For individual '"+ind.id+"' the value for transition '"+trans[tr].name+"' in column "+to_string(data_column+1)+" is expected to be 'no'.");
            }
          }           
        }
        else{
          timer[tr].tmin = group[g].inference_range.tmin; // Default if no observation is made on transition  
          if(trans[tr].type == INFECTION) timer[tr].tmax = group[g].inference_range.tmax;
          else timer[tr].tmax = LARGE;
          
          if(data_column != UNSET){                       // The data column provides information about this transition
            auto val = tab.ele[r][data_column];
            if(val == "no"){                              // The transition is not observered during the observation range
              if(trans[tr].type == INFECTION){
                timer[tr].tmin = group[g].observation_range.tmax; 
                timer[tr].tmax = group[g].inference_range.tmax;
              }
              else{
                timer[tr].tmin = group[g].observation_range.tmax; 
                timer[tr].tmax = LARGE;
              }
            
              if(timer[tr].tmin == timer[tr].tmax) change_status(ind,NOT_INFECTED);   
            }
            else{                                          // The transition time is observed exactly     
              if(val == "."){                              // Observation time is missing
              }
              else{
                change_status(ind,INFECTED);
                
                auto t = number(val);
                if(t <= group[g].observation_range.tmin || t > group[g].observation_range.tmax){
                  emsg("Observed transition time of '"+val+"' must lie within the observation range");
                }
              
                timer[tr].tmin = t;
                timer[tr].tmax = t;
              }
            } 
          }
        }
      }
    }
  
		// In information comes about one transition this can limit the possible times for other transitions
    for(auto tr = 0u; tr < ntrans; tr++){
      if(ind.inside_group == true){ 
        if(ind.status == INFECTED){     
          auto tmin = timer[tr].tmin;
          if(tmin != UNSET){
            for(auto tr2 = tr+1; tr2 < ntrans; tr2++){           // Subsequent transitions must occur after this time
              if(timer[tr2].tmin < tmin) timer[tr2].tmin = tmin;
            }
          }
          
          auto tmax = timer[tr].tmax;
          if(tmax != UNSET){
            for(auto tr2 = 0u; tr2 < tr; tr2++){                 // Prior transitions must occur before this time
              if(timer[tr2].tmax > tmax) timer[tr2].tmax = tmax;
            }
          }
        }
      }
    }                 
    
		// This take and information provides about compartmental state (from the 'comp_status' column)
    if(comp_status_column != UNSET){
      auto val = tab.ele[r][comp_status_column];
      if(ind.inside_group == true){  
        bool infected_flag = false;
        bool not_infected_flag = false;
        
        auto time_max = -LARGE;
        auto vec = split(val,',');
        for(const auto &v : vec){
          if(v.length() == 0) emsg("There was a problem with the format of '"+val+"'");
          if(v.substr(0,1) != "[") emsg("There was a problem with the format of '"+val+"'");
          if(v.substr(v.length()-1,1) != "]") emsg("There was a problem with the format of '"+val+"'");
          auto str = v.substr(1,v.length()-2);
          
          auto spl = split(str,':');
          if(spl.size() != 2) emsg("There was a problem with the format of '"+val+"'");
          
          auto pos = split(spl[0],'|');
          auto cmin = ncomp-1, cmax = 0;
          for(auto p : pos){
            auto c = 0u; while(c < ncomp && p != comp[c].name) c++;
            if(c < cmin) cmin = c;
            if(c > cmax) cmax = c;
          }
          if(cmin > 0) infected_flag = true;
          if(cmax == 0) not_infected_flag = true; else not_infected_flag = false;
          
          auto t = number(spl[1]);
          if(t > time_max) time_max = t;
          
          for(auto tr = cmax; tr < ntrans; tr++){                       // Ensures that the transition occur after observation
            if(timer[tr].tmin < t){
              timer[tr].tmin = t; 
              if(timer[tr].tmin > timer[tr].tmax){
                emsg("Inconsistent data for individual '"+ind.id+"'");
              }
            }
          }
          
          for(auto tr = 0u; tr < cmin; tr++){                           // Ensures that the transition occurs before observation
            if(timer[tr].tmax > t){
              timer[tr].tmax = t; 
              if(timer[tr].tmin > timer[tr].tmax){
                emsg("Inconsistent data for individual '"+ind.id+"'");
              }
            }
          }
        }
        
        /// Set the infection status of the individual
        if(infected_flag == true) change_status(ind,INFECTED);  
        else{
          if(not_infected_flag == true) change_status(ind,NOT_INFECTED);
        }
      }
      else{
        if(val != "NA") emsg("Individuals not in a group must have 'comp_status' set to 'NA'.");
      }     
    }
    
		 // Gets information about individual effects (e.g. true breeding values)
    ind.ind_effect_value.resize(nind_effect);                      
    for(auto e = 0u; e < nind_effect; e++){
      auto col = ind_effect_column[e];
      if(col == UNSET) ind.ind_effect_value[e] = UNSET;
      else{
        ind.ind_effect_value[e] = number(tab.ele[r][col]);
      }
    }
    
		// Gets the design matrix for fixed effects
    ind.fixed_effect_X.resize(nfixed_effect);                
    for(auto fe = 0u; fe < nfixed_effect; fe++){
      ind.fixed_effect_X[fe] = number(tab.ele[r][fixed_effect_column[fe]]);
    }
    ind.fixed_effect_X_unshifted = ind.fixed_effect_X;
    
		// Gets the genotypes for SNPs
    ind.SNP_genotype.resize(nsnp_effect);                
    for(auto se = 0u; se < nsnp_effect; se++){
			auto geno = tab.ele[r][snp_effect_column[se]];
			if(geno == "AA") ind.SNP_genotype[se] = AA;
			else{
				if(geno == "AB") ind.SNP_genotype[se] = AB;
				else{
					if(geno == "BB") ind.SNP_genotype[se] = BB;
					else emsg("Genotype '"+geno+"' not recognised");
				}
			}
    }
	
		// Reads in any disease diagnostic test results
		ind.diag_test_result.resize(ndiag_test);
		for(auto dt = 0u; dt < ndiag_test; dt++){
			auto val = tab.ele[r][diag_test[dt].data_column];
			
			if(ind.inside_group == true){  	
				auto vec = split(val,',');
				for(const auto &v : vec){
					if(v.length() == 0) emsg("There was a problem with the format of '"+val+"'");
					if(v.substr(0,1) != "[") emsg("There was a problem with the format of '"+val+"'");
					if(v.substr(v.length()-1,1) != "]") emsg("There was a problem with the format of '"+val+"'");
					auto str = v.substr(1,v.length()-2);
						
					auto spl = split(str,':');
					if(spl.size() != 2) emsg("There was a problem with the format of '"+val+"'");
					DiagTestResult dtr;
					if(spl[0] == "+") dtr.positive = true;
					else{
						if(spl[0] == "-") dtr.positive = false;
						else emsg("The diagnostic test result '"+spl[0]+"' should be '+' or '-'");
					}
					dtr.time = number(spl[1]);
					
					ind.diag_test_result[dt].push_back(dtr);
				}
			}
			 else{
        if(val != "NA") emsg("Individuals not in a group must have test result '"+diag_test[dt].name+"' set to 'NA'.");
      }     
		}
		
		individual.push_back(ind);
  }
  nindividual = individual.size();
}


/// This shifts the design matrix of fixed effects so they have an average of zero
void Model::shift_fixed_effects()
{
  for(auto fe = 0u; fe < nfixed_effect; fe++){
    auto sum = 0.0;
    for(auto &indiv : individual) sum += indiv.fixed_effect_X[fe];
    sum /= N;
    
    for(auto &indiv : individual) indiv.fixed_effect_X[fe] -= sum;
  }
}


/// Adds an individual effect to the model
int Model::add_ind_effect(string st)
{
  auto i = 0u; while(i < ind_effect.size() && ind_effect[i].name != st) i++;
  if(i == ind_effect.size()){
    IndEffect indeff;
    indeff.name = st;
    indeff.covar_ref = UNSET;
    indeff.covar_num = UNSET;
    ind_effect.push_back(indeff);
  }
  return i;
}


/// Adds a fixed effect to the model
int Model::add_fixed_effect(string st, string type)
{
  auto i = 0u; while(i < fixed_effect.size() && fixed_effect[i].name != st) i++;
  if(i == fixed_effect.size()){
    FixedEffect fixeff;
    fixeff.name = st;
    fixeff.param = find_param(st,FIXED_EFFECT);
    fixeff.L_inf_update = false;
    fixeff.L_trans_update = false;
    fixed_effect.push_back(fixeff);
  }
  
  if(type == "infection" || type == "infectivity") fixed_effect[i].L_inf_update = true;
  else{
    if(type == "gamma") fixed_effect[i].L_trans_update = true;
    else emsg("Fixed effect type '"+type+"' not recognised");
  }
  
  return i;
}


/// Adds a snp effect to the model
int Model::add_snp_effect(string mag, string dom, string type)
{
  auto i = 0u; while(i < snp_effect.size() && snp_effect[i].name != mag) i++;
  if(i == snp_effect.size()){
    SNPEffect snpeff;
    snpeff.name = mag;
    snpeff.param_mag = find_param(mag,SNP_EFFECT_MAG);
		snpeff.param_dom = find_param(dom,SNP_EFFECT_DOM);
    snpeff.L_inf_update = false;
    snpeff.L_trans_update = false;
    snp_effect.push_back(snpeff);
  }
  
  if(type == "infection" || type == "infectivity") snp_effect[i].L_inf_update = true;
  else{
    if(type == "gamma") snp_effect[i].L_trans_update = true;
    else emsg("SNP effect type '"+type+"' not recognised");
  }
  
  return i;
}


/// Adds groups of indivduals for which we want to calculate prediction accuracies
void Model::add_pred_acc(XMLNode *child)
{
  PredAcc pa;
  pa.name = get(child,"name");
  auto vec = split(get(child,"ind"),',');

  for(const auto &id : vec){
    auto i = 0u; while(i < nindividual && individual[i].id != id) i++;
    if(i == nindividual) emsg("Individual '"+id+"' cannot be found.");
    pa.ind.push_back(i);
  }
  
  pred_acc.push_back(pa);
}


/// This initialises ev_change which provides a list of all the changes to the events which need to be made in a group
void Model::ev_change_initialise()
{
  // Calculates a characteristic timescale from the data
  auto av = 0.0, nav = 0.0;
  for(const auto &indiv : individual){
    if(indiv.status == INFECTED){
      auto g = indiv.group;
      auto tmin = group[g].inference_range.tmin;
      auto tmax = group[g].inference_range.tmax;
      
      for(auto tr = 0u; tr < ntrans; tr++){
        auto tmi = indiv.trans_time_range[tr].tmin;
        auto tma = indiv.trans_time_range[tr].tmax;
        if(tmi != UNSET && tma != UNSET){
          if(tmi > tmin && tma < tmax){
            av += (tmi+tma)/2 - tmin; nav++;
          }
        }         
      }
    }
  }
	
	if(nav == 0 && ndiag_test > 0){   // Attempts to get from disease diagnostic test results
		for(const auto &indiv : individual){
			if(indiv.inside_group == true){
				auto g = indiv.group;
				auto tmin = group[g].inference_range.tmin;
				auto tmax = group[g].inference_range.tmax;
				for(auto dt = 0u; dt < ndiag_test; dt++){
					for(const auto &dtr : indiv.diag_test_result[dt]){
						if(dtr.positive == true){ av += dtr.time - tmin; nav++;}
					}
				}
			}
		}
	}
	
  if(nav == 0) emsg("Cannot find characteristic time for epidemics");
  av /= nav;
  
	/// Depending on whether the time gaps are small / large compared to av we set either peaked or uniform proposals
  for(auto &gr : group){
    for(auto i : gr.ind_ref){
      const auto &indiv = individual[i];
      if(indiv.status != NOT_INFECTED){
        for(auto tr = 0u; tr < ntrans; tr++){
          const auto &tran = indiv.trans_time_range[tr];
          if(tran.tmin != tran.tmax){
            EvChange evch; evch.ind = i; evch.trans = tr;
            if(tran.tmax -  tran.tmin > av/3 && ntrans > 1) evch.ev_samp_type = PEAKED_SAMPLE;
            else evch.ev_samp_type = UNIFORM_SAMPLE;
            
            gr.ev_change.push_back(evch);
          }
        }
      }
    }
    gr.nev_change = gr.ev_change.size();
  }
}


/// This initialises the group property "unknown", which provides a list of all individuals with unknown disease status
void Model::unknown_initialise()
{
  for(auto &gr : group){
    for(auto i : gr.ind_ref){
      if(individual[i].status == UNKNOWN) gr.unknown.push_back(i);
    }
    gr.nunknown = gr.unknown.size();
  }   
}


/// Initialises how individual effects are updated
// For the infectious individual effects there are two types of update:
// 1) We go through each compartment and update all individual effects affecting that the infectiousness of that compartment
// 2) We go through each individual effect (much acts on multiple compartments) and update
// For the transition mean individual effects there are two types of update:
// 1) We go through each transition and update all individual effects affecting that the mean of that transition
// 2) We go through each individual effect (much acts on multiple transition) and update
// NOTE: individual effects cannot act on combinations of susceptibility/infectivity/transition mean	
void Model::ind_effect_proposal_initialise()
{
	for(auto &ie : ind_effect){                                       // Checks model specified correctly
		auto flag_sus = false; if(ie.susceptibility.size() > 0) flag_sus = true;
		auto flag_inf = false; if(ie.infectivity_comp.size() > 0) flag_inf = true;
		auto flag_trans = false; if(ie.trans_mean.size() > 0) flag_trans = true;
		if(flag_sus == true && flag_inf == true && flag_trans == true){
			emsg("Individual effect '"+ie.name+"' cannot act on susceptibility, infectivity and transition mean at the same time");
		}
		
		if(flag_sus == true && flag_inf == true){
			emsg("Individual effect '"+ie.name+"' cannot act on susceptibility and infectivity at the same time");
		}
		
		if(flag_sus == true && flag_trans == true){
			emsg("Individual effect '"+ie.name+"' cannot act on susceptibility and transition mean at the same time");
		}
		
		if(flag_inf == true && flag_trans == true){
			emsg("Individual effect '"+ie.name+"' cannot act on infectivity  and transition mean at the same time");
		}
		
		check_not_repeated(ie.susceptibility,ie.name);
		check_not_repeated(ie.infectivity_comp,ie.name);
		check_not_repeated(ie.trans_mean,ie.name);
	}
	
	for(auto &ie : ind_effect) ie.single_ind_effect_proposal = true;
	
  for(auto &co : comp){
		co.all_ind_effect_proposal = false;
		auto nie = co.ind_effect.size();
		if(nie > 0){
			auto ie = 0u; while(ie < nie && ind_effect[ie].infectivity_comp.size() == 1) ie++;
			if(ie == nie){
				co.all_ind_effect_proposal = true;
				for(auto ie :  co.ind_effect) ind_effect[ie].single_ind_effect_proposal = false;
			}
		}
	}
	
	for(auto &tra : trans){
		tra.all_ind_effect_proposal = false;
		auto nie = tra.ind_effect.size(); 
		if(nie > 0){
			if(tra.type == INFECTION) tra.all_ind_effect_proposal = true;
			else{
				auto ie = 0u; while(ie < nie && ind_effect[ie].trans_mean.size() == 1) ie++;
				if(ie == nie) tra.all_ind_effect_proposal = true;
			}
		}
		
		if(tra.all_ind_effect_proposal == true){
			for(auto ie :  tra.ind_effect) ind_effect[ie].single_ind_effect_proposal = false;
		}
	}
	
  if(false){	
		for(auto &ie : ind_effect){
			cout << ie.name << " ";
			if(ie.single_ind_effect_proposal == true) cout << "on" << endl;
			else cout << "off" << endl;
		}
		
		for(auto &co : comp){
			cout << co.name << " ";
			if(co.all_ind_effect_proposal == true) cout << "on" << endl; 
			else cout << "off" << endl;
		}
		 
		for(auto &tra : trans){
			cout << tra.name << " "; 
			if(tra.all_ind_effect_proposal == true) cout << "on" << endl; 
			else cout << "off" << endl;
		}
	}
}
