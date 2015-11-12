//////////////////////////////////////////////////
//
//
//////////////////////////////////////////////////
#include "MyXAODTools/EventCounter.h"
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <ctype.h>
#include <fstream>
#include <sstream>
#include <memory>
#include <exception>

using namespace std;

EventCounter::EventCounter()
{
    is_dc14_ = false;
    is_mc_ = true;
    zee_mc_id_vec = NULL;
    GetCrossSection();
}

EventCounter::EventCounter(bool is_dc14):is_dc14_(is_dc14)
{
    zee_mc_id_vec = NULL;
    is_mc_ = false;
    GetCrossSection();
}

EventCounter::EventCounter(bool is_dc14, bool is_mc):
    is_dc14_(is_dc14),
    is_mc_(is_mc)
{
    zee_mc_id_vec = NULL;
    GetCrossSection();
}

void EventCounter::GetCrossSection()
{
    // cross section in *pb*
    if(cross_section_dic_.size() > 0) return;
    if(!is_mc_) return;

    root_core_path = string(getenv("ROOTCOREBIN"));
    if(root_core_path.compare("") == 0){
        cout<<"ERROR: ROOTCOREBIN is not setup!"<<endl;
    }
    string xs_name = "susy_crosssections_13TeV.txt";
    if(is_dc14_){
        xs_name = "susy_crosssections_13TeV_dc14.txt";
    }
    string full_path(root_core_path+"/data/MyXAODTools/"+xs_name);
    cout<<"Loading cross section from: "<< full_path<<endl;
    // background 
    this->initial_cs(full_path.c_str());
    //manually add mono-jet DM signal
    cross_section_dic_[191040] = 5.4513E-03 * 2.7437E-01 * 1e3; //D5, 400, MET>100
    cross_section_dic_[191041] = 5.4498E-03 * 6.1546E-02 * 1e3; //D5, 400, MET>300
    cross_section_dic_[191042] = 5.4552E-03 * 2.0209E-02 * 1e3; //D5, 400, MET>500
    cross_section_dic_[191043] = 1.3001E-02 * 1.8771E-01 * 1e3; //D5, 50, MET>100
    // add mono-Higgs signal
    LoadXS(Form("%s/data/MyXAODTools/xs_zprime.txt", 
                root_core_path.c_str()));
    LoadXS(Form("%s/data/MyXAODTools/xs_scalar.txt", 
                root_core_path.c_str()));
    LoadXS(Form("%s/data/MyXAODTools/xs_sm.txt", 
                root_core_path.c_str()));
    LoadXS(Form("%s/data/MyXAODTools/xs_monojet_signal.txt", 
                root_core_path.c_str()));
    cout<<"Cross section DB is loaded"<<endl;
}

EventCounter::~EventCounter()
{
    if(zee_mc_id_vec != NULL) delete zee_mc_id_vec;
}

TChain* EventCounter::getChain(const char* filename, const char* chainName)
{
    std::string fileName(filename);
	TChain* fc = new TChain(chainName);
    if(fileName.find("root") != string::npos){
        fc->AddFile(filename);
    }else{
        std::ifstream input(filename, std::ifstream::in);
        TString name;
        while ( input >> name){
            fc->AddFile(name.Data());
        }
    }
    return fc;
}

void EventCounter::GetTotalEventsDic(const char* file_name, const char* tree_name)
{
    if(! total_events_dic_.empty())  total_events_dic_.clear();
    if(! all_events_noweight_dic_.empty()) all_events_noweight_dic_.clear();

    TChain* associate = getChain(file_name, tree_name);
    if(is_mc_){
        int mc_channel_number; Long64_t total_entries;
        float mcweight ;
        uint64_t n_events_processed;
        double sum_of_weights;
        associate->SetBranchAddress("mc_channel_number", &mc_channel_number);
        associate->SetBranchAddress("MCWeight", &mcweight);
        associate->SetBranchAddress("nEventsProcessed", &n_events_processed);
        associate->SetBranchAddress("nSumEventWeights", &sum_of_weights);

        total_entries = associate->GetEntries();
        for(Long64_t ientry = 0; ientry < total_entries; ientry++){
            associate->LoadTree(ientry);
            associate->GetEntry(ientry);
            try{
                total_events_dic_.at(mc_channel_number) += sum_of_weights;
                all_events_noweight_dic_.at(mc_channel_number) += n_events_processed;
            }catch (const std::out_of_range& oor){
                total_events_dic_[mc_channel_number] = sum_of_weights;
                all_events_noweight_dic_[mc_channel_number] = n_events_processed;
            }
        }
    } else {
        int run_number; 
        Long64_t total_entries;
        associate->SetBranchAddress("RunNumber", &run_number);
        total_entries = associate->GetEntries();
        for(Long64_t ientry = 0; ientry < total_entries; ientry++){
            associate->LoadTree(ientry);
            associate->GetEntry(ientry);
            try{
                all_events_noweight_dic_.at(run_number) += 1.0;
                total_events_dic_.at(run_number) += 1.0;
            }catch (const std::out_of_range& oor){
                all_events_noweight_dic_[run_number] = 1.0;
                total_events_dic_[run_number] = 1.0;
            }
        }
    }
    delete associate;
}

double EventCounter::getCrossSection(int mc_channel_number)
{
    return this->getUtil(mc_channel_number, this->cross_section_dic_, "cross section DB");
}

void EventCounter::initial_cs(const char* filename)
{
    string line;

    std::ifstream in(filename);
    if (!in) return;
    while ( getline(in,line) ) 
    {
        // skip leading blanks (in case there are some in front of a comment)
        if ( !line.empty() )
        {
            while ( line[0] == ' ' ) line.erase(0,1);
        }
        // skip lines that do not start with a number, they are comments
        if ( !line.empty() && isdigit(line[0]) )
        {
            std::stringstream is(line);
            int id;
            string name;
            float xsect, kfactor, efficiency, relunc;
            is >> id >> name >> xsect >> kfactor >> efficiency >> relunc;
            cross_section_dic_[id] = xsect*kfactor*efficiency*relunc;
        }
    }
}

double EventCounter::getTotalEvents(int mc_channel_number){
    return this->getUtil(mc_channel_number, this->total_events_dic_, "total events DB");
}

double EventCounter::getTotalEventsNoWeight(int mc_channel_number){
    return this->getUtil(mc_channel_number, this->all_events_noweight_dic_, "total events(no Weight) DB");
}

double EventCounter::getUtil(int mc_channel_number, std::map<int, double>& dic, const char* name)
{
    
    double total = 0.0;
    try{
        total = dic.at(mc_channel_number);
    }catch(const std::out_of_range& oor){
        cout<<"Error: "<< mc_channel_number<< " not exist in "<<name<<endl;
        throw std::out_of_range("CHECK DICTIONARY");
    }
    return total;
}

void EventCounter::printTotalEvents(){
    for(auto& x: total_events_dic_){
        if(is_mc_) cout<< x.first <<" : "<< x.second <<" : "<<this->getTotalEventsNoWeight(x.first)<<" : "<<this->getCrossSection(x.first)<< "\n";
        else cout<< x.first <<" : "<< x.second <<" : "<<this->getTotalEventsNoWeight(x.first)<<endl;
    }
}

void EventCounter::SaveTotalEvents(const char* file_dir){
    string in_dir;
    if (file_dir == NULL){
        in_dir = root_core_path+"/data/StackProcess/mc_weights.txt";
    } else {
        in_dir = string(file_dir);
    }
    fstream output_text(in_dir.c_str(), fstream::out);
    // output_text << "mc_id : totalEventsWithWeight : totalEventsNoWeight : crossSection" <<endl;
    for(auto& x: total_events_dic_){
        output_text << x.first <<" : "<< x.second <<" : "<<this->getTotalEventsNoWeight(x.first)<<" : "<<this->getCrossSection(x.first)<< "\n";
    }
    output_text.close();
}

bool EventCounter::isZeeSample(int mc_id){
    if(zee_mc_id_vec == NULL){
        zee_mc_id_vec = new vector<int>();
        // load the sample list when needed.
        string _path(root_core_path+"/data/StackProcess/Zee.list");
        std::cout<<"read Zee info "<< _path.c_str() << std::endl;
        std::ifstream input(_path, std::ifstream::in);
        int mc_id;
        while (!input.eof() && input.good()){
            input >> mc_id;
            zee_mc_id_vec ->push_back(mc_id);
        }
    }
    bool has_id = false;
    for(unsigned int i = 0; i < zee_mc_id_vec->size(); i++){
        if(zee_mc_id_vec->at(i) == mc_id){
            has_id = true;
            break;
        }
    }
    return has_id;
}

void EventCounter::ReadEventInfo(const char* input_name)
{
    if (input_name == NULL){
        printf("[ERROR] (%s:%d:) Input point is NULL ", __FILE__, __LINE__);
        throw std::invalid_argument("Give a Weight File!");
    }
    fstream in_file(input_name, fstream::in);
    do {
        int id;
        char c_tmp;
        double events, events_no_weight, xs;
        in_file >> id >> c_tmp >> events >> c_tmp >> events_no_weight >> c_tmp >> xs;
        total_events_dic_[id] = events;
        all_events_noweight_dic_[id] = events_no_weight;
        cross_section_dic_[id] = xs;
    } while (!in_file.eof() && in_file.good());
    in_file.close();
}

void EventCounter::LoadXS(const char* in_name)
{
    cout << "reading: " << in_name << endl;
    fstream in_file(in_name, fstream::in);
    string line;
    while( getline(in_file, line) ){
        if( line[0] == '#' ){
            continue;
        } else {
            istringstream foo(line);
            int id; 
            float tmp;
            double xs;
            foo >> id >> tmp >> tmp >> xs;
            cross_section_dic_[id] = xs;
        }
    }
    in_file.close();
}
