#include <stdlib.h>

#include "MyXAODTools/ElectronBranch.h"
#include "CPAnalysisExamples/errorcheck.h"

const char* ElectronBranch::APP_NAME = "ElectronBranch";

ElectronBranch::ElectronBranch(){
    CreateBranch();
}

bool ElectronBranch::CreateBranch()
{
    p4_ = new vector<TLorentzVector>();
    return true;
}

ElectronBranch::~ElectronBranch(){
    delete p4_;
}

void ElectronBranch::ClearBranch(){
    total_ = 0;
    p4_->clear();
}

void ElectronBranch::AttachBranchToTree(TTree& tree){
    tree.Branch("n_ele", &total_, "n_ele/I");
    tree.Branch("ele_p4", &p4_);
}

void ElectronBranch::Fill(const xAOD::Electron& ele) 
{
    total_ ++;
    p4_->push_back(ele.p4());
}
