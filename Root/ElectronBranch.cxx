#include <stdlib.h>

#include "MyXAODTools/ElectronBranch.h"
#include "CPAnalysisExamples/errorcheck.h"

const char* ElectronBranch::APP_NAME = "ElectronBranch";

ElectronBranch::ElectronBranch(){
    m_isBranchCreated = false;
}

bool ElectronBranch::CreateBranch()
{
    m_isBranchCreated = true;
    p4_ = new vector<TLorentzVector>();
    return true;
}

ElectronBranch::~ElectronBranch(){
    if(m_isBranchCreated){
        delete p4_;
    }
}

void ElectronBranch::ClearBranch(){
    total_ = 0;
    if(m_isBranchCreated){
        p4_->clear();
    }
}

void ElectronBranch::AttachBranchToTree(TTree& tree){
    if(!m_isBranchCreated){
        CreateBranch();
    }
    tree.Branch("n_ele", &total_, "n_ele/I");
    tree.Branch("ele_p4", &p4_);
}

void ElectronBranch::Fill(const xAOD::Electron& ele) 
{
    total_ ++;
    p4_->push_back(ele.p4());
}
