#include <stdlib.h>

#include "MyXAODTools/MuonBranch.h"
#include "CPAnalysisExamples/errorcheck.h"

const char* MuonBranch::APP_NAME = "MuonBranch";

MuonBranch::MuonBranch(){
    CreateBranch();
}

bool MuonBranch::CreateBranch()
{
    p4_ = new vector<TLorentzVector>();
    return true;
}

MuonBranch::~MuonBranch(){
    delete p4_;
}

void MuonBranch::ClearBranch(){
    total_ = 0;
    p4_->clear();
}

void MuonBranch::AttachBranchToTree(TTree& tree){
    tree.Branch("n_muon", &total_, "n_muon/I");
    tree.Branch("mu_p4", &p4_);
}

void MuonBranch::Fill(const xAOD::Muon& muon)
{
    total_ ++;
    p4_->push_back(muon.p4());
}
