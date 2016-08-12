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
    id_tag_ = new vector<int> ();
    return true;
}

MuonBranch::~MuonBranch(){
    delete p4_;
    delete id_tag_;
}

void MuonBranch::ClearBranch(){
    total_ = 0;
    p4_->clear();
    id_tag_->clear();
}

void MuonBranch::AttachBranchToTree(TTree& tree){
    tree.Branch("n_muon", &total_, "n_muon/I");
    tree.Branch("mu_p4", &p4_);
    tree.Branch("mu_id", &id_tag_);
}

void MuonBranch::Fill(const xAOD::Muon& muon)
{
    total_ ++;
    p4_->push_back(muon.p4());
}
