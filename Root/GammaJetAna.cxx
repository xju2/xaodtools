#include <stdlib.h>

#include "MyXAODTools/GammaJetAna.h"
#include "CPAnalysisExamples/errorcheck.h"


GammaJetAna::GammaJetAna():AnalysisBase()
{
    if(APP_NAME==NULL) APP_NAME = "GammaJetAna";
    string maindir(getenv("ROOTCOREBIN"));
    m_susy_config = Form("%s/data/MyXAODTools/gamma_jet.conf", maindir.c_str());

    trigger_map_ = {
        {"HLT_g120_loose", false},
        {"HLT_g200_etcut", false}
    };
    CreateBranch();
    AttachBranchToTree();
}


GammaJetAna::~GammaJetAna(){
}

void GammaJetAna::CreateBranch()
{
    return ;
}

void GammaJetAna::ClearBranch(){
    AnalysisBase::ClearBranch();
}


void GammaJetAna::AttachBranchToTree()
{
    AnalysisBase::AttachBranchToTree();

    event_br->AttachBranchToTree(*physics);
    // jet_br->AttachBranchToTree(*physics);
}

int GammaJetAna::process(Long64_t ientry)
{
    int sc = AnalysisBase::process(ientry);
    if(m_debug) {
        Info(APP_NAME, " GammaJetAna: processing");
    }
    event_br->Fill(*ei);
    if(sc==0){
        tree->Fill();
        physics->Fill();
    }
    return sc;
}
