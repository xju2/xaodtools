#include <stdlib.h>

#include "MyXAODTools/UpsilonBranch.h"
#include "CPAnalysisExamples/errorcheck.h"

const char* UpsilonBranch::APP_NAME = "UpsilonBranch";

UpsilonBranch::UpsilonBranch(){
    trigger_map_ = {
        // single electron
        {"HLT_e24_lhmedium_L1EM18VH", false},
        {"HLT_e24_lhmedium_L1EM20VH", false},
        {"HLT_e24_lhtight_nod0_ivarloose", false},
        {"HLT_e26_lhtight_nod0_ivarloose", false},
        // di-electrons
        {"HLT_2e12_lhloose_L12EM10VH", false},
        {"HLT_2e15_lhvloose_nod0_L12EM13VH", false},
        {"HLT_2e17_lhvloose_nod0", false},
        // Tri-electrons
        {"HLT_e17_lhloose_2e9_lhloose", false},
        {"HLT_e17_lhloose_nod0_2e9_lhloose_nod0", false},
        {"HLT_e17_lhmedium_nod0_2e9_lhmedium_nod0", false},
        // single Muon
        {"HLT_mu20_iloose_L1MU15", false},
        {"HLT_mu24_ivarloose_L1MU15", false},
        {"HLT_mu24_ivarmedium", false},
        {"HLT_mu24_imedium", false},
        {"HLT_mu26_ivarmedium", false},
        {"HLT_mu26_imedium", false},
        // Di muon
        {"HLT_2mu10", false},
        {"HLT_mu18_mu8noL1", false},
        {"HLT_2mu10_nomucomb", false},
        {"HLT_mu20_mu8noL1", false},
        {"HLT_mu20_nomucomb_mu6noL1_nscan03", false},
        {"HLT_2mu14_nomucomb", false},
        {"HLT_mu22_mu8noL1", false},
        {"HLT_2mu14", false},
        // Trig Muons
        {"HLT_3mu6", false},
        {"HLT_3mu6_msonly", false},
        {"HLT_mu18_2mu4noL1", false},
        {"HLT_mu20_2mu4noL1", false},
        {"HLT_3mu4", false},
        {"HLT_mu6_2mu4", false},
        {"HLT_mu11_nomucomb_2mu4noL1_nscan03_L1MU11_2MU6", false},
        {"HLT_mu20_msonly_mu10noL1_msonly_nscan05_noComb", false},
        {"HLT_3mu6", false},
        {"HLT_3mu6_msonly", false},
        // elctron-muon
        {"HLT_e17_lhloose_mu14", false},
        {"HLT_2e12_lhloose_mu10", false},
        {"HLT_e12_lhloose_2mu10", false},
        {"HLT_e24_medium_L1EM20VHI_mu8noL1", false},
        {"HLT_e7_medium_mu24", false},
        {"HLT_e17_lhloose_nod0_mu14", false},
        {"HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1", false},
        {"HLT_e7_lhmedium_nod0_mu24", false},
        {"HLT_e12_lhloose_nod0_2mu10", false},
        {"HLT_2e12_lhloose_nod0_mu10", false}
        // 
    };
    CreateBranch();
}

bool UpsilonBranch::CreateBranch()
{
    return true;
}

UpsilonBranch::~UpsilonBranch(){
}

void UpsilonBranch::ClearBranch(){
    for(auto keys: trigger_map_) {
        keys.second = false;
    }
    pass_trigger_ = false;
    m_upsilon_ = 0;
    m_4l_ = 0;
    event_type_ = -1;
}

void UpsilonBranch::AttachBranchToTree(TTree& MyTree){

    // Trigger Info
    for(auto kv : trigger_map_){
        TString key(kv.first);
        TString br_name(kv.first);
        br_name.ReplaceAll("HLT", "trig");
        MyTree.Branch(br_name.Data(), &trigger_map_[key.Data()],
                Form("%s/O",br_name.Data()));
    }

    MyTree.Branch("passTrigger", &pass_trigger_, "passTrigger/O");
    MyTree.Branch("m_upsilon", &m_upsilon_, "m_upsilon/F");
    MyTree.Branch("m_4l", &m_4l_, "m_4l/F");
    MyTree.Branch("event_type", &event_type_, "event_type/I");
}

void UpsilonBranch::Fill()
{
    // Just fill variables in main code.
    return;
}
