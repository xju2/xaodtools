#include "MyXAODTools/AnalysisBase.h"

AnalysisBase::AnalysisBase(
        const char* file_name,
        const char* tree_name,
        const char* associate_tree_name)
{
    m_debug = false;
    APP_NAME = NULL;
    trigger_map_.clear();
    event = NULL;
    ei = NULL;
    vertice = NULL;
    pv = NULL;
    m_susy_config = "";

    f_out = TFile::Open(file_name, "recreate");
    tree = new TTree(associate_tree_name, associate_tree_name);
    physics = new TTree(tree_name, tree_name);

    // initiate tools
    cp_tools = NULL;
    event_br = new EventInfoCreator();
    muon_br = new MuonBranch();
    el_br = new ElectronBranch();
    jet_br = new JetBranch();
    ph_br = new PhotonBranch();
    m_objTool = NULL;

    CreateBranch();
    AttachBranchToTree();
}

AnalysisBase::~AnalysisBase()
{
    if (f_out) {
        f_out->cd();
        tree->Write();
        physics->Write();
        f_out->Close();
    }
    if (cp_tools) delete cp_tools;
    if (event_br) delete event_br;
    if (muon_br) delete muon_br;
    if (el_br) delete el_br;
    if (jet_br) delete jet_br;
    if (ph_br) delete ph_br;
}

void AnalysisBase::CreateBranch()
{
    return ;
}

void AnalysisBase::SetEvent(xAOD::TEvent* event){
    this->event = event;
}


bool AnalysisBase::SaveProcessedInfo(
        uint64_t total_evts, double sum_of_weight, double sum_of_w_sq)
{
    event->getEntry(0);
    CHECK( event->retrieve( ei, "EventInfo" ) );
    if(ei) {
        m_isData = CPToolsHelper::SaveProcessedEvents(
                *tree, *ei, total_evts, sum_of_weight, sum_of_w_sq);
    }
    return true;
}

void AnalysisBase::GetSUSYTool(const char* config)
{
    if(m_susy_config != ""){
        m_objTool = CPToolsHelper::GetSUSYTools(m_isData, m_susy_config.c_str());
    } else if(!config){
        m_objTool = CPToolsHelper::GetSUSYTools(m_isData, config);
    } else {
    }
    return;
}

void AnalysisBase::AttachBranchToTree()
{
    // Trigger Info
    for(auto& kv : trigger_map_){
        TString key(kv.first);
        TString br_name(kv.first);
        br_name.ReplaceAll("HLT", "trig");
        physics->Branch(br_name.Data(), &trigger_map_[key.Data()],
                Form("%s/O",br_name.Data()));
    }

    if(physics->GetBranch("passTrigger") == NULL) {
        physics->Branch("passTrigger", &pass_trigger_, "passTrigger/O");
    }
}

void AnalysisBase::ClearBranch()
{
    for(auto keys: trigger_map_) {
        keys.second = false;
    }
    pass_trigger_ = false;

    if(event_br)    event_br->ClearBranch();
    if(muon_br)     muon_br->ClearBranch();
    if(el_br)       el_br->ClearBranch();
    if(jet_br)      jet_br->ClearBranch();
    if(ph_br)       ph_br   ->ClearBranch();
}

int AnalysisBase::process(Long64_t ientry)
{
    event->getEntry( ientry );
    CHECK( event->retrieve( ei, "EventInfo" ) );
    //Hack...
    // if(ei->eventNumber() != 1773037184) return 1;

    if(m_debug) Info(APP_NAME, " AnalysisBase:processing: %d %llu", (int) ei->runNumber(), ei->eventNumber());

    CHECK(m_objTool->ApplyPRWTool());


    for(auto& kv : trigger_map_)
    {
        if(m_objTool->IsTrigPassed(kv.first.c_str())) 
        {
            if(m_debug) Info(APP_NAME, "fired the trigger %s", kv.first.c_str());
            pass_trigger_ = kv.second = true;
        } else {
            kv.second = false;
        }
        if(m_debug){
            Info(APP_NAME, "%s trigger: %d %d", kv.first.c_str(),
                    (int) kv.second, (int) pass_trigger_);
        }
    }

    CHECK( event->retrieve(vertice, "PrimaryVertices") );

    if(cp_tools == NULL){
        cp_tools = new CPToolsHelper();
    }

    if(! cp_tools->HasPrimaryVertex(*vertice, 1))
        return 1;

    pv = CPToolsHelper::GetPrimVtx( *vertice );

    if(! cp_tools->PassGRL(*ei) ){
        return 1;
    }
    return 0;
}

void AnalysisBase::SetVerbose(){
    m_debug = true;
    if(m_objTool) {
        m_objTool->setProperty("DebugMode", m_debug).ignore();
    }
}
