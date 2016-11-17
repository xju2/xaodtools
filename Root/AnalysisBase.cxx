#include "MyXAODTools/AnalysisBase.h"

AnalysisBase::AnalysisBase(
        const char* file_name,
        const char* tree_name,
        const char* associate_tree_name)
{
    m_debug = false;
    m_withGRL = true;
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


}

int AnalysisBase::initializeBasicTools()
{
    // initiate tools
    cp_tools =  unique_ptr<CPToolsHelper>( new CPToolsHelper() );
    if(! GetSUSYTool()){
        return 1;
    }
    return 0;
}

void AnalysisBase::setSUSYConfig(const string& config){
    m_susy_config = config;
}

AnalysisBase::~AnalysisBase()
{
    if (f_out) {
        f_out->cd();
        tree->Write();
        physics->Write();
        f_out->Close();
    }
}

void AnalysisBase::CreateBasicBranch()
{
    event_br =  unique_ptr<EventInfoCreator>( new EventInfoCreator() );
    muon_br =   unique_ptr<MuonBranch>( new MuonBranch() );
    el_br =     unique_ptr<ElectronBranch>( new ElectronBranch() );
    jet_br =    unique_ptr<JetBranch>( new JetBranch() );
    ph_br =     unique_ptr<PhotonBranch>( new PhotonBranch() );
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
    // m_totalEvents = total_evts;
    return true;
}

bool AnalysisBase::GetSUSYTool(const char* config)
{
    if(m_susy_config != ""){
        m_objTool.reset( CPToolsHelper::GetSUSYTools(m_isData, m_susy_config.c_str()) );
    } else if(!config){
        m_objTool.reset( CPToolsHelper::GetSUSYTools(m_isData, config) );
    } else {
        Error(APP_NAME, "failed to get SUSYTool");
        return false;
    }
    return true;
}

void AnalysisBase::AttachBasicToTree()
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

void AnalysisBase::ClearBasicBranch()
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

int AnalysisBase::Start(Long64_t ientry)
{
    event->getEntry( ientry );

    // check if isData flag is set true when running data, or program will complain!
    CHECK(m_objTool->ApplyPRWTool());

    CHECK( event->retrieve( ei, "EventInfo" ) );
    //Hack...
    // if(ei->eventNumber() != 1773037184) return 1;

    if(m_debug) Info(APP_NAME, " AnalysisBase:processing: %d %llu", (int) ei->runNumber(), ei->eventNumber());
    if( ientry%1000 == 0) {
        Info(APP_NAME, "processed %lld, %.2f %%", ientry, 100.*ientry/m_totalEvents);
    }

    for(auto& kv : trigger_map_)
    {
        if(m_objTool->IsTrigPassed(kv.first.c_str()))
        {
            if(m_debug) Info(APP_NAME, "fired the trigger %s", kv.first.c_str());
            pass_trigger_ = kv.second = true;
        } else {
            kv.second = false;
        }
    }

    CHECK( event->retrieve(vertice, "PrimaryVertices") );

    if(! cp_tools->HasPrimaryVertex(*vertice, 1))
        return 1;

    pv = CPToolsHelper::GetPrimVtx( *vertice );

    if(m_withGRL && ! cp_tools->PassGRL(*ei) ){
        return 1;
    }
    return 0;
}

void AnalysisBase::SetVerbose(){
    m_debug = true;
    if(m_objTool.get()) {
        m_objTool->setProperty("DebugMode", m_debug).ignore();
    }
}
void AnalysisBase::setGRLTag(bool tag){
    m_withGRL = tag;
}
void AnalysisBase::SetTotalEventsToProcess(Long64_t nentries)
{
    m_totalEvents = nentries;
}
