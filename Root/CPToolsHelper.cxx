#include <stdlib.h>
#include <TTree.h>

#include "CPAnalysisExamples/errorcheck.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackParticlexAODHelpers.h"

#ifdef ROOTCORE 
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#endif

#include "MyXAODTools/CPToolsHelper.h"
#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

const char* CPToolsHelper::APP_NAME = "CPToolsHelper";

CPToolsHelper::CPToolsHelper(){
    initialize();
    iso_wp_ = "Loose";
}

CPToolsHelper::~CPToolsHelper(){
    delete grl_tool_; 
    if(iso_tool_) delete iso_tool_;
    delete ele_medium_LLH_tool_;
}

bool CPToolsHelper::initialize(){
    // GRL
    grl_tool_ = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    std::vector<std::string> myvals;
    string maindir(getenv("ROOTCOREBIN"));
    myvals.push_back(maindir+"/data/MyXAODTools/data15_13TeV_25ns_all2015_3_32fb.xml");
    myvals.push_back(maindir+"/data/MyXAODTools/data16_13TeV_all.xml");
    CHECK( grl_tool_->setProperty( "GoodRunsListVec", myvals) );
    CHECK( grl_tool_->setProperty("PassThrough", false) );
    CHECK( grl_tool_->initialize() );
    Info( APP_NAME, "GRL tool initialized... " );

    // Isolation tool
    iso_tool_ = new CP::IsolationSelectionTool("iso_tool");
    // CHECK(iso_tool_->setProperty("ElectronWP", "FixedCutLoose")); 
    // CHECK(iso_tool_->setProperty("MuonWP", "FixedCutLoose")); 
    CHECK(iso_tool_->setProperty("ElectronWP", "LooseTrackOnly")); 
    CHECK(iso_tool_->setProperty("MuonWP", "LooseTrackOnly")); 
    CHECK(iso_tool_->initialize());
    Info( APP_NAME, "Isolation tool initialized..." );

    // Medium Electron Likelihood
    ele_medium_LLH_tool_ = new AsgElectronLikelihoodTool("medium_ele_LLH"); 
    ele_medium_LLH_tool_->setProperty("primaryVertexContainer", "PrimaryVertices");
    ele_medium_LLH_tool_->setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150429/ElectronLikelihoodMediumOfflineConfig2015.conf");
    CHECK(ele_medium_LLH_tool_->initialize());

    // Trigger
    /**
    m_trigConfigTool_ = new TrigConf::xAODConfigTool("xAODConfigTool"); // gives us access to the meta-data
    CHECK( m_trigConfigTool_->initialize() );
    ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle( m_trigConfigTool_ );
    m_trigDecisionTool_ = new Trig::TrigDecisionTool("TrigDecisionTool");
    m_trigDecisionTool_->setProperty( "ConfigTool", trigConfigHandle );
    m_trigDecisionTool_->setProperty( "TrigDecisionKey", "xTrigDecision" );
    CHECK( m_trigDecisionTool_->initialize() );
    Info( APP_NAME, "Trigger decision tool initialized..." );
    **/

    return true;
}

bool CPToolsHelper::PassGRL(const xAOD::EventInfo& ei){
    if(!(ei.eventType(xAOD::EventInfo::IS_SIMULATION)) && 
            !(grl_tool_->passRunLB(ei.runNumber(), ei.lumiBlock()))
      ){
        return false;
    } else {
        return true;
    }
}

bool CPToolsHelper::PassGRL(int run_number, int lumi_block){
    return grl_tool_->passRunLB(run_number, lumi_block);
}

bool CPToolsHelper::PassEventCleaning(const xAOD::EventInfo& ei)
{
    // Instructions
    // https://twiki.cern.ch/twiki/bin/viewauth/Atlas/DataPreparationCheckListForPhysicsAnalysis
    bool bad_tag = (ei.errorState(xAOD::EventInfo::LAr) == xAOD::EventInfo::Error ||
            ei.errorState(xAOD::EventInfo::Tile) == xAOD::EventInfo::Error ||
            ei.errorState(xAOD::EventInfo::SCT) == xAOD::EventInfo::Error ||
            // ei.isEventFlagBitSet(xAOD::EventInfo::Core, 18) ||
            (ei.eventFlags(xAOD::EventInfo::Core) & 0x40000)  // incomplete events
            );
    return !bad_tag;
}

bool CPToolsHelper::HasPrimaryVertex(const xAOD::VertexContainer& vertice, unsigned int n_trks)
{
    bool passVertex = false;
    for(const auto& vxp : vertice){  
        if ( vxp->nTrackParticles() >= n_trks ){
            passVertex = true;
            break;
        }
    }
    return passVertex;
}


bool CPToolsHelper::PassIsolation(const xAOD::Muon& muon){
    return iso_tool_->accept(muon);
}

bool CPToolsHelper::PassIsolation(const xAOD::Electron& electron){
    return iso_tool_->accept(electron);
}

bool CPToolsHelper::PassEleMediumLLH(const xAOD::Electron& electron){
    return ele_medium_LLH_tool_->accept(electron);
}

bool CPToolsHelper::GetTrackSumPt(const xAOD::Vertex& vertex, 
        float& sum_px, float& sum_py)
{
    int n_tracks = (int) vertex.nTrackParticles();
    if(n_tracks <= 0) return false;

    float px = 0, py = 0;
    for(int i = 0; i < n_tracks; i++){
        const xAOD::TrackParticle* track = vertex.trackParticle(i);
        if(!track) continue;
        // add possible track quality cuts ?
        // float weight = vertex.trackWeight(i);
        px += track->p4().Px();
        py += track->p4().Py();
    }
    sum_px = px;
    sum_py = py;
    return true;
}

const xAOD::Vertex* CPToolsHelper::GetPrimVtx(const xAOD::VertexContainer& vertices)
{
    for ( const auto& vx : vertices ) {
      if (vx->vertexType() == xAOD::VxType::PriVtx) {
        return vx;
      }
    }
    return nullptr;
}

void CPToolsHelper::GetTrackQuality(const xAOD::TrackParticle* track, 
        const xAOD::EventInfo& evtInfo, 
        const xAOD::VertexContainer& vertices,
        float& d0, float& z0, float& zp)
{
  // d0
  d0 = xAOD::TrackingHelpers::d0significance( track , 
          evtInfo.beamPosSigmaX(), 
          evtInfo.beamPosSigmaY(), 
          evtInfo.beamPosSigmaXY() );
  // z0 
  const xAOD::Vertex* pv = GetPrimVtx(vertices);
  if(pv != nullptr){
      double primvertex_z = pv ? pv->z() : 0;
      z0 = (track->z0() + track->vz() - primvertex_z);
      zp = z0*TMath::Sin(track->theta());
  }
}

void CPToolsHelper::GetTrackQuality(const xAOD::Electron& input, 
        const xAOD::EventInfo& evtInfo, 
        const xAOD::VertexContainer& vertices,
        float& d0, float& z0, float& zp)
{
  const xAOD::TrackParticle* track =  input.trackParticle();
  GetTrackQuality(track, evtInfo, vertices, d0, z0, zp);
}

void CPToolsHelper::GetTrackQuality(const xAOD::Muon& input, 
        const xAOD::EventInfo& evtInfo, 
        const xAOD::VertexContainer& vertices,
        float& d0, float& z0, float& zp)
{
  const xAOD::TrackParticle* track =  input.primaryTrackParticle();
  GetTrackQuality(track, evtInfo, vertices, d0, z0, zp);
}

bool CPToolsHelper::GetProcessEventsInfo(xAOD::TEvent& event, 
        uint64_t& n_events_processed,
        double& sum_of_weights,
        double& sum_of_weights_squared)
{
    //Read the CutBookkeeper container
    const xAOD::CutBookkeeperContainer* completeCBC = 0;
    if (!event.retrieveMetaInput(completeCBC, "CutBookkeepers").isSuccess()) 
    {
        Error( APP_NAME, "Failed to retrieve CutBookkeepers from MetaData! Exiting.");
        return false;
    } else {
        Info( APP_NAME, "retrieve CutBookkeepers from MetaData!");
    }

    const xAOD::CutBookkeeper* allEventsCBK = 0;
    int maxCycle = -1;
    for ( auto cbk :  *completeCBC ) {
        if(!cbk) continue;
        if ( cbk->name() == "AllExecutedEvents" &&
             cbk->inputStream() == "StreamAOD" && 
             cbk->cycle() > maxCycle 
           )
        {
            allEventsCBK = cbk;
            maxCycle = cbk->cycle();
        }
    }
    if(allEventsCBK) {
        n_events_processed = allEventsCBK->nAcceptedEvents();
        sum_of_weights = allEventsCBK->sumOfEventWeights();
        sum_of_weights_squared = allEventsCBK->sumOfEventWeightsSquared();
        // Info(APP_NAME, "input stream: %s", allEventsCBK->inputStream().c_str());
    } else { 
        Info( APP_NAME, "No relevent CutBookKeepers found" ); 
        return false; 
    }	

    return true;
}

bool CPToolsHelper::GetProcessEventsInfo(const char* file_name, 
        uint64_t& n_events_processed,
        double& sum_of_weights,
        double& sum_of_weights_squared)
{
    TFile* ifile = TFile::Open(file_name, "READ");
    xAOD::TEvent event( xAOD::TEvent::kClassAccess );
    CHECK( event.readFrom( ifile ) );
    CHECK( GetProcessEventsInfo(event, n_events_processed, sum_of_weights, sum_of_weights_squared));
    ifile->Close();
    return true;
}

bool CPToolsHelper::SaveProcessedEvents(
        TTree& tree, const xAOD::EventInfo& ei,
        uint64_t total_evts_pro, double sum_of_evt_w,
        double sum_of_evt_w_sq)
{
    tree.Branch("nEventsProcessed", &total_evts_pro, "nEventsProcessed/l");
    tree.Branch("nSumEventWeights", &sum_of_evt_w, "nSumEventWeights/D");
    tree.Branch("nSumEventWeightsSquared", &sum_of_evt_w_sq, 
            "nSumEventWeightsSquared/D");
    int run_number_ = -1;
    int event_number_ = -1;
    int mc_channel_number_ = -1;
    tree.Branch("RunNumber", &run_number_, "RunNumber/I");
    tree.Branch("EventNumber", &event_number_, "EventNumber/I");
    tree.Branch("mc_channel_number", &mc_channel_number_, "mc_channel_number/I");

    run_number_ = ei.runNumber();
    event_number_ = ei.eventNumber();
    bool is_data = true;
    if(ei.eventType(xAOD::EventInfo::IS_SIMULATION)) 
    {
        mc_channel_number_= ei.mcChannelNumber();
        is_data = false;
    }
    tree.Fill();
    return is_data;
}

ST::SUSYObjDef_xAOD* CPToolsHelper::GetSUSYTools(bool isData, const char* config_name)
{
    Info( APP_NAME, "Creating SUSY Tools initialized... " );
    // create SUSYTools and config it
    ST::SUSYObjDef_xAOD* objTool = new ST::SUSYObjDef_xAOD("SUSYObjDef_Upsilon");
    objTool->msg().setLevel(MSG::ERROR);   // MSG::VERBOSE

    // Configure the SUSYObjDef instance
    ST::ISUSYObjDef_xAODTool::DataSource data_source = isData ? ST::ISUSYObjDef_xAODTool::Data : ST::ISUSYObjDef_xAODTool::FullSim;
    if(! objTool->setProperty("DataSource", data_source) ) return NULL;

    // general configuration
    string maindir(getenv("ROOTCOREBIN"));
    // string config_file = Form("%s/data/MyXAODTools/upsilon.conf", maindir.c_str());
    if(! objTool->setProperty("ConfigFile", config_name) ) return NULL;

    // pileup reweight
    vector<string> prw_conf;
    prw_conf.push_back(maindir+"/data/MyXAODTools/mc15c.prw.root");
    if(! objTool->setProperty("PRWConfigFiles", prw_conf) ) return NULL;

    vector<string> prw_lumicalc;
    prw_lumicalc.push_back(maindir+"/data/MyXAODTools/ilumicalc_histograms_None_276262-284484_final_20.7.root");
    prw_lumicalc.push_back(maindir+"/data/MyXAODTools/ilumicalc_histograms_None_297730-303892.root");
    if(! objTool->setProperty("PRWLumiCalcFiles", prw_lumicalc) ) return NULL;

    if( objTool->initialize() != StatusCode::SUCCESS){
        Error( APP_NAME, "Cannot intialize SUSYObjDef_xAOD..." );
        Error( APP_NAME, "Exiting... " );
        return NULL;
    }else{
        Info( APP_NAME, "SUSYObjDef_xAOD initialized... " );
    }
    return objTool;
}

bool CPToolsHelper::PassTrigger(const string& trig_name)
{
    return m_trigDecisionTool_->isPassed(trig_name); 
}
