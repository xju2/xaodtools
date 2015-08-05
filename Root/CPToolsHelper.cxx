#include <stdlib.h>

#include "CPAnalysisExamples/errorcheck.h"
#include "xAODTracking/TrackParticle.h"

#include "MyXAODTools/CPToolsHelper.h"

CPToolsHelper::CPToolsHelper(){
    initialize();
    iso_wp_ = "Loose";
}

CPToolsHelper::~CPToolsHelper(){
    delete grl_tool_; 
    delete pjvtag_tool_;
    delete m_trigConfigTool_;
    delete m_trigDecisionTool_;
    delete iso_tool_;
    delete ele_medium_LLH_tool_;
}

bool CPToolsHelper::initialize(){
    // GRL
    grl_tool_ = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    std::vector<std::string> myvals;
    string maindir(getenv("ROOTCOREBIN"));
    myvals.push_back(maindir+"/data/MyXAODTools/data15_13TeV.CombinedPerf.xml");
    CHECK( grl_tool_->setProperty( "GoodRunsListVec", myvals) );
    CHECK( grl_tool_->setProperty("PassThrough", false) );
    CHECK( grl_tool_->initialize() );
    Info( APP_NAME, "GRL tool initialized... " );

    // JVT
    pjvtag_tool_ = new JetVertexTaggerTool("jvtag");
    CHECK(pjvtag_tool_->setProperty("JVTFileName", "JetMomentTools/JVTlikelihood_20140805.root"));
    CHECK(pjvtag_tool_->initialize());
    Info( APP_NAME, "JVT tool initialized... " );

    // Trigger
    m_trigConfigTool_ = new TrigConf::xAODConfigTool("xAODConfigTool"); // gives us access to the meta-data
    CHECK( m_trigConfigTool_->initialize() );
    ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle( m_trigConfigTool_ );
    m_trigDecisionTool_ = new Trig::TrigDecisionTool("TrigDecisionTool");
    m_trigDecisionTool_->setProperty( "ConfigTool", trigConfigHandle );
    m_trigDecisionTool_->setProperty( "TrigDecisionKey", "xTrigDecision" );
    CHECK( m_trigDecisionTool_->initialize() );
    Info( APP_NAME, "Trigger decision tool initialized..." );

    // Isolation tool
    iso_tool_ = new CP::IsolationSelectionTool("iso_tool");
    CHECK(iso_tool_->setProperty("ElectronWP", "Loose")); 
    CHECK(iso_tool_->setProperty("MuonWP", "Loose")); 
    CHECK(iso_tool_->initialize());
    Info( APP_NAME, "Isolation tool initialized..." );

    // Medium Electron Likelihood
    ele_medium_LLH_tool_ = new AsgElectronLikelihoodTool("medium_ele_LLH"); 
    ele_medium_LLH_tool_->setProperty("primaryVertexContainer", "PrimaryVertices");
    ele_medium_LLH_tool_->setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150429/ElectronLikelihoodMediumOfflineConfig2015.conf");
    CHECK(ele_medium_LLH_tool_->initialize());
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

float CPToolsHelper::NewJVT(const xAOD::Jet& jet)
{
    return pjvtag_tool_->updateJvt(jet);
}

bool CPToolsHelper::PassEventCleaning(const xAOD::EventInfo& ei)
{
    bool bad_tag = (ei.errorState(xAOD::EventInfo::LAr) == xAOD::EventInfo::Error ||
            ei.errorState(xAOD::EventInfo::Tile) == xAOD::EventInfo::Error ||
            ei.isEventFlagBitSet(xAOD::EventInfo::Core, 18));
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

bool CPToolsHelper::PassTrigger(const string& trig_name)
{
    return m_trigDecisionTool_->isPassed(trig_name); 
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
