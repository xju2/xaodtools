#ifndef __MYXAODTOOLS_CPTOOLSHELPER_H__
#define __MYXAODTOOLS_CPTOOLSHELPER_H__

#include "xAODJet/Jet.h"
#include "xAODEgamma/Electron.h"
#include "xAODMuon/Muon.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/Vertex.h"
#include "xAODTracking/TrackParticle.h"

#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "JetMomentTools/JetVertexTaggerTool.h"
#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "IsolationSelection/IsolationSelectionTool.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"

#include <string>

using namespace std;
class CPToolsHelper{
public:
    CPToolsHelper();
    virtual ~CPToolsHelper();

    bool PassGRL(const xAOD::EventInfo& ei);
    float NewJVT(const xAOD::Jet& jet);
    bool PassTrigger(const string& trig_name);
    bool PassIsolation(const xAOD::Muon& muon);
    bool PassIsolation(const xAOD::Electron& electron);
    bool PassEleMediumLLH(const xAOD::Electron& electron);

    static bool PassEventCleaning(const xAOD::EventInfo& ei);
    static bool HasPrimaryVertex(const xAOD::VertexContainer& vertice, 
            unsigned int n_trks = 2);
    static bool GetTrackSumPt(const xAOD::Vertex& vertex, 
            float& sum_px, float& sum_py);
    // Track quality
    static void GetTrackQuality(const xAOD::Electron& el, 
            const xAOD::EventInfo& ei, 
            const xAOD::VertexContainer& vertices,
            float& d0, float& z0, float& zp);
    static void GetTrackQuality(const xAOD::Muon& muon,
            const xAOD::EventInfo& ei, 
            const xAOD::VertexContainer& vertices,
            float& d0, float& z0, float& zp);
    static void GetTrackQuality(const xAOD::TrackParticle* track, 
            const xAOD::EventInfo& ei, 
            const xAOD::VertexContainer& vertices,
            float& d0, float& z0, float& zp);
    static const xAOD::Vertex* GetPrimVtx(const xAOD::VertexContainer& vertices);
    static bool GetProcessEventsInfo(const char* file_name, 
            uint64_t& n_events_processed,
            double& sum_of_weights,
            double& sum_of_weights_squared);

private:    
    static const char* APP_NAME; 
    string iso_wp_ ;

    GoodRunsListSelectionTool *grl_tool_;
    JetVertexTaggerTool* pjvtag_tool_; 
    TrigConf::xAODConfigTool *m_trigConfigTool_;
    Trig::TrigDecisionTool *m_trigDecisionTool_;
    CP::IsolationSelectionTool* iso_tool_;
    AsgElectronLikelihoodTool* ele_medium_LLH_tool_;


    bool initialize();
};
#endif
