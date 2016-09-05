#ifndef __MYXAODTOOLS_CPTOOLSHELPER_H__
#define __MYXAODTOOLS_CPTOOLSHELPER_H__

#include "xAODJet/Jet.h"
#include "xAODEgamma/Electron.h"
#include "xAODMuon/Muon.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/Vertex.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODBase/IParticle.h"

#include "CPAnalysisExamples/errorcheck.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"

#include "GoodRunsLists/GoodRunsListSelectionTool.h"
// #include "JetMomentTools/JetVertexTaggerTool.h"
#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "IsolationSelection/IsolationSelectionTool.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "PileupReweighting/PileupReweightingTool.h"
#include "AsgTools/ToolHandle.h"

#include <string>

using namespace std;
class CPToolsHelper{
public:
    CPToolsHelper();
    virtual ~CPToolsHelper();

    bool PassGRL(const xAOD::EventInfo& ei);
    bool PassGRL(int run_number, int lumi_block);
    float NewJVT(const xAOD::Jet& jet);

    bool PassIsolation(const xAOD::Muon& muon);
    bool PassIsolation(const xAOD::Electron& electron);
    bool PassEleMediumLLH(const xAOD::Electron& electron);

    static bool PassEventCleaning(const xAOD::EventInfo& ei);
    static bool HasPrimaryVertex(const xAOD::VertexContainer& vertice, 
            unsigned int n_trks = 2);

    // Track quality
    static bool GetTrackSumPt(const xAOD::Vertex& vertex,
            float& sum_px, float& sum_py);
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

    static bool GetProcessEventsInfo(xAOD::TEvent& event,
            uint64_t& n_events_processed,
            double& sum_of_weights,
            double& sum_of_weights_squared);

    static bool SaveProcessedEvents(
            TTree& tree, const xAOD::EventInfo& ei,
            uint64_t total_evts_pro, double sum_of_evt_w,
            double sum_of_evt_w_sq);

    static ST::SUSYObjDef_xAOD* GetSUSYTools(bool is_data, const char* config_name);

    bool PassTrigger(const string& trig_name);
private:
    static const char* APP_NAME;
    string iso_wp_ ;

    GoodRunsListSelectionTool *grl_tool_;
    CP::IsolationSelectionTool* iso_tool_;
    AsgElectronLikelihoodTool* ele_medium_LLH_tool_;
    // ToolHandle<CP::IPileupReweightingTool> m_prw_tool;
    TrigConf::xAODConfigTool *m_trigConfigTool_;
    Trig::TrigDecisionTool *m_trigDecisionTool_;


    bool initialize();
};
#endif
