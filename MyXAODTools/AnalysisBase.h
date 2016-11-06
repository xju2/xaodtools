#ifndef __MYXAODTOOLS_ANALYSISBASE_H__
#define __MYXAODTOOLS_ANALYSISBASE_H__

// Infrastructure include(s):
#ifdef ROOTCORE
#   include "xAODRootAccess/Init.h"
#   include "xAODRootAccess/TEvent.h"
#endif // ROOTCORE
#include <memory>

#include <TFile.h>
#include <TTree.h>

#include "CPAnalysisExamples/errorcheck.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"

#include "MyXAODTools/CPToolsHelper.h"
#include "MyXAODTools/EventInfoCreator.h"
#include "MyXAODTools/ElectronBranch.h"
#include "MyXAODTools/MuonBranch.h"
#include "MyXAODTools/JetBranch.h"
#include "MyXAODTools/PhotonBranch.h"

static SG::AuxElement::Decorator<char> dec_baseline("baseline");
static SG::AuxElement::Decorator<char> dec_signal("signal");
static SG::AuxElement::Decorator<char> dec_passOR("passOR");
static SG::AuxElement::Decorator<char> dec_bad("bad");
static SG::AuxElement::Decorator<char> dec_bjet("bjet");
static SG::AuxElement::Decorator<char> dec_cosmic("cosmic");
static SG::AuxElement::Decorator<char> dec_tightBad("tightBad");

using namespace std;
class AnalysisBase{

public:
    AnalysisBase(
            const char* file_name="reduced_ntup.root",
            const char* tree_name="physics",
            const char* associate_tree_name="associate");
    virtual ~AnalysisBase();

    virtual int initialize() = 0;
    virtual void ClearBranch() = 0;
    virtual int process(Long64_t ientry) = 0;

    void setSUSYConfig(const string& config);
    void SetEvent(xAOD::TEvent* event);
    void SetVerbose();

    bool SaveProcessedInfo(uint64_t total_evts, double sum_of_weight, double sum_of_w_sq);
    static bool descend_on_pt(xAOD::IParticle* p1, xAOD::IParticle* p2){
        return p1->pt() > p2->pt();
    }


    bool isPassGRL();
    void setGRLTag(bool);
protected: // methods

    bool GetSUSYTool(const char* config=NULL);
    // bookings of tree
    int initializeBasicTools();
    void AttachBasicToTree();
    void CreateBasicBranch();
    void ClearBasicBranch();
    int Start(Long64_t ientry);

protected:
    bool m_isData;
    bool m_debug;
    bool m_withGRL;
    string  m_susy_config;
    const char* APP_NAME;
    map<string, bool> trigger_map_;
    uint64_t m_totalEvents;

    // points to physics objects
    xAOD::TEvent* event;
    const xAOD::EventInfo* ei;
    const xAOD::VertexContainer* vertice;
    const xAOD::Vertex* pv;

    // tools for branch booking
    unique_ptr<EventInfoCreator> event_br;
    unique_ptr<MuonBranch> muon_br;
    unique_ptr<ElectronBranch> el_br;
    unique_ptr<JetBranch> jet_br;
    unique_ptr<PhotonBranch> ph_br;

    // tools for CP recommendations
    unique_ptr<CPToolsHelper> cp_tools;
    ST::SUSYObjDef_xAOD* m_objTool;

    // output setup
    TFile* f_out;
    TTree* tree;
    TTree* physics;

    // output branches
    bool pass_trigger_;

};
#endif
