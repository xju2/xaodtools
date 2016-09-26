#ifndef __MYXAODTOOLS_ANALYSISBASE_H__
#define __MYXAODTOOLS_ANALYSISBASE_H__

// Infrastructure include(s):
#ifdef ROOTCORE
#   include "xAODRootAccess/Init.h"
#   include "xAODRootAccess/TEvent.h"
#endif // ROOTCORE

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

class AnalysisBase{

public:
    AnalysisBase(
            const char* file_name="reduced_ntup.root",
            const char* tree_name="physics",
            const char* associate_tree_name="associate");
    virtual ~AnalysisBase();

    void SetEvent(xAOD::TEvent* event);
    void SetVerbose();

    bool SaveProcessedInfo(uint64_t total_evts, double sum_of_weight, double sum_of_w_sq);
    void GetSUSYTool(const char* config=NULL);
    static bool descend_on_pt(xAOD::IParticle* p1, xAOD::IParticle* p2){
        return p1->pt() > p2->pt();
    }

    // bookings of tree
    virtual void AttachBranchToTree();
    virtual void CreateBranch();
    virtual void ClearBranch();

    virtual int process(Long64_t ientry);
    bool isPassGRL();

protected:
    bool m_isData;
    bool m_debug;
    string  m_susy_config;
    const char* APP_NAME;
    map<string, bool> trigger_map_;

    // points to physics objects
    xAOD::TEvent* event;
    const xAOD::EventInfo* ei;
    const xAOD::VertexContainer* vertice;
    const xAOD::Vertex* pv;

    // tools for branch booking
    EventInfoCreator* event_br;
    MuonBranch* muon_br;
    ElectronBranch* el_br;
    JetBranch* jet_br;
    PhotonBranch* ph_br;

    // tools for CP recommendations
    CPToolsHelper* cp_tools;
    ST::SUSYObjDef_xAOD* m_objTool;

    // output setup
    TFile* f_out;
    TTree* tree;
    TTree* physics;

    // output branches
    bool pass_trigger_;

};
#endif
