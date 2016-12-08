#!/usr/bin/evn python

head_file="""
/*
 * description of the analysis
 *
 */
#ifndef __MYXAODTOOLS_ANALYSISNAME_H__
#define __MYXAODTOOLS_ANALYSISNAME_H__

#include <vector>
#include <string>

#include "MyXAODTools/AnalysisBase.h"

#include "AsgTools/ToolHandle.h"

using namespace std;

class ANALYSISNAME : public AnalysisBase
{
public:
    ANALYSISNAME();
    virtual ~ANALYSISNAME();

    int initialize();
    void ClearBranch();
    int process(Long64_t ientry); // main program

private:
    /** private methods */
    void AttachBranchToTree();
    void CreateBranch();

private:
    /* specific branches used in this analysis */

private:
    /* specific Tools used in this analysis */
};

#endif

"""

src_file = """
#include <stdlib.h>

#include <TFile.h>

#include "MyXAODTools/ANALYSISNAME.h"
#include "MyXAODTools/Helper.h"
#include "CPAnalysisExamples/errorcheck.h"
#include "xAODBase/IParticleHelpers.h"

#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODMissingET/MissingETContainer.h"

ANALYSISNAME::ANALYSISNAME():
    AnalysisBase()
{
    if(APP_NAME==NULL) APP_NAME = "ANALYSISNAME";
    string maindir(getenv("ROOTCOREBIN"));

    // don't forget to change your SUSY configuration!
    m_susy_config = Form("%s/data/MyXAODTools/monojet.conf", maindir.c_str());

    trigger_map_ = { // your triggers go here.
    // {"HLT_xe70", false},
    };
}

int ANALYSISNAME::initialize()
{
    if( initializeBasicTools() != 0 ){
        return 1;
    }
    CreateBranch();
    AttachBranchToTree();

    // initiate tools
    return 0;
}

ANALYSISNAME::~ANALYSISNAME(){
}

void ANALYSISNAME::CreateBranch()
{
    CreateBasicBranch();
    return ;
}

void ANALYSISNAME::ClearBranch(){
    ClearBasicBranch();

}


void ANALYSISNAME::AttachBranchToTree()
{
    AttachBasicToTree();

    // event_br->AttachBranchToTree(*physics);
    // muon_br, el_br, jet_br, ph_br

    // set your own branches
    // physics->Branch("has_bad_muon", &m_hasBadMuon, "has_bad_muon/O");
}

int ANALYSISNAME::process(Long64_t ientry)
{
    int sc = Start(ientry);
    if(m_debug) {
        Info(APP_NAME, " ANALYSISNAME: processing");
    }
    if(sc != 0) return sc;
    // event_br->Fill(*ei);

    // Start ANALYSISNAME 

    /*get physics objects*/
    // Electrons
    // xAOD::ElectronContainer* electrons_copy = NULL;
    // xAOD::ShallowAuxContainer* electrons_copyaux = NULL;
    // CHECK( m_objTool->GetElectrons(electrons_copy, electrons_copyaux, true) );

    // Muons
    // xAOD::MuonContainer* muons_copy = NULL;
    // xAOD::ShallowAuxContainer* muons_copyaux = NULL;
    // CHECK( m_objTool->GetMuons(muons_copy, muons_copyaux, true) );

    // Jets
    // xAOD::JetContainer* jets_copy = NULL;
    // xAOD::ShallowAuxContainer* jets_copyaux = NULL;
    // CHECK( m_objTool->GetJets(jets_copy,jets_copyaux, true) );

    // Photons
    // xAOD::PhotonContainer* ph_copy = nullptr;
    // xAOD::ShallowAuxContainer* phAux_copy = nullptr;
    // CHECK(m_objTool->GetPhotons(ph_copy, phAux_copy, true));

    ///////////////////////
    // do overlap removal before object selection
    // turn off the harmonization
    ///////////////////////
    // bool doHarmonization = false;
    // CHECK( m_objTool->OverlapRemoval(
    //             electrons_copy, muons_copy,
    //             jets_copy, ph_copy) );

    //discard the event if any jets is labelled as 'bad'
    // bool passJetCleaning = true;
    // for(const auto& jet : *jets_copy){
    //     m_objTool->IsBJet(*jet) ;
    //     if ( jet->pt() > 20e3 )
    //     {
    //         if( dec_bad(*jet) && dec_passOR(*jet)){
    //             passJetCleaning = false;
    //             break;
    //         }
    //     }
    // }
    // if ( !passJetCleaning ) return 1;

    // electron selections
    // for(const auto&  el : *electrons_copy){
    //     if( !(bool) dec_baseline(*el) || !(bool) dec_passOR(*el)){
    //         continue;
    //     }
    //     m_nBaseEl ++;
    // }
    // // muon selections
    // for(const auto& mu : *muons_copy){
    //     if( !(bool) dec_baseline(*mu) || !(bool) dec_passOR(*mu) ){
    //         continue;
    //     }
    //     if( dec_bad(*mu) ) m_hasBadMuon = true;
    //     if( dec_cosmic(*mu) ) m_hasCosmicMuon = true;
    //     muon_br->Fill(*mu, ei, vertice);
    //     m_nBaseMu ++;
    // }
    // // photon selections
    // for(const auto& ph : *ph_copy) {
    //     if( !(bool) dec_baseline(*ph) || !(bool) dec_passOR(*ph) ){
    //         continue;
    //     }
    //     m_nBasePhoton ++;
    // }

    // Fill your tree!!!!
    physics->Fill();
    return 0;
}

"""
