
#include <stdlib.h>

#include <TFile.h>

#include "MyXAODTools/FakeMuonAna.h"
#include "MyXAODTools/Helper.h"
#include "CPAnalysisExamples/errorcheck.h"
#include "xAODBase/IParticleHelpers.h"

#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "MuonSelectorTools/MuonSelectionTool.h"

FakeMuonAna::FakeMuonAna():
    AnalysisBase()
{
    if(APP_NAME==NULL) APP_NAME = "FakeMuonAna";
    string maindir(getenv("ROOTCOREBIN"));

    // don't forget to change your SUSY configuration!
    m_susy_config = Form("%s/data/MyXAODTools/fakemuon.conf", maindir.c_str());

    trigger_map_ = { // your triggers go here.
        {"HLT_3mu4_nomucomb_delayed", false},
        {"HLT_3mu4", false},
        {"HLT_3mu4_bDimu", false},
        {"HLT_3mu4_bUpsi", false},
        {"HLT_3mu4_bUpsi_delayed", false}
    };
}

int FakeMuonAna::initialize()
{
    if( initializeBasicTools() != 0 ){
        return 1;
    }
    CreateBranch();
    AttachBranchToTree();

    // initiate tools
    return 0;
}

FakeMuonAna::~FakeMuonAna(){
}

void FakeMuonAna::CreateBranch()
{
    CreateBasicBranch();
    m_isPromptMuon.reset( new vector<bool>() );
    m_isPassTrack.reset( new vector<bool>() );
    m_isPassPresel.reset( new vector<bool>() );
    m_isPassLoose.reset( new vector<bool>() );
    m_isPassMedium.reset( new vector<bool>() );

    return ;
}

void FakeMuonAna::ClearBranch(){
    ClearBasicBranch();

    m_isPromptMuon  ->clear();
    m_isPassTrack   ->clear();
    m_isPassPresel  ->clear();
    m_isPassLoose   ->clear();
    m_isPassMedium  ->clear();
}


void FakeMuonAna::AttachBranchToTree()
{
    AttachBasicToTree();

    event_br->AttachBranchToTree(*physics);
    muon_br ->AttachBranchToTree(*physics);

    // set your own branches
    physics->Branch("mu_isPrompt", m_isPromptMuon.get());
    physics->Branch("mu_isPassTrack", m_isPassTrack.get());
    physics->Branch("mu_isPassPresel", m_isPassPresel.get());
    physics->Branch("mu_isPassLoose", m_isPassLoose.get());
    physics->Branch("mu_isPassMedium", m_isPassMedium.get());
}

int FakeMuonAna::process(Long64_t ientry)
{
    int sc = Start(ientry);
    if(m_debug) {
        Info(APP_NAME, " FakeMuonAna: processing");
    }
    if(sc != 0) return sc;
    event_br->Fill(*ei);

    // Start FakeMuonAna 

    /*get physics objects*/

    // Muons
    xAOD::MuonContainer* muons_copy = NULL;
    xAOD::ShallowAuxContainer* muons_copyaux = NULL;
    CHECK( m_objTool->GetMuons(muons_copy, muons_copyaux, true) );


    // muon selections
    for(const auto& mu : *muons_copy){
        m_isPromptMuon->push_back( isPrompt(*mu) );
        
        m_isPassTrack->push_back( muon_br->m_muonSelectionTool->passedIDCuts(*mu) );
        m_isPassPresel->push_back( muon_br->m_muonSelectionTool->passedMuonCuts(*mu) );
        auto quality = muon_br->m_muonSelectionTool->getQuality(*mu);
        m_isPassLoose->push_back( (bool) (quality <= xAOD::Muon::Quality::Loose) );
        m_isPassMedium->push_back( (bool) (quality <= xAOD::Muon::Quality::Medium) );
        muon_br->Fill(*mu, ei, vertice);
    }

    // Fill your tree!!!!
    physics->Fill();
    return 0;
}

bool FakeMuonAna::isPrompt(const xAOD::Muon& muon) const
{
    return true;
}
