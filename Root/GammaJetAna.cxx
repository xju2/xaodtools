#include <stdlib.h>


#include "MyXAODTools/GammaJetAna.h"
#include "MyXAODTools/Helper.h"
#include "CPAnalysisExamples/errorcheck.h"

#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/ElectronAuxContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"


GammaJetAna::GammaJetAna():
    AnalysisBase(),
    LEADING_PHOTON_CUT(150E3),
    LEADING_JET_CUT(150E3)
{
    if(APP_NAME==NULL) APP_NAME = "GammaJetAna";
    string maindir(getenv("ROOTCOREBIN"));
    m_susy_config = Form("%s/data/MyXAODTools/gamma_jet.conf", maindir.c_str());

    trigger_map_ = {
        {"HLT_g120_loose", false},
        {"HLT_g200_etcut", false}
    };
    CreateBranch();
    AttachBranchToTree();
    h_cutflow = new TH1F("h_cutflow", "cut flow", 101, -0.5, 100.5);
}


GammaJetAna::~GammaJetAna(){
    if(f_out){
        f_out->cd();
        h_cutflow->Write();
        tree->Write();
        physics->Write();
        f_out->Close();
    }
}

void GammaJetAna::CreateBranch()
{
    return ;
}

void GammaJetAna::ClearBranch(){
    AnalysisBase::ClearBranch();
    m_mass = -999;
}


void GammaJetAna::AttachBranchToTree()
{
    AnalysisBase::AttachBranchToTree();

    event_br->AttachBranchToTree(*physics);
    jet_br  ->AttachBranchToTree(*physics);
    ph_br   ->AttachBranchToTree(*physics);
    el_br   ->AttachBranchToTree(*physics);
    physics->Branch("mass", &m_mass, "mass/F"); 
}

int GammaJetAna::process(Long64_t ientry)
{
    h_cutflow->Fill(0);
    int sc = AnalysisBase::process(ientry);
    if(m_debug) {
        Info(APP_NAME, " GammaJetAna: processing");
    }

    if(sc != 0) return sc;
    event_br->Fill(*ei);

    // Get Photons
    xAOD::PhotonContainer* ph_copy = NULL;
    xAOD::ShallowAuxContainer* ph_copyaxu = NULL;
    CHECK( m_objTool->GetPhotons(ph_copy, ph_copyaxu, true) );
    sort(ph_copy->begin(), ph_copy->end(), descend_on_pt);
    if(m_debug) Info(APP_NAME, "Got photons");

    // Get Electrons
    xAOD::ElectronContainer* el_copy = NULL;
    xAOD::ShallowAuxContainer* el_copyaux = NULL;
    CHECK( m_objTool->GetElectrons(el_copy, el_copyaux, true) );
    if(m_debug) Info(APP_NAME, "Got Electrons");
    // Fill baseline electrons
    for(auto el_itr = el_copy->begin(); el_itr != el_copy->end(); ++el_itr){
        if(! dec_baseline(**el_itr) ) continue;
        el_br->Fill( **el_itr );
    }

    // Get Jets
    xAOD::JetContainer* jets_copy = NULL;
    xAOD::ShallowAuxContainer* jets_copyaux = NULL;
    CHECK( m_objTool->GetJets(jets_copy, jets_copyaux, true) );
    sort(jets_copy->begin(), jets_copy->end(), descend_on_pt);
    if(m_debug) Info(APP_NAME, "Got Jets");

    // find leading good photon
    vector<xAOD::Photon*> good_photons;
    auto leading_ph_id = ph_copy->end();
    for(auto ph_itr = ph_copy->begin(); ph_itr != ph_copy->end(); ++ph_itr) {
        if(! dec_signal(**ph_itr) ) continue;
        xAOD::Photon* photon = (*ph_itr);
        if(photon->pt() > LEADING_PHOTON_CUT){
            leading_ph_id = ph_itr;
        }
        // Save all good photons
        ph_br->Fill(**ph_itr);
    }
    if(leading_ph_id == ph_copy->end()) return 2;
    if(m_debug) Info(APP_NAME, "Got leading photon");

    // find leading good jet
    auto leading_jet_id = jets_copy->end();
    for(auto jet_itr = jets_copy->begin(); jet_itr != jets_copy->end(); ++jet_itr){
        if(! dec_signal(**jet_itr)) continue;
        // Save all good jets 
        jet_br->Fill(**jet_itr);

        // check overlap
        bool overlap_w_ph = false;
        for(auto ph_itr = ph_copy->begin(); ph_itr != ph_copy->end(); ++ph_itr){
            if(! dec_baseline(**ph_itr) ) continue;
            if( (*ph_itr)->p4().DeltaR( (*jet_itr)->p4() ) < 0.4 ){
                overlap_w_ph = true;
                break;
            }
        }
        if(overlap_w_ph) continue;
        bool overlap_w_el = false;
        for(auto el_itr = el_copy->begin(); el_itr != el_copy->end(); ++el_itr){
            if(! dec_baseline(**el_itr) ) continue;
            if( (*el_itr)->p4().DeltaR( (*el_itr)->p4() ) < 0.2) {
               overlap_w_el = true;
               break;
            }
        }
        if(overlap_w_el) continue;

        if( (*jet_itr)->pt() > LEADING_JET_CUT){
            leading_jet_id = jet_itr;
        }
    }
    if(leading_jet_id == jets_copy->end()) return 3;

    //photon and jet separated
    float delta_eta = fabs((*leading_ph_id)->eta() - (*leading_jet_id)->eta());
    if(delta_eta <= 1.6) return 4;


    // leading photon is well isolation from other jets.
    bool has_overlap_jet = false;
    float ph_eta = (*leading_jet_id)->eta();
    float ph_phi = (*leading_jet_id)->phi();
    for(auto jet_itr = jets_copy->begin(); jet_itr != jets_copy->end(); ++jet_itr)
    {
        if(! dec_signal(**jet_itr) || ! dec_passOR(**jet_itr)) continue;
        if( (*jet_itr)->pt() < 30E3 ) continue;

        float jet_eta = (*jet_itr)->eta();
        float jet_phi = (*jet_itr)->phi();
        float delta_R = MyXAODTools::delta_r(jet_eta, jet_phi, ph_eta, ph_phi);
        if(delta_R < 0.8){
            has_overlap_jet = true;
            break;
        }
    }
    if(has_overlap_jet) return 5;

    // invariant mass of the leading jet and photon
    m_mass = ((*leading_ph_id)->p4() + (*leading_jet_id)->p4()).M();

    physics->Fill();
    return 0;
}
