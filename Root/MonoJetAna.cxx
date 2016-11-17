#include <stdlib.h>

#include <TFile.h>

#include "MyXAODTools/MonoJetAna.h"
#include "MyXAODTools/Helper.h"
#include "CPAnalysisExamples/errorcheck.h"
#include "xAODBase/IParticleHelpers.h"

#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/ElectronAuxContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"



MonoJetAna::MonoJetAna():
    AnalysisBase(),
    smearJet("smearjet"),
    m_singleJetTrigger{
        // Order matters, don't change.
        // Don't forget a "j", otherwise you would waste your whole weekend.
        "HLT_j400", "HLT_j380", "HLT_j360", "HLT_j320",
        "HLT_j300", "HLT_j260", "HLT_j200", "HLT_j175",
        "HLT_j150", "HLT_j110", "HLT_j100", "HLT_j85",
        "HLT_j60",  "HLT_j55",  "HLT_j25",  "HLT_j15"
    },
    m_JET_PT_CUT(200),
    m_MET_ET_CUT(200)
{
    m_doSmearing = false;
    if(APP_NAME==NULL) APP_NAME = "MonoJetAna";
    string maindir(getenv("ROOTCOREBIN"));
    m_susy_config = Form("%s/data/MyXAODTools/monojet.conf", maindir.c_str());

    trigger_map_ = {
        {"HLT_xe70", false},
        {"HLT_xe80_tc_lcw_L1XE50",  false},
        {"HLT_xe90_mht_L1XE50",     false},
        {"HLT_xe100_mht_L1XE50",    false},
        {"HLT_xe110_mht_L1XE50",    false},
        // electron
        // {"HLT_e24_lhmedium_L1EM18VH", false}, // 2015
        // {"HLT_e24_lhmedium_L1EM20VH", false}, // 2015
        // {"HLT_e24_lhtight_nod0_ivarloose", false}, // 2016
        // {"HLT_e26_lhtight_nod0_ivarloose", false}, // 2016
        // {"HLT_e60_lhmedium_nod0", false}, // 2016
        // {"HLT_e140_lhloose_nod0", false}, // 2016
        // {"HLT_e60_medium", false}, // 2016
        // // photon
        // {"HLT_g140_loose", false},
        // pre-scaled jet
        {"HLT_j15", false},
        {"HLT_j25", false},
        {"HLT_j55", false},
        {"HLT_j60", false},
        {"HLT_j85", false},
        {"HLT_j100", false},
        {"HLT_j110", false},
        {"HLT_j150", false},
        {"HLT_j175", false},
        {"HLT_j200", false},
        {"HLT_j260", false},
        {"HLT_j300", false},
        {"HLT_j320", false},
        {"HLT_j360", false},
        {"HLT_j380", false},
        {"HLT_j400", false}
    };

}

int MonoJetAna::initialize()
{
    if( initializeBasicTools() != 0 ){
        return 1;
    }
    lightJetResponse = NULL;
    bJetResponse = NULL;
    CreateBranch();
    AttachBranchToTree();

    // initiate tools
    m_jetCleaningTool.reset( new JetCleaningTool("JetCleaningToolTight") );
    CHECK(m_jetCleaningTool->setProperty("CutLevel", "TightBad"));
    CHECK(m_jetCleaningTool->initialize());

    if (m_doSmearing){
        unsigned int test = 1000;
        m_mySmearingTool.reset( new JetSmearing::JetMCSmearingTool("MonoJetSmearingTool") );
        m_mySmearingTool->setProperty("NumberOfSmearedEvents",test);
        m_mySmearingTool->setProperty("DoPhiSmearing", true);
        m_mySmearingTool->initialize();

        string maindir(getenv("ROOTCOREBIN"));
        std::string input_light_jet(maindir+"/data/JetSmearing/MC15/R_map2016_bveto_OP77_EJES_p2666.root");
        TFile* lightJetFile = TFile::Open(input_light_jet.c_str(), "read");
        lightJetResponse = (TH2F*)lightJetFile->Get("responseEJES_p2666");
        lightJetResponse->SetDirectory(0);
        lightJetFile->Close();

        std::string input_bjet(maindir+"/data/JetSmearing/MC15/R_map2016_btag_OP77_EJES_p2666.root");
        TFile* bJetFile = TFile::Open(input_bjet.c_str(), "read");
        bJetResponse = (TH2F*)bJetFile->Get("responseEJES_p2666");
        bJetResponse->SetDirectory(0);
        bJetFile->Close();

        m_mySmearingTool->SetResponseMaps(lightJetResponse, bJetResponse);


        // setup pre-scale tool
        m_prescaleTool.reset( new JetSmearing::PreScaleTool("PreScaleTool") );
        // m_prescaleTool = new JetSmearing::PreScaleTool("PreScaleTool");
        CHECK(m_prescaleTool->setProperty("HistoPath", maindir+"/data/JetSmearing/PreScales/prescale_histos_combined_2015-2016_v81-pro20-10.root"));
        CHECK(m_prescaleTool->initialize());
    }
    return 0;
}

MonoJetAna::~MonoJetAna(){
    if(lightJetResponse) delete lightJetResponse;
    if(bJetResponse) delete bJetResponse;

    if(m_jetP4) delete m_jetP4;
    if(m_doSmearing && m_smearedData){
        delete m_smearedData;
    }
}

void MonoJetAna::CreateBranch()
{
    CreateBasicBranch();
    m_jetP4 = new vector<TLorentzVector>;
    if (m_doSmearing){
        m_smearedData = new vector<SmearedInfo>();
    }
    return ;
}

void MonoJetAna::ClearBranch(){
    ClearBasicBranch();

    m_jetP4->clear();
    m_nGoodJets = 0;
    m_nJetsBtagged = 0;
    m_nBaseEl = 0;
    m_nBaseMu = 0;
    m_nBasePhoton = 0;
    m_hasCosmicMuon = false;
    m_hasBadMuon = false;
    if (m_doSmearing){
        m_smearedData->clear();
    }

}


void MonoJetAna::AttachBranchToTree()
{
    AttachBasicToTree();

    event_br->AttachBranchToTree(*physics);

    physics->Branch("n_good_jet", &m_nGoodJets, "n_good_jet/I");
    physics->Branch("n_jet_btagged", &m_nJetsBtagged, "n_jet_btagged/I");
    physics->Branch("jet_p4", &m_jetP4);
    physics->Branch("triggerWeight", &m_triggerWeight, "triggerWeight/F");
    // physics->Branch("n_base_el",  &m_nBaseEl,  "n_base_el/I");
    // physics->Branch("n_base_mu",  &m_nBaseMu,  "n_base_mu/I");
    // physics->Branch("n_base_ph",  &m_nBasePhoton,  "n_base_ph/I");
    physics->Branch("min_dphi_jetMET", &m_minDphiJetsMET, "min_dphi_jetMET/F");
    physics->Branch("MET_etx", &m_metEtx, "MET_etx/F");
    physics->Branch("MET_ety", &m_metEty, "MET_ety/F");
    physics->Branch("MET_et", &m_metEt, "MET_et/F");
    physics->Branch("MET_sumet", &m_metSumet, "MET_sumet/F");
    physics->Branch("MET_et_soft", &m_metEtSoft, "MET_et_soft/F");
    physics->Branch("has_comic_muon", &m_hasCosmicMuon, "has_comic_muon/O");
    physics->Branch("has_bad_muon", &m_hasBadMuon, "has_bad_muon/O");
    if(m_doSmearing){
        physics->Branch("pseudoData", m_smearedData);
    }
}

int MonoJetAna::process(Long64_t ientry)
{
    int sc = Start(ientry);
    if(m_debug) {
        Info(APP_NAME, " MonoJetAna: processing");
    }
    if(sc != 0) return sc;
    event_br->Fill(*ei);

    // Start mono-jet analysis

    /*get physics objects*/
    // Electrons
    xAOD::ElectronContainer* electrons_copy = NULL;
    xAOD::ShallowAuxContainer* electrons_copyaux = NULL;
    CHECK( m_objTool->GetElectrons(electrons_copy, electrons_copyaux, true) );

    // Muons
    xAOD::MuonContainer* muons_copy = NULL;
    xAOD::ShallowAuxContainer* muons_copyaux = NULL;
    CHECK( m_objTool->GetMuons(muons_copy, muons_copyaux, true) );

    // Jets
    xAOD::JetContainer* jets_copy = NULL;
    xAOD::ShallowAuxContainer* jets_copyaux = NULL;
    CHECK( m_objTool->GetJets(jets_copy,jets_copyaux, true) );

    // Photons
    xAOD::PhotonContainer* ph_copy = nullptr;
    xAOD::ShallowAuxContainer* phAux_copy = nullptr;
    CHECK(m_objTool->GetPhotons(ph_copy, phAux_copy, true));

    ///////////////////////
    // do overlap removal before object selection
    // turn off the harmonization
    ///////////////////////
    // bool doHarmonization = false;
    CHECK( m_objTool->OverlapRemoval(
                electrons_copy, muons_copy,
                jets_copy, ph_copy) );

    //discard the event if any jets is labelled as 'bad'
    bool passJetCleaning = true;
    for(const auto& jet : *jets_copy){
        m_objTool->IsBJet(*jet) ;
        if ( jet->pt() > 20e3 )
        {
            if( dec_bad(*jet) && dec_passOR(*jet)){
                passJetCleaning = false;
                break;
            }
        }
    }
    if ( !passJetCleaning ) return 1;

    // electron selections
    for(const auto&  el : *electrons_copy){
        if( !(bool) dec_baseline(*el) || !(bool) dec_passOR(*el)){
            continue;
        }
        m_nBaseEl ++;
    }
    // muon selections
    for(const auto& mu : *muons_copy){
        if( !(bool) dec_baseline(*mu) || !(bool) dec_passOR(*mu) ){
            continue;
        }
        if( dec_bad(*mu) ) m_hasBadMuon = true;
        if( dec_cosmic(*mu) ) m_hasCosmicMuon = true;
        m_nBaseMu ++;
    }
    // photon selections
    for(const auto& ph : *ph_copy) {
        if( !(bool) dec_baseline(*ph) || !(bool) dec_passOR(*ph) ){
            continue;
        }
        m_nBasePhoton ++;
    }

    // require no baseline electron/muon/photon
    if (m_nBaseEl > 0 || m_nBaseMu > 0 || m_nBasePhoton > 0) {
        if(m_debug){
            Info(APP_NAME, "has lepton/photon in the event %d %d %d",
                    m_nBaseEl, m_nBaseMu, m_nBasePhoton);
        }
        return 2;
    }

    // first calculate MET so that later one cat obtain minimum Dphi
    /* ********
     * MET with muon invisible. 
     * ********/
    auto met = new xAOD::MissingETContainer;
    auto metAux = new xAOD::MissingETAuxContainer;
    met->setStore(metAux);
    // CHECK( store.record( met, "MET_MyRefFinal" ) );
    // CHECK( store.record( metAux, "MET_MyRefFinalAux." ) );
    CHECK( m_objTool->GetMET(*met,
                jets_copy,
                electrons_copy,
                muons_copy, // muon term
                ph_copy, // photon
                nullptr,//taus_copy,
                true,  // TST
                true, // JVT
                0// 0 
                ) );
    auto met_it = met->find("Final");

    if (met_it == met->end() )
    {
        Error( APP_NAME, "No RefFinal inside MET container" );
    }
    m_metEtx = (float) (*met_it)->mpx();
    m_metEty = (float) (*met_it)->mpy();
    double met_phi= (*met_it)->phi();
    m_metEt  = (float) (*met_it)->met();
    m_metSumet  = (float) (*met_it)->sumet();
    m_metEtSoft = (float) (*(met->find("PVSoftTrk")))->met();

    delete met;
    delete metAux;

    // jet selections
    sort(jets_copy->begin(), jets_copy->end(), descend_on_pt);
    double min_dphi = 99999;
    const xAOD::Jet* leading_jet = NULL; 
    for(const auto& jet : *jets_copy)
    {
        if ( !(bool) dec_baseline(*jet) ||  !(bool)dec_passOR(*jet) ){
            continue;
        }
        jet->auxdata<char>(smearJet) = true;

        if(! dec_signal(*jet) ) continue;

        m_jetP4->push_back( jet->p4() );
        m_nGoodJets ++ ;
        if(leading_jet == NULL){
            leading_jet = jet;
        }

        double dphi = fabs(TVector2::Phi_mpi_pi(met_phi - jet->phi()));
        if (dphi < min_dphi) {
            min_dphi = dphi;
        }
        if( dec_bjet(*jet)) m_nJetsBtagged ++;

        // decorate if jet pass tight jet-cleaning
        dec_tightBad(*jet) = (bool) m_jetCleaningTool->keep(*jet);
    }
    m_minDphiJetsMET = min_dphi;

    // at least one leading jet
    if( m_nGoodJets < 1 || m_jetP4->at(0).Pt()/1E3 < m_JET_PT_CUT ||
        dec_tightBad(*leading_jet) != 1 ||
        fabs(m_jetP4->at(0).Eta()) >= 2.4 )
    { 

        if(m_debug){
            if(m_nGoodJets < 1){
                Info(APP_NAME, "no good leading jet");
            } else {
                Info(APP_NAME, "leading-jet pT is low %.2f, pass tightBad %d, eta: %.2f", m_jetP4->at(0).Pt()/1E3, (int)dec_tightBad(*leading_jet), m_jetP4->at(0).Eta() );
            }
        }
        return 3;
    }

    // obtain prescale factor for the event
    if(m_doSmearing){
        // obtain prescale factor
        vector<pair<string, bool> > triggerPass;
        for(const auto sjt: m_singleJetTrigger) {
            triggerPass.push_back(make_pair(sjt, m_objTool->IsTrigPassed(sjt)));
        }
        const xAOD::JetContainer * hlt_jet = 0;
        TString mc_name="HLT_xAOD__JetContainer_a4tcemsubjesFS";
        if( ! event->retrieve( hlt_jet, mc_name.Data()).isSuccess() ) {
            Error("execute()", Form("failed to retrieve %s", mc_name.Data()));
        }
        //use the prescale tool to retrieve the prescale weight for the event
        if(hlt_jet->size()>0){
            m_triggerWeight = m_prescaleTool->getTriggerPrescale(triggerPass,*(hlt_jet->at(0)), ei->runNumber());
        } else {
            m_triggerWeight = 0.0;
        }

        // !!!start to smear the data!!!!
        std::vector<std::unique_ptr<SmearData > > smrMc;
        m_mySmearingTool->DoSmearing(smrMc,*jets_copy);
        for (auto& SmearedEvent : smrMc){
            // clear up smeared data, or the program will consume infinite memory!
            xAOD::JetContainer* theJetContainer =  SmearedEvent->jetContainer;
            SmearedInfo smeared;
            if (get_smeared_info(theJetContainer, muons_copy, electrons_copy, ph_copy, smeared)
               ){
                m_smearedData->push_back(smeared);
            }
        }
    }
    // Fill your tree!!!!
    physics->Fill();
    return 0;
}

bool MonoJetAna::get_smeared_info(
        xAOD::JetContainer* jets,
        xAOD::MuonContainer* muons,
        xAOD::ElectronContainer* electrons,
        xAOD::PhotonContainer* photons,
        SmearedInfo& smeared_info)
{
    sort(jets->begin(), jets->end(), descend_on_pt);
    smeared_info.leading_jet_pt_ = (float)jets->at(0)->p4().Pt();
    smeared_info.leading_jet_eta_ = (float)jets->at(0)->p4().Eta();
    smeared_info.leading_jet_phi_ = (float)jets->at(0)->p4().Phi();
    if( smeared_info.leading_jet_pt_/1E3 < m_JET_PT_CUT ||
        fabs(smeared_info.leading_jet_eta_) >= 2.4 ||
        dec_tightBad(*(jets->at(0))) != 1 ||
        jets->at(0)->auxdata< char >(smearJet) == false
    ) return false;

    unique_ptr<xAOD::MissingETContainer> met(new xAOD::MissingETContainer() );
    unique_ptr<xAOD::MissingETAuxContainer> metAux(new xAOD::MissingETAuxContainer() );
    met->setStore(metAux.get());
    m_objTool->GetMET(*met,
                jets,
                electrons,
                muons, // muon term
                photons, // photon
                nullptr,//taus_copy,
                true,  // TST
                true, // JVT
                0 // 0 
                );
    xAOD::MissingETContainer::const_iterator met_it = met->find("Final");

    smeared_info.met_ =(float) (*met_it)->met();
    smeared_info.sum_et_ =(float) (*met_it)->sumet();
    // since xe80 used for the analysis, cut on 80 GeV to reduce size of pseudo-data.
    if(smeared_info.met_/1E3 < m_MET_ET_CUT){
        return false;
    }
    float min_dphi_jetMET  = 9999;
    int n_good_jets = 0;
    smeared_info.HT_ = 0.0;
    smeared_info.sub_leading_jet_pt_ = -999.0;
    smeared_info.sub_leading_jet_eta_ = -999.0;
    smeared_info.sub_leading_jet_phi_ = -999.0;
    smeared_info.l3rd_jet_pt_ = -999.0;
    smeared_info.l3rd_jet_eta_ = -999.0;
    smeared_info.l3rd_jet_phi_ = -999.0;
    smeared_info.l4th_jet_pt_ = -999.0;
    smeared_info.l4th_jet_eta_ = -999.0;
    smeared_info.l4th_jet_phi_ = -999.0;
    for(auto jet: *jets) {
        if ( jet->auxdata< char >(smearJet) == false ) continue;
        if(jet->pt() <= 30E3 || fabs(jet->eta()) >= 2.8)  continue;
        smeared_info.HT_ += jet->pt();
        float dphi = (float) fabs(TVector2::Phi_mpi_pi((*met_it)->phi() - jet->phi()));
        n_good_jets ++;
        if (n_good_jets == 2) {
            smeared_info.sub_leading_jet_pt_ = (float)jet->p4().Pt();
            smeared_info.sub_leading_jet_eta_ = (float)jet->p4().Eta();
            smeared_info.sub_leading_jet_phi_ = (float)jet->p4().Phi();
        } else if(n_good_jets == 3){
            smeared_info.l3rd_jet_pt_ = (float)jet->p4().Pt();
            smeared_info.l3rd_jet_eta_ = (float)jet->p4().Eta();
            smeared_info.l3rd_jet_phi_ = (float)jet->p4().Phi();
        } else if (n_good_jets == 4){
            smeared_info.l4th_jet_pt_ = (float)jet->p4().Pt();
            smeared_info.l4th_jet_eta_ = (float)jet->p4().Eta();
            smeared_info.l4th_jet_phi_ = (float)jet->p4().Phi();
        } else {}
        if(dphi < min_dphi_jetMET) min_dphi_jetMET = dphi;
    }
    smeared_info.min_jets_met_ = min_dphi_jetMET;
    smeared_info.n_good_jets_ = n_good_jets;

    /* Track Missing Et */
    unique_ptr<xAOD::MissingETContainer> met_track(new xAOD::MissingETContainer);
    unique_ptr<xAOD::MissingETAuxContainer> metAux_track(new xAOD::MissingETAuxContainer);
    met_track->setStore(metAux_track.get());
    m_objTool->GetTrackMET(*met_track,
                jets,
                electrons,
                muons
                );
    xAOD::MissingETContainer::const_iterator met_track_it = met_track->find("Track");
    smeared_info.dphi_EP_ = fabs(TVector2::Phi_mpi_pi((*met_track_it)->phi() - (*met_it)->phi()));

    return true;
}

void MonoJetAna::setSmear(bool smear_){
    m_doSmearing = smear_;
}
