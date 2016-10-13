#include <stdlib.h>

#include "MyXAODTools/UpsilonAna.h"
#include "MyXAODTools/Helper.h"
#include "CPAnalysisExamples/errorcheck.h"

#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/ElectronAuxContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"


UpsilonAna::UpsilonAna(): AnalysisBase()
{
    if(APP_NAME==NULL) APP_NAME = "UpsilonAna";
    string maindir(getenv("ROOTCOREBIN"));
    m_susy_config = Form("%s/data/MyXAODTools/upsilon.conf", maindir.c_str());

    trigger_map_ = {
        // single Muon
        {"HLT_mu20_iloose_L1MU15", false},
        {"HLT_mu24_ivarloose_L1MU15", false},
        {"HLT_mu24_ivarmedium", false},
        {"HLT_mu24_imedium", false},
        {"HLT_mu26_ivarmedium", false},
        {"HLT_mu26_imedium", false},
        // Di muon
        {"HLT_2mu10", false},
        {"HLT_mu18_mu8noL1", false},
        {"HLT_2mu10_nomucomb", false},
        {"HLT_mu20_mu8noL1", false},
        {"HLT_mu20_nomucomb_mu6noL1_nscan03", false},
        {"HLT_2mu14_nomucomb", false},
        {"HLT_mu22_mu8noL1", false},
        {"HLT_2mu14", false},
        // Trig Muons
        {"HLT_3mu6", false},
        {"HLT_3mu6_msonly", false},
        {"HLT_mu18_2mu4noL1", false},
        {"HLT_mu20_2mu4noL1", false},
        {"HLT_3mu4", false},
        {"HLT_mu6_2mu4", false},
        {"HLT_mu11_nomucomb_2mu4noL1_nscan03_L1MU11_2MU6", false},
        {"HLT_mu20_msonly_mu10noL1_msonly_nscan05_noComb", false},
        {"HLT_3mu6", false},
        {"HLT_3mu6_msonly", false}
    };
    CreateBranch();
    AttachBranchToTree();
}


UpsilonAna::~UpsilonAna(){
}

void UpsilonAna::CreateBranch()
{
    // Onia information
    m_onia_muon1id = new std::vector<int>;
    m_onia_muon2id = new std::vector<int>;
    m_onia_charge = new std::vector<float>;

    m_onia_pt_fitted = new std::vector<float>;
    m_onia_eta_fitted = new std::vector<float>;
    m_onia_phi_fitted = new std::vector<float>;
    m_onia_mass_fitted = new std::vector<float>;
    m_onia_x = new std::vector<float>;
    m_onia_y = new std::vector<float>;
    m_onia_z = new std::vector<float>;
    m_onia_chi2 = new std::vector<float>;

    m_onia_mass   = new std::vector<float>;
    m_onia_pt     = new std::vector<float>;
    m_onia_eta    = new std::vector<float>;
    m_onia_phi    = new std::vector<float>;

    m_onia_track_mass   = new std::vector<float>;
    m_onia_track_pt     = new std::vector<float>;
    m_onia_track_eta    = new std::vector<float>;
    m_onia_track_phi    = new std::vector<float>;

    // upsilon information
    m_quad_charge = new std::vector<float>;
    m_quad_chi2 = new std::vector<float>;
    m_quad_x      = new std::vector<float>;
    m_quad_y      = new std::vector<float>;
    m_quad_z      = new std::vector<float>;

    m_quad_nCombined = new std::vector<int>;
    m_quad_id1 = new std::vector<int>;
    m_quad_id2 = new std::vector<int>;
    m_quad_id3 = new std::vector<int>;
    m_quad_id4 = new std::vector<int>;

    m_quad_mass   = new std::vector<float>;
    m_quad_pt     = new std::vector<float>;
    m_quad_eta    = new std::vector<float>;
    m_quad_phi    = new std::vector<float>;

    m_quad_track_mass   = new std::vector<float>;
    m_quad_track_pt     = new std::vector<float>;
    m_quad_track_eta    = new std::vector<float>;
    m_quad_track_phi    = new std::vector<float>;

    m_quad_fitted_mass   = new std::vector<float>;
    m_quad_fitted_pt     = new std::vector<float>;
    m_quad_fitted_eta    = new std::vector<float>;
    m_quad_fitted_phi    = new std::vector<float>;

    return ;
}

void UpsilonAna::ClearBranch(){
    AnalysisBase::ClearBranch();

    // Onia information
    m_n_onia = 0;
    m_onia_muon1id->clear();
    m_onia_muon2id->clear();
    m_onia_charge ->clear();

    m_onia_pt_fitted->clear();
    m_onia_eta_fitted->clear();
    m_onia_phi_fitted->clear();
    m_onia_mass_fitted->clear();
    m_onia_x->clear();
    m_onia_y->clear();
    m_onia_z->clear();
    m_onia_chi2->clear();

    m_onia_mass->clear();
    m_onia_pt->clear();
    m_onia_eta->clear();
    m_onia_phi->clear();
    m_onia_track_mass->clear();
    m_onia_track_pt->clear();
    m_onia_track_eta->clear();
    m_onia_track_phi->clear();

    // 
    m_n_quad = 0;
    m_quad_charge->clear();
    m_quad_chi2->clear();
    m_quad_x->clear();
    m_quad_y->clear();
    m_quad_z->clear();

    m_quad_nCombined->clear();
    m_quad_id1->clear();
    m_quad_id2->clear();
    m_quad_id3->clear();
    m_quad_id4->clear();

    m_quad_mass->clear();
    m_quad_pt->clear();
    m_quad_eta->clear();
    m_quad_phi->clear();

    m_quad_track_mass->clear();
    m_quad_track_pt->clear();
    m_quad_track_eta->clear();
    m_quad_track_phi->clear();

    m_quad_fitted_mass->clear();
    m_quad_fitted_pt->clear();
    m_quad_fitted_eta->clear();
    m_quad_fitted_phi->clear();
}


void UpsilonAna::AttachBranchToTree()
{
    AnalysisBase::AttachBranchToTree();

    event_br->AttachBranchToTree(*physics);
    muon_br ->AttachBranchToTree(*physics);

    // Onia information
    physics->Branch("n_onia", &m_n_onia, "n_onia/I");
    physics->Branch("onia_id1", &m_onia_muon1id);
    physics->Branch("onia_id2", &m_onia_muon2id);
    physics->Branch("onia_charge", &m_onia_charge);

    physics->Branch("onia_fitted_pt", &m_onia_pt_fitted);
    physics->Branch("onia_fitted_eta", &m_onia_eta_fitted);
    physics->Branch("onia_fitted_phi", &m_onia_phi_fitted);
    physics->Branch("onia_fitted_mass", &m_onia_mass_fitted);

    physics->Branch("onia_x", &m_onia_x);
    physics->Branch("onia_y", &m_onia_y);
    physics->Branch("onia_z", &m_onia_z);
    physics->Branch("onia_chi2", &m_onia_chi2);

    physics->Branch("onia_mass", &m_onia_mass);
    physics->Branch("onia_pt", &m_onia_pt);
    physics->Branch("onia_eta", &m_onia_eta);
    physics->Branch("onia_phi", &m_onia_phi);

    physics->Branch("onia_track_mass", &m_onia_track_mass);
    physics->Branch("onia_track_pt", &m_onia_track_pt);
    physics->Branch("onia_track_eta", &m_onia_track_eta);
    physics->Branch("onia_track_phi", &m_onia_track_phi);

    //physics->Branch("has_upsilon", &has_upsilon, "has_upsilon/O");

    // upsilon 
    physics->Branch("n_quad", &m_n_quad, "n_quad/I");
    physics->Branch("quad_charge", &m_quad_charge);
    physics->Branch("quad_chi2", &m_quad_chi2);
    physics->Branch("quad_x", &m_quad_x);
    physics->Branch("quad_y", &m_quad_y);
    physics->Branch("quad_z", &m_quad_z);
    physics->Branch("quad_nCombined", &m_quad_nCombined);
    physics->Branch("quad_id1", &m_quad_id1);
    physics->Branch("quad_id2", &m_quad_id2);
    physics->Branch("quad_id3", &m_quad_id3);
    physics->Branch("quad_id4", &m_quad_id4);

    physics->Branch("quad_mass", &m_quad_mass);
    physics->Branch("quad_pt", &m_quad_pt);
    physics->Branch("quad_eta", &m_quad_eta);
    physics->Branch("quad_phi", &m_quad_phi);
    physics->Branch("quad_track_mass", &m_quad_track_mass);
    physics->Branch("quad_track_pt", &m_quad_track_pt);
    physics->Branch("quad_track_eta", &m_quad_track_eta);
    physics->Branch("quad_track_phi", &m_quad_track_phi);
    physics->Branch("quad_fitted_mass", &m_quad_fitted_mass);
    physics->Branch("quad_fitted_pt", &m_quad_fitted_pt);
    physics->Branch("quad_fitted_eta", &m_quad_fitted_eta);
    physics->Branch("quad_fitted_phi", &m_quad_fitted_phi);
}

int UpsilonAna::process(Long64_t ientry)
{
    int sc = AnalysisBase::process(ientry);
    if(m_debug) {
        Info(APP_NAME, " UpsilonAna: processing");
    }
    if(sc != 0) return sc;
    event_br->Fill(*ei);
   
    m_bphy4_quad = NULL;
    m_bphy4_pair = NULL;
    // obtain the fitted BPHY4
    if(event->contains<xAOD::VertexContainer>("BPHY4Quads")) {
        CHECK(event->retrieve(m_bphy4_quad, "BPHY4Quads"));
    }
    if(event->contains<xAOD::VertexContainer>("BPHY4Pairs")){
        CHECK(event->retrieve(m_bphy4_pair, "BPHY4Pairs"));
    }

    // get muons
    xAOD::MuonContainer* muons_copy = NULL;
    xAOD::ShallowAuxContainer* muons_copyaux = NULL;
    CHECK( m_objTool->GetMuons(muons_copy, muons_copyaux, true) );
    sort(muons_copy->begin(), muons_copy->end(), descend_on_pt);

    int n_muon = 0;
    int imuon = -1;
    MuonVect* good_muons = new MuonVect();
    int n_pos = 0;
    int n_neg = 0;

    for(auto mu_itr = muons_copy->begin(); mu_itr != muons_copy->end(); ++mu_itr)
    {
        imuon ++;
        if( (*mu_itr)->muonType() != xAOD::Muon::Combined &&
            (*mu_itr)->muonType() != xAOD::Muon::SegmentTagged ) continue;

        const xAOD::TrackParticle* id_track = MuonBranch::getTrack( (**mu_itr) );
        if(id_track){
            if(fabs(id_track->eta()) > 2.5) continue;
            if(id_track->p4().Pt() < 3) continue;
        }
        float charge = (*mu_itr)->charge();
        bool consitent_charge = (id_track->charge() == charge);
        if(!consitent_charge) continue;

        if( (bool) dec_baseline(**mu_itr) ){
            good_muons->push_back( (*mu_itr) );
            n_muon ++;
            if(charge < 0) n_neg ++;
            else n_pos ++;

            muon_br->Fill(**mu_itr, ei, vertice);
            int muIndex = (*mu_itr)->auxdataConst<int>("BPHY4MuonIndex");

            if(m_debug){
                cout << "Type: " << (*mu_itr)->muonType() << endl;
                cout << "Index: " << muIndex << endl;
                cout << "pT: " << (*mu_itr)->p4().Pt()/1E3 << endl;
            }
            muon_br->eloss_->push_back( (*mu_itr)->auxdataConst<float>("EnergyLoss") );
            muon_br->etcone30_->push_back( (*mu_itr)->auxdataConst<float>("etcone30") );
            muon_br->ptvarcone30_->push_back( (*mu_itr)->auxdataConst<float>("ptvarcone30") );
        }
    }
    //if (n_pos >= 2 && n_neg >= 2)
    if (n_muon >= 4)
    {
        this->buildTwoMuons( *good_muons );
        this->buildFourMuons( *good_muons );
        physics->Fill();
    }
    delete good_muons;

    return 0;
}

void UpsilonAna::buildTwoMuons(const MuonVect& muons)
{
    if(muons.size() < 2) return;

    for(int i = 0; i < (int) muons.size(); i++) {
        const xAOD::Muon* muon1 = dynamic_cast<const xAOD::Muon*>( muons.at(i) );
        float mu_charge_1 = muon1->charge();

        for(int j = i+1; j < (int) muons.size(); ++j) {
            const xAOD::Muon* muon2 = dynamic_cast<const xAOD::Muon*>( muons.at(j) );
            float mu_charge_2 = muon2->charge();
            // if( (mu_charge_1 + mu_charge_2 ) != 0) continue;

            this->fillOniaInfo(*muon1, *muon2);
            m_onia_muon1id->push_back(i);
            m_onia_muon2id->push_back(j);
        }
    }
}

void UpsilonAna::buildFourMuons(const MuonVect& muons)
{
    if(muons.size() < 4) return;
    for(int i = 0; i < (int) muons.size(); ++i) {
        const xAOD::Muon* muon1 = dynamic_cast<const xAOD::Muon*>( muons.at(i) );

        for(int j=i+1; j < (int) muons.size(); j++){
            const xAOD::Muon* muon2 = dynamic_cast<const xAOD::Muon*>( muons.at(j) );

            for( int k=j+1; k < (int) muons.size(); ++k){
                const xAOD::Muon* muon3 = dynamic_cast<const xAOD::Muon*>( muons.at(k) );

                for(int l=k+1; l < (int) muons.size(); ++l){
                    const xAOD::Muon* muon4 = dynamic_cast<const xAOD::Muon*>( muons.at(l) );
                    // if( (muon1->charge() + muon2->charge() + muon3->charge() + muon4->charge()) != 0) continue;
                    // so far we have four neutral tracks

                    // require at least three combined muons
                    int n_combined = 0;
                    if( muon1->muonType() == xAOD::Muon::Combined ) n_combined ++;
                    if( muon2->muonType() == xAOD::Muon::Combined ) n_combined ++;
                    if( muon3->muonType() == xAOD::Muon::Combined ) n_combined ++;
                    if( muon4->muonType() == xAOD::Muon::Combined ) n_combined ++;
                    // if(n_combined < 3) continue;

                    m_quad_nCombined->push_back(n_combined);
                    this->fillQuadInfo(*muon1, *muon2, *muon3, *muon4);
                    m_quad_id1->push_back(i);
                    m_quad_id2->push_back(j);
                    m_quad_id3->push_back(k);
                    m_quad_id4->push_back(l);
                }
            }
        }
    }
}

void UpsilonAna::fillOniaInfo(const xAOD::Muon& muon1, const xAOD::Muon& muon2)
{
    m_n_onia ++;
    m_onia_charge->push_back(muon1.charge() + muon2.charge());

    MuonVect* muons = new MuonVect();
    muons->push_back( &muon1 );
    muons->push_back( &muon2 );
    
    float onia_mass_fitted = -999;
    float onia_x = -999;
    float onia_y = -999;
    float onia_z = -999;
    float onia_chi2 = -999;

    if(m_bphy4_pair) {
        const xAOD::Vertex* fitted_v = matchFittedVertex(*m_bphy4_pair, *muons);
        if(fitted_v){
            onia_mass_fitted = fitted_v->auxdataConst<float>("PAIR_mass");
            onia_x = fitted_v->x();
            onia_y = fitted_v->y();
            onia_z = fitted_v->z();
            onia_chi2 = fitted_v->chiSquared() / fitted_v->numberDoF();
        }
    }

    m_onia_x->push_back( onia_x );
    m_onia_y->push_back( onia_y );
    m_onia_z->push_back( onia_z );
    m_onia_chi2->push_back( onia_chi2 );
    m_onia_mass_fitted->push_back( onia_mass_fitted );

    // four-momentum from combined measurement 
    TLorentzVector tlv_total = (muon1.p4() + muon2.p4());
    m_onia_pt   ->push_back( tlv_total.Pt() );
    m_onia_eta  ->push_back( tlv_total.Eta() );
    m_onia_phi  ->push_back( tlv_total.Phi() );
    m_onia_mass ->push_back( tlv_total.M() );

    // four-momentum from track particle
    const xAOD::TrackParticle* track1 = MuonBranch::getTrack(muon1);
    const xAOD::TrackParticle* track2 = MuonBranch::getTrack(muon2);
    TLorentzVector track_tlv_total = track1->p4() + track2->p4();
    m_onia_track_pt->push_back( track_tlv_total.Pt() );
    m_onia_track_eta->push_back( track_tlv_total.Eta() );
    m_onia_track_phi->push_back( track_tlv_total.Phi() );
    m_onia_track_mass->push_back( track_tlv_total.M() );
    
    delete muons;
}

void UpsilonAna::fillQuadInfo(
        const xAOD::Muon& muon1, const xAOD::Muon& muon2,
        const xAOD::Muon& muon3, const xAOD::Muon& muon4)
{
    m_n_quad ++;
    m_quad_charge->push_back(muon1.charge() + muon2.charge() + muon3.charge() + muon4.charge());

    MuonVect* muons_can = new MuonVect;
    muons_can->push_back( &muon1 );
    muons_can->push_back( &muon2 );
    muons_can->push_back( &muon3 );
    muons_can->push_back( &muon4 );

    float quad_mass_fitted = -999;
    float quad_x = -999;
    float quad_y = -999;
    float quad_z = -999;
    float quad_chi2 = -999;

    if(m_bphy4_quad) {
        const xAOD::Vertex* fitted_v = matchFittedVertex(*m_bphy4_quad, *muons_can);
        if(fitted_v){
            quad_mass_fitted = fitted_v->auxdataConst<float>("QUAD_mass");
            quad_x = fitted_v->x();
            quad_y = fitted_v->y();
            quad_z = fitted_v->z();
            quad_chi2 = fitted_v->chiSquared() / fitted_v->numberDoF();
        }
    }
    m_quad_x    ->push_back( quad_x );
    m_quad_y    ->push_back( quad_y );
    m_quad_z    ->push_back( quad_z );
    m_quad_chi2 ->push_back( quad_chi2 );
    m_quad_fitted_mass->push_back( quad_mass_fitted );

    delete muons_can;

    // four-momentum from combined measurement 
    TLorentzVector tlv_total = (muon1.p4() + muon2.p4() + muon3.p4() + muon4.p4());
    m_quad_pt   ->push_back( tlv_total.Pt() );
    m_quad_eta  ->push_back( tlv_total.Eta() );
    m_quad_phi  ->push_back( tlv_total.Phi() );
    m_quad_mass ->push_back( tlv_total.M() );

    // four-momentum from track particle
    const xAOD::TrackParticle* track1 = MuonBranch::getTrack(muon1);
    const xAOD::TrackParticle* track2 = MuonBranch::getTrack(muon2);
    const xAOD::TrackParticle* track3 = MuonBranch::getTrack(muon3);
    const xAOD::TrackParticle* track4 = MuonBranch::getTrack(muon4);

    TLorentzVector track_tlv_total = (track1->p4()+ track2->p4() + track3->p4() + track4->p4());
    m_quad_track_pt ->push_back( track_tlv_total.Pt() );
    m_quad_track_eta->push_back( track_tlv_total.Eta() );
    m_quad_track_phi->push_back( track_tlv_total.Phi() );
    m_quad_track_mass->push_back( track_tlv_total.M() );
}

const xAOD::Vertex* UpsilonAna::matchFittedVertex(
        const xAOD::VertexContainer& muonVertexCont,
        const MuonVect& muons)
{
    if(m_debug){
        Info(APP_NAME, "in matchFittedVertex");
    }
    const xAOD::Vertex* res_vertex = NULL;
    vector<const xAOD::TrackParticle*> inputTracks(0);
    // vector<ElementLink<xAOD::TrackParticleContainer> > inputTrackLinks(0);
    for(const auto& muon : muons) {
        // const xAOD::TrackParticle* tp = MuonBranch::getTrack( *muon );
        const xAOD::TrackParticle* tp = muon->trackParticle(xAOD::Muon::InnerDetectorTrackParticle);
        if(tp) { inputTracks.push_back(tp); }
    }
    unsigned int nTracks = inputTracks.size();

    for(const auto& v: muonVertexCont) {
        if(!v) continue;

        set<const xAOD::TrackParticle*> vtxTrks;
        for(unsigned int i = 0; i < v->nTrackParticles(); ++i) {
            vtxTrks.insert(v->trackParticle(i));
        }

        unsigned int nMatch = 0;
        for(unsigned int i = 0; i < inputTracks.size(); i++){
            if(vtxTrks.find(inputTracks.at(i)) != vtxTrks.end()) {
                nMatch ++;
            }
        }
        if(nMatch == nTracks) {
            res_vertex = v;
            break;
        }
    }
    if(m_debug && !res_vertex){
        Info(APP_NAME, "cannot find a mattched vertex!");
    }
    return res_vertex;
}
