#include <stdlib.h>

#include "MyXAODTools/MuonBranch.h"
#include "CPAnalysisExamples/errorcheck.h"
#include "xAODTracking/TrackParticlexAODHelpers.h"
#include "MyXAODTools/CPToolsHelper.h"
#include "MuonSelectorTools/MuonSelectionTool.h"

const char* MuonBranch::APP_NAME = "MuonBranch";

MuonBranch::MuonBranch():
    m_track(NULL)
{
    m_isBranchCreated = false;
    m_addDetails = true;
}

int MuonBranch::initial_tools(){
    string toolName("Muon_4Branch");
    m_muonSelectionTool =  unique_ptr<CP::MuonSelectionTool>(new CP::MuonSelectionTool(toolName));
    m_muonSelectionTool->setProperty( "MaxEta", 2.7);
    m_muonSelectionTool->setProperty( "MuQuality", (int) xAOD::Muon::Quality::Medium );
    // m_muonSelectionTool->setProperty( "TurnOffMomCorr", true);
    m_muonSelectionTool->setProperty( "TrtCutOff", true);
    m_muonSelectionTool->msg().setLevel(MSG::ERROR);
    CHECK( m_muonSelectionTool->initialize().isSuccess() );
    return 0;
}

bool MuonBranch::CreateBranch()
{
    m_isBranchCreated = true;
    author_ = new vector<int>();
    type_   = new vector<int>();
    quality_ = new vector<int>;

    pt_     = new vector<float>();
    eta_    = new vector<float>();
    phi_    = new vector<float>();
    e_      = new vector<float>();

    track_charge_     = new vector<float>();
    track_pt_     = new vector<float>();
    track_eta_    = new vector<float>();
    track_phi_    = new vector<float>();
    track_e_      = new vector<float>();

    charge_ = new vector<float>();
    d0_     = new vector<float>();
    z0_sintheta_ = new vector<float>();
    d0_sig_ = new vector<float>();

    eloss_  = new vector<float>();
    etcone30_ = new vector<float>();
    ptvarcone30_ = new vector<float>();

    pvID_ = new vector<int>;

    // variables used in MuonSelectorTools
    if(m_addDetails) {
        m_qOverpSignif.reset( new vector<float>() );
        m_nPrecLayer.reset( new vector<uint8_t>() );
        m_nPrecHoleLayer.reset( new vector<uint8_t>() );
    }

    return true;
}

MuonBranch::~MuonBranch()
{
    if(m_isBranchCreated){
        delete author_;

        delete pt_;
        delete eta_;
        delete phi_;
        delete e_;

        delete track_charge_;
        delete track_pt_;
        delete track_eta_;
        delete track_phi_;
        delete track_e_;

        delete charge_;
        delete type_;
        delete d0_;
        delete z0_sintheta_;
        delete d0_sig_;

        delete eloss_;
        delete etcone30_;
        delete ptvarcone30_;

        delete pvID_;
        delete quality_;
    }
}

void MuonBranch::ClearBranch(){
    total_ = 0;
    if (m_isBranchCreated){
        author_ ->clear();

        pt_ ->clear();
        eta_ ->clear();
        phi_ ->clear();
        e_ ->clear();

        track_charge_ ->clear();
        track_pt_ ->clear();
        track_eta_ ->clear();
        track_phi_ ->clear();
        track_e_ ->clear();

        charge_->clear();
        type_->clear();
        d0_->clear();
        z0_sintheta_->clear();
        d0_sig_->clear();

        eloss_      ->clear();
        etcone30_   ->clear();
        ptvarcone30_->clear();

        pvID_   ->clear();
        quality_->clear();

        if(m_addDetails){
            m_qOverpSignif  ->clear();
            m_nPrecLayer    ->clear();
            m_nPrecHoleLayer->clear();
        }
    }
}

void MuonBranch::AttachBranchToTree(TTree& tree)
{
    if(! m_isBranchCreated) {
        CreateBranch();
        initial_tools();
    }
    tree.Branch("n_muon", &total_, "n_muon/I");

    tree.Branch("mu_author", &author_);
    tree.Branch("mu_pt", &pt_);
    tree.Branch("mu_eta", &eta_);
    tree.Branch("mu_phi", &phi_);
    tree.Branch("mu_e", &e_);


    tree.Branch("mu_track_charge", &track_charge_);
    tree.Branch("mu_track_pt", &track_pt_);
    tree.Branch("mu_track_eta", &track_eta_);
    tree.Branch("mu_track_phi", &track_phi_);
    tree.Branch("mu_track_e", &track_e_);

    tree.Branch("mu_charge", &charge_);
    tree.Branch("mu_type", &type_);
    tree.Branch("mu_d0", &d0_);
    tree.Branch("mu_z0_sintheta", &z0_sintheta_);
    tree.Branch("mu_d0_sig", &d0_sig_);

    tree.Branch("mu_eloss", &eloss_);
    tree.Branch("mu_etcone30", &etcone30_);
    tree.Branch("mu_ptvarcone30", &ptvarcone30_);

    tree.Branch("mu_pvID", &pvID_);
    tree.Branch("mu_quality", &quality_);

    if(m_addDetails){
        tree.Branch("mu_qOverPsig",     (m_qOverpSignif.get()) );
        tree.Branch("mu_nPrecLayer",    (m_nPrecLayer.get()) );
        tree.Branch("mu_nPrecHoleLayer",    (m_nPrecHoleLayer.get()) );
    }
}

void MuonBranch::Fill(const xAOD::Muon& muon,
        const xAOD::EventInfo* evtInfo, const xAOD::Vertex* pv)
{
    total_ ++;
    author_->push_back( (int) muon.author() );
    pt_     ->push_back( muon.pt() );
    eta_    ->push_back( muon.eta() );
    phi_    ->push_back( muon.phi() );
    e_      ->push_back( muon.e() );

    m_track = getTrack(muon);
    if (m_track){
        track_charge_   ->push_back( m_track->charge() );
        track_pt_   ->push_back( m_track->pt() );
        track_eta_  ->push_back( m_track->eta() );
        track_phi_  ->push_back( m_track->phi() );
        track_e_    ->push_back( m_track->e() );
    }

    float charge = (float) muon.charge();
    charge_->push_back(charge);
    type_->push_back((int) muon.muonType());

    double primvertex_z = pv? pv->z(): 0;
    float d0 = -9999;
    float z0_sintheta = -9999; 
    float d0_sig = -9999;
    if(m_track) {
        d0 = m_track->d0();
        // z0_sintheta = (m_track->z0() + m_track->vz() - primvertex_z) * TMath::Sin(muon.p4().Theta());
        z0_sintheta = (m_track->z0() + m_track->vz() - primvertex_z) * TMath::Sin(m_track->p4().Theta());
        d0_sig = xAOD::TrackingHelpers::d0significance(m_track, evtInfo->beamPosSigmaX(), evtInfo->beamPosSigmaY(), evtInfo->beamPosSigmaXY());
    }
    d0_->push_back(d0);
    z0_sintheta_->push_back(z0_sintheta);
    d0_sig_->push_back(d0_sig);

    quality_->push_back( (int)m_muonSelectionTool->getQuality(muon) );
    if (m_addDetails) addDetailedInfo(muon);
}

void MuonBranch::addDetailedInfo(const xAOD::Muon& muon)
{
    m_qOverpSignif  ->push_back( getQoverPsig(muon) );

    uint8_t nprecisionLayers, nprecisionHoleLayers;
    getPrecisionLayer(muon, nprecisionLayers, nprecisionHoleLayers);
    m_nPrecLayer        ->push_back( nprecisionLayers );
    m_nPrecHoleLayer    ->push_back( nprecisionHoleLayers );
}

float MuonBranch::getQoverPsig(const xAOD::Muon& muon) const
{
    const xAOD::TrackParticle* idtrack = muon.trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
    const xAOD::TrackParticle* metrack = muon.trackParticle( xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle );
    if(!idtrack || !metrack) {
        return -999999.;
    }

    float mePt = -999999., idPt = -999999.;
    try{
        static SG::AuxElement::Accessor<float> mePt_acc("MuonSpectrometerPt");
        static SG::AuxElement::Accessor<float> idPt_acc("InnerDetectorPt");
        mePt = mePt_acc(muon);
        idPt = idPt_acc(muon);
    } catch ( SG::ExcNoAuxStore b ) {
        mePt = metrack->pt();
        idPt = idtrack->pt();
    }

    float meP  = 1.0 / ( sin(metrack->theta()) / mePt);
    float idP  = 1.0 / ( sin(idtrack->theta()) / idPt);
    float qOverPsigma  = sqrt( idtrack->definingParametersCovMatrix()(4,4) + metrack->definingParametersCovMatrix()(4,4) );
    float qOverPsignif  = fabs( (metrack->charge() / meP) - (idtrack->charge() / idP) ) / qOverPsigma;

    return qOverPsignif;
}

void MuonBranch::getPrecisionLayer(const xAOD::Muon& muon, uint8_t& precLayer, uint8_t& precHoleLayer) const
{
    if ( !muon.summaryValue(precLayer, xAOD::SummaryType::numberOfPrecisionLayers)){
        precLayer = -1;
    }
    if ( !muon.summaryValue(precHoleLayer, xAOD::SummaryType::numberOfPrecisionHoleLayers) ) {
        precHoleLayer = -1;
    }
    // Fix for CSC
    if( fabs(muon.eta()) > 2.0 ) {
        precLayer = 0;
        uint8_t innerSmallHits, innerLargeHits, middleSmallHits, middleLargeHits, outerSmallHits, outerLargeHits;
        if ( !muon.summaryValue(innerSmallHits, xAOD::MuonSummaryType::innerSmallHits) ||
                !muon.summaryValue(innerLargeHits, xAOD::MuonSummaryType::innerLargeHits) ||
                !muon.summaryValue(middleSmallHits, xAOD::MuonSummaryType::middleSmallHits) ||
                !muon.summaryValue(middleLargeHits, xAOD::MuonSummaryType::middleLargeHits) ||
                !muon.summaryValue(outerSmallHits, xAOD::MuonSummaryType::outerSmallHits) ||
                !muon.summaryValue(outerLargeHits, xAOD::MuonSummaryType::outerLargeHits) )
        {
            precLayer = -1;
        } else{
            if( innerSmallHits>1  || innerLargeHits>1  ) precLayer += 1;
            if( middleSmallHits>2 || middleLargeHits>2 ) precLayer += 1;
            if( outerSmallHits>2  || outerLargeHits>2  ) precLayer += 1;
        }
    }
}

const xAOD::TrackParticle* MuonBranch::getTrack(const xAOD::Muon& muon)
{
    /*** return primary track, not used in 4mu analysis
      const xAOD::TrackParticle* track = NULL;
      if (muon.muonType() == xAOD::Muon::SiliconAssociatedForwardMuon) {
      track = muon.trackParticle(xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle);
      } else {
      track = muon.primaryTrackParticle();
      }
      return track;
     ***/
    // return tracks from inner detector
    return muon.trackParticle(xAOD::Muon::InnerDetectorTrackParticle);
}

int MuonBranch::matchPV(const xAOD::Muon& muon, const xAOD::VertexContainer& vertice)
{
    int res = -1;
    m_track = getTrack(muon);
    if(!m_track){ return res; }

    float mu_vx = m_track->vx();
    float mu_vy = m_track->vy();
    float mu_vz = m_track->vz();

    int index_pv = -1;
    for(const auto& vxp: vertice) {
        index_pv ++;
        float vxt_x = vxp->x();
        float vxt_y = vxp->y();
        float vxt_z = vxp->z();
        if(vxt_x == mu_vx && vxt_y == mu_vy && vxt_z == mu_vz) {
            res = index_pv;
            break;
        }
    }
    return res;
}

void MuonBranch::Fill(const xAOD::Muon& muon,
        const xAOD::EventInfo* evtInfo, const xAOD::VertexContainer* vertice)
{
    const xAOD::Vertex* pv = CPToolsHelper::GetPrimVtx(*vertice);
    this->Fill(muon, evtInfo, pv);

    pvID_->push_back( matchPV(muon, *vertice) );
}
