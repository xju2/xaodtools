#include <stdlib.h>

#include "MyXAODTools/MuonBranch.h"
#include "CPAnalysisExamples/errorcheck.h"
#include "xAODTracking/TrackParticlexAODHelpers.h"

const char* MuonBranch::APP_NAME = "MuonBranch";

MuonBranch::MuonBranch(){
    CreateBranch();
}

bool MuonBranch::CreateBranch()
{
    p4_ = new vector<TLorentzVector>();
    charge_ = new vector<float>();
    type_ = new vector<int>();
    d0_ = new vector<float>();
    z0_sintheta_ = new vector<float>();
    d0_sig_ = new vector<float>();

    eloss_ = new vector<float>();
    etcone30_ = new vector<float>();
    ptvarcone30_ = new vector<float>();

    return true;
}

MuonBranch::~MuonBranch(){
    delete p4_;
    delete charge_;
    delete type_;
    delete d0_;
    delete z0_sintheta_;
    delete d0_sig_;

    delete eloss_;
    delete etcone30_;
    delete ptvarcone30_;
}

void MuonBranch::ClearBranch(){
    total_ = 0;
    p4_->clear();
    charge_->clear();
    type_->clear();
    d0_->clear();
    z0_sintheta_->clear();
    d0_sig_->clear();

    eloss_->clear();
    etcone30_->clear();
    ptvarcone30_->clear();
}

void MuonBranch::AttachBranchToTree(TTree& tree){
    tree.Branch("n_muon", &total_, "n_muon/I");
    tree.Branch("mu_p4", &p4_);
    tree.Branch("mu_charge", &charge_);
    tree.Branch("mu_type", &type_);
    tree.Branch("mu_d0", &d0_);
    tree.Branch("mu_z0_sintheta", &z0_sintheta_);
    tree.Branch("mu_d0_sig", &d0_sig_);

    tree.Branch("mu_eloss", &eloss_);
    tree.Branch("mu_etcone30", &etcone30_);
    tree.Branch("mu_ptvarcone30", &ptvarcone30_);
}

void MuonBranch::Fill(const xAOD::Muon& muon,
        const xAOD::EventInfo* evtInfo, const xAOD::Vertex* pv)
{
    total_ ++;
    p4_->push_back(muon.p4());
    float charge = (float) muon.charge();
    charge_->push_back(charge);
    type_->push_back((int) muon.muonType());

    const xAOD::TrackParticle* track;
    if (muon.muonType() == xAOD::Muon::SiliconAssociatedForwardMuon) {
        track = muon.trackParticle(xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle);
        // if (!track) return StatusCode::SUCCESS; // don't treat SAF muons without ME track further
    } else {
        track = muon.primaryTrackParticle();
    }
    
    double primvertex_z = pv? pv->z(): 0;
    float d0 = -9999;
    float z0_sintheta = -9999; 
    float d0_sig = -9999;
    if(track) {
        d0 = track->d0();
        z0_sintheta = (track->z0() + track->vz() - primvertex_z) * TMath::Sin(muon.p4().Theta());
        d0_sig = xAOD::TrackingHelpers::d0significance(track, evtInfo->beamPosSigmaX(), evtInfo->beamPosSigmaY(), evtInfo->beamPosSigmaXY());
    }
    d0_->push_back(d0);
    z0_sintheta_->push_back(z0_sintheta);
    d0_sig_->push_back(d0_sig);
}
