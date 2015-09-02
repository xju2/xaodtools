#include <stdlib.h>
#include "MyXAODTools/PhotonBranch.h"

PhotonBranch::PhotonBranch()
{
    CreateBranch();
}

PhotonBranch::~PhotonBranch()
{
    if(ph_is_tight_) delete ph_is_tight_;
    if(ph_is_medium_) delete ph_is_medium_;
    if(ph_topoetcone40_) delete ph_topoetcone40_;
    if(ph_p4_) delete ph_p4_;
    if(ph_etaBE_) delete ph_etaBE_;
}

void PhotonBranch::CreateBranch(){
    ph_is_tight_ = new vector<bool>();
    ph_is_medium_ = new vector<bool>();
    ph_topoetcone40_ = new vector<float>();
    ph_p4_ = new vector<TLorentzVector>();
    ph_etaBE_ = new vector<float>();
}

void PhotonBranch::AttachBranchToTree(TTree& tree)
{
    tree.Branch("n_base_ph", &n_base_ph_, "n_base_ph/I");
    tree.Branch("ph_is_tight", &ph_is_tight_);
    tree.Branch("ph_is_medium", &ph_is_medium_);
    tree.Branch("ph_topoetcone40", &ph_topoetcone40_);
    tree.Branch("ph_p4", &ph_p4_);
    tree.Branch("ph_etaBE", &ph_etaBE_);
}

void PhotonBranch::Fill(const xAOD::Photon& photon)
{
    n_base_ph_ ++;
    bool ph_istight = false;
    if ( photon.passSelection(ph_istight,"Tight") )
        ph_is_tight_->push_back( ph_istight );
    /****
    bool ph_isMedium = true;
    if ( photon.passSelection(ph_isMedium,"Medium") )
        ph_is_medium_->push_back( ph_isMedium );
     ***/
    ph_is_medium_->push_back( true );
    float topoetcone40 = 0;
    if ( photon.isolationValue(topoetcone40,xAOD::Iso::topoetcone40) )
        ph_topoetcone40_->push_back( topoetcone40 );
    ph_p4_->push_back(photon.p4());
    ph_etaBE_->push_back(photon.caloCluster()->etaBE(2));
}

void PhotonBranch::ClearBranch()
{
    n_base_ph_ = 0;
    ph_is_tight_->clear();
    ph_is_medium_->clear();
    ph_topoetcone40_->clear();
    ph_p4_->clear();
    ph_etaBE_->clear();
}
