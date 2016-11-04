#ifndef __MYXAODTOOLS_PHOTONBRANCH_H__
#define __MYXAODTOOLS_PHOTONBRANCH_H__

#include <TLorentzVector.h>
#include "xAODEgamma/Photon.h"
#include "MyXAODTools/BranchCreatorBase.h"

#include <vector>
using namespace std;
class PhotonBranch : public BranchCreatorBase
{
public:
    explicit PhotonBranch ();
    ~PhotonBranch();
    void AttachBranchToTree(TTree& );
    void Fill(const xAOD::Photon& photon, bool signal=false);
    void ClearBranch();
    void CreateBranch();

private:
    int n_base_ph_;
    vector<bool>* ph_is_tight_;
    vector<bool>* ph_is_medium_;
    vector<bool>* ph_is_signal_;
    vector<float>* ph_topoetcone40_;
    vector<TLorentzVector>* ph_p4_;
    vector<float>* ph_etaBE_;
};
#endif
