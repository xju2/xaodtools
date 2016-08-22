#ifndef __MYXAODTOOLS_MUONBRANCH_H__
#define __MYXAODTOOLS_MUONBRANCH_H__

#include <vector>
#include <TLorentzVector.h>

#include "xAODMuon/Muon.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/Vertex.h"


#include "MyXAODTools/BranchCreatorBase.h"
using namespace std;

class MuonBranch : public BranchCreatorBase
{
public:
    MuonBranch();
    virtual ~MuonBranch();

    void AttachBranchToTree(TTree& );
    void ClearBranch();
    bool CreateBranch();
    void Fill(const xAOD::Muon& muon, const xAOD::EventInfo* evtInfo,
            const xAOD::Vertex* pv = NULL);

public:
    static const char* APP_NAME;

    int total_;
    vector<TLorentzVector>* p4_;
    vector<float>* charge_;
    vector<int>* type_;
    vector<float>* d0_;
    vector<float>* z0_sintheta_;
    vector<float>* d0_sig_;

    // These variables have to be filled seperately
    vector<float>* eloss_;
    vector<float>* etcone30_;
    // TODO: add topoetcone30
    vector<float>* ptvarcone30_;
};

#endif
