#ifndef __MYXAODTOOLS_MUONBRANCH_H__
#define __MYXAODTOOLS_MUONBRANCH_H__

#include <vector>
#include <TLorentzVector.h>

#include "xAODMuon/Muon.h"
#include "xAODMuon/MuonContainer.h"

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
    void Fill(const xAOD::Muon& muon);
public:
    static const char* APP_NAME;
    int total_;
    vector<TLorentzVector>* p4_;
};

#endif
