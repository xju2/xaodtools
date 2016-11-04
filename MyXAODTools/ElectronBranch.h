#ifndef __MYXAODTOOLS_ELECTRONBRANCH_H__
#define __MYXAODTOOLS_ELECTRONBRANCH_H__

#include <vector>
#include <TLorentzVector.h>

#include "xAODEgamma/Electron.h"
#include "xAODEgamma/ElectronContainer.h"

#include "MyXAODTools/BranchCreatorBase.h"
using namespace std;

class ElectronBranch : public BranchCreatorBase
{
public:
    ElectronBranch();
    virtual ~ElectronBranch();

    void AttachBranchToTree(TTree& );
    void ClearBranch();
    bool CreateBranch();
    void Fill(const xAOD::Electron& ele);
public:
    static const char* APP_NAME;
    int total_;
    vector<TLorentzVector>* p4_;
};

#endif
