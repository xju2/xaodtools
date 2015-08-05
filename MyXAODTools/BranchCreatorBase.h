#ifndef __MYXAODTOOLS_BRANCHCREATORBASE_H__
#define __MYXAODTOOLS_BRANCHCREATORBASE_H__

#include "TTree.h"
#include "xAODBase/IParticle.h"
#include <string>

using namespace std;
class BranchCreatorBase{
public:
    BranchCreatorBase();
    virtual ~BranchCreatorBase();
    virtual void AttachBranchToTree(TTree& ){}
    // virtual void Fill(const xAOD::IParticle& ){}
    virtual void ClearBranch(){}
    void CreateBranch() {}

private:
    // string name_;
};
#endif
