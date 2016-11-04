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
    void AttachBranchToTree(TTree& ){}
    void ClearBranch(){}
    void CreateBranch() { m_isBranchCreated = true; }

protected:
    bool m_isBranchCreated;
};
#endif
