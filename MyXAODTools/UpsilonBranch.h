#ifndef __MYXAODTOOLS_UPSILONBRANCH_H__
#define __MYXAODTOOLS_UPSILONBRANCH_H__

#include <vector>
#include <string>
#include <map>

#include "MyXAODTools/BranchCreatorBase.h"

using namespace std;

class UpsilonBranch : public BranchCreatorBase
{
public:
    UpsilonBranch();
    virtual ~UpsilonBranch();

    void AttachBranchToTree(TTree& );
    void ClearBranch();
    bool CreateBranch();
    void Fill();

public:
    static const char* APP_NAME;
    map<string, bool> trigger_map_;

    bool pass_trigger_;

    int event_type_;
    float m_upsilon_;
    float m_4l_;
    float vtx4l_chi2ndf_;
    float m34_;
};

#endif
