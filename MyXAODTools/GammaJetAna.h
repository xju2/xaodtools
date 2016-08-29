#ifndef __MYXAODTOOLS_GAMMAJETANA_H__
#define __MYXAODTOOLS_GAMMAJETANA_H__

#include <vector>
#include <string>
#include <map>

#include "MyXAODTools/AnalysisBase.h"

using namespace std;

class GammaJetAna : public AnalysisBase
{
public:
    GammaJetAna();
    virtual ~GammaJetAna();

    void AttachBranchToTree();
    void CreateBranch();
    void ClearBranch();

    int process(Long64_t ientry); // main program

private:
    int event_type_;
};

#endif
