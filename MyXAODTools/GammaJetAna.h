#ifndef __MYXAODTOOLS_GAMMAJETANA_H__
#define __MYXAODTOOLS_GAMMAJETANA_H__

#include <vector>
#include <string>
#include <map>

#include <TH1F.h>

#include "MyXAODTools/AnalysisBase.h"

using namespace std;

class GammaJetAna : public AnalysisBase
{
public:
    GammaJetAna();
    virtual ~GammaJetAna();

    int initialize();
    void ClearBranch();
    int process(Long64_t ientry); // main program

private:
    void AttachBranchToTree();
    void CreateBranch();

private:
    const float LEADING_PHOTON_CUT;
    const float LEADING_JET_CUT;

    TH1F* h_cutflow;
    float m_mass;
};

#endif
