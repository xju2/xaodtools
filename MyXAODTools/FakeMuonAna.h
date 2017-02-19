
/*
 * Study the fake rate of good tracks 
 *
 */
#ifndef __MYXAODTOOLS_FakeMuonAna_H__
#define __MYXAODTOOLS_FakeMuonAna_H__

#include <vector>
#include <string>

#include "MyXAODTools/AnalysisBase.h"
#include "AsgTools/ToolHandle.h"

#include "xAODMuon/MuonContainer.h"

using namespace std;

class FakeMuonAna : public AnalysisBase
{
public:
    FakeMuonAna();
    virtual ~FakeMuonAna();

    int initialize();
    void ClearBranch();
    int process(Long64_t ientry); // main program

private:
    /** private methods */
    void AttachBranchToTree();
    void CreateBranch();
    bool isPrompt(const xAOD::Muon& muon) const;

private:
    /* specific branches used in this analysis */
    unique_ptr< vector<bool> > m_isPromptMuon;
    unique_ptr< vector<bool> > m_isPassTrack;
    unique_ptr< vector<bool> > m_isPassPresel;
    unique_ptr< vector<bool> > m_isPassLoose;
    unique_ptr< vector<bool> > m_isPassMedium;


private:
    /* specific Tools used in this analysis */
};

#endif

