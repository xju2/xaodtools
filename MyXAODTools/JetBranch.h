#ifndef __MYXAODTOOLS_JETBRANCH_H__
#define __MYXAODTOOLS_JETBRANCH_H__

#include <vector>

#include "xAODJet/Jet.h"
#include "MyXAODTools/BranchCreatorBase.h"

using namespace std;

class JetBranch : public BranchCreatorBase 
{
public:
    JetBranch();
    virtual ~JetBranch();

    void AttachBranchToTree(TTree& );
    void ClearBranch();
    void CreateBranch();
    void Fill(const xAOD::Jet& jet);
private:
    vector<float>* emF_;
    vector<float>* hecF_;
    vector<float>* larQ_;
    vector<float>* hecQ_;
    vector<float>* sumpttrk_;
    vector<float>* frac_sampling_max_;
    vector<float>* negE_;
    vector<float>* avg_larQF_;
    vector<int>*  frac_sampling_max_index_;
};

#endif
