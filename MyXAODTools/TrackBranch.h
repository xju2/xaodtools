#ifndef __MYXAODTOOLS_TRACKBRANCH_H__
#define __MYXAODTOOLS_TRACKBRANCH_H__

#include "xAODTracking/Vertex.h"
// #include "xAODTracking/VertexContainer.h"

#include <TTree.h>

#include <vector>

#include "MyXAODTools/BranchCreatorBase.h"
using namespace std;
class TrackBranch : public BranchCreatorBase
{
public: 
    explicit TrackBranch();
    ~TrackBranch();
    void AttachBranchToTree(TTree& );
    void Fill(const xAOD::Vertex& vertex);
    void ClearBranch();
    void CreateBranch();

private:
    vector<float>* sum_pt_;
    vector<int>* n_tracks_;
    vector<float>* sum_px_;
    vector<float>* sum_py_;
    vector<float>* sum_phi_;
};

#endif
