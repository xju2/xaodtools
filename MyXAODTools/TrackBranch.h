#ifndef __MYXAODTOOLS_TRACKBRANCH_H__
#define __MYXAODTOOLS_TRACKBRANCH_H__

#include "xAODTracking/Vertex.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticle.h"
#include "TrackVertexAssociationTool/LooseTrackVertexAssociationTool.h"
#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"


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
    void AttachSumOfTracksToTree(TTree& );
    void AttachTrackToTree(TTree& );
    /* 
     * fill the sumation information of the tracks that 
     * are associated with the vertex
     */
    void Fill(const xAOD::Vertex& vertex);
    /* only save the tracks passed Loose cut! */
    void FillTrack(const xAOD::TrackParticle& track, const xAOD::VertexContainer&);
    void ClearBranch();
    void CreateBranch();

private:
    bool initialize_tools();
private:
    static const char* APP_NAME; 
    ///////////////////////////////
    // tracks associated with a Vertex
    ///////////////////////////////
    vector<float>* sum_pt_;
    vector<int>* n_tracks_;
    vector<float>* sum_px_;
    vector<float>* sum_py_;
    vector<float>* sum_phi_;

    ///////////////////////////////
    // tracks from TrackParticleContainer
    ///////////////////////////////
    vector<float>* track_pt_;
    vector<float>* track_eta_;
    vector<float>* track_phi_;
    vector<float>* track_theta_;
    vector<float>* track_d0_;
    vector<float>* track_z0_;
    vector<float>* track_beamspot_z0_;
    vector<float>* track_match_vz_;
    vector<bool>* track_pass_LooseP_;
    vector<bool>* track_pass_TightP_;

    // Utility tools
    CP::LooseTrackVertexAssociationTool *trktovxtool_;
    InDet::InDetTrackSelectionTool* trkSelTool_Loose_;
    InDet::InDetTrackSelectionTool* trkSelTool_LooseP_;
    InDet::InDetTrackSelectionTool* trkSelTool_TightP_;
};

#endif
