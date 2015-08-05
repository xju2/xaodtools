#include <stdlib.h>
#include <TVector2.h>
#include "MyXAODTools/CPToolsHelper.h"

#include "MyXAODTools/TrackBranch.h"

TrackBranch::TrackBranch(){
    CreateBranch();
}

TrackBranch::~TrackBranch()
{
    if(sum_pt_) delete sum_pt_;
    if(n_tracks_) delete n_tracks_;
    if(sum_px_) delete sum_px_;
    if(sum_py_) delete sum_py_;
    if(sum_phi_) delete sum_phi_;
}

void TrackBranch::CreateBranch(){
    sum_pt_ = new vector<float>();
    n_tracks_ = new vector<int>();
    sum_px_ = new vector<float>();
    sum_py_ = new vector<float>();
    sum_phi_ = new vector<float>();
}

void TrackBranch::AttachBranchToTree(TTree& tree)
{
    tree.Branch("track_sum_pt", &sum_pt_);
    tree.Branch("track_n", &n_tracks_);
    tree.Branch("track_sum_px", &sum_px_);
    tree.Branch("track_sum_py", &sum_py_);
    tree.Branch("track_sum_phi", &sum_phi_);
}

void TrackBranch::Fill(const xAOD::Vertex& vertex)
{
    float sum_px, sum_py;
    if(! CPToolsHelper::GetTrackSumPt(vertex, sum_px, sum_py)) {
        return;
    }
    sum_px_->push_back(-1 * sum_px);
    sum_py_->push_back(-1 * sum_py);
    sum_pt_->push_back(sqrt(sum_px*sum_px + sum_py*sum_py));
    TVector2 vec2(-1*sum_px, -1*sum_py);
    float delta_phi = TVector2::Phi_mpi_pi(vec2.Phi());
    sum_phi_->push_back(delta_phi);
    n_tracks_->push_back(vertex.nTrackParticles());
}

void TrackBranch::ClearBranch()
{
    sum_pt_->clear();
    n_tracks_->clear();
    sum_px_->clear();
    sum_py_->clear();
    sum_phi_->clear();
}
