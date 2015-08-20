#include <stdlib.h>
#include <exception>

#include <TVector2.h>
#include "CPAnalysisExamples/errorcheck.h"

#include "MyXAODTools/CPToolsHelper.h"
#include "MyXAODTools/TrackBranch.h"

const char* TrackBranch::APP_NAME = "TrackBranch";

TrackBranch::TrackBranch(){
    CreateBranch();
    if(!initialize_tools())
    {
        throw runtime_error("Cannot initialize TrackBranch");
    }
}


TrackBranch::~TrackBranch()
{
    if(sum_pt_) delete sum_pt_;
    if(n_tracks_) delete n_tracks_;
    if(sum_px_) delete sum_px_;
    if(sum_py_) delete sum_py_;
    if(sum_phi_) delete sum_phi_;

    if(track_pt_) delete track_pt_;
    if(track_eta_) delete track_eta_;
    if(track_phi_) delete track_phi_;
    if(track_theta_) delete track_theta_;
    if(track_d0_) delete track_d0_;
    if(track_z0_) delete track_z0_;
    if(track_beamspot_z0_) delete track_beamspot_z0_;
    if(track_match_vz_) delete track_match_vz_;
    if(track_pass_LooseP_) delete track_pass_LooseP_;
    if(track_pass_TightP_) delete track_pass_TightP_;
    
    if(trktovxtool_) delete trktovxtool_;
    if(trkSelTool_Loose_) delete trkSelTool_Loose_;
    if(trkSelTool_LooseP_) delete trkSelTool_LooseP_;
    if(trkSelTool_TightP_) delete trkSelTool_TightP_;
}

void TrackBranch::CreateBranch(){
    sum_pt_ = new vector<float>();
    n_tracks_ = new vector<int>();
    sum_px_ = new vector<float>();
    sum_py_ = new vector<float>();
    sum_phi_ = new vector<float>();

    track_pt_ = new vector<float>();
    track_eta_ = new vector<float>();
    track_phi_ = new vector<float>();
    track_theta_ = new vector<float>();
    track_d0_ = new vector<float>();
    track_z0_ = new vector<float>();
    track_beamspot_z0_ = new vector<float>();
    track_match_vz_ = new vector<float>();
    track_pass_LooseP_ = new vector<bool>();
    track_pass_TightP_ = new vector<bool>();
}

bool TrackBranch::initialize_tools()
{
    trktovxtool_ = new CP::LooseTrackVertexAssociationTool("LooseTrackVertexAssociationTool");
    trktovxtool_->msg().setLevel( MSG::ERROR);
    CHECK( trktovxtool_->setProperty("dzSinTheta_cut", 100.) );

    trkSelTool_Loose_ = new InDet::InDetTrackSelectionTool( "TrackSelectionLoose" );
    CHECK( trkSelTool_Loose_->setProperty("CutLevel", "Loose") );
    CHECK( trkSelTool_Loose_->setProperty("maxZ0SinTheta", 3.) );
    CHECK( trkSelTool_Loose_->setProperty("maxD0overSigmaD0", 1.5) );
    CHECK( trkSelTool_Loose_->initialize() );

    trkSelTool_LooseP_ = new InDet::InDetTrackSelectionTool( "TrackSelectionLooseP" );
    CHECK( trkSelTool_LooseP_->setProperty("CutLevel", "LoosePrimary") );
    CHECK( trkSelTool_LooseP_->setProperty("maxZ0SinTheta", 3.) );
    CHECK( trkSelTool_LooseP_->setProperty("maxD0overSigmaD0", 1.5) );
    CHECK( trkSelTool_LooseP_->initialize() );

    trkSelTool_TightP_ = new InDet::InDetTrackSelectionTool( "TrackSelectionTightP" );
    CHECK( trkSelTool_TightP_->setProperty("CutLevel", "TightPrimary") );
    CHECK( trkSelTool_TightP_->setProperty("maxZ0SinTheta", 3.) );
    CHECK( trkSelTool_TightP_->setProperty("maxD0overSigmaD0", 1.5) );
    CHECK( trkSelTool_TightP_->initialize() ); 
    return true;
}

void TrackBranch::AttachBranchToTree(TTree& tree)
{
    this->AttachSumOfTracksToTree(tree);
    this->AttachTrackToTree(tree);
}

void TrackBranch::AttachSumOfTracksToTree(TTree& tree){
    tree.Branch("track_sum_pt", &sum_pt_);
    tree.Branch("track_n", &n_tracks_);
    tree.Branch("track_sum_px", &sum_px_);
    tree.Branch("track_sum_py", &sum_py_);
    tree.Branch("track_sum_phi", &sum_phi_);
}

void TrackBranch::AttachTrackToTree(TTree& tree){
    tree.Branch("track_pt", &track_pt_);
    tree.Branch("track_eta", &track_eta_);
    tree.Branch("track_phi", &track_phi_);
    tree.Branch("track_theta", &track_theta_);
    tree.Branch("track_d0", &track_d0_);
    tree.Branch("track_z0", &track_z0_);
    tree.Branch("track_beamspot_z0", &track_beamspot_z0_);
    tree.Branch("track_match_vz", &track_match_vz_);
    tree.Branch("track_LooseP", &track_pass_LooseP_);
    tree.Branch("track_TightP", &track_pass_TightP_);
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

    track_pt_->clear();
    track_eta_->clear();
    track_phi_->clear();
    track_theta_->clear();
    track_d0_->clear();
    track_z0_->clear();
    track_beamspot_z0_->clear();
    track_match_vz_->clear();
    track_pass_LooseP_->clear();
    track_pass_TightP_->clear();
}

void TrackBranch::FillTrack(const xAOD::TrackParticle& track, const xAOD::VertexContainer& vxCont)
{
    ElementLink< xAOD::VertexContainer> match_vx = trktovxtool_->getUniqueMatchVertexLink(track,vxCont);
    if ( match_vx.isValid() )
    {
        if ( !trkSelTool_Loose_->accept(track, *match_vx) ) return;
        track_match_vz_->push_back( (*match_vx)->z() );
        track_pass_LooseP_->push_back( trkSelTool_LooseP_->accept(track, *match_vx) );
        track_pass_TightP_->push_back( trkSelTool_TightP_->accept(track, *match_vx) );
    } else {
        if ( !trkSelTool_Loose_->accept(track) ) return;
        track_match_vz_->push_back( 0. );
        track_pass_LooseP_->push_back( trkSelTool_LooseP_->accept(track) );
        track_pass_TightP_->push_back( trkSelTool_TightP_->accept(track) );
    }
    track_pt_->push_back( track.pt() );
    track_eta_->push_back( track.eta() );
    track_phi_->push_back( track.phi() );
    track_theta_->push_back( track.theta() );
    track_d0_->push_back( track.d0() );
    track_z0_->push_back( track.z0() );
    track_beamspot_z0_->push_back( track.vz() );
}
