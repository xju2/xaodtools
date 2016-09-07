#ifndef __MYXAODTOOLS_MUONBRANCH_H__
#define __MYXAODTOOLS_MUONBRANCH_H__

#include <vector>
#include <TLorentzVector.h>

#include "xAODMuon/Muon.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/Vertex.h"
#include "xAODTracking/VertexContainer.h"


#include "MyXAODTools/BranchCreatorBase.h"
using namespace std;

class MuonBranch : public BranchCreatorBase
{
public:
    MuonBranch();
    virtual ~MuonBranch();

    void AttachBranchToTree(TTree& );
    void ClearBranch();
    bool CreateBranch();
    void Fill(const xAOD::Muon& muon, const xAOD::EventInfo* evtInfo,
            const xAOD::Vertex* pv = NULL);

    void Fill(const xAOD::Muon& muon,
            const xAOD::EventInfo* evtInfo, const xAOD::VertexContainer* vertice);

    int matchPV(const xAOD::Muon& muon, const xAOD::VertexContainer& vertice);

public:
    static const xAOD::TrackParticle* getTrack(const xAOD::Muon& muon);

public:
    static const char* APP_NAME;

private:
    const xAOD::TrackParticle* m_track;
private:
    // used for booking trees
    int total_;
    vector<int>* author_;

    vector<float>* pt_;
    vector<float>* eta_;
    vector<float>* phi_;
    vector<float>* e_;

    vector<float>* track_charge_;
    vector<float>* track_pt_;
    vector<float>* track_eta_;
    vector<float>* track_phi_;
    vector<float>* track_e_;


    vector<float>* charge_;
    vector<int>* type_;
    vector<float>* d0_;
    vector<float>* z0_sintheta_;
    vector<float>* d0_sig_;

    vector<int>* pvID_;

public:
    // These variables have to be filled seperately
    vector<float>* eloss_;
    vector<float>* etcone30_;
    // TODO: add topoetcone30
    vector<float>* ptvarcone30_;

};

#endif
