#ifndef __MYXAODTOOLS_MUONBRANCH_H__
#define __MYXAODTOOLS_MUONBRANCH_H__

#include <vector>
#include <memory>
#include <TLorentzVector.h>

#include "xAODMuon/Muon.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/Vertex.h"
#include "xAODTracking/VertexContainer.h"

#include "MyXAODTools/BranchCreatorBase.h"

using namespace std;

namespace CP{
class MuonSelectionTool;
}

class MuonBranch: public BranchCreatorBase
{
public:
    MuonBranch();
    virtual ~MuonBranch();

    void AttachBranchToTree(TTree& );
    void ClearBranch();
    void Fill(const xAOD::Muon& muon, const xAOD::EventInfo* evtInfo,
            const xAOD::Vertex* pv = NULL);

    void Fill(const xAOD::Muon& muon,
            const xAOD::EventInfo* evtInfo, const xAOD::VertexContainer* vertice);

    int matchPV(const xAOD::Muon& muon, const xAOD::VertexContainer& vertice);
private:
    bool CreateBranch();
    int initial_tools();
    void addDetailedInfo(const xAOD::Muon& muon);
    void getPrecisionLayer(const xAOD::Muon& muon, uint8_t& precLayer, uint8_t& precHoleLayer) const;
    float getQoverPsig(const xAOD::Muon& muon) const;

public:
    static const xAOD::TrackParticle* getTrack(const xAOD::Muon& muon);

public:
    static const char* APP_NAME;
    unique_ptr<CP::MuonSelectionTool> m_muonSelectionTool;

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
    vector<int>* quality_;

public:
    // These variables have to be filled seperately
    vector<float>* eloss_;
    vector<float>* etcone30_;
    // TODO: add topoetcone30
    vector<float>* ptvarcone30_;

    // variables used in MuonSelectorTools
    unique_ptr< vector<float> > m_qOverpSignif;
    unique_ptr< vector<uint8_t> > m_nPrecLayer;
    unique_ptr< vector<uint8_t> > m_nPrecHoleLayer;
};

#endif
