#ifndef __MYXAODTOOLS_UPSILONANA_H__
#define __MYXAODTOOLS_UPSILONANA_H__

#include <vector>
#include <string>
#include <map>

#include "MyXAODTools/AnalysisBase.h"

using namespace std;
typedef vector<const xAOD::MuonContainer::base_value_type*> MuonVect;

class UpsilonAna : public AnalysisBase
{
public:
    UpsilonAna();
    virtual ~UpsilonAna();

    void AttachBranchToTree();
    void CreateBranch();
    void ClearBranch();

    int process(Long64_t ientry); // main program
private:
    void buildTwoMuons(const MuonVect& muons);
    void buildFourMuons(const MuonVect& muons);
    void fillOniaInfo(const xAOD::Muon& muon1, const xAOD::Muon& muon2);
    void fillQuadInfo(
            const xAOD::Muon& muon1, const xAOD::Muon& muon2,
            const xAOD::Muon& muon3, const xAOD::Muon& muon4
            );

    const xAOD::Vertex* matchFittedVertex(
            const xAOD::VertexContainer& muonVertexCont,
            const MuonVect& muons);
private:
    const xAOD::VertexContainer* m_bphy4_quad;
    const xAOD::VertexContainer* m_bphy4_pair;

private:
   // Onia information
   int m_n_onia;
   std::vector<int>* m_onia_muon1id;
   std::vector<int>* m_onia_muon2id;
   std::vector<float>* m_onia_charge;

   std::vector<float>* m_onia_pt_fitted;
   std::vector<float>* m_onia_eta_fitted;
   std::vector<float>* m_onia_phi_fitted;
   std::vector<float>* m_onia_mass_fitted;
   std::vector<float>* m_onia_x;
   std::vector<float>* m_onia_y;
   std::vector<float>* m_onia_z;
   std::vector<float>* m_onia_chi2;

   std::vector<float>* m_onia_mass;
   std::vector<float>* m_onia_pt;
   std::vector<float>* m_onia_eta;
   std::vector<float>* m_onia_phi;

   std::vector<float>* m_onia_track_mass;
   std::vector<float>* m_onia_track_pt;
   std::vector<float>* m_onia_track_eta;
   std::vector<float>* m_onia_track_phi;

   // upsilon information
   int m_n_quad;
   std::vector<float>* m_quad_charge;
   std::vector<float>* m_quad_chi2;
   std::vector<float>* m_quad_x;
   std::vector<float>* m_quad_y;
   std::vector<float>* m_quad_z;

   std::vector<int>* m_quad_nCombined;
   std::vector<int>* m_quad_id1;
   std::vector<int>* m_quad_id2;
   std::vector<int>* m_quad_id3;
   std::vector<int>* m_quad_id4;

   std::vector<float>* m_quad_mass;
   std::vector<float>* m_quad_pt;
   std::vector<float>* m_quad_eta;
   std::vector<float>* m_quad_phi;

   std::vector<float>* m_quad_track_mass;
   std::vector<float>* m_quad_track_pt;
   std::vector<float>* m_quad_track_eta;
   std::vector<float>* m_quad_track_phi;

   std::vector<float>* m_quad_fitted_mass;
   std::vector<float>* m_quad_fitted_pt;
   std::vector<float>* m_quad_fitted_eta;
   std::vector<float>* m_quad_fitted_phi;

};

#endif
