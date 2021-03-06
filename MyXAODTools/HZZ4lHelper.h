#ifndef __MYXAODTOOLS_HZZ4LHELPER_H__
#define __MYXAODTOOLS_HZZ4LHELPER_H__

#include "MyXAODTools/Candidate.h"
#include "xAODTruth/TruthParticleContainer.h"

class TTree;

#include <vector>
using namespace std;

class HZZ4lHelper
{
private:
    const float kZMASS;
    TTree* m_tree;
    float m_m4l;
    float m_Hpt_;
    float m_Hphi_;
    float m_mZ1;
    float m_mZ2;
    float m_Z1_lepplus_pt;
    float m_Z1_lepminus_pt;
    float m_Z2_lepplus_pt;
    float m_Z2_lepminus_pt;
    float m_Z1_lepplus_eta;
    float m_Z1_lepminus_eta;
    float m_Z2_lepplus_eta;
    float m_Z2_lepminus_eta;

    float truth_h_mass;
    float truth_z1_mass;
    float truth_z2_mass;
    float truth_h_pt;
    float truth_z1_pt;
    float truth_z2_pt;

    float truth_l1_pt;
    float truth_l2_pt;
    float truth_l3_pt;
    float truth_l4_pt;
    float truth_l1_eta;
    float truth_l2_eta;
    float truth_l3_eta;
    float truth_l4_eta;

    int event_type;
    int event_type_truth;
private:
    static bool descend_sort(int i,int j) { return (i>j); }

    bool passElectronFid(float pt, float eta);
    bool passMuonFid(float pt, float eta);
    bool passLepPt(vector<float>& lep_pt);
    bool passZMass(float z1, float z2, float m4l);
    bool passKinemitic(vector<float>& lep_pt, vector<float>& lep_eta, int event_type);
public:
    explicit HZZ4lHelper();
    ~HZZ4lHelper();
    bool Is_close2Z(vector<Candidate*>* lep_4vec, double& m12, int& idL1, int& idL2);
    bool Is_4mu4e(vector<Candidate*>* lep_4vec );
    bool Is_2e2mu(vector<Candidate*>* ele_4vec, vector<Candidate*>* muon_4vec, int& type);

    bool MakeOutputTree(TTree& MyTree);
    bool MakeTruthTree(TTree& MyTree);
    bool GetTruthInfo(const xAOD::TruthParticleContainer& particles);

    void Clear();
    int passFiducial();
    int passFiducialTruth();
    void setType(int type){ event_type = type; }
    int getTruthType(){return event_type_truth;}
};

#endif
