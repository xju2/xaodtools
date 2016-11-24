#ifndef __MYXAODTOOLS_HZZ4LHELPER_H__
#define __MYXAODTOOLS_HZZ4LHELPER_H__

#include "MyXAODTools/Candidate.h"
class TTree;

#include <vector>
using namespace std;

class HZZ4lHelper{
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

public:
    explicit HZZ4lHelper();
    bool Is_close2Z(vector<Candidate*>* lep_4vec, double& m12, int& idL1, int& idL2);
    bool Is_4mu4e(vector<Candidate*>* lep_4vec );
    bool Is_2e2mu(vector<Candidate*>* ele_4vec, vector<Candidate*>* muon_4vec, int& type);

    bool MakeOutputTree(TTree& MyTree);
};

#endif
