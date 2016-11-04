/*
 *
 *
 */
#ifndef __MYXAODTOOLS_MONOJETANA_H__
#define __MYXAODTOOLS_MONOJETANA_H__

#include <vector>
#include <string>
#include <memory>

#include "MyXAODTools/AnalysisBase.h"
#include "MyXAODTools/SmearedInfo.h"

#include "AsgTools/ToolHandle.h"
#include "JetSelectorTools/JetCleaningTool.h"

#include "JetSmearing/JetMCSmearingTool.h"
#include "JetSmearing/SmearData.h"
#include "JetSmearing/PreScaleTool.h"

using namespace std;

class MonoJetAna : public AnalysisBase
{
public:
    MonoJetAna();
    virtual ~MonoJetAna();

    void AttachBranchToTree();
    void CreateBranch();
    void ClearBranch();

    int process(Long64_t ientry); // main program
private:
    int initial_tools();
    bool get_smeared_info(
            xAOD::JetContainer* jets,
            xAOD::MuonContainer* muons,
            xAOD::ElectronContainer* electrons,
            xAOD::PhotonContainer* photons,
            SmearedInfo& smeared_info);
private:
    bool m_doSmearing;
    const std::string smearJet;

    int m_nGoodJets;
    int m_nJetsBtagged;
    vector<TLorentzVector>* m_jetP4;
    float m_triggerWeight;
    int m_nBaseEl;
    int m_nBaseMu;
    int m_nBasePhoton;
    float m_minDphiJetsMET;
    float m_metEtx;
    float m_metEty;
    float m_metEt;
    float m_metSumet;
    float m_metEtSoft;
    bool m_hasCosmicMuon;
    bool m_hasBadMuon;
    // smeared data
    vector<SmearedInfo>* m_smearedData;

    // for jet-smearing
private:
    unique_ptr<JetCleaningTool> m_jetCleaningTool;
    unique_ptr<JetSmearing::JetMCSmearingTool> m_mySmearingTool;
    unique_ptr<JetSmearing::PreScaleTool> m_prescaleTool;
    TH2F* lightJetResponse;
    TH2F* bJetResponse;
    const string m_singleJetTrigger[16];
};

#endif
