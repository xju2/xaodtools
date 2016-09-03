// dear emacs, this is -*- C++ -*-
#ifndef ANALYSISEXAMPLES_VFITZMMONAOD_H
#define ANALYSISEXAMPLES_VFITZMMONAOD_H 1

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Name    : VFitZmmOnAOD.h
// Package : offline/PhysicsAnalysis/AnalysisCommon/AnalysisExamples
// Author  : Shuwei YE
// Created : September 2009
// Modifed from ZeeZmmOnAOD in AnalysisExamples
//
// DESCRIPTION:
//
// Example of Vertex Fitting on Z->mumu with AOD
//  using TrkVKalVrtFitter
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

#include "GaudiKernel/ToolHandle.h"
//#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ITHistSvc.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "StoreGate/StoreGateSvc.h"
#include "AnalysisTools/AnalysisTools.h"
#include "AthenaBaseComps/AthAlgorithm.h"

#include "muonEvent/MuonContainer.h"
#include "muonEvent/Muon.h"

#include "TrkVKalVrtFitter/TrkVKalVrtFitter.h"
#include "VxVertex/VxCandidate.h"   

#include <string>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TLorentzVector.h"

typedef std::vector<const Analysis::MuonContainer::base_value_type*> MuonVect;

// forwards
class VFitZmmOnAOD;

// the selection function for muons
bool selectMuon(VFitZmmOnAOD *self, const MuonVect &ll);

class VFitZmmOnAOD : public AthAlgorithm {

 public:

   VFitZmmOnAOD(const std::string& name, ISvcLocator* pSvcLocator);
   ~VFitZmmOnAOD();

   StatusCode initialize();
   StatusCode execute();
   StatusCode initEvent();
   StatusCode finalize();

 private:

   /// the selection function for electrons
   friend bool selectMuon(VFitZmmOnAOD *self, const MuonVect &ll);

   /// Z->mumu reconstruction with Z as composite particles
   StatusCode zmm_on_aod();

  /// add event info to ntuple
  StatusCode addEventInfo();

  bool passMuon(const Analysis::Muon& muon, bool* blayer=NULL);
  TLorentzVector getLorentzVector(const Analysis::Muon& muon);
  TLorentzVector* getTrackLorentzV(const Analysis::Muon& muon);

  Trk::VxCandidate* VkVrtFit(const MuonVect& muons,
          HepLorentzVector* fitted_momentum=NULL
          );
  Trk::VxCandidate* VkVrtFit(
          const Analysis::Muon& muon1,
          const Analysis::Muon& muon2,
          HepLorentzVector* fitted_momentum=NULL
          );

  void buildFourMuons(const MuonVect& muons);
  void buildTwoMuons(const MuonVect& muons);
  void fillOniaInfo(const Analysis::Muon& muon1, const Analysis::Muon& muon2);
  void fillQuadInfo(
          const Analysis::Muon& muon1,
          const Analysis::Muon& muon2,
          const Analysis::Muon& muon3,
          const Analysis::Muon& muon4
          );

  bool passOniaCuts(
    const Analysis::Muon& muon1,
    const Analysis::Muon& muon2,
    const Analysis::Muon& muon3,
    const Analysis::Muon& muon4,
    float& m_light, float& m_heavy,
    float& chi2_onia1, float& chi2_onia2
    );

  int matchPV(const Analysis::Muon& muon1);
  float getChi2(Trk::VxCandidate* vx_can);

 private:

   /// get a handle to the tool helper
   ToolHandle<AnalysisTools> m_analysisTools;

   /// a handle on Store Gate 
   StoreGateSvc* m_storeGate;

   /// a handle on the Hist/TTree registration service
   ITHistSvc * m_thistSvc;

   // for handling Trk::TrkVKalVrtFitter
   ToolHandle<Trk::IVertexFitter> m_ToolIVrtFitter;
   ToolHandle<Trk::ITrkVKalVrtFitter> m_VKVrtFitter;
   // ToolHandle<Trk::TrkVKalVrtFitter> m_VKVrtFitter;

   /// the AOD muon container to retrieve
   std::string m_muonContainerName;

   // the vertex container
   std::string m_VxContainerName; 
   /// use selection cuts - for muons 
   /// to be modified thru job options
   double m_etMuonCut;
   double m_etaMuonCut;

   // variables for "tree_Zll"
   TTree* m_tree_Zll;
   double m_v0_x, m_v0_y, m_v0_z;

   // muons
   int n_muon;
   std::vector<int>* mu_author_;
   std::vector<float>* mu_pt_;
   std::vector<float>* mu_eta_;
   std::vector<float>* mu_phi_;
   std::vector<float>* mu_e_;

   std::vector<float>* mu_track_pt_;
   std::vector<float>* mu_track_eta_;
   std::vector<float>* mu_track_phi_;
   std::vector<float>* mu_track_e_;

   std::vector<float>* mu_charge_;
   std::vector<int>* mu_type_;
   std::vector<float>* mu_d0_;
   std::vector<float>* mu_d0_pv;
   std::vector<float>* mu_z0_sintheta_;
   std::vector<float>* mu_z0_;
   std::vector<float>* mu_z0_pv;
   std::vector<float>* mu_d0_sig_;
   std::vector<float>* mu_eloss_;
   std::vector<float>* mu_etcone30_;
   std::vector<float>* mu_ptcone30_;
   std::vector<int>* mu_pvID;
   std::vector<bool>* mu_blayer_;

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

   //std::vector<float>* m_onia_pt_cst;
   //std::vector<float>* m_onia_eta_cst;
   //std::vector<float>* m_onia_phi_cst;
   //std::vector<float>* m_onia_mass_cst;

   std::vector<float>* m_onia_mass;
   std::vector<float>* m_onia_pt;
   std::vector<float>* m_onia_eta;
   std::vector<float>* m_onia_phi;

   std::vector<float>* m_onia_track_mass;
   std::vector<float>* m_onia_track_pt;
   std::vector<float>* m_onia_track_eta;
   std::vector<float>* m_onia_track_phi;

   bool has_upsilon;

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

  /// variables to store Event Info stuff

  unsigned int    m_runNumber;
  unsigned int    m_eventNumber;
  unsigned int    m_eventTime;
  unsigned int    m_lumiBlock;
  unsigned int    m_bCID;
  unsigned int    m_lVL1ID;
  double  m_eventWeight;
  unsigned int    m_statusElement;
  unsigned int    m_lvl1TriggerType;
  std::vector<unsigned int>* m_lvl1TriggerInfo;
  std::vector<unsigned int>* m_lvl2TriggerInfo;
  std::vector<unsigned int>* m_evtFilterInfo;
  std::vector<std::string>* m_streamTagName;
  std::vector<std::string>* m_streamTagType;

  /// the muon histograms
  TH1F* m_aod_muon_pt;
  TH1F* m_aod_muon_eta;
  TH1F* m_aod_muon_chi2;
  TH1F* m_aod_zmm_mass_hist;
  TH1F* m_aod_muon_charge;

  TH1F* m_chi2_fit_12;
  TH1F* m_chi2_fit_34;
  TH1F* m_cutFlow;
  TH1F* m_muons_cutFlow;
};

#endif // VFITZMMONAOD_H

