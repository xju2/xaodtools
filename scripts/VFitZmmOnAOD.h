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

  bool passMuon(const Analysis::Muon& muon);
  TLorentzVector getLorentzVector(const Analysis::Muon& muon);
  float VkVrtFit(const MuonVect& muons);
  float VkVrtFit(const Analysis::Muon& muon1, const Analysis::Muon& muon2);
  void buildFourMuons(const MuonVect& muons);
  bool passOniaCuts(
          const Analysis::Muon& muon1,
          const Analysis::Muon& muon2,
          const Analysis::Muon& muon3,
          const Analysis::Muon& muon4,
          float& m_light, float& m_heavy
          );

  int matchPV(const Analysis::Muon& muon1);

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
   std::vector<float>* mu_track_pt_;
   std::vector<float>* mu_eta_;
   std::vector<float>* mu_phi_;
   std::vector<float>* mu_e_;

   std::vector<float>* mu_charge_;
   std::vector<int>* mu_type_;
   std::vector<float>* mu_d0_;
   std::vector<float>* mu_z0_sintheta_;
   std::vector<float>* mu_d0_sig_;
   std::vector<float>* mu_eloss_;
   std::vector<float>* mu_etcone30_;
   std::vector<float>* mu_ptcone30_;

   // upsilon information
   float m_upsilon_;
   float m_4l_;
   float vtx4l_chi2ndf_;
   float m34_;
   bool same_vertex_;
   float m_4l_fitted_;
   int n_combined_muons_;
   int m_index_1_;
   int m_index_2_;
   int m_index_3_;
   int m_index_4_;
   int m_pvID_1_;
   int m_pvID_2_;
   int m_pvID_3_;
   int m_pvID_4_;


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

