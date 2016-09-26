#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IToolSvc.h"


#include "TH1.h"
#include "TH2.h"

// The Muon
#include "muonEvent/Muon.h"

// #include "TrkParticleBase/TrackParticleBase.h"

// PDT particle property
#include "GaudiKernel/IPartPropSvc.h"
#include "HepPDT/ParticleDataTable.hh"
#include "HepPDT/ParticleData.hh"

// Constituent navigation
#include "Navigation/NavigationToken.h"

// common implementation of all particles
#include "ParticleEvent/ParticleBaseContainer.h"

// vertex
#include "VxVertex/VxContainer.h"   
#include "VxVertex/VxCandidate.h"   
#include "VxVertex/RecVertex.h"
#include "VxVertex/VxTrackAtVertex.h"

// the composite particle
#include "ParticleEvent/CompositeParticle.h"

// particle jets
#include "JetEvent/JetCollection.h"

// analysis tools
#include "AnalysisUtils/AnalysisCombination.h"

// Event Info
#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"
#include "EventInfo/EventType.h"
#include "EventInfo/TriggerInfo.h"

// Electron Info
#include "egammaEvent/egammaPIDdefs.h"

// the header file
#include "VFitZmmOnAOD.h"

#include <stdint.h>
#include <algorithm>
#include <math.h>
#include <functional>

#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#endif

// static const double mMuon = 105.658369;
static double mMuon = 105.66; // MeV
const float ONIA_ONE_PT_CUT = 4000;
const float ONIA_TWO_PT_CUT = 3000;
const float UPSILON_MASS = 9460.30;

const Trk::RecVertex* primVtx = NULL;
const VxContainer* vxCont = NULL;

using namespace Analysis;

//////////////////////////////////////////////////////////////////////////////////////
// Constructor

VFitZmmOnAOD::VFitZmmOnAOD(const std::string& name,ISvcLocator* pSvcLocator): 
  AthAlgorithm(name, pSvcLocator),
  m_analysisTools( "AnalysisTools", this ),
  m_VKVrtFitter("Trk::TrkVKalVrtFitter", this)
{
  // switches to control the analysis through job options :: these are the default
  // to changed in the job options

  // The Electron AOD container name & selection cuts
  declareProperty( "AnalysisTools", m_analysisTools );

  // The Muon AOD container name & selection cuts
  declareProperty("MuonContainer", m_muonContainerName = "MuonCollection");
  declareProperty("MuonEtCut", m_etMuonCut = 6.0*GeV);
  declareProperty("MuonEtaCut", m_etaMuonCut = 2.5);

  // Name of the Vertex container
  declareProperty("VxContainer", m_VxContainerName= "VxPrimaryCandidate",
                                 "Name of the Vertex container");

  // vertex fitter
  declareProperty("TrkVKalVrtFitterTool", m_ToolIVrtFitter);
  declareProperty("TrackToVertexTool", m_trackToVertexTool);

  // if do mass constraint
  declareProperty("doMassConstraint", m_doMassConstraint = false, "Perform mass constraint");
  declareProperty("doPVConstraint", m_doPVConstraint = false, "Perform PV constraint");
  }


/////////////////////////////////////////////////////////////////////////////////////
// Destructor - check up memory allocation
// delete any memory allocation on the heap

VFitZmmOnAOD::~VFitZmmOnAOD() {}

////////////////////////////////////////////////////////////////////////////////////
// Initialize
// get a handle on the analysis tools
// book histograms

StatusCode VFitZmmOnAOD::initialize() {

  ATH_MSG_INFO("Initializing VFitZmmOnAOD");

  // get tool svc
  //
  IToolSvc* toolSvc;
  StatusCode sc = service("ToolSvc",toolSvc); 
  if (StatusCode::SUCCESS != sc) {
    ATH_MSG_ERROR("Unable to retrieve ToolSvc");
    return StatusCode::FAILURE;
  }
  

//  VKalVrt vertex fitter
    if (m_VKVrtFitter.retrieve().isFailure()) {
      ATH_MSG_WARNING("Can't find Trk::TrkVKalVrtFitter");
      return StatusCode::SUCCESS; 
    } else {
      ATH_MSG_INFO("Trk::TrkVKalVrtFitter found");
    }

    // Retrieve TrackToVertex extrapolator tool
    if (m_trackToVertexTool.retrieve().isFailure()) {
        ATH_MSG_WARNING("Cannot find TrackToVertex tool");
        return StatusCode::SUCCESS;
    }else {
        ATH_MSG_INFO("TrackToVertex Tool");
    }

  // get Particle Property service
  IPartPropSvc* p_PartPropSvc = 0;
  sc = service("PartPropSvc", p_PartPropSvc, true);
  if ( sc.isFailure() || 0 == p_PartPropSvc) {
    ATH_MSG_ERROR("Unable to initialize Particle Properties Service");
    return StatusCode::FAILURE;;
  }
  const HepPDT::ParticleDataTable* particleTable = p_PartPropSvc->PDT();
  const HepPDT::ParticleData* muonPDT = particleTable->particle(13);
  ATH_MSG_DEBUG(":muonPDT: mass=" << muonPDT->mass() << ", id=" << muonPDT->pid() << ", name=" << muonPDT->name() );
  mMuon = muonPDT->mass();

  // get a handle on the analysis tools
  sc = m_analysisTools.retrieve();
  if ( sc.isFailure() ) {
    ATH_MSG_ERROR("Can't get handle on analysis tools");
    return sc;
  }

  // Return a pointer to THistSvc
  sc = service("THistSvc", m_thistSvc);
  if(sc.isFailure() ){
    ATH_MSG_ERROR("Unable to retrieve pointer to THistSvc");
    return sc;
  }

  // create vector for good-leptons
  m_good_electrons = new ElecVect();
  m_good_muons = new MuonVect();

  // muons
  mu_author_ = new std::vector<int>;
  mu_pt_ = new std::vector<float>;
  mu_eta_ = new std::vector<float>;
  mu_phi_ = new std::vector<float>;
  mu_e_ = new std::vector<float>;

  mu_track_pt_ = new std::vector<float>();
  mu_track_eta_ = new std::vector<float>();
  mu_track_phi_ = new std::vector<float>();
  mu_track_e_ = new std::vector<float>();

  mu_charge_ = new std::vector<float>;
  mu_type_ = new std::vector<int>;
  mu_d0_ = new std::vector<float>;
  mu_d0_pv = new std::vector<float>;
  mu_z0_sintheta_ = new std::vector<float>;
  mu_z0_pv_sintheta_ = new std::vector<float>;
  mu_z0_ = new std::vector<float>;
  mu_z0_pv = new std::vector<float>;
  mu_d0_sig_ = new std::vector<float>;
  mu_d0_pv_sig_ = new std::vector<float>;
  mu_eloss_ = new std::vector<float>;
  mu_etcone30_ = new std::vector<float>;
  mu_ptcone30_ = new std::vector<float>;
  mu_pvID = new std::vector<int>;
  mu_blayer_ = new std::vector<bool>;

  // Onia information
  m_onia_type = new std::vector<int>;
  m_onia_muon1id = new std::vector<int>;
  m_onia_muon2id = new std::vector<int>;
  m_onia_charge = new std::vector<float>;

  m_onia_pt_fitted = new std::vector<float>;
  m_onia_eta_fitted = new std::vector<float>;
  m_onia_phi_fitted = new std::vector<float>;
  m_onia_mass_fitted = new std::vector<float>;
  m_onia_x = new std::vector<float>;
  m_onia_y = new std::vector<float>;
  m_onia_z = new std::vector<float>;
  m_onia_chi2 = new std::vector<float>;

  m_onia_mass   = new std::vector<float>;
  m_onia_pt     = new std::vector<float>;
  m_onia_eta    = new std::vector<float>;
  m_onia_phi    = new std::vector<float>;

  m_onia_track_mass   = new std::vector<float>;
  m_onia_track_pt     = new std::vector<float>;
  m_onia_track_eta    = new std::vector<float>;
  m_onia_track_phi    = new std::vector<float>;

  // has_upsilon = false;

  // upsilon information
  m_quad_type = new std::vector<int>;
  m_quad_charge = new std::vector<float>;
  m_quad_chi2 = new std::vector<float>;
  m_quad_x      = new std::vector<float>;
  m_quad_y      = new std::vector<float>;
  m_quad_z      = new std::vector<float>;

  m_quad_nCombined = new std::vector<int>;
  m_quad_id1 = new std::vector<int>;
  m_quad_id2 = new std::vector<int>;
  m_quad_id3 = new std::vector<int>;
  m_quad_id4 = new std::vector<int>;

  m_quad_mass   = new std::vector<float>;
  m_quad_pt     = new std::vector<float>;
  m_quad_eta    = new std::vector<float>;
  m_quad_phi    = new std::vector<float>;

  m_quad_track_mass   = new std::vector<float>;
  m_quad_track_pt     = new std::vector<float>;
  m_quad_track_eta    = new std::vector<float>;
  m_quad_track_phi    = new std::vector<float>;

  m_quad_fitted_mass   = new std::vector<float>;
  m_quad_fitted_pt     = new std::vector<float>;
  m_quad_fitted_eta    = new std::vector<float>;
  m_quad_fitted_phi    = new std::vector<float>;

  //
  // event info variables
  m_lvl1TriggerInfo = new std::vector<unsigned int>; 
  m_lvl2TriggerInfo = new std::vector<unsigned int>; 
  m_evtFilterInfo   = new std::vector<unsigned int>;   
  m_streamTagName   = new std::vector<std::string>;   
  m_streamTagType   = new std::vector<std::string>;   

  // the TTree
  m_tree_Zll = new TTree("physics","physics");
  sc = m_thistSvc->regTree("/AANT/physics", m_tree_Zll);

  // first add Event info stuff
  m_tree_Zll->Branch("Run",  &m_runNumber,   "Run/I");    // run number
  m_tree_Zll->Branch("Event",&m_eventNumber, "Event/I");  // event number
  m_tree_Zll->Branch("Time", &m_eventTime,   "Time/I");   // time stamp
  m_tree_Zll->Branch("LumiBlock", &m_lumiBlock,"LumiBlock/I"); // lum block num 
  m_tree_Zll->Branch("BCID", &m_bCID,"BCID/I"); // bunch crossing ID
  m_tree_Zll->Branch("LVL1ID", &m_lVL1ID,"LVL1ID/I"); // trigger LVL1 ID
  m_tree_Zll->Branch("Weight", &m_eventWeight, "Weight/D"); // weight
  m_tree_Zll->Branch("StatusElement",  &m_statusElement, "StatusElement/I");
  m_tree_Zll->Branch("LVL1TriggerType",  &m_lvl1TriggerType, "LVL1TriggerType/I");
  m_tree_Zll->Branch("LVL1TriggerInfo",&m_lvl1TriggerInfo);
  m_tree_Zll->Branch("LVL2TriggerInfo",&m_lvl2TriggerInfo);
  m_tree_Zll->Branch("EventFilterInfo",&m_evtFilterInfo);
  m_tree_Zll->Branch("StreamTagName",&m_streamTagName);
  m_tree_Zll->Branch("StreamTagType",&m_streamTagType);
  
  // now variables from algorithm
  m_v0_x = new std::vector<float>;
  m_v0_y = new std::vector<float>;
  m_v0_z = new std::vector<float>;
  m_tree_Zll->Branch("v0_x", &m_v0_x);
  m_tree_Zll->Branch("v0_y", &m_v0_y);
  m_tree_Zll->Branch("v0_z", &m_v0_z);

  // muon variables
  m_tree_Zll->Branch("n_muon", &n_muon, "n_muon/I");
  m_tree_Zll->Branch("mu_author", &mu_author_);
  m_tree_Zll->Branch("mu_pt", &mu_pt_);
  m_tree_Zll->Branch("mu_eta", &mu_eta_);
  m_tree_Zll->Branch("mu_phi", &mu_phi_);
  m_tree_Zll->Branch("mu_e", &mu_e_);

  m_tree_Zll->Branch("mu_track_pt", &mu_track_pt_);
  m_tree_Zll->Branch("mu_track_eta", &mu_track_eta_);
  m_tree_Zll->Branch("mu_track_phi", &mu_track_phi_);
  m_tree_Zll->Branch("mu_track_e", &mu_track_e_);

  m_tree_Zll->Branch("mu_charge", &mu_charge_);
  m_tree_Zll->Branch("mu_type", &mu_type_);
  m_tree_Zll->Branch("mu_d0", &mu_d0_);
  m_tree_Zll->Branch("mu_d0_pv", &mu_d0_pv);
  m_tree_Zll->Branch("mu_z0_sintheta", &mu_z0_sintheta_);
  m_tree_Zll->Branch("mu_z0_pv_sintheta", &mu_z0_pv_sintheta_);
  m_tree_Zll->Branch("mu_z0", &mu_z0_);
  m_tree_Zll->Branch("mu_z0_pv", &mu_z0_pv);
  m_tree_Zll->Branch("mu_d0_sig", &mu_d0_sig_);
  m_tree_Zll->Branch("mu_d0_pv_sig", &mu_d0_pv_sig_);
  m_tree_Zll->Branch("mu_eloss", &mu_eloss_);
  m_tree_Zll->Branch("mu_etcone30", &mu_etcone30_);
  m_tree_Zll->Branch("mu_ptvarcone30", &mu_ptcone30_);
  m_tree_Zll->Branch("mu_pvID", &mu_pvID);
  m_tree_Zll->Branch("mu_blayer", &mu_blayer_);

  // Onia information
  m_tree_Zll->Branch("n_onia", &m_n_onia, "n_onia/I");
  m_tree_Zll->Branch("onia_type", &m_onia_type);
  m_tree_Zll->Branch("onia_id1", &m_onia_muon1id);
  m_tree_Zll->Branch("onia_id2", &m_onia_muon2id);
  m_tree_Zll->Branch("onia_charge", &m_onia_charge);

  m_tree_Zll->Branch("onia_fitted_pt", &m_onia_pt_fitted);
  m_tree_Zll->Branch("onia_fitted_eta", &m_onia_eta_fitted);
  m_tree_Zll->Branch("onia_fitted_phi", &m_onia_phi_fitted);
  m_tree_Zll->Branch("onia_fitted_mass", &m_onia_mass_fitted);

  m_tree_Zll->Branch("onia_x", &m_onia_x);
  m_tree_Zll->Branch("onia_y", &m_onia_y);
  m_tree_Zll->Branch("onia_z", &m_onia_z);
  m_tree_Zll->Branch("onia_chi2", &m_onia_chi2);

  m_tree_Zll->Branch("onia_mass", &m_onia_mass);
  m_tree_Zll->Branch("onia_pt", &m_onia_pt);
  m_tree_Zll->Branch("onia_eta", &m_onia_eta);
  m_tree_Zll->Branch("onia_phi", &m_onia_phi);

  m_tree_Zll->Branch("onia_track_mass", &m_onia_track_mass);
  m_tree_Zll->Branch("onia_track_pt", &m_onia_track_pt);
  m_tree_Zll->Branch("onia_track_eta", &m_onia_track_eta);
  m_tree_Zll->Branch("onia_track_phi", &m_onia_track_phi);

  //m_tree_Zll->Branch("has_upsilon", &has_upsilon, "has_upsilon/O");

  // upsilon 
  m_tree_Zll->Branch("n_quad", &m_n_quad, "n_quad/I");
  m_tree_Zll->Branch("quad_type", &m_quad_type);
  m_tree_Zll->Branch("quad_charge", &m_quad_charge);
  m_tree_Zll->Branch("quad_chi2", &m_quad_chi2);
  m_tree_Zll->Branch("quad_x", &m_quad_x);
  m_tree_Zll->Branch("quad_y", &m_quad_y);
  m_tree_Zll->Branch("quad_z", &m_quad_z);
  m_tree_Zll->Branch("quad_nCombined", &m_quad_nCombined);
  m_tree_Zll->Branch("quad_id1", &m_quad_id1);
  m_tree_Zll->Branch("quad_id2", &m_quad_id2);
  m_tree_Zll->Branch("quad_id3", &m_quad_id3);
  m_tree_Zll->Branch("quad_id4", &m_quad_id4);

  m_tree_Zll->Branch("quad_mass", &m_quad_mass);
  m_tree_Zll->Branch("quad_pt", &m_quad_pt);
  m_tree_Zll->Branch("quad_eta", &m_quad_eta);
  m_tree_Zll->Branch("quad_phi", &m_quad_phi);
  m_tree_Zll->Branch("quad_track_mass", &m_quad_track_mass);
  m_tree_Zll->Branch("quad_track_pt", &m_quad_track_pt);
  m_tree_Zll->Branch("quad_track_eta", &m_quad_track_eta);
  m_tree_Zll->Branch("quad_track_phi", &m_quad_track_phi);
  m_tree_Zll->Branch("quad_fitted_mass", &m_quad_fitted_mass);
  m_tree_Zll->Branch("quad_fitted_pt", &m_quad_fitted_pt);
  m_tree_Zll->Branch("quad_fitted_eta", &m_quad_fitted_eta);
  m_tree_Zll->Branch("quad_fitted_phi", &m_quad_fitted_phi);

  // Add information for electrons
  m_elec_eta    = new std::vector<float>;
  m_elec_phi    = new std::vector<float>;
  m_elec_pt     = new std::vector<float>;
  m_elec_auth   = new std::vector<int>;
  m_elec_charge = new std::vector<float>;
  m_elec_id     = new std::vector<bool>;
  m_elec_ISEM   = new std::vector<double>;
 
  m_elec_etCl   = new std::vector<float>; 
  m_elec_etas2  = new std::vector<float>; 
  m_elec_rhad   = new std::vector<double>; 
  m_elec_rhad1  = new std::vector<double>; 
  m_elec_reta   = new std::vector<float>; 
  m_elec_w2     = new std::vector<float>; 

  m_elec_f1     = new std::vector<double>; 
  m_elec_f3     = new std::vector<double>; 
  m_elec_wstot  = new std::vector<double>; 
  m_elec_DEmaxs1= new std::vector<double>; 
  m_elec_deltaEta = new std::vector<double>; 
  m_elec_deltaPhiRescaled   = new std::vector<double>; 

  m_elec_nSi    = new std::vector<int>; 
  m_elec_nSiDeadSensors = new std::vector<int>; 
  m_elec_nPix   = new std::vector<int>; 
  m_elec_nPixDeadSensors    = new std::vector<int>; 
  m_elec_nTRThigh   = new std::vector<int>; 
  m_elec_nTRThighOutliers   = new std::vector<int>; 
  m_elec_nTRT   = new std::vector<int>; 
  m_elec_nTRTOutliers   = new std::vector<int>; 

  m_elec_dpOverp    = new std::vector<double>;
  m_elec_rTRT       = new std::vector<double>;
  m_elec_expBLayer  = new std::vector<bool>;
  m_elec_trackd0    = new std::vector<double>;

  m_tree_Zll->Branch("el_eta", &m_elec_eta);
  m_tree_Zll->Branch("el_phi", &m_elec_phi);
  m_tree_Zll->Branch("el_pt", &m_elec_pt);
  m_tree_Zll->Branch("el_auth", &m_elec_auth);
  m_tree_Zll->Branch("el_charge", &m_elec_charge);
  m_tree_Zll->Branch("el_id", &m_elec_id);
  m_tree_Zll->Branch("el_isEM", &m_elec_ISEM);

  m_tree_Zll->Branch("el_etCluster", &m_elec_etCl);
  m_tree_Zll->Branch("el_etas2", &m_elec_etas2);
  m_tree_Zll->Branch("el_rhad", &m_elec_rhad);
  m_tree_Zll->Branch("el_rhad1", &m_elec_rhad1);
  m_tree_Zll->Branch("el_reta", &m_elec_reta);
  m_tree_Zll->Branch("el_w2", &m_elec_w2);

  m_tree_Zll->Branch("el_f1", &m_elec_f1);
  m_tree_Zll->Branch("el_f3", &m_elec_f3);
  m_tree_Zll->Branch("el_wstot", &m_elec_wstot);
  m_tree_Zll->Branch("el_DEmaxs1", &m_elec_DEmaxs1);
  m_tree_Zll->Branch("el_deltaEta", &m_elec_deltaEta);
  m_tree_Zll->Branch("el_deltaPhi", &m_elec_deltaPhiRescaled);

  m_tree_Zll->Branch("el_nSi", &m_elec_nSi);
  m_tree_Zll->Branch("el_nSiDeadSensors", &m_elec_nSiDeadSensors);
  m_tree_Zll->Branch("el_nPix", &m_elec_nPix);
  m_tree_Zll->Branch("el_nPixDeadSensors", &m_elec_nPixDeadSensors);
  m_tree_Zll->Branch("el_nTRThigh", &m_elec_nTRThigh);
  m_tree_Zll->Branch("el_nTRThighOutliers", &m_elec_nTRThighOutliers);
  m_tree_Zll->Branch("el_nTRT", &m_elec_nTRT);
  m_tree_Zll->Branch("el_nTRTOutliers", &m_elec_nTRTOutliers);

  m_tree_Zll->Branch("el_EoverP", &m_elec_dpOverp);
  m_tree_Zll->Branch("el_rTRT", &m_elec_rTRT);
  m_tree_Zll->Branch("el_expectedBLayer", &m_elec_expBLayer);
  m_tree_Zll->Branch("el_trackd0", &m_elec_trackd0);

  // the histograms
  m_cutFlow = new TH1F("cutFlow", "cut flow", 20, 0.5, 20.5);
  sc = m_thistSvc->regHist("/AANT/Muon/cutFlow", m_cutFlow);

  m_muons_cutFlow = new TH1F("muons_cutFlow", "cut flow", 20, 0.5, 20.5);
  sc = m_thistSvc->regHist("/AANT/Muon/muons_cutFlow", m_muons_cutFlow);
  m_muons_cutFlow->GetXaxis()->SetBinLabel(1, "all");
  m_muons_cutFlow->GetXaxis()->SetBinLabel(2, "kinematic");
  m_muons_cutFlow->GetXaxis()->SetBinLabel(3, "ID Hits");
  m_muons_cutFlow->GetXaxis()->SetBinLabel(4, "muonType");
  m_muons_cutFlow->GetXaxis()->SetBinLabel(5, "chargeConsitent");

  // Muon histogram booking 
  m_aod_muon_pt        = new TH1F("aod_muon_pt","aod pt mu",50,0,250.*GeV);
  sc = m_thistSvc->regHist("/AANT/Muon/aod_muon_pt", m_aod_muon_pt);
  m_aod_muon_eta       = new TH1F("aod_muon_eta","aod eta mu",70,-3.5,3.5);
  sc = m_thistSvc->regHist("/AANT/Muon/aod_muon_eta", m_aod_muon_eta);
  m_aod_muon_chi2      = new TH1F("aod_muon_chi2","aod chi2 mu",200,0.0,1000.0);
  sc = m_thistSvc->regHist("/AANT/Muon/aod_muon_chi2", m_aod_muon_chi2);
  m_aod_zmm_mass_hist  = new TH1F("Mmm","Mmm",50,0,250.*GeV);
  sc = m_thistSvc->regHist("/AANT/Muon/Mmm", m_aod_zmm_mass_hist);
  m_aod_muon_charge = new TH1F("Muon_charge","Muon_charge",10,-5,5);
  sc = m_thistSvc->regHist("/AANT/Muon/Muon_charge", m_aod_muon_charge);

  m_chi2_fit_12 = new TH1F("chi2_fit_12", "chi2 fit for 1 and 2",
          100, -3, 17);
  sc = m_thistSvc->regHist("/AANT/Muon/chi2_fit_12", m_chi2_fit_12);
  m_chi2_fit_34 = new TH1F("chi2_fit_34", "chi2 fit for 3 and 4",
          100, -3, 17);
  sc = m_thistSvc->regHist("/AANT/Muon/chi2_fit_34", m_chi2_fit_34);


  if( sc.isFailure() ){
    ATH_MSG_ERROR("ROOT Hist aod_zmm_mass_hist registration failed");
    return sc;
  }
  return StatusCode::SUCCESS;
}		 

///////////////////////////////////////////////////////////////////////////////////
// Finalize - delete any memory allocation from the heap

StatusCode VFitZmmOnAOD::finalize() {
  return StatusCode::SUCCESS;

}

///////////////////////////////////////////////////////////////////////////////////
// initialization before processing a new event
StatusCode VFitZmmOnAOD::initEvent() {

  ATH_MSG_DEBUG("initEvent()");

  m_good_electrons->clear();
  m_good_muons->clear();

  m_v0_x ->clear();
  m_v0_y ->clear();
  m_v0_z ->clear();
  
  n_muon = 0;
  mu_author_->clear();
  mu_pt_->clear();
  mu_eta_->clear();
  mu_phi_->clear();
  mu_e_->clear();

  mu_track_pt_->clear();
  mu_track_eta_->clear();
  mu_track_phi_->clear();
  mu_track_e_->clear();

  mu_charge_->clear();
  mu_type_->clear();
  mu_d0_->clear();
  mu_d0_pv->clear();
  mu_z0_sintheta_->clear();
  mu_z0_pv_sintheta_->clear();
  mu_z0_->clear();
  mu_z0_pv->clear();
  mu_d0_sig_->clear();
  mu_d0_pv_sig_->clear();
  mu_eloss_->clear();
  mu_etcone30_->clear();
  mu_ptcone30_->clear();
  mu_pvID->clear();
  mu_blayer_->clear();

  // Onia information
  m_n_onia = 0;
  m_onia_type->clear();
  m_onia_muon1id->clear();
  m_onia_muon2id->clear();
  m_onia_charge ->clear();

  m_onia_pt_fitted->clear();
  m_onia_eta_fitted->clear();
  m_onia_phi_fitted->clear();
  m_onia_mass_fitted->clear();
  m_onia_x->clear();
  m_onia_y->clear();
  m_onia_z->clear();
  m_onia_chi2->clear();

  m_onia_mass->clear();
  m_onia_pt->clear();
  m_onia_eta->clear();
  m_onia_phi->clear();
  m_onia_track_mass->clear();
  m_onia_track_pt->clear();
  m_onia_track_eta->clear();
  m_onia_track_phi->clear();

  // has_upsilon = false;

  // 
  m_n_quad = 0;
  m_quad_type->clear();
  m_quad_charge->clear();
  m_quad_chi2->clear();
  m_quad_x->clear();
  m_quad_y->clear();
  m_quad_z->clear();

  m_quad_nCombined->clear();
  m_quad_id1->clear();
  m_quad_id2->clear();
  m_quad_id3->clear();
  m_quad_id4->clear();

  m_quad_mass->clear();
  m_quad_pt->clear();
  m_quad_eta->clear();
  m_quad_phi->clear();

  m_quad_track_mass->clear();
  m_quad_track_pt->clear();
  m_quad_track_eta->clear();
  m_quad_track_phi->clear();

  m_quad_fitted_mass->clear();
  m_quad_fitted_pt->clear();
  m_quad_fitted_eta->clear();
  m_quad_fitted_phi->clear();

  // clear electron information

  m_elec_eta    ->clear();
  m_elec_phi    ->clear();
  m_elec_pt     ->clear();
  m_elec_auth   ->clear();
  m_elec_charge ->clear();
  m_elec_id     ->clear();
  m_elec_ISEM   ->clear();

  m_elec_etCl   ->clear();
  m_elec_etas2  ->clear();
  m_elec_rhad   ->clear();
  m_elec_rhad1  ->clear();
  m_elec_reta   ->clear();
  m_elec_w2     ->clear();

  m_elec_f1     ->clear();
  m_elec_f3     ->clear();
  m_elec_wstot  ->clear();
  m_elec_DEmaxs1    ->clear();
  m_elec_deltaEta   ->clear();
  m_elec_deltaPhiRescaled   ->clear();

  m_elec_nSi    ->clear();
  m_elec_nSiDeadSensors ->clear();
  m_elec_nPix   ->clear();
  m_elec_nPixDeadSensors    ->clear();
  m_elec_nTRThigh   ->clear();
  m_elec_nTRThighOutliers   ->clear();
  m_elec_nTRT   ->clear();
  m_elec_nTRTOutliers   ->clear();

  m_elec_dpOverp    ->clear();
  m_elec_rTRT   ->clear();
  m_elec_expBLayer  ->clear();
  m_elec_trackd0    ->clear();

  //
  m_runNumber=0;
  m_eventNumber=0;
  m_eventTime=0;
  m_lumiBlock=0;
  m_bCID=0;
  m_lVL1ID=0;
  m_eventWeight=0;
  m_statusElement=0;
  m_lvl1TriggerType=0;
  m_lvl1TriggerInfo->clear();
  m_lvl2TriggerInfo->clear();
  m_evtFilterInfo->clear();
  m_streamTagName->clear();
  m_streamTagType->clear();
  //

  return StatusCode::SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////////
// Execute - called by the event loop on event by event

StatusCode VFitZmmOnAOD::execute() {

  ATH_MSG_DEBUG("execute()");
  m_cutFlow->Fill(1);

  StatusCode sc = StatusCode::SUCCESS;

  // initialize first before processing each event
  sc = initEvent();
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("initEvent failed. Continue");
  }

  // add event info to ntuple
  sc = addEventInfo();
  if ( sc.isFailure() ) {
    ATH_MSG_WARNING("could not get EventInfo. Continue");
  }

  // Get primary vertex from StoreGate
  sc = evtStore()->retrieve(vxCont,m_VxContainerName);
  if (sc.isFailure() ) {
    ATH_MSG_WARNING("No Vertex conainter with key "
		    << m_VxContainerName << " found in StoreGate. Continue");
  } else {
    ATH_MSG_VERBOSE("Found primary vertex info: " <<  m_VxContainerName);
    if(vxCont) {
      // int npvx = vxCont->size();
      VxContainer::const_iterator fz = vxCont->begin();
      const Trk::RecVertex& primaryVertex = (*fz)->recVertex();
      primVtx = &primaryVertex;
      
      int count_vtx = 0;
      for(VxContainer::const_iterator vtxItr = vxCont->begin();
              vtxItr != vxCont->end(); ++ vtxItr) 
      {
          if(count_vtx++ > 10) break;
          m_v0_x->push_back( (float) (*vtxItr)->recVertex().position().x());
          m_v0_y->push_back( (float) (*vtxItr)->recVertex().position().y());
          m_v0_z->push_back( (float) (*vtxItr)->recVertex().position().z());
      }
    }
  }

  // do electrons, if has one oppositely charged electron-pair
  StatusCode has_2OC_ele = zee_on_aod();
  if(has_2OC_ele.isSuccess()){
    ATH_MSG_DEBUG("has two electrons");
  }
  // do muons, if has two oppositely charge muon-pairs
  StatusCode has_2OC_muons = zmm_on_aod();

  if ( has_2OC_muons.isFailure() ) {
    ATH_MSG_DEBUG("no four neutral tracks available in the event");
    // check 2mu-2e case
    if(has_2OC_ele.isSuccess()){
        if( this->buildTwoMuonTwoEle(*m_good_muons, *m_good_electrons) ) {
            m_tree_Zll->Fill();
        }
    }
    return StatusCode::SUCCESS;
  } else {
      m_tree_Zll->Fill();
  }

  if(sc.isFailure()) sc = StatusCode::SUCCESS;
  return sc;

}


//////////////////////////////////////////////////////////////////////////////////
// zmm on aod: called by execute()

StatusCode VFitZmmOnAOD::zmm_on_aod() {

  ATH_MSG_DEBUG("zmm_on_aod()");

  StatusCode sc = StatusCode::SUCCESS;

  // retrieve the user container of pre selected muon from TDS
  const MuonContainer* muonTDS=0;
  sc= evtStore()->retrieve( muonTDS, m_muonContainerName);
  if( sc.isFailure()  ||  !muonTDS ) {
    ATH_MSG_WARNING("No AOD muon container of muons found in TDS");
    return StatusCode::SUCCESS;
  }  
  ATH_MSG_DEBUG("MuonContainer successfully retrieved");

  // iterators over the container
  MuonContainer::const_iterator muonItr  = muonTDS->begin();
  MuonContainer::const_iterator muonItrE = muonTDS->end();
  int imuon=0;

  int n_pos = 0;
  int n_neg = 0;
  for (; muonItr != muonItrE; ++muonItr) 
  {
     imuon++;
     TLorentzVector tlv = this->getLorentzVector( (**muonItr) );

     float charge = (*muonItr)->charge();
     bool passBLayer = false;

     // apply selections
     if( !this->passMuon(**muonItr, &passBLayer) ) continue;
     m_muons_cutFlow->Fill(3);
     if( !(*muonItr)->isCombinedMuon() and\
         !(*muonItr)->isSegmentTaggedMuon()){
        continue;
     }
     m_muons_cutFlow->Fill(4);
     // charge of ID track and muon
     const Rec::TrackParticle* id_track = (*muonItr)->inDetTrackParticle();

     bool consitent_charge = (id_track->charge() == charge);
     if(!consitent_charge) continue;
     m_muons_cutFlow->Fill(5);
     if(charge < 0) {
         n_neg ++;
     } else {
         n_pos ++;
     }

     m_good_muons->push_back( (*muonItr) );
     // Fill muon variables
     n_muon ++;
     mu_author_->push_back( (int) (*muonItr)->author() );
     mu_pt_->push_back(tlv.Pt());
     mu_eta_->push_back(tlv.Eta());
     mu_phi_->push_back(tlv.Phi());
     mu_e_->push_back(tlv.E());

     mu_charge_->push_back(charge);
     int type = 0;
     if( (*muonItr)->isCombinedMuon() ) type = 1;
     mu_type_->push_back(type);

     const Rec::TrackParticle* track = (*muonItr)->track();
     if (track){
         const Trk::MeasuredPerigee* aMeasPer = track->measuredPerigee();
         if(aMeasPer){
             float d0 = aMeasPer->parameters()[Trk::d0];
             mu_d0_->push_back(d0);
             float z0 = aMeasPer->parameters()[Trk::z0];
             mu_z0_->push_back(z0);
             float z0_sintheta = z0 * tlv.Theta();
             mu_z0_sintheta_->push_back(z0_sintheta);
             const Trk::ErrorMatrix& theErrorMatrix = aMeasPer->localErrorMatrix();
             float d0_error = (theErrorMatrix.covariance())[Trk::d0][Trk::d0];
             mu_d0_sig_->push_back( d0/sqrt(d0_error) );
         }
         // extrapolate to primary vertex
         // const Trk::TrackParameters* newPerigee = m_trackToVertexTool->perigeeAtVertex(*track, primVtx->position());
         const Trk::MeasuredPerigee* newPerigee = m_trackToVertexTool->perigeeAtVertex(*track, primVtx->position());
         if(newPerigee){
             mu_d0_pv->push_back( newPerigee->parameters()[Trk::d0] );
             mu_z0_pv->push_back( newPerigee->parameters()[Trk::z0] );
             const Trk::ErrorMatrix& theErrorMatrix = newPerigee->localErrorMatrix();
             float d0_error = (theErrorMatrix.covariance())[Trk::d0][Trk::d0];
             mu_d0_pv_sig_->push_back(newPerigee->parameters()[Trk::d0]/sqrt(d0_error));
             float z0_sintheta =  newPerigee->parameters()[Trk::z0] * tlv.Theta();
             mu_z0_pv_sintheta_->push_back(z0_sintheta);
         }
     }

     // d0 w.r.t primary vertex
     float   track_pt   = id_track->pt();
     float   track_eta  = id_track->eta();
     float   track_phi  = id_track->phi();
     float   track_e    = id_track->e();
     mu_track_pt_->push_back(track_pt);
     mu_track_eta_->push_back(track_eta);
     mu_track_phi_->push_back(track_phi);
     mu_track_e_->push_back(track_e);

     mu_eloss_->push_back( (*muonItr)->energyLoss().first );
     mu_etcone30_->push_back( (*muonItr)->parameter(MuonParameters::etcone30) );
     mu_ptcone30_->push_back( (*muonItr)->parameter(MuonParameters::ptcone30) );
     mu_pvID->push_back( this->matchPV( (**muonItr) ) );
     mu_blayer_->push_back(passBLayer);

     m_aod_muon_charge->Fill( (*muonItr)->charge() );
     m_aod_muon_pt->Fill( (*muonItr)->pt(), 1.);
     m_aod_muon_eta->Fill( (*muonItr)->eta(), 1.);
     m_aod_muon_chi2->Fill( (*muonItr)->matchChi2OverDoF(), 1.);
     ATH_MSG_DEBUG("imuon=" <<imuon<< ", AthenaBarCode=" << (*muonItr)->getAthenaBarCode());
  }

  // remove event with less than four muons
  if(n_muon < 4) return StatusCode::FAILURE;
  bool has_neutral_track = (n_pos >= 2 && n_neg >=2);
  if(!has_neutral_track) return StatusCode::FAILURE;
  m_cutFlow->Fill(2);

  this->buildTwoMuons( *m_good_muons );
  this->buildFourMuons( *m_good_muons );

  return StatusCode::SUCCESS;
}

// muon selection function
bool selectMuon(VFitZmmOnAOD *self, const MuonVect &ll) 
{
    bool test1 = ll[0]->charge() == -(ll[1]->charge());
    bool test2 = (ll[0]->pt() > self->m_etMuonCut) &&
        (ll[1]->pt() > self->m_etMuonCut);
    bool test3 = ( fabs(ll[0]->eta()) < self->m_etaMuonCut ) &&
        ( fabs(ll[1]->eta()) < self->m_etaMuonCut );
    return (test1 && test2 && test3);
}

bool VFitZmmOnAOD::passMuon(const Analysis::Muon& muon, bool* bLayer)
{
    m_muons_cutFlow->Fill(1);
    bool result = muon.pt() > this->m_etMuonCut && fabs(muon.eta()) < this->m_etaMuonCut;
    if(!result) return false;

    m_muons_cutFlow->Fill(2);
    // add track quality cuts
    const Rec::TrackParticle* id_track = muon.inDetTrackParticle();
    if(!id_track) return false;
    if(fabs(id_track->eta()) > 2.5) return false;

    const Trk::TrackSummary* trkSummary = id_track->trackSummary();
    if(!trkSummary) return false;
    int n_blayerHits = trkSummary->get(Trk::numberOfBLayerHits);
    int expectedBLayerHits = trkSummary->get(Trk::expectBLayerHit);

    int n_pixHits = trkSummary->get(Trk::numberOfPixelHits);
    int n_sctHits = trkSummary->get(Trk::numberOfSCTHits);
    int n_trtHits = trkSummary->get(Trk::numberOfTRTHits);
    int n_trtOutlier = trkSummary->get(Trk::numberOfTRTOutliers);
    int n_pixelHoles = trkSummary->get(Trk::numberOfPixelHoles);
    int n_sctHoles = trkSummary->get(Trk::numberOfSCTHoles);
    int n_pixelDeadSensor = trkSummary->get(Trk::numberOfPixelDeadSensors);
    int n_sctDeadSenesor = trkSummary->get(Trk::numberOfSCTDeadSensors);

    bool passBLayer = (!expectedBLayerHits || n_blayerHits > 0); // donot apply blayer fornow
    if(!bLayer) *bLayer = passBLayer;
    // result &= (n_pixHits + n_pixelDeadSensor) > 1;
    // result &= (n_sctHits + n_sctDeadSenesor) >= 6;
    result &= (n_pixHits + n_pixelDeadSensor) >= 1;
    result &= (n_sctHits + n_sctDeadSenesor) >= 5;
    result &= (n_pixelHoles + n_sctHoles) < 3;

    int n_trt = n_trtHits + n_trtOutlier;
    float eta = muon.eta();
    if (abs(eta) > 0.1 && abs(eta) < 1.9) {
        result &= (n_trt > 5 && n_trtOutlier < n_trt*0.9);
    } else if (n_trt > 5) {
        result &= (n_trtOutlier < n_trt * 0.9);
    }else {
    }

    return result;
}

TLorentzVector VFitZmmOnAOD::getLorentzVector(const Analysis::Muon& muon) {
    float pT  = muon.pt();
    float eta = muon.eta();
    float phi = muon.phi();
    float energy = muon.e();
    TLorentzVector tlv;
    tlv.SetPtEtaPhiE(pT, eta, phi, energy);
    return tlv;
}

TLorentzVector VFitZmmOnAOD::getLorentzVector(const Analysis::Electron& ele) {
    float pT  = ele.pt();
    float eta = ele.eta();
    float phi = ele.phi();
    float energy = ele.e();
    TLorentzVector tlv;
    tlv.SetPtEtaPhiE(pT, eta, phi, energy);
    return tlv;
}

TLorentzVector* VFitZmmOnAOD::getTrackLorentzV(const Analysis::Muon& muon)
{
    const Rec::TrackParticle* id_track = muon.inDetTrackParticle();
    if(id_track){
        float pT  = id_track->pt();
        float eta = id_track->eta();
        float phi = id_track->phi();
        float energy = id_track->e();
        TLorentzVector* tlv = new TLorentzVector;
        tlv->SetPtEtaPhiE(pT, eta, phi, energy);
        return tlv;
    } else {
        return NULL;
    }
}

TLorentzVector* VFitZmmOnAOD::getTrackLorentzV(const Analysis::Electron& electrons)
{
    const Rec::TrackParticle* id_track = electrons.trackParticle();
    if(id_track){
        float pT  = id_track->pt();
        float eta = id_track->eta();
        float phi = id_track->phi();
        float energy = id_track->e();
        TLorentzVector* tlv = new TLorentzVector;
        tlv->SetPtEtaPhiE(pT, eta, phi, energy);
        return tlv;
    } else {
        return NULL;
    }
}

///
StatusCode VFitZmmOnAOD::addEventInfo() {

    // this code has been taken from CBNT_execute
    // Reconstruction/CBNT_Athena/src/CBNTAA_EventInfo.cxx
    // Need this in the ntuple, but don't want to invoke the CBNT class
    // I have the actual EventNumber, but skipped the sequential count of event #
    // 

    //get EventInfo for run and event number

    const EventInfo* eventInfo;
    StatusCode sc = evtStore()->retrieve(eventInfo);

    if (sc.isFailure())
    {
        ATH_MSG_WARNING("Could not retrieve event info");
        return sc;
    }

    const EventID* myEventID=eventInfo->event_ID();
    //
    m_runNumber=myEventID->run_number();
    m_eventNumber=myEventID->event_number();
    m_eventTime= myEventID->time_stamp() ; 
    m_lumiBlock=myEventID->lumi_block() ;
    m_bCID=myEventID->bunch_crossing_id() ;

    const EventType* myEventType=eventInfo->event_type();
    if (myEventType!=0) {
        m_eventWeight=myEventType->mc_event_weight();
    }else
    {
        m_eventWeight=-999;
    }

    const TriggerInfo* myTriggerInfo=eventInfo->trigger_info();
    if (myTriggerInfo!=0) {
        m_lVL1ID=myTriggerInfo->extendedLevel1ID();
        m_statusElement=myTriggerInfo->statusElement();
        m_lvl1TriggerType=myTriggerInfo->level1TriggerType();

        std::vector<TriggerInfo::number_type>::const_iterator lvl1TrigIt=myTriggerInfo->level1TriggerInfo().begin();
        std::vector<TriggerInfo::number_type>::const_iterator lvl1TrigIt_e=myTriggerInfo->level1TriggerInfo().end();
        for (;lvl1TrigIt!=lvl1TrigIt_e;lvl1TrigIt++)
            m_lvl1TriggerInfo->push_back(*lvl1TrigIt);


        std::vector<TriggerInfo::number_type>::const_iterator lvl2TrigIt=myTriggerInfo->level2TriggerInfo().begin();
        std::vector<TriggerInfo::number_type>::const_iterator lvl2TrigIt_e=myTriggerInfo->level2TriggerInfo().end();
        for (;lvl2TrigIt!=lvl2TrigIt_e;lvl2TrigIt++)
            m_lvl2TriggerInfo->push_back(*lvl2TrigIt);

        std::vector<TriggerInfo::number_type>::const_iterator evtFilterIt=myTriggerInfo->eventFilterInfo().begin();
        std::vector<TriggerInfo::number_type>::const_iterator evtFilterIt_e=myTriggerInfo->eventFilterInfo().end();
        for (;evtFilterIt!=evtFilterIt_e;evtFilterIt++)
            m_evtFilterInfo->push_back(*evtFilterIt);


        std::vector<TriggerInfo::StreamTag>::const_iterator streamInfoIt=myTriggerInfo->streamTags().begin();
        std::vector<TriggerInfo::StreamTag>::const_iterator streamInfoIt_e=myTriggerInfo->streamTags().end();
        for (;streamInfoIt!=streamInfoIt_e;streamInfoIt++) { 
            m_streamTagName->push_back(streamInfoIt->name());
            m_streamTagType->push_back(streamInfoIt->type());
        }

    }else
    {
        m_lVL1ID=0;
        m_statusElement=0;
        m_lvl1TriggerType=0;
    }


    return StatusCode::SUCCESS;

}

Trk::VxCandidate* VFitZmmOnAOD::VkVrtFit(
        const MuonVect& muons, HepLorentzVector* momentum)
{
    m_VKVrtFitter->setDefault();
    std::vector<const Rec::TrackParticle *> myTracks;
    std::vector<const Trk::TrackParticleBase *> myTrackBases;
    std::vector<double> myMuonMasses;

    Trk::RecVertex primaryVtx = *primVtx;

    Hep3Vector  primVtxPos(primVtx->position().x(), primVtx->position().y(), primVtx->position().z());

    // int index = 1;
    for(MuonVect::const_iterator mu_itr = muons.begin();
            mu_itr != muons.end(); ++ mu_itr) {
        // const Rec::TrackParticle* tp = (*mu_itr)->track(); //Use primary track
        const Rec::TrackParticle* tp = (*mu_itr)->inDetTrackParticle(); //Use ID track
        if(tp) {
            myTracks.push_back(tp);
            myTrackBases.push_back(tp);
            myMuonMasses.push_back(mMuon);
            // indices.push_back(index);
            // index++;
        }
    }

    m_VKVrtFitter->setMassInputParticles(myMuonMasses);

    if (m_doMassConstraint && muons.size() == 4) {
        int upsilon_id1 = -1;
        int upsilon_id2 = -1;

        // only for four muons, try to find indices of muon pair coming from upsilon
        for(int i=0; i < (int) muons.size(); i++){
            const Analysis::Muon* muon1 = dynamic_cast<const Analysis::Muon*>(muons.at(i));
            float charge_muon1 = muon1->charge();
            float track_pt_1 = muon1->inDetTrackParticle()->pt();
            if(track_pt_1 < 4E3) continue;

            for(int j=i+1; j < (int)muons.size(); j++){
                const Analysis::Muon* muon2 = dynamic_cast<const Analysis::Muon*>(muons.at(j));
                float charge_muon2 = muon2->charge();
                float track_pt_2 = muon2->inDetTrackParticle()->pt();
                if(track_pt_2 < 4E3) continue;

                if( (charge_muon1 + charge_muon2) != 0) continue;
                Trk::VxCandidate* vx_can = VkVrtFit(*muon1, *muon2);
                float chi2 = getChi2( vx_can );
                delete vx_can;
                if (chi2 > 3) continue;

                TLorentzVector* tlv1 = this->getTrackLorentzV(*muon1);
                TLorentzVector* tlv2 = this->getTrackLorentzV(*muon2);
                float inv_m  = (*tlv1 + *tlv2).M();
                if(tlv1) delete  tlv1;
                if(tlv2) delete  tlv2;
                if(inv_m > 9.2E3 && inv_m < 9.7E3){
                    upsilon_id1 = i;
                    upsilon_id2 = j;
                    break;
                }
            }
            if(upsilon_id1 >= 0) break;
        }
        std::vector<int> indices;
        if (upsilon_id1 >= 0) {
            indices.push_back(upsilon_id1+1);
            indices.push_back(upsilon_id2+1);
        }
        m_VKVrtFitter->setMassForConstraint(UPSILON_MASS, indices);
        m_VKVrtFitter->setCnstType(1);
        m_VKVrtFitter->setRobustness(0);
    }

    return VkVrtFit(myTrackBases, momentum);
}

Trk::VxCandidate* VFitZmmOnAOD::VkVrtFit(
        std::vector<const Trk::TrackParticleBase*>& myTrackBases,
        HepLorentzVector* momentum)
{
    // m_VKVrtFitter->setDefault();
    // std::vector<const Rec::TrackParticle *> myTracks;
    // std::vector<const Trk::TrackParticleBase *> myTrackBases;

    Trk::RecVertex primaryVtx = *primVtx;

    Hep3Vector  primVtxPos(primVtx->position().x(), primVtx->position().y(), primVtx->position().z());

    // vertex contrained to the primary vertex
    if(m_doPVConstraint){
        m_VKVrtFitter->setVertexForConstraint(primaryVtx);
    }

    // define variables returned by the vertex fit
    Hep3Vector appVertex;
    m_VKVrtFitter->setMomCovCalc(0);  // No total momentum and its covariance matrix
    StatusCode sc = m_VKVrtFitter->VKalVrtFitFast(myTrackBases, appVertex);
    if (sc.isFailure()) {
        std::cout << "Warning from VKaVrt - fast fit failed!" << std::endl;
    }
    m_VKVrtFitter->setApproximateVertex(appVertex.x(),appVertex.y(),appVertex.z()); /*Use as starting point*/

    Hep3Vector finalVertex;
    // HepLorentzVector   momentum;
    long int charge = 0;
    double  fitChi2_vk=0.;
    std::vector<double> errMatrix,chi2PerTrk; // Fit error matrix and chi2 per track
    std::vector< std::vector<double> > trkAtVrt; // "True" tracks passing through vertex [phi, theta, q/pt*sin(theta)]
    m_VKVrtFitter->setMomCovCalc(1);  /* Total momentum and its covariance matrix are calculated*/

    if (momentum == NULL) {
        HepLorentzVector   momentum;
        sc = m_VKVrtFitter->VKalVrtFit(myTrackBases, finalVertex, momentum, charge,errMatrix,chi2PerTrk,trkAtVrt,fitChi2_vk);
    } else {
        sc = m_VKVrtFitter->VKalVrtFit(myTrackBases, finalVertex, *momentum, charge,errMatrix,chi2PerTrk,trkAtVrt,fitChi2_vk);
    }

    if (sc.isSuccess()) {
        Trk::TrkVKalVrtFitter* m_VKVFitter = dynamic_cast<Trk::TrkVKalVrtFitter*>(& (*m_VKVrtFitter));
        if(m_VKVFitter){
            return m_VKVFitter->makeVxCandidate(myTrackBases, finalVertex, errMatrix, chi2PerTrk, trkAtVrt, fitChi2_vk);
        } else {
            return NULL;
        }
    } else {
        return NULL;
    }
}

Trk::VxCandidate* VFitZmmOnAOD::VkVrtFit(const Analysis::Muon& muon1, const Analysis::Muon& muon2,
        HepLorentzVector* momentum)
{
    MuonVect* muons = new MuonVect;
    muons->push_back(&muon1);
    muons->push_back(&muon2);
    Trk::VxCandidate* results = this->VkVrtFit(*muons, momentum);
    delete muons;

    return results;
}

void VFitZmmOnAOD::buildFourMuons(const MuonVect& muons)
{
    for(int i = 0; i < (int) muons.size(); ++i) {
        const Analysis::Muon* muon1 = dynamic_cast<const Analysis::Muon*>( muons.at(i) );

        for(int j=i+1; j < (int) muons.size(); j++){
            const Analysis::Muon* muon2 = dynamic_cast<const Analysis::Muon*>( muons.at(j) );

            for( int k=j+1; k < (int) muons.size(); ++k){
                const Analysis::Muon* muon3 = dynamic_cast<const Analysis::Muon*>( muons.at(k) );

                for(int l=k+1; l < (int) muons.size(); ++l){
                    const Analysis::Muon* muon4 = dynamic_cast<const Analysis::Muon*>( muons.at(l) );
                    if( (muon1->charge() + muon2->charge() + muon3->charge() + muon4->charge()) != 0) continue;
                    // so far we have four neutral tracks

                    // require at least three combined muons
                    int n_combined = 0;
                    if( muon1->isCombinedMuon() ) n_combined ++;
                    if( muon2->isCombinedMuon() ) n_combined ++;
                    if( muon3->isCombinedMuon() ) n_combined ++;
                    if( muon4->isCombinedMuon() ) n_combined ++;
                    // if(n_combined < 3) continue;

                    m_quad_nCombined->push_back(n_combined);

                    MuonVect muons_can;
                    muons_can.push_back( muon1 );
                    muons_can.push_back( muon2 );
                    muons_can.push_back( muon3 );
                    muons_can.push_back( muon4 );

                    HepLorentzVector* momentum = new HepLorentzVector();
                    Trk::VxCandidate* vx_can = this->VkVrtFit( muons_can, momentum);

                    this->fillQuadInfo(vx_can, momentum, 0);
                    if(momentum) delete momentum;
                    if(vx_can) delete vx_can;

                    // four-momentum from combined measurement 
                    TLorentzVector tlv1 = this->getLorentzVector(*muon1) ;
                    TLorentzVector tlv2 = this->getLorentzVector(*muon2) ;
                    TLorentzVector tlv3 = this->getLorentzVector(*muon3) ;
                    TLorentzVector tlv4 = this->getLorentzVector(*muon4) ;
                    TLorentzVector tlv_total = (tlv1 + tlv2 + tlv3 + tlv4);
                    m_quad_pt->push_back( tlv_total.Pt() );
                    m_quad_eta->push_back( tlv_total.Eta() );
                    m_quad_phi->push_back( tlv_total.Phi() );
                    m_quad_mass->push_back( tlv_total.M() );

                    // four-momentum from track particle
                    TLorentzVector* track_tlv1 = this->getTrackLorentzV(*muon1);
                    TLorentzVector* track_tlv2 = this->getTrackLorentzV(*muon2);
                    TLorentzVector* track_tlv3 = this->getTrackLorentzV(*muon3);
                    TLorentzVector* track_tlv4 = this->getTrackLorentzV(*muon4);
                    TLorentzVector track_tlv_total = (*track_tlv1 + *track_tlv2 + *track_tlv3 + *track_tlv4);
                    m_quad_track_pt->push_back( track_tlv_total.Pt() );
                    m_quad_track_eta->push_back( track_tlv_total.Eta() );
                    m_quad_track_phi->push_back( track_tlv_total.Phi() );
                    m_quad_track_mass->push_back( track_tlv_total.M() );
                    delete track_tlv1;
                    delete track_tlv2;
                    delete track_tlv3;
                    delete track_tlv4;

                    m_quad_id1->push_back(i);
                    m_quad_id2->push_back(j);
                    m_quad_id3->push_back(k);
                    m_quad_id4->push_back(l);
                }
            }
        }
    }
}

bool VFitZmmOnAOD::passOniaCuts(
        const Analysis::Muon& muon1, const Analysis::Muon& muon2,
        const Analysis::Muon& muon3, const Analysis::Muon& muon4,
        float& m_light, float& m_heavy,
        float& chi2_ndf_1, float& chi2_ndf_2 
        )
{
    Trk::VxCandidate* vx_can1 = VkVrtFit(muon1, muon2);
    chi2_ndf_1= getChi2( vx_can1 );
    if (vx_can1) delete vx_can1;

    Trk::VxCandidate* vx_can2 = VkVrtFit(muon3, muon4);
    chi2_ndf_2= getChi2( vx_can2 );
    if (vx_can2) delete vx_can2;

    m_chi2_fit_12->Fill(chi2_ndf_1);
    m_chi2_fit_34->Fill(chi2_ndf_2);

    if(chi2_ndf_1 == -1 || chi2_ndf_2 == -1 || chi2_ndf_1 >= 3 || chi2_ndf_2 >=3 ) 
        return false;

    const Rec::TrackParticle* id_track_1 = muon1.inDetTrackParticle();
    const Rec::TrackParticle* id_track_2 = muon2.inDetTrackParticle();
    const Rec::TrackParticle* id_track_3 = muon3.inDetTrackParticle();
    const Rec::TrackParticle* id_track_4 = muon4.inDetTrackParticle();

    TLorentzVector tlv_1 = this->getLorentzVector( muon1 );
    TLorentzVector tlv_2 = this->getLorentzVector( muon2 );
    TLorentzVector tlv_3 = this->getLorentzVector( muon3 );
    TLorentzVector tlv_4 = this->getLorentzVector( muon4 );

    float m_12 = (tlv_1 + tlv_2).M();
    float m_34 = (tlv_3 + tlv_4).M();

    // check if first pair pass Onia-1 cut: track pT > 4 GeV
    float m_onia1 = -999;
    float m_onia2 = -999;

    if(id_track_1 && id_track_1->pt() > ONIA_ONE_PT_CUT &&\
            id_track_2 && id_track_2->pt() > ONIA_ONE_PT_CUT)
    {
        m_onia1 = m_12;
        if (id_track_3 && id_track_3->pt() > ONIA_TWO_PT_CUT &&\
                id_track_4 && id_track_4->pt() > ONIA_TWO_PT_CUT)
        {
            m_onia2 = m_34;
        }
    }

    if (id_track_3 && id_track_3->pt() > ONIA_ONE_PT_CUT &&\
            id_track_4 && id_track_4->pt() > ONIA_ONE_PT_CUT && m_34 > m_onia1) 
    {
        if(m_onia1 < 0){
            if(id_track_1 && id_track_1->pt() > ONIA_TWO_PT_CUT &&\
                    id_track_2 && id_track_2->pt() > ONIA_TWO_PT_CUT){
                m_onia2 = m_12;
            }
        } else {
            m_onia2 = m_12;
        }
        m_onia1 = m_34;
    }
    m_heavy = m_onia1;
    m_light = m_onia2;
    if( m_onia1 > 0 && m_onia2 > 0) return true;
    else return false;
}

int VFitZmmOnAOD::matchPV(const Analysis::Muon& muon)
{
    // return the index of the primary vertex that matches the vertex associated with the ID track of the muon
    int res = -1;

    const Rec::TrackParticle* track = muon.inDetTrackParticle();
    if(!track) return res;
    const Trk::VxCandidate* muon_vxCan = track->reconstructedVertex();
    if(!muon_vxCan) return res;

    float mu_d0 = muon_vxCan->recVertex().position()[Trk::d0];
    float mu_z0 = muon_vxCan->recVertex().position()[Trk::z0];

    int index_pv = -1;
    for(VxContainer::const_iterator vtxItr = vxCont->begin();
            vtxItr != vxCont->end(); ++ vtxItr)
    {
        index_pv ++;
        float vxt_d0 = (*vtxItr)->recVertex().position()[Trk::d0];
        float vxt_z0 = (*vtxItr)->recVertex().position()[Trk::z0];
        if (fabs(vxt_d0 - mu_d0) < 1E-5 && fabs(vxt_z0 - mu_z0) < 1E-5) {
            // ATH_MSG_INFO("found a match: " << index_pv);
            res = index_pv;
            break;
        }
    }
    return res;
}

void VFitZmmOnAOD::buildTwoMuons(const MuonVect& muons)
{
    for(int i = 0; i < (int) muons.size(); i++) {
        const Analysis::Muon* muon1 = dynamic_cast<const Analysis::Muon*>( muons.at(i) );
        float mu_charge_1 = muon1->charge();
        for(int j = i+1; j < (int) muons.size(); ++j) {
            const Analysis::Muon* muon2 = dynamic_cast<const Analysis::Muon*>( muons.at(j) );
            float mu_charge_2 = muon2->charge();
            if( (mu_charge_1 + mu_charge_2 ) != 0) continue;

            // perform the fit
            HepLorentzVector* momentum = new HepLorentzVector();
            Trk::VxCandidate* vx_can = this->VkVrtFit(*muon1, *muon2, momentum);

            if(! this->fillOniaInfo(vx_can, momentum, 0)) continue;
            delete momentum;
            delete vx_can;
            m_onia_muon1id->push_back(i);
            m_onia_muon2id->push_back(j);

            // four-momentum from combined measurement 
            TLorentzVector tlv1 = this->getLorentzVector(*muon1);
            TLorentzVector tlv2 = this->getLorentzVector(*muon2);
            TLorentzVector tlv_total = (tlv1 + tlv2);
            m_onia_pt->push_back( tlv_total.Pt() );
            m_onia_eta->push_back( tlv_total.Eta() );
            m_onia_phi->push_back( tlv_total.Phi() );
            m_onia_mass->push_back( tlv_total.M() );

            // four-momentum from track particle
            TLorentzVector* track_tlv1 = this->getTrackLorentzV(*muon1);
            TLorentzVector* track_tlv2 = this->getTrackLorentzV(*muon2);
            TLorentzVector track_tlv_total = (*track_tlv1 + *track_tlv2);
            m_onia_track_pt->push_back( track_tlv_total.Pt() );
            m_onia_track_eta->push_back( track_tlv_total.Eta() );
            m_onia_track_phi->push_back( track_tlv_total.Phi() );
            m_onia_track_mass->push_back( track_tlv_total.M() );
            delete track_tlv1;
            delete track_tlv2;

        }
    }
}

bool VFitZmmOnAOD::buildTwoElectrons(const ElecVect& electrons)
{
    bool has_2OC = false;
    for(int i = 0; i < (int) electrons.size(); i++){
        const Analysis::Electron* ele1 = dynamic_cast<const Analysis::Electron*>( electrons.at(i) );
        float ele_charge_1 = ele1->charge();
        for(int j = i+1; j < (int) electrons.size(); ++j){
            const Analysis::Electron* ele2 = dynamic_cast<const Analysis::Electron*>( electrons.at(j) );
            float ele_charge_2 = ele2->charge();

            if( (ele_charge_1 + ele_charge_2) != 0) continue;

            std::vector<const Trk::TrackParticleBase*>  myTracks;
            myTracks.push_back( ele1->trackParticle() );
            myTracks.push_back( ele2->trackParticle() );

            HepLorentzVector* momentum = new HepLorentzVector();
            Trk::VxCandidate* vx_can = this->VkVrtFit(myTracks, momentum);
            if( ! this->fillOniaInfo(vx_can, momentum, 1) ){ 
                delete momentum;
                continue;
            }
            delete momentum;
            delete vx_can;
            m_onia_muon1id->push_back(i);
            m_onia_muon2id->push_back(j);
            has_2OC = true;

            // four-momentum from combined measurement 
            TLorentzVector tlv1 = this->getLorentzVector(*ele1);
            TLorentzVector tlv2 = this->getLorentzVector(*ele2);
            TLorentzVector tlv_total = (tlv1 + tlv2);
            m_onia_pt->push_back( tlv_total.Pt() );
            m_onia_eta->push_back( tlv_total.Eta() );
            m_onia_phi->push_back( tlv_total.Phi() );
            m_onia_mass->push_back( tlv_total.M() );

            // four-momentum from track particle
            TLorentzVector* track_tlv1 = this->getTrackLorentzV(*ele1);
            TLorentzVector* track_tlv2 = this->getTrackLorentzV(*ele2);
            TLorentzVector track_tlv_total = (*track_tlv1 + *track_tlv2);
            m_onia_track_pt->push_back( track_tlv_total.Pt() );
            m_onia_track_eta->push_back( track_tlv_total.Eta() );
            m_onia_track_phi->push_back( track_tlv_total.Phi() );
            m_onia_track_mass->push_back( track_tlv_total.M() );
            delete track_tlv1;
            delete track_tlv2;
        }
    }
    return has_2OC;
}

bool VFitZmmOnAOD::fillOniaInfo( Trk::VxCandidate* vx_can, HepLorentzVector* momentum, int type)
{
    if(!vx_can || !momentum) return false;
    m_n_onia ++;
    m_onia_type->push_back(type);

    m_onia_chi2->push_back( getChi2( vx_can ) );
    m_onia_x->push_back( vx_can->recVertex().position()[0] );
    m_onia_y->push_back( vx_can->recVertex().position()[1] );
    m_onia_z->push_back( vx_can->recVertex().position()[2] );

    m_onia_pt_fitted->push_back(momentum->perp());
    m_onia_eta_fitted->push_back(momentum->pseudoRapidity());
    m_onia_phi_fitted->push_back(momentum->phi());
    m_onia_mass_fitted->push_back(momentum->m());

    return true;
}

void VFitZmmOnAOD::fillOniaInfo(
        const Analysis::Muon& muon1,
        const Analysis::Muon& muon2
        )
{
    m_n_onia ++;
    m_onia_charge->push_back(muon1.charge() + muon2.charge());

    HepLorentzVector* momentum = new HepLorentzVector();
    Trk::VxCandidate* vx_can = this->VkVrtFit(muon1, muon2, momentum);
    m_onia_chi2->push_back( getChi2( vx_can ) );
    if(vx_can){ 
        m_onia_x->push_back( vx_can->recVertex().position()[0] );
        m_onia_y->push_back( vx_can->recVertex().position()[1] );
        m_onia_z->push_back( vx_can->recVertex().position()[2] );

        m_onia_pt_fitted->push_back(momentum->perp());
        m_onia_eta_fitted->push_back(momentum->pseudoRapidity());
        m_onia_phi_fitted->push_back(momentum->phi());
        m_onia_mass_fitted->push_back(momentum->m());
        delete vx_can;
    }
    delete momentum;

    // four-momentum from combined measurement 
    TLorentzVector tlv1 = this->getLorentzVector(muon1);
    TLorentzVector tlv2 = this->getLorentzVector(muon2);
    TLorentzVector tlv_total = (tlv1 + tlv2);
    m_onia_pt->push_back( tlv_total.Pt() );
    m_onia_eta->push_back( tlv_total.Eta() );
    m_onia_phi->push_back( tlv_total.Phi() );
    m_onia_mass->push_back( tlv_total.M() );

    // four-momentum from track particle
    TLorentzVector* track_tlv1 = this->getTrackLorentzV(muon1);
    TLorentzVector* track_tlv2 = this->getTrackLorentzV(muon2);
    TLorentzVector track_tlv_total = (*track_tlv1 + *track_tlv2);
    m_onia_track_pt->push_back( track_tlv_total.Pt() );
    m_onia_track_eta->push_back( track_tlv_total.Eta() );
    m_onia_track_phi->push_back( track_tlv_total.Phi() );
    m_onia_track_mass->push_back( track_tlv_total.M() );
    delete track_tlv1;
    delete track_tlv2;

}

float VFitZmmOnAOD::getChi2(Trk::VxCandidate* vx_can)
{
    float res = -1;
    if(vx_can) {
        const Trk::FitQuality& fitQuality = vx_can->recVertex().fitQuality();
        res = fitQuality.chiSquared()/fitQuality.numberDoF();
    }
    return res;
}

void VFitZmmOnAOD::fillQuadInfo(
        const Analysis::Muon& muon1,
        const Analysis::Muon& muon2,
        const Analysis::Muon& muon3,
        const Analysis::Muon& muon4)
{
    m_n_quad ++;
    m_quad_type->push_back(0);
    m_quad_charge->push_back(muon1.charge() + muon2.charge() + muon3.charge() + muon4.charge());
    HepLorentzVector* momentum = new HepLorentzVector();
    MuonVect* muons_can = new MuonVect;
    muons_can->push_back( &muon1 );
    muons_can->push_back( &muon2 );
    muons_can->push_back( &muon3 );
    muons_can->push_back( &muon4 );

    Trk::VxCandidate* vx_can = this->VkVrtFit( *muons_can, momentum);
    m_quad_chi2->push_back( getChi2( vx_can ) );
    if(vx_can){
        m_quad_x->push_back( vx_can->recVertex().position()[0] );
        m_quad_y->push_back( vx_can->recVertex().position()[1] );
        m_quad_z->push_back( vx_can->recVertex().position()[2] );

        m_quad_fitted_pt->push_back(momentum->perp());
        m_quad_fitted_eta->push_back(momentum->pseudoRapidity());
        m_quad_fitted_phi->push_back(momentum->phi());
        m_quad_fitted_mass->push_back(momentum->m());
        delete vx_can;
    }
    delete momentum;
    delete muons_can;

    // four-momentum from combined measurement 
    TLorentzVector tlv1 = this->getLorentzVector(muon1) ;
    TLorentzVector tlv2 = this->getLorentzVector(muon2) ;
    TLorentzVector tlv3 = this->getLorentzVector(muon3) ;
    TLorentzVector tlv4 = this->getLorentzVector(muon4) ;
    TLorentzVector tlv_total = (tlv1 + tlv2 + tlv3 + tlv4);
    m_quad_pt->push_back( tlv_total.Pt() );
    m_quad_eta->push_back( tlv_total.Eta() );
    m_quad_phi->push_back( tlv_total.Phi() );
    m_quad_mass->push_back( tlv_total.M() );

    // four-momentum from track particle
    TLorentzVector* track_tlv1 = this->getTrackLorentzV(muon1);
    TLorentzVector* track_tlv2 = this->getTrackLorentzV(muon2);
    TLorentzVector* track_tlv3 = this->getTrackLorentzV(muon3);
    TLorentzVector* track_tlv4 = this->getTrackLorentzV(muon4);
    TLorentzVector track_tlv_total = (*track_tlv1 + *track_tlv2 + *track_tlv3 + *track_tlv4);
    m_quad_track_pt->push_back( track_tlv_total.Pt() );
    m_quad_track_eta->push_back( track_tlv_total.Eta() );
    m_quad_track_phi->push_back( track_tlv_total.Phi() );
    m_quad_track_mass->push_back( track_tlv_total.M() );
    delete track_tlv1;
    delete track_tlv2;
    delete track_tlv3;
    delete track_tlv4;
}

bool VFitZmmOnAOD::fillQuadInfo(Trk::VxCandidate* vx_can, HepLorentzVector* momentum, int type)
{
    m_n_quad ++;
    m_quad_type->push_back(type);
    m_quad_chi2->push_back( getChi2( vx_can ) );

    if(!vx_can || !momentum){ 
        m_quad_x->push_back(-9999);
        m_quad_y->push_back(-9999);
        m_quad_z->push_back(-9999);
        m_quad_fitted_pt->push_back(-9999);
        m_quad_fitted_phi->push_back(-9999);
        m_quad_fitted_eta->push_back(-9999);
        m_quad_fitted_mass->push_back(-9999);
        return false;
    }

    // m_quad_charge->push_back(muon1.charge() + muon2.charge() + muon3.charge() + muon4.charge());
    m_quad_x->push_back( vx_can->recVertex().position()[0] );
    m_quad_y->push_back( vx_can->recVertex().position()[1] );
    m_quad_z->push_back( vx_can->recVertex().position()[2] );

    m_quad_fitted_pt->push_back(momentum->perp());
    m_quad_fitted_eta->push_back(momentum->pseudoRapidity());
    m_quad_fitted_phi->push_back(momentum->phi());
    m_quad_fitted_mass->push_back(momentum->m());

    return true;
}

StatusCode VFitZmmOnAOD::zee_on_aod() {
    ATH_MSG_DEBUG("zee_on_aod()");

    const ElectronContainer* elecTES;
    StatusCode sc = StatusCode::SUCCESS;

    sc = evtStore()->retrieve(elecTES, m_electronContainerName);
    if( sc.isFailure() || !elecTES) {
        ATH_MSG_WARNING("No AOD muon container of muons found in TDS");
        return StatusCode::SUCCESS;
    }
    ATH_MSG_DEBUG("ElectronContainer successfully retrieved");


    // iterator over the electrons
    ElectronContainer::const_iterator elecItr = elecTES->begin();
    ElectronContainer::const_iterator elecItrE = elecTES->end();

    int ielec = 0;
    int n_pos = 0;
    int n_neg = 0;
    for (; elecItr != elecItrE; ++elecItr){
        ielec ++;

        int author = (*elecItr)->author();
        if(author != 1 && author != 3) continue;
        if(! (*elecItr)->trackParticle() ) continue;

        const CaloCluster* cluster = (*elecItr)->cluster();
        if(!cluster) continue;

        float track_eta = fabs( (*elecItr)->trackParticle()->eta() );
        float eta = fabs( cluster->etaBE(2) );
        if(track_eta > 2.47) continue;

        float charge = (*elecItr)->charge();
        if(charge < 0) n_neg ++;
        else n_pos ++;

        n_el ++;
        m_good_electrons->push_back( (*elecItr) );
        bool pass_ID = (*elecItr)->passID(egammaPID::ElectronIDLoosePP);

        m_elec_eta->push_back( (*elecItr)->eta() );
        m_elec_phi->push_back( (*elecItr)->phi() );
        m_elec_pt->push_back( (*elecItr)->pt() );
        m_elec_auth->push_back( (*elecItr)->author() );
        m_elec_charge->push_back(charge);
        m_elec_id->push_back(pass_ID);
        
        // eID
        float et = cluster->energy()/cosh(eta);
        const Electron* eg = dynamic_cast<const Electron*>( *elecItr );
        float e237 = eg->detailValue(egammaParameters::e237);
        float e277   = eg->detailValue(egammaParameters::e277);
        float Reta = e277 != 0 ? e237/e277 : 0.;
        float weta2c = eg->detailValue(egammaParameters::weta2);

        double ethad1 = eg->detailValue(egammaParameters::ethad1);
        double ethad = eg->detailValue(egammaParameters::ethad);
        double raphad1 = fabs(et) != 0. ? ethad1/et : 0.;
        double raphad  = fabs(et) != 0. ? ethad/et : 0.;

        // E of 2nd max between max and min in strips
        double emax2  = eg->detailValue(egammaParameters::e2tsts1);
        // E of 1st max in strips
        double emax   = eg->detailValue(egammaParameters::emaxs1);
        // fraction of energy reconstructed in the 1st sampling
        double f1     = eg->detailValue(egammaParameters::f1);
        double f3     = eg->detailValue(egammaParameters::f3);
        // E of 1st max in strips
        double wtot   = eg->detailValue(egammaParameters::wtots1);
        // (Emax1-Emax2)/(Emax1+Emax2)
        double demaxs1 = (emax+emax2)==0. ? 0 : (emax-emax2)/(emax+emax2);

        // number of precision hits (Pixels+SCT)
        int nSi = 0;
        int nPi = 0;
        int nSiOutliers = 0; 
        int nPiOutliers = 0; 
        int nTRThigh          = 0;
        int nTRThighOutliers  = 0;
        int nTRT         = 0;
        int nTRTOutliers = 0;

        // retrieve associated track
        const Rec::TrackParticle* track  = eg->trackParticle();    

        const Trk::TrackSummary* summary = track->trackSummary();

        if (summary) {
            nPi = summary->get(Trk::numberOfPixelHits);
            nSi = summary->get(Trk::numberOfPixelHits)+summary->get(Trk::numberOfSCTHits);
            nPiOutliers = summary->get(Trk::numberOfPixelOutliers);   
            nSiOutliers = summary->get(Trk::numberOfPixelOutliers)+summary->get(Trk::numberOfSCTOutliers);

            nTRThigh          = summary->get(Trk::numberOfTRTHighThresholdHits);
            nTRThighOutliers  = summary->get(Trk::numberOfTRTHighThresholdOutliers);
            nTRT         = summary->get(Trk::numberOfTRTHits);
            nTRTOutliers = summary->get(Trk::numberOfTRTOutliers); 

        }
        // delta eta
        double deltaEta = fabs(eg->detailValue(egammaParameters::deltaEta1));
        double deltaphi = eg->detailValue(egammaParameters::deltaPhi2);

        bool expectHitInBLayer = eg->detailValue(egammaParameters::expectHitInBLayer);
        double trackd0 = fabs(eg->detailValue(egammaParameters::trackd0_physics));
        double rTRT = (nTRT+nTRTOutliers) > 0 ? ((double)(nTRThigh+nTRThighOutliers)/(nTRT+nTRTOutliers) ) : 0.;
        // E/p
        double trackp = track->p();
        double ep = 0.;
        if(fabs(trackp) > 0) ep = cluster->energy()/trackp;

        m_elec_etCl ->push_back(et);
        m_elec_etas2    ->push_back( cluster->etaBE(2)  );
        m_elec_rhad ->push_back(raphad);
        m_elec_rhad1    ->push_back(raphad1);
        m_elec_reta ->push_back(Reta);
        m_elec_w2   ->push_back(weta2c);

        m_elec_f1   ->push_back(f1);
        m_elec_f3   ->push_back(f3);
        m_elec_wstot    ->push_back(wtot);
        m_elec_DEmaxs1  ->push_back(demaxs1);
        m_elec_deltaEta ->push_back(deltaEta);
        m_elec_deltaPhiRescaled ->push_back(deltaphi);

        m_elec_nSi  ->push_back(nSi);
        m_elec_nSiDeadSensors   ->push_back( summary->get(Trk::numberOfSCTDeadSensors) + summary->get(Trk::numberOfPixelDeadSensors) );
        m_elec_nPix ->push_back(nPi);
        m_elec_nPixDeadSensors  ->push_back( summary->get(Trk::numberOfPixelDeadSensors) );
        m_elec_nTRThigh ->push_back(nTRThigh);
        m_elec_nTRThighOutliers ->push_back(nTRThighOutliers);
        m_elec_nTRT ->push_back(nTRT);
        m_elec_nTRTOutliers ->push_back(nTRTOutliers);

        m_elec_dpOverp  ->push_back(ep);
        m_elec_rTRT ->push_back(rTRT);
        m_elec_expBLayer    ->push_back(expectHitInBLayer);
        m_elec_trackd0  ->push_back(trackd0);
    }

    if (n_el < 2) return StatusCode::FAILURE;
    bool has_neutral_track = (n_pos >= 1 && n_neg >=1);
    if(!has_neutral_track) return StatusCode::FAILURE;

    // Fill onia
    this->buildTwoElectrons( *m_good_electrons );
    return StatusCode::SUCCESS;
}

bool VFitZmmOnAOD::buildTwoMuonTwoEle(const MuonVect& muons, const ElecVect& electrons)
{
    bool has_quadruplet = false;
    for(int i = 0; i < (int) muons.size(); i++){
        const Analysis::Muon* muon1 = dynamic_cast<const Analysis::Muon*>( muons.at(i) );
        if(! muon1->isCombinedMuon()) continue;

        for(int j = i+1; j < (int) muons.size(); j++){
            const Analysis::Muon* muon2 = dynamic_cast<const Analysis::Muon*>( muons.at(j) );
            if(! muon2->isCombinedMuon()) continue;
            if( (muon1->charge() + muon2->charge()) != 0) continue;

            for(int k = 0; k < (int) electrons.size(); k++){
                const Analysis::Electron* ele1 = dynamic_cast<const Analysis::Electron*>( electrons.at(k) );

                for(int l = k+1; l < (int) electrons.size(); l++){
                    const Analysis::Electron* ele2 = dynamic_cast<const Analysis::Electron*>( electrons.at(l) );
                    if((ele1->charge() + ele2->charge()) != 0) continue;

                    // start to build
                    std::vector<const Trk::TrackParticleBase*> myTracks;
                    myTracks.push_back( muon1->inDetTrackParticle() );
                    myTracks.push_back( muon2->inDetTrackParticle() );
                    myTracks.push_back( ele1->trackParticle() );
                    myTracks.push_back( ele2->trackParticle() );

                    HepLorentzVector* momentum = new HepLorentzVector();
                    Trk::VxCandidate* vx_can = this->VkVrtFit(myTracks, momentum);
                    this->fillQuadInfo(vx_can, momentum, 1);
                    if (momentum) delete momentum;
                    if (vx_can) delete vx_can;
                    has_quadruplet = true;

                    // four-momentum from combined measurement
                    TLorentzVector tlv1 = this->getLorentzVector(*muon1) ;
                    TLorentzVector tlv2 = this->getLorentzVector(*muon2) ;
                    TLorentzVector tlv3 = this->getLorentzVector(*ele1) ;
                    TLorentzVector tlv4 = this->getLorentzVector(*ele2) ;
                    TLorentzVector tlv_total = (tlv1 + tlv2 + tlv3 + tlv4);
                    m_quad_pt->push_back( tlv_total.Pt() );
                    m_quad_eta->push_back( tlv_total.Eta() );
                    m_quad_phi->push_back( tlv_total.Phi() );
                    m_quad_mass->push_back( tlv_total.M() );

                    // four-momentum from track particle
                    TLorentzVector* track_tlv1 = this->getTrackLorentzV(*muon1);
                    TLorentzVector* track_tlv2 = this->getTrackLorentzV(*muon2);
                    TLorentzVector* track_tlv3 = this->getTrackLorentzV(*ele1);
                    TLorentzVector* track_tlv4 = this->getTrackLorentzV(*ele2);
                    TLorentzVector track_tlv_total = (*track_tlv1 + *track_tlv2 + *track_tlv3 + *track_tlv4);
                    m_quad_track_pt->push_back( track_tlv_total.Pt() );
                    m_quad_track_eta->push_back( track_tlv_total.Eta() );
                    m_quad_track_phi->push_back( track_tlv_total.Phi() );
                    m_quad_track_mass->push_back( track_tlv_total.M() );
                    delete track_tlv1;
                    delete track_tlv2;
                    delete track_tlv3;
                    delete track_tlv4;

                    m_quad_id1->push_back(i);
                    m_quad_id2->push_back(j);
                    m_quad_id3->push_back(k);
                    m_quad_id4->push_back(l);
                }
            }
        }
    }
    return has_quadruplet;
}
