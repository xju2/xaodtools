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
  
  // muons
  mu_author_ = new std::vector<int>;
  mu_pt_ = new std::vector<float>;
  mu_eta_ = new std::vector<float>;
  mu_phi_ = new std::vector<float>;
  mu_e_ = new std::vector<float>;
  mu_track_pt_ = new std::vector<float>();

  mu_charge_ = new std::vector<float>;
  mu_type_ = new std::vector<int>;
  mu_d0_ = new std::vector<float>;
  mu_z0_sintheta_ = new std::vector<float>;
  mu_d0_sig_ = new std::vector<float>;
  mu_eloss_ = new std::vector<float>;
  mu_etcone30_ = new std::vector<float>;
  mu_ptcone30_ = new std::vector<float>;

  m_upsilon_ = 0;
  m_4l_ = -9999;
  vtx4l_chi2ndf_ = -9999;
  m34_ = -9999;
  same_vertex_ = false;
  m_4l_fitted_ = -9999;
  n_combined_muons_ = 0;
  m_index_1_ = -1;
  m_index_2_ = -1;
  m_index_3_ = -1;
  m_index_4_ = -1;

  m_pvID_1_ = -1;
  m_pvID_2_ = -1;
  m_pvID_3_ = -1;
  m_pvID_4_ = -1;


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
  m_tree_Zll->Branch("v0_x", &m_v0_x, "v0_x/D");
  m_tree_Zll->Branch("v0_y", &m_v0_y, "v0_y/D");
  m_tree_Zll->Branch("v0_z", &m_v0_z, "v0_z/D");

  // muon variables
  m_tree_Zll->Branch("n_muon", &n_muon, "n_muon/I");
  m_tree_Zll->Branch("mu_author", &mu_author_);
  m_tree_Zll->Branch("mu_pt", &mu_pt_);
  m_tree_Zll->Branch("mu_track_pt", &mu_track_pt_);
  m_tree_Zll->Branch("mu_eta", &mu_eta_);
  m_tree_Zll->Branch("mu_phi", &mu_phi_);
  m_tree_Zll->Branch("mu_e", &mu_e_);
  m_tree_Zll->Branch("mu_charge", &mu_charge_);
  m_tree_Zll->Branch("mu_type", &mu_type_);
  m_tree_Zll->Branch("mu_d0", &mu_d0_);
  m_tree_Zll->Branch("mu_z0_sintheta", &mu_z0_sintheta_);
  m_tree_Zll->Branch("mu_d0_sig", &mu_d0_sig_);
  m_tree_Zll->Branch("mu_eloss", &mu_eloss_);
  m_tree_Zll->Branch("mu_etcone30", &mu_etcone30_);
  m_tree_Zll->Branch("mu_ptvarcone30", &mu_ptcone30_);


  m_tree_Zll->Branch("mUpsilon", &m_upsilon_, "mUpsilon/F");
  m_tree_Zll->Branch("m4l", &m_4l_, "m4l/F");
  m_tree_Zll->Branch("m4l_fitted", &m_4l_fitted_, "m4l_fitted/F");
  m_tree_Zll->Branch("vtx4l_chi2ndf", &vtx4l_chi2ndf_, "vtx4l_chi2ndf/F");
  m_tree_Zll->Branch("m34", &m34_, "m34/F");
  m_tree_Zll->Branch("same_vertex", &same_vertex_, "same_vertex/O");
  m_tree_Zll->Branch("n_combined_muons", &n_combined_muons_, "n_combined_muons/I");

  m_tree_Zll->Branch("index_1", &m_index_1_, "index_1/I");
  m_tree_Zll->Branch("index_2", &m_index_2_, "index_2/I");
  m_tree_Zll->Branch("index_3", &m_index_3_, "index_3/I");
  m_tree_Zll->Branch("index_4", &m_index_4_, "index_4/I");
  m_tree_Zll->Branch("pvID_1", &m_pvID_1_, "pvID_1/I");
  m_tree_Zll->Branch("pvID_2", &m_pvID_2_, "pvID_2/I");
  m_tree_Zll->Branch("pvID_3", &m_pvID_3_, "pvID_3/I");
  m_tree_Zll->Branch("pvID_4", &m_pvID_4_, "pvID_4/I");

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
          50, 0, 10);
  sc = m_thistSvc->regHist("/AANT/Muon/chi2_fit_12", m_chi2_fit_12);
  m_chi2_fit_34 = new TH1F("chi2_fit_34", "chi2 fit for 3 and 4",
          50, 0, 10);
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

  m_v0_x = -99;
  m_v0_y = -99;
  m_v0_z = -99;
  
  n_muon = 0;
  mu_author_->clear();
  mu_pt_->clear();
  mu_track_pt_->clear();
  mu_eta_->clear();
  mu_phi_->clear();
  mu_e_->clear();
  mu_charge_->clear();
  mu_type_->clear();
  mu_d0_->clear();
  mu_z0_sintheta_->clear();
  mu_d0_sig_->clear();
  mu_eloss_->clear();
  mu_etcone30_->clear();
  mu_ptcone30_->clear();


  m_upsilon_ = 0;
  m_4l_ = -9999;
  vtx4l_chi2ndf_ = -9999;
  m34_ = -9999;
  same_vertex_ = false;
  m_4l_fitted_ = -9999;


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
      m_v0_x = primaryVertex.position().x();
      m_v0_y = primaryVertex.position().y();
      m_v0_z = primaryVertex.position().z();
    }
  }

  // do the Z->mm reconstruction on AOD
  sc = zmm_on_aod();
  if ( sc.isFailure() ) {
    // ATH_MSG_WARNING("Z->mm reconstruction on AOD failed");
    return StatusCode::SUCCESS;
  } else {
      m_tree_Zll->Fill();
  }
  //
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

  MuonVect*  good_muons = new MuonVect();
  for (; muonItr != muonItrE; ++muonItr) {
     imuon++;
     TLorentzVector tlv = this->getLorentzVector( (**muonItr) );

     float charge = (*muonItr)->charge();

     // apply selections
     if( !this->passMuon(**muonItr) ) continue;
     m_muons_cutFlow->Fill(3);
     if( !(*muonItr)->isCombinedMuon() and\
         !(*muonItr)->isSegmentTaggedMuon()){
        continue;
     }
     m_muons_cutFlow->Fill(4);
     // charge of ID track and muon
     const Rec::TrackParticle* id_track = (*muonItr)->inDetTrackParticle();
     bool consitent_charge = false;
     if(id_track){
        consitent_charge = id_track->charge() == charge;
     }
     if(!consitent_charge) continue;
     m_muons_cutFlow->Fill(5);

     good_muons->push_back( (*muonItr) );
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
             float z0_sintheta = aMeasPer->parameters()[Trk::z0] * tlv.Theta();
             mu_z0_sintheta_->push_back(z0_sintheta);
             const Trk::ErrorMatrix& theErrorMatrix = aMeasPer->localErrorMatrix();
             float d0_error = (theErrorMatrix.covariance())[Trk::d0][Trk::d0];
             mu_d0_sig_->push_back( d0/sqrt(d0_error) );
         }
     }
     float track_pt = -999;
     if(id_track) track_pt = id_track->pt();
     mu_track_pt_->push_back(track_pt);

     mu_eloss_->push_back( (*muonItr)->energyLoss().first );
     mu_etcone30_->push_back( (*muonItr)->parameter(MuonParameters::etcone30) );
     mu_ptcone30_->push_back( (*muonItr)->parameter(MuonParameters::ptcone30) );

     m_aod_muon_charge->Fill( (*muonItr)->charge() );
     m_aod_muon_pt->Fill( (*muonItr)->pt(), 1.);
     m_aod_muon_eta->Fill( (*muonItr)->eta(), 1.);
     m_aod_muon_chi2->Fill( (*muonItr)->matchChi2OverDoF(), 1.);
     ATH_MSG_DEBUG("imuon=" <<imuon<< ", AthenaBarCode=" << (*muonItr)->getAthenaBarCode());

  }
  // remove event with less than four muons
  if(n_muon < 4) return StatusCode::FAILURE;
  m_cutFlow->Fill(2);

  this->buildFourMuons( *good_muons );

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

bool VFitZmmOnAOD::passMuon(const Analysis::Muon& muon)
{
    m_muons_cutFlow->Fill(1);
    bool result = muon.pt() > this->m_etMuonCut && fabs(muon.eta()) < this->m_etaMuonCut;
    if(!result) return false;

    m_muons_cutFlow->Fill(2);
    // add track quality cuts
    const Rec::TrackParticle* id_track = muon.inDetTrackParticle();
    if(!id_track) return false;

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

    // passBLayer = (!expectedBLayerHits || n_blayerHits > 0); // donot apply blayer fornow
    result &= (n_pixHits + n_pixelDeadSensor) > 1;
    result &= (n_sctHits + n_sctDeadSenesor) >= 6;
    result &= (n_pixelHoles + n_sctHoles) < 3;

    int n_trt = n_trtHits + n_trtOutlier;
    float eta = muon.eta();
    if (abs(eta) < 1.9) {
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

float VFitZmmOnAOD::VkVrtFit(const MuonVect& muons)
{
    std::vector<const Rec::TrackParticle *> myTracks;
    std::vector<const Trk::TrackParticleBase *> myTrackBases;
    std::vector<double> myMuonMasses;

    Trk::RecVertex primaryVtx = *primVtx;

    Hep3Vector  primVtxPos(primVtx->position().x(), primVtx->position().y(), primVtx->position().z());
    for(MuonVect::const_iterator mu_itr = muons.begin();
            mu_itr != muons.end(); ++ mu_itr) {
        const Rec::TrackParticle* tp = (*mu_itr)->track();
        if(tp) {
            myTracks.push_back(tp);
            myTrackBases.push_back(tp);
            myMuonMasses.push_back(mMuon);
        }
    }

    // vertex contrained to the primary vertex
    m_VKVrtFitter->setVertexForConstraint(primaryVtx);
    m_VKVrtFitter->setMassInputParticles(myMuonMasses);

    // define variables returned by the vertex fit
    Hep3Vector retVtxPos;
    HepLorentzVector   retMom4;
    long int retQ = 0;
    std::vector<double> errMatrix,chi2PerTrk; // Fit error matrix and chi2 per track
    std::vector< std::vector<double> > trkAtVrt; // "True" tracks passing through vertex [phi, theta, q/pt*sin(theta)]
    double  fitChi2_vk=0.;

    StatusCode sc = m_VKVrtFitter->VKalVrtFit(myTrackBases, retVtxPos,retMom4,retQ,errMatrix,chi2PerTrk,trkAtVrt,fitChi2_vk);
    int dof = m_VKVrtFitter->VKalGetNDOF();

    if (sc.isSuccess()) {
        return fitChi2_vk/dof;
    } else {
        return -1;
    }
}

float VFitZmmOnAOD::VkVrtFit(const Analysis::Muon& muon1, const Analysis::Muon& muon2)
{
    MuonVect* muons = new MuonVect;
    muons->push_back(&muon1);
    muons->push_back(&muon2);
    float results = this->VkVrtFit(*muons);
    delete muons;

    return results;
}

void VFitZmmOnAOD::buildFourMuons(const MuonVect& muons)
{
    // ATH_MSG_INFO("Building four muons");

    // require the charge of muons to be: - + + -
    bool has_neutral_muons = false;
    bool has_two_onia = false;
    float min_chi2_ndf_4muons = 9E6;

    int id_1 = -1;
    for( MuonVect::const_iterator muonItr1  = muons.begin();
            muonItr1 != muons.end(); ++muonItr1) {
        id_1 ++;
        float mu_charge_1 = (*muonItr1)->charge();
        if(mu_charge_1 > 0) continue;

        int id_2 = -1;
        for( MuonVect::const_iterator muonItr2  = muons.begin();
                muonItr2 != muons.end(); ++muonItr2) {
            id_2 ++;
            if(muonItr2 == muonItr1) continue;
            float mu_charge_2 = (*muonItr2)->charge();
            if(mu_charge_2 < 0) continue;

            int id_3 = -1;
            for( MuonVect::const_iterator muonItr3  = muons.begin();
                    muonItr3 != muons.end(); ++muonItr3) {

                id_3 ++;
                if(muonItr3 == muonItr1 || muonItr3 == muonItr2) continue;
                if( (*muonItr3)->charge() < 0) continue;

                int id_4 = -1;
                for( MuonVect::const_iterator muonItr4  = muons.begin();
                        muonItr4 != muons.end(); ++muonItr4) {

                    id_4 ++;
                    if(muonItr4 == muonItr1 || muonItr4 == muonItr2 || muonItr4 == muonItr3) continue;
                    if( (*muonItr4)->charge() > 0) continue;
                    // so far we have four neutral tracks

                    // require at least three combined muons
                    int n_combined = 0;
                    if( (*muonItr1)->isCombinedMuon() ) n_combined ++;
                    if( (*muonItr2)->isCombinedMuon() ) n_combined ++;
                    if( (*muonItr3)->isCombinedMuon() ) n_combined ++;
                    if( (*muonItr4)->isCombinedMuon() ) n_combined ++;
                    if(n_combined < 3) continue;
                    has_neutral_muons = true;
                    n_combined_muons_ = n_combined;

                    // start to build Onia(1,2), Onia(3,4), Onia(1,3) and Onia(2, 4)
                    float m_light_can1 = -1;
                    float m_light_can2 = -1;
                    float m_heavy_can1 = -1;
                    float m_heavy_can2 = -1;
                    bool pass_comb1 = this->passOniaCuts( (**muonItr1), (**muonItr2), (**muonItr3), (**muonItr4), m_light_can1, m_heavy_can1);
                    bool pass_comb2 = this->passOniaCuts( (**muonItr1), (**muonItr3), (**muonItr2), (**muonItr4), m_light_can2, m_heavy_can2);
                    if(!pass_comb1 && !pass_comb2) continue;

                    TLorentzVector tlv_1 = this->getLorentzVector( (**muonItr1) );
                    TLorentzVector tlv_2 = this->getLorentzVector( (**muonItr2) );
                    TLorentzVector tlv_3 = this->getLorentzVector( (**muonItr3) );
                    TLorentzVector tlv_4 = this->getLorentzVector( (**muonItr4) );
                    m_4l_ = (tlv_1 + tlv_2 + tlv_3 + tlv_4).M();

                    has_two_onia = true;

                    MuonVect* muons_can = new MuonVect;
                    muons_can->push_back( (*muonItr1) );
                    muons_can->push_back( (*muonItr2) );
                    muons_can->push_back( (*muonItr3) );
                    muons_can->push_back( (*muonItr4) );

                    float chi2_ndf_4muons = this->VkVrtFit( *muons_can );
                    if(chi2_ndf_4muons >= min_chi2_ndf_4muons) {
                        continue;
                    } else {
                        vtx4l_chi2ndf_ = min_chi2_ndf_4muons = chi2_ndf_4muons;
                    }

                    if(pass_comb1 && pass_comb2) {
                        m_upsilon_ = (m_heavy_can1 > m_heavy_can2)?m_heavy_can1:m_heavy_can2;
                        m34_ = (m_heavy_can1 > m_heavy_can2)?m_light_can1:m_light_can2;
                        m_index_1_ = (m_heavy_can1 > m_heavy_can2)?id_1:id_3;
                        m_index_2_ = (m_heavy_can1 > m_heavy_can2)?id_2:id_4;
                        m_index_3_ = (m_heavy_can1 > m_heavy_can2)?id_3:id_1;
                        m_index_4_ = (m_heavy_can1 > m_heavy_can2)?id_4:id_2;
                    }else if(pass_comb1) {
                        m_upsilon_ = m_heavy_can1;
                        m34_ = m_light_can1;
                        m_index_1_ = id_1;
                        m_index_2_ = id_2;
                        m_index_3_ = id_3;
                        m_index_4_ = id_4;
                    }else if(pass_comb2) {
                        m_upsilon_ = m_heavy_can2;
                        m34_ = m_light_can2;
                        m_index_1_ = id_3;
                        m_index_2_ = id_4;
                        m_index_3_ = id_1;
                        m_index_4_ = id_2;
                    }else {
                        ;
                    }
                }
            }

        }
    }
    vtx4l_chi2ndf_ = min_chi2_ndf_4muons;
    if(has_neutral_muons) m_cutFlow->Fill(3);
    if(has_two_onia) m_cutFlow->Fill(4);
    if(m_index_1_ >= 0){ 
        m_pvID_1_ = this->matchPV( *(muons.at(m_index_1_)) );
        m_pvID_2_ = this->matchPV( *(muons.at(m_index_2_)) );
        m_pvID_3_ = this->matchPV( *(muons.at(m_index_3_)) );
        m_pvID_4_ = this->matchPV( *(muons.at(m_index_4_)) );
    }
}

bool VFitZmmOnAOD::passOniaCuts(const Analysis::Muon& muon1, const Analysis::Muon& muon2, const Analysis::Muon& muon3, const Analysis::Muon& muon4,
        float& m_light, float& m_heavy)
{
    float chi2_ndf_1 = VkVrtFit( muon1, muon2 );
    float chi2_ndf_2 = VkVrtFit( muon3, muon4 );
    m_chi2_fit_12->Fill(chi2_ndf_1);
    m_chi2_fit_34->Fill(chi2_ndf_2);
    if(chi2_ndf_1 >= 3 || chi2_ndf_2 >=3 ) return false;

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
    // require  PV association.
    int index_pv = -1;
    int res = -1;
    float min_dr = 1E6;
    for(VxContainer::const_iterator vtxItr = vxCont->begin();
            vtxItr != vxCont->end(); ++ vtxItr)
    {
        index_pv ++;
        const std::vector<Trk::VxTrackAtVertex*>* vtxTracks = (*vtxItr)->vxTrackAtVertex();
        unsigned int nTracks = vtxTracks->size();

        const Rec::TrackParticle* track = muon.inDetTrackParticle();
        // const Rec::TrackParticle* track = muon.track();
        const Trk::MeasuredPerigee* aMeasPer = track->measuredPerigee();
        
        float min_dr_track = 9E6;
        for(unsigned int j = 0; j < nTracks; j++){
            // const Trk::TrackParameters* para = dynamic_cast<const Trk::TrackParameters*>(vtxTracks->at(j)->perigeeAtVertex());
            const Trk::ParametersBase* para = vtxTracks->at(j)->perigeeAtVertex();
            // const Trk::ParametersBase* para = vtxTracks->at(j)->origTrack();
            if (para && aMeasPer){
                //float d0_diff = para->parameters()[Trk::d0] - aMeasPer->parameters()[Trk::d0];
                //float z0_diff = para->parameters()[Trk::z0] - aMeasPer->parameters()[Trk::z0];
                float d0_diff = para->position()[Trk::d0] - aMeasPer->parameters()[Trk::d0];
                float z0_diff = para->position()[Trk::z0] - aMeasPer->parameters()[Trk::z0];
                float dr = sqrt(d0_diff*d0_diff + z0_diff*z0_diff);
                if( dr < min_dr_track){
                    min_dr_track = dr;
                }
            }
        }
        if(min_dr_track < min_dr){
            min_dr = min_dr_track;
            res = index_pv;
        }
    }
    return res;
}
