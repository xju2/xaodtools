// System include(s):
#include <memory>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>

// ROOT include(s):
#include <TFile.h>
#include <TError.h>
#include <TString.h>
#include <TChain.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TVector2.h>

// Infrastructure include(s):
#ifdef ROOTCORE
#   include "xAODRootAccess/Init.h"
#   include "xAODRootAccess/TEvent.h"
#endif // ROOTCORE

#include "AsgTools/AsgTool.h"
#include "AsgTools/ToolHandle.h"

// EDM include(s):
#include "xAODEventInfo/EventInfo.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODCaloEvent/CaloCluster.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthEvent.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODBTaggingEfficiency/BTaggingEfficiencyTool.h"
#include "xAODBTaggingEfficiency/BTaggingSelectionTool.h"

#include "xAODBase/IParticleHelpers.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"

// Local include(s):
#include "CPAnalysisExamples/errorcheck.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"

// Other includes
#include "PATInterfaces/SystematicVariation.h"
#include "PATInterfaces/SystematicRegistry.h"
#include "PATInterfaces/SystematicCode.h"
#include "boost/unordered_map.hpp"

#include "MyXAODTools/EventCounter.h"
#include "MyXAODTools/HZZ4lHelper.h"

using namespace std;

bool descend_on_pt(xAOD::IParticle* p1, xAOD::IParticle* p2){
    return p1->pt() > p2->pt();
}

int main( int argc, char* argv[] ) {

    if(argc > 1 && string(argv[1]) == "help"){
        cout << argv[0] << " toberun.txt num_evts" << endl;
        return 0;
    }
	StatusCode::enableFailure();
	CP::SystematicCode::enableFailure();
	CP::CorrectionCode::enableFailure();

	// The application's name:
	const char* APP_NAME = argv[ 0 ];

	// Initialise the application:
	CHECK( xAOD::Init( APP_NAME ) );
    string inputFileName = "toberun.txt";
    if(argc > 1){
        inputFileName = string(argv[1]);
    }

	TChain* fc = EventCounter::getChain(inputFileName.c_str(), "CollectionTree"); 

	// Create a TEvent object:
	xAOD::TEvent event( xAOD::TEvent::kClassAccess );
	CHECK( event.readFrom( fc ) );
	Info( APP_NAME, "Number of events in the chain: %i",
			static_cast< int >( event.getEntries() ) );


	// Decide how many events to run over:
	Long64_t entries = event.getEntries();
	if( argc > 2 ) {
		const Long64_t e = atoll( argv[ 2 ] );
		if( e>0 && e < entries ) {
			entries = e;
		}
	}
    Info(APP_NAME, "Will process: %d", (int)entries);

	int RunNumber;
	int EventNumber;
	int mc_channel_number;
	double MCWeight;
    double mc_weight_up;
    double mc_weight_down;
    vector<double>* mc_weights = new vector<double>();

	TFile *fOutputFile = new TFile( "reduced_ntuple.root", "recreate" );
	TTree MyTree( "physics", "physics" );

	MyTree.Branch("RunNumber", &RunNumber, "RunNumber/I");
	MyTree.Branch("EventNumber", &EventNumber, "EventNumber/I");
	MyTree.Branch("mc_channel_number", &mc_channel_number, "mc_channel_number/I");

    unique_ptr<HZZ4lHelper> h4l_helper(new HZZ4lHelper());
    h4l_helper->MakeOutputTree(MyTree);
    h4l_helper->MakeTruthTree(MyTree);

    int type;
    MyTree.Branch("truth_event_type", &type, "truth_event_type/I");
    int truth_type;
    MyTree.Branch("type_truth", &truth_type, "type_truth/I");

	MyTree.Branch("MCWeight", &MCWeight, "MCWeight/D");
	MyTree.Branch("MCWeightUp", &mc_weight_up, "MCWeightUp/D");
	MyTree.Branch("MCWeightDown", &mc_weight_down, "MCWeightDown/D");
	MyTree.Branch("MCWeights", &mc_weights);

    // add jet information
    int n_jets_30 = 0;
    double dijet_invmass = -999;
    double dijet_deltaeta = -999;
	MyTree.Branch("n_jets", &n_jets_30, "n_jets/I");
	MyTree.Branch("dijet_deltaeta_fid", &dijet_deltaeta, "dijet_deltaeta_fid/D");
	MyTree.Branch("dijet_invmass_fid", &dijet_invmass, "dijet_invmass_fid/D");

    int pass_fid = -1;
    MyTree.Branch("pass_fid_cut", &pass_fid, "pass_fid_cut/I");
    int pass_fid_truth = -1;
    MyTree.Branch("pass_fid_cut_truth", &pass_fid_truth, "pass_fid_cut_truth/I");

	for( Long64_t entry = 0; entry < entries; ++entry ) {

        mc_weights->clear();
        h4l_helper->Clear();
        dijet_invmass = -999;
        dijet_deltaeta = -999;
        pass_fid = -1;
        type = -1;
        truth_type = -1;

		// Tell the object which entry to look at:
		event.getEntry( entry );

		const xAOD::EventInfo* ei = 0;
		CHECK( event.retrieve( ei, "EventInfo" ) );

		RunNumber = ei->runNumber();
		EventNumber = ei->eventNumber();
		mc_channel_number = ei-> mcChannelNumber();

		MCWeight = ei->mcEventWeight();

        const vector<float>& weights = ei->mcEventWeights();
        if(weights.size() > 2) {
            mc_weight_up = weights.at(1);
            mc_weight_down = weights.at(2);
        }
        for(int i = 0; i < (int)weights.size(); i ++){
            mc_weights->push_back(weights.at(i));
        }

		// Get the Jets from the event:
		const xAOD::JetContainer* jets_origin = 0;
        CHECK( event.retrieve( jets_origin, "AntiKt4TruthJets" ) );
        std::pair<xAOD::JetContainer*, xAOD::ShallowAuxContainer*> jet_shallowcopy = xAOD::shallowCopyContainer(*jets_origin);
        xAOD::JetContainer* jets = jet_shallowcopy.first;
       
        // count number of jets with pt > 30 GeV and |eta| < 4.5
        n_jets_30 = 0;
        for(const auto& jet : *jets){
            if( jet->p4().Pt() > 30E3 && fabs(jet->p4().Eta()) < 4.5 ) {
                n_jets_30 ++;
            }
        }

        // calculate dijet invariant mass and delta-eta
        if (jets->size() >= 2){
            sort(jets->begin(), jets->end(), descend_on_pt);
            if(jets->at(1)->p4().Pt()/1E3 > 30
               && fabs(jets->at(0)->p4().Eta()) < 4.5
               && fabs(jets->at(1)->p4().Eta()) < 4.5
               ){
                dijet_invmass = (jets->at(0)->p4() + jets->at(1)->p4()).M();
                dijet_deltaeta = fabs(jets->at(0)->eta() - jets->at(1)->eta());
            }
        }

        // get truth info
       const xAOD::TruthParticleContainer* particles(0);
       if( event.retrieve(particles, "TruthParticles").isSuccess() )
       {
           // get truth information
           h4l_helper->GetTruthInfo( *particles );
           truth_type = h4l_helper->getTruthType();
       } else {
           cout << "No truth particles" << endl;
       }

        ////get electrons and muons, then to mimic our fiducial selection
        const xAOD::TruthParticleContainer* electrons(0);
        const xAOD::TruthParticleContainer* muons(0);
        CHECK( event.retrieve(muons, "TruthMuons") );
        if( ! event.retrieve(electrons,"TruthElectrons").isSuccess() ||
            ! event.retrieve(muons, "TruthMuons").isSuccess() )
        {
            if (entry == 0){
                cout << "Cannot find TruthMuons or TruthElectrons" << endl;
            }
            continue;
        }
        xAOD::TruthParticleContainer::const_iterator ele_itr = electrons->begin();
        xAOD::TruthParticleContainer::const_iterator ele_end = electrons->end();

        vector<Candidate*>* ele_4vec = new vector<Candidate*>();
        int ele_index = 0;
        for(; ele_itr != ele_end; ele_itr ++){
            if( (*ele_itr)->status() == 1 && (*ele_itr)->absPdgId() == 11)
            {
                TLorentzVector ele;
                ele.SetPxPyPzE((*ele_itr)->px(),(*ele_itr)->py(),
                        (*ele_itr)->pz(), (*ele_itr)->e());
                if(ele.Pt() > 7e3 && ele.Eta() < 2.47){
                    float charge = (*ele_itr)->pdgId()/11.;
                    Candidate* can = new Candidate(ele_index, ele);
                    can->setCharge(charge);
                    ele_4vec ->push_back(can);
                }
            }
            ele_index ++;
        }
        // Info(APP_NAME, "Number of electrons: %d", (int)ele_4vec->size());

        xAOD::TruthParticleContainer::const_iterator muon_itr = muons->begin();
        xAOD::TruthParticleContainer::const_iterator muon_end = muons->end();
        // cout << "total muons: " << muons->size() << endl;
        vector<Candidate*>* muon_4vec = new vector<Candidate*>();
        int index = 0;
        for(; muon_itr != muon_end; muon_itr++){
            if( (*muon_itr)->status() == 1 && (*muon_itr)->absPdgId() == 13)
            {
                TLorentzVector mu ;
                mu.SetPxPyPzE((*muon_itr)->px(), (*muon_itr)->py(),
                        (*muon_itr)->pz(), (*muon_itr)->e());
                if(mu.Pt() > 5e3 && mu.Eta() < 2.7){
                    float charge = (*muon_itr)->pdgId()/13.;
                    Candidate* can = new Candidate(index, mu);
                    can->setCharge(charge);
                    muon_4vec->push_back(can);
                }
            }
            index ++;
        }
        // Info(APP_NAME, "Number of muons: %d", (int)muon_4vec->size());

        if(ele_4vec->size() >= 4 && h4l_helper->Is_4mu4e(ele_4vec)){
            type = 1; //4electrons
        }else if(muon_4vec->size() >= 4 && h4l_helper->Is_4mu4e(muon_4vec)){
            type = 0; //4muons
        }else if (ele_4vec->size() >=2 && muon_4vec->size() >= 2){ 
            h4l_helper->Is_2e2mu(ele_4vec, muon_4vec, type);
        }else{
            type = -1;
        }
        h4l_helper->setType(type);

        pass_fid = h4l_helper->passFiducial();
        pass_fid_truth = h4l_helper->passFiducialTruth();

        MyTree.Fill();
        for(int i =0; i < (int)ele_4vec->size(); i++){
            delete ele_4vec->at(i);
        }
        delete ele_4vec;
        for(int i = 0; i < (int)muon_4vec->size(); i++){
            delete muon_4vec->at(i);
        }
        delete muon_4vec;
    }
    delete mc_weights;

    fOutputFile->cd();
    MyTree.Write();
    fOutputFile->Close();

    Info( APP_NAME, "Successfully finished analysis; Exitting..." );
    return 1;
}
