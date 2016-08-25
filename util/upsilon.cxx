// System include(s):
#include <memory>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <typeinfo>
#include <math.h>

// ROOT include(s):
#include <TFile.h>
#include <TError.h>
#include <TString.h>
#include <TChain.h>
#include <TH1F.h>

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
#include "xAODEgamma/ElectronAuxContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"

#include "xAODCore/ShallowCopy.h"

// Local include(s):
#include "CPAnalysisExamples/errorcheck.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"

// Other includes
// #include "PATInterfaces/SystematicList.h"
#include "PATInterfaces/SystematicVariation.h"
#include "PATInterfaces/SystematicRegistry.h"
#include "PATInterfaces/SystematicCode.h"

/** Vertexing
#include "TrkVertexFitterInterfaces/IVertexFitter.h"
#include "TrkVKalVrtFitter/TrkVKalVrtFitter.h"
#include "CLHEP/Units/SystemOfUnits.h"
***/

#include "MyXAODTools/CPToolsHelper.h"
#include "MyXAODTools/TrackBranch.h"
#include "MyXAODTools/EventInfoCreator.h"
#include "MyXAODTools/ElectronBranch.h"
#include "MyXAODTools/MuonBranch.h"
#include "MyXAODTools/UpsilonBranch.h"

#include <TError.h>
using namespace std;

static SG::AuxElement::Decorator<double> dec_dphi_MET("dphi_MET");
static SG::AuxElement::Decorator<char> dec_baseline("baseline");
static SG::AuxElement::Decorator<char> dec_signal("signal");
static SG::AuxElement::Decorator<char> dec_passOR("passOR");
static SG::AuxElement::Decorator<char> dec_bad("bad");
static SG::AuxElement::Decorator<char> dec_bjet("bjet");
static SG::AuxElement::Decorator<char> dec_cosmic("cosmic");
static SG::AuxElement::Decorator<double> dec_effscalefact("effscalefact");
static SG::AuxElement::Decorator<char> dec_isol("isol");
static SG::AuxElement::Decorator<char> dec_tightBad("tightBad");
static SG::AuxElement::Decorator<int> dec_muonIndex("BPHY4MuonIndex");

bool descend_on_pt(xAOD::IParticle* p1, xAOD::IParticle* p2){
    return p1->pt() > p2->pt();
}

int main( int argc, char* argv[] ) 
{
    // gErrorIgnoreLevel = kError;

    if (argc > 1 && string(argv[1]) == "help") {
        cout << argv[0] << " toberun.txt number_evts isData=1 debug=0" << endl;
        exit(0);
    }
    // The application's name:
    const char* APP_NAME = argv[ 0 ];
    string inputFileName = "toberun.txt";
    if( argc > 1 ){
        inputFileName = string( argv[1] ); 
    }

    // Initialise the application:
    CHECK( xAOD::Init( APP_NAME ) );

    ifstream input_file( inputFileName.c_str() );
    Info( APP_NAME, "Input file name: %s", inputFileName.c_str() );
    TString name_file;
    TChain* fc = new TChain("CollectionTree");
    uint64_t total_evts_pro = 0;
    double sum_of_evt_w = 0;
    double sum_of_evt_w_sq = 0;
    while ( input_file>>name_file ) {
        uint64_t n_events_process = 0;
        double sum_of_evt_weights = 0;
        double sum_of_evt_weight_sqd = 0;
        CHECK(CPToolsHelper::GetProcessEventsInfo(name_file.Data(), 
                    n_events_process,
                    sum_of_evt_weights,
                    sum_of_evt_weight_sqd));
        total_evts_pro += n_events_process;
        sum_of_evt_w += sum_of_evt_weights;
        sum_of_evt_w_sq += sum_of_evt_weight_sqd;

        fc->Add(name_file);
    }

    Info(APP_NAME, "Total events: %lu", total_evts_pro);
    Info(APP_NAME, "Sum of evt weights: %f", sum_of_evt_w);
    Info(APP_NAME, "Sum of evt weights sq: %f", sum_of_evt_w_sq);

    // Create a TEvent object:
    xAOD::TEvent event( xAOD::TEvent::kClassAccess );
    CHECK( event.readFrom( fc ) );
    Info( APP_NAME, "Number of events in the chain: %i",
            static_cast< int >( event.getEntries() ) );

    // StatusCode::enableFailure();
    // CP::SystematicCode::enableFailure();
    // CP::CorrectionCode::enableFailure();
    xAOD::TStore store;

    // Decide how many events to run over:
    Long64_t entries = event.getEntries();
    if( argc > 2 ) {
        const Long64_t e = atoll( argv[ 2 ] );
        if( e>0 && e < entries ) {
            entries = e;
        }
    }

    int isData = 0;
    bool do_debug = false;

    for (int i= 3 ; i<argc ; i++) {
        const char* key = strtok(argv[i],"=") ;
        const char* val = strtok(0," ") ;
        Info( APP_NAME,  "processing key %s  with value %s", key, val );
        if (strcmp(key,"isData")==0) isData = atoi(val);
        if (strcmp(key,"debug")==0) do_debug = (bool)atoi(val);
    }

    Info( APP_NAME, "Number of events to process: %i", static_cast<int>( entries ) );
    if (isData) {
        cout <<"You are running data, congratuations" << endl;
    }

    ////////////////////////////////////////////////////////
    // create SUSYTools and config it
    ST::SUSYObjDef_xAOD objTool("SUSYObjDef_Upsilon");
    std::cout << " ABOUT TO INITIALIZE SUSYTOOLS " << std::endl;
    if(do_debug) objTool.msg().setLevel(MSG::VERBOSE);
    else objTool.msg().setLevel(MSG::ERROR);

    // Configure the SUSYObjDef instance
    ST::ISUSYObjDef_xAODTool::DataSource data_source = isData ? ST::ISUSYObjDef_xAODTool::Data : ST::ISUSYObjDef_xAODTool::FullSim;
    CHECK( objTool.setProperty("DataSource", data_source) );

    // general configuration
    string maindir(getenv("ROOTCOREBIN"));
    string config_file = Form("%s/data/MyXAODTools/upsilon.conf", maindir.c_str());
    CHECK( objTool.setProperty("ConfigFile", config_file) );

    // pileup reweight
    vector<string> prw_conf;
    prw_conf.push_back(maindir+"/data/MyXAODTools/mc15c.prw.root");
    CHECK( objTool.setProperty("PRWConfigFiles", prw_conf) );

    vector<string> prw_lumicalc;
    prw_lumicalc.push_back(maindir+"/data/MyXAODTools/ilumicalc_histograms_None_276262-284484_final_20.7.root");
    prw_lumicalc.push_back(maindir+"/data/MyXAODTools/ilumicalc_histograms_None_297730-303892.root");
    CHECK( objTool.setProperty("PRWLumiCalcFiles", prw_lumicalc) );

    if( objTool.initialize() != StatusCode::SUCCESS){
        Error( APP_NAME, "Cannot intialize SUSYObjDef_xAOD..." );
        Error( APP_NAME, "Exiting... " );
        exit(-1);
    }else{
        Info( APP_NAME, "SUSYObjDef_xAOD initialized... " );
    }
    std::cout << " INITIALIZED SUSYTOOLS " << std::endl;
    // end of SUSYTools
    ////////////////////////////////////////////
    
    // vertexing stuff
    // ToolHandle < Trk::IVertexFitter > m_iVertexFitter;
    // Trk::TrkVKalVrtFitter* m_VKVFitter = new Trk::TrkVKalVrtFitter();


    TFile *fOutputFile = new TFile( "reduced_ntup.root", "recreate" );

    // GRL, JVT and other tools are moved to MyXAODTools
    CPToolsHelper* cp_tools = new CPToolsHelper();

    EventInfoCreator* event_br = new EventInfoCreator(isData);


    /*record the number of processed events*/
    /*****
    auto tree = new TTree("associate","associate");
    event_br->AttachMiniToTree(*tree);
    tree->Branch("nEventsProcessed", &total_evts_pro, "nEventsProcessed/l");
    tree->Branch("nSumEventWeights", &sum_of_evt_w, "nSumEventWeights/D");
    tree->Branch("nSumEventWeightsSquared", &sum_of_evt_w_sq, 
            "nSumEventWeightsSquared/D");

    for (Long64_t entry = 0; entry < 1; ++entry) {
        event_br->ClearBranch();

        event.getEntry( entry );
        const xAOD::EventInfo* ei = 0;
        CHECK( event.retrieve( ei, "EventInfo" ) );

        event_br->Fill(*ei);
        tree ->Fill();
    }
    fOutputFile->cd();
    tree ->Write();
    *****/

    // Fill physics tree
    string tree_name = "physics";

    TTree MyTree(tree_name.c_str(), tree_name.c_str());

    UpsilonBranch* output = new UpsilonBranch();
    MuonBranch* muon_br = new MuonBranch(); 
    ElectronBranch* el_br = new ElectronBranch();

    event_br->AttachBranchToTree(MyTree);
    output->AttachBranchToTree(MyTree);
    el_br->AttachBranchToTree(MyTree);
    muon_br->AttachBranchToTree(MyTree);
    
    TH1F* h_cutflow = new TH1F("h_cutflow", "cut flow", 101, -0.5, 100.5);

    // Set Branches
    vector<float>* br_muon_energyloss = NULL;
    vector<float>* br_muon_etcone30 = NULL;
    vector<float>* br_muon_ptvarcone30 = NULL;
    vector<string>* br_bphy4quads_combinationCode = NULL;
    fc->SetBranchAddress("MuonsAuxDyn.EnergyLoss", &br_muon_energyloss);
    fc->SetBranchAddress("MuonsAuxDyn.etcone30", &br_muon_etcone30);
    fc->SetBranchAddress("MuonsAuxDyn.ptvarcone30", &br_muon_ptvarcone30);

    // fc->SetBranchAddress("BPHY4QuadsAuxDyn.CombinationCode", &br_bphy4quads_combinationCode);
    
    const float UPSILON_MASS = 9.46E3;
    for( Long64_t entry = 0; entry < entries; ++entry ) 
    {
        output->ClearBranch();
        event_br->ClearBranch();
        el_br->ClearBranch();
        muon_br->ClearBranch();

        h_cutflow->Fill(0);

        // Tell the object which entry to look at:
        fc->GetEntry( entry );
        event.getEntry( entry );

        const xAOD::EventInfo* ei = 0;
        CHECK( event.retrieve( ei, "EventInfo" ) );
        CHECK(objTool.ApplyPRWTool());

        event_br->Fill(*ei);

        // if(entry == 0){
        //     cout << "Dumping trigger info" << endl;
        //     auto chainGroup = objTool.GetTrigChainGroup("HLT_.*");
        //     for (auto& trig : chainGroup->getListOfTriggers())
        //         cout << " " << trig << endl;
        // }


        if(! cp_tools->PassGRL(*ei)) continue;
        h_cutflow->Fill(1);
        if(! cp_tools->PassEventCleaning(*ei)) continue;
        h_cutflow->Fill(2);
        // Info(APP_NAME, "In event: %d", (int)ei->eventNumber());

        // if(ei->eventNumber() != 620) continue;
        // cout<< "reading 620!"<<endl;

        for(auto& kv : output->trigger_map_)
        {
            if(objTool.IsTrigPassed(kv.first.c_str())) 
            {
                // Info(APP_NAME, "Trigger %s is fired", kv.first.c_str());
                output->pass_trigger_ = kv.second = true;
            }
        }

        //////////////////
        //primary vertex
        //////////////////
        const xAOD::VertexContainer* vertice = 0;
        CHECK( event.retrieve(vertice, "PrimaryVertices") );
        if(! cp_tools->HasPrimaryVertex(*vertice, 1)) continue;
        h_cutflow->Fill(3);

        const xAOD::Vertex* pv = CPToolsHelper::GetPrimVtx(*vertice);

        /*get physics objects*/
        // Electrons
        /****
        xAOD::ElectronContainer* electrons_copy = NULL;
        xAOD::ShallowAuxContainer* electrons_copyaux = NULL;
        CHECK( objTool.GetElectrons(electrons_copy, electrons_copyaux, true) );
        */

        // Muons
        xAOD::MuonContainer* muons_copy = NULL;
        xAOD::ShallowAuxContainer* muons_copyaux = NULL;
        CHECK( objTool.GetMuons(muons_copy, muons_copyaux, true) );

        // Jets
        /***
        xAOD::JetContainer* jets_copy = NULL;
        xAOD::ShallowAuxContainer* jets_copyaux = NULL;
        CHECK( objTool.GetJets(jets_copy, jets_copyaux, true) );
        ***/
        // no jet cleaning
        // no overlap removal
        /**
        CHECK( objTool.OverlapRemoval(electrons_copy, muons_copy,
                    jets_copy) );
         ***/

        //////////////////////
        // Electron's selection
        ///////////////////////
        /**
        int n_ele = 0;
        for(auto el = electrons_copy->begin(); el != electrons_copy->end(); el++)
        {
            if( (bool) dec_signal(**el) && (bool) dec_passOR(**el) ){
                n_ele ++;
                el_br->Fill(**el);
            }
        }
        ***/

        ///////////////////
        // Muon's Selection
        ///////////////////
        int n_muon = 0;
        int imuon = 0;
        if(do_debug)
            cout <<"processing: "<< ei->runNumber() << " " << ei->eventNumber() << endl;
        for(auto mu_itr = muons_copy->begin();
                mu_itr != muons_copy->end(); ++mu_itr)
        {
            // if( (bool) dec_signal(**mu_itr) && (bool) dec_passOR(**mu_itr))
            if( (bool) dec_baseline(**mu_itr) )
            {
                n_muon ++;
                muon_br->Fill(**mu_itr, ei, pv);
                int muIndex = (*mu_itr)->auxdataConst<int>("BPHY4MuonIndex");
                
                if(do_debug){
                    cout << "Type: " << (*mu_itr)->muonType() << endl;
                    cout << "Index: " << muIndex << endl;
                    cout << "pT: " << (*mu_itr)->p4().Pt()/1E3 << endl;
                    cout << "Energy loss: " << br_muon_energyloss->at(imuon) << endl;
                    cout << "ETcone30: " << br_muon_etcone30->at(imuon) << endl;
                    cout << "pTvarcone30: " << br_muon_ptvarcone30->at(imuon) << endl;
                }
                muon_br->eloss_->push_back(br_muon_energyloss->at(imuon));
                muon_br->etcone30_->push_back(br_muon_etcone30->at(imuon));
                muon_br->ptvarcone30_->push_back(br_muon_ptvarcone30->at(imuon));
            }
            imuon ++;
        }
        // reject events with less than four muons!
        // otherwise too many events to handle with
        if (n_muon < 4) {
            store.clear();
            continue;
        }

        const float UPSILON_LOW = 8E3;
        const float UPSILON_HI = 12E3;
        if( n_muon >= 2 ) {
            h_cutflow->Fill(4);
            output->event_type_ = 0;
            // upsilon decays to two muons
            // loop over muons and reconstruct upsilon
            auto upsilon_mu1 = muons_copy->end();
            auto upsilon_mu2 = muons_copy->end();
            for(auto mu_itr1 = muons_copy->begin();
                    mu_itr1 != muons_copy->end(); ++mu_itr1)
            {
                // if(!dec_signal(**mu_itr1) || !dec_passOR(**mu_itr1)) continue;
                if(!dec_baseline(**mu_itr1)) continue;
                const TLorentzVector& mu_tlv_1 = (*mu_itr1)->p4();
                // if(mu_tlv_1.Pt() < 5E3) continue;

                float mu_charge_1 = (*mu_itr1)->charge();

                for(auto mu_itr2 = mu_itr1 + 1; mu_itr2 != muons_copy->end(); ++mu_itr2){
                    // if(!dec_signal(**mu_itr2) || !dec_passOR(**mu_itr2)) continue;
                    if(!dec_baseline(**mu_itr2)) continue;
                    const TLorentzVector& mu_tlv_2 = (*mu_itr2)->p4();

                    // if(mu_tlv_2.Pt() < 5E3) continue;
                    float mu_charge_2 = (*mu_itr2)->charge();
                    // opposite charge
                    if(mu_charge_2 * mu_charge_1 > 0) continue;

                    auto inv_mass = (mu_tlv_1 + mu_tlv_2).M();
                    if (inv_mass > UPSILON_LOW && inv_mass < UPSILON_HI &&
                            (fabs(inv_mass - UPSILON_MASS) < fabs(output->m_upsilon_ - UPSILON_MASS))) {
                        output->m_upsilon_ = inv_mass;
                        upsilon_mu1 = mu_itr1;
                        upsilon_mu2 = mu_itr2;
                    }
                }
            }

            if (upsilon_mu1 != muons_copy->end()) {
                // only when there's a upsilon candidate, save the information!
                if (n_muon >= 4) {
                    
                    // require the second two muons to be opposite charge
                    vector<xAOD::Muon*> muon_quad;
                    muon_quad.push_back(*upsilon_mu1);
                    muon_quad.push_back(*upsilon_mu2);

                    float mass34_can = -900E3;
                    for(auto mu_itr3 = muons_copy->begin(); mu_itr3 != muons_copy->end(); ++mu_itr3)
                    {
                        if(mu_itr3 == upsilon_mu1 || mu_itr3 == upsilon_mu2) continue;
                        if(!dec_baseline(**mu_itr3)) continue;
                        float mu_charge_3 = (*mu_itr3)->charge();

                        for(auto mu_itr4 = mu_itr3+1; mu_itr4 != muons_copy->end(); ++mu_itr4){
                            if(mu_itr4 == upsilon_mu1 || mu_itr4 == upsilon_mu2 || mu_itr4 == mu_itr3) continue;
                            if(!dec_baseline(**mu_itr4)) continue;
                            float mu_charge_4 = (*mu_itr4)->charge();

                            if(mu_charge_3 * mu_charge_4 > 0) continue;
                        
                            float m4l = (float) ((*upsilon_mu1)->p4() + (*upsilon_mu2)->p4() + (*mu_itr3)->p4() + (*mu_itr4)->p4()).M()/1E3;
                            float m34 = (float) ((*mu_itr3)->p4() + (*mu_itr4)->p4()).M()/1E3;
                            if (m34 > mass34_can) 
                            { // select the most energetic ones
                                output->m_4l_ = m4l;
                                output->m34_ = mass34_can = m34;
                                if(muon_quad.size() > 2){
                                    muon_quad[2] = *mu_itr3;
                                    muon_quad[3] = *mu_itr4;
                                } else {
                                    muon_quad.push_back(*mu_itr3);
                                    muon_quad.push_back(*mu_itr4);
                                }
                            }

                        }
                    }

                    if(muon_quad.size() == 4)
                    {
                        // vertex fitting
                        vector<const xAOD::TrackParticle*> inputTracks(0);
                        // vector<ElementLink<xAOD::TrackParticleContainer> > inputTrackLinks(0);
                        for(auto& muon : muon_quad){
                            const xAOD::TrackParticle* tp = muon->trackParticle(xAOD::Muon::InnerDetectorTrackParticle);
                            if(tp) { inputTracks.push_back(tp); }
                        }
                        unsigned int nTracks = inputTracks.size();

                        const xAOD::VertexContainer* fourMuonsVertexCont = 0;
                        if (event.contains<xAOD::VertexContainer>("BPHY4Quads")) {
                            CHECK(event.retrieve(fourMuonsVertexCont, "BPHY4Quads"));

                            const xAOD::Vertex* vert = 0;
                            output->n_bphy4_quad_ = fourMuonsVertexCont->size();

                            int iquads = 0;
                            float min_chi2 = 900E3;
                            for (const auto& v: *fourMuonsVertexCont) 
                            {
                                if(!v) continue;
                                float chi2 = (float) v->chiSquared()/v->numberDoF();
                                // if(v->nTrackParticles() != nTracks) continue;
                                // if(do_debug) cout << "combination code: " << br_bphy4quads_combinationCode->at(iquads) << endl;

                                set<const xAOD::TrackParticle*> vtxTrks;
                                for(unsigned int i = 0; i < v->nTrackParticles(); ++i) {
                                    vtxTrks.insert(v->trackParticle(i));
                                }

                                unsigned int nMatch = 0;
                                for(unsigned int i = 0; i < inputTracks.size(); i++){
                                    if(vtxTrks.find(inputTracks.at(i)) != vtxTrks.end()) {
                                        nMatch ++;
                                    }
                                }
                                if(nMatch == nTracks && chi2 < min_chi2){
                                    min_chi2 = chi2;
                                    vert = v;
                                    // break;
                                }
                            }
                            if(vert) {
                                if(do_debug) cout<<"vertex: " << vert->chiSquared() << " NDOF:" << vert->numberDoF() << endl;
                                output->vtx4l_chi2ndf_  = (float) vert->chiSquared()/vert->numberDoF();
                            }
                        }
                        // Loop over the primary vertices
                        for(const auto& v: *vertice){
                            if(!v) continue;

                            set<const xAOD::TrackParticle*> vtxTrks;
                            for(unsigned int i = 0; i < v->nTrackParticles(); ++i) {
                                vtxTrks.insert(v->trackParticle(i));
                            }

                            unsigned int nMatch = 0;
                            for(unsigned int i = 0; i < inputTracks.size(); i++){
                                if(vtxTrks.find(inputTracks.at(i)) != vtxTrks.end()) nMatch ++;
                            }
                            if(nMatch == nTracks){
                                output->same_vertex_ = true;
                                break;
                            }
                        }
                    }
                }
                h_cutflow->Fill(5);
                MyTree.Fill();
            }
        }

        // The containers created by the shallow copy are owned by you.
        // Remember to delete them
        store.clear();
    }

    fOutputFile->cd();
    h_cutflow->Write();
    MyTree.Write();

    fOutputFile->cd();
    fOutputFile->Close();
    Info( APP_NAME, "finished analysis; Cleaning..." );

    delete cp_tools;

    Info( APP_NAME, "Successfully finished analysis; Exitting..." );
    return 1;
}

