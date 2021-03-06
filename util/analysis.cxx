// System include(s):
#include <TString.h>
#include <TChain.h>
// #include <TH1F.h>

// Infrastructure include(s):
#ifdef ROOTCORE
#   include "xAODRootAccess/Init.h"
#   include "xAODRootAccess/TEvent.h"
#endif // ROOTCORE

#include "AsgTools/AsgTool.h"
#include "AsgTools/ToolHandle.h"

// EDM include(s):
#include "xAODEventInfo/EventInfo.h"
#include "CPAnalysisExamples/errorcheck.h"

#include "MyXAODTools/GammaJetAna.h"
#include "MyXAODTools/UpsilonAna.h"
#include "MyXAODTools/MonoJetAna.h"
#include "MyXAODTools/FakeMuonAna.h"

using namespace std;

int main( int argc, char* argv[] )
{
    // gErrorIgnoreLevel = kError;

    if ((argc > 1 && string(argv[1]) == "help") ||(argc < 3))
    {
        cout << argv[0] << " analysisName toberun.txt number_evts isData=1 debug=0 noGRL=0 doSmear=1 useBphy1=1" << endl;
        cout << "analysisName: gammajet, upsilon, monojet, fakemuon" << endl;
        exit(1);
    }

    // The application's name:
    const char* APP_NAME = argv[ 1 ];

    // Initialise the application:
    CHECK( xAOD::Init( APP_NAME ) );

    string inputFileName(argv[2]);
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

    xAOD::TStore store;

    TString anaName(argv[1]);
    anaName.ToLower();
    // setup the analysis.
    AnalysisBase* ana = NULL;
    if(anaName == "gammajet") {
        ana = new GammaJetAna();
    } else if(anaName == "upsilon") {
        ana = new UpsilonAna();
    } else if(anaName == "monojet") {
        ana = new MonoJetAna();
    } else if(anaName == "fakemuon") {
        ana = new FakeMuonAna();
    }else{
        Error(APP_NAME, "cannot find algorithm: %s", anaName.Data());
        return 1;
    }


    // Decide how many events to run over:
    Long64_t entries = event.getEntries();
    if( argc > 3 ) {
        const Long64_t e = atoll( argv[ 3 ] );
        if( e>0 && e < entries ) {
            entries = e;
        }
    }

    int isData = 0;
    bool do_debug = false;
    bool no_grl = false;
    bool do_smear = true;
    bool use_bphy1 = false;

    for (int i= 4 ; i<argc ; i++) {
        const char* key = strtok(argv[i],"=") ;
        const char* val = strtok(0," ") ;
        Info( APP_NAME,  "processing key %s  with value %s", key, val );
        if (strcmp(key,"isData")==0) isData = atoi(val);
        if (strcmp(key,"debug")==0) do_debug = (bool)atoi(val);
        if (strcmp(key,"noGRL")==0) no_grl = (bool)atoi(val);
        if (strcmp(key,"doSmear")==0) do_smear = (bool)atoi(val);
        if (strcmp(key,"useBphy1")==0) use_bphy1 = (bool)atoi(val);
    }

    // setup each analysis!
    if (no_grl){
        ana->setGRLTag(false);
    }

    // smearing for monojet only
    if( anaName == "monojet"){
        MonoJetAna* monojet_ana = dynamic_cast<MonoJetAna*>(ana);
        monojet_ana->setSmear(do_smear);
    }
    if( anaName == "upsilon" && use_bphy1 ){
        UpsilonAna* upsilon_ana = dynamic_cast<UpsilonAna*>(ana);
        upsilon_ana->UseBPHY1();
    }

    Info( APP_NAME, "Number of events to process: %i", static_cast<int>( entries ) );
    if (isData) {
        cout <<"You are running data, congratuations" << endl;
    }
    event.getEntry(0);

    ana->SetEvent(&event); // don't change the order with following commands

    // save total number of events first, which also tells if it's Data/MC
    // SUSYTool needs it's information to be properly initialized!
    ana->SaveProcessedInfo(total_evts_pro, sum_of_evt_w, sum_of_evt_w_sq);
    ana->SetTotalEventsToProcess(entries);

    if(ana->initialize() != 0) {
        Error(APP_NAME, "cannot initialize %s", anaName.Data());
        delete ana;
        exit(1);
    }
    if(do_debug) ana->SetVerbose();

    // ana->GetSUSYTool();


    for( Long64_t entry = 0; entry < entries; ++entry )
    {
        ana->ClearBranch();
        ana->process(entry);

        // Remember to delete them
        store.clear();
    }

    Info( APP_NAME, "finished; Cleaning..." );

    delete ana;

    Info( APP_NAME, "Successfully finished analysis; Exitting..." );
    return 0;
}

