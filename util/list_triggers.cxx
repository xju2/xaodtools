#include <stdlib.h>
#include <utility>
#include <sstream>
#include <iostream>
#include <vector>

// ROOT include(s):
#include <TFile.h>
#include <TChain.h>

#ifdef ROOTCORE 
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#endif

// EDM include(s):
#include "xAODEventInfo/EventInfo.h"

#include "CPAnalysisExamples/errorcheck.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"

#include "MyXAODTools/CPToolsHelper.h"
#include "MyXAODTools/Helper.h"

using namespace std;
pair<long, long> get_run_event(const string& input){
    istringstream iss(input);
    long run, event;
    char dummy;
    iss >> run >> dummy >> event;
    return make_pair(run, event);
}
bool match_candiate(const vector<pair<long, long> >& candidates, long run, long event){
    bool res = false;
    for(auto& can : candidates){
        if(can.first == run && can.second == event) {
            res = true;
            break;
        }
    }
    return res;
}

int main( int argc, char* argv[]) 
{
    const char* APP_NAME = "list_triggers";
    if ((argc > 1 && string(argv[1]) == "help") ||(argc < 2))
    {
        cout << argv[0] << " toberun.txt [run,event;]" << endl;
        exit(1);
    }

    CHECK( xAOD::Init( APP_NAME ) );
    string inputFileName(argv[1]);
    Info( APP_NAME, "Input file name: %s", inputFileName.c_str() );

    vector<pair<long, long> > candidates;
    if(argc > 2) {
        string input(argv[2]);
        if (input.find(";") != string::npos){
            vector<string> tokens;
            MyXAODTools::tokenizeString(input, ';', tokens);
            for(auto& input_str : tokens){
                candidates.push_back( get_run_event(input_str) );
            }
        } else {
            candidates.push_back( get_run_event(input) );
        }
    }
    cout <<"# of candidates: " << candidates.size() << endl;


    TChain* fc = new TChain("CollectionTree");
    fc->Add(inputFileName.c_str());

    // Create a TEvent object:
    xAOD::TEvent event( xAOD::TEvent::kClassAccess );
    CHECK( event.readFrom( fc ) );
    Info( APP_NAME, "Number of events in the chain: %i",
            static_cast< int >( event.getEntries() ) );

    xAOD::TStore store;

    bool is_data = true;


    // Tell the object which entry to look at:
    event.getEntry(0);
    const xAOD::EventInfo* ei = 0;
    CHECK( event.retrieve( ei, "EventInfo" ) );
    if(ei->eventType(xAOD::EventInfo::IS_SIMULATION)) {
        is_data = false;
    }

    string maindir(getenv("ROOTCOREBIN"));
    string susy_config = Form("%s/data/MyXAODTools/upsilon.conf", maindir.c_str());
    ST::SUSYObjDef_xAOD* objTool = CPToolsHelper::GetSUSYTools(is_data, susy_config.c_str());

    int n_found = 0;
    Long64_t nentries = event.getEntries();
    if(candidates.size() == 0){
        nentries = 1;
    }

    for(Long64_t ientry = 0; ientry < nentries; ientry++){

        event.getEntry(ientry);
        CHECK( event.retrieve( ei, "EventInfo" ) );
        long run_number = (long) ei->runNumber();
        long event_number = (long) ei->eventNumber();

        bool found = false;
        if( candidates.size() == 0 ||
                (candidates.size() > 0 &&
                n_found < (int)candidates.size() &&
                match_candiate(candidates, run_number, event_number))
          )
        {
            found = true;
        } else if (found == (int)candidates.size()){
            break;
        }
        if(!found) continue;
        n_found ++;

        cout <<"Processing: " << ei->runNumber()<<" " << ei->eventNumber() << endl;
        auto chainGroup = objTool->GetTrigChainGroup("HLT_.*");
        for (auto& trig : chainGroup->getListOfTriggers())
        {
            cout << " " << trig << " is passed: " << objTool->IsTrigPassed(trig) << endl;
        }
    }

    delete objTool;
    delete fc;

    return 0;
}
