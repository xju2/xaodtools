#include "MyXAODTools/EventInfoCreator.h"
#include <stdlib.h>

EventInfoCreator::EventInfoCreator(){
    is_data_ = false;
}

EventInfoCreator::EventInfoCreator(bool is_data):
    is_data_(is_data)
{
}

EventInfoCreator::~EventInfoCreator(){

}

void EventInfoCreator::AttachMiniToTree(TTree& tree){
    tree.Branch("RunNumber", &run_number_, "RunNumber/I");
    tree.Branch("EventNumber", &event_number_, "EventNumber/I");
    tree.Branch("mc_channel_number", &mc_channel_number_, "mc_channel_number/I");
    tree.Branch("MCWeight", &mc_weight_, "MCWeight/F");
}

void EventInfoCreator::AttachBranchToTree(TTree& tree)
{
    AttachMiniToTree(tree);
    tree.Branch("actualIPC", &actualIPC_, "actualIPC/F");
    tree.Branch("averageIPC", &averageIPC_, "averageIPC/F");
    tree.Branch("bcid", &bcid_, "bcid/i");
    tree.Branch("lumiblock", &lumiBlock_, "lumiblock/i");
}

void EventInfoCreator::ClearBranch(){
    run_number_ = -1;
    event_number_ = -1;
}

void EventInfoCreator::Fill(const xAOD::EventInfo& ei){
    run_number_ = ei.runNumber();
    event_number_ = ei.eventNumber();
    if(ei.eventType(xAOD::EventInfo::IS_SIMULATION)) 
    {
        mc_channel_number_= ei.mcChannelNumber();
        const vector<float>& weights = ei.mcEventWeights();
        if(weights.size() > 0) mc_weight_ = weights[0];
        else mc_weight_ = 1.0;
    }
    actualIPC_ = ei.actualInteractionsPerCrossing();
    averageIPC_ = ei.averageInteractionsPerCrossing();
    bcid_ = ei.bcid();
    lumiBlock_ = ei.lumiBlock();
}
