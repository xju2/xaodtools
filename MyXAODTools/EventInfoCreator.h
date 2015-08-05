#ifndef __MYXAODTOOLS_EVENTINFOCREATOR_H__
#define __MYXAODTOOLS_EVENTINFOCREATOR_H__

#include <TTree.h>

#include "xAODEventInfo/EventInfo.h"

#include "MyXAODTools/BranchCreatorBase.h"
class EventInfoCreator
{

public: 
    EventInfoCreator();
    EventInfoCreator(bool is_data);
    virtual ~EventInfoCreator();
    // add full set of branches
    void AttachBranchToTree(TTree& tree);
    // add minimum set of branches
    void AttachMiniToTree(TTree& tree);
    void Fill(const xAOD::EventInfo& ei);
    void ClearBranch();

private:
    bool is_data_;

    int run_number_;
    int event_number_;
    int mc_channel_number_;
    float mc_weight_;

    float actualIPC_; // actual interaction per crossing for the current BCID -- for in-time pile-up
    float averageIPC_;// average interaction per corssing for all BCIDs -- for out-of-time pile-up
    uint32_t bcid_;
    uint32_t lumiBlock_;
};
#endif
