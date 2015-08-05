#ifndef __MYXAODTOOLS_EVENTCOUNTER_H__
#define __MYXAODTOOLS_EVENTCOUNTER_H__

#include <string>
#include <map>
#include <utility>

#include "TChain.h"


class EventCounter
{

private:
    bool is_dc14_;
    bool is_mc_;
    std::string root_core_path;
    std::vector<int>* zee_mc_id_vec;
    std::map<int, double> total_events_dic_; //key:mc_channel_number, value: total_events
    std::map<int, double> all_events_noweight_dic_; //key:mc_channel_number, value: total_events
    std::map<int, double> cross_section_dic_;  //key:mc_channel_number, value: cross section

public:
    EventCounter(bool is_dc14);
    EventCounter(bool is_dc14, bool is_mc);
    EventCounter();
    ~EventCounter();
    void GetTotalEventsDic(const char* file_name, const char* tree_name);
    double getCrossSection(int mc_channel_number);
    double getTotalEvents(int mc_channel_number);
    double getTotalEventsNoWeight(int mc_channel_number);
    static TChain* getChain(const char* filename, const char* chainName);
    void printTotalEvents();
    void SaveTotalEvents(const char* file_dir = NULL);
    bool isZeeSample(int mc_id);
    void ReadEventInfo(const char* input_name);

private:
    void initial_cs(const char* filename);
    void GetCrossSection();
    double getUtil(int mc_channel_number, std::map<int, double>&, const char*);
    void LoadXS(const char* in_name);
};
#endif
