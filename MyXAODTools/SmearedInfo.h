#ifndef __MYXAODTOOLS_SMEAREDINFO_H__
#define __MYXAODTOOLS_SMEAREDINFO_H__

typedef struct SmearedInfo 
{
    float leading_jet_pt_;
    float leading_jet_eta_;
    float leading_jet_phi_;
    float sub_leading_jet_pt_;
    float sub_leading_jet_eta_;
    float sub_leading_jet_phi_;
    float HT_;
    float met_;
    float sum_et_;
    float min_jets_met_;
    float dphi_EP_;
    uint32_t n_good_jets_;
    uint32_t n_vertices_;
    float l3rd_jet_pt_;
    float l3rd_jet_eta_;
    float l3rd_jet_phi_;
    float l4th_jet_pt_;
    float l4th_jet_eta_;
    float l4th_jet_phi_;
} SmearedInfo;

#endif
