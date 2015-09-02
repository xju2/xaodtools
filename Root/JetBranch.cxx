#include <stdlib.h>

#include "MyXAODTools/JetBranch.h"
#include "CPAnalysisExamples/errorcheck.h"

const char* JetBranch::APP_NAME = "JetBranch";

JetBranch::JetBranch(){
    CreateBranch();
}

bool JetBranch::CreateBranch()
{
    emF_ = new vector<float>();
    hecF_ = new vector<float>();
    larQ_ = new vector<float>();
    hecQ_ = new vector<float>();
    sumpttrk_ = new vector<float>();
    frac_sampling_max_ = new vector<float>();
    negE_ = new vector<float>();
    avg_larQF_ = new vector<float>();
    frac_sampling_max_index_ = new vector<int>();
    jet_isBadTight_ = new vector<bool>();
    jet_timing_ = new vector<float>();

    jetCleaningTool_ = new JetCleaningTool("JetCleaningToolTight");
    CHECK(jetCleaningTool_->setProperty("CutLevel", "TightBad"));
    CHECK(jetCleaningTool_->initialize());
    return true;
}

JetBranch::~JetBranch(){
    delete emF_;
    delete hecF_;
    delete larQ_;
    delete hecQ_;
    delete sumpttrk_;
    delete frac_sampling_max_;
    delete negE_;
    delete avg_larQF_;
    delete frac_sampling_max_index_;
    delete jet_isBadTight_;
    delete jet_timing_;
    delete jetCleaningTool_;
}

void JetBranch::ClearBranch(){
    emF_->clear();
    hecF_->clear();
    larQ_->clear();
    hecQ_->clear();
    sumpttrk_->clear();
    frac_sampling_max_->clear();
    negE_->clear();
    avg_larQF_->clear();
    frac_sampling_max_index_->clear();
    jet_isBadTight_->clear();
    jet_timing_->clear();
}

void JetBranch::AttachBranchToTree(TTree& tree){
    tree.Branch("jet_emf", &emF_);
    tree.Branch("jet_hecf", &hecF_);
    tree.Branch("jet_LArQ", &larQ_);
    tree.Branch("jet_hecQ", &hecQ_);
    tree.Branch("jet_sumpt_trk", &sumpttrk_);
    tree.Branch("jet_fracSamplingMax", &frac_sampling_max_);
    tree.Branch("jet_negE", &negE_);
    tree.Branch("jet_avgLArQF", &avg_larQF_);
    tree.Branch("jet_fracSamplingMaxIndex", &frac_sampling_max_index_);
    tree.Branch("jet_isBadTight", &jet_isBadTight_);
    tree.Branch("jet_timing", &jet_timing_);
}

void JetBranch::Fill(const xAOD::Jet& jet) 
{
    std::vector<float> sumPtTrkvec;
    jet.getAttribute( xAOD::JetAttribute::SumPtTrkPt500, sumPtTrkvec );
    float sumpttrk = 0;
    if( ! sumPtTrkvec.empty() ) sumpttrk = sumPtTrkvec[0];

    emF_->push_back(jet.getAttribute<float>(xAOD::JetAttribute::EMFrac));
    hecF_->push_back(jet.getAttribute<float>(xAOD::JetAttribute::HECFrac));
    larQ_->push_back(jet.getAttribute<float>(xAOD::JetAttribute::LArQuality));
    hecQ_->push_back(jet.getAttribute<float>(xAOD::JetAttribute::HECQuality));
    sumpttrk_->push_back(sumpttrk);
    frac_sampling_max_->push_back(jet.getAttribute<float>(xAOD::JetAttribute::FracSamplingMax));
    negE_->push_back(jet.getAttribute<float>(xAOD::JetAttribute::NegativeE));
    avg_larQF_->push_back(jet.getAttribute<float>(xAOD::JetAttribute::AverageLArQF));
    frac_sampling_max_index_->push_back(jet.getAttribute<int>(xAOD::JetAttribute::FracSamplingMaxIndex));
    jet_timing_->push_back(jet.getAttribute<float>(xAOD::JetAttribute::Timing));

    jet_isBadTight_->push_back((bool)jetCleaningTool_->keep(jet));
}
