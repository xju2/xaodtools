/*
 * =====================================================================================
 *
 *       Filename:  SmearedInfo.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/27/2015 07:58:39 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xiangyang Ju (), xiangyang.ju@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef __MYXAODTOOLS_SMEAREDINFO_H_
#define __MYXAODTOOLS_SMEAREDINFO_H_

#include <vector>
#include "MyXAODTools/BranchCreatorBase.h"

class SmearedInfoBranch : public BranchCreatorBase 
{
public:
    SmearedInfoBranch();
    virtual ~SmearedInfoBranch();

    typedef struct SmearedInfo 
    {
        double leading_jet_pt_;
        double met_;
        double min_jets_met_;
        double dphi_EP_;
        int n_good_jets_;
    } SmearedInfo;

    void AttachBranchToTree(TTree& );
    void ClearBranch();
    bool CreateBranch();
    void Fill(const SmearedInfo& smeared);

private:
    double leading_jet_pt_;
    double met_;
    double min_jets_met_;
    double dphi_EP_;
    int n_good_jets_;
};

#endif
