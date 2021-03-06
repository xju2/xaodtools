#ifndef __MYXAODTOOLS_LINKDEF_H__
#define __MYXAODTOOLS_LINKDEF_H__

#include <MyXAODTools/BranchCreatorBase.h>
#include <MyXAODTools/CPToolsHelper.h>
#include <MyXAODTools/EventCounter.h>
#include <MyXAODTools/EventInfoCreator.h>
#include <MyXAODTools/PhotonBranch.h>
#include <MyXAODTools/TrackBranch.h>
#include <MyXAODTools/ElectronBranch.h>
#include <MyXAODTools/MuonBranch.h>
#include <MyXAODTools/SmearedInfo.h>
#include <MyXAODTools/Candidate.h>
#include <MyXAODTools/HZZ4lHelper.h>

#ifdef __MAKECINT__
#pragma link C++ class SmearedInfo+;
#pragma link C++ class vector<SmearedInfo>+;
#endif

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;

#pragma link C++ class vector<TLorentzVector>;
#pragma link C++ class vector<SmearedInfo>+;
#pragma link C++ class BranchCreatorBase;
#pragma link C++ class CPToolsHelper;
#pragma link C++ class EventCounter;
#pragma link C++ class EventInfoCreator;
#pragma link C++ class PhotonBranch;
#pragma link C++ class TrackBranch;
#pragma link C++ class ElectronBranch;
#pragma link C++ class MuonBranch;
#pragma link C++ class Candidate;
#pragma link C++ class HZZ4lHelper;

#endif
#endif
