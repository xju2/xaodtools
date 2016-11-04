#!/usr/bin/env python

import ROOT
import sys
ROOT.gROOT.SetBatch()
import string
import math
from array import array
from sets import Set
import time
from optparse import OptionParser

if not hasattr(ROOT, "myText"):
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/AtlasUtils.C")
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c")

if not hasattr(ROOT, "passMultiLepton"):
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/work/upsilon/code/MyXAODTools/scripts/MultiLeptonDefs_new.cxx")

interested_event = [
    (195847, 9013)
]

CHI2_DIMUON_CUT = 3
CHI2_DIELE_CUT = 8
CHI2_4MUON_CUT = 50

LEADING_MUON_PT = 4E3
SUBLEADING_MUON_PT = 3E3

M12CUT_LOW = 9.2E3
M12CUT_HI = 9.7E3

M34CUT_LOW = 2E3
M34CUT_HI = 50E3

# author information
author_cuts_reco = {
    1 : "MuonBoy", # 0
    2 : "STACO",   # 1
    4 : "MuTag",   # 2
    8 : "MuidSA",  # 3
    16 : "MuidCo", # 4
    32 : "MuGirl", # 5
    64 : "MuGirlLowBeta", # 6
    128 : "CaloMuonId",   # 7
    256 : "CaloTag",      # 8
    512 : "CaloLikelihood", # 9
    1024 : "MuTagIMO",      # 10
    2048 : "MuonCombinedRefit",  # 11
    4096 : "ExtrapolateMuonToIP", # 12
}
MUON_AUTHORS = [
    "unknown", # 0
    "MuidCo",  # 1
    "STACO",   # 2
    "MuTag",   # 3
    "MuTagIMO",# 4
    "MuidSA",  # 5
    "MuGirl",  # 6
    "MuGirlLowBeta", # 7
    "CaloTag",       # 8
    "CaloLikelihood", # 9
    "ExtrapolatMuonToIP", #10
]

class BLSana:
    """
    select 4muons
    """
    def __init__(self, out_name):
        print "start to work"
        self.out_name = out_name

        # save output
        self.do_save = True

        # don't apply PV matching cut
        self.ignore_PV_cut = True

        # select upsilon first or quadruplet first
        self.upsilon_first = False

        # do 13 TeV
        self.do_13TeV = True

        # do debug
        self.m_debug = False

        self.cut_flow = ROOT.TH1F("cut_flow", "cut flow", 100, 0.5, 100.5)
        self.cut_flow_no3mu4 = ROOT.TH1F("cut_flow_no3mu4", "cut flow", 100, 0.5, 100.5)
        self.cut_flow_with_3mu4 = ROOT.TH1F("cut_flow_with3mu4", "cut flow", 100, 0.5, 100.5)

        self.out_events = ""

        # use new type
        self.use_new_type = True

        # use new pT
        self.use_new_pT = True

    def book_tree(self):
        self.out_tree = ROOT.TTree("bls", "bls")

        self.m4l = array('f', [0])
        self.m4l_track = array('f', [0])
        self.m4l_fitted = array('f', [0])
        self.m12 = array('f', [0])
        self.m34 = array('f', [0])
        self.dis_v1_x = array('f', [0])
        self.dis_v2_x = array('f', [0])
        self.dis_v1_v2 = array('f', [0])
        self.z_v1_v2 = array('f', [0])
        self.x_dis = array('f', [0])
        self.x_chi2 = array('f', [0])
        self.x_lxy = array('f', [0])
        self.x_track_pt = array('f', [0])
        self.x_fitted_pt = array('f', [0])
        self.x_type = array('f', [0])
        self.x_nDummyPV = array('f', [0])

        self.m1_pt = array('f', [0])
        self.m2_pt = array('f', [0])
        self.m3_pt = array('f', [0])
        self.m4_pt = array('f', [0])

        self.m1_track_pt = array('f', [0])
        self.m2_track_pt = array('f', [0])
        self.m3_track_pt = array('f', [0])
        self.m4_track_pt = array('f', [0])

        self.m1_track_eta = array('f', [0])
        self.m2_track_eta = array('f', [0])
        self.m3_track_eta = array('f', [0])
        self.m4_track_eta = array('f', [0])

        self.m1_track_p = array('f', [0])
        self.m2_track_p = array('f', [0])
        self.m3_track_p = array('f', [0])
        self.m4_track_p = array('f', [0])

        self.m1_trackD0PV = array('f', [0])
        self.m2_trackD0PV = array('f', [0])
        self.m3_trackD0PV = array('f', [0])
        self.m4_trackD0PV = array('f', [0])

        self.m1_trackD0SigPV = array('f', [0])
        self.m2_trackD0SigPV = array('f', [0])
        self.m3_trackD0SigPV = array('f', [0])
        self.m4_trackD0SigPV = array('f', [0])

        self.m1_trackZ0PV = array('f', [0])
        self.m2_trackZ0PV = array('f', [0])
        self.m3_trackZ0PV = array('f', [0])
        self.m4_trackZ0PV = array('f', [0])

        self.m1_type = array('i', [0])
        self.m2_type = array('i', [0])
        self.m3_type = array('i', [0])
        self.m4_type = array('i', [0])


        self.u1_chi2 = array('f', [0])
        self.u2_chi2 = array('f', [0])


        self.out_tree.Branch("m4l_fitted", self.m4l_fitted, "m4l_fitted/F")
        self.out_tree.Branch("m4l_track", self.m4l_track, "m4l_track/F")
        self.out_tree.Branch("m4l", self.m4l, "m4l/F")
        self.out_tree.Branch("m12", self.m12, "m12/F")
        self.out_tree.Branch("m34", self.m34, "m34/F")
        self.out_tree.Branch("dis_v1_x", self.dis_v1_x, "dis_v1_x/F")
        self.out_tree.Branch("dis_v2_x", self.dis_v2_x, "dis_v2_x/F")
        self.out_tree.Branch("dis_v1_v2", self.dis_v1_v2, "dis_v1_v2/F")
        self.out_tree.Branch("z_v1_v2", self.z_v1_v2, "z_v1_v2/F")
        self.out_tree.Branch("x_dis", self.x_dis, "x_dis/F")
        self.out_tree.Branch("x_chi2", self.x_chi2, "x_chi2/F")
        self.out_tree.Branch("x_lxy", self.x_lxy, "x_lxy/F")
        self.out_tree.Branch("x_track_pt", self.x_track_pt, "x_track_pt/F")
        self.out_tree.Branch("x_fitted_pt", self.x_fitted_pt, "x_fitted_pt/F")
        self.out_tree.Branch("x_type", self.x_type, "x_type/F")
        self.out_tree.Branch("x_nDummyPV", self.x_nDummyPV, "x_nDummyPV/F")

        self.out_tree.Branch("m1_pt", self.m1_pt, "m1_pt/F")
        self.out_tree.Branch("m2_pt", self.m2_pt, "m2_pt/F")
        self.out_tree.Branch("m3_pt", self.m3_pt, "m3_pt/F")
        self.out_tree.Branch("m4_pt", self.m4_pt, "m4_pt/F")

        self.out_tree.Branch("m1_track_pt", self.m1_track_pt, "m1_track_pt/F")
        self.out_tree.Branch("m2_track_pt", self.m2_track_pt, "m2_track_pt/F")
        self.out_tree.Branch("m3_track_pt", self.m3_track_pt, "m3_track_pt/F")
        self.out_tree.Branch("m4_track_pt", self.m4_track_pt, "m4_track_pt/F")

        self.out_tree.Branch("m1_track_eta", self.m1_track_eta, "m1_track_eta/F")
        self.out_tree.Branch("m2_track_eta", self.m2_track_eta, "m2_track_eta/F")
        self.out_tree.Branch("m3_track_eta", self.m3_track_eta, "m3_track_eta/F")
        self.out_tree.Branch("m4_track_eta", self.m4_track_eta, "m4_track_eta/F")

        self.out_tree.Branch("m1_track_p", self.m1_track_p, "m1_track_p/F")
        self.out_tree.Branch("m2_track_p", self.m2_track_p, "m2_track_p/F")
        self.out_tree.Branch("m3_track_p", self.m3_track_p, "m3_track_p/F")
        self.out_tree.Branch("m4_track_p", self.m4_track_p, "m4_track_p/F")

        self.out_tree.Branch("m1_trackD0PV", self.m1_trackD0PV, "m1_trackD0PV/F")
        self.out_tree.Branch("m2_trackD0PV", self.m2_trackD0PV, "m2_trackD0PV/F")
        self.out_tree.Branch("m3_trackD0PV", self.m3_trackD0PV, "m3_trackD0PV/F")
        self.out_tree.Branch("m4_trackD0PV", self.m4_trackD0PV, "m4_trackD0PV/F")
        self.out_tree.Branch("m1_trackD0SigPV", self.m1_trackD0SigPV, "m1_trackD0SigPV/F")
        self.out_tree.Branch("m2_trackD0SigPV", self.m2_trackD0SigPV, "m2_trackD0SigPV/F")
        self.out_tree.Branch("m3_trackD0SigPV", self.m3_trackD0SigPV, "m3_trackD0SigPV/F")
        self.out_tree.Branch("m4_trackD0SigPV", self.m4_trackD0SigPV, "m4_trackD0SigPV/F")
        self.out_tree.Branch("m1_trackZ0PV", self.m1_trackZ0PV, "m1_trackZ0PV/F")
        self.out_tree.Branch("m2_trackZ0PV", self.m2_trackZ0PV, "m2_trackZ0PV/F")
        self.out_tree.Branch("m3_trackZ0PV", self.m3_trackZ0PV, "m3_trackZ0PV/F")
        self.out_tree.Branch("m4_trackZ0PV", self.m4_trackZ0PV, "m4_trackZ0PV/F")

        self.out_tree.Branch("m1_type", self.m1_type, "m1_type/I")
        self.out_tree.Branch("m2_type", self.m2_type, "m2_type/I")
        self.out_tree.Branch("m3_type", self.m3_type, "m3_type/I")
        self.out_tree.Branch("m4_type", self.m4_type, "m4_type/I")

        # authors of each muon
        self.m1_author = array('f', [0])
        self.m2_author = array('f', [0])
        self.m3_author = array('f', [0])
        self.m4_author = array('f', [0])
        self.out_tree.Branch("m1_author", self.m1_author, "m1_author/F")
        self.out_tree.Branch("m2_author", self.m2_author, "m2_author/F")
        self.out_tree.Branch("m3_author", self.m3_author, "m3_author/F")
        self.out_tree.Branch("m4_author", self.m4_author, "m4_author/F")

        self.out_tree.Branch("u1_chi2", self.u1_chi2, "u1_chi2/F")
        self.out_tree.Branch("u2_chi2", self.u2_chi2, "u2_chi2/F")

        self.run = array('i', [0])
        self.event = array('i', [0])
        self.out_tree.Branch("run", self.run, "run/I")
        self.out_tree.Branch("event", self.event, "event/I")

        self.trig_3mu4 = array('i', [0])
        self.out_tree.Branch("trig_3mu4", self.trig_3mu4, "trig_3mu4/I")

        self.charge = array('f', [0])
        self.out_tree.Branch("charge", self.charge, "charge/F")

        # quality of muons
        self.m1_quality = array('i', [0])
        self.m2_quality = array('i', [0])
        self.m3_quality = array('i', [0])
        self.m4_quality = array('i', [0])
        self.out_tree.Branch("m1_quality", self.m1_quality, "m1_quality/I")
        self.out_tree.Branch("m2_quality", self.m2_quality, "m2_quality/I")
        self.out_tree.Branch("m3_quality", self.m3_quality, "m3_quality/I")
        self.out_tree.Branch("m4_quality", self.m4_quality, "m4_quality/I")

    def clear_tree(self):
        pass

    def book_upsilon(self):
        self.tree_onia = ROOT.TTree("upsilon", "upsilon")

        self.up_mass = ROOT.vector('float')()
        self.up_pt = ROOT.vector('float')()
        self.up_chi2 = ROOT.vector('float')()
        self.up_d0_1 = ROOT.vector('float')()
        self.up_d0_2 = ROOT.vector('float')()
        self.up_z0_1 = ROOT.vector('float')()
        self.up_z0_2 = ROOT.vector('float')()

        self.tree_onia.Branch("run", self.run, "run/I")
        self.tree_onia.Branch("event", self.event, "event/I")
        self.tree_onia.Branch("mass", self.up_mass)
        self.tree_onia.Branch("pt", self.up_pt)
        self.tree_onia.Branch("chi2", self.up_chi2)
        self.tree_onia.Branch("mu_d0_1", self.up_d0_1)
        self.tree_onia.Branch("mu_d0_2", self.up_d0_2)
        self.tree_onia.Branch("mu_z0_1", self.up_z0_1)
        self.tree_onia.Branch("mu_z0_2", self.up_z0_2)
        self.tree_onia.Branch("trig_3mu4", self.trig_3mu4, "trig_3mu4/I")

        self.up_pass_dionia = array('i', [0])
        self.tree_onia.Branch("pass_diOnia", self.up_pass_dionia, "pass_diOnia/I")

    def clear_upsilon(self):
        self.up_mass.clear()
        self.up_pt.clear()
        self.up_chi2.clear()
        self.up_d0_1.clear()
        self.up_d0_2.clear()
        self.up_z0_1.clear()
        self.up_z0_2.clear()
        self.up_pass_dionia[0] = 0

    def book_ss(self):
        self.tree_ss = ROOT.TTree("ss", "ss")

        self.ss_m4l = array('f', [0])
        self.ss_m12 = array('f', [0])
        self.ss_m34 = array('f', [0])
        self.tree_ss.Branch("m4l", self.ss_m4l, "m4l/F")
        self.tree_ss.Branch("m12", self.ss_m12, "m12/F")
        self.tree_ss.Branch("m34", self.ss_m34, "m34/F")

    def clear_ss(self):
        self.ss_m4l.clear()
        self.ss_m12.clear()
        self.ss_m34.clear()

    def book_hists(self):
        self.h_upsilon = ROOT.TH1F("h_upsilon", "upsilon mass;m_{#varUpsilon} [GeV];Events / 50 MeV", 400, 1, 21)
        self.h_jpsi = ROOT.TH1F("h_jpsi", "jpsi mass;m_{#mu+#mu} [GeV];Events / 100 MeV", 200, 1, 21)

        nbins_4l = 400
        low_4l = 10
        hi_4l = 30
        self.h_m4l = ROOT.TH1F("h_m4l", "m4l;m_{#varUpsilon+#mu+#mu} [GeV];Events / 50 MeV", nbins_4l, low_4l, hi_4l)
        self.h_m4l_LHS = ROOT.TH1F("h_m4l_LHS", "m4l LHS;LHS m_{#varUpsilon+#mu+#mu} [GeV];Events / 200 MeV", nbins_4l, low_4l, hi_4l)
        self.h_m4l_RHS_U1 = ROOT.TH1F("h_m4l_RHS_U1", "m4l RHS;RHS U1 m_{#varUpsilon+#mu+#mu} [GeV];Events / 200 MeV", nbins_4l, low_4l, hi_4l)
        self.h_m4l_RHS_U2 = ROOT.TH1F("h_m4l_RHS_U2", "m4l RHS;RHS U2 m_{#varUpsilon+#mu+#mu} [GeV];Events / 200 MeV", nbins_4l, low_4l, hi_4l)

        # quad properties
        self.h_dis_v1_quad = ROOT.TH1F("h_dis_v1_quad", ";D(v1,x);", 20, 0, 20)
        self.h_dis_v2_quad = ROOT.TH1F("h_dis_v2_quad", ";D(v2,x);", 20, 0, 20)
        self.h_dis_quad = ROOT.TH1F("h_dis_quad", ";D(4l);", 20, 0, 20)
        self.h_chi_quad = ROOT.TH1F("h_chi2_quad",  ";Quad #chi^{2};", 900, 0, 9E6)

    def fill_hists(self, file_names):
        print "debug:", self.m_debug
        tree = ROOT.TChain("physics", "physics")
        for file_ in file_names:
            print "adding.. ", file_
            tree.Add(file_)

        nentries = tree.GetEntries()
        print "total entries:",nentries
        if nentries == 0:
            self.do_save = False
            return

        imatched = 0
        start_time = time.time()
        #if hasattr(tree, "mu_z0_pv_sintheta"):
        #    self.do_13TeV = False

        print "do 13TeV:", self.do_13TeV

        for ientry in xrange(nentries):
            tree.GetEntry(ientry)
            self.clear_upsilon()
            self.clear_tree()

            if hasattr(tree, "Run"):
                run = tree.Run
                event = tree.Event
            else:
                run = tree.RunNumber
                event = tree.EventNumber

            self.run[0] = run
            self.event[0] = event

            self.fill_cut_flow(1)


            if self.m_debug and imatched == len(interested_event):
                break
            if self.m_debug and (run, event) not in interested_event:
                continue

            #if self.run[0] != 307710:
            #    continue

            if self.m_debug:
                imatched += 1

            if self.m_debug:
                print "process: ",run, event

            if ientry%50000 == 0:
                print "processed:", ientry,"with time: {:.2f} min".format((time.time()-start_time)/60.)


            if tree.n_muon < 2: # don't forget mu-mu-e-e
                continue

            # add trigger info
            if hasattr(tree, "trig_3mu4"):
                self.trig_3mu4[0] = int(tree.trig_3mu4)

            if self.upsilon_first and tree.n_muon < 4:
                results = self.select_upsilon(tree)
            else:
                results = self.select_v04(tree)

            if results is None:
                continue

            mU, m34, mass_id, onia1_id, onia2_id, quad_id, quad_type = results

            self.h_upsilon.Fill(mU)
            self.h_jpsi.Fill(m34)
            self.h_m4l.Fill(tree.quad_fitted_mass[mass_id]/1E3)

            onia1_x = tree.onia_x[onia1_id]
            onia1_y = tree.onia_y[onia1_id]
            onia1_z = tree.onia_z[onia1_id]
            onia2_x = tree.onia_x[onia2_id]
            onia2_y = tree.onia_y[onia2_id]
            onia2_z = tree.onia_z[onia2_id]
            quad_x = tree.quad_x[mass_id]
            quad_y = tree.quad_y[mass_id]
            quad_z = tree.quad_z[mass_id]

            dis_v1_quad = self.get_dis(onia1_x, onia1_y, onia1_z, quad_x, quad_y, quad_z)
            dis_v2_quad = self.get_dis(onia2_x, onia2_y, onia2_z, quad_x, quad_y, quad_z)
            dis_v1_v2 = self.get_dis(onia1_x, onia1_y, onia1_z, onia2_x, onia2_y, onia2_z)
            z_v1_v2 = abs(onia1_z-onia2_z)

            dis_quad = math.sqrt(quad_x**2 + quad_y**2 + quad_z**2)
            chi_quad = tree.quad_chi2[quad_id]

            self.h_dis_v1_quad.Fill(dis_v1_quad)
            self.h_dis_v2_quad.Fill(dis_v2_quad)
            self.h_dis_quad.Fill(dis_quad)
            self.h_chi_quad.Fill(chi_quad)

            self.m4l_fitted[0] =  tree.quad_fitted_mass[mass_id]/1E3
            self.m4l_track[0] =  tree.quad_track_mass[quad_id]/1E3
            self.m4l[0] =  tree.quad_mass[quad_id]/1E3
            self.m12[0] = mU
            self.m34[0] = m34
            self.dis_v1_x[0] = dis_v1_quad
            self.dis_v2_x[0] = dis_v2_quad
            self.dis_v1_v2[0] = dis_v1_v2
            self.z_v1_v2[0] = z_v1_v2
            self.x_dis[0] = dis_quad
            self.x_chi2[0] = chi_quad
            self.x_lxy[0] = math.sqrt(quad_x**2 + quad_y**2)
            self.x_track_pt[0] = tree.quad_track_pt[quad_id]
            #self.x_fitted_pt[0] = tree.quad_fitted_pt[mass_id]
            self.x_type[0] = quad_type

            mu1_id =  tree.quad_id1[quad_id]
            mu2_id =  tree.quad_id2[quad_id]
            mu3_id =  tree.quad_id3[quad_id]
            mu4_id =  tree.quad_id4[quad_id]

            self.m1_track_pt[0] = tree.mu_track_pt[ mu1_id ]
            self.m2_track_pt[0] = tree.mu_track_pt[ mu2_id ]
            if abs(quad_type) < 1E-6:
                self.m3_track_pt[0] = tree.mu_track_pt[ mu3_id ]
                self.m4_track_pt[0] = tree.mu_track_pt[ mu4_id ]
            else:
                self.m3_track_pt[0] = tree.el_pt[ mu3_id ]
                self.m4_track_pt[0] = tree.el_pt[ mu4_id ]

            self.m1_track_eta[0] = tree.mu_track_eta[ mu1_id ]
            self.m2_track_eta[0] = tree.mu_track_eta[ mu2_id ]
            if abs(quad_type) < 1E-6:
                self.m3_track_eta[0] = tree.mu_track_eta[ mu3_id ]
                self.m4_track_eta[0] = tree.mu_track_eta[ mu4_id ]
            else:
                self.m3_track_eta[0] = tree.el_eta[ mu3_id ]
                self.m4_track_eta[0] = tree.el_eta[ mu4_id ]

            self.m1_pt[0] = tree.mu_pt[ mu1_id ]
            self.m2_pt[0] = tree.mu_pt[ mu2_id ]
            if abs(quad_type) < 1E-6:
                self.m3_pt[0] = tree.mu_pt[ mu3_id ]
                self.m4_pt[0] = tree.mu_pt[ mu4_id ]
            else:
                self.m3_pt[0] = tree.el_pt[ mu3_id ]
                self.m4_pt[0] = tree.el_pt[ mu4_id ]

            # get momentum
            self.m1_track_p[0] = tree.mu_track_pt[mu1_id]*ROOT.TMath.CosH(tree.mu_track_eta[mu1_id])
            self.m2_track_p[0] = tree.mu_track_pt[mu2_id]*ROOT.TMath.CosH(tree.mu_track_eta[mu2_id])
            self.m3_track_p[0] = tree.mu_track_pt[mu3_id]*ROOT.TMath.CosH(tree.mu_track_eta[mu3_id])
            self.m4_track_p[0] = tree.mu_track_pt[mu4_id]*ROOT.TMath.CosH(tree.mu_track_eta[mu4_id])

            if self.do_13TeV:
                self.m1_trackD0PV[0] = tree.mu_d0[ tree.quad_id1[quad_id] ]
                self.m2_trackD0PV[0] = tree.mu_d0[ tree.quad_id2[quad_id] ]
                self.m3_trackD0PV[0] = tree.mu_d0[ tree.quad_id3[quad_id] ]
                self.m4_trackD0PV[0] = tree.mu_d0[ tree.quad_id4[quad_id] ]
                self.m1_trackD0SigPV[0] = tree.mu_d0_sig[ tree.quad_id1[quad_id] ]
                self.m2_trackD0SigPV[0] = tree.mu_d0_sig[ tree.quad_id2[quad_id] ]
                self.m3_trackD0SigPV[0] = tree.mu_d0_sig[ tree.quad_id3[quad_id] ]
                self.m4_trackD0SigPV[0] = tree.mu_d0_sig[ tree.quad_id4[quad_id] ]
                self.m1_trackZ0PV[0] = tree.mu_z0_sintheta[ tree.quad_id1[quad_id] ]
                self.m2_trackZ0PV[0] = tree.mu_z0_sintheta[ tree.quad_id2[quad_id] ]
                self.m3_trackZ0PV[0] = tree.mu_z0_sintheta[ tree.quad_id3[quad_id] ]
                self.m4_trackZ0PV[0] = tree.mu_z0_sintheta[ tree.quad_id4[quad_id] ]

            elif abs(quad_type) == 0 and hasattr(tree, "mu_d0_pv_sig"):
                self.m1_trackD0PV[0] = tree.mu_d0_pv[ tree.quad_id1[quad_id] ]
                self.m2_trackD0PV[0] = tree.mu_d0_pv[ tree.quad_id2[quad_id] ]
                self.m3_trackD0PV[0] = tree.mu_d0_pv[ tree.quad_id3[quad_id] ]
                self.m4_trackD0PV[0] = tree.mu_d0_pv[ tree.quad_id4[quad_id] ]
                self.m1_trackD0SigPV[0] = tree.mu_d0_pv_sig[ tree.quad_id1[quad_id] ]
                self.m2_trackD0SigPV[0] = tree.mu_d0_pv_sig[ tree.quad_id2[quad_id] ]
                self.m3_trackD0SigPV[0] = tree.mu_d0_pv_sig[ tree.quad_id3[quad_id] ]
                self.m4_trackD0SigPV[0] = tree.mu_d0_pv_sig[ tree.quad_id4[quad_id] ]
                self.m1_trackZ0PV[0] = tree.mu_z0_pv[ tree.quad_id1[quad_id] ]
                self.m2_trackZ0PV[0] = tree.mu_z0_pv[ tree.quad_id2[quad_id] ]
                self.m3_trackZ0PV[0] = tree.mu_z0_pv[ tree.quad_id3[quad_id] ]
                self.m4_trackZ0PV[0] = tree.mu_z0_pv[ tree.quad_id4[quad_id] ]

            self.u1_chi2[0] = tree.onia_chi2[ onia1_id ]
            self.u2_chi2[0] = tree.onia_chi2[ onia2_id ]

            # get total charge
            if quad_type == 0:
                self.charge[0] = tree.mu_charge[mu1_id]+tree.mu_charge[mu2_id]+tree.mu_charge[mu3_id]+tree.mu_charge[mu4_id]
            else:
                self.charge[0] = tree.mu_charge[mu1_id]+tree.mu_charge[mu2_id]+tree.el_charge[mu3_id]+tree.el_charge[mu4_id]

            # charge of muons
            if abs(quad_type) == 0:
                if self.use_new_type:
                    self.m1_type[0] = self.get_muon_type(tree, mu1_id)
                    self.m2_type[0] = self.get_muon_type(tree, mu2_id)
                    self.m3_type[0] = self.get_muon_type(tree, mu3_id)
                    self.m4_type[0] = self.get_muon_type(tree, mu4_id)
                else:
                    self.m1_type[0] = tree.mu_type[ mu1_id ]
                    self.m2_type[0] = tree.mu_type[ mu2_id ]
                    self.m3_type[0] = tree.mu_type[ mu3_id ]
                    self.m4_type[0] = tree.mu_type[ mu4_id ]

                # add author information
                self.m1_author[0] = tree.mu_author[mu1_id]
                self.m2_author[0] = tree.mu_author[mu2_id]
                self.m3_author[0] = tree.mu_author[mu3_id]
                self.m4_author[0] = tree.mu_author[mu4_id]

                if hasattr(tree, "mu_quality"):
                    self.m1_quality[0] = tree.mu_quality[mu1_id]
                    self.m2_quality[0] = tree.mu_quality[mu1_id]
                    self.m3_quality[0] = tree.mu_quality[mu1_id]
                    self.m4_quality[0] = tree.mu_quality[mu1_id]

            self.out_tree.Fill()

    def select_upsilon(self, tree):
        has_upsilon = False
        mU = -1
        can_id = []
        onia1_id = 0
        upsilon_chi2 = 9999
        for i,onia_mass in enumerate(tree.onia_fitted_mass):
            id1 = tree.onia_id1[i]
            id2 = tree.onia_id2[i]
            #print id1, id2

            # require to be muon-muon pair
            if hasattr(tree, "onia_type") and tree.onia_type[i] != 0:
                if self.m_debug:
                    print "failed onia type cut",tree.onia_type[i]
                continue

            if tree.mu_track_pt[id1] > LEADING_MUON_PT and\
               tree.mu_track_pt[id2] > LEADING_MUON_PT and\
               tree.onia_chi2[i] < CHI2_DIMUON_CUT and\
               onia_mass > M12CUT_LOW and\
               onia_mass < M12CUT_HI and\
               tree.onia_chi2[i] < upsilon_chi2:

                has_upsilon = True
                upsilon_chi2 = tree.onia_chi2[i]
                mU = onia_mass/1E3
                if len(can_id) < 2:
                    can_id.append(id1)
                    can_id.append(id2)
                else:
                    can_id[0] = id1
                    can_id[1] = id2

                onia1_id = i
            else:
                if self.m_debug:
                    print "failed pT cuts"

        if not has_upsilon:
            if self.m_debug:
                print "cannot find upsilon"
            return None

        if self.m_debug:
            print "Upsilon mass: ", mU, self.event[0], self.run[0]

        m34 = -1
        onia2_chi2 = 99999
        hadd_add_ee = False
        quad_type = -1
        chi2_cuts = -1
        has_4muon = False
        #first find a 4muon, if not find mumu-ee
        for i,onia_mass in enumerate(tree.onia_fitted_mass):
            id1 = tree.onia_id1[i]
            id2 = tree.onia_id2[i]
            if i== onia1_id:
                continue

            if id1 in can_id or id2 in can_id:
                continue

            if hasattr(tree, "onia_type") and tree.onia_type[i] == 1:
                continue
            else:
                quad_type = 0
                ## 4mu channel
                track_pt1 = tree.mu_track_pt[id1]
                track_pt2 = tree.mu_track_pt[id2]
                chi2_cut = CHI2_DIMUON_CUT

            if track_pt1 > SUBLEADING_MUON_PT and\
               track_pt2 > SUBLEADING_MUON_PT and\
               tree.onia_chi2[i] < chi2_cut and\
               onia_mass > M34CUT_LOW and\
               onia_mass < M34CUT_HI and\
               tree.onia_chi2[i] < onia2_chi2:

                has_4muon = True
                onia2_chi2 = tree.onia_chi2[i]
                m34 = onia_mass/1E3
                onia2_id = i
                if len(can_id) < 4:
                    can_id.append(id1)
                    can_id.append(id2)
                else:
                    can_id[2] = id1
                    can_id[3] = id2

        if not has_4muon and hasattr(tree, "onia_type"):
            for i,onia_mass in enumerate(tree.onia_fitted_mass):
                id1 = tree.onia_id1[i]
                id2 = tree.onia_id2[i]
                if i == onia1_id:
                    continue
                if tree.onia_type[i] == 0:
                    continue

                if not self.passElectronID(tree, id1) or\
                   not self.passElectronID(tree, id2):
                    continue

                quad_type = 1
                track_pt1 = tree.el_pt[id1]
                track_pt2 = tree.el_pt[id2]
                chi2_cut = CHI2_DIELE_CUT
                if track_pt1 > SUBLEADING_MUON_PT and\
                   track_pt2 > SUBLEADING_MUON_PT and\
                   tree.onia_chi2[i] < chi2_cut and\
                   onia_mass > M34CUT_LOW and\
                   onia_mass < M34CUT_HI and\
                   tree.onia_chi2[i] < onia2_chi2:
                    onia2_chi2 = tree.onia_chi2[i]
                    m34 = onia_mass/1E3
                    onia2_id = i
                    if len(can_id) < 4:
                        can_id.append(id1)
                        can_id.append(id2)
                    else:
                        can_id[2] = id1
                        can_id[3] = id2

        if m34 < 0:
            if self.m_debug:
                print "cannot second pair"
            return None

        if self.m_debug:
            print "muonID:"," ".join([str(x) for x in can_id])

        mass_id = -1
        quad_id = -1
        j = -1
        for i,x_chi2 in enumerate(tree.quad_chi2):
            if x_chi2 < 0:
                if self.m_debug:
                    print self.run[0],self.event[0],"Fit failed",x_chi2
                continue

            j += 1
            id1 = tree.quad_id1[i]
            id2 = tree.quad_id2[i]
            id3 = tree.quad_id3[i]
            id4 = tree.quad_id4[i]
            if id1 not in can_id or\
               id2 not in can_id or\
               id3 not in can_id or\
               id4 not in can_id:
                continue

            if tree.quad_nCombined[i] < 3:
                if self.m_debug:
                    print "less than 3 combined muons"
                continue

            ## vertex association cut
            pvID1 = tree.mu_pvID[can_id[0]]
            pvID2 = tree.mu_pvID[can_id[1]]
            if abs(quad_type) < 1E-6:
                pvID3 = tree.mu_pvID[can_id[2]]
                pvID4 = tree.mu_pvID[can_id[3]]
            else:
                ### not available for electrons, just let it pass!
                pvID3 = pvID1
                pvID4 = pvID1

            pass_vertex = (pvID1 == pvID2) and (pvID3 == pvID2 or pvID4 == pvID2)
            #pass_vertex = pvID1 == 0 and pvID1 == pvID2 and pvID3 == pvID2 and pvID4 == pvID3

            if not self.ignore_PV_cut and not pass_vertex:
                if self.m_debug:
                    print "not from same vertex",pvID1,pvID2,pvID3,pvID4
                continue

            mass_id = j
            quad_mass = tree.quad_fitted_mass[mass_id]
            if self.m_debug:
                print "mass-> ",quad_mass,mass_id
            quad_id = i

        if quad_id < 0:
            return None

        if self.m_debug:
            print "quadID:", quad_id, "mass:", quad_mass," massID:", mass_id,"quadType:",quad_type

        return (mU, m34, mass_id, onia1_id, onia2_id, quad_id, quad_type)

    def select_v04(self, tree):
        """
        corresponding to Note v0p7, apply pT > 4 GeV on leading two muons
        """
        ##first select four good muons
        good_cb_muons = []
        good_st_muons = []
        combined_type_cut = 0
        if not self.do_13TeV:
            combined_type_cut = 1

        for i in sorted(range(len(tree.mu_track_pt)), key=lambda k:tree.mu_track_pt[k], reverse=True):
            if not self.passMuonID(tree, i):
                continue

            if self.use_new_type:
                type_muon = self.get_muon_type(tree, i)
            else:
                type_muon = tree.mu_type[i]

            if type_muon == combined_type_cut:
                good_cb_muons.append(i)
            elif type_muon == 2:
                good_st_muons.append(i)
            else:
                continue

        if len(good_cb_muons) + len(good_st_muons) > 3:
            self.fill_cut_flow(2)

        good_muons = []
        if len(good_cb_muons) >= 4:
            good_muons = good_cb_muons[0:4]
        elif len(good_cb_muons) == 3 and len(good_st_muons) >= 1:
            good_muons = good_cb_muons + [good_st_muons[0]]
        else:
            return None
        self.fill_cut_flow(3)

        # for chi2 studies
        self.fill_onia(tree, good_muons)

        # find quadruplet pair
        quad_id = -1
        for i in range(len(tree.quad_id1)):
            id1 = tree.quad_id1[i]
            id2 = tree.quad_id2[i]
            id3 = tree.quad_id3[i]
            id4 = tree.quad_id4[i]
            if id1 not in good_muons or\
               id2 not in good_muons or\
               id3 not in good_muons or\
               id4 not in good_muons:
                continue

            quad_id = i
            break

        pass_dionia = False
        if quad_id >= 0:
            m4l = tree.quad_fitted_mass[quad_id]
            if m4l < 50E3 and m4l > 0:
                # onia cuts
                onia_pair_index = self.find_onia_pair(tree, good_muons, self.passOniaCuts)
                if len(onia_pair_index) > 0:
                    pass_dionia = True

        # neutral charge
        total_charge = 0
        n_muon_with_pT_gt_4GeV = 0
        for i in good_muons:
            total_charge += tree.mu_charge[i]
            if tree.mu_track_pt[i] > 4E3:
                n_muon_with_pT_gt_4GeV += 1

        self.up_pass_dionia[0] = int(pass_dionia and total_charge == 0)
        self.tree_onia.Fill()
        if not pass_dionia:
            return None
        self.fill_cut_flow(4)

        if self.m_debug:
            print "onia_pars: "
            print onia_pair_index

        charge_weight = 0
        if abs(total_charge) == 2:
            charge_weight = 10
        elif abs(total_charge) == 0:
            charge_weight = 0
            self.print_event(tree)
        else:
            return None

        self.fill_cut_flow(5+charge_weight)

        # if pass 4 GeV cut
        if self.use_new_pT and n_muon_with_pT_gt_4GeV < 2:
            return None
        # if has upsilon
        chi2_upsilon = 9E9
        id_upsilon = [-1, -1, -1, -1, -1, -1]
        for i,j in onia_pair_index:
            mu_id1 = tree.onia_id1[i]
            mu_id2 = tree.onia_id2[i]
            mu_id3 = tree.onia_id1[j]
            mu_id4 = tree.onia_id2[j]

            if self.passUpsilon(tree, i):
                id_upsilon[0] = i
                id_upsilon[1] = mu_id1
                id_upsilon[2] = mu_id2
                id_upsilon[3] = j
                id_upsilon[4] = mu_id3
                id_upsilon[5] = mu_id4
            elif self.passUpsilon(tree, j):
                id_upsilon[0] = j
                id_upsilon[1] = mu_id3
                id_upsilon[2] = mu_id4
                id_upsilon[3] = i
                id_upsilon[4] = mu_id1
                id_upsilon[5] = mu_id2
            else:
                pass

        if id_upsilon[0] < 0:
            return None
        self.fill_cut_flow(6+charge_weight)

        # count number of combined muons
        ncombined = 0
        for i in good_muons:
            if self.use_new_type:
                type_muon = self.get_muon_type(tree, i)
            else:
                type_muon = tree.mu_type[i]

            if type_muon == 0:
                ncombined += 1

        if ncombined == 4:
            self.fill_cut_flow(7+charge_weight)

        onia1_id = id_upsilon[0]
        onia2_id = id_upsilon[3]
        m12 = tree.onia_fitted_mass[onia1_id]
        m34 = tree.onia_fitted_mass[onia2_id]
        if self.m_debug:
            print "mass: ", m12, m34, m4l

        if hasattr(tree, 'trig_3mu4') and tree.trig_3mu4:
            self.fill_cut_flow(8+charge_weight)

        if m4l < 50E3 and m4l > 0:
            self.fill_cut_flow(9+charge_weight)

        return (m12, m34, quad_id, onia1_id, onia2_id, quad_id, 0)

    def get_dis(self, x1, y1, z1, x2, y2, z2):
        return math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

    def save_hists(self):
        if not self.do_save:
            return

        # save output
        fout = ROOT.TFile.Open(out_name, "recreate")
        self.cut_flow.Write()
        self.cut_flow_with_3mu4.Write()
        self.cut_flow_no3mu4.Write()
        self.h_upsilon.Write()
        self.h_m4l.Write()
        self.h_m4l_LHS.Write()
        self.h_m4l_RHS_U1.Write()
        self.h_m4l_RHS_U2.Write()
        self.h_jpsi.Write()

        self.h_dis_v1_quad.Write()
        self.h_dis_v2_quad.Write()
        self.h_dis_quad.Write()
        self.h_chi_quad.Write()

        self.out_tree.Write()
        self.tree_onia.Write()

        with open("print_cutflow.txt", 'w') as f:
            f.write(self.out_events)
        fout.Close()

    def passElectronID(self, tree, el_id):
        et = tree.el_etCluster[el_id]
        eta = tree.el_etas2[el_id]
        rhad = tree.el_rhad[el_id]
        rhad1 = tree.el_rhad1[el_id]
        reta = tree.el_reta[el_id]
        w2 = tree.el_w2[el_id]
        f1 = tree.el_f1[el_id]
        f3 = tree.el_f3[el_id]
        wstot = tree.el_wstot[el_id]
        DEmaxs1 = tree.el_DEmaxs1[el_id]
        deltaEta = tree.el_deltaEta[el_id]
        deltaPhi = tree.el_deltaPhi[el_id]
        nSi = tree.el_nSi[el_id]
        nSiDeadSensors = tree.el_nSiDeadSensors[el_id]
        nPix = tree.el_nPix[el_id]
        nPixDeadSensors = tree.el_nPixDeadSensors[el_id]
        nTRThigh = tree.el_nTRThigh[el_id]
        nTRThighOutliers = tree.el_nTRThighOutliers[el_id]
        nTRT = tree.el_nTRT[el_id]
        nTRTOutliers = tree.el_nTRTOutliers[el_id]
        eoverp = tree.el_EoverP[el_id]
        rTRT = tree.el_rTRT[el_id]
        expectedBLayer = tree.el_expectedBLayer[el_id]
        trackd0 = tree.el_trackd0[el_id]
        nTRTtotal = nTRT + nTRTOutliers

        res = ROOT.passMultiLepton(eta, et, rhad, rhad1, reta, w2, f1, f3, wstot, DEmaxs1, deltaEta, nSi,nSiDeadSensors, nPix, nPixDeadSensors, deltaPhi, eoverp, rTRT, nTRTtotal, 0, False, False)
        return True

    def passMuonID(self, tree, mu_id):
        """
            pass muon selections
        """
        mu_pt = tree.mu_track_pt[mu_id]
        res = mu_pt > SUBLEADING_MUON_PT
        mu_eta = tree.mu_track_eta[mu_id]
        res = res and abs(mu_eta) < 2.5

        if self.do_13TeV:
            d0_sig = tree.mu_d0_sig[mu_id]
            z0_ = tree.mu_z0_sintheta[mu_id]
        #elif hasattr(tree, "mu_d0_pv_sig"):
        #    d0_sig = tree.mu_d0_pv_sig[mu_id]
        #    z0_ = tree.mu_z0_pv_sintheta[mu_id]
        else:
            return True

        res = res and abs(d0_sig) < 6
        res = res and abs(z0_) < 1

        if self.m_debug:
            print "d0: ", d0_sig
            print "z0: ", z0_
            print "passed: ", res

        return res

    def passOniaCuts(self, tree, onia_id):
        """
            pass basic onia selections
        """
        #mu_id1 = tree.onia_id1[id_onia2]
        #mu_id2 = tree.onia_id2[id_onia2]
        #if tree.mu_charge[mu_id1] + tree.mu_charge[mu_id2] != 0:
        #    return False
        if hasattr(tree, "onia_type") and tree.onia_type[onia_id]==1:
            return False

        res = tree.onia_chi2[onia_id] < CHI2_DIMUON_CUT
        mass_onia = tree.onia_fitted_mass[onia_id]
        res = res and mass_onia > M34CUT_LOW and mass_onia < M34CUT_HI
        return res

    def passUpsilon(self, tree, onia_id):
        #if not self.passOniaCuts(tree, onia_id):
        #    return False

        if not self.is_neutral_onia(tree, onia_id):
            return False

        mass_onia = tree.onia_fitted_mass[onia_id]
        return (mass_onia > 9.2E3 and mass_onia < 9.7E3)

    def is_neutral_onia(self, tree, onia_id):
        """
        neutral onia, with muon pT > 4 GeV
        """

        mu_id1 = tree.onia_id1[onia_id]
        mu_id2 = tree.onia_id2[onia_id]
        if tree.mu_charge[mu_id1] + tree.mu_charge[mu_id2] != 0:
            return False

        mu_pt1 = tree.mu_track_pt[mu_id1]
        mu_pt2 = tree.mu_track_pt[mu_id2]
        if not self.use_new_pT and (mu_pt1 <= 4E3 or mu_pt2 <= 4E3):
            return False

        return True

    def notZero(self, x):
        return x!=0

    def fill_cut_flow(self, ncut):
        self.cut_flow.Fill(ncut)
        if self.run[0] >= 307619:
            self.cut_flow_with_3mu4.Fill(ncut)
        else:
            self.cut_flow_no3mu4.Fill(ncut)

    def fill_onia(self, tree, good_muons):
        """
         Fill in tree_onia, for chi2 performance studies
         Use all the onia that passed the criteria
        """
        if self.m_debug:
            print "in fill_onia"
            print "good muons:"+",".join([str(x) for x in good_muons])
        used_muons = []
        #for i in sorted(range(len(tree.onia_chi2)), key=lambda k:tree.onia_chi2[k], reverse=False):

        for i in range(len(tree.onia_chi2)):
            mass_onia = tree.onia_fitted_mass[i]
            if self.m_debug:
                print "read: {:.0f} {:.2f} {:.2f}".format(i, tree.onia_chi2[i],mass_onia)
            mu_id1 = tree.onia_id1[i]
            mu_id2 = tree.onia_id2[i]
            if self.m_debug:
                print "muon id:",mu_id1, mu_id2

            if mu_id1 not in good_muons or\
               mu_id2 not in good_muons:
                continue

            # charge
            if tree.mu_charge[mu_id1] + tree.mu_charge[mu_id2] != 0:
                continue
            if self.m_debug:
                print "passed charge"
            # pT
            mu_pt1 = tree.mu_track_pt[mu_id1]
            mu_pt2 = tree.mu_track_pt[mu_id2]
            if mu_pt1 <= SUBLEADING_MUON_PT or mu_pt2 <= SUBLEADING_MUON_PT:
                continue
            if self.m_debug:
                print "passed pT"
            # mass
            if mass_onia < M34CUT_LOW or mass_onia > M34CUT_HI:
                continue

            if self.m_debug:
                print "passed mass"

            if self.m_debug:
                print "onia: ", i, "accepted:", tree.onia_chi2[i]
            used_muons.append(mu_id1)
            used_muons.append(mu_id2)

            self.up_mass.push_back(mass_onia/1E3)
            self.up_pt.push_back(tree.onia_track_pt[i]/1E3)
            self.up_chi2.push_back(tree.onia_chi2[i])
            self.up_d0_1.push_back(tree.mu_d0_sig[mu_id1])
            self.up_d0_2.push_back(tree.mu_d0_sig[mu_id2])
            self.up_z0_1.push_back(tree.mu_z0_sintheta[mu_id1])
            self.up_z0_2.push_back(tree.mu_z0_sintheta[mu_id2])

        #self.tree_onia.Fill()

    def find_onia_pair(self, tree, good_muons, pass_cuts):
        """
        return the onia-pair index build from good_muons and pass cuts.
        """
        onia_pair_index = []
        for i in range(len(tree.onia_fitted_mass)):
            mu_id1 = tree.onia_id1[i]
            mu_id2 = tree.onia_id2[i]
            if mu_id1 not in good_muons or\
               mu_id2 not in good_muons:
                continue

            if not pass_cuts(tree, i):
                continue

            for j in range(i+1, len(tree.onia_fitted_mass)):
                mu_id3 = tree.onia_id1[j]
                mu_id4 = tree.onia_id2[j]
                if mu_id3 not in good_muons or\
                   mu_id4 not in good_muons:
                    continue

                if mu_id3 in [mu_id1, mu_id2] or\
                   mu_id4 in [mu_id1, mu_id2]:
                    continue

                if not pass_cuts(tree, j):
                    continue

                onia_pair_index.append( (i,j) )
        return onia_pair_index

    def print_event(self, tree):
        self.out_events += "{:.0f} {:.0f}\n".format(self.run[0], self.event[0])

    @staticmethod
    def is_author(author_val):
        """ check the type of the author,
        return a list of its authors."""
        has_author = False
        authors_list = []
        for i,key in enumerate( sorted(author_cuts.keys()) ):
            if author_val&key:
                authors_list.append(i)
                has_author = True

        if has_author:
            return authors_list
        else:
            return [-1]

    def get_muon_type(self, tree, mu_id):
        """
        define muon type as the following:
            ST muons: author = 4 or 6
            Combined muons: 1 or 2
        """
        author = tree.mu_author[mu_id]
        if author == 1 or author == 2:
            return 0
        elif author == 4 or author == 6:
            return 2
        else:
            print "I don't know this author:", author


def draw(file_name, post_fix):
    f1 = ROOT.TFile.Open(file_name)

    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)

    x_pos = 0.70
    y_pos = 0.85
    n_rebin = 1

    upsilon_low = 8
    upsilon_hi = 12
    h1 = f1.Get("h_upsilon")
    h1.Rebin(n_rebin)
    h1.GetYaxis().SetTitle("Events / {:.0f} MeV".format(50*n_rebin))
    h1.GetXaxis().SetRangeUser(upsilon_low, upsilon_hi)
    h1.Draw("EP")
    ROOT.myText(x_pos, y_pos, 1, "Entries: {:.0f}".format(h1.Integral()))
    canvas.SaveAs("mUpsilon_"+post_fix+".pdf")

    jspi_low = 1
    jspi_hi = 15
    #orignal
    h1.GetXaxis().SetRangeUser(jspi_low, jspi_hi)
    h1.Draw("EP")
    ROOT.myText(x_pos, y_pos, 1, "Entries: {:.0f}".format(h1.Integral()))
    canvas.SaveAs("mUpsilon_Full_"+post_fix+".pdf")

    h6 = f1.Get("h_jpsi")
    h6.Rebin(n_rebin)
    h6.GetYaxis().SetTitle("Events / {:.0f} MeV".format(100*n_rebin))
    h6.GetXaxis().SetRangeUser(jspi_low, jspi_hi)
    h6.Draw("EP")
    ROOT.myText(x_pos, y_pos, 1, "Entries: {:.0f}".format(h6.Integral()))
    canvas.SaveAs("mJpsi_"+post_fix+".pdf")

    m4l_bin_width= 50
    n_rebin_4l = 4
    m4l_xlow = 15
    m4l_xhi = 25
    m4l_yaxis = m4l_bin_width * n_rebin_4l

    h2 = f1.Get("h_m4l_LHS")
    h2.Rebin(n_rebin_4l)
    h2.GetYaxis().SetTitle("Events / {:.0f} MeV".format(m4l_yaxis))
    h2.GetXaxis().SetRangeUser(m4l_xlow, m4l_xhi)
    h2.Draw("EP")
    ROOT.myText(x_pos, y_pos, 1, "Entries: {:.0f}".format(h2.Integral()))
    canvas.SaveAs("m4l_LHS_"+post_fix+".pdf")

    h3 = f1.Get("h_m4l_RHS_U1")
    h3.Rebin(n_rebin_4l)
    h3.GetYaxis().SetTitle("Events / {:.0f} MeV".format(m4l_yaxis))
    h3.GetXaxis().SetRangeUser(m4l_xlow, m4l_xhi)
    h3.Draw("EP")
    ROOT.myText(x_pos, y_pos, 1, "Entries: {:.0f}".format(h3.Integral()))
    canvas.SaveAs("m4l_RHS_U1_"+post_fix+".pdf")

    h4 = f1.Get("h_m4l_RHS_U2")
    h4.Rebin(n_rebin_4l)
    h4.GetYaxis().SetTitle("Events / {:.0f} MeV".format(m4l_yaxis))
    h4.GetXaxis().SetRangeUser(m4l_xlow, m4l_xhi)
    h4.Draw("EP")
    ROOT.myText(x_pos, y_pos, 1, "Entries: {:.0f}".format(h4.Integral()))
    canvas.SaveAs("m4l_RHS_U2_"+post_fix+".pdf")

    h5 = f1.Get("h_m4l")
    h5.Rebin(n_rebin_4l)
    h5.GetYaxis().SetTitle("Events / {:.0f} MeV".format(m4l_yaxis))
    h5.GetXaxis().SetRangeUser(m4l_xlow, m4l_xhi)
    h5.Draw("EP")
    ROOT.myText(x_pos, y_pos, 1, "Entries: {:.0f}".format(h5.Integral()))
    canvas.SaveAs("m4l_"+post_fix+".pdf")

def compare_sideband(file_name, br_name="x_chi2"):

    f1 = ROOT.TFile.Open(file_name)
    tree = f1.Get("bls")
    h_sig = ROOT.TH1F("h_sig", "signal", 10, 0, 50)
    h_left = h_sig.Clone("h_left")
    h_right = h_sig.Clone("h_right")

    sig_cut = "m4l_track > 17.8 && m4l_track < 18.8"
    left_cut = "m4l_track < 17.8"
    right_cut = "m4l_track > 18.8 && m4l_track < 25"

    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)

    tree.Draw(br_name+">>h_sig", br_name+" < 100 && "+sig_cut)
    tree.Draw(br_name+">>h_left", br_name+" < 100 && "+left_cut)
    tree.Draw(br_name+">>h_right", br_name+" < 100 && "+right_cut)

    n_sig = h_sig.Integral()
    h_left.Scale(n_sig / h_left.Integral())
    h_right.Scale(n_sig / h_right.Integral())

    #canvas.SetLogy()
    #canvas.SetLogx()
    h_sig.SetMarkerColor(1)
    h_sig.Draw("EP")
    h_left.SetLineColor(2)
    h_right.SetLineColor(4)
    h_left.SetMarkerColor(2)
    h_right.SetMarkerColor(4)
    h_left.SetMarkerStyle(24)
    h_right.SetMarkerStyle(27)

    h_left.Draw("same EP")
    h_right.Draw("same EP")
    legend = ROOT.myLegend(0.7, 0.8, 0.9, 0.9)
    legend.AddEntry(h_sig, "SR", "EP")
    legend.AddEntry(h_left, "LHS", "EP")
    legend.AddEntry(h_right, "RHS", "EP")
    legend.Draw("same")
    canvas.SaveAs("cmp_"+br_name+".pdf")

if __name__ == "__main__":

    usage = sys.argv[0]+" make/draw/cmp file_index,file2 out_name"
    parser = OptionParser(usage=usage, description="apply final cuts for upsilon analysis")
    parser.add_option("--do8TeV", action="store_true", dest="do8TeV",help="Perform 8 TeV analysis",default=False)
    parser.add_option("--uf", action="store_true", dest="uf",help="select upsilon first",default=False)
    parser.add_option("-v", action="store_true", dest="verbose",help="debug mode",default=False)
    parser.add_option("--oldMuonType", action="store_true", dest="oldMuonType", help="make: use old defition of muon type", default=False)
    parser.add_option("--oldPtCut", action="store_true", dest="oldPtCut", help="make: use old pT cut", default=False)

    (options, args) = parser.parse_args()
    if len(args) < 3:
        parser.print_help()
        exit(0)

    option = args[0]
    out_name = args[2]
    print option
    print out_name

    import AtlasStyle

    if option == "make":
        input_name = args[1].split(',')
        bls_ana = BLSana(out_name)
        if options.do8TeV:
            bls_ana.do_13TeV = False

        if options.uf:
            bls_ana.upsilon_first = True

        if options.verbose:
            bls_ana.m_debug = True

        if options.oldMuonType:
            bls_ana.use_new_type = False

        if options.oldPtCut:
            bls_ana.use_new_pT = False

        bls_ana.book_hists()
        bls_ana.book_tree()
        bls_ana.book_upsilon()
        bls_ana.fill_hists(input_name)
        bls_ana.save_hists()

    elif option == "draw":
        input_name = args[1]
        draw(input_name, out_name)
    else:
        input_name = args[1]
        compare_sideband(input_name, out_name)
