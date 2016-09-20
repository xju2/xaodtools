#!/usr/bin/env python

import ROOT
import sys
ROOT.gROOT.SetBatch()
import string
import math
from array import array

import AtlasStyle
if not hasattr(ROOT, "myText"):
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/AtlasUtils.C")
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c")

m_debug = False
CHI2_DIMUON_CUT = 6
CHI2_4MUON_CUT = 50
LEADING_MUON_PT = 3E3
SUBLEADING_MUON_PT = 3E3

class BLSana:
    """
    select upsilon+upsilon
    """
    def __init__(self, out_name):
        print "start to work"
        self.out_name = out_name
        self.do_save = True
        self.ignore_PV_cut = True

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

        self.u1_chi2 = array('f', [0])
        self.u2_chi2 = array('f', [0])

        self.run = array('i', [0])
        self.event = array('i', [0])


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

        self.out_tree.Branch("u1_chi2", self.u1_chi2, "u1_chi2/F")
        self.out_tree.Branch("u2_chi2", self.u2_chi2, "u2_chi2/F")

        self.out_tree.Branch("run", self.run, "run/I")
        self.out_tree.Branch("event", self.event, "event/I")

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
        tree = ROOT.TChain("physics", "physics")
        for file_ in file_names:
            print "adding.. ", file_
            tree.Add(file_)

        nentries = tree.GetEntries()
        print "total entries:",nentries
        if nentries == 0:
            self.do_save = False
            return

        ## disable branches
        tree.SetBranchStatus("*", 0)
        tree.SetBranchStatus("onia_type", 1)
        tree.SetBranchStatus("onia_fitted_mass", 1)
        tree.SetBranchStatus("onia_chi2", 1)
        tree.SetBranchStatus("onia_id1", 1)
        tree.SetBranchStatus("onia_id2", 1)

        tree.SetBranchStatus("quad_type", 1)
        tree.SetBranchStatus("quad_fitted_mass", 1)
        tree.SetBranchStatus("quad_fitted_pt", 1)
        tree.SetBranchStatus("quad_track_mass", 1)
        tree.SetBranchStatus("quad_track_pt", 1)
        tree.SetBranchStatus("quad_mass", 1)
        tree.SetBranchStatus("quad_nCombined", 1)

        tree.SetBranchStatus("n_muon", 1)
        tree.SetBranchStatus("mu_track_pt", 1)
        tree.SetBranchStatus("mu_track_eta", 1)
        tree.SetBranchStatus("mu_pt", 1)
        tree.SetBranchStatus("mu_pvID", 1)

        tree.SetBranchStatus("onia_x", 1)
        tree.SetBranchStatus("onia_y", 1)
        tree.SetBranchStatus("onia_z", 1)
        tree.SetBranchStatus("quad_x", 1)
        tree.SetBranchStatus("quad_y", 1)
        tree.SetBranchStatus("quad_z", 1)
        tree.SetBranchStatus("quad_chi2", 1)
        tree.SetBranchStatus("quad_id1", 1)
        tree.SetBranchStatus("quad_id2", 1)
        tree.SetBranchStatus("quad_id3", 1)
        tree.SetBranchStatus("quad_id4", 1)
        tree.SetBranchStatus("Event", 1)
        tree.SetBranchStatus("Run", 1)
        tree.SetBranchStatus("EventNumber", 1)
        tree.SetBranchStatus("RunNumber", 1)

        tree.SetBranchStatus("v0_x", 1)
        tree.SetBranchStatus("v0_y", 1)
        tree.SetBranchStatus("v0_z", 1)

        interested_event = [
            #(208662, 168839095),
            (202712, 14264115),
            (205071, 120286323),
                           ]

        for ientry in xrange(nentries):
            tree.GetEntry(ientry)

            #if m_debug and (tree.Run != 203636 or tree.Event != 50171978):
            #if m_debug and (tree.Run != 202712 or tree.Event != 4744031):
            #if m_debug and (tree.Run != 202712 or tree.Event != 42230077):

            if hasattr(tree, "Run"):
                run = tree.Run
                event = tree.Event
            else:
                run = tree.RunNumber
                event = tree.EventNumber

            if m_debug and (run, event) not in interested_event:
                continue

            if m_debug:
                print "process: ",run, event

            if tree.n_muon != 4:
                continue

            self.run[0] = run
            self.event[0] = event
            ##
            results = self.select_upsilon(tree)
            #results = self.select_quadruplet(tree)
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

            self.m1_track_pt[0] = tree.mu_track_pt[ tree.quad_id1[quad_id] ]
            self.m2_track_pt[0] = tree.mu_track_pt[ tree.quad_id2[quad_id] ]
            self.m3_track_pt[0] = tree.mu_track_pt[ tree.quad_id3[quad_id] ]
            self.m4_track_pt[0] = tree.mu_track_pt[ tree.quad_id4[quad_id] ]

            self.m1_track_eta[0] = tree.mu_track_eta[ tree.quad_id1[quad_id] ]
            self.m2_track_eta[0] = tree.mu_track_eta[ tree.quad_id2[quad_id] ]
            self.m3_track_eta[0] = tree.mu_track_eta[ tree.quad_id3[quad_id] ]
            self.m4_track_eta[0] = tree.mu_track_eta[ tree.quad_id4[quad_id] ]

            self.m1_pt[0] = tree.mu_pt[ tree.quad_id1[quad_id] ]
            self.m2_pt[0] = tree.mu_pt[ tree.quad_id2[quad_id] ]
            self.m3_pt[0] = tree.mu_pt[ tree.quad_id3[quad_id] ]
            self.m4_pt[0] = tree.mu_pt[ tree.quad_id4[quad_id] ]

            self.u1_chi2[0] = tree.onia_chi2[ onia1_id ]
            self.u2_chi2[0] = tree.onia_chi2[ onia2_id ]

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

            if tree.mu_track_pt[id1] > LEADING_MUON_PT and\
               tree.mu_track_pt[id2] > LEADING_MUON_PT and\
               tree.onia_chi2[i] < CHI2_DIMUON_CUT and\
               onia_mass > 9.2E3 and\
               onia_mass < 9.7E3 and\
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

        if not has_upsilon:
            if m_debug:
                print "cannot find upsilon"
            return None

        if m_debug:
            print "Upsilon mass: ", mU

        has_add_mu = False
        m34 = -1
        onia2_chi2 = 99999
        for i,onia_mass in enumerate(tree.onia_fitted_mass):
            id1 = tree.onia_id1[i]
            id2 = tree.onia_id2[i]
            if i== onia1_id:
                continue

            if id1 in can_id or id2 in can_id:
                continue

            if tree.mu_track_pt[id1] > SUBLEADING_MUON_PT and\
               tree.mu_track_pt[id2] > SUBLEADING_MUON_PT and\
               tree.onia_chi2[i] < CHI2_DIMUON_CUT and\
               onia_mass > 9.2E3 and\
               onia_mass < 9.7E3 and\
               tree.onia_chi2[i] < onia2_chi2:

                has_add_mu = True
                onia2_chi2 = tree.onia_chi2[i]
                m34 = onia_mass/1E3
                onia2_id = i
                if len(can_id) < 4:
                    can_id.append(id1)
                    can_id.append(id2)
                else:
                    can_id[2] = id1
                    can_id[3] = id2

        if not has_add_mu:
            if m_debug:
                print "cannot second pair"
            return None

        if m_debug:
            print "muonID:"," ".join([str(x) for x in can_id])

        mass_id = -1
        quad_id = -1
        j = -1
        for i,x_chi2 in enumerate(tree.quad_chi2):
            if x_chi2 < 0:
                if m_debug:
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
                if m_debug:
                    print "less than 3 combined muons"
                continue

            ## vertex association cut
            pvID1 = tree.mu_pvID[can_id[0]]
            pvID2 = tree.mu_pvID[can_id[1]]
            pvID3 = tree.mu_pvID[can_id[2]]
            pvID4 = tree.mu_pvID[can_id[3]]
            pass_vertex = (pvID1 == pvID2) and (pvID3 == pvID2 or pvID4 == pvID2)
            #pass_vertex = pvID1 == 0 and pvID1 == pvID2 and pvID3 == pvID2 and pvID4 == pvID3
            if not self.ignore_PV_cut and not pass_vertex:
                if m_debug:
                    print "not from same vertex",pvID1,pvID2,pvID3,pvID4
                continue

            mass_id = j
            quad_mass = tree.quad_fitted_mass[mass_id]
            if m_debug:
                print "mass-> ",quad_mass,mass_id
            quad_id = i

        if quad_id < 0:
            return None

        if m_debug:
            print "quadID:", quad_id, "mass:", quad_mass," massID:", mass_id

        return (mU, m34, mass_id, onia1_id, onia2_id, quad_id)

    def select_quadruplet(self, tree):

        m12 = [-1, False]
        m34 = [-1, False]
        mass_id = -1
        onia1_id = -1
        onia2_id = -1
        quad_id =  -1
        min_quad_chi2 = 999999
        has_upsilon = False
        has_jpsi = False
        for id_onia1 in range(len(tree.onia_fitted_mass)):
            mu_id1 = tree.onia_id1[id_onia1]
            mu_id2 = tree.onia_id2[id_onia1]

            mu_pt1 = tree.mu_track_pt[mu_id1]
            mu_pt2 = tree.mu_track_pt[mu_id2]

            if tree.onia_chi2[id_onia1] > CHI2_DIMUON_CUT:
                continue

            if mu_pt1 <= 3E3 or mu_pt2 <= 3E3:
                continue

            mass_onia1 = tree.onia_fitted_mass[id_onia1]
            if mass_onia1 < 9.2E3 or\
               mass_onia1 > 9.7E3:
                continue
            if m_debug:
                print "find a good onia"

            for id_onia2 in range(id_onia1+1, len(tree.onia_fitted_mass)):
                mu_id3 = tree.onia_id1[id_onia2]
                mu_id4 = tree.onia_id2[id_onia2]

                mu_pt3 = tree.mu_track_pt[mu_id3]
                mu_pt4 = tree.mu_track_pt[mu_id4]
                if tree.onia_chi2[id_onia2] > CHI2_DIMUON_CUT:
                    continue

                if mu_pt3 <= 3E3 or mu_pt4 <= 3E3:
                    continue

                mass_onia2 = tree.onia_fitted_mass[id_onia2]
                if mass_onia2 < 9.2E3 or\
                   mass_onia2 > 9.7E3:
                    continue

                if mu_id3 in [mu_id1, mu_id2] or\
                   mu_id4 in [mu_id1, mu_id2]:
                    continue
                if m_debug:
                    print "find a second onia"

                # pT cut, either 44,33
                #assign mU and m34
                m12[0] = mass_onia1
                m12[1] = mu_pt1 > 4E3 and mu_pt2 > 4E3
                m34[0] = mass_onia2
                m34[1] = mu_pt3 > 4E3 and mu_pt4 > 4E3
                #if not m12[1] and not m34[1]:
                #    continue

                if m_debug:
                    print "find onia pass 4GeV cut"
                # find the chi2 of the four muons
                muon_cans = [mu_id1, mu_id2, mu_id3, mu_id4]
                if m_debug:
                    print "muonID:"," ".join([str(x) for x in muon_cans])
                j = -1
                for iquad, x_chi2 in enumerate(tree.quad_chi2):
                    if x_chi2 < 0:
                        continue

                    j += 1
                    id1 = tree.quad_id1[iquad]
                    id2 = tree.quad_id2[iquad]
                    id3 = tree.quad_id3[iquad]
                    id4 = tree.quad_id4[iquad]
                    if m_debug:
                        print id1,id2,id3,id4
                    if id1 not in muon_cans or\
                       id2 not in muon_cans or\
                       id3 not in muon_cans or\
                       id4 not in muon_cans:
                        continue

                    if tree.quad_nCombined[iquad] < 3:
                        if m_debug:
                            print "less than 3 combined muons",tree.quad_nCombined[iquad]
                        continue

                    ## vertex association cut
                    pvID1 = tree.mu_pvID[muon_cans[0]]
                    pvID2 = tree.mu_pvID[muon_cans[1]]
                    pvID3 = tree.mu_pvID[muon_cans[2]]
                    pvID4 = tree.mu_pvID[muon_cans[3]]
                    pass_vertex = (pvID1 == pvID2) and (pvID3 == pvID2 or pvID4 == pvID2)
                    if not pass_vertex and not self.ignore_PV_cut:
                        if m_debug:
                            print "not from same vertex",pvID1,pvID2,pvID3,pvID4
                        continue

                    # min chi2
                    if tree.quad_chi2[iquad] > min_quad_chi2:
                        if m_debug:
                            print "too large chi2"
                        continue

                    mass_id = j
                    quad_id = iquad
                    min_quad_chi2 = tree.quad_chi2[iquad]
                    onia1_id = id_onia1
                    onia2_id = id_onia2

        # check if there's quadruplet
        if quad_id < 0:
            if m_debug:
                print self.run[0],self.event[0],"no quadruplet"
            return None

        return (m12[0], m34[0], mass_id, onia1_id, onia2_id, quad_id)
        # apply the upsilon cut!
        #if (m12[1] and m12[0] > 9.2E3 and m12[0] < 9.7E3) or\
        #   (m34[1] and m34[0] > 9.2E3 and m34[0] < 9.7E3):
        #    return (m12[0], m34[0], mass_id, onia1_id, onia2_id, quad_id)
        #else:
        #    if m_debug:
        #        print self.run[0],self.event[0],"no upsilon"
        #    return None

    def get_dis(self, x1, y1, z1, x2, y2, z2):
        return math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

    def save_hists(self):
        if not self.do_save:
            return

        # save output
        fout = ROOT.TFile.Open(out_name, "recreate")
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
        fout.Close()


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
    if len(sys.argv) < 4:
        print sys.argv[0], "make/draw/cmp file_index,file2 out_name"
        exit(1)

    base_name = "/afs/cern.ch/user/x/xju/work/upsilon/run/data12_v1/split_and_merge/merged_"
    #input_files = [base_name+x+".root" for x in string.ascii_lowercase]
    option = sys.argv[1]
    out_name = sys.argv[3]
    print sys.argv[2]
    print out_name

    if option == "make":
        input_name = sys.argv[2].split(',')
        bls_ana = BLSana(out_name)
        bls_ana.book_hists()
        bls_ana.book_tree()
        bls_ana.fill_hists(input_name)
        bls_ana.save_hists()
    elif option == "draw":
        input_name = sys.argv[2]
        draw(input_name, out_name)
    else:
        input_name = sys.argv[2]
        compare_sideband(input_name, out_name)
