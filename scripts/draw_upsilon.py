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

class BLSana:
    def __init__(self, out_name):
        print "start to work"
        self.out_name = out_name

    def book_tree(self):
        self.out_tree = ROOT.TTree("bls", "bls")

        self.m4l = array('f', [0])
        self.m4l_track = array('f', [0])
        self.m4l_fitted = array('f', [0])
        self.m12 = array('f', [0])
        self.m34 = array('f', [0])
        self.dis_v1_x = array('f', [0])
        self.dis_v2_x = array('f', [0])
        self.x_dis = array('f', [0])
        self.x_chi2 = array('f', [0])

        self.out_tree.Branch("m4l_fitted", self.m4l_fitted, "m4l_fitted/F")
        self.out_tree.Branch("m4l_track", self.m4l_track, "m4l_track/F")
        self.out_tree.Branch("m4l", self.m4l, "m4l/F")
        self.out_tree.Branch("m12", self.m12, "m12/F")
        self.out_tree.Branch("m34", self.m34, "m34/F")
        self.out_tree.Branch("dis_v1_x", self.dis_v1_x, "dis_v1_x/F")
        self.out_tree.Branch("dis_v2_x", self.dis_v2_x, "dis_v2_x/F")
        self.out_tree.Branch("x_dis", self.x_dis, "x_dis/F")
        self.out_tree.Branch("x_chi2", self.x_chi2, "x_chi2/F")

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
            tree.Add(file_)

        nentries = tree.GetEntries()
        print "total entries:",nentries

        ## disable branches
        tree.SetBranchStatus("*", 0)
        tree.SetBranchStatus("onia_fitted_mass", 1);
        tree.SetBranchStatus("onia_chi2", 1);
        tree.SetBranchStatus("onia_id1", 1);
        tree.SetBranchStatus("onia_id2", 1);
        tree.SetBranchStatus("quad_fitted_mass", 1);
        tree.SetBranchStatus("quad_track_mass", 1);
        tree.SetBranchStatus("quad_mass", 1);
        tree.SetBranchStatus("mu_track_pt", 1);
        tree.SetBranchStatus("n_muon", 1);

        tree.SetBranchStatus("onia_x", 1);
        tree.SetBranchStatus("onia_y", 1);
        tree.SetBranchStatus("onia_z", 1);
        tree.SetBranchStatus("quad_x", 1);
        tree.SetBranchStatus("quad_y", 1);
        tree.SetBranchStatus("quad_z", 1);
        tree.SetBranchStatus("quad_chi2", 1);
        tree.SetBranchStatus("quad_id1", 1);
        tree.SetBranchStatus("quad_id2", 1);
        tree.SetBranchStatus("quad_id3", 1);
        tree.SetBranchStatus("quad_id4", 1);
        tree.SetBranchStatus("Event", 1);

        tree.SetBranchStatus("mu_pvID", 1);
        tree.SetBranchStatus("Event", 1);

        for ientry in xrange(nentries):
            tree.GetEntry(ientry)

            if tree.n_muon < 4:
                continue
            ##
            has_upsilon = False
            mU = -1
            can_id = []
            onia1_id = 0
            upsilon_chi2 = 9999
            for i,onia_mass in enumerate(tree.onia_fitted_mass):
                id1 = tree.onia_id1[i]
                id2 = tree.onia_id2[i]
                #print id1, id2
                if tree.mu_track_pt[id1] > 4E3 and\
                   tree.mu_track_pt[id2] > 4E3 and\
                   tree.onia_chi2[id1] < 3 and\
                   tree.onia_chi2[id2] < 3 and\
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
                continue

            has_add_mu = False
            m34 = -1
            onia2_chi2 = 99999
            for i,onia_mass in enumerate(tree.onia_fitted_mass):
                id1 = tree.onia_id1[i]
                id2 = tree.onia_id2[i]
                if id1 in can_id or id2 in can_id:
                    continue

                if tree.mu_track_pt[id1] > 3E3 and\
                   tree.mu_track_pt[id2] > 3E3 and\
                   tree.onia_chi2[id1] < 3 and\
                   tree.onia_chi2[id2] < 3 and\
                   onia_mass > 2E3 and\
                   onia_mass < 20E3 and\
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
                continue

            quad_id = -1
            if len(can_id) != 4:
                print "ERROR: ", tree.Event

            for i,quad_mass in enumerate(tree.quad_fitted_mass):
                id1 = tree.quad_id1[i]
                id2 = tree.quad_id2[i]
                id3 = tree.quad_id3[i]
                id4 = tree.quad_id4[i]
                if id1 not in can_id or\
                   id2 not in can_id or\
                   id3 not in can_id or\
                   id4 not in can_id or\
                   quad_mass > 50E3:
                    continue

                ## vertex association cut
                pvID1 = tree.mu_pvID[can_id[0]]
                pvID2 = tree.mu_pvID[can_id[1]]
                pvID3 = tree.mu_pvID[can_id[2]]
                pvID4 = tree.mu_pvID[can_id[3]]
                pass_vertex = (pvID1 == pvID2) and (pvID3 == pvID2 or pvID4 == pvID2)
                if not pass_vertex:
                    continue

                quad_id = i

            if quad_id < 0:
                continue
                #print tree.Event

            self.h_upsilon.Fill(mU)
            self.h_jpsi.Fill(m34)
            self.h_m4l.Fill(tree.quad_fitted_mass[quad_id]/1E3)

            onia1_x = tree.onia_x[onia1_id]
            onia1_y = tree.onia_y[onia1_id]
            onia1_z = tree.onia_z[onia1_id]
            onia2_x = tree.onia_x[onia2_id]
            onia2_y = tree.onia_y[onia2_id]
            onia2_z = tree.onia_z[onia2_id]
            quad_x = tree.quad_x[quad_id]
            quad_y = tree.quad_y[quad_id]
            quad_z = tree.quad_z[quad_id]

            dis_v1_quad = self.get_dis(onia1_x, onia1_y, onia1_z, quad_x, quad_y, quad_z)
            dis_v2_quad = self.get_dis(onia2_x, onia2_y, onia2_z, quad_x, quad_y, quad_z)
            dis_quad = math.sqrt(quad_x**2 + quad_y**2 + quad_z**2)
            chi_quad = tree.quad_chi2[quad_id]

            self.h_dis_v1_quad.Fill(dis_v1_quad)
            self.h_dis_v2_quad.Fill(dis_v2_quad)
            self.h_dis_quad.Fill(dis_quad)
            self.h_chi_quad.Fill(chi_quad)

            self.m4l_fitted[0] =  tree.quad_fitted_mass[quad_id]/1E3
            self.m4l_track[0] =  tree.quad_track_mass[quad_id]/1E3
            self.m4l[0] =  tree.quad_mass[quad_id]/1E3
            self.m12[0] = mU
            self.m34[0] = m34
            self.dis_v1_x[0] = dis_v1_quad
            self.dis_v2_x[0] = dis_v2_quad
            self.x_dis[0] = dis_quad
            self.x_chi2[0] = chi_quad
            self.out_tree.Fill()



    def get_dis(self, x1, y1, z1, x2, y2, z2):
        return math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

    def save_hists(self):
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


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print sys.argv[0], "make/draw file_index,file2 out_name"
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
    else:
        input_name = sys.argv[2]
        draw(input_name, out_name)
