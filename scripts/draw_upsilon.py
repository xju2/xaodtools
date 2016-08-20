#!/usr/bin/env python

import ROOT
import sys
ROOT.gROOT.SetBatch()

import AtlasStyle
if not hasattr(ROOT, "myText"):
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/AtlasUtils.C")
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c")

def make_hists():
    tree = ROOT.TChain("physics", "physics")
    tree.Add("data16.merged.p2689.root")
    tree.Add("data16.merged.p2667.root")
    nentries = tree.GetEntries()
    print "total entries:",nentries

    h_upsilon = ROOT.TH1F("h_upsilon", "upsilon mass;m_{#varUpsilon} [GeV];Events / 50 MeV", 80, 8, 12)
    h_m4l = ROOT.TH1F("h_m4l", "m4l;m_{#varUpsilon+#mu+#mu} [GeV];Events / 200 MeV", 50, 14, 24)
    h_m4l_LHS = ROOT.TH1F("h_m4l_LHS", "m4l LHS;LHS m_{#varUpsilon+#mu+#mu} [GeV];Events / 200 MeV", 50, 14, 24)
    h_m4l_RHS_U1 = ROOT.TH1F("h_m4l_RHS_U1", "m4l RHS;RHS U1 m_{#varUpsilon+#mu+#mu} [GeV];Events / 200 MeV", 50, 14, 24)
    h_m4l_RHS_U2 = ROOT.TH1F("h_m4l_RHS_U2", "m4l RHS;RHS U2 m_{#varUpsilon+#mu+#mu} [GeV];Events / 200 MeV", 50, 14, 24)
    h_jpsi = ROOT.TH1F("h_jpsi", "jpsi mass;m_{#mu+#mu} [GeV];Events / 100 MeV", 100, 1, 11)

    ## disable branches
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("mUpsilon", 1)
    tree.SetBranchStatus("m4l", 1)
    tree.SetBranchStatus("n_muon", 1)
    tree.SetBranchStatus("mu_p4", 1)

    for ientry in xrange(nentries):
        tree.GetEntry(ientry)
        mU = tree.mUpsilon/1E3
        m4l = tree.m4l/1E3
        if tree.n_muon != 4:
            continue
        # require a jpsi
        mu_id1 = -1
        mu_id2 = -1
        for i in range(4):
            p4_1 = tree.mu_p4[i]
            for j in range(4):
                if j <= i:
                    continue
                p4_2 = tree.mu_p4[j]
                inv_mass = (p4_1+p4_2).M()
                if inv_mass > 8E3 and inv_mass < 12E3:
                    mU = inv_mass/1E3
                    mu_id1 = i
                    mu_id2 = j

        mu_tlv = ROOT.TLorentzVector()
        for i in range(4):
            if i != mu_id1 and i != mu_id2:
                mu_tlv += tree.mu_p4[i]
        m_mu2 = mu_tlv.M()/1E3
        h_jpsi.Fill(m_mu2)
        h_upsilon.Fill(mU)

        if m_mu2 < 2.5 or m_mu2 > 3.5:
            continue

        if tree.n_muon > 3:
            if mU < 9:
                h_m4l_LHS.Fill(m4l)
            elif mU >= 9 and mU < 9.8:
                h_m4l.Fill(m4l)
            elif mU >= 9.8 and mU < 10.8:
                h_m4l_RHS_U1.Fill(m4l)
            else:
                h_m4l_RHS_U2.Fill(m4l)

    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    h_upsilon.Draw()
    canvas.SaveAs("mUpsilon.pdf")
    h_m4l_LHS.Draw()
    canvas.SaveAs("m4l_LHS.pdf")
    h_m4l_RHS_U1.Draw()
    canvas.SaveAs("m4l_RHS_U1.pdf")
    h_m4l_RHS_U2.Draw()
    canvas.SaveAs("m4l_RHS_U2.pdf")
    h_m4l.Draw()
    canvas.SaveAs("m4l.pdf")
    h_jpsi.Draw()
    canvas.SaveAs("mJpsi.pdf")

    # save output
    fout = ROOT.TFile.Open("hist_merged.root", "recreate")
    h_upsilon.Write()
    h_m4l.Write()
    h_m4l_LHS.Write()
    h_m4l_RHS_U1.Write()
    h_m4l_RHS_U2.Write()
    h_jpsi.Write()
    fout.Close()

    #f1.Close()

def draw(file_name):
    f1 = ROOT.TFile.Open(file_name)
    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    f1.Get("h_upsilon").Draw()
    canvas.SaveAs("mUpsilon.pdf")
    f1.Get("h_m4l_LHS").Draw()
    canvas.SaveAs("m4l_LHS.pdf")
    f1.Get("h_m4l_RHS_U1").Draw()
    canvas.SaveAs("m4l_RHS_U1.pdf")
    f1.Get("h_m4l_RHS_U2").Draw()
    canvas.SaveAs("m4l_RHS_U2.pdf")
    f1.Get("h_m4l").Draw()
    canvas.SaveAs("m4l.pdf")

if __name__ == "__main__":
    #file_name = "data16.merged.root"
    #make_hists("data16.merged.p2689.root")
    make_hists()
    #draw("hist.merged.root")
