#!/usr/bin/env python

import ROOT
ROOT.gROOT.SetBatch()

import AtlasStyle
if not hasattr(ROOT, "myText"):
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/AtlasUtils.C")
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c")

def draw():
    file_name = "data16.merged.root"
    f1 = ROOT.TFile.Open(file_name)
    tree = f1.Get("physics")
    nentries = tree.GetEntries()
    print "total entries:",nentries

    h_upsilon = ROOT.TH1F("h_upsilon", "upsilon mass;m_{#varUpsilon} [GeV];Events / 50 MeV", 80, 8, 12)
    h_m4l = ROOT.TH1F("h_m4l", "m4l;m_{#varUpsilon+#mu+#mu} [GeV];Events / 200 MeV", 50, 14, 24)

    ## disable branches
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("mUpsilon", 1)
    tree.SetBranchStatus("m4l", 1)
    tree.SetBranchStatus("n_muon", 1)
    for ientry in xrange(nentries):
        tree.GetEntry(ientry)
        h_upsilon.Fill(tree.mUpsilon/1E3)
        if tree.n_muon > 3:
            h_m4l.Fill(tree.m4l/1E3)

    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    h_upsilon.Draw()
    canvas.SaveAs("mUpsilon.pdf")
    h_m4l.Draw()
    canvas.SaveAs("m4l.pdf")
    f1.Close()

if __name__ == "__main__":
    draw()
