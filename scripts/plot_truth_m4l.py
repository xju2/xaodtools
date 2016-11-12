#!/usr/bin/env python
"""
read ntuple produced by Truth_JETMET...
make some validation plots
"""

import ROOT
ROOT.gROOT.SetBatch()

def make_plots(file_name):
    f1 = ROOT.TFile.Open(file_name)
    tree = f1.Get("physics")

    h_m4l = ROOT.TH1F("h_m4l", "m4l;m_{4l} [GeV]", 80, 80, 160)
    h_mZ1 = ROOT.TH1F("h_mZ1", "mZ1;m_{Z1} [GeV]", 70, 30, 110)
    h_mZ2 = ROOT.TH1F("h_mZ2", "mZ2:m_{Z2} [GeV]", 90, 0,  90)
    h_Z1_lepplus_pt = ROOT.TH1F("h_Z1_lepplus_pt", "Z1_lepplus_pt;l^{+} of Z1 p_{T} [GeV]", 70, 0,  140)
    h_Z2_lepplus_pt = ROOT.TH1F("h_Z2_lepplus_pt", "Z2_lepplus_pt;l^{+} of Z2 p_{T} [GeV]", 70, 0,  140)
    h_Z1_lepminus_pt = ROOT.TH1F("h_Z1_lepminus_pt", "Z1_lepminus_pt;l^{-} of Z1 p_{T} [GeV]", 70, 0,  140)
    h_Z2_lepminus_pt = ROOT.TH1F("h_Z2_lepminus_pt", "Z2_lepminus_pt;l^{-} of Z2 p_{T} [GeV]", 70, 0,  140)

    tree.Draw("m4l/1E3>>"+h_m4l.GetName(), "")
    tree.Draw("mZ1/1E3>>"+h_mZ1.GetName(), "")
    tree.Draw("mZ2/1E3>>"+h_mZ2.GetName(), "")
    tree.Draw("Z1_lepplus_pt/1E3>>"+h_Z1_lepplus_pt.GetName(), "")
    tree.Draw("Z2_lepplus_pt/1E3>>"+h_Z2_lepplus_pt.GetName(), "")
    tree.Draw("Z1_lepminus_pt/1E3>>"+h_Z1_lepminus_pt.GetName(), "")
    tree.Draw("Z2_lepminus_pt/1E3>>"+h_Z2_lepminus_pt.GetName(), "")

    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    hists = [h_m4l, h_mZ1, h_mZ2, h_Z1_lepplus_pt, h_Z2_lepplus_pt, h_Z1_lepminus_pt, h_Z2_lepminus_pt] 
    for hist in hists:
        hist.Draw()
        canvas.SaveAs(hist.GetName()+".pdf")

if __name__ == "__main__":
    make_plots("reduced_ntuple.root")
