#!/usr/bin/env python
"""
read ntuple produced by Truth_JETMET...
make some validation plots
"""

import ROOT
ROOT.gROOT.SetBatch()
from optparse import OptionParser

def make_plots(file_name, post_fix):
    import AtlasStyle
    f1 = ROOT.TFile.Open(file_name)
    tree = f1.Get("physics")

    h_m4l = ROOT.TH1F("h_m4l", "m4l;m_{4l} [GeV];Events/2 GeV", 100, 100, 600)
    h_mZ1 = ROOT.TH1F("h_mZ1", "mZ1;m_{Z1} [GeV];Events/1 GeV", 70, 30, 110)
    h_mZ2 = ROOT.TH1F("h_mZ2", "mZ2;m_{Z2} [GeV];Events/2 GeV", 60, 0,  120)
    h_Z1_lepplus_pt = ROOT.TH1F("h_Z1_lepplus_pt", "Z1_lepplus_pt;l^{+} of Z1 p_{T} [GeV];Events/4 GeV", 35, 0,  140)
    h_Z2_lepplus_pt = ROOT.TH1F("h_Z2_lepplus_pt", "Z2_lepplus_pt;l^{+} of Z2 p_{T} [GeV];Events/4 GeV", 35, 0,  140)
    h_Z1_lepminus_pt = ROOT.TH1F("h_Z1_lepminus_pt", "Z1_lepminus_pt;l^{-} of Z1 p_{T} [GeV];Events/ 4 GeV", 35, 0,  140)
    h_Z2_lepminus_pt = ROOT.TH1F("h_Z2_lepminus_pt", "Z2_lepminus_pt;l^{-} of Z2 p_{T} [GeV];Events/ 4 GeV", 35, 0,  140)

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
        canvas.SaveAs(post_fix+"_"+hist.GetName()+".pdf")

if __name__ == "__main__":
    usage = "%prog file_name out_tag"
    parser = OptionParser(usage=usage, description="read truth file, plot basic variables")
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print parser.print_help()
        exit(1)

    file_ = args[0]
    out_ = args[1]
    make_plots(file_, out_)
