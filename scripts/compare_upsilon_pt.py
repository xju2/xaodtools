#!/usr/bin/env python
import sys
import ROOT
import AtlasStyle

ROOT.gROOT.SetBatch()

if len(sys.argv) < 3:
    print sys.argv[0]," file_8TeV, file_13TeV"
    exit(1)

tree_name = "upsilon"
f8tev_name = sys.argv[1]
f13tev_name = sys.argv[2]

f8tev = ROOT.TFile.Open(f8tev_name)
f13tev = ROOT.TFile.Open(f13tev_name)

h_upsilon_pt = ROOT.TH1F("h_upsilon_pt", "upsilon pT;onia p_{T} [GeV]", 30, 0, 30)
h_upsilon_pt_8TeV = h_upsilon_pt.Clone("h_upsilon_pt_8TeV")
h_upsilon_pt_13TeV = h_upsilon_pt.Clone("h_upsilon_pt_13TeV")
h_upsilon_pt_13TeV_3mu4 = h_upsilon_pt.Clone("h_upsilon_pt_13TeV_3mu4")

tree8 = f8tev.Get(tree_name)
tree13 = f13tev.Get(tree_name)

cut = ROOT.TCut("mass > 8 && mass < 12 && pt <= 30 && chi2 < 3 && pass_diOnia == 1")
tree8.Draw("pt>>"+h_upsilon_pt_8TeV.GetName(), cut)
tree13.Draw("pt>>"+h_upsilon_pt_13TeV.GetName(), cut)
tree13.Draw("pt>>"+h_upsilon_pt_13TeV_3mu4.GetName(), cut+ROOT.TCut("trig_3mu4"))
print h_upsilon_pt_8TeV.Integral(), h_upsilon_pt_13TeV.Integral(), h_upsilon_pt_13TeV_3mu4.Integral()

canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
h_upsilon_pt_8TeV.Scale(1.0/h_upsilon_pt_8TeV.Integral())
h_upsilon_pt_13TeV.Scale(1.0/h_upsilon_pt_13TeV.Integral())
h_upsilon_pt_13TeV_3mu4.Scale(1.0/h_upsilon_pt_13TeV_3mu4.Integral())

h_upsilon_pt_8TeV.SetLineColor(2)
h_upsilon_pt_13TeV.SetLineColor(4)
h_upsilon_pt_13TeV_3mu4.SetLineColor(3)

h_upsilon_pt_13TeV.Draw()
h_upsilon_pt_8TeV.Draw("same")
h_upsilon_pt_13TeV_3mu4.Draw("same")

legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.SetFillColor(0)
legend.SetBorderSize(0)
legend.SetTextFont(42)
legend.SetTextSize(0.04)
legend.AddEntry(h_upsilon_pt_8TeV, "8TeV", "L")
legend.AddEntry(h_upsilon_pt_13TeV, "13TeV", "L")
legend.AddEntry(h_upsilon_pt_13TeV_3mu4, "13TeV 3mu4", "L")
legend.Draw()

canvas.SaveAs("cmp_upsilon_pT.pdf")
canvas.SaveAs("cmp_upsilon_pT.eps")
