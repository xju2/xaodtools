#!/usr/bin/env python

import ROOT
ROOT.gROOT.SetBatch()

import AtlasStyle
if not hasattr(ROOT, "myText"):
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/AtlasUtils.C")
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c")

def make_hist():
    f1 = ROOT.TFile.Open("all.root")
    f2 = ROOT.TFile.Open("can_ts.root")

    tree1 = f1.Get("bls")
    tree2 = f2.Get("physics")

    chi2_sig_xy = ROOT.TH1F("chi2_sig_xy", "chi2 signal;quad #chi^{2};Events/2", 50, 0, 100)
    chi2_sig_ts = chi2_sig_xy.Clone("chi2_sig_ts")

    chi2_left_xy = chi2_sig_xy.Clone("chi2_left_xy")
    chi2_left_ts = chi2_sig_xy.Clone("chi2_left_ts")

    chi2_right_xy = chi2_sig_xy.Clone("chi2_right_xy")
    chi2_right_ts = chi2_sig_xy.Clone("chi2_right_ts")

    sig_cut_xy = "m4l_fitted > 17.8 && m4l_fitted < 18.8"
    left_cut_xy = "m4l_fitted < 17.8"
    right_cut_xy = "m4l_fitted > 18.8 && m4l_fitted < 25"

    sig_cut_ts = "pX.m > 17.8 && pX.m < 18.8"
    left_cut_ts = "pX.m < 17.8"
    right_cut_ts = "pX.m > 18.8 && pX.m < 25"

    # times NDof
    tree1.Draw("x_chi2*5>>chi2_sig_xy", sig_cut_xy)
    tree1.Draw("x_chi2*5>>chi2_left_xy", left_cut_xy)
    tree1.Draw("x_chi2*5>>chi2_right_xy", right_cut_xy)

    tree2.Draw("pX.chi2>>chi2_sig_ts", sig_cut_ts)
    tree2.Draw("pX.chi2>>chi2_left_ts", left_cut_ts)
    tree2.Draw("pX.chi2>>chi2_right_ts", right_cut_ts)

    #zoom-in range
    m4l_xy = ROOT.TH1F("m4l_xy", "m4l;m_{4l} [GeV]; Events/0.2GeV", 60, 14, 26)
    m4l_ts = m4l_xy.Clone("m4l_ts")
    tree1.Draw("m4l_fitted>>m4l_xy")
    tree2.Draw("pX.m >> m4l_ts")

    # full range
    m4l_full_xy = ROOT.TH1F("m4l_full_xy", "m4l;m_{4l} [GeV]; Events/0.4GeV", 100, 10, 50)
    m4l_full_ts = m4l_full_xy.Clone("m4l_full_ts")
    tree1.Draw("m4l_fitted>>m4l_full_xy")
    tree2.Draw("pX.m >> m4l_full_ts")

    #with chi2-cut
    m4l_chi2_xy = ROOT.TH1F("m4l_chi2_xy", "m4l;m_{4l} [GeV]; Events/0.2GeV", 60, 14, 26)
    m4l_chi2_ts = m4l_chi2_xy.Clone("m4l_chi2_ts")
    chi2_cut = 10
    tree1.Draw("m4l_fitted>>m4l_chi2_xy", "x_chi2 < "+str(chi2_cut))
    tree2.Draw("pX.m >> m4l_chi2_ts", "pX.chi2 < "+str(chi2_cut))

    fout = ROOT.TFile.Open("chi2_hist.root", "recreate")
    chi2_sig_xy.Write()
    chi2_left_xy.Write()
    chi2_right_xy.Write()
    chi2_sig_ts.Write()
    chi2_left_ts.Write()
    chi2_right_ts.Write()

    m4l_xy.Write()
    m4l_ts.Write()
    m4l_full_xy.Write()
    m4l_full_ts.Write()
    m4l_chi2_xy.Write()
    m4l_chi2_ts.Write()
    fout.Close()

def save_compare(h_xy, h_ts, out_name):
    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    h_xy.SetLineColor(4)
    h_ts.SetLineColor(2)
    h_xy.SetMarkerColor(4)
    h_ts.SetMarkerColor(2)
    h_xy.SetMarkerStyle(20)
    h_ts.SetMarkerStyle(24)

    max_y = h_xy.GetMaximum()
    if h_ts.GetMaximum() > max_y:
        max_y = h_ts.GetMaximum()

    h_xy.GetYaxis().SetRangeUser(0, max_y*1.1)
    h_xy.Draw("EP")
    h_ts.Draw("same EP")
    legend = ROOT.myLegend(0.7, 0.8, 0.9, 0.9)
    legend.AddEntry(h_xy, "XY", "EP")
    legend.AddEntry(h_ts, "TS", "EP")
    legend.Draw("same")
    canvas.SaveAs(out_name+".pdf")

def norm_hist(hist):
    hist.Sumw2()
    hist.Scale(1./hist.Integral())
    return hist

def compare_hists():
    fin = ROOT.TFile.Open("chi2_hist.root")
    sig_xy = fin.Get("chi2_sig_xy")
    left_xy = fin.Get("chi2_left_xy")
    right_xy = fin.Get("chi2_right_xy")

    norm_hist(sig_xy)
    norm_hist(left_xy)
    norm_hist(right_xy)

    sig_ts = fin.Get("chi2_sig_ts")
    left_ts = fin.Get("chi2_left_ts")
    right_ts = fin.Get("chi2_right_ts")

    norm_hist(sig_ts)
    norm_hist(left_ts)
    norm_hist(right_ts)

    save_compare(sig_xy, sig_ts, "chi2_sig")
    save_compare(left_xy, left_ts, "chi2_left")
    save_compare(right_xy, right_ts, "chi2_right")

    save_compare(norm_hist(fin.Get("m4l_xy")), norm_hist(fin.Get("m4l_ts")), "m4l")
    save_compare(norm_hist(fin.Get("m4l_full_xy")), norm_hist(fin.Get("m4l_full_ts")), "m4l_full")
    save_compare(norm_hist(fin.Get("m4l_chi2_xy")), norm_hist(fin.Get("m4l_chi2_ts")), "m4l_chi2")

    fin.Close()

if __name__ == "__main__":
    make_hist()
    compare_hists()
