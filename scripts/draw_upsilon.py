#!/usr/bin/env python

import ROOT
import sys
ROOT.gROOT.SetBatch()
import string

import AtlasStyle
if not hasattr(ROOT, "myText"):
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/AtlasUtils.C")
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c")


chi2_cut = 10
def make_hists(file_names, out_name):
    tree = ROOT.TChain("physics", "physics")
    #tree.Add("data16.merged.p2689.root")
    #tree.Add("data16.merged.p2667.root")
    for file_ in file_names:
        tree.Add(file_)
    #tree.Add("test_data12.root")
    nentries = tree.GetEntries()
    print "total entries:",nentries

    h_upsilon = ROOT.TH1F("h_upsilon", "upsilon mass;m_{#varUpsilon} [GeV];Events / 50 MeV", 80, 8, 12)
    h_jpsi = ROOT.TH1F("h_jpsi", "jpsi mass;m_{#mu+#mu} [GeV];Events / 100 MeV", 100, 1, 11)

    nbins_4l = 100
    low_4l = 10
    hi_4l = 30
    h_m4l = ROOT.TH1F("h_m4l", "m4l;m_{#varUpsilon+#mu+#mu} [GeV];Events / 200 MeV", nbins_4l, low_4l, hi_4l)
    h_m4l_LHS = ROOT.TH1F("h_m4l_LHS", "m4l LHS;LHS m_{#varUpsilon+#mu+#mu} [GeV];Events / 200 MeV", nbins_4l, low_4l, hi_4l)
    h_m4l_RHS_U1 = ROOT.TH1F("h_m4l_RHS_U1", "m4l RHS;RHS U1 m_{#varUpsilon+#mu+#mu} [GeV];Events / 200 MeV", nbins_4l, low_4l, hi_4l)
    h_m4l_RHS_U2 = ROOT.TH1F("h_m4l_RHS_U2", "m4l RHS;RHS U2 m_{#varUpsilon+#mu+#mu} [GeV];Events / 200 MeV", nbins_4l, low_4l, hi_4l)

    ## disable branches
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("mUpsilon", 1)
    tree.SetBranchStatus("m4l", 1)
    tree.SetBranchStatus("m34", 1)
    tree.SetBranchStatus("n_muon", 1)
    tree.SetBranchStatus("vtx4l_chi2ndf", 1);
    #tree.SetBranchStatus("passTrigger", 1);
    #tree.SetBranchStatus("mu_p4", 1)


    for ientry in xrange(nentries):
        tree.GetEntry(ientry)
        has_4l = tree.n_muon > 3 and tree.m4l > 0 and tree.vtx4l_chi2ndf < chi2_cut and\
                tree.vtx4l_chi2ndf > 0
        if not has_4l:
            continue
        mU = tree.mUpsilon/1E3

        h_upsilon.Fill(mU)
        m4l = tree.m4l

        if has_4l:
            h_jpsi.Fill(tree.m34)
            if mU < 9:
                h_m4l_LHS.Fill(m4l)
            elif mU >= 9 and mU < 9.8:
                h_m4l.Fill(m4l)
            elif mU >= 9.8 and mU < 10.8:
                h_m4l_RHS_U1.Fill(m4l)
            else:
                h_m4l_RHS_U2.Fill(m4l)

    # save output
    #out_name = "hist_merged"+post_fix+"_"+str(chi2_cut)+".root"
    fout = ROOT.TFile.Open(out_name, "recreate")
    h_upsilon.Write()
    h_m4l.Write()
    h_m4l_LHS.Write()
    h_m4l_RHS_U1.Write()
    h_m4l_RHS_U2.Write()
    h_jpsi.Write()
    fout.Close()

    #return out_name

def draw(file_name, post_fix):
    f1 = ROOT.TFile.Open(file_name)

    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)

    x_pos = 0.70
    y_pos = 0.85
    n_rebin = 2


    h1 = f1.Get("h_upsilon")
    h1.Rebin(n_rebin)
    h1.GetYaxis().SetTitle("Events / {:.0f} MeV".format(50*n_rebin))
    h1.Draw("EP")
    ROOT.myText(x_pos, y_pos, 1, "Entries: {:.0f}".format(h1.Integral()))
    canvas.SaveAs("mUpsilon"+post_fix+"_"+str(chi2_cut)+".pdf")

    m4l_bin_width = 200

    h2 = f1.Get("h_m4l_LHS")
    h2.Rebin(n_rebin)
    h2.GetYaxis().SetTitle("Events / {:.0f} MeV".format(m4l_bin_width*n_rebin))
    h2.Draw("EP")
    ROOT.myText(x_pos, y_pos, 1, "Entries: {:.0f}".format(h2.Integral()))
    canvas.SaveAs("m4l_LHS"+post_fix+"_"+str(chi2_cut)+".pdf")

    h3 = f1.Get("h_m4l_RHS_U1")
    h3.Rebin(n_rebin)
    h3.GetYaxis().SetTitle("Events / {:.0f} MeV".format(m4l_bin_width*n_rebin))
    h3.Draw("EP")
    ROOT.myText(x_pos, y_pos, 1, "Entries: {:.0f}".format(h3.Integral()))
    canvas.SaveAs("m4l_RHS_U1"+post_fix+"_"+str(chi2_cut)+".pdf")

    h4 = f1.Get("h_m4l_RHS_U2")
    h4.Rebin(n_rebin)
    h4.GetYaxis().SetTitle("Events / {:.0f} MeV".format(m4l_bin_width*n_rebin))
    h4.Draw("EP")
    ROOT.myText(x_pos, y_pos, 1, "Entries: {:.0f}".format(h4.Integral()))
    canvas.SaveAs("m4l_RHS_U2"+post_fix+"_"+str(chi2_cut)+".pdf")

    h5 = f1.Get("h_m4l")
    h5.Rebin(n_rebin)
    h5.GetYaxis().SetTitle("Events / {:.0f} MeV".format(m4l_bin_width*n_rebin))
    h5.Draw("EP")
    ROOT.myText(x_pos, y_pos, 1, "Entries: {:.0f}".format(h5.Integral()))
    canvas.SaveAs("m4l"+post_fix+"_"+str(chi2_cut)+".pdf")

    h6 = f1.Get("h_jpsi")
    h6.Rebin(n_rebin)
    h6.GetYaxis().SetTitle("Events / {:.0f} MeV".format(100*n_rebin))
    h6.Draw("EP")
    ROOT.myText(x_pos, y_pos, 1, "Entries: {:.0f}".format(h6.Integral()))
    canvas.SaveAs("mJpsi"+post_fix+"_"+str(chi2_cut)+".pdf")

if __name__ == "__main__":

    base_name = "split_and_merge/merged_xa"
    input_tags = "abcdefghijklmnopqr"
    #input_files = [base_name+x for x in string.ascii_lowercase]
    input_files = [base_name+x+".root" for x in input_tags]
    out_name = "hist_a2r_100bins.root"

    make_hists(input_files, out_name)
    draw(out_name, "chi2Cut_100bins")
