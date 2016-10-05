#!/usr/bin/env python

import ROOT 
ROOT.gROOT.SetBatch()

import AtlasStyle
if not hasattr(ROOT, "my,Text"):
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/AtlasUtils.C")
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c")

import compare_mine_TS as cpt
def make():
    f1 = ROOT.TFile.Open("input.root")
    tree = f1.Get("bls")
    cut = ROOT.TCut("m34 > 9.2 && m34 < 9.7")
    #cut = ROOT.TCut("1==1")

    h_chi2 = ROOT.TH1F("h_chi2", "chi2;Log10(#chi^{2});Events/0.2", 40, -2, 6)
    tree.Draw("TMath::Log10(x_chi2*5) >> h_chi2", cut)

    h_chi2_mu1_pt = ROOT.TH2F("h_chi2_mu1_pt", "chi2 vs mu1;Log10(#chi^{2});pT muon-1[GeV]", 40, -2, 6, 40, 3E3, 13E3)
    tree.Draw("m1_track_pt:TMath::Log10(x_chi2*5) >> h_chi2_mu1_pt", cut)
    h_chi2_mu2_pt = ROOT.TH2F("h_chi2_mu2_pt", "chi2 vs mu2;Log10(#chi^{2});pT muon-2[GeV]", 40, -2, 6, 40, 3E3, 13E3)
    tree.Draw("m2_track_pt:TMath::Log10(x_chi2*5) >> h_chi2_mu2_pt", cut)
    h_chi2_mu3_pt = ROOT.TH2F("h_chi2_mu3_pt", "chi2 vs mu3;Log10(#chi^{2});pT muon-3[GeV]", 40, -2, 6, 40, 3E3, 13E3)
    tree.Draw("m3_track_pt:TMath::Log10(x_chi2*5) >> h_chi2_mu3_pt", cut)
    h_chi2_mu4_pt = ROOT.TH2F("h_chi2_mu4_pt", "chi2 vs mu4;Log10(#chi^{2});pT muon-4[GeV]", 40, -2, 6, 40, 3E3, 13E3)
    tree.Draw("m4_track_pt:TMath::Log10(x_chi2*5) >> h_chi2_mu4_pt", cut)

    h_chi2_mu1_eta = ROOT.TH2F("h_chi2_mu1_eta", "chi2 vs mu1;Log10(#chi^{2});#eta muon-1[GeV]", 40, -2, 6, 54, -2.7, 2.7)
    tree.Draw("m1_track_eta:TMath::Log10(x_chi2*5) >> h_chi2_mu1_eta", cut)
    h_chi2_mu2_eta = ROOT.TH2F("h_chi2_mu2_eta", "chi2 vs mu2;Log10(#chi^{2});#eta muon-2[GeV]", 40, -2, 6, 54, -2.7, 2.7)
    tree.Draw("m2_track_eta:TMath::Log10(x_chi2*5) >> h_chi2_mu2_eta", cut)
    h_chi2_mu3_eta = ROOT.TH2F("h_chi2_mu3_eta", "chi2 vs mu3;Log10(#chi^{2});#eta muon-3[GeV]", 40, -2, 6, 54, -2.7, 2.7)
    tree.Draw("m3_track_eta:TMath::Log10(x_chi2*5) >> h_chi2_mu3_eta", cut)
    h_chi2_mu4_eta = ROOT.TH2F("h_chi2_mu4_eta", "chi2 vs mu4;Log10(#chi^{2});#eta muon-4[GeV]", 40, -2, 6, 54, -2.7, 2.7)
    tree.Draw("m4_track_eta:TMath::Log10(x_chi2*5) >> h_chi2_mu4_eta", cut)

    h_chi2_mu1_d0 = ROOT.TH2F("h_chi2_mu1_d0", "chi2 vs mu1;Log10(#chi^{2});d0 muon-1[GeV]", 40, -2, 6, 50, -1, 1.)
    tree.Draw("m1_trackD0PV:TMath::Log10(x_chi2*5) >> h_chi2_mu1_d0", cut)
    h_chi2_mu2_d0 = ROOT.TH2F("h_chi2_mu2_d0", "chi2 vs mu2;Log10(#chi^{2});d0 muon-2[GeV]", 40, -2, 6, 50, -1., 1.)
    tree.Draw("m2_trackD0PV:TMath::Log10(x_chi2*5) >> h_chi2_mu2_d0", cut)
    h_chi2_mu3_d0 = ROOT.TH2F("h_chi2_mu3_d0", "chi2 vs mu3;Log10(#chi^{2});d0 muon-3[GeV]", 40, -2, 6, 50, -1., 1.)
    tree.Draw("m3_trackD0PV:TMath::Log10(x_chi2*5) >> h_chi2_mu3_d0", cut)
    h_chi2_mu4_d0 = ROOT.TH2F("h_chi2_mu4_d0", "chi2 vs mu4;Log10(#chi^{2});d0 muon-4[GeV]", 40, -2, 6, 50, -1., 1.)
    tree.Draw("m4_trackD0PV:TMath::Log10(x_chi2*5) >> h_chi2_mu4_d0", cut)

    h_chi2_mu1_d0sig = ROOT.TH2F("h_chi2_mu1_d0sig", "chi2 vs mu1;Log10(#chi^{2});d0sig muon-1[GeV]", 40, -2, 6, 50, -1, 1.)
    tree.Draw("m1_trackD0PV:TMath::Log10(x_chi2*5) >> h_chi2_mu1_d0sig", cut)
    h_chi2_mu2_d0sig = ROOT.TH2F("h_chi2_mu2_d0sig", "chi2 vs mu2;Log10(#chi^{2});d0sig muon-2[GeV]", 40, -2, 6, 50, -1., 1.)
    tree.Draw("m2_trackD0PV:TMath::Log10(x_chi2*5) >> h_chi2_mu2_d0sig", cut)
    h_chi2_mu3_d0sig = ROOT.TH2F("h_chi2_mu3_d0sig", "chi2 vs mu3;Log10(#chi^{2});d0sig muon-3[GeV]", 40, -2, 6, 50, -1., 1.)
    tree.Draw("m3_trackD0PV:TMath::Log10(x_chi2*5) >> h_chi2_mu3_d0sig", cut)
    h_chi2_mu4_d0sig = ROOT.TH2F("h_chi2_mu4_d0sig", "chi2 vs mu4;Log10(#chi^{2});d0sig muon-4[GeV]", 40, -2, 6, 50, -1., 1.)
    tree.Draw("m4_trackD0PV:TMath::Log10(x_chi2*5) >> h_chi2_mu4_d0sig", cut)

    h_chi2_mu1_z0 = ROOT.TH2F("h_chi2_mu1_z0", "chi2 vs mu1;Log10(#chi^{2});z0 muon-1[GeV]", 40, -2, 6, 50, -2., 2.)
    tree.Draw("m1_trackZ0PV:TMath::Log10(x_chi2*5) >> h_chi2_mu1_z0", cut)
    h_chi2_mu2_z0 = ROOT.TH2F("h_chi2_mu2_z0", "chi2 vs mu2;Log10(#chi^{2});z0 muon-2[GeV]", 40, -2, 6, 50, -2., 2.)
    tree.Draw("m2_trackZ0PV:TMath::Log10(x_chi2*5) >> h_chi2_mu2_z0", cut)
    h_chi2_mu3_z0 = ROOT.TH2F("h_chi2_mu3_z0", "chi2 vs mu3;Log10(#chi^{2});z0 muon-3[GeV]", 40, -2, 6, 50, -2., 2.)
    tree.Draw("m3_trackZ0PV:TMath::Log10(x_chi2*5) >> h_chi2_mu3_z0", cut)
    h_chi2_mu4_z0 = ROOT.TH2F("h_chi2_mu4_z0", "chi2 vs mu4;Log10(#chi^{2});z0 muon-4[GeV]", 40, -2, 6, 50, -2., 2.)
    tree.Draw("m4_trackZ0PV:TMath::Log10(x_chi2*5) >> h_chi2_mu4_z0", cut)

    ## p vs pT
    hists = []
    for i in range(1,5):
        hist = ROOT.TH2F("h_mu"+str(i)+"_p_pT", ";track p"+str(i)+"[GeV];track pT"+str(i)+" [GeV]", 100, 0, 50, 100, 0., 50.)
        tree.Draw("m"+str(i)+"_track_pt/1E3:m"+str(i)+"_track_p/1E3 >> h_mu"+str(i)+"_p_pT", cut)
        hists.append(hist)

    fout = ROOT.TFile.Open("output_chi2.root", "recreate")
    h_chi2.Write()
    h_chi2_mu1_pt.Write()
    h_chi2_mu2_pt.Write()
    h_chi2_mu3_pt.Write()
    h_chi2_mu4_pt.Write()
    h_chi2_mu1_eta.Write()
    h_chi2_mu2_eta.Write()
    h_chi2_mu3_eta.Write()
    h_chi2_mu4_eta.Write()

    h_chi2_mu1_d0.Write()
    h_chi2_mu2_d0.Write()
    h_chi2_mu3_d0.Write()
    h_chi2_mu4_d0.Write()
    h_chi2_mu1_z0.Write()
    h_chi2_mu2_z0.Write()
    h_chi2_mu3_z0.Write()
    h_chi2_mu4_z0.Write()

    for hist in hists:
        hist.Write()

    fout.Close()
    f1.Close()

def make_ts():
    f1 = ROOT.TFile.Open("can_ts.root")
    tree = f1.Get("physics")
    h_chi2_ts = ROOT.TH1F("h_chi2_ts", "chi2;Log10(#chi^{2});Events/0.2", 40, -2, 6)
    tree.Draw("TMath::Log10(pX.chi2) >> h_chi2_ts", "")

    fout = ROOT.TFile.Open("output_chi2.root", "update")
    h_chi2_ts.Write()
    fout.Close()


class MakeHists:
    def __init__(self, file_name, tree_name):
        self.f1 = ROOT.TFile.Open(file_name)
        self.tree = self.f1.Get(tree_name)

        self.cut_left = ROOT.TCut("mass > 8.5 && mass < 9.0")
        self.cut_mid = ROOT.TCut("mass > 9.2 && mass < 9.7")
        self.cut_right = ROOT.TCut("mass > 11 && mass < 11.5")
        self.bands = {
            "LHS":self.cut_left,
            "SR":self.cut_mid,
            "RHS":self.cut_right
        }
        print "total entries:",self.tree.GetEntries()

    def make_hists(self, h_temp, var):
        hists = []
        for name,cut in self.bands.iteritems():
            hist = h_temp.Clone(h_temp.GetName()+"_"+name)
            self.tree.Draw(var+name, cut)
            hists.append(hist)

        return hists

    def make_upsilon(self, out_name):
        # plot chi2 in upsilon side-band

        hists = []
        h_chi2_mass = ROOT.TH2F("h_chi2_mass", "chi2 vs mass;Log10(#chi^{2});m_{onia} [GeV]", 40, -2, 6, 40, 8, 12)
        self.tree.Draw("mass:TMath::Log10(chi2)>>h_chi2_mass")
        hists.append(h_chi2_mass)

        h_chi2 = ROOT.TH1F("h_chi2", "chi2;Log10(#chi^{2});Events", 40, -2, 6)
        hists += self.make_hists(h_chi2, "TMath::Log10(chi2)>>h_chi2_")

        h_chi2_linear = ROOT.TH1F("h_chi2_linear", "chi2;#chi^{2};Events", 10, 0, 10)
        hists += self.make_hists(h_chi2_linear, "chi2>>h_chi2_linear_")
        # save output
        self.out_name = out_name
        fout = ROOT.TFile.Open(out_name, "recreate")
        for hist in hists:
            hist.Write()
        fout.Close()

    def compare_hists(self):
        if not hasattr(self, "out_name"):
            self.out_name = "hist_sideband.root"

        fin = ROOT.TFile.Open(self.out_name)
        base_name = "h_chi2_linear_"
        hists = {}
        for name,cut in self.bands.iteritems():
            hists[name] = fin.Get(base_name+name)
            print name,":",hists[name].Integral()

        #save_compare([cpt.norm_hist(x) for x in hists.values()], hists.keys(), "chi2_upsilon")
        save_compare(hists.values(), hists.keys(), "chi2_upsilon_linear", False)
        save_compare(hists.values(), hists.keys(), "chi2_upsilon_linear", True)
        fin.Close()

def get_chi2_slice(hist):
    # [<10, 10-100, 100]
    xaxis = hist.GetXaxis()
    bin_10 = xaxis.FindBin(1)
    bin_100 = xaxis.FindBin(2)
    h1 = hist.ProjectionY(hist.GetName()+"_s1", 1, bin_10)
    h2 = hist.ProjectionY(hist.GetName()+"_s2", bin_10+1, bin_100)
    h3 = hist.ProjectionY(hist.GetName()+"_s3", bin_100+1, xaxis.GetNbins())
    return [h1, h2, h3]

def save_compare(hists, labels, out_name, is_log=False):
    if len(hists) != len(labels):
        return None
    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    if is_log:
        canvas.SetLogy()
    colors = [4, 2, 8]
    styles = [20, 24, 34]

    max_y = -1
    for i,hist in enumerate(hists):
        hist.SetLineColor( colors[i] )
        hist.SetMarkerColor( colors[i] )
        hist.SetMarkerStyle( styles[i] )
        if hist.GetMaximum() > max_y:
            max_y = hist.GetMaximum()

    h1 = hists[0]
    if is_log:
        h1.GetYaxis().SetRangeUser(1E-3, max_y*100)
    else:
        h1.GetYaxis().SetRangeUser(0, max_y*1.1)

    h1.Draw("EP")
    for i in range(1, len(hists)):
        hists[i].Draw("same EP")

    legend = ROOT.myLegend(0.7, 0.9-0.05*len(hists), 0.9, 0.9)
    for hist,label in zip(hists, labels):
        legend.AddEntry(hist, label, "EP")

    legend.Draw("same")
    if is_log:
        canvas.SaveAs(out_name+"_log.pdf")
        canvas.SaveAs(out_name+"_log.eps")
    else:
        canvas.SaveAs(out_name+".pdf")
        canvas.SaveAs(out_name+".eps")

def plot():
    fin = ROOT.TFile.Open("output_chi2.root")
    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    for key in fin.GetListOfKeys():
        h1 = key.ReadObj()
        name = h1.GetName()
        if "pt" in name:
            h1.GetYaxis().SetLabelSize(0.035)
        if h1.ClassName() == "TH2F":
            h1.Draw("box")
        else:
            h1.Draw("EP")
        canvas.SaveAs(name+".pdf")

        if h1.ClassName() == "TH2F":
            py = h1.ProfileY()
            py.Draw()
            canvas.SaveAs(name+"_profileY.pdf")

    cpt.save_compare(cpt.norm_hist(fin.Get("h_chi2")), cpt.norm_hist(fin.Get("h_chi2_ts")), "chi2_cmp")
    save_compare([cpt.norm_hist(fin.Get("h_chi2")), cpt.norm_hist(fin.Get("h_chi2_ts"))],["#Upsilon(1S)#Upsilon(1S)","baseline"],"chi2_cmp")

    # compare the slices
    labels = ["#chi^{2} < 10", "#chi^{2} [10, 100]", "#chi^{2} > 100"]
    for i in range(4):
        hists = get_chi2_slice( fin.Get("h_chi2_mu"+str(i+1)+"_d0") )
        save_compare([cpt.norm_hist(x) for x in hists], labels, "cmp_mu"+str(i+1)+"_d0")

    for i in range(4):
        hists = get_chi2_slice( fin.Get("h_chi2_mu"+str(i+1)+"_z0") )
        save_compare([cpt.norm_hist(x) for x in hists], labels, "cmp_mu"+str(i+1)+"_z0")

    fin.Close()

def test():
    #hist_maker = MakeHists("all_v5.root", "upsilon")
    #hist_maker = MakeHists("all_v5_withGRL.root", "upsilon") ##
    hist_maker = MakeHists("all_v7_chi2Studies.root", "upsilon") ##
    hist_maker.make_upsilon("hist_sideband.root")
    hist_maker.compare_hists()

if __name__ == "__main__":
    #make()
    #make_ts()
    #plot()
    test()
