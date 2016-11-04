#!/usr/bin/env python

import ROOT
from ROOT import RooRealVar
from ROOT import RooArgSet
from ROOT import RooArgList

import random
from array import array

import sys
import math
from optparse import OptionParser

import AtlasStyle
if not hasattr(ROOT, "my,Text"):
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/AtlasUtils.C")
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c")

ROOT.gROOT.SetBatch()
class Fit2L:
    def __init__(self, file_name, br_name="mass"):
        self.file_name = file_name
        self.br_name = br_name
        min_x = 8
        max_x = 12
        self.obs = RooRealVar("obs", "m4l", min_x, max_x)
        self.nbins = int((max_x - min_x) * 20)
        self.obs.setBin(self.nbins)
        #self.ws = ROOT.RooWorkspace("combined", "combined")
        self.dummy_hists = []
        self.chi2_cut = 10000

        # dionia selection or onia selection
        self.dionia_selection = False
        # different onia_pt cuts, 0:noCut, 1:<10, 2:10,20, 3:>20
        self.onia_pt_cut = 1
        # with 3mu4 trigger or not
        self.with_3mu4 = True

    def pass_onia_pt_cut(self, pT):
        pT_cut = int(self.onia_pt_cut)
        if pT_cut == 1:
            return (pT < 10)
        elif pT_cut == 2:
            return (pT >= 10 and pT < 20)
        elif pT_cut == 3:
            return pT > 20
        else:
            return True

    def print_onia_pt_cut(self):
        pT_cut = int(self.onia_pt_cut)
        if pT_cut == 1:
            return "p_{T}^{onia} in [0, 10) GeV"
        elif pT_cut == 2:
            return "p_{T}^{onia} in [10, 20) GeV"
        elif pT_cut == 3:
            return "p_{T}^{onia} in [20, #infty) GeV"
        else:
            return "p_{T}^{onia} in [-#infty, #infty] GeV"

    def print_dionia_selection(self):
        if self.dionia_selection:
            return "di-onia selection"
        else:
            return "onia selection"

    def print_trigger(self):
        if self.with_3mu4:
            return "with trigger"
        else:
            return "no trigger"

    def get_ws_name(self):
        return "ws_chi2Cut"+str(self.chi2_cut)+"_DiOnia"+str(self.dionia_selection)+"_pTCut"+str(self.onia_pt_cut)+"_withTrigger"+str(self.with_3mu4)+".root"

    def build_model(self):
        #mean = RooRealVar("mean", "mass of 1S", 9.46, 9.2, 9.7)
        #sigma = RooRealVar("sigma", "sigma of gaussian", 0.14, 0.09, 0.3)
        mean = RooRealVar("mean", "mass of 1S", 9.46)
        sigma = RooRealVar("sigma", "sigma of gaussian", 0.18)
        gaussian = ROOT.RooGaussian("gauss", "gauss", self.obs, mean, sigma)

        ## try Crystal Ball
        #alpha = RooRealVar("alpha", "alpha of CB", 5.9, 0, 100)
        #cb_n = RooRealVar("cb_n", "n of CB",  55.4, 0, 100)
        #gaussian = ROOT.RooCBShape("gauss", "gauss", self.obs, mean, sigma, alpha, cb_n)

        n_sig = RooRealVar("n_sig", "number of signal" , 5000, 0, 1E6)
        esig = ROOT.RooExtendPdf("esig", "esig", gaussian, n_sig)

        m2_shift = RooRealVar("m2_shift", "m2 shift", 0.56296)
        m2 = ROOT.RooFormulaVar("m2", "mass of 2S", "@0+ @1", RooArgList(mean, m2_shift))
        s2 = ROOT.RooFormulaVar("s2", "sigma*(unit+m2_shift)/mean", "@0*(1+@1/9.46)", RooArgList(sigma, m2_shift))
        g2 = ROOT.RooGaussian("g2", "gauss", self.obs, m2, s2)
        n2 = RooRealVar("n2", "number of 2S" , 800, 4, 100000)
        #n2 = ROOT.RooFormulaVar("n2", "number of 2S" ,"n_sig*0.26", RooArgList(n_sig))
        esig2 = ROOT.RooExtendPdf("esig2", "esig2", g2, n2)

        m3_shift = RooRealVar("m3_shift", "m3 shift", 0.8949)
        m3 = ROOT.RooFormulaVar("m3", "mass of 3S", "mean + m3_shift", RooArgList(mean, m3_shift))
        s3 = ROOT.RooFormulaVar("s3", "sigma*(unit+m2_shift)/mean", "@0*(1+@1/9.46)", RooArgList(sigma, m3_shift))
        g3 = ROOT.RooGaussian("g3", "gauss", self.obs, m3, s3)
        n3 = RooRealVar("n3", "number of 3S" , 1, 0, 1000)
        #n3 = ROOT.RooFormulaVar("n3", "number of 3S", "n2*0.45", RooArgList(n2))
        esig3 = ROOT.RooExtendPdf("esig3", "esig3", g3, n3)

        n_bkg = RooRealVar("n_bkg", "number of bkg" , 1000, 0, 1E7)
        p0 = RooRealVar("p0", "p0", -1E6, 1E6)
        p1 = RooRealVar("p1", "p1", -1E6, 1E6)
        p2 = RooRealVar("p2", "p2", -1E6, 1E6)
        p3 = RooRealVar("p3", "p3", -1E6, 1E6)
        p4 = RooRealVar("p4", "p4", -1E6, 1E6)
        #p0 = RooRealVar("p0", "p0", -4.56380e+03)
        #p1 = RooRealVar("p1", "p1", 1.48976e+03)
        #p2 = RooRealVar("p2", "p2", -7.60820e+01)
        #p3 = RooRealVar("p3", "p3", -1.35696e-02)
        #p4 = RooRealVar("p4", "p4", -1.46353e-02)
        bkg = ROOT.RooChebychev("bkg", "bkg", self.obs, RooArgList(p0, p1, p2, p3, p4))
        #bkg = ROOT.RooPolynomial("bkg", "bkg", self.obs, RooArgList(p0, p1, p2))
        ebkg = ROOT.RooExtendPdf("ebkg", "ebkg", bkg, n_bkg)
        model = ROOT.RooAddPdf("model", "model", RooArgList(esig, esig2, esig3, ebkg))
        getattr(self.ws, "import")(model)

    def get_data(self):
        fin = ROOT.TFile.Open(self.file_name)
        tree = fin.Get("upsilon")
        nentries = tree.GetEntries()
        print "total: ", nentries
        print "onia pT cut: ", self.onia_pt_cut
        obs_set = RooArgSet(self.obs)
        data = ROOT.RooDataSet("data", "data", obs_set)
        for ientry in xrange(nentries):
            tree.GetEntry(ientry)
            for i,m4l in enumerate(getattr(tree, self.br_name)):
                #m4l = m4l/1000
                if m4l > self.obs.getMax() or m4l < self.obs.getMin():
                    continue
                if tree.chi2[i] > self.chi2_cut:
                    continue
                # apply trigger requirement
                if self.with_3mu4 and not tree.trig_3mu4:
                    continue
                # apply dionia selection
                if self.dionia_selection and tree.pass_diOnia != 1:
                    continue
                # apply onia pT cut
                pT = tree.pt[i]
                if not self.pass_onia_pt_cut(pT):
                    continue

                self.obs.setVal(m4l)
                data.add(obs_set)

        getattr(self.ws, "import")(data)
        print "selected events: ", data.sumEntries()
        fin.Close()

    def fit(self):
        print "my configuration:",self.get_ws_name()
        if not hasattr(self, "ws"):
            self.ws = ROOT.RooWorkspace("combined", "combined")

        self.build_model()
        self.get_data()
        data = self.ws.obj("data")
        model = self.ws.obj("model")
        print "total data:", data.sumEntries()
        nll = model.createNLL(data)

        do_bkg_only = False
        n_sig = self.ws.var("n_sig")
        n2 = self.ws.var("n2")
        n3 = self.ws.var("n3")
        if do_bkg_only:
            n_sig.setVal(0.0)
            n2.setVal(0.0)
            n3.setVal(0.0)
            n_sig.setConstant()
            n2.setConstant()
            n3.setConstant()
            model.fitTo(data)
            nll_condition = nll.getVal()
            model.plotOn(frame, ROOT.RooFit.LineColor(4))
        else:
            nll_condition = 0.0

        #n_sig.setVal(50000)
        n_sig.setConstant(False)
        #n2.setVal(50000)
        #n3.setVal(50000)
        #n2.setConstant(False)
        #n3.setConstant(False)
        model.fitTo(data)
        nll_uncondition = nll.getVal()
        self.ws.saveSnapshot("splusb", self.ws.allVars())

        try:
            sig = math.sqrt(2*(nll_condition - nll_uncondition))
            print "significance: ", sig
        except:
            print nll_condition,nll_uncondition

        self.ws.writeToFile(self.get_ws_name())

    def plot(self, title):
        ## plot
        if not hasattr(self, "ws"):
            f1 = ROOT.TFile.Open(self.get_ws_name())
            self.ws = f1.Get("combined")

        self.ws.loadSnapshot("splusb")

        ## get number of signal and background in signal region
        obs = self.ws.var('obs')
        obs.setRange("signal", 9.2, 9.7)
        nsig = self.ws.obj("n_sig").getVal()
        nbkg = self.ws.obj("n_bkg").getVal()
        nall = self.ws.obj('model').expectedEvents(RooArgSet(obs))
        print "full range:", nsig, nbkg, nall
        frac_sig_full = self.ws.obj("gauss").createIntegral(
            RooArgSet(obs)
        ).getVal()
        frac_sig = self.ws.obj("gauss").createIntegral(
            RooArgSet(obs),
            ROOT.RooFit.Range("signal")
        ).getVal()
        frac_bkg_full = self.ws.obj('bkg').createIntegral(
            RooArgSet(obs)
        ).getVal()
        frac_bkg = self.ws.obj("bkg").createIntegral(
            RooArgSet(obs),
            ROOT.RooFit.Range("signal")
        ).getVal()
        frac_all = self.ws.obj("model").createIntegral(
            RooArgSet(obs),
            ROOT.RooFit.Range("signal")
        ).getVal()
        print "fraction: ", frac_sig, frac_sig_full, frac_bkg, frac_bkg_full, frac_all
        nsig *= frac_sig/frac_sig_full
        nbkg *= frac_bkg/frac_bkg_full
        print "signal region:", nsig, nbkg

        # prepare for the plot
        frame = obs.frame(self.nbins)
        self.ws.obj("data").plotOn(frame)
        self.ws.obj("model").plotOn(frame, ROOT.RooFit.LineColor(2))
        self.ws.obj("bkg").plotOn(
            frame, ROOT.RooFit.LineColor(4),
            ROOT.RooFit.Normalization(self.ws.obj("n_bkg").getVal(), ROOT.RooAbsReal.NumEvent)
        )

        self.ws.obj("gauss").plotOn(
            frame, ROOT.RooFit.LineColor(3),
            ROOT.RooFit.Normalization(self.ws.obj("n_sig").getVal(), ROOT.RooAbsReal.NumEvent)
        )
        self.ws.obj("g2").plotOn(
            frame, ROOT.RooFit.LineColor(6),
            ROOT.RooFit.Normalization(self.ws.obj("n2").getVal(), ROOT.RooAbsReal.NumEvent)
        )
        self.ws.obj("g3").plotOn(
            frame, ROOT.RooFit.LineColor(8),
            ROOT.RooFit.Normalization(self.ws.obj("n3").getVal(), ROOT.RooAbsReal.NumEvent)
        )
        canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
        canvas.SetTopMargin(0.07)
        canvas.SetBottomMargin(0.14)
        frame.Draw()
        max_y = frame.GetMaximum()
        frame.GetXaxis().SetTitle(title)
        #frame.SetMinimum(max_y*0.2)
        #frame.SetMaximum(max_y*1.7)
        ## add legend
        y_start = 0.43
        x_start = 0.60
        # title and mass and sigma
        ROOT.myText(0.2, 0.955, 1, "{},{},{}".format(self.print_onia_pt_cut(), self.print_dionia_selection(), self.print_trigger()))
        ROOT.myText(0.18, 0.87, 1, "m = {:.2f} GeV".format(self.ws.var("mean").getVal()))
        ROOT.myText(0.18, 0.87-0.05, 1, "#sigma = {:.2f} GeV".format(self.ws.var("sigma").getVal()))
        # left side
        ROOT.myText(0.2, y_start, 1, "Total: {:.0f}".format(self.ws.obj("data").sumEntries()))
        ROOT.myText(0.2, y_start-0.05, 1, "N(1S) = {:.1f}".format(self.ws.obj("n_sig").getVal()) )
        ROOT.myText(0.2, y_start-0.05*2, 1, "N(2S) = {:.1f}".format(self.ws.obj("n2").getVal()) )
        ROOT.myText(0.2, y_start-0.05*3, 1, "N(3S) = {:.1f}".format(self.ws.obj("n3").getVal()) )
        ROOT.myText(0.2, y_start-0.05*4, 1, "N(bkg) = {:.1f}".format(self.ws.obj("n_bkg").getVal()) )
        # right side
        ROOT.myText(x_start, y_start, 1, "In [9.2, 9.7] GeV")
        ROOT.myText(x_start, y_start-0.05, 1, "with #chi^{2} < "+str(self.chi2_cut))
        ROOT.myText(x_start, y_start-0.05*2, 1, "N(1S) = {:.1f}".format(nsig))
        ROOT.myText(x_start, y_start-0.05*3, 1, "N(bkg) = {:.1f}".format(nbkg))
        ROOT.myText(x_start, y_start-0.05*4, 1, "S/sqrt(B):{:.1f}".format(nsig/math.sqrt(nbkg)))

        canvas.SaveAs(self.get_ws_name().replace("root", "pdf").replace("ws_", "fit_"))
        self.ws.obj("s2").Print()
        self.ws.obj("s3").Print()

    def add_dummy_entry(self, legend, color, label):
        h1 = ROOT.TH1F(label, label, 10, 0, 10)
        h1.SetLineColor(color)
        legend.AddEntry(h1, label, "L")
        self.dummy_hists.append(h1)

if __name__ == "__main__":
    usage = sys.argv[0]+" file_name br_name"
    parser = OptionParser(usage=usage)
    parser.add_option('--dionia', action="store_true", dest="dionia", help="apply dionia selection", default=False)
    parser.add_option('--oniaPt', dest='oniapt', help="onia pT cut, 0/1/2/3", default=0)
    parser.add_option('--noTrigger', action="store_true", dest='notrigger', help="don't apply trigger", default=False)
    parser.add_option('--title', dest='title', help="title of x-axis", default="m_{#mu#mu} [GeV]")


    (options, args) = parser.parse_args()
    if len(args) < 2:
        parser.print_help()
        exit(1)

    file_name = args[0]
    br_name = "mass"
    if len(sys.argv) > 2:
        br_name = args[1]

    #cuts = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 8, 10, 10000]
    cuts = [3.]
    for cut in cuts:
        fit_4l = Fit2L(file_name, br_name)
        if options.dionia:
            fit_4l.dionia_selection = True

        if options.notrigger:
            fit_4l.with_3mu4 = False

        fit_4l.onia_pt_cut = options.oniapt
        fit_4l.chi2_cut = cut
        fit_4l.fit()
        fit_4l.plot(options.title)
