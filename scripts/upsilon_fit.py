#!/usr/bin/env python

import ROOT
from ROOT import RooRealVar

import sys
import math

ROOT.gROOT.SetBatch()

class Fit4L:
    def __init__(self, file_name, br_name):
        self.fin = ROOT.TFile.Open(file_name)
        self.br_name = br_name
        self.tree = self.fin.Get("bls")
        self.obs = RooRealVar("obs", "m4l", 14, 21)
        self.nbins = 35
        self.obs.setBin(self.nbins)
        self.ws = ROOT.RooWorkspace("combined", "combined")

    def build_model(self):
        #mean = RooRealVar("mean", "mean of gaussian", 18.1, 17, 19)
        #sigma = RooRealVar("sigma", "sigma of gaussian", 0.21, 0.15, 1.0)
        mean = RooRealVar("mean", "mean of gaussian", 18.15)
        sigma = RooRealVar("sigma", "sigma of gaussian", 0.19)
        gaussian = ROOT.RooGaussian("gauss", "gauss", self.obs, mean, sigma)
        n_sig = RooRealVar("n_sig", "number of signal" , 50, 0, 1000)
        esig = ROOT.RooExtendPdf("esig", "esig", gaussian, n_sig)

        n_bkg = RooRealVar("n_bkg", "number of bkg" , 100, 0, 5000)
        p0 = RooRealVar("p0", "p0", -1E6, 1E6)
        p1 = RooRealVar("p1", "p1", -1E6, 1E6)
        p2 = RooRealVar("p2", "p2", -1E6, 1E6)
        p3 = RooRealVar("p3", "p3", -1E6, 1E6)
        p4 = RooRealVar("p4", "p4", -1E6, 1E6)
        #bkg = ROOT.RooChebychev("bkg", "bkg", self.obs, ROOT.RooArgList(p0, p1, p2, p3, p4))
        bkg = ROOT.RooChebychev("bkg", "bkg", self.obs, ROOT.RooArgList(p0, p1))
        ebkg = ROOT.RooExtendPdf("ebkg", "ebkg", bkg, n_bkg)
        model = ROOT.RooAddPdf("model", "model", ROOT.RooArgList(esig, ebkg))
        getattr(self.ws, "import")(model)

    def get_data(self):
        nentries = self.tree.GetEntries()
        print "total: ", nentries
        obs_set = ROOT.RooArgSet(self.obs)
        data = ROOT.RooDataSet("data", "data", obs_set)
        for ientry in xrange(nentries):
            self.tree.GetEntry(ientry)
            m4l = getattr(self.tree, self.br_name)
            if m4l > self.obs.getMax() or m4l < self.obs.getMin():
                continue
            # add chi2 cut
            #if self.tree.x_chi2 > 5:
            #    continue

            self.obs.setVal(m4l)
            data.add(obs_set)
        getattr(self.ws, "import")(data)

    def fit(self):
        self.build_model()
        self.get_data()
        data = self.ws.obj("data")
        model = self.ws.obj("model")
        nll = model.createNLL(data)

        n_sig = self.ws.var("n_sig")
        n_sig.setVal(0.0)
        n_sig.setConstant()
        #self.minimize(nll)
        model.fitTo(data)
        nll_condition = nll.getVal()

        n_sig.setVal(50)
        n_sig.setConstant(False)
        model.fitTo(data)
        #self.minimize(nll)
        nll_uncondition = nll.getVal()


        print "significance", math.sqrt(2*(nll_condition - nll_uncondition))

        ## plot
        frame = self.obs.frame(self.nbins)
        data.plotOn(frame)
        model.plotOn(frame, ROOT.RooFit.LineColor(2))
        canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
        frame.Draw()
        frame.GetYaxis().SetRangeUser(0, 100)
        canvas.SaveAs("fit.pdf")

        self.save()

    def save(self):
        self.ws.writeToFile("combined.root")
        self.fin.Close()

    def minimize(self, nll):
        minim = ROOT.RooMinimizer(nll)
        minim.setStrategy(1)
        minim.setProfile()
        minim.setPrintLevel(1);
        status = minim.minimize("Minuit2", ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo());
        print "fit stats:", status


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print sys.argv[0], " file_name [br_name]"
        exit(1)

    file_name = sys.argv[1]
    br_name = "m4l_track"
    if len(sys.argv) > 2:
        br_name = sys.argv[2]

    fit_4l = Fit4L(file_name, br_name)
    fit_4l.fit()
