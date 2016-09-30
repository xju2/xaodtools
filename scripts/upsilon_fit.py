#!/usr/bin/env python

import ROOT
from ROOT import RooRealVar
from ROOT import RooArgSet
from ROOT import RooArgList

import random
from array import array

import sys
import math

ROOT.gROOT.SetBatch()

class Fit4L:
    def __init__(self, file_name, br_name):
        self.file_name = file_name
        self.br_name = br_name
        min_x = 15
        max_x = 25
        self.obs = RooRealVar("obs", "m4l", min_x, max_x)
        self.nbins = (max_x - min_x) * 5
        self.obs.setBin(self.nbins)
        self.ws = ROOT.RooWorkspace("combined", "combined")

    def build_model(self):
        #mean = RooRealVar("mean", "mean of gaussian", 18.1, 17, 19)
        #sigma = RooRealVar("sigma", "sigma of gaussian", 0.21, 0.01, 1.0)
        mean = RooRealVar("mean", "mean of gaussian", 18.06)
        sigma = RooRealVar("sigma", "sigma of gaussian", 0.19)
        #sigma = RooRealVar("sigma", "sigma of gaussian", 0.08)
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
        fin = ROOT.TFile.Open(self.file_name)
        tree = fin.Get("bls")
        nentries = tree.GetEntries()
        print "total: ", nentries
        obs_set = RooArgSet(self.obs)
        data = ROOT.RooDataSet("data", "data", obs_set)
        for ientry in xrange(nentries):
            tree.GetEntry(ientry)
            m4l = getattr(tree, self.br_name)
            if m4l > self.obs.getMax() or m4l < self.obs.getMin():
                continue
            # add chi2 cut
            #if tree.x_chi2 > 10:
            #    continue
            #if tree.m34 < 4:
            #    continue

            self.obs.setVal(m4l)
            data.add(obs_set)
        getattr(self.ws, "import")(data)
        fin.Close()

    def get_data_ts(self):
        """
        get data provided by TieSheng
        """
        fin = ROOT.TFile.Open(self.file_name)
        tree = fin.Get("physics")
        nentries = tree.GetEntries()
        print "total: ", nentries
        obs_set = RooArgSet(self.obs)
        data = ROOT.RooDataSet("data", "data", obs_set)
        for ientry in xrange(nentries):
            tree.GetEntry(ientry)
            m4l = getattr(tree, "pX.m")
            chi2 = getattr(tree, "pX.chi2")
            if m4l > self.obs.getMax() or m4l < self.obs.getMin():
                continue
            # add chi2 cut
            if chi2 <= 50:
                continue

            self.obs.setVal(m4l)
            data.add(obs_set)
        getattr(self.ws, "import")(data)
        fin.Close()

    def get_random_data(self):
        NSIGNAL = 86
        fin = ROOT.TFile.Open(self.file_name)
        tree = fin.Get("bls")
        nentries = tree.GetEntries()
        print "total: ", nentries
        obs_set = RooArgSet(self.obs)
        data = ROOT.RooDataSet("data", "data", obs_set)
        sample = random.sample([x for x in range(nentries)], NSIGNAL)
        #for ientry in xrange(nentries):
        for ientry in sample:
            tree.GetEntry(ientry)
            m4l = getattr(tree, self.br_name)
            if m4l > self.obs.getMax() or m4l < self.obs.getMin():
                continue
            # add chi2 cut
            #if tree.x_chi2 > 10:
            #    continue
            #if tree.m34 < 4:
            #    continue

            self.obs.setVal(m4l)
            data.add(obs_set)

        fin.Close()
        return data

    def fit(self):
        self.build_model()
        self.get_data()
        #self.get_data_ts()
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

        try:
            sig = math.sqrt(2*(nll_condition - nll_uncondition))
            print "significance: ", sig
        except:
            print nll_condition,nll_uncondition

        ## plot
        frame = self.obs.frame(self.nbins)
        data.plotOn(frame)
        model.plotOn(frame, ROOT.RooFit.LineColor(2))
        canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
        frame.Draw()
        #frame.GetYaxis().SetRangeUser(0, 120)
        canvas.SaveAs("fit.pdf")

        self.save()

    def save(self):
        self.ws.writeToFile("combined.root")

    def minimize(self, nll):
        minim = ROOT.RooMinimizer(nll)
        minim.setStrategy(1)
        minim.setProfile()
        minim.setPrintLevel(1);
        status = minim.minimize("Minuit2", ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo());
        print "fit stats:", status

    def fit_signal_random(self):
        mean = RooRealVar("mean", "mean of gaussian", 18.1, 17, 19)
        sigma = RooRealVar("sigma", "sigma of gaussian", 0.21, 0.01, 1.0)
        gaussian = ROOT.RooGaussian("gauss", "gauss", self.obs, mean, sigma)
        NTRIES = 10000
        fout = ROOT.TFile.Open("ensamble.root", "recreate")
        br_mean = array('f', [0])
        br_sigma = array('f', [0])
        out_tree = ROOT.TTree("physics", "physics")
        out_tree.Branch("mean", br_mean, "mean/F")
        out_tree.Branch("sigma", br_sigma, "sigma/F")
        hists = []
        self.obs.setRange("peak", 17.7, 18.3);
        for itry in range(NTRIES):
            data = self.get_random_data()
            gaussian.fitTo(data, ROOT.RooFit.Range("peak"))
            br_mean[0] = mean.getVal()
            br_sigma[0] = sigma.getVal()
            out_tree.Fill()
            hists.append(super(ROOT.RooDataSet, data).createHistogram("mass_"+str(itry), self.obs))

        fout.cd()
        for hist in hists:
            hist.Write()
        out_tree.Write()
        fout.Close()

class Fit2L:

    def __init__(self, file_name, br_name="mass"):
        self.file_name = file_name
        self.br_name = br_name
        min_x = 8
        max_x = 12
        self.obs = RooRealVar("obs", "m4l", min_x, max_x)
        self.nbins = int((max_x - min_x) * 10)
        self.obs.setBin(self.nbins)
        self.ws = ROOT.RooWorkspace("combined", "combined")

    def build_model(self):
        mean = RooRealVar("mean", "mass of 1S", 9.46, 9.2, 9.7)
        sigma = RooRealVar("sigma", "sigma of gaussian", 0.12, 0.09, 0.3)
        gaussian = ROOT.RooGaussian("gauss", "gauss", self.obs, mean, sigma)
        n_sig = RooRealVar("n_sig", "number of signal" , 5000, 0, 100000)
        esig = ROOT.RooExtendPdf("esig", "esig", gaussian, n_sig)

        m2_shift = RooRealVar("m2_shift", "m2 shift", 0.56296)
        m2 = ROOT.RooFormulaVar("m2", "mass of 2S", "mean + m2_shift", RooArgList(mean,m2_shift))
        s2 = ROOT.RooFormulaVar("s2", "sigma of 2S", "sigma*(1+m2_shift)/mean", RooArgList(sigma, m2_shift, mean))
        g2 = ROOT.RooGaussian("g2", "gauss", self.obs, m2, s2)
        n2 = RooRealVar("n2", "number of 2S" , 5000, 0, 100000)
        esig2 = ROOT.RooExtendPdf("esig2", "esig2", g2, n2)

        m3_shift = RooRealVar("m3_shift", "m3 shift", 0.8949)
        m3 = ROOT.RooFormulaVar("m3", "mass of 3S", "mean + m3_shift", RooArgList(mean, m3_shift))
        s3 = ROOT.RooFormulaVar("s3", "sigma of 3S", "sigma*(1+m3_shift)/mean", RooArgList(sigma, m3_shift, mean))
        g3 = ROOT.RooGaussian("g3", "gauss", self.obs, m3, s3)
        n3 = RooRealVar("n3", "number of 3S" , 5000, 0, 100000)
        esig3 = ROOT.RooExtendPdf("esig3", "esig3", g3, n3)

        n_bkg = RooRealVar("n_bkg", "number of bkg" , 1000, 0, 100000)
        p0 = RooRealVar("p0", "p0", -1E6, 1E6)
        p1 = RooRealVar("p1", "p1", -1E6, 1E6)
        p2 = RooRealVar("p2", "p2", -1E6, 1E6)
        p3 = RooRealVar("p3", "p3", -1E6, 1E6)
        p4 = RooRealVar("p4", "p4", -1E6, 1E6)
        bkg = ROOT.RooChebychev("bkg", "bkg", self.obs, RooArgList(p0, p1, p2, p3, p4))
        ebkg = ROOT.RooExtendPdf("ebkg", "ebkg", bkg, n_bkg)
        model = ROOT.RooAddPdf("model", "model", RooArgList(esig, esig2, esig3, ebkg))
        getattr(self.ws, "import")(model)

    def get_data(self):
        fin = ROOT.TFile.Open(self.file_name)
        tree = fin.Get("upsilon")
        nentries = tree.GetEntries()
        print "total: ", nentries
        obs_set = RooArgSet(self.obs)
        data = ROOT.RooDataSet("data", "data", obs_set)
        for ientry in xrange(nentries):
            tree.GetEntry(ientry)
            for i,m4l in enumerate(getattr(tree, self.br_name)):
                m4l = m4l/1000
                if m4l > self.obs.getMax() or m4l < self.obs.getMin():
                    continue
                if tree.chi2[i] > 10000:
                    continue

                self.obs.setVal(m4l)
                data.add(obs_set)

        getattr(self.ws, "import")(data)
        fin.Close()

    def fit(self):
        self.build_model()
        self.get_data()
        data = self.ws.obj("data")
        model = self.ws.obj("model")
        print "total data:", data.sumEntries()
        nll = model.createNLL(data)

        ## plot
        frame = self.obs.frame(self.nbins)
        data.plotOn(frame)

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

        n_sig.setVal(50000)
        n_sig.setConstant(False)
        n2.setVal(50000)
        n3.setVal(50000)
        n2.setConstant(False)
        n3.setConstant(False)
        model.fitTo(data)
        nll_uncondition = nll.getVal()

        try:
            sig = math.sqrt(2*(nll_condition - nll_uncondition))
            print "significance: ", sig
        except:
            print nll_condition,nll_uncondition

        model.plotOn(frame, ROOT.RooFit.LineColor(2))
        canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
        frame.Draw()
        #frame.GetYaxis().SetRangeUser(0, 120)
        canvas.SaveAs("fit.pdf")
        self.ws.writeToFile("combined.root")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print sys.argv[0], " file_name [br_name]"
        exit(1)

    file_name = sys.argv[1]
    br_name = "mass"
    if len(sys.argv) > 2:
        br_name = sys.argv[2]

    #fit_4l = Fit4L(file_name, br_name)
    #fit_4l.fit()
    #fit_4l.fit_signal_random()

    fit_4l = Fit2L(file_name, br_name)
    fit_4l.fit()
