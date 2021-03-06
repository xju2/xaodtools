#!/usr/bin/env python
"""
It's likely the fit cannot converge out of box,
need to tune the model in build_model()
"""

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
if not hasattr(ROOT, "myText"):
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
        # with 3mu4 trigger or not
        self.with_3mu4 = True
        # only unprescaled runs
        self.only_unprescaled = True

        # different onia_pt cuts, 0:noCut, 1:<5, 2:5-10, 3:10,20, 4:>20
        self.onia_pt_cut = 0
        self.onia_pt_cuts = [0, 5, 10, 20]
        #self.onia_pt_cuts = [0, 10, 20]

        # different muon pT cuts
        # suggested by Terry to look at onia mass that have one muon with pT (3, 4) GeV.
        self.use_low_pt_muon = False

        # it is found that these low pT muons contributes a lot the background 40%
        # while only gain 15% signal
        self.no_low_pt = False

        # if build new model
        self.new_model = False
        self.model_name = "model"

        # workspace components
        self.data = None
        self.model = None

    def pass_onia_pt_cut(self, pT):
        pT_cut = int(self.onia_pt_cut)
        if pT_cut >= len(self.onia_pt_cuts):
            return pT >= self.onia_pt_cuts[ len(self.onia_pt_cuts)-1 ]
        elif pT_cut > 0:
            return pT >= self.onia_pt_cuts[pT_cut-1] and pT < self.onia_pt_cuts[pT_cut]
        else:
            return True

    def print_onia_pt_cut(self):
        pT_cut = int(self.onia_pt_cut)
        if pT_cut >= len(self.onia_pt_cuts):
            return "p_{T}^{onia} in ["+str(self.onia_pt_cuts[ len(self.onia_pt_cuts)-1 ])+", #infty) GeV"
        elif pT_cut > 0:
            return "p_{T}^{onia} in ["+str(self.onia_pt_cuts[pT_cut-1])+", "+str(self.onia_pt_cuts[pT_cut])+") GeV"
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
        res = "ws_chi2Cut"+str(self.chi2_cut)+"_DiOnia"+str(self.dionia_selection)+"_pTCut"+str(self.onia_pt_cut)+"_withTrigger"+str(self.with_3mu4)
        if self.use_low_pt_muon:
            res += "_LowPtMuon"
        if self.no_low_pt:
            res += "_noLowPtMuon"

        res += ".root"
        return res

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
        n3 = RooRealVar("n3", "number of 3S" , 1, 0, 2000)
        #n3 = ROOT.RooFormulaVar("n3", "number of 3S", "n2*0.45", RooArgList(n2))
        esig3 = ROOT.RooExtendPdf("esig3", "esig3", g3, n3)

        n_bkg = RooRealVar("n_bkg", "number of bkg" , 1000, 0, 1E7)
        p0 = RooRealVar("p0", "p0", -1E6, 1E6)
        p1 = RooRealVar("p1", "p1", -1E6, 1E6)
        p2 = RooRealVar("p2", "p2", -1E6, 1E6)
        p3 = RooRealVar("p3", "p3", -1E6, 1E6)
        p4 = RooRealVar("p4", "p4", -1E6, 1E6)
        #p0 = RooRealVar("p0", "p0", 9.93810e-02)
        #p1 = RooRealVar("p1", "p1", -3.51406e-02)
        #p2 = RooRealVar("p2", "p2", 5.92968e-03)
        #p3 = RooRealVar("p3", "p3", -6.53710e-03)
        #p4 = RooRealVar("p4", "p4", 9.76852e-03)

        bkg = ROOT.RooChebychev("bkg", "bkg", self.obs, RooArgList(p0, p1, p2, p3, p4))
        #bkg = ROOT.RooChebychev("bkg", "bkg", self.obs, RooArgList(p0, p1, p2))
        #bkg = ROOT.RooPolynomial("bkg", "bkg", self.obs, RooArgList(p0, p1, p2))
        ebkg = ROOT.RooExtendPdf("ebkg", "ebkg", bkg, n_bkg)
        model = ROOT.RooAddPdf(self.model_name, self.model_name, RooArgList(esig, esig2, esig3, ebkg))
        getattr(self.ws, "import")(model, ROOT.RooFit.RecycleConflictNodes())
        #getattr(self.ws, "import")(model, ROOT.RooFit.RenameConflictNodes("new"))

    def is_unprescaled_runs(self, tree):
        run_ = tree.run
        res = True
        if run_ >= 297730 and run_ < 307619:
            res = False
        return res

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
            if self.only_unprescaled and not self.is_unprescaled_runs(tree):
                continue

            for i,m4l in enumerate(getattr(tree, self.br_name)):
                #m4l = m4l/1000
                if m4l > self.obs.getMax() or m4l < self.obs.getMin():
                    continue
                # apply chi2 
                if tree.chi2[i] > self.chi2_cut:
                    continue
                # apply trigger 
                if self.with_3mu4 and not tree.trig_3mu4:
                    continue
                # apply dionia selection
                if self.dionia_selection and tree.pass_diOnia != 1:
                    continue

                has_one_lowPt = (tree.mu_pt_1[i] > 3E3 and tree.mu_pt_1[i] < 4E3) or (tree.mu_pt_2[i] > 3E3 and tree.mu_pt_2[i] < 4E3)
                if self.use_low_pt_muon and not has_one_lowPt:
                    continue

                if self.no_low_pt and has_one_lowPt:
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
            # try to look for the workspace and reused the data
            f1 = ROOT.TFile.Open(self.get_ws_name())
            self.ws = ROOT.RooWorkspace("combined", "combined")
            if f1:
                self.f1 = f1
                ws = f1.Get("combined")
                if ws:
                    getattr(self.ws, "import")(ws.obj('data'))
                if self.new_model:
                    print "building new models"
                    self.build_model()
                else:
                    getattr(self.ws, "import")(ws.obj(self.model_name))
            else:
                self.build_model()
                self.get_data()

        data = self.ws.obj("data")
        self.obs.setRange("fit_range", 8.2, 11.7)
        model = self.ws.obj(self.model_name)
        model.Print()
        print "total data:", data.sumEntries()
        nll = model.createNLL(data)

        #self.change_model()

        model.fitTo(data, ROOT.RooFit.Range("fit_range"))
        #model.fitTo(data)

        nll_uncondition = nll.getVal()
        self.ws.saveSnapshot("splusb", self.ws.allVars())

        if hasattr(self, "f1"):
            self.f1.Close()
        else:
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
        mass_offset = 0.18
        if self.use_low_pt_muon:
            mass_offset = 0.7
        ROOT.myText(mass_offset, 0.87, 1, "m = {:.2f} GeV".format(self.ws.var("mean").getVal()))
        ROOT.myText(mass_offset, 0.87-0.05, 1, "#sigma = {:.2f} GeV".format(self.ws.var("sigma").getVal()))
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
        canvas.SaveAs(self.get_ws_name().replace("root", "eps").replace("ws_", "fit_"))
        self.ws.obj("s2").Print()
        self.ws.obj("s3").Print()

    def change_model(self):
        self.ws.loadSnapshot("splusb")
        # free signal
        mean_ = self.ws.var("mean")
        mean_.setMin(9.2)
        mean_.setMax(9.7)
        mean_.setConstant(False)
        sigma_ = self.ws.var("sigma")
        sigma_.setMin(0.09)
        sigma_.setMax(0.3)
        sigma_.setConstant(False)
        # fix background
        self.ws.var("p0").setConstant(True)
        self.ws.var("p1").setConstant(True)
        self.ws.var("p2").setConstant(True)
        if self.ws.var('p3'):
            self.ws.var('p3').setConstant(True)
        if self.ws.var('p4'):
            self.ws.var('p4').setConstant(True)


if __name__ == "__main__":
    usage = sys.argv[0]+" file_name br_name"
    parser = OptionParser(usage=usage)
    parser.add_option('--dionia', action="store_true", dest="dionia", help="apply dionia selection", default=False)
    parser.add_option('--oniaPt', dest='oniapt', help="onia pT cut, 0/1/2/3/4", default=0)
    parser.add_option('--noTrigger', action="store_true", dest='notrigger', help="don't apply trigger", default=False)
    parser.add_option('--title', dest='title', help="title of x-axis", default="m_{#mu#mu} [GeV]")
    parser.add_option('--allRuns', dest='allruns', help="use all runs", default=False, action="store_true")
    parser.add_option('--newModel', dest='model', help="build new model, regardless of exiting one in the workspace", default=False, action="store_true")
    parser.add_option('--lowPtMuon', dest='lowPt', help="use the onia with at least one muon with pT of (3, 4) GeV", default=False, action="store_true")
    parser.add_option('--noLowPtMuon', dest='noLowPt', help="use the onia that both muon with pT > 4 GeV", default=False, action="store_true")
    #parser.add_option('--chi2', dest='chi2', help="chi2 cut applied in onium", default=None)


    (options, args) = parser.parse_args()
    if len(args) < 2:
        parser.print_help()
        exit(1)

    file_name = args[0]
    br_name = "mass"
    if len(sys.argv) > 2:
        br_name = args[1]

    #cuts = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 8, 10, 10000]
    cuts = [3.0]
    for cut in cuts:
        fit_4l = Fit2L(file_name, br_name)
        if options.dionia:
            fit_4l.dionia_selection = True

        if options.notrigger:
            fit_4l.with_3mu4 = False

        if options.allruns:
            fit_4l.only_unprescaled = False

        if options.model:
            fit_4l.new_model = True

        if options.lowPt and options.noLowPt:
            print "what do you want?? with low pT or without? chose one!"
            exit(2)

        if options.lowPt:
            fit_4l.use_low_pt_muon = True

        if options.noLowPt:
            fit_4l.no_low_pt = True

        fit_4l.onia_pt_cut = options.oniapt
        fit_4l.chi2_cut = cut
        fit_4l.fit()
        fit_4l.plot(options.title)
