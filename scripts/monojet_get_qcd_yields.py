#!/usr/bin/env python
"""
from smeared data to get fraction of yields in each MET bin,
print out the yields in each MET bins
"""
import tools
import ROOT
from optparse import OptionParser
import os

ROOT.gROOT.SetBatch()
class MonoJetQCDYields:
    def __init__(self):
        self.hist_input = "met_signal_fraction.root"
        self.hist_name = "h_met"
        pass

    def get_fractions(self, file_name):
        if os.path.exists(self.hist_input):
            return
        f1 = ROOT.TFile.Open(file_name)
        tree = f1.Get("smeared")

        met_xbins = [
            250, 300, 350, 400, 500, 600, 700,
            800, 900, 2000
        ]
        h_met = tools.make_hist(self.hist_name, met_xbins)
        cuts = ROOT.TCut("weight*(njets < 5 && met_et > 250 && leading_jet_pt > 250 && min_dphi > 0.4)")
        tree.Draw("met_et>>"+h_met.GetName(), cuts)

        fout = ROOT.TFile.Open(self.hist_input, "recreate")
        h_met.Scale(1./h_met.Integral())
        h_met.Write()
        fout.Close()
        f1.Close()

    def print_yileds(self, n_total):
        fin = ROOT.TFile.Open(self.hist_input)
        if not fin:
            print "please run get_fraction first"
            return
        hist = fin.Get(self.hist_name)
        for ibin in range(hist.GetXaxis().GetNbins()):
            low = hist.GetBinLowEdge(ibin+1)
            hi = hist.GetBinLowEdge(ibin+2)
            print ibin+1,low,hi,hist.GetBinContent(ibin+1)*n_total

        fin.Close()


if __name__ == "__main__":
    usage = "%prog n_total file_name "
    parser = OptionParser(usage=usage, description=__doc__)

    (options, args) = parser.parse_args()
    if len(args) < 1:
        print parser.print_help()
        exit(1)

    n_total = float(args[0])

    qcd_yields = MonoJetQCDYields()
    if len(args) > 1:
        qcd_yields.get_fractions(args[1])

    qcd_yields.print_yileds(n_total)
