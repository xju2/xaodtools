#!/usr/bin/env python
__doc__="""
This script reads the Ntuple produced by giuliano.gustavino@cern.ch
or that produced by jetsmearing.py
make some plots for each background, then plotted by plot_monojet.py
"""
import ROOT
import tools
from optparse import OptionParser

ROOT.gROOT.SetBatch()

class MonoJetNtupleReader:
    def __init__(self, file_name, tree_name, out_name):
        self.fin = ROOT.TFile.Open(file_name)
        self.tree = self.fin.Get(tree_name)
        self.out_name = out_name

    def book_hists(self):
        met_xbins = [
            250, 275, 300, 350, 400,
            450, 500, 600, 800, 1000, 1200, 2000
        ]
        # first create a template then make a clone
        # so that one can use the same template for the varaibles that share the same binning
        h_met_temp = tools.make_hist("h_met_temp", met_xbins)
        h_met_temp.SetXTitle("E_{T}^{miss} [GeV]")
        self.h_met = h_met_temp.Clone('h_met')

        h_met_temp.SetXTitle("leading jet p_{T} [GeV]")
        self.h_jetPt = h_met_temp.Clone('h_jetPt')

        h_njets_temp = ROOT.TH1F("h_njets_temp", "h_njets_temp;n_{jets}", 11, -0.5, 10.5)
        self.h_njets = h_njets_temp.Clone("h_njets")

        h_minDphi_temp = ROOT.TH1F("h_minDphi_temp", "h_minDphi_temp;min #Delta#phi(jets, E_{T}^{miss})", 16, 0, 3.2 )
        self.h_minDphi = h_minDphi_temp.Clone("h_minDphi")

    def loop_events(self, is_giuliano):
        if is_giuliano:
            print "reading Giuliano's input"
        else:
            print "reading Smearing input"

        # loop over the events and fill histograms
        nentries = self.tree.GetEntries()
        print "total events:",nentries
        tree = self.tree
        for ientry in xrange(nentries):
            tree.GetEntry(ientry)
            weight = tree.weight

            # Get variables
            if is_giuliano:
                njets = tree.n_jet
                min_dphi = tree.jet_met_dphi_min

                # apply QCD selections as suggested.
                pass_sel = tree.isQCD and tree.applyTrigger
                pass_sel = pass_sel and min_dphi < 0.4 and njets < 5
                met = tree.met_nomuon_tst_et/1E3
                pT = tree.jet_pt[0]/1E3
            else:
                njets = tree.njets
                min_dphi = tree.min_dphi
                pass_sel = njets < 5 and min_dphi < 0.4
                met = tree.met_et
                pT = tree.leading_jet_pt

            pass_sel = pass_sel and met > 250 and pT > 250
            if not pass_sel:
                continue


            self.h_met.Fill(met, weight)
            self.h_jetPt.Fill(pT, weight)
            self.h_njets.Fill(njets, weight)
            self.h_minDphi.Fill(min_dphi, weight)

    def save_hists(self):
        fout = ROOT.TFile.Open(self.out_name, "recreate")
        self.h_met.Write()
        self.h_jetPt.Write()
        self.h_njets.Write()
        self.h_minDphi.Write()

        fout.Close()
        if self.fin:
            self.fin.Close()

    def process(self, is_giuliano):
        self.book_hists()
        self.loop_events(is_giuliano)
        self.save_hists()

if __name__ == "__main__":
    usage = "%prog file_name tree_name out_name"
    parser = OptionParser(usage=usage, description="read Ntuple")
    parser.add_option("-s", "--smearing", dest="smear", default=False, action="store_true", help="input from jetsmearing!")

    (options, args) = parser.parse_args()
    if len(args) < 3:
        parser.print_help()
        exit(1)

    file_name = args[0]
    tree_name = args[1]
    out_name = args[2]

    ntuple = MonoJetNtupleReader(file_name, tree_name, out_name)
    ntuple.process(not options.smear)
