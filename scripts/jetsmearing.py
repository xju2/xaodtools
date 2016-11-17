#!/usr/bin/env python

import ROOT
import AtlasStyle
from optparse import OptionParser

import os
import sys
import threading
import math
from array import array
import tools
from ploter import Ploter

ROOT.gROOT.SetBatch()

class MiniTree:
    def __init__(self):
        self.m_debug = False
        self.save_seed = False

    @staticmethod
    def get_met_sig(chain):
        return (chain.MET_et/1E3 - 8)/math.sqrt(chain.MET_sumet/1E3)
        #return (chain.MET_et/1E3)/math.sqrt(chain.MET_sumet/1E3)

    def change_file(self, file_name, out_name):
        if not hasattr(ROOT, "loader"):
            ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c")

        if not hasattr(ROOT, "SmearedInfo"):
            ROOT.gROOT.Macro(os.getenv("ROOTCOREDIR")+"/scripts/load_packages.C")

        chain = ROOT.loader(file_name, "physics")
        nentries = chain.GetEntries()

        # book branches for output trees
        outfile = ROOT.TFile.Open(out_name, "recreate")
        outtree = ROOT.TTree("smeared", "smeared")
        seedtree = ROOT.TTree("seed", "seed")
        met = array('f', [0])
        sumet = array('f', [0])
        jetpt = array('f', [0])
        jeteta = array('f', [0])
        jetphi = array('f', [0])
        subjetpt = array('f', [0])
        subjeteta = array('f', [0])
        subjetphi = array('f', [0])

        njets = array('i', [0])
        dphi = array('f', [0])
        dphiEP = array('f', [0])
        rmet_pt = array('f', [0])
        run_number = array('i', [0])
        event_number = array('i', [0])
        lb = array('i', [0])
        weight = array('f', [0])

        ht = array('f', [0])
        l3rd_jet_pt = array('f', [0])
        l3rd_jet_eta = array('f', [0])
        l3rd_jet_phi = array('f', [0])
        l4th_jet_pt = array('f', [0])
        l4th_jet_eta = array('f', [0])
        l4th_jet_phi = array('f', [0])
        n_vertices = array('i', [0])
        mass_eff = array('f', [0])


        outtree.Branch('met_et', met, 'met_et/F')
        outtree.Branch('sumet', sumet, 'sumet/F')
        outtree.Branch('njets', njets, 'njets/I')
        outtree.Branch('leading_jet_pt', jetpt, 'leading_jet_pt/F')
        outtree.Branch('leading_jet_eta', jeteta, 'leading_jet_eta/F')
        outtree.Branch('sub_leading_jet_pt', subjetpt, 'sub_leading_jet_pt/F')
        outtree.Branch('sub_leading_jet_eta', subjeteta, 'sub_leading_jet_eta/F')
        outtree.Branch('min_dphi', dphi, 'min_dphi/F')
        outtree.Branch("dphi_ep", dphiEP, 'dphi_ep/F')
        outtree.Branch('rmet_pt', rmet_pt, 'rmet_pt/F')
        outtree.Branch('run', run_number, 'run/I')
        outtree.Branch('event', event_number, 'event/I')
        outtree.Branch('lb', lb, 'lb/I')
        outtree.Branch('weight', weight, 'weight/F')

        outtree.Branch('Ht', ht, 'Ht/F')
        outtree.Branch('leading_jet_phi', jetphi, 'leading_jet_phi/F')
        outtree.Branch('sub_leading_jet_phi', subjetphi, 'sub_leading_jet_phi/F')
        outtree.Branch('l3rd_jet_pt', l3rd_jet_pt, 'l3rd_jet_pt/F')
        outtree.Branch('l3rd_jet_eta', l3rd_jet_eta, 'l3rd_jet_eta/F')
        outtree.Branch('l3rd_jet_phi', l3rd_jet_phi, 'l3rd_jet_phi/F')
        outtree.Branch('l4th_jet_pt', l4th_jet_pt, 'l4th_jet_pt/F')
        outtree.Branch('l4th_jet_eta', l4th_jet_eta, 'l4th_jet_eta/F')
        outtree.Branch('l4th_jet_phi', l4th_jet_phi, 'l4th_jet_phi/F')
        outtree.Branch('n_vertices', n_vertices, 'n_vertices/I')
        outtree.Branch('mass_eff', mass_eff, 'mass_eff/F')

        seedtree.Branch('met_et', met, 'met_et/F')
        seedtree.Branch('njets', njets, 'njets/I')
        seedtree.Branch('leading_jet_pt', jetpt, 'leading_jet_pt/F')
        seedtree.Branch('leading_jet_eta', jeteta, 'leading_jet_eta/F')
        seedtree.Branch('min_dphi', dphi, 'min_dphi/F')
        #seedtree.Branch("dphi_ep", dphiEP, 'dphi_ep/F')
        #seedtree.Branch('rmet_pt', rmet_pt, 'rmet_pt/F')
        seedtree.Branch('run', run_number, 'run/I')
        seedtree.Branch('event', event_number, 'event/I')
        seedtree.Branch('lb', lb, 'lb/I')
        seedtree.Branch('weight', weight, 'weight/F')

        met_sig_br = array('f', [0])
        frac_soft_br = array('f', [0])
        seedtree.Branch("met_sig", met_sig_br, 'met_sig/F')
        seedtree.Branch("frac_soft", frac_soft_br, 'frac_soft/F')

        # book histograms for all-events, seed events (i.e. after METsig cut) and smeared events..
        # Then compare them to see the "bias" of the seed selection.
        met_xbins = [
            # 0, 50, 100, 125, 150, 175,
            # 200, 225, 
            250, 275, 300, 350, 400,
            450, 500, 600, 800, 1000
        ]
        h_met_temp = tools.make_hist("h_met_temp", met_xbins)
        h_met_temp.SetXTitle("E_{T}^{miss} [GeV]")
        h_met_seed = h_met_temp.Clone('h_met_seed')
        h_met_smeared = h_met_temp.Clone('h_met_smeared')
        h_met_all = h_met_temp.Clone('h_met_all')

        h_met_temp.SetXTitle("leading jet p_{T} [GeV]")
        h_jetPt_seed = h_met_temp.Clone('h_jetPt_seed')
        h_jetPt_smeared = h_met_temp.Clone('h_jetPt_smeared')
        h_jetPt_all = h_met_temp.Clone('h_jetPt_all')

        h_njets_temp = ROOT.TH1F("h_njets_temp", "h_njets_temp;n_{jets}", 11, -0.5, 10.5)
        h_njets_seed = h_njets_temp.Clone("h_njets_seed")
        h_njets_smeared = h_njets_temp.Clone("h_njets_smeared")
        h_njets_all = h_njets_temp.Clone("h_njets_all")

        h_minDphi_temp = ROOT.TH1F("h_minDphi_temp", "h_minDphi_temp;min #Delta#phi(jets, E_{T}^{miss})", 16, 0, 3.2 )
        h_minDphi_seed = h_minDphi_temp.Clone("h_minDphi_seed")
        h_minDphi_smeared = h_minDphi_temp.Clone("h_minDphi_smeared")
        h_minDphi_all = h_minDphi_temp.Clone("h_minDphi_all")
        #########################

        for ientry in xrange(nentries):
            if chain.LoadTree(ientry) < 0:
                break
            chain.GetEntry(ientry)

            jet_pt_origin = chain.jet_p4[0].Pt()

            # get trigger weight from the branch
            weight[0] = chain.triggerWeight

            # lepton/photon veto has been applied in mini-tree

            #select the seed events
            met_sig =  self.get_met_sig(chain)
            frac_soft = chain.MET_et_soft/chain.MET_et
            frac_soft_br[0] = frac_soft

            met_sig_br[0] = met_sig
            run_number[0] = chain.RunNumber
            event_number[0] = chain.EventNumber
            lb[0] = chain.lumiblock
            met[0] = chain.MET_et/1E3
            jetpt[0] = jet_pt_origin/1E3
            jeteta[0] = chain.jet_p4[0].Eta()
            njets[0] = chain.n_good_jet
            dphi[0] = chain.min_dphi_jetMET

            # don't save seed tree unless it's really needed!
            # because there are too many events...
            if self.save_seed:
                seedtree.Fill()

            # Fill all these events into all (leading pT > 250 GeV, already applied)
            # that's before the METsig cut. (unbiased)
            h_met_all.Fill(met[0], weight[0])
            h_jetPt_all.Fill(jetpt[0], weight[0])
            h_njets_all.Fill(njets[0], weight[0])
            h_minDphi_all.Fill(dphi[0], weight[0])

            #n_vertices[0] = chain.n_vertices

            if met_sig < 0.5 + 0.1*chain.n_jet_btagged:
                # METsig cut is recommended at
                # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetSmearing


                # events that passed METsig cuts, treated as seed events 
                h_met_seed.Fill(met[0], weight[0])
                h_jetPt_seed.Fill(jetpt[0], weight[0])
                h_njets_seed.Fill(njets[0], weight[0])
                h_minDphi_seed.Fill(dphi[0], weight[0])

                for data in chain.pseudoData:
                    if self.m_debug:
                        print data.leading_jet_pt_
                    met[0] = data.met_/1E3
                    jetpt[0] = data.leading_jet_pt_/1E3
                    #if met[0] < 200:
                    #    continue
                    if jetpt[0] < 200 or met[0] < 80:
                        continue

                    # these are smeared data...
                    h_met_smeared.Fill(met[0], weight[0])
                    h_jetPt_smeared.Fill(jetpt[0], weight[0])
                    h_njets_smeared.Fill(njets[0], weight[0])
                    h_minDphi_smeared.Fill(dphi[0], weight[0])

                    sumet[0] = data.sum_et_/1E3
                    jeteta[0] = data.leading_jet_eta_
                    subjetpt[0] = data.sub_leading_jet_pt_/1E3
                    subjeteta[0] = data.sub_leading_jet_eta_
                    njets[0] = data.n_good_jets_
                    dphi[0] = data.min_jets_met_
                    dphiEP[0] = data.dphi_EP_
                    rmet_pt[0] = data.met_/data.leading_jet_pt_

                    ht[0] = data.HT_/1E3
                    jetphi[0] = data.leading_jet_phi_
                    subjetphi[0] = data.sub_leading_jet_phi_
                    l3rd_jet_pt[0] = data.l3rd_jet_pt_/1E3
                    l3rd_jet_eta[0] = data.l3rd_jet_eta_
                    l3rd_jet_phi[0] = data.l3rd_jet_phi_
                    l4th_jet_pt[0] = data.l4th_jet_pt_/1E3
                    l4th_jet_eta[0] = data.l4th_jet_eta_
                    l4th_jet_phi[0] = data.l4th_jet_phi_
                    mass_eff[0] = jetpt[0]+subjetpt[0]+l3rd_jet_pt[0]+l4th_jet_pt[0]+met[0]
                    outtree.Fill()

        outfile.cd()
        outtree.Write()
        if self.save_seed:
            seedtree.Write()

        # write out histograms
        h_met_seed.Write()
        h_met_smeared.Write()
        h_met_all.Write()

        h_jetPt_seed.Write()
        h_jetPt_smeared.Write()
        h_jetPt_all.Write()

        h_njets_seed.Write()
        h_njets_smeared.Write()
        h_njets_all.Write()

        h_minDphi_seed.Write()
        h_minDphi_smeared.Write()
        h_minDphi_all.Write()

        outfile.Close()

    def compare_all_seed_smeared(self, file_name):
        # use the hitogram produced from change_file()
        fin = ROOT.TFile.Open(file_name)
        plot = Ploter()
        hist_names = ["met", "jetPt", "njets", "minDphi"]
        tag_names = ["all", "seed", "smeared"]
        for hist_name in hist_names:
            h_all = fin.Get("h_"+hist_name+"_all")
            h_seed = fin.Get("h_"+hist_name+"_seed")
            h_smeared = fin.Get("h_"+hist_name+"_smeared")
            plot.compare_hists([h_all, h_seed, h_smeared], tag_names, "check_"+hist_name)



if __name__ == "__main__":
    usage = "%prog file_name [out_name]"
    parser = OptionParser(usage=usage, description="analyze jet smearing data")
    parser.add_option("-v", "--verbose", dest='debug', default=False, help="swich to debug mode", action="store_true")
    parser.add_option("--saveSeed", dest='save_seed', default=False, help="save seed events", action="store_true")
    parser.add_option("--get_hist", dest='change', default=False, help="get smeared histograms and tree, need out_name", action="store_true")
    parser.add_option("--plot_hist", dest='plot', default=False, help="plot histograms", action="store_true")

    (options, args) = parser.parse_args()

    if options.change and len(args) < 2:
        parser.print_help()
        exit(1)

    if options.plot and len(args) < 1:
        parser.print_help()
        exit(2)

    if not options.change and not options.plot:
        print "provide a option!"
        parser.print_help()
        exit(3)


    in_name = args[0]

    minitree = MiniTree()
    if options.save_seed:
        minitree.save_seed = True

    if options.debug:
        minitree.m_debug = True

    if options.change:
        out_name = args[1]
        minitree.change_file(in_name, out_name)

    if options.plot:
        minitree.compare_all_seed_smeared(in_name)
