#!/usr/bin/env python

import ROOT
ROOT.gROOT.SetBatch()

import AtlasStyle
from ploter import Ploter

import os
from array import array

class Comparison:
    def __init__(self):
        self.out_root_name = "histograms.root"
        self.ps = Ploter()

        if os.path.exists(self.out_root_name):
            self.fout = ROOT.TFile.Open(self.out_root_name)
        else:
            self.fout = None

        self.cuts_on_weight = {}
        self.cuts_on_weight[364250] = 20
        self.cuts_on_weight[364256] = 30
        self.cuts_on_weight[364257] = 20
        self.cuts_on_weight[364258] = 40
        self.cuts_on_weight[364259] = 100
        self.cuts_on_weight[364260] = 25

    def has_large_weight(self, tree):
        try:
            return tree.w_MCw > self.cuts_on_weight[tree.run]
        except KeyError:
            return False

    def get_hist(self, tree, name, cut):
        if self.fout:
            return self.fout.Get(name)

        #h1 = ROOT.TH1F(name, name, 70, 130, 1530)
        # define a histogram with various bin-width
        bin_list = []
        mass = 200
        while mass <= 1500:
            bin_list.append(mass)
            if mass < 300:
                mass += 20
            elif mass < 600:
                mass += 25
            elif mass < 100:
                mass += 50
            else:
                mass += 100
        h1 = ROOT.TH1F(name, name, len(bin_list)-1, array('f', bin_list))

        tree.SetBranchStatus("*", 0)
        tree.SetBranchStatus("run", 1)
        tree.SetBranchStatus("w_MCw", 1)
        tree.SetBranchStatus("higgs_m_fidBorn_4lsel", 1)
        tree.SetBranchStatus("dijet_m_fidBorn_4lsel", 1)
        tree.SetBranchStatus("dijet_deltaeta_fidBorn_4lsel", 1)
        tree.SetBranchStatus("event_type_fidBorn_4lsel", 1)
        total_weight = 0
        for ientry in xrange(tree.GetEntries()):
            tree.GetEntry(ientry)
            if self.has_large_weight(tree):
                continue

            total_weight += tree.w_MCw
            if tree.higgs_m_fidBorn_4lsel == -999:
                continue

            pass_VBF = tree.dijet_m_fidBorn_4lsel > 400 and abs(tree.dijet_deltaeta_fidBorn_4lsel) > 3.3
            pass_cuts = False
            event_type = tree.event_type_fidBorn_4lsel
            #if cut == 1 and not pass_VBF and event_type == 0:
            #    # ggF 4mu
            #    pass_cuts = True
            #elif cut == 2 and not pass_VBF and event_type == 1:
            #    # ggF 4e
            #    pass_cuts = True
            #elif cut == 3 and not pass_VBF and (event_type == 2 or event_type == 3):
            #    # ggF 2mu2e
            #    pass_cuts = True
            #elif cut == 4 and pass_VBF:
            #    pass_cuts = True
            #elif cut == -1:
            #    pass_cuts = True
            #else:
            #    pass

            if pass_VBF:
                h1.Fill(tree.higgs_m_fidBorn_4lsel, tree.w_MCw)

        print "total weight:", total_weight
        if h1.GetEntries() == 0 or h1.GetIntegral() == 0:
            print h1.GetName(),"is empty!"
            exit(1)

        h1.Scale(1./total_weight)
        return h1

    def go(self):
        base_dir = "/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/MiniTrees/Prod_v12/mc/TheorySystematics/fidBorn/"
        tree_name = "tree_incl_all"

        f1_name = "mc15_13TeV.364250.Sherpa_222_NNPDF30NNLO_llll.root"
        f1_1_name = "mc15_13TeV.364251.Sherpa_222_NNPDF30NNLO_llll_m4l100_300_filt100_150.root"
        f1_2_name = "mc15_13TeV.364252.Sherpa_222_NNPDF30NNLO_llll_m4l300.root"

        f2_name = "mc15_13TeV.364256.Sherpa_222_NNPDF30NNLO_llll_CKKW15.root"
        f3_name = "group.phys-higgs.user.sabidi.364257.Sherpa_222_NNPDF30NNLO_llll_CKKW30.TRUTH4.root"

        f4_name = "group.phys-higgs.user.sabidi.364258.Sherpa_222_NNPDF30NNLO_llll_QSF025.TRUTH4.root"
        f5_name = "mc15_13TeV.364259.Sherpa_222_NNPDF30NNLO_llll_QSF4.root"

        f6_name = "group.phys-higgs.user.sabidi.364260.Sherpa_222_NNPDF30NNLO_llll_CSSKIN.TRUTH4.root"

        #f1 = ROOT.TFile.Open(base_dir+f1_name)
        #t1 = f1.Get(tree_name)
        t1 = ROOT.TChain(tree_name, tree_name)
        t1.Add(base_dir+f1_name)
        #t1.Add(base_dir+f1_1_name)
        #t1.Add(base_dir+f1_2_name)

        f2 = ROOT.TFile.Open(base_dir+f2_name)
        t2 = f2.Get(tree_name)

        f3 = ROOT.TFile.Open(base_dir+f3_name)
        t3 = f3.Get(tree_name)

        f4 = ROOT.TFile.Open(base_dir+f4_name)
        t4 = f4.Get(tree_name)

        f5 = ROOT.TFile.Open(base_dir+f5_name)
        t5 = f5.Get(tree_name)

        f6 = ROOT.TFile.Open(base_dir+f6_name)
        t6 = f6.Get(tree_name)

        cuts_dic = {}
        #cuts_dic["inc"] = -1
        #cuts_dic["4mu"] = 1
        #cuts_dic["4e"] = 2
        #cuts_dic["2e2mu"] = 3
        cuts_dic["VBF"] = 4

        all_histograms = []
        for key,value in cuts_dic.iteritems():
            print key,value
            h1 = self.get_hist(t1, 'h_nominal_'+key, value)
            h2 = self.get_hist(t2, 'h_ckkw15_'+key, value)
            h3 = self.get_hist(t3, 'h_ckkw30_'+key, value)
            h4 = self.get_hist(t4, 'h_qsf0p25_'+key, value)
            h5 = self.get_hist(t5, 'h_qsf4_'+key, value)
            h6 = self.get_hist(t6, 'h_csskin_'+key, value)

            all_histograms.append(h1)
            all_histograms.append(h2)
            all_histograms.append(h3)
            all_histograms.append(h4)
            all_histograms.append(h5)
            all_histograms.append(h6)
            for hist in all_histograms:
                hist.Rebin(4)
                hist.Scale(1./hist.Integral())

            list_ckkw = [h1, h2, h3]
            tag_ckkw = ["Nominal", "ckkw 15", "ckkw 30"]
            opt = {}
            opt['out_name'] = "sys_ckkw_"+key
            opt['add_yields'] = True
            opt['no_fill'] = True
            opt['ratio_title'] = "Variation/Nominal"
            self.ps.compare_hists(list_ckkw, tag_ckkw, **opt)
            #Ressummation scale
            list_qfs = [h1, h4, h5]
            tag_qfs = ["Nominal", "QSF 1/4", "QSF 4"]
            opt['out_name']= 'sys_qsf_'+key
            self.ps.compare_hists(list_qfs, tag_qfs, **opt)
            #Catani-Seymour shower option
            list_csskin = [h1, h6]
            tag_csskin = ["Nominal", "CSSKIN"]
            opt['out_name']= 'sys_csskin_'+key
            self.ps.compare_hists(list_csskin, tag_csskin, **opt)
            #everything
            list_all = [h1, h2, h3, h4, h5, h6]
            tag_all = ["Nominal", "ckkw 15", "ckkw 30", "QSF 1/4", "QSF 4", "CSSKIN"]
            opt['out_name']= 'sys_all_'+key
            self.ps.compare_hists(list_all, tag_all, **opt)

        if not self.fout:
                self.fout = ROOT.TFile.Open(self.out_root_name, "recreate")
                for hist in all_histograms:
                    hist.Write()
                self.fout.Close()

        #f1.Close()
        f2.Close()
        f3.Close()
        f4.Close()
        f5.Close()
        f6.Close()


if __name__ == "__main__":
    cmp = Comparison()
    cmp.go()
