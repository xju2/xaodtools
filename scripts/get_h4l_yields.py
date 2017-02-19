#!/usr/bin/env python
"""
get the yields for H4l analysis,
inputs are minitress
"""

import ROOT
from optparse import OptionParser

ROOT.gROOT.SetBatch()

class MinitreeReader():
    def __init__(self, file_name, options):
        print "file name:", file_name
        #self.fin = ROOT.TFile.Open(file_name)
        TREE_NAME = "tree_incl_all"
        self.tree = ROOT.TChain(TREE_NAME, TREE_NAME)
        self.tree.Add(file_name)
        self.options = options
        #self.lumi = 36.5

    def get_highmass_cuts(self):
        dic = {
            "ggF_2e2mu":"pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<1300&&(event_type==3||event_type==2) && prod_type_HM==0",
            "ggF_4e": "pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<1300&&(event_type==1) && prod_type_HM==0",
            "ggF_4mu": "pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<1300&&(event_type==0) && prod_type_HM==0",
            "VBF_incl": "pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<1300&&prod_type_HM==1",
        }
        return dic

    def get_lowmass_cuts(self):
        dic = {
            "ggF_2e2mu_13TeV":"pass_vtx4lCut==1 && 115<m4l_constrained && m4l_constrained<130 && event_type==3",
            "ggF_2mu2e_13TeV":"pass_vtx4lCut==1 && 115<m4l_constrained && m4l_constrained<130 && event_type==2",
            "ggF_4e_13TeV"   :"pass_vtx4lCut==1 && 115<m4l_constrained && m4l_constrained<130 && event_type==1",
            "ggF_4mu_13TeV"  :"pass_vtx4lCut==1 && 115<m4l_constrained && m4l_constrained<130 && event_type==0",
        }
        return dic

    def get_yield(self, dic):
        self.tree.GetEntry(0)
        print "Luminosity:", round(self.tree.w_lumi, 2),"fb-1"
        # use the one in the minitree
        self.lumi = self.tree.w_lumi
        for channel, cuts in dic.iteritems():
            # use weight_jet when using jet information!
            cut = ROOT.TCut(self.weight_name+"*("+cuts+")")
            self.tree.Draw("m4l_constrained_HM>>h1", cut)
            yields = ROOT.gDirectory.Get("h1").Integral()
            print channel,"{:.4f}  {:.4f}".format(yields, yields/self.lumi)

    def process(self):
        cuts_dir = None
        if self.options.analysis == "HighMass":
            cuts_dir = self.get_highmass_cuts()
            self.weight_name = "weight_jet"
        elif self.options.analysis == "LowMass":
            cuts_dir = self.get_lowmass_cuts()
            self.weight_name = "weight"
        else:
            pass

        if cuts_dir is not None:
            self.get_yield(cuts_dir)
        else:
            print "I don't know"

if __name__ == "__main__":
    usage = "%prog [options] file_name"
    version="%prog 1.0"
    parser = OptionParser(usage=usage, description="get yields for WS", version=version)
    parser.add_option("--analysis", dest='analysis', default='HighMass')

    (options,args) = parser.parse_args()
    if len(args) < 1:
        print parser.print_help()
        exit(1)

    reader = MinitreeReader(args[0], options)
    reader.process()
