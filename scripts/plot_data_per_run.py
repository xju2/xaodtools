#!/usr/bin/env python
"""
plot visible cross section for data as function of run-number
inputs:
    @minitree for data
    @XML that contains luminiosity per run
"""
BASE_NAME = "/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/MiniTrees/Prod_v11/data/Nominal/data_13TeV.root"

import sys
import math
from array import array
from sets import Set

import AtlasStyle
import ROOT

from get_lumi import LumiCal
from data_period import PERIOD2016
import helper as Helper

ROOT.gROOT.SetBatch()

class EventPerRun:
    def __init__(self, file_name, out_name):
        self.f1 = ROOT.TFile.Open(file_name)

        # setup the luminosity calculator
        lumi_table_2015 = "/afs/cern.ch/user/x/xju/work/myRooCoreTools/MyXAODTools/scripts/lumitable_data15_all.csv"
        lumi_table_2016 = "/afs/cern.ch/user/x/xju/work/myRooCoreTools/MyXAODTools/" \
                          "scripts/lumitable_data16_297730_311481_v88_pro20-21.csv"
        self.cal = LumiCal()
        self.cal.load_lumi(lumi_table_2015)
        self.cal.load_lumi(lumi_table_2016)

        # number of events observed for each run.
        self.num_events_dic = {}
        for key in self.cal.lumi_dic.keys():
            self.num_events_dic[key] = 0


    def loop_events(self):
        runs_not_in_GRL = []
        ntotal = 0

        tree = self.f1.Get("tree_incl_all")
        nentries = tree.GetEntries()
        for ientry in range(nentries):
            tree.GetEntry(ientry)
            run_ = tree.run
            event_ = tree.event
            m4l = tree.m4l_constrained_HM
            event_type = tree.event_type

            if not tree.pass_vtx4lCut: 
                print run_,event_,m4l,event_type,"not pass vertex"
                continue
            # only look at 4e candidates in [230, 250] GeV
            if not event_type == 1 or not tree.weight:
                continue
            if not (230 < m4l < 250):
                continue
        
            ntotal += 1
            if run_ in self.cal.lumi_dic:
                self.num_events_dic[run_] += 1
            else:
                runs_not_in_GRL.append(run_)

        print "Total number of events in 4e", ntotal

    def events(self, start, end):
        """
        total number of events for run in range of (start, end)
        """
        total = 0
        for key, value in self.num_events_dic.iteritems():
            if start <= key <= end:
                total += value
        return total

    def per_period(self):
        period_ = [x+1 for x in range(len(PERIOD2016)-1)]
        xs_ = []
        xs_E_ = []

        lumi_list = []
        n_event_list = []
        curr_lumi = 0
        curr_nEvt = 0
        for start, end in PERIOD2016:
            # print start,end
            lumi, total_runs = self.cal.get_lumi_for_range(start, end)
            lumi /= 1E3
            n_event = self.events(start, end)

            # count accumulated number of events and luminosity
            curr_lumi += lumi
            curr_nEvt += n_event

            if lumi == 0:
                continue

            xs, xs_E = Helper.get_eff_error(
                n_event, math.sqrt(n_event),
                lumi, 0.032*lumi)
            xs_.append(xs)
            xs_E_.append(xs_E)

            lumi_list.append(curr_lumi)
            n_event_list.append(curr_nEvt)
            # print lumi, n_event,
            # print curr_lumi, curr_nEvt

        gr = Helper.make_graphError(
            "XS", period_, [0.]*len(period_),
            xs_, xs_E_)

        canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
        gr.SetMarkerStyle(20)
        gr.Draw('A*')
        gr.GetXaxis().SetTitle("Data Period")
        gr.GetYaxis().SetTitle("# of observed / Int #it{L} [fb]")
        # add expected XS
        exp_bkg_xs = 0.675
        line = ROOT.TLine(0, exp_bkg_xs, 12, exp_bkg_xs)
        line.SetLineWidth(2)
        line.SetLineColor(4)
        line.Draw("same")
        canvas.SaveAs(out_name+"_XS.pdf")
        canvas.SaveAs(out_name+"_XS.eps")

        gr_total = Helper.make_graphError(
            "Total", lumi_list, [0.]*len(lumi_list),
            n_event_list, [math.sqrt(x) for x in n_event_list]
        )
        gr_total.SetMarkerStyle(20)
        gr_total.Draw('A*')
        gr_total.GetXaxis().SetTitle("Integrated luminosity [fb^{-1}]")
        gr_total.GetYaxis().SetTitle("# of observed in [230, 250] GeV")
        # add expected and Fit
        #fun_p1 = ROOT.TF1("fun_p1", "{:.3f}*x".format(exp_bkg_xs), 0, 36.5)
        #fun_p1.Draw("same")
        fun_p2 = ROOT.TF1("fun_p2", "[0]+[1]*x+[2]*x*x", 0, 36.5)
        gr_total.Fit("fun_p2")
        fun_p2.Draw("same")
        fun_p2.SetLineColor(4)
        canvas.SaveAs(out_name+"_nEvt.pdf")
        canvas.SaveAs(out_name+"_nEvt.eps")

        ## get significance as function of luminosity
        signal_list = [x-exp_bkg_xs*y for x,y in zip(n_event_list, lumi_list)]
        sigma_list = [Helper.significance(x, exp_bkg_xs*y) for x,y in zip(signal_list, lumi_list)]
        gr_sigma = Helper.make_graphError(
            "Total", lumi_list, [0.]*len(lumi_list),
            sigma_list, [0.]*len(sigma_list)
        )
        gr_sigma.SetMarkerStyle(20)
        gr_sigma.Draw('A*')
        gr_sigma.GetXaxis().SetTitle("Integrated luminosity [fb^{-1}]")
        gr_sigma.GetYaxis().SetTitle("significance")
        canvas.SaveAs(out_name+"_Sigma.pdf")
        canvas.SaveAs(out_name+"_Sigma.eps")

        fout = ROOT.TFile.Open(out_name+".root", "recreate")
        gr.Write()
        gr_total.Write()
        fout.Close()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print sys.argv[0]," file_name out_name"
        sys.exit(1)

    input_name = sys.argv[1]
    out_name = sys.argv[2]
    handle = EventPerRun(input_name, out_name)
    handle.loop_events()
    handle.per_period()

