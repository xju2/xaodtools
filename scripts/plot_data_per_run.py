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
        lumi_table_2016 = "/afs/cern.ch/user/x/xju/work/myRooCoreTools/MyXAODTools/scripts/lumitable_data16_297730_311481_v88_pro20-21.csv"
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
            ## only look at 4e candidates in [230, 250] GeV
            if not event_type == 1:
                continue
            if not (m4l > 240 and m4l < 245):
                continue
        
            ntotal += 1
            if run_ in self.cal.lumi_dic:
                self.num_events_dic[run_] += 1
            else:
                runs_not_in_GRL.append(run_)

        print "Total number of events in 4e",ntotal


    def events(self, start, end):
        """
        total number of events for run in range of (start, end)
        """
        total = 0
        for key, value in self.num_events_dic.iteritems():
            if key >= start and key <= end:
                total += value
        return total

    def per_period(self):
        period_ = [x+1 for x in range(len(PERIOD2016))]
        xs_ = []
        xs_E_ = []
        for start,end in PERIOD2016:
            #print start,end
            lumi,total_runs = self.cal.get_lumi_for_range(start, end)
            if lumi == 0:
                continue
            #print lumi,total_runs
            n_event = self.events(start, end)
            xs, xs_E = Helper.get_eff_error(
                n_event, math.sqrt(n_event),
                lumi, 0.032*lumi)
            xs_.append(xs)
            xs_E_.append(xs_E)

        gr = Helper.make_graphError(
            "XS", period_, [0.]*len(period_),
            xs_, xs_E_)

        canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
        gr.SetMarkerStyle(20)
        gr.Draw('A*')
        canvas.SaveAs("TEST_XS.pdf")


def junk():
    #save visible cross to histogram
    run_list = []
    xs_list = []
    evt_list = []

    # number of events as funcition of integrated luminosity
    lumi_list = []
    n_event_list = []
    curr_lumi = 0
    curr_nEvt = 0
    print "Total runs",len(num_events_dic.keys())
    for key in sorted(num_events_dic.keys()):
        lumi = lumi_dic[key]
        nEvt = num_events_dic[key]
        xs = nEvt/lumi
        if xs > 0.02:
            print key,round(xs,3),round(lumi,2),nEvt

        curr_lumi += lumi/1E3
        curr_nEvt += nEvt
        if nEvt < 1:
            continue

        run_list.append(key)
        xs_list.append(xs)
        evt_list.append(nEvt)
        lumi_list.append(curr_lumi)
        n_event_list.append(curr_nEvt)

    gr = ROOT.TGraph(len(run_list), array('f', run_list), array('f', xs_list))
    fout = ROOT.TFile.Open(out_name+".root", "recreate")
    gr.SetName("visualXS")
    gr.Write()

    #gr.SetMarkerStyle(21)
    gr.Draw("A*")
    canvas.SaveAs(out_name+".pdf")

    gr_lumi = ROOT.TGraph(len(lumi_list), array('f', lumi_list), array('f', n_event_list))
    gr_lumi.SetName("totalEvts")
    gr_lumi.Write()
    #gr_lumi.SetMarkerStyle(21)
    gr_lumi.Draw("A*")
    canvas.SaveAs(out_name+"_lumi.pdf")
    ##

    gr_evt = ROOT.TGraph(len(run_list), array('f', run_list), array('f', evt_list))
    gr_evt.SetName("nEvts")
    gr_evt.Write()
    #gr_evt.SetMarkerStyle(21)
    gr_evt.Draw("A*")
    canvas.SaveAs(out_name+"_nEvt.pdf")
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

