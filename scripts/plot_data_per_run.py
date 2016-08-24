#!/usr/bin/env python
"""
plot visible cross section for data as function of run-number
inputs:
    @minitree for data
    @XML that contains luminiosity per run
"""
BASE_NAME = "/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/MiniTrees/Prod_v05/data/Nominal/data_13TeV.root"

import ROOT
import sys
from array import array
from sets import Set

ROOT.gROOT.SetBatch()

if not hasattr(ROOT, "myText"):
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/AtlasUtils.C")
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c")

def load_lumi(csv_name="lumitable_MD1.csv"):
    lumi_dic = {}
    with open(csv_name, 'r') as f:
        for line in f:
            if "Run" in line: continue
            if "Total" in line: continue
            items = line[:-1].split(',')
            run = float(items[0])
            lumi = float(items[5])
            lumi_dic[run] = lumi

    return lumi_dic

def plot(file_name, out_name):
    f1 = ROOT.TFile.Open(file_name)
    tree = f1.Get("tree_incl_all")
    nentries = tree.GetEntries()
    num_events_dic = {}
    lumi_dic = load_lumi()
    for key in lumi_dic.keys():
        num_events_dic[key] = 0

    runs_not_in_GRL = []
    masses = []
    for ientry in range(nentries):
        tree.GetEntry(ientry)
        run_ = tree.run
        event_ = tree.event
        if not tree.pass_vtx4lCut: 
            print "not pass vertex"
            continue
        if run_ in lumi_dic:
            num_events_dic[run_] += 1
        elif run_ >= 297730:
            runs_not_in_GRL.append(run_)
        else:
            pass

    # runs not in GRL
    print "events should not there: ", len(runs_not_in_GRL)
    print "runs not in GRL: ",
    for run in sorted(Set(runs_not_in_GRL)):
        print run

    #save visible cross to histogram
    run_list = []
    xs_list = []
    for key in sorted(num_events_dic.keys()):
        run_list.append(key)
        xs = num_events_dic[key]/lumi_dic[key]
        xs_list.append(xs)

    gr = ROOT.TGraph(len(run_list), array('f', run_list), array('f', xs_list))
    fout = ROOT.TFile.Open(out_name, "recreate")
    gr.Write()

    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    gr.Draw("ACL*")
    canvas.SaveAs("visible_xs.pdf")
    fout.Close()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print sys.argv[0]," file_name out_name"
        sys.exit(1)

    input_name = sys.argv[1]
    out_name = sys.argv[2]
    plot(input_name, out_name)
