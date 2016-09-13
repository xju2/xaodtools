#!/usr/bin/env python

import ROOT
from array import array
ROOT.gROOT.SetBatch()

import AtlasStyle

if not hasattr(ROOT, "myText"):
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/AtlasUtils.C")
    ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c")

import compare_mine_TS as cpt

def get_dic():
    file_name = "note_xevt_up_to_50GeV_PVassoc_bak.txt"
    can_michigan = {}
    with open(file_name, 'r') as f:
        icount = 0
        for line in f:
            icount += 1
            if icount == 1:
                continue
            items = line[:-1].split('*')
            #print items
            run = int(items[2])
            event = int(items[3])
            mass = float(items[4])
            chi2 = float(items[5])
            can_michigan[(run,event)] = (mass, chi2)

    return can_michigan
    #print can_michigan

def process():
    can_michigan = get_dic()
    #print can_michigan
    f1 = ROOT.TFile.Open("all.root")
    tree = f1.Get("bls")
    nentries = tree.GetEntries()

    fout = ROOT.TFile.Open("diff.root", "recreate")
    out_tree = ROOT.TTree("physics", "")

    b_event = array('i', [0])
    b_run = array('i', [0])
    xy_chi2 = array('f', [0])
    xy_m4l = array('f', [0])
    ts_chi2 = array('f', [0])
    ts_m4l = array('f', [0])

    out_tree.Branch("xy_chi2", xy_chi2, "xy_chi2/F")
    out_tree.Branch("xy_mass", xy_m4l,  "xy_mass/F")
    out_tree.Branch("ts_chi2", ts_chi2, "ts_chi2/F")
    out_tree.Branch("ts_mass", ts_m4l,  "ts_mass/F")
    out_tree.Branch("event", b_event, "event/I")
    out_tree.Branch("run", b_run, "run/I")

    chi2_sig_xy = ROOT.TH1F("chi2_sig_xy", "chi2 signal;quad #chi^{2};Events/2", 50, 0, 100)
    chi2_sig_ts = chi2_sig_xy.Clone("chi2_sig_ts")
    m4l_xy = ROOT.TH1F("m4l_xy", "m4l;m_{4l} [GeV]; Events/0.2GeV", 60, 14, 26)
    m4l_ts = m4l_xy.Clone("m4l_ts")

    m4l_full_xy = ROOT.TH1F("m4l_full_xy", "m4l;m_{4l} [GeV]; Events/0.4GeV", 100, 10, 50)
    m4l_full_ts = m4l_full_xy.Clone("m4l_full_ts")

    for ientry in xrange(nentries):
        tree.GetEntry(ientry)
        run = int(tree.run)
        event = int(tree.event)
        #print run,event,tree.run,tree.event
        #print run,event
        if (run, event) not in can_michigan:
            print run,event,"not in michigan"
            continue

        michigan_mass,michigan_chi2 = can_michigan[(run, event)]
        xy_m4l[0] = tree.m4l_fitted
        xy_chi2[0] = tree.x_chi2*5.0

        ts_chi2[0] = michigan_chi2
        ts_m4l[0] = michigan_mass
        b_event[0] = event
        b_run [0] = run
        out_tree.Fill()

        chi2_sig_xy.Fill(tree.x_chi2*5.0)
        chi2_sig_ts.Fill(michigan_chi2)
        m4l_xy.Fill(tree.m4l_fitted)
        m4l_ts.Fill(michigan_mass)
        m4l_full_xy.Fill(tree.m4l_fitted)
        m4l_full_ts.Fill(michigan_mass)


    fout.cd()
    out_tree.Write()
    chi2_sig_ts.Write()
    chi2_sig_xy.Write()
    m4l_xy.Write()
    m4l_ts.Write()
    m4l_full_xy.Write()
    m4l_full_ts.Write()

    fout.Close()
    f1.Close()

def compare_hist():
    f1 = ROOT.TFile.Open("diff.root")

    cpt.save_compare(cpt.norm_hist(f1.Get("chi2_sig_xy")), cpt.norm_hist(f1.Get("chi2_sig_ts")), "chi2_sameevents")
    cpt.save_compare(cpt.norm_hist(f1.Get("m4l_xy")), cpt.norm_hist(f1.Get("m4l_ts")), "mass_sameevents")
    cpt.save_compare(cpt.norm_hist(f1.Get("m4l_full_xy")), cpt.norm_hist(f1.Get("m4l_full_ts")), "mass_fullRange_sameevents")

    f1.Close()

if __name__ == "__main__":
    process()
    compare_hist()
