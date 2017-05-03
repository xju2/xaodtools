#!/usr/bin/env python
__doc__="""
compare the histograms made from read_monojet_minitree.py
"""
from ploter import Ploter
import ROOT
ROOT.gROOT.SetBatch()

from optparse import OptionParser

def process( dir_name="histograms", scale_qcd=True, is_debug=False):
    import AtlasStyle
    SAMPLES=[
        # defines name of sample
        #(tag_name, file_name)
        ("data", "hist_data.root"),
        ("Z(#rightarrow #nu#nu) + jets", "hist_Znunu.root"),
        ("W(#rightarrow l#nu) + jets", "hist_wlnu.root"),
        ("Z(#rightarrow ll) + jets", "hist_zll.root"),
        ("Dibosons", "hist_VV.root"),
        ("#bar{t}t + single t", "hist_ttbar.root"),
        ("#gamma + jets", "hist_gamma.root"),
        # !!!!put the QCD component in the end!!!!
        ("Multi-jets", "hist_qcd.root"),
    ]
    file_list = []
    tag_list =[]
    for tag, sample_name in SAMPLES:
        file_list.append( ROOT.TFile.Open(dir_name+"/"+sample_name) )
        tag_list.append(tag)

    hist_names = [
        # (hist_name, is_log)
        ("met", True),
        ("jetPt", True),
        ("njets", False),
        ("minDphi", False)
    ]
    has_data = True
    n_qcd = 0
    total_lumi = 36.4
    for hist_name, is_log in hist_names:
        hist_list = []
        scale_by_width = False
        if "met" in hist_name or "jetPt" in hist_name:
            scale_by_width = True

        for ifile, file_ in enumerate(file_list):
            hist = file_.Get("h_"+hist_name)
            if has_data and ifile == 0:
                scale_factor = 1.
            else:
                scale_factor = total_lumi

            if scale_by_width:
                hist.Scale(scale_factor, "width")
            else:
                hist.Scale(scale_factor)

            hist.SetName(hist.GetName()+"_"+file_.GetName())

            # check if QCD background
            if "qcd" in file_.GetName() and scale_qcd and has_data:
                n_data = hist_list[0].Integral()
                n_bkg = 0
                for i,hist_tmp in enumerate(hist_list):
                    if i==0:
                        continue
                    n_yield = hist_tmp.Integral()
                    if is_debug:
                        print hist_tmp.GetName(),n_yield
                    n_bkg += n_yield
                if is_debug:
                    print "data: ",n_data,"bkg:",n_bkg
                # scale QCD background
                n_qcd = n_data-n_bkg
                hist.Scale( (n_qcd)/hist.Integral() )

            hist_list.append(hist)

        out_dir = "plots/"
        out_name = out_dir+hist_name+"_log"+str(is_log)+"_scaled"+str(scale_qcd)+".pdf"

        has_ratio = True
        plot_helper = Ploter()
        plot_helper.stack_hists(
            hist_list, tag_list, out_name,
            hist.GetXaxis().GetTitle(),
            hist.GetYaxis().GetTitle(),
            is_log, has_data)

    print "QCD:",n_qcd
    # close all the files
    for file_ in file_list:
        file_.Close()

if __name__ == "__main__":
    usage = "%prog dir_name"
    parser = OptionParser(usage=usage, description="make plots")
    parser.add_option("-v", "--verbose", dest="verbose", default=False, action="store_true", help="debug mode")
    parser.add_option("--noScale", dest="scale", default=True, action="store_false", help="scale QCD to data")
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        exit(0)
    else:
        dir_name = args[0]

    process(dir_name, options.scale, options.verbose)
