#!/usr/bin/env python
__author__ = "Xiangyang Ju"
__version__ = "0.1"
__doc__ = """
auxillary tools
"""

import math
import ROOT
from array import array
ROOT.gROOT.SetBatch()

def get_ratio_and_error(a, b):
    f = a/b
    error = f*math.sqrt((a+b)/(a*b))
    print f,error
    return (f, error)

def make_hist(hist_name, bin_list):
    """
    create TH1F using a list as x-axis
    """
    nbins = len(bin_list) -1
    h1 = ROOT.TH1F(hist_name, hist_name, nbins, array('f', bin_list))
    return h1

def add_ratio_pad(hist_list):
    """
    treat first histogram as denominator, plot the ratio of other histograms over the first histogram.
    """

def compare_hists(
    hist_list, tag_list, out_name,
    x_title, y_title,
    only_shape=True, is_log=False, add_ratio=False
):
    """
    compare histograms
    """
    colors = [1, 2, 4, 8, 6]
    if len(hist_list) != len(tag_list):
        print "size of histograms' list different from size of tags"
        return False

    if len(hist_list) > len(colors):
        print "I don't enough colors.."
        return False

    # clone current histograms, so that inputs are untouched.
    hist_list_clone = [x.Clone(x.GetName()+"_clone") for x in hist_list]

    hist_id = 0
    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    for hist, tag in zip(hist_list_clone, tag_list):
        color = colors[hist_id]
        if only_shape:
            hist.Scale(1./hist.Integral())

        # set general properties
        hist.SetLineColor(color)
        if hist_id == 0:
            # set axis title, range and etc..
            hist.GetXaxis().SetTitle(x_title)
            hist.GetYaxis().SetTitle(y_title)
            hist.SetMarkerColor(color)
            hist.SetMarkerSize(0.5)
            hist.Draw("EP")
            legend.AddEntry(hist, tag, "LP")
        else:
            hist.Draw("same")
            legend.AddEntry(hist, tag, "L")

        hist_id += 1

    legend.Draw("same")
    canvas.SaveAs(out_name)
