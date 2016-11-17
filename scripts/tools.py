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

def add_ratio_pad(hist_list, y_title, y_min, y_max, reverse=False):
    """
    hist_list = [Data, MC1, MC2]
    plot Data/MC1, Data/MC2
    @para=reverse, plot MC1/Data
    """
    if len(hist_list) < 2:
        print "less than 2 histograms, kidding?"
        return None

    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.15, 1.0, 1.0)
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.01, 1, 0.286)
    pad1.Draw()
    pad2.Draw()

    pad2.cd()
    hist_list_cp = [x.Clone(x.GetName()+"_clone") for x in hist_list]
    h_refer = hist_list_cp[0]
    for i, hist in enumerate(hist_list_cp):
        if i==0:
            hist.Sumw2()
            hist.Divide(h_refer)

            hist.GetYaxis().SetTitle(y_title)
            hist.GetYaxis().SetRangeUser(y_min, y_max)

            hist.Draw("E3")
        else:
            # start to calculate the ratio
            if reverse: # MC/Data
                this_hist = hist.Clone(hist.GetName()+"_cp")
                this_hist.Divide(h_refer)
            else: # Data/MC
                this_hist = h_refer.Clone(h_refer.GetName()+"_cp")
                this_hist.Divide(hist)

            this_hist.Draw("HIST SAME")

    pad1.cd()
    return pad1


def compare_hists(
    hist_list, tag_list, out_name,
    x_title, y_title,
    only_shape=True, is_log=False, add_ratio=False
):
    """
    compare histograms in hist_list, named with tag_list
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

def stack_hists(
    hist_list, tag_list, out_name,
    x_title, y_title, is_log=False, has_data=True, debug=False
):
    # put data to be the first one, plz.
    colors = [1, 206, 64, 95, 28, 29, 209, 5]
    if len(hist_list) > len(colors):
        print "I don't have enough colors", len(hist_list), len(colors)
        return

    # clone current histograms, so that inputs are untouched.
    hist_list_all = []
    hist_list_cp = []
    h_data = None  # the first element is assumed to be data
    for i, hist in enumerate(hist_list):
        new_hist = hist.Clone(hist.GetName()+"_clone")
        color = colors[i]
        new_hist.SetLineColor(color)
        new_hist.SetFillColor(color)
        if i==0 and has_data:
            # decorate data points
            new_hist.SetMarkerStyle(20)
            new_hist.SetMarkerSize(1.2)
            h_data = new_hist
        else:
            hist_list_cp.append(new_hist)

    # always plot the smallest component in the bottom
    hist_sorted_list = sorted(hist_list_cp, key=lambda k:k.Integral())
    hs = ROOT.THStack("hs", "")
    for hist in hist_sorted_list:
        hs.Add(hist)

    # start to plot them
    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)

    y_max = hs.GetMaximum()
    y_min = hs.GetMinimum()
    if has_data and h_data is not None:
        if y_max < h_data.GetMaximum():
            y_max = h_data.GetMaximum()
        if y_min > h_data.GetMinimum():
            y_min = h_data.GetMinimum()

    if has_data:
        this_hist = h_data
    else:
        this_hist = hs

    if is_log:
        canvas.SetLogy()
        this_hist.GetYaxis().SetRangeUser(4E-3, y_max*1e3)
    else:
        this_hist.GetYaxis().SetRangeUser(1E-3, y_max*1.1)

    this_hist.SetNdivisions(8, "X")
    this_hist.SetXTitle(x_title)
    this_hist.SetYTitle(y_title)

    if has_data:
        h_data.Draw("EP")
        hs.Draw("HISTsame")
        h_data.Draw("AXISsame")
        h_data.Draw("EPsame")
    else:
        hs.Draw("HIST")

    # add legend
    legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.94)
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    hist_id = 0
    if has_data:
        hist_all = [h_data] + hist_list_cp
    else:
        hist_all = hist_list_cp
    for hist, tag in zip(hist_all, tag_list):
        if has_data and hist_id == 0:
            legend.AddEntry(hist, tag+" {:.0f}".format(hist.Integral()), "LP")
        else:
            legend.AddEntry(hist, tag+" {:.1f}".format(hist.Integral()), "F")
        hist_id += 1

    legend.Draw("same")
    canvas.SaveAs(out_name)
