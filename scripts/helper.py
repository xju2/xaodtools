#!/usr/bin/env python

import ROOT
import math
from array import array

def loader(infile, tree_name):
    chain = ROOT.TChain(tree_name)
    if "root" in infile:
        chain.Add(infile)
        return chain

    # infile is a list of input root files
    ncounter = 0
    with open(infile, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            chain.Add(line[:-1])
            ncounter += 1

    print "total events:",chain.GetEntries(),"in",ncounter,"files"
    return chain

def get_eff(a, b):
    eff = a*1.0/b
    error = eff * math.sqrt(1./a + 1./b)
    print "{:.4f} {:.4f}".format(eff,error)
    return eff,error

def get_eff_error(a, ae, b, be):
    eff = a*1.0/b
    error = eff*math.sqrt((ae/a)**2 + (be/b)**2)
    return eff,error

def make_hist(hist_name, bin_list):
    """
    create TH1F using a list as x-axis
    """
    nbins = len(bin_list) -1
    h1 = ROOT.TH1F(hist_name, hist_name, nbins, array('f', bin_list))
    return h1

def make_graph(name_, x_, y_):
    """
    Graph only
    """
    gr = ROOT.TGraph(len(x_), array('f', x_), array('f', y_))
    gr.SetName(name_)
    return gr

def make_graphError(name_, x_, xe_, y_, ye_):
    """
    Symmetric errors
    """
    gr = ROOT.TGraphErrors(
        len(x_), array('f', x_), array('f', y_),
        array('f', xe_), array('f', ye_)
    )
    gr.SetName(name_)
    return gr

def column(matrix, i):
    return [row[i] for row in matrix]

def flatten(matrix):
    return [x for y in matrix for x in y]

def significance(s, b):
    if s < 0 or b < 0:
        return 0
    if abs(b) < 1E-6:
        return 0
    return math.sqrt(2*((s+b)*math.log(1+s/b) - s))
