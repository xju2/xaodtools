#!/usr/bin/env python
"""
obtain global pvalue according to
- number of upper crosssing, w.r.t 0.7sigma
- local pvalue
ref: http://arxiv.org/abs/1005.1891
"""
import ROOT
import math
import sys

def get_global_sigma(input_, upper_xs, is_pvalue):
    if is_pvalue:
        local_p = input_
        local_sig = ROOT.RooStats.PValueToSignificance(input_)
    else:
        local_p = ROOT.RooStats.SignificanceToPValue(input_)
        local_sig = input_

    trial_factor = 1 + math.sqrt(math.pi) / 2. * upper_xs * local_sig
    global_pvalue = local_p * trial_factor
    global_significance = ROOT.RooStats.PValueToSignificance(global_pvalue)
    return global_significance

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print sys.argv[0]," local_p upper_xs is_pvalue"
        sys.exit(1)

    local_p = float(sys.argv[1])
    upper_xs = float(sys.argv[2])
    is_pvalue = bool(int(sys.argv[3]))
    global_sigma = get_global_sigma(local_p, upper_xs, is_pvalue)
    if is_pvalue:
        print "local sigma: {:.2f}".format(ROOT.RooStats.PValueToSignificance(local_p))
    else:
        print "local sigma: {:.2f}".format(local_p)

    print "global sigma: {:.2f}".format(global_sigma)

