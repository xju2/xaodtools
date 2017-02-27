#!/usr/bin/env python

import ROOT

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
