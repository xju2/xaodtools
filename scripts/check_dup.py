#!/usr/bin/env python

import ROOT
import sys
from sets import Set

def check(file_name, tree_name, run_name, event_name):
    fin = ROOT.TFile.Open(file_name)
    tree = fin.Get(tree_name)
    nentries = tree.GetEntries()

    seen = Set([])
    for ientry in range(nentries):
        tree.GetEntry(ientry)
        run = getattr(tree, run_name)
        event = getattr(tree, event_name)
        if (run,event) in seen:
            return True
        seen.add( (run,event) )

    return False

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print sys.argv[0],"file_name tree_name run_name event_name"
        exit(0)

    if check(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]):
        print "has duplicated events!"
    else:
        print "no duplication"
    
