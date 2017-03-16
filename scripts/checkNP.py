#!/usr/bin/env python

import ROOT
import os
from optparse import OptionParser

def run(config, path, options):
    ROOT.gROOT.LoadMacro("/Users/xju/work/myRooCoreTools/MyXAODTools/scripts/checkNP.C")
    print config
    print path
    pathTS = ROOT.TString(path)
    ROOT.checkNP(config, pathTS, options.nNPs)


if __name__ == "__main__":
    usage = "%prog [option] configFile pathToFiles"
    parser = OptionParser(description="plot NPs from a text file", usage=usage)
    parser.add_option('-n', help="number of NPs", default=27, type='int', dest='nNPs')
    options,args = parser.parse_args()

    if len(args) < 1:
        print parser.print_help()
        exit(1)

    configFile = args[0]
    if len(args) > 1:
        pathToFiles = args[1]
    else:
        pathToFiles = os.path.dirname(os.path.abspath(configFile))
    run(configFile, pathToFiles, options)
