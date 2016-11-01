#!/usr/bin/env python
import glob
import os
import ROOT
from optparse import OptionParser

usage = "%prog [option] "
parser = OptionParser(description="check if any files missing", usage=usage)
parser.add_option("-a", "--aod", dest='daod', default=False, help="check output of DAOD", action="store_true")
#parser.add_option("-m", "--mini", dest='mini', default=False, help="check output of mini-tree", action="store_true")

options,args = parser.parse_args()

has_bad = False
tag="../split_and_merge/x*"
for file_ in glob.glob(tag):
    if options.daod:
        file_name = "mini_"+os.path.basename(file_)+".root"
        check_size = False
    else:
        file_name = "merged_"+os.path.basename(file_)+"_hist.root"
        check_size = True
    if not os.path.exists(file_name):
        has_bad = True

        if check_size:
            # check if the input is empty
            tag_name = "data16_v11_onlyTrack_addQuality"
            base_name="root://eosatlas//eos/atlas/unpledged/group-wisc/users/xju/bphys/merged/"+tag_name+"/split_and_merge/"
            root_name = base_name + "merged_"+os.path.basename(file_)+".root"
            f1 = ROOT.TFile.Open(root_name)
            tree = f1.Get("physics")
            print file_name, tree.GetEntries()
            f1.Close()
        else:
            print file_name

if not has_bad:
    print "all there"
