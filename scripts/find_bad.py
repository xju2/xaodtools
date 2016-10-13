#!/usr/bin/env python
import glob
import os

has_bad = False
for file_ in glob.glob("../split_and_merge/x*"):
    #file_name = "merged_"+os.path.basename(file_)+"_hist.root"
    file_name = "mini_"+os.path.basename(file_)+".root"
    if not os.path.exists(file_name):
        print file_name
        has_bad = True

if not has_bad:
    print "all there"
