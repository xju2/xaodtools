#!/usr/bin/env python
import glob
import os
import ROOT
from optparse import OptionParser
import commands

usage = "%prog [option] "
parser = OptionParser(description="check if any files missing", usage=usage)
parser.add_option("-a", "--aod", dest='daod', default=False, help="check output of DAOD", action="store_true")
parser.add_option("--eos", dest='eos', default=False, help="put reprocessed DAOD to EOS", action="store_true")
parser.add_option("--noProcess", dest='process', default=True, help="reprocess DAOD?", action="store_false")
#parser.add_option("-m", "--mini", dest='mini', default=False, help="check output of mini-tree", action="store_true")

options,args = parser.parse_args()

pwd = os.getcwd()
#tag_name = "data16_v11_onlyTrack_addQuality"
tag_name = pwd.split("/")[8]
print tag_name


has_bad = False
tag="../split_and_merge/x*"
for file_ in glob.glob(tag):
    split_name = os.path.basename(file_)
    if options.daod:
        file_name = "mini_"+split_name+".root"
        check_size = False
    else:
        file_name = "merged_"+split_name+"_hist.root"
        check_size = True
    if not os.path.exists(file_name):
        has_bad = True

        if check_size:
            # check if the input is empty
            base_name="root://eosatlas//eos/atlas/unpledged/group-wisc/users/xju/bphys/merged/"+tag_name+"/split_and_merge/"
            root_name = base_name + "merged_"+os.path.basename(file_)+".root"
            f1 = ROOT.TFile.Open(root_name)
            tree = f1.Get("physics")
            print file_name, tree.GetEntries()
            f1.Close()
        else:
            print file_name, file_
            # re-process these runs.
            if options.daod and options.process:
                cmd = "/afs/cern.ch/user/x/xju/work/upsilon/run/run_analysis.sh "+file_+" "+file_name
                status = 0
                status,output=commands.getstatusoutput(cmd)
                #print cmd
                if status != 0:
                    print cmd,"failed"
                if options.eos:
                    current_tag = os.getcwd().split('/')[-2]
                    print current_tag
                    eos_dir = "root://eosatlas//eos/atlas/unpledged/group-wisc/users/xju/bphys/merged/"+current_tag+"/split_and_merge/"
                    status = 0
                    cmd_eos = "xrdcp "+file_name+" "+eos_dir+"merged_"+split_name+".root"
                    #print cmd_eos
                    status,output=commands.getstatusoutput(cmd_eos)
                    if status != 0:
                        print cmd_eos,"failed"


if not has_bad:
    print "all there"
