 #!/usr/bin/env python
import sys
import string
import commands
import os
from optparse import OptionParser

import glob

base_dir = os.getcwd()
exe_base = "/afs/cern.ch/user/x/xju/work/upsilon/code/MyXAODTools/bsubs/monojet/"

def check_dir(dir_):
    if not os.path.exists(dir_):
        os.mkdir(dir_)

def submit(exe, out_log_name):
    print "executable:", exe
    print "log file:", out_log_name
    bad_jobs = 0
    good_jobs = 0
    input_dir = base_dir + "/split_and_merge/"
    input_all = glob.glob(input_dir+"x*")

    check_dir(base_dir+"/histograms/")
    for input_name in input_all:

        out_name = base_dir+"/histograms/merged_"+os.path.basename(input_name)+"_hist.root"
        run_cmd = exe + " " +input_name+" "+out_name

        bsubs_cmd = "bsub -q wisc  -R 'pool>4000' -C 0 -o " + \
                base_dir+ "/"+ out_log_name+" "+run_cmd

        #print bsubs_cmd
        status,output=commands.getstatusoutput(bsubs_cmd)
        status = 0
        if status != 0:
            bad_jobs += 1
        else:
            good_jobs += 1

    print "Good jobs: "+ str(good_jobs)+", "+str(bad_jobs)+" failed!"

if __name__ == "__main__":

    usage = "%prog log_name"
    parser = OptionParser(description="submit jobs for monojet", usage=usage)
    parser.add_option("--new_file", dest="new_file", default=False, action="store_true", help="create new file")
    (options,args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        exit(1)

    out_log_name = args[0]

    if options.new_file:
        exe = exe_base+"run_jetsmearing.sh"
    else:
        parser.print_help()
        exit(2)

    submit(exe, out_log_name)
