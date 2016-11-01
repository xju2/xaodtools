 #!/usr/bin/env python
import sys
import string
import commands
import os

import glob

base_dir = os.getcwd()
exe_base = "/afs/cern.ch/user/x/xju/work/upsilon/code/MyXAODTools/bsubs/"

def check_dir(dir_):
    if not os.path.exists(dir_):
        os.mkdir(dir_)

def submit(exe, out_log_name):
    print "executable:", exe
    print "log file:", out_log_name
    bad_jobs = 0
    good_jobs = 0
    input_dir = base_dir + "/split_and_merge/"
    #input_dir = base_dir + "/split_and_merge_Jpsi/"
    input_all = glob.glob(input_dir+"x*")

    for input_name in input_all:
        run_cmd = exe + " " +input_dir+" "+os.path.basename(input_name)

        bsubs_cmd = "bsub -q wisc  -R 'pool>4000' -C 0 -o " + \
                base_dir+ "/"+ out_log_name+" "+run_cmd

        #print bsubs_cmd

        status,output=commands.getstatusoutput(bsubs_cmd)
        if status != 0:
            bad_jobs += 1
        else:
            good_jobs += 1

    print "Good jobs: "+ str(good_jobs)+", "+str(bad_jobs)+" failed!"

def submit_tree(exe, out_log_name):
    print "executable:", exe
    print "log file:", out_log_name
    bad_jobs = 0
    good_jobs = 0
    input_dir = base_dir + "/split_and_merge/"
    input_all = glob.glob(input_dir+"x*")

    check_dir(base_dir+"/histograms/")
    for input_name in input_all:
        input_basename = "merged_"+os.path.basename(input_name)+".root"
        out_name = base_dir+"/histograms/"+os.path.basename(input_basename).replace(".root", "_hist.root")
        input_name = os.path.dirname(input_name)+"/"+input_basename
        run_cmd = exe + " " +input_name+" "+out_name

        bsubs_cmd = "bsub -q wisc  -R 'pool>4000' -C 0 -o " + \
                base_dir+ "/"+ out_log_name+" "+run_cmd

        #print bsubs_cmd
        status,output=commands.getstatusoutput(bsubs_cmd)
        if status != 0:
            bad_jobs += 1
        else:
            good_jobs += 1

    print "Good jobs: "+ str(good_jobs)+", "+str(bad_jobs)+" failed!"

def submit_daod(exe, out_log_name, directory):
    print "executable:", exe
    print "log file:", out_log_name
    print "directory:", directory
    bad_jobs = 0
    good_jobs = 0

    input_dir = base_dir + "/"+directory+"/"
    input_all = glob.glob(input_dir+"x*")
    print "total files:",len(input_all)
    check_dir(base_dir+"/"+directory)
    for input_name in input_all:
        out_name = input_dir+"mini_"+os.path.basename(input_name)+".root"
        run_cmd = exe + " " +input_name+" "+out_name

        bsubs_cmd = "bsub -q wisc  -R 'pool>4000' -C 0 -o " + \
                base_dir+ "/"+ out_log_name+" "+run_cmd

        #print bsubs_cmd
        status,output=commands.getstatusoutput(bsubs_cmd)
        if status != 0:
            bad_jobs += 1
        else:
            good_jobs += 1

    print "Good jobs: "+ str(good_jobs)+", "+str(bad_jobs)+" failed!"

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print sys.argv[0]," merge_tree/make_aod/merge_aod/make_tree/read_daod log_name/directory"
        exit(1)

    option = sys.argv[1]
    out_log_name = sys.argv[2]

    if option == "merge_tree":
        exe = exe_base+"run_merge.sh"
    elif option == "merge_aod":
        exe = exe_base+"run_merge_aod.sh"
    elif option == "make_aod":
        exe = exe_base+"run_aod.sh"
    elif option == "make_tree":
        exe = exe_base+"run_draw_upsilon.sh"
        submit_tree(exe, out_log_name)
        exit(0)
    elif option == "read_daod":
        exe = exe_base+"run_daod.sh"
        if len(sys.argv) < 4:
            print "add directory"
            exit(1)

        directory= sys.argv[3]
        submit_daod(exe, out_log_name, directory)
        exit(0)
    else:
        print option," not supported!"
        exit(2)

    submit(exe, out_log_name)
