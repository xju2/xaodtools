#!/usr/bin/env python
import commands
import os

class BsubHandle:
    def __init__(self):
        self.out_log = "log.submit"
        self.no_submit = True
        self.bad_jobs = 0
        self.good_jobs = 0

    def submit(self, run_cmd):
        bsubs_cmd = "bsub -q wisc  -R 'pool>4000' -C 0 -o " + \
                os.getcwd()+ "/"+ self.out_log+" "+ run_cmd
        if self.no_submit:
            print bsubs_cmd
            return

        status,output=commands.getstatusoutput(bsubs_cmd)
        if status != 0:
            self.bad_jobs += 1
        else:
            self.good_jobs += 1

    def print_summary(self):
        if self.no_submit:
            print "-------not submitted!!---------"
        else:
            print "Good jobs: "+ str(self.good_jobs)+", "+str(self.bad_jobs)+" failed!"

