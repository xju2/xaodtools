#!/usr/bin/env python
import ROOT
import sys

def load_lumi(csv_name):
    lumi_dic = {}
    with open(csv_name, 'r') as f:
        for line in f:
            if "Run" in line: continue
            if "Total" in line: continue
            items = line[:-1].split(',')
            run = float(items[0])
            lumi = float(items[5])
            lumi_dic[run] = lumi
    return lumi_dic

def get_lumi(table, first_run, last_run): 
    lumi_dic = load_lumi(table)
    total_lumi = 0
    total_runs = 0
    for key,value in lumi_dic.iteritems():
        if key >= first_run and key <= last_run:
            total_lumi += value
            total_runs += 1

    print "total runs:", total_runs
    print "total lumi:", total_lumi

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print sys.argv[0]," table first end"
        sys.exit(1)

    table = sys.argv[1]
    first = float(sys.argv[2])
    last = float(sys.argv[3])
    get_lumi(table, first, last)
