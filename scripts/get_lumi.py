#!/usr/bin/env python
from optparse import OptionParser
from sets import Set

total_summary = ""
def load_lumi(csv_name):
    print "loading luminosity from", csv_name
    lumi_dic = {}
    with open(csv_name, 'r') as f:
        for line in f:
            if "Run" in line: continue
            if "Total" in line:
                total_summary = line
                continue
            items = line[:-1].split(',')
            run = int(items[0])
            lumi = float(items[6])
            lumi_dic[run] = lumi
    return lumi_dic

def check_file(lumi_dic, file_name):
    available_runs = []
    with open(file_name, 'r') as f:
        for lines in f:
            this_run = int(lines.split('.')[1])
            available_runs.append(this_run)

    available_runs_set = Set(available_runs)
    runs_in_GRL = Set(lumi_dic.keys())
    missing_runs = runs_in_GRL.difference(available_runs_set)
    more_runs = available_runs_set.difference(runs_in_GRL)

    if len(missing_runs) == 0:
        print "all runs are there!"
    else:
        print "runs in GRL not in the file:"
        for run in sorted(missing_runs):
            print run,
    print
    if len(more_runs) == 0:
        print "no extra runs"
    else:
        print "Extra runs not in GRL:"
        for run in sorted(more_runs):
            print run,

usage = "%prog [option] table"
parser = OptionParser(description="calculate lumonisity", usage=usage)
parser.add_option("-r", "--range", dest='range', default="None", help="set range of runs: 310015,311481")
parser.add_option("-l", "--list", dest='list', default="None", help="give a list of runs  310969,311170")
parser.add_option("--check", dest='check', default="None", help="check if all the good runs available in a container list")
options,args = parser.parse_args()

if len(args) < 1:
    print parser.print_help()
    exit(1)

table = args[0]
lumi_dic = load_lumi(table)
total_lumi = 0
total_runs = 0

if options.range != "None":
    first, end = options.range.split(',')
    first_run = int(first)
    end_run = int(end)
    for key,value in lumi_dic.iteritems():
        if key >= first_run and key <= end_run:
            total_lumi += value
            total_runs += 1

elif options.list != "None":
    run_list = options.list.split(',')
    total_lumi = 0
    total_runs = len(run_list)
    for run in run_list:
        run = int(run)
        try:
            lumi = lumi_dic[run]
        except KeyError:
            print run,"does not in GRL"
            continue

        total_lumi += lumi

else:
    total_runs = len(lumi_dic)
    total_lumi = total_summary
    #total_lumi = float(total_summary.split(',')[5])

print "total runs:", total_runs
print "total lumi:", total_lumi

if options.check != "None":
    file_name = options.check
    check_file(lumi_dic, file_name)

