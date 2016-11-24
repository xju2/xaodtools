#!/usr/bin/env python
from optparse import OptionParser

total_summary = ""
def load_lumi(csv_name):
    lumi_dic = {}
    with open(csv_name, 'r') as f:
        for line in f:
            if "Run" in line: continue
            if "Total" in line:
                total_summary = line
                continue
            items = line[:-1].split(',')
            run = int(items[0])
            lumi = float(items[5])
            lumi_dic[run] = lumi
    return lumi_dic

usage = "%prog [option] table"
parser = OptionParser(description="calculate lumonisity", usage=usage)
parser.add_option("-r", "--range", dest='range', default="None", help="set range of runs: 310015,311481")
parser.add_option("-l", "--list", dest='list', default="None", help="give a list of runs  310969,311170")
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
        total_lumi += lumi_dic[run]

else:
    total_runs = len(lumi_dic)
    total_lumi = total_summary
    #total_lumi = float(total_summary.split(',')[5])

print "total runs:", total_runs
print "total lumi:", total_lumi
