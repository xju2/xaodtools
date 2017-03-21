#!/usr/bin/env python
from optparse import OptionParser
from sets import Set

total_summary = ""

class LumiCal:
    def __init__(self):
        self.lumi_dic = {}

    def load_lumi(self, csv_name):
        print "loading luminosity from", csv_name
        with open(csv_name, 'r') as f:
            for line in f:
                if "Run" in line:
                    continue
                if "Total" in line:
                    continue
                items = line[:-1].split(',')
                try:
                    run = int(items[0])
                except:
                    print items[0]
                    exit(1)

                lumi = float(items[6])
                self.lumi_dic[run] = lumi

        print "Runs:", len(self.lumi_dic)

    def check_file(self, file_name):
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

    def get_lumi_for_runs(self, run_list):
        total_lumi = 0
        for run in run_list:
            run = int(run)
            try:
                lumi = self.lumi_dic[run]
            except KeyError:
                print run,"no in GRL"
                continue
            total_lumi += lumi
        return total_lumi

    def get_lumi_for_range(self, start, end):
        total_runs = 0
        total_lumi = 0
        for key,value in self.lumi_dic.iteritems():
            if key >= start and key <= end:
                total_lumi += value
                total_runs += 1
        return (total_lumi, total_runs)

    def total_lumi(self):
        total_runs = len(self.lumi_dic)
        total_lumi = sum(self.lumi_dic.values())
        return (total_lumi, total_runs)


if __name__ == "__main__":
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
    cal = LumiCal()
    cal.load_lumi(table)

    total_lumi = 0
    total_runs = 0

    if options.range != "None":
        first, end = options.range.split(',')
        first_run = int(first)
        end_run = int(end)
        total_lumi,total_runs = cal.get_lumi_for_range(first_run, end_run)

    elif options.list != "None":
        run_list = options.list.split(',')
        total_lumi, total_runs = cal.get_lumi_for_runs(run_list)

    else:
        total_lumi,total_runs =cal.total_lumi()

    print "total runs:", total_runs
    print "total lumi:", total_lumi

    if options.check != "None":
        file_name = options.check
        check_file(lumi_dic, file_name)

