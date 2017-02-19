#!/usr/bin/env python
import subprocess
import os
from optparse import OptionParser
import pyAMI.client
import pyAMI.atlas.api as AtlasAPI
import re

# NEED to setup following environment
# setupATLAS
# lsetup pyami --- for get Info from AMI
# lsetup rucio --- for rucio


# links to help AMI:
# http://ami.in2p3.fr/pyAMI/command_line.html
# http://ami.in2p3.fr/pyAMI/pyAMI5_atlas_api.html

class DataSetManager:
    ATLFAST_PAT_STR = ".*e[0-9]+_a[0-9]+_a[0-9]+_r[0-9]+"
    atlfast_pat = re.compile(ATLFAST_PAT_STR)
    def __init__(self, mc_tag, daod, ptag, only_fs=True):
        """
        DAOD indicate derivation name, such as: DAOD_SUSY1, DAOD_HSG2
        ptag indicate the desired tag for the derivation
        """
        self.mc_tag = mc_tag
        self.daod = daod
        self.ptag = ptag
        self.usr_name = 'xju'
        self.client = pyAMI.client.Client('atlas')
        self.only_fs = only_fs
        print 'DataSet Manager is Created!'

    def get_pattern(self, number):
        res = self.mc_tag+"."+str(number)+".%"+self.daod+".%"+self.ptag
        return res

    def find_dataset(self, mc_id):
        out = ""
        pattern = self.get_pattern(mc_id)
        print pattern
        #output = subprocess.Popen(['rucio','list-dids', pattern],
        #                          stdout=subprocess.PIPE).communicate()
        dic_datasets = AtlasAPI.list_datasets(self.client,
                                              patterns=[pattern])

        out = None
        if len(dic_datasets) < 1:
            print pattern," not found"
            return None
        else:
            if len(dic_datasets) > 2:
                print pattern," more than 3 found"
            for dataset in dic_datasets:
                ds_name = dataset['ldn']
                #print re.match(DataSetManager.atlfast_pat, ds_name)
                if re.match(DataSetManager.atlfast_pat, ds_name) and self.only_fs:
                    print "Skip fast simulation"
                    continue
                else:
                    out = ds_name
                    break
        return out

    def create_container(self, ctn_name, input_ds_list):
        """
        ctn_name, container name
        input_ds_list, a text file containing the datasets
        """
        final_name = 'user.'+self.usr_name+'.'+self.mc_tag+'.'+ctn_name+'.'+self.daod+'.'+self.ptag+'/'
        print "Container: ",final_name
        arg = ['rucio', 'add-container', final_name]
        arg.append(final_name)
        q = subprocess.Popen(arg,subprocess.PIPE).communicate()
        print final_name,"is created"

        # add dataset list
        arg_add_list = ['rucio', 'attach', final_name, '-f', input_ds_list]
        q = subprocess.Popen(arg,subprocess.PIPE).communicate()
        print q

    def info(self, dataset_name):
        AtlasAPI.init()
        files_info = AtlasAPI.get_dataset_info(self.client, dataset_name)[0]
        #print files_info
        total_N = int(files_info["totalEvents"])
        try:
            xs = float(files_info['crossSection'])
        except:
            xs = -999.0
        try:
            filter_eff = float(files_info['GenFiltEff_mean'])
        except:
            print "no GenFiltEff_mean, set to 1"
            filter_eff = 1.0
        eff_xs = xs* filter_eff
        return (total_N, xs, filter_eff, eff_xs)

    def get_info(self, dsid):
        ds = self.find_dataset(dsid)
        if ds is not None:
            print ds
            return (ds,self.info(ds))
        else:
            return (None,None)


if __name__ == "__main__":
    # example to run...
    # python rucio_mgr.py -r 341274,341274 --daod EVNT -p e3940
    usage = "%prog [options]"
    version="%prog 1.0"
    parser = OptionParser(usage=usage, description="get XS info for DSID", version=version)
    parser.add_option("-r", "--range", dest='range', default=None, help="set range of runs: 310015,311481")
    parser.add_option("-l", "--list", dest='list', default=None, help="give a list of runs  310969,311170")
    parser.add_option("--mc", dest='mc', default="mc15_13TeV", help="MC production version")
    parser.add_option("--daod", dest='daod', default="DAOD_HIGG2D1", help="Derivation name")
    parser.add_option("-p", "--ptag", dest='ptag', default="p2666", help="tag for derivation")
    (options,args) = parser.parse_args()

    if options.range is None and options.list is None:
        print parser.print_help()
        exit(1)

    if options.range:
        first, end = options.range.split(',')
        dsid_list = range(int(first), int(end)+1)

    if options.list:
        dsid_list = [int(x) for x in options.list.split(',')]

    ds_mgr = DataSetManager(options.mc, options.daod, options.ptag)
    out_text = "DSID & Name & Total Events & Cross Section [$nb$] & Efficiency of Generator Filter & Effective Cross section [$nb$] \\\\ \\hline \n"
    for dsid in dsid_list:
        print dsid
        ds_name,info_ = ds_mgr.get_info(dsid)
        event_type = ds_name.split('.')[2].replace('_', '\_')
        if info_ is not None:
            out_text += str(dsid) + " & "+ event_type +" & "\
                    " & ".join([str(x) for x in info_])+" \\\\ \n"
    print out_text
