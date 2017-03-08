#!/usr/bin/env python
"""
get the yields for H4l analysis,
The samples, input minitress are encoded!
Change the function if it does fit you.
"""

import ROOT
from optparse import OptionParser
from collections import OrderedDict
import math
import re

ROOT.gROOT.SetBatch()

class MinitreeReader():
    def __init__(self, options):
        self.TREE_NAME = "tree_incl_all"
        self.weight_name = "weight"
        self.options = options

    def get_cuts(self):
        current_ana = self.options.analysis
        m4l_ = self.options.poi
        dic = OrderedDict()
        if current_ana == "HighMass":

            if self.options.no_VBF:
                dic["2mu2e"] = "pass_vtx4lCut==1 && 130<"+m4l_+"&&"+m4l_+"<1500&&(event_type==3||event_type==2)"
                dic["4e"] = "pass_vtx4lCut==1 && 130<"+m4l_+"&&"+m4l_+"<1500&&(event_type==1)"
                dic["4mu"] = "pass_vtx4lCut==1 && 130<"+m4l_+"&&"+m4l_+"<1500&&(event_type==0)"
            else:
                dic["2mu2e"] = "pass_vtx4lCut==1 && 130<"+m4l_+"&&"+m4l_+"<1500&&(event_type==3||event_type==2) && prod_type_HM==0"
                dic["4e"] = "pass_vtx4lCut==1 && 130<"+m4l_+"&&"+m4l_+"<1500&&(event_type==1) && prod_type_HM==0"
                dic["4mu"] = "pass_vtx4lCut==1 && 130<"+m4l_+"&&"+m4l_+"<1500&&(event_type==0) && prod_type_HM==0"
                dic["VBF"] = "pass_vtx4lCut==1 && 130<"+m4l_+"&&"+m4l_+"<1500&&prod_type_HM==1"

        elif current_ana == "LowMass":
            dic = {
                "ggF_2e2mu_13TeV":"pass_vtx4lCut==1 && 115<m4l_constrained && m4l_constrained<130 && event_type==3",
                "ggF_2mu2e_13TeV":"pass_vtx4lCut==1 && 115<m4l_constrained && m4l_constrained<130 && event_type==2",
                "ggF_4e_13TeV"   :"pass_vtx4lCut==1 && 115<m4l_constrained && m4l_constrained<130 && event_type==1",
                "ggF_4mu_13TeV"  :"pass_vtx4lCut==1 && 115<m4l_constrained && m4l_constrained<130 && event_type==0",
            }
        else:
            pass
        return dic

    def get_reducible(self):
        """return a dictionary for pre-defined background, the keys should match ones in get_cuts"""
        current_ana = self.options.analysis
        dic = OrderedDict()
        norm = stats = sys = 1
        if current_ana == "HighMass":
            if self.options.no_VBF:
                dic["2mu2e"] = (7.95706, 0.2926057, 1.1205216)
                dic["4e"] = (4.48462, 0.241686, 0.789469)
                dic["4mu"] = (3.78732, 0.159662, 0.884825)
            else:
                dic["2mu2e"] = (7.77571, 0, 1.05572)
                dic["4e"] = (4.37282, 0, 0.759606)
                dic["4mu"] = (3.71372, 0, 0.808304)
                dic["VBF"] = (0.3724366, 0, 0.049943)
        else:
            print "I don't know"

        return dic


    def get_samples(self):
        analysis = self.options.analysis
        mc_dir = self.options.mcDir
        sample_list = OrderedDict()
        if analysis == "HighMass":

            # qqZZ
            if self.options.powheg:
                # Powheg
                print "use PowHeg for qqZZ"
                sample_list['qqZZ']  = mc_dir + 'mc15_13TeV.361603.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZllll_mll4.root,'
                sample_list['qqZZ'] += mc_dir + 'mc15_13TeV.342556.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZllll_mll4_m4l_100_150.root,'
                sample_list['qqZZ'] += mc_dir + 'mc15_13TeV.343232.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZllll_mll4_m4l_500_13000.root,'
            elif self.options.sherpa == 2.1:
                print "use Sherpa 2.1 for qqZZ"
                # qqZZ, Sherpa, 2.1
                sample_list['qqZZ'] = mc_dir + 'mc15_13TeV.361090.Sherpa_CT10_llll_M4l100.root,'
            elif self.options.sherpa == 2.2:
                print "use Sherpa 2.2 for qqZZ"
                # qqZZ, Sherpa, 2.2.1,
                # It's likely the overlap is not removed...
                # don't use 345107, 345108, for now. 2017-03-01
                sample_list['qqZZ'] = mc_dir + 'mc15_13TeV.363490.Sherpa_221_NNPDF30NNLO_llll.root,' # full mass range,
                #sample_list['qqZZ'] = mc_dir + 'mc15_13TeV.345107.Sherpa_221_NNPDF30NNLO_llll_m4l100_300_filt100_150.root,' # 100, 130GeV
                #sample_list['qqZZ'] += mc_dir + 'mc15_13TeV.345108.Sherpa_221_NNPDF30NNLO_llll_m4l300.root,' # >= 300 GeV


            # ggZZ
            sample_list['ggZZ'] = mc_dir + 'mc15_13TeV.361073.Sherpa_CT10_ggllll.root,'

            # qqZZjj
            #sample_list['qqZZjj'] = mc_dir + 'mc15_13TeV.361072.Sherpa_CT10_lllljj_EW6.root,'
            sample_list['qqZZ'] += mc_dir + 'mc15_13TeV.361072.Sherpa_CT10_lllljj_EW6.root,'

            # reducible
            sample_list['reducible'] = self.get_reducible()

            # ttV
            sample_list['ttV']     = mc_dir  + 'mc15_13TeV.361621.Sherpa_CT10_WWZ_4l2v.root,'
            sample_list['ttV']    += mc_dir  + 'mc15_13TeV.361623.Sherpa_CT10_WZZ_5l1v.root,'
            sample_list['ttV']    += mc_dir  + 'mc15_13TeV.361625.Sherpa_CT10_ZZZ_6l0v.root,'
            sample_list['ttV']    += mc_dir  + 'mc15_13TeV.361626.Sherpa_CT10_ZZZ_4l2v.root,'
            #sample_list['ttV']    += mc_dir  + 'mc15_13TeV.410069.MadGraphPythia8EvtGen_A14NNPDF23LO_ttZllonshell_Np0.root' + ','
            #sample_list['ttV']    += mc_dir  + 'mc15_13TeV.410070.MadGraphPythia8EvtGen_A14NNPDF23LO_ttZllonshell_Np1.root' + ','
            sample_list['ttV']    += mc_dir  + 'mc15_13TeV.410144.Sherpa_NNPDF30NNLO_ttW.root,'
            sample_list['ttV']    += mc_dir  + 'mc15_13TeV.410142.Sherpa_NNPDF30NNLO_ttll_mll5.root,'

            # data
            sample_list['data']   = self.options.dataDir + 'data_13TeV.root,'
        else:
            print "I don't know"

        return sample_list

    def get_sys_list(self):
        analysis = self.options.analysis
        sys_dir = self.options.sysDir
        sys_list = OrderedDict()
        if analysis == "HighMass":
            # these text files follows the structure of the inputs for workspace
            sys_list['qqZZ'] = self.get_sys( sys_dir + 'norm_qqZZ.txt' )
            sys_list['ggZZ'] = self.get_sys( sys_dir + 'norm_ggllll.txt' )
            sys_list['qqZZjj'] = self.get_sys( sys_dir + 'norm_qqZZ.txt' )
            sys_list['ttV'] = self.get_sys( sys_dir + 'norm_qqZZ.txt' )
        else:
            pass

        return sys_list

    def get_sys(self, file_name):
        """This will read one of the workspace file"""
        sysMap = {}
        total_sys = 0.
        currSection = ''
        prev_section = None
        fileObj = open(file_name)
        for line in fileObj:
            line = line.strip()
            if(len(line) <= 0):
                continue

            if(line.startswith('#')):
                continue

            #print line, len(line)
            # check if the current line is a heading
            #isSection = re.match('^\[.*\]$', line)

            # update the current section and add a section to the map
            if '[' in line:
                currSection = line[1:-1].strip()
                if prev_section is not None:
                    sysMap[prev_section] = math.sqrt(total_sys)

                prev_section = currSection
                total_sys = 0
                continue


            # now these are options for section
            # split on '=' but not on '=='.. this uses lookForward and lookBehind regex search
            #print line[:-1].split('=')[1].split()
            try:
                low, high = line[:-1].split('=')[1].split()
                mean = (abs(float(low)-1) + abs(float(high)-1))/2.
                total_sys += mean**2
            except :
                pass
                #print line,"cannot be recognised"


        sysMap[prev_section] = math.sqrt(total_sys)
        return sysMap

    def get_yield(self, sample, cut):

        tree = ROOT.TChain(self.TREE_NAME, self.TREE_NAME)
        w_name = self.weight_name
        if "NNPDF30NNLO_llll" in sample:
            #print sample
            w_name += '*w_sherpaLep'

        ## use , to separate the samples.
        n_files = 0
        if ',' in sample:
            for s in sample.split(','):
                if s != "":
                    n_files += 1
                    tree.Add(s)
        else:
            tree.Add(sample)
            n_files += 1

        #print n_files,"Files with entries:", tree.GetEntries()
        if 'data' in sample:
            cut_t = ROOT.TCut(cut)
        else:
            if self.options.lumi > 0:
                lumi = self.options.lumi
            else:
                tree.GetEntry(0)
                lumi = tree.w_lumi

            cut_t = ROOT.TCut(w_name+"/w_lumi*"+str(lumi)+"*("+cut+")")
            #print "Luminosity:", round(lumi, 2),"fb-1"

        tree.Draw(self.options.poi+">>h1", cut_t)
        hist = ROOT.gDirectory.Get("h1")
        yields = hist.Integral()
        stats_error = yields/math.sqrt(hist.GetEntries())
        del hist

        return yields,stats_error

    def get_str(self, nominal, stats, sys):
        """
        return a str: "10 $\pm$ 10 $\pm$ 10"
        """
        split_sys = self.options.split
        dd = self.options.digits
        if split_sys:
            res = str(round(nominal, dd))+' $\pm$ '+str(round(stats, dd))+' $\pm$ '+str(round(sys,dd))
        else:
            total_ = math.sqrt(stats**2 + sys**2)
            res = str(round(nominal, dd))+' $\pm$ '+str(round(total_, dd))

        return res

    def process(self):
        samples = self.get_samples()
        cuts = self.get_cuts()
        sys = self.get_sys_list()

        out_text = ""
        ic = 0
        for sample_name in samples.keys():
            if ic == 0:
                out_text += sample_name
            else:
                out_text += " & " + sample_name
            ic += 1

        out_text += " \\\\ \\hline \n"

        combined = [0.]*len(samples)
        comb_stats = [0.]*len(samples)
        comb_sys = [0.]*len(samples)

        # summation of non-data for each category
        comb_chs = [0.]*len(cuts)
        comb_chs_stats = [0.]*len(cuts)
        comb_chs_sys = [0.]*len(cuts)

        if cuts is not None and len(samples) > 0:
            ichan = 0
            for ch_name,cut in cuts.iteritems():
                out_text += ch_name
                ic = 0
                for sample_name, sample_dir in samples.iteritems():

                    exp_ = stat_ = sys_ = 0
                    if type(sample_dir) is str:
                        exp_, stat_ = self.get_yield(sample_dir, cut)
                        if 'data' in sample_name:
                            sys_ = 0
                        else:
                            sys_ = sys[sample_name][ch_name]*exp_
                    else:
                        # for the pre-defined backgrounds, such as reducible-bkg
                        exp_, stat_, sys_ = sample_dir[ch_name]

                    combined[ic] += exp_
                    comb_stats[ic] += stat_**2
                    comb_sys[ic] += sys_**2

                    ic += 1
                    if 'data' in sample_name:
                        stats_ = math.sqrt(comb_chs_stats[ichan])
                        sys_ = math.sqrt(comb_chs_sys[ichan])
                        sum_bkg = self.get_str(comb_chs[ichan], stats_, sys_)
                        out_text += ' & ' + sum_bkg + ' & {:.0f}'.format(exp_)
                    else:
                        # sum of the samples, except data.
                        comb_chs[ichan] += exp_
                        comb_chs_stats[ichan] += stat_**2
                        comb_chs_sys[ichan] += sys_**2

                        # print out info
                        out_text += ' & ' + self.get_str(exp_, stat_, sys_)

                out_text += " \\\\ \n"
                ichan += 1
        else:
            print "I don't know"

        comb_tex = "Total"
        for icat, sample_name in enumerate(samples.keys()):

            if "data" in sample_name:
                # add total of non-data samples.
                sum_ = sum(comb_chs)
                stats_ = math.sqrt(sum(comb_chs_stats))
                sys_ = math.sqrt(sum(comb_chs_sys))
                sum_bkg = self.get_str(sum_, stats_, sys_)
                comb_tex += ' & '+sum_bkg+' & {:.0f}'.format(combined[icat])
            else:
                comb_tex += ' & '+self.get_str(combined[icat], math.sqrt(comb_stats[icat]), math.sqrt(comb_sys[icat]))

        comb_tex += "\\\\ \\hline \n"

        out_text += comb_tex
        print out_text

if __name__ == "__main__":
    usage = "%prog [options]"
    version="%prog 1.1"
    parser = OptionParser(usage=usage, description="get yields for WS", version=version)
    parser.add_option("--analysis", dest='analysis', default='HighMass', help='which analysis, affecting the built-in cuts')
    parser.add_option("--poi", dest='poi', default='m4l_constrained_HM', help='which variable used for counting')
    parser.add_option("--mcDir", dest='mcDir', default='/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/MiniTrees/Prod_v11/mc/Nominal/', help="directory for MC")
    parser.add_option("--dataDir", dest='dataDir', default='/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/MiniTrees/Prod_v10/data/Nominal/', help="directory for data")
    parser.add_option("--sysDir", dest='sysDir', help="directory for data", default="/Users/xju/Documents/Higgs/H4l/highmass/yields/")

    parser.add_option("--lumi", dest='lumi', default=-1, type='float', help='final luminosity')
    parser.add_option("--digits", dest='digits', default=2, type='int', help="digits in final numbers")
    parser.add_option("--split", dest='split', default=False, action='store_true', help="split stats and sys")

    parser.add_option("--test", dest='test', default=False, action='store_true', help="no VBF in highmass")
    ## change qqZZ samples
    parser.add_option("--powHeg", dest='powheg', default=False, action='store_true', help="use PowHeg for qqZZ")
    parser.add_option("--sherpa", dest='sherpa', default=2.2, type='float', help="Sherpa version")
    ## no VBF-like category in HighMass
    parser.add_option("--noVBF", dest='no_VBF', default=False, action='store_true', help="no VBF-like category")


    (options,args) = parser.parse_args()

    reader = MinitreeReader(options)
    reader.process()
