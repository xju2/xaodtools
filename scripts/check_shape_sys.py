#/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import helper
from optparse import OptionParser
from ploter import Ploter
import AtlasStyle

class ShapeSysReader(object):
    """
    compare the input shape systematic uncertainty for each signal
    Input is: 1) a list of nuisances, 2) a systematic input generated by HZZ framework 
    """
    def __init__(self, input_name, np_list_name, options):
        self.fin = ROOT.TFile.Open(input_name)
        self.np_list = helper.read_np_list(np_list_name)
        self.obs_name = "m4l"
        self.categories = ["ggF_2mu2e_13TeV", "ggF_4e_13TeV", "ggF_4mu_13TeV"]
        self.sample_name = options.sample_name
        self.x_min, self.x_max = map(lambda x: float(x), options.range.split(','))

        self.ps = Ploter("Internal", 36.1)
        
        self.plot_option = {}
        self.plot_option["ratio_title"] = "Variation/Nominal"
        self.plot_option["ratio_range"] = (0.55, 1.42)
        self.plot_option["logY"] = False
        self.plot_option["out_folder"] = "shape_sys/"+self.sample_name
        #self.plot_option['x_offset'] = 0.20

    def loop_category(self, np_name):
        """
        compare the nominal and variations
        """
        for category in self.categories:
            norm_name = "-".join([self.obs_name, self.sample_name, category])
            up_name = "-".join([self.obs_name, np_name, category, "up"])
            down_name = "-".join([self.obs_name, np_name, category, "down"])

            h_norm = self.fin.Get(norm_name)
            h_up = self.fin.Get(up_name)
            h_down = self.fin.Get(down_name)

            if h_norm and h_up and h_down:
                hist_list = [h_norm, h_up, h_down]
                for hist in hist_list:
                    hist.GetXaxis().SetRangeUser(self.x_min, self.x_max)
                    hist.GetXaxis().SetTitle("m_{4l} [GeV]")
                
                cat_name = category.replace("ATLAS_", "")
                tag_list = ["Nominal", np_name+" Up", np_name+" Down"]
                self.compare(hist_list, tag_list, "_".join([np_name, category]))
            else:
                print norm_name,up_name,down_name

    def compare(self, hist_list, tag_list, out_name):
        self.plot_option["out_name"] = out_name

        self.ps.compare_hists(hist_list, tag_list, **self.plot_option)


    def process(self):
        """
        main function, loop over nuisance parameters
        """
        for np_name in self.np_list:
            self.loop_category(np_name)

if __name__ == "__main__":
    usage = "%prog [option] input_file np_list"
    parser = OptionParser(description="compare shape variations", usage=usage)
    parser.add_option('--obs', dest="obs", help="observable name", default="m4l")
    parser.add_option('-s', dest="sample_name", help="Name of the sample", default="ggF200")
    parser.add_option('-r', dest="range", help="range of observable", default="0,240")

    options,args = parser.parse_args()

    if len(args) < 2:
        print parser.print_help()
        exit(1)

    reader = ShapeSysReader(args[0], args[1], options)
    reader.process()
