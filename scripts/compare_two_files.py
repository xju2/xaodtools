#!/usr/bin/env python
import ROOT
import helper
from optparse import OptionParser
from ploter import Ploter

class FileCompare:
    def __init__(self, f1_name, t1_name, f2_name, t2_name, options):
        import AtlasStyle
        self.f1_name = f1_name
        self.t1_name = t1_name
        self.f2_name = f2_name
        self.t2_name = t2_name

        self.out_text = ""
        hist_temp = ROOT.TH1F("hist_temp", "diff", 100, -10, 10)
        self.hist_list = []

        self.ch_names = ["4mu", "4e", "2mu2e", "2e2mu"]
        self.mass_bins = ["115", "115_130", "130_200", "200_500", "500_"]
        self.mass_bin_boundies = [115, 130, 200, 500]

        self.n_ch = len(self.ch_names)
        self.n_mass_bins = len(self.mass_bins)
        for ch_name in self.ch_names:
            for mass_bin in self.mass_bins:
                self.hist_list.append(hist_temp.Clone("hist_"+\
                                                      ch_name+"_"+mass_bin))

        # book 2D histograms
        h2d_temp = ROOT.TH2F("h2d_temp", "m4l vs diff", 50, 60, 1200, 100, -10, 10)
        self.h2d_list = []
        for ch_name in self.ch_names:
            self.h2d_list.append( h2d_temp.Clone("h2d_"+ch_name) )

        self.lumi = 36.1
        self.ps = Ploter("Internal", self.lumi)

        self.options = options

    def get_bin_index(self, event_type, m4l):
        mass_id = -1
        for imass in range(len(self.mass_bin_boundies)):
            mass_cut = self.mass_bin_boundies[imass]
            if m4l < mass_cut:
                mass_id = imass
                break

            if imass == len(self.mass_bin_boundies)-1 and m4l > mass_cut:
                mass_id = imass + 1

        final_id = (mass_id + event_type*self.n_mass_bins)
        #print event_type,m4l,mass_id,final_id
        return final_id

    def compare_ebye(self):
        ch1 = helper.loader(self.f1_name, self.t1_name)
        ch2 = helper.loader(self.f2_name, self.t2_name)
        ch2.BuildIndex("run", "event")

        n_f1 = ch1.GetEntries()
        n_f2 = ch2.GetEntries()

        self.out_text = self.f1_name+" contains "+str(n_f1)+" events\n"
        self.out_text += self.f2_name+" contains "+str(n_f2)+" events\n"

        n_matched = 0
        n_missed = 0

        for iEvt in range(n_f1): # loop file1
            #ch1.LoadTree(iEvt)
            ch1.GetEntry(iEvt)

            #found = ch2.GetEntryWithIndex(ch1.run, ch1.event)
            run = ch1.run
            event = ch1.event
            # the match scheme is stupid, not necessary to have precise match!
            found_entry = ch2.GetEntryNumberWithBestIndex(run, event)
            if found_entry < 0:
                n_missed += 1
                print run,event,"found no matches"
                continue
            else:
                ch2.GetEntry(found_entry)

            if ch2.run != run or ch2.event != event:
                print "NO!--",run,event,ch1.event_type,ch1.m4l_unconstrained,ch1.m4l_constrained_HM
                n_missed += 1
                continue
            else:
                n_matched += 1
                bin_id = self.get_bin_index(ch1.event_type,
                                            ch1.m4l_unconstrained)
                hist = self.hist_list[bin_id]

                m4l_diff = ch1.m4l_unconstrained - ch2.m4l_unconstrained
                if m4l_diff > 100:
                    print ch2.run,ch2.event,ch1.run,ch1.event,ch1.event_type,ch1.m4l_unconstrained,ch2.m4l_unconstrained,m4l_diff,m4l_diff/ch1.m4l_unconstrained
                hist.Fill(m4l_diff)

                h2d = self.h2d_list[ch1.event_type]
                h2d.Fill(ch1.m4l_unconstrained, m4l_diff)

        self.out_text += str(n_matched)+" "+str(round(n_matched*1./n_f1,3))+\
                " matches found in "+self.f2_name+"\n"
        self.out_text += str(n_missed)+" "+str(round(n_missed*1./n_f1,3))+\
                " no matches found in "+self.f2_name+"\n"

        fout = ROOT.TFile.Open("diff_mass.root", "recreate")
        for hist in self.hist_list:
            hist.Write()

        for h2d in self.h2d_list:
            h2d.Write()

        fout.Close()
        print self.out_text

    def compare_spectrum(self, add_cuts, out_name):
        ch1 = helper.loader(self.f1_name, self.t1_name)
        ch2 = helper.loader(self.f2_name, self.t2_name)

        h4l_temp = ROOT.TH1F("h4l_temp", "temp", 70, 130, 1530)
        h4l_ch1 = h4l_temp.Clone("h4l_ch1")
        h4l_ch2 = h4l_temp.Clone("h4l_ch2")

    
        cut = ROOT.TCut("weight*(pass_vtx4lCut==1 && m4l_constrained_HM > 130 && m4l_constrained_HM < 1500 && "+add_cuts+")")
        ch1.Draw("m4l_constrained_HM>>h4l_ch1", cut)
        ch2.Draw("m4l_constrained_HM>>h4l_ch2", cut)

        h4l_ch2.SetXTitle("m_{4l} [GeV]")
        h4l_ch1.GetXaxis().SetTitle("m_{4l} [GeV]")
        h4l_ch2.SetYTitle("Events/5 GeV")

        h4l_ch1.Sumw2()
        h4l_ch2.Sumw2()
        h4l_ch1.SetMarkerSize(0.5)
        h4l_ch1.SetMarkerStyle(4)
        h4l_ch2.SetMarkerSize(0.5)

        hist_list = [h4l_ch1, h4l_ch2]
        self.ps.color(hist_list)

        self.ps.prepare_2pad_canvas('canvas', 600, 600)
        self.ps.pad2.cd()
        self.ps.add_ratio_panel(hist_list, "v10/v09", 0.55, 1.42, True)
        self.ps.pad1.cd()
        self.ps.get_offset(h4l_ch1)
        is_logy = self.options.logy

        legend = self.ps.get_legend(len(hist_list) + 2)

        this_hist = self.ps.set_y_range(h4l_ch1, h4l_ch2, is_logy)
        y_axis = h4l_ch2.GetMaximum()
        h4l_ch2.GetYaxis().SetRangeUser(0, y_axis*1.5)
        h4l_ch2.Draw("EP")
        h4l_ch1.Draw("EP SAME")

        legend.AddEntry(h4l_ch1, "v09: {:.0f}".format(h4l_ch1.Integral()), "EP")
        legend.AddEntry(h4l_ch2, "v10: {:.0f}".format(h4l_ch2.Integral()), "EP")

        legend.Draw("same")
        self.ps.add_atlas()
        self.ps.add_lumi()

        if is_logy:
            self.ps.can.SaveAs(out_name+"_Log.eps")
        else:
            self.ps.can.SaveAs(out_name+".eps")


if __name__ == "__main__":
    usage = "%prog [option] f1 t1 f2 t2 out_name"
    parser = OptionParser(description="compare two TTree", usage=usage)
    parser.add_option('--logY', dest="logy", help="Log-Y", default=False, action='store_true')
    parser.add_option('--cmp', dest="cmp", help="compare the shapes", default=False, action='store_true')

    options,args = parser.parse_args()

    if len(args) < 4:
        print parser.print_help()
        exit(1)

    fc = FileCompare(args[0], args[1], args[2], args[3], options)
    out_name = args[4]
    if options.cmp:
        fc.compare_spectrum("1", out_name+"_v09_v10_4l")
        fc.compare_spectrum("event_type==0", out_name+"_v09_v10_4mu")
        fc.compare_spectrum("event_type==1", out_name+"_v09_v10_4e")
        fc.compare_spectrum("(event_type==2||event_type==3)", out_name+"_v09_v10_2e2mu")
    else:
        fc.compare_ebye()
