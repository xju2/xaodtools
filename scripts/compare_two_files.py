#!/usr/bin/env python
import ROOT
import helper
from optparse import OptionParser

class FileCompare:
    def __init__(self, f1_name, t1_name, f2_name, t2_name):
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

    def compare(self):
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
                print "NO!--",run,event,ch1.m4l_unconstrained
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

        self.out_text += str(n_matched)+" "+str(n_matched/n_f1)+\
                " matches found in "+self.f2_name+"\n"
        self.out_text += str(n_missed)+" "+str(n_missed/n_f1)+\
                " no matches found in "+self.f2_name+"\n"

        fout = ROOT.TFile.Open("diff_mass.root", "recreate")
        for hist in self.hist_list:
            hist.Write()

        for h2d in self.h2d_list:
            h2d.Write()

        fout.Close()
        print self.out_text


if __name__ == "__main__":
    usage = "%prog [option] f1 t1 f2 t2"
    parser = OptionParser(description="compare two TTree", usage=usage)

    options,args = parser.parse_args()
    if len(args) < 4:
        print parser.print_help()
        exit(1)

    fc = FileCompare(args[0], args[1], args[2], args[3])
    fc.compare()
