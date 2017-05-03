#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = "Xiangyang Ju"
__version__ = "0.1"
import ROOT
import helper

class Ploter:
    def __init__(self, status="Internal", lumi=36.1):
        self.status = status
        self.lumi = lumi

        # few options
        self.add_ratio = True

        self.has_data = False

        # predefined colors
        self.COLORS = [64, 95, 28, 29, 209, 5, 432, 433, 434, 435, 436, 8, 6]
        self.LINE_STYLE = [1]*10
        self.VerticalCanvasSplit = 0.4

        # parameters for label/title size
        self.t_size = 0.05
        self.x_title_size = 0.05
        self.text_size = 0.05

        # my canvas
        self.can = None
        self.pad1 = None
        self.pad2 = None

        # legend
        self.legend = None

        # atlas and legend offset
        self.x_offset = 0.20
        self.y_offset = 0.80

        # show sum of background for x-check
        self.show_sum_bkg = True

        self.totalObj = []


    def prepare_2pad_canvas(self, cname, width=600, height=600):

        self.can = ROOT.TCanvas(cname,cname,width,height)
        self.pad1 = ROOT.TPad("p1_"+cname,cname,0.0, self.VerticalCanvasSplit,1.0,1.0)
        self.pad2 = ROOT.TPad("p2_"+cname,cname,0.0,0.0, 1.0, self.VerticalCanvasSplit)
        self.pad1.SetBottomMargin(0)
        self.pad1.SetTopMargin(0.09)
        self.pad1.SetLeftMargin(0.17)
        #self.pad2.Draw()
        self.pad2.SetTopMargin(0)
        self.pad2.SetBottomMargin(0.4)
        self.pad2.SetLeftMargin(0.17)
        self.pad2.SetGridy()
        self.pad2.SetGridx()
        self.can.cd()
        self.pad1.Draw()
        self.pad2.Draw()

    def add_ratio_panel(self, hist_list, y_title,
                        y_min, y_max, reverse=False):
        """
        hist_list = [Data, MC1, MC2]
        plot Data/MC1, Data/MC2
        @para=reverse, plot MC1/Data
        """
        self.add_ratio = True
        #y_title = "Data/MC"
        #if reverse:
        #    y_title = "MC/Data"

        if len(hist_list) < 2:
            print "less than 2 histograms, kidding?"
            return None

        hist_list_cp = [x.Clone(x.GetName()+"_cloneRatio") for x in hist_list]
        self.totalObj.append(hist_list_cp)

        h_refer = hist_list_cp[0].Clone("Histreference")
        h_refer.Sumw2()
        self.totalObj.append(h_refer)
        # print "REFER:", h_refer.Integral()
        for i, hist in enumerate(hist_list_cp):
            if i==0:

                hist.Sumw2()
                hist.Divide(h_refer)
                hist.SetFillColor(1)
                hist.SetFillStyle(3010)
                hist.SetMarkerSize(0.001)

                labelscalefact = 1. / (1. - self.VerticalCanvasSplit)
                hist.GetYaxis().SetTitle(y_title)
                hist.GetYaxis().SetTitleSize(self.t_size*labelscalefact)
                hist.GetYaxis().SetLabelSize(self.t_size*labelscalefact)
                hist.GetYaxis().SetRangeUser(y_min, y_max)
                hist.GetYaxis().SetTitleOffset(0.8)

                hist.GetXaxis().SetLabelSize(self.t_size*labelscalefact)
                hist.GetXaxis().SetTitleSize(self.x_title_size*labelscalefact)
                hist.GetXaxis().SetTitleOffset(1.4)

                #hist.Draw("AXIS")
                hist.Draw("E2")
            else:
                # start to calculate the ratio
                if reverse: # MC/Data
                    this_hist = hist.Clone(hist.GetName()+"_cp")
                    this_hist.Divide(h_refer)
                else: # Data/MC
                    this_hist = h_refer.Clone(hist.GetName()+"_cpDI")
                    this_hist.SetLineColor(hist.GetLineColor())
                    this_hist.Divide(hist)
                    print "Yields:",hist.Integral(), h_refer.Integral()

                self.totalObj.append(this_hist)
                this_hist.Draw("EP SAME")
        helper.add_line(hist_list[0], 1.0)


    def stack_hists(self,
        hist_list, tag_list, out_name,
        x_title, y_title,
        is_log=False, has_data=True,
    ):
        ## not saving the plots..

        # In hist_list, the data should be first element, if has_data
        if len(hist_list) > len(self.COLORS):
            print "I don't have enough colors", len(hist_list), len(self.COLORS)
            return

        # clone current histograms, so that inputs are untouched.
        hist_list_cp = [] # a list of non-data histograms
        h_data = None  # the first element is assumed to be data

        hist_sum = None
        for i, hist in enumerate(hist_list):
            new_hist = hist.Clone(hist.GetName()+"_clone")
            color = self.COLORS[i]
            new_hist.SetLineColor(color)

            if i==0 and has_data:
                # decorate data points
                new_hist.SetMarkerStyle(20)
                new_hist.SetMarkerSize(1.2)
                h_data = new_hist
                self.get_offset(h_data)
                continue
            elif i==0 or hist_sum is None:
                hist_sum = hist
            else:
                hist_sum.Add(hist)

            new_hist.SetFillColor(color)
            hist_list_cp.append(new_hist)

        # always plot the smallest component in the bottom
        hist_sorted_list = sorted(hist_list_cp, key=lambda k:k.Integral())
        hs = ROOT.THStack("hs", "")
        for hist in hist_sorted_list:
            hs.Add(hist)

        # start to plot them
        if self.add_ratio and has_data:
            self.prepare_2pad_canvas("canvas", 600, 600)
            self.pad2.cd()
            self.pad2.SetGridy()
            hist_sum.SetLineColor(8)
            new_data_copy = h_data.Clone("data_copy")
            #self.add_ratio_panel([new_data_copy, hist_sum], y_title, 0.5, 1.52)
            self.add_ratio_panel([new_data_copy, hist_sum], "Data/MC", 0.55, 1.42)
            self.pad1.cd()
        else:
            self.can = ROOT.TCanvas("canvas", "canvas", 600, 600)

        y_max = hs.GetMaximum()
        y_min = hs.GetMinimum()
        if has_data and h_data is not None:
            if y_max < h_data.GetMaximum():
                y_max = h_data.GetMaximum()
            if y_min > h_data.GetMinimum():
                y_min = h_data.GetMinimum()

        if has_data:
            this_hist = h_data
        else:
            this_hist = hs

        if is_log:
            if self.add_ratio:
                self.pad1.SetLogy()
            else:
                self.can.SetLogy()
            this_hist.GetYaxis().SetRangeUser(4E-3, y_max*1e3)
        else:
            this_hist.GetYaxis().SetRangeUser(1E-3, y_max*1.1)

        #this_hist.SetNdivisions(8, "X")
        this_hist.SetXTitle(x_title)
        this_hist.SetYTitle(y_title)

        if has_data:
            h_data.Draw("EP")
            hs.Draw("HISTsame")
            h_data.Draw("AXISsame")
            h_data.Draw("EPsame")
        else:
            hs.Draw("HIST")

        # add legend
        hist_id = 0
        if has_data:
            hist_all = [h_data] + hist_list_cp
        else:
            hist_all = hist_list_cp

        if self.show_sum_bkg:
            legend = self.get_legend(len(hist_all)+1)
        else:
            legend = self.get_legend(len(hist_all))

        for hist, tag in zip(hist_all, tag_list):
            if has_data and hist_id == 0:
                legend.AddEntry(hist, tag+" {:.0f}".format(hist.Integral()), "LP")
                # add sum of background...
                if self.show_sum_bkg:
                    hist_sum.SetFillColor(0)
                    hist_sum.SetLineColor(0)
                    legend.AddEntry(hist_sum, "Total Bkg {:.0f}".format(hist_sum.Integral()), "F")
            else:
                legend.AddEntry(hist, tag+" {:.1f}".format(hist.Integral()), "F")
            hist_id += 1

        legend.Draw("same")
        self.add_atlas()
        self.add_lumi()

        if out_name is not None:
            self.can.SaveAs(out_name)


    def get_legend(self, nentries):
        x_min = self.x_offset
        x_max = x_min + 0.3
        y_max = self.y_offset-self.text_size*2-0.001
        y_min = y_max - self.t_size*nentries

        legend = ROOT.TLegend(x_min, y_min, x_max, y_max)
        legend.SetFillColor(0)
        legend.SetBorderSize(0)
        legend.SetTextFont(42)

        return legend


    def add_text(self, x, y, color, text, font=42):
        l = ROOT.TLatex()
        l.SetTextSize(self.text_size)
        l.SetNDC()
        l.SetTextColor(color)
        l.SetTextFont(font)
        l.DrawLatex(x, y, text)

    def add_atlas(self):
        self.add_text(self.x_offset, self.y_offset,
                      1, "#bf{#it{ATLAS}} "+self.status)

    def add_lumi(self):
        self.add_text(self.x_offset,
                      self.y_offset - self.text_size - 0.007,
                      1, "13 TeV, "+str(self.lumi)+" fb^{-1}")

    def get_offset(self, hist):
        max_bin = hist.GetMaximumBin()
        # nbins = hist.GetXaxis().GetNbins()
        last_bin = hist.GetXaxis().GetLast()
        first_bin = hist.GetXaxis().GetFirst()
        if max_bin < first_bin + (last_bin - first_bin)/2.:
            self.x_offset = 0.60


    def stack(self, hist_list):
        hist_list_cp = [] # a list of non-data histograms
        for hist in hist_list:
            hist_list_cp.append( hist.Clone(hist.GetName()+"stackClone"))

        self.totalObj.append(hist_list_cp)

        hist_sorted_list = sorted(hist_list_cp, key=lambda k:k.Integral())
        hs = ROOT.THStack("hs", "")
        hist_sum = None
        for hist in hist_sorted_list:
            hs.Add(hist)
            if hist_sum is None:
                hist_sum = hist.Clone(hist.GetName()+"sumClone")
            else:
                hist_sum.Add(hist)

        self.totalObj.append(hist_sum)
        return hist_sum, hs

    def color(self, hist_list, no_fill=False):
        for i, hist in enumerate(hist_list):
            color_ = self.COLORS[i]
            hist.SetLineColor(color_)
            if not no_fill:
                hist.SetFillColor(color_)
            hist.SetMarkerColor(color_)
            hist.SetLineStyle(self.LINE_STYLE[i])

    def set_y_range(self, hist_data, hist_splusb, is_logY):
        y_max = hist_splusb.GetMaximum()
        y_min = hist_splusb.GetMinimum()
        if hist_data:
            this_hist = hist_data
            if y_max < hist_data.GetMaximum():
                y_max = hist_data.GetMaximum()

            if y_min > hist_data.GetMinimum():
                y_min = hist_data.GetMinimum()

        else:
            this_hist = hist_splusb


        if is_logY:
            if self.add_ratio:
                self.pad1.SetLogy()
            else:
                self.can.SetLogy()

            # this_hist.GetYaxis().SetRangeUser(4E-3, y_max*1e2)
        else:
            this_hist.GetYaxis().SetRangeUser(1E-3, y_max*1.10)
            pass

        return this_hist

    def compare_hists(self, hist_list, tag_list, **kwargs):
        """
        a list of histograms,
        Key words include:  
            ratio_title, ratio_range, logY, out_name
            no_fill, x_offset, draw_option, add_yields, 
            out_folder
        """
        self.del_obj()

        if len(hist_list) < 2:
            print "not enough hitograms for comparison"
            return 
        try:
            no_fill = kwargs["no_fill"]
        except KeyError:
            no_fill = False

        self.color(hist_list, no_fill)
      
        if self.add_ratio:
            self.prepare_2pad_canvas('canvas', 600, 600)
            self.pad2.cd()
            try:
                ratio_title = kwargs["ratio_title"]
            except KeyError:
                ratio_title = "MC/Data"

            try:
                ratio_x, ratio_y = kwargs["ratio_range"]
            except KeyError:
                ratio_x, ratio_y = 0.55, 1.42
        
            self.add_ratio_panel(hist_list, ratio_title, ratio_x, ratio_y, True)
            self.pad1.cd()
        else:
            self.text_size = 0.04
            self.can = ROOT.TCanvas("canvas", "canvas", 600, 600)

        self.set_y_offset()
        try:
            self.x_offset = kwargs["x_offset"]
        except KeyError:
            self.get_offset(hist_list[0])

        try:
            is_logy = kwargs["logY"]
        except KeyError:
            is_logy = False

        legend = self.get_legend(len(hist_list))

        this_hist = self.set_y_range(hist_list[0], hist_list[1], is_logy)
        # y_axis = this_hist.GetMaximum()
        # this_hist.GetYaxis().SetRangeUser(0, y_axis*1.5)
        try:
            draw_option = kwargs["draw_option"]
        except KeyError:
            draw_option = "HIST"

        try:
            add_yield = kwargs['add_yields']
        except KeyError:
            add_yield = False

        for i, hist in enumerate(hist_list):
            if add_yield:
                legend.AddEntry(hist, "{}: {:.0f}".format(tag_list[i], hist.Integral()))
            else:
                legend.AddEntry(hist, tag_list[i])

            if i==0:
                hist.Draw(draw_option)
            else:
                hist.Draw(draw_option+" SAME")

        legend.Draw("same")
        self.add_atlas()
        self.add_lumi()

        try:
            out_name = kwargs["out_name"]
        except KeyError:
            out_name = "TEST"

        try:
            out_folder = kwargs["out_folder"]
        except KeyError:
            out_folder = "./"

        helper.mkdir_p(out_folder)


        if is_logy:
            self.can.SaveAs(out_folder+"/"+out_name+"_Log.eps")
        else:
            self.can.SaveAs(out_folder+"/"+out_name+".eps")


    def set_y_offset(self):
        if not self.add_ratio:
            self.y_offset = 0.88

    def del_obj(self):
        for obj in self.totalObj:
            del obj

        if self.can:
            del self.can

        if self.pad1:
            del self.pad1

        if self.pad2:
            del self.pad2
