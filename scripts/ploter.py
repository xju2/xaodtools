#!/usr/bin/env python
__author__ = "Xiangyang Ju"
__version__ = "0.1"
__doc__ = """
to make comparison plots
"""
import ROOT

class Ploter:
    def __init__(self, status="Internal", lumi=36.47):
        self.status = status
        self.lumi = lumi
        self.totalObj = []

        # predefined colors
        self.COLORS = [1, 206, 64, 95, 28, 29, 209, 5]
        self.VerticalCanvasSplit = 0.4

        # parameters for label/title size
        self.t_size = 0.04
        self.x_title_size = 0.05
        self.text_size = 0.05

        # my canvas
        self.can = None
        self.pad1 = None
        self.pad2 = None

    def prepare_2pad_canvas(self, cname, width=600, height=600):

        self.can = ROOT.TCanvas(cname,cname,width,height)
        self.pad1 = ROOT.TPad("p1_"+cname,cname,0.0, self.VerticalCanvasSplit,1.0,1.0)
        self.pad2 = ROOT.TPad("p2_"+cname,cname,0.0,0.0, 1.0, self.VerticalCanvasSplit)
        self.pad1.SetBottomMargin(0)
        self.pad1.SetTopMargin(0.09)
        self.pad1.SetLeftMargin(0.17)
        self.pad2.Draw()
        self.pad2.SetTopMargin(0)
        self.pad2.SetBottomMargin(0.4)
        self.pad2.SetLeftMargin(0.17)
        self.pad2.SetGridy()
        self.can.cd()
        self.pad1.Draw()
        self.pad2.Draw()

    def add_ratio_panel(self, hist_list_cp, y_title,
                        y_min, y_max, reverse=False):
        """
        hist_list = [Data, MC1, MC2]
        plot Data/MC1, Data/MC2
        @para=reverse, plot MC1/Data
        """
        if len(hist_list_cp) < 2:
            print "less than 2 histograms, kidding?"
            return None

        #hist_list_cp = [x.Clone(x.GetName()+"_clone") for x in hist_list]
        h_refer = hist_list_cp[0].Clone("Histreference")
        self.totalObj.append(h_refer)
        print "REFER:", h_refer.Integral()
        for i, hist in enumerate(hist_list_cp):
            if i==0:

                hist.Sumw2()
                hist.Divide(h_refer)
                hist.SetFillColor(1)
                hist.SetFillStyle(3010)
                hist.SetMarkerSize(0.001)

                labelscalefact = 1. / (1. - self.VerticalCanvasSplit)
                hist.GetYaxis().SetTitleSize(self.t_size*labelscalefact)
                hist.GetYaxis().SetLabelSize(self.t_size*labelscalefact)
                hist.GetXaxis().SetLabelSize(self.t_size*labelscalefact)

                hist.GetXaxis().SetTitleSize(self.x_title_size*labelscalefact)

                hist.GetYaxis().SetTitle(y_title)
                hist.GetYaxis().SetRangeUser(y_min, y_max)

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
                this_hist.Draw("HIST SAME")


    def stack_hists(self,
        hist_list, tag_list, out_name,
        x_title, y_title,
        is_log=False, has_data=True,
        add_ratio=False
    ):
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
        if add_ratio and has_data:
            self.prepare_2pad_canvas("canvas", 600, 600)
            canvas = self.can
            self.pad2.cd()
            hist_sum.SetLineColor(4)
            new_data_copy = h_data.Clone("data_copy")
            self.add_ratio_panel([new_data_copy, hist_sum], y_title, 0.5, 1.52)
            self.pad1.cd()
        else:
            canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)

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
            if add_ratio:
                self.pad1.SetLogy()
            else:
                canvas.SetLogy()
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

        legend = self.get_legend(len(hist_all))
        for hist, tag in zip(hist_all, tag_list):
            if has_data and hist_id == 0:
                legend.AddEntry(hist, tag+" {:.0f}".format(hist.Integral()), "LP")
            else:
                legend.AddEntry(hist, tag+" {:.1f}".format(hist.Integral()), "F")
            hist_id += 1

        legend.Draw("same")
        x_offset = 0.17
        self.add_atlas(x_offset, 0.80)
        self.add_lumi(x_offset, 0.80 - self.text_size - 0.007)

        canvas.SaveAs(out_name)



    def get_legend(self, nentries, corner="RT"):
        # corner = LT, RT, LB, RB
        # LT = left-hand top
        # LB = left-hand bottom
        x_min = 0.6
        x_max = 0.9
        y_min = 0.55
        y_max = y_min + self.t_size*nentries

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

    def add_atlas(self, x, y):
        self.add_text(x, y, 1, "#bf{#it{ATLAS}} "+self.status)

    def add_lumi(self, x, y):
        self.add_text(x, y, 1, "13 TeV, "+str(self.lumi)+" fb^{-1}")
