#ifndef __loader_c__
#define __loader_c__

#include "/afs/cern.ch/user/x/xju/tool/AtlasStyle.C"
#include "/afs/cern.ch/user/x/xju/tool/AtlasUtils.C"
#include <TString.h>
#include <exception>
#include <TH1F.h>
#include <TVector2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TChain.h>
#include <TCut.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TMath.h>

#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <RooTFnPdfBinding.h>
#include <RooBinning.h>

#include <iostream>
#include <stdio.h>
#include <string>

using namespace std;
using namespace RooFit;
int my_color_list[] = {kRed, kBlue, kGreen+1, kOrange+6, kAzure+2, kViolet-3, kAzure+5, kOrange+6, kViolet-9, kAzure+6, kViolet-4};

TChain* loader(const char* inFile_name, const char* chain_name = "physics")
{
    TChain* chain = new TChain(chain_name);
    TString in_name(inFile_name);
    if(in_name.Contains("root")) {
        chain->Add(inFile_name);
        return chain;
    }
    fstream input(inFile_name, fstream::in);
    string file_name;
    int ncounter = 0;
    while (input >> file_name){
        // cout << "adding: " << file_name << endl;
        if (file_name.at(0) == '#') continue;
        chain->Add(file_name.c_str());
        ncounter ++;
    }
    cout << "total events: " << chain->GetEntries() << " in " << ncounter << " files." << endl;
    input.close();
    return chain;
}

void SetAtlasHist(TH1* h_copy){
    h_copy->GetYaxis()->SetRangeUser(0.6, 1.8);
    h_copy->GetYaxis()->SetNdivisions(7);
    h_copy->GetYaxis()->SetLabelSize(0.1);
    h_copy->GetYaxis()->SetTitleSize(0.15);
    h_copy->GetYaxis()->SetTitleOffset(0.3);
    h_copy->GetXaxis()->SetLabelSize(0.1);
    h_copy->GetXaxis()->SetTitleSize(0.15);
    h_copy->GetXaxis()->SetTitleOffset(0.78);
}

typedef struct BranchInfo {
    string name_;
    int n_;
    float low_;
    float high_;
    BranchInfo(){
        name_ = "";
        n_ = 100;
        low_ = 0;
        high_ = 1000;
    }
} BranchInfo;

TH1F* draw_hist_from_chain(TChain* chain, const char* branch_name, 
        const TCut& cut, const char* hist_name, 
        int n_bins, float low_value, float high_value)
{
    TH1F* h1 = new TH1F(hist_name, "h1", n_bins, low_value, high_value);
    chain->Draw(Form("%s >> %s", branch_name, hist_name), cut);
    h1->SetDirectory(0);
    return h1;
}

TH1F* draw_hist_from_chain(TChain* chain, const TCut& cut, 
        const char* hist_name, const BranchInfo& br)
{
    return draw_hist_from_chain(chain, br.name_.c_str(), cut, hist_name, 
            br.n_, br.low_, br.high_);
}

TH1F* draw_hist_from_file(const char* file_name, const char* chain_name, 
        const TCut& cut, const char* hist_name, const BranchInfo& br){
    TChain* chain = loader(file_name, chain_name);
    return draw_hist_from_chain(chain, br.name_.c_str(), cut, hist_name, 
            br.n_, br.low_, br.high_);
}

TH1F* create_hist(const char* file_name, const char* chain_name,
        const char* branch_name, const TCut& cut, const char* hist_name,
        int n_bins, float low_value, float high_value, int color)
{
    TChain* chain = loader(file_name, chain_name);

    TH1F* h1 = new TH1F(hist_name, "h1", n_bins, low_value, high_value);
    chain->Draw(Form("%s >> %s", branch_name, hist_name), cut);

    h1->SetDirectory(0);
    h1->SetLineColor(color);
    h1->SetMarkerSize(0.03);
    delete chain;
    return h1;
}

void save_hist(TH1F* h1, const char* out_file_name, 
        const char* hist_name = "met_all")
{
    TFile* file = TFile::Open(out_file_name, "UPDATE");
    TString oldName(h1->GetName());
    h1->SetName(hist_name);
    h1->Write();
    file->Close();
    h1->SetName(oldName.Data());
}

TPad* add_ratio_pad(TH1* h_signal, const TList& h_bkgs)
{
    TPad* pad1 = new TPad("pad1","pad1",0, 0.15,1.00,1.00);
    pad1->SetBottomMargin(0.161);
    TPad* pad2 = new TPad("pad2","pad2",0., 0.010, 1.00, 0.286);
    pad2->SetBottomMargin(0.28);
    pad1->Draw();
    pad2->Draw();

    // plot the ratio
    pad2->cd();
    // add error of signal
    TH1F* h_signal_copy = (TH1F*) h_signal->Clone("signal_copy");
    h_signal_copy->Sumw2();
    h_signal_copy->Divide(h_signal);
    h_signal_copy->SetFillColor(h_signal->GetFillColor());

    // cout << "background size: " << h_bkgs.GetSize() << endl;
    TIter next(h_bkgs.MakeIterator());
    TH1F* h_bkg = NULL;
    int icount = 0;
    while ((h_bkg = (TH1F*) next())){
        TH1F* h_copy = (TH1F*) h_bkg->Clone(Form("r_%s", h_bkg->GetName()));
        if(!h_copy->GetDefaultSumw2())  h_copy->Sumw2();
        h_copy->Divide(h_signal);
        h_copy->SetMarkerSize(0.02);
        h_copy->SetMarkerColor(h_bkg->GetMarkerColor());
        h_copy->SetLineColor(h_bkg->GetLineColor());
        if(icount == 0){
            h_copy->GetYaxis()->SetRangeUser(0.8, 1.5);
            h_copy->GetYaxis()->SetNdivisions(7);
            h_copy->GetYaxis()->SetLabelSize(0.1);
            h_copy->GetYaxis()->SetTitleSize(0.15);
            h_copy->GetYaxis()->SetTitle("Ratio");
            h_copy->GetYaxis()->SetTitleOffset(0.3);
            h_copy->GetXaxis()->SetLabelSize(0.1);
            h_copy->GetXaxis()->SetTitleSize(0.15);
            h_copy->GetXaxis()->SetTitleOffset(0.78);
            h_copy->Draw("HIST");
            h_signal_copy->Draw("SAME E3");
            h_copy->Draw("SAME HIST");
            icount ++;
        } else {
            h_copy->Draw("HIST SAME");
        }
    }

    AddLine(h_signal, 1, h_signal->GetLineColor(), h_signal->GetLineStyle());


    // plot the comparison
    pad1->cd();
    return pad1;
}

TPad* add_ratio_pad(TH1* h_signal, TH1* h_bkg)
{
    TList* list = new TList();
    list->Add(h_bkg);
    TPad* res = add_ratio_pad(h_signal, *list);
    delete list;
    return res;
}

void norm_hist(TH1* h1){
    h1->Scale(1./h1->Integral());
}

void CheckNull(TObject* obj)
{
    if(obj == NULL) {
        throw invalid_argument(Form("%s does not exist", obj->GetName()));
    }
}

void SetAtlasStyleCanvas(TCanvas* canvas, bool for_2d = false)
{
    if(!canvas) return;
    canvas->SetTopMargin(0.05);
    canvas->SetRightMargin(0.05);
    canvas->SetBottomMargin(0.16);
    canvas->SetLeftMargin(0.16);
    if(for_2d){
        canvas->SetLeftMargin(0.15);
        canvas->SetRightMargin(0.14);
    }

    // use large fonts
    //Int_t font=72; // Helvetica italics
    Int_t font=42; // Helvetica
    Double_t tsize=0.05;
    Int_t icol = 0;
    canvas->SetFrameBorderMode(icol);
    canvas->SetFrameFillColor(icol);
    canvas->SetBorderMode(icol);

    canvas->SetTickx(1);
    canvas->SetTicky(1);
}

void SetAtlasStyleHist(TH1* h1){
    if(!h1) return;
    TString className(h1->ClassName());
    int font = 42;
    double tsize = 0.05;
    h1->SetLineWidth(2.);
    if(className.Contains("TH1") || className.Contains("TH2")){
        TAxis* xaxis = h1->GetXaxis();
        xaxis->SetTitleFont(font);
        xaxis->SetTitleSize(tsize);
    } 
    if (className.Contains("TH2")){
        TAxis* yaxis = h1->GetYaxis();
        yaxis->SetTitleFont(font);
        yaxis->SetTitleSize(tsize);
    }
}

void SetAtlasOpt(){
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
}

void add_hist(TList* objarray, TH1F* h1)
{
    //cout << "adding: " << h1->GetName() << endl;
    // items in TList are sorted by their area
    int total = objarray->GetEntries();
    if(total < 1){
        objarray->Add(h1);
        return;
    } else{
        int index = total-1;
        int pos = 0;
        while(index >= 0){
            TH1F* h2 = (TH1F*) objarray->At(index);
            if(h2->Integral() < h1->Integral()){
                objarray->AddAfter(h2, h1);
                break;
            }
            index --;
        }
        if(index < 0){
            TH1F* h2 = (TH1F*) objarray->At(0);
            objarray->AddBefore(h2, h1);
        }
        return;
    }
}

TH1* SumHistsWithSysUncertainties(const TList& hists, 
        const vector<double>& systematics, bool only_sys) 
{
    TH1* h_all = NULL;
    TIter next(hists.MakeIterator());
    TH1* hist = NULL;
    int icount = 0;
    while ((hist = (TH1*) next())) {
        double sys = systematics.at(icount);
        if(icount == 0){
            h_all = (TH1*) hist->Clone("h_all");
            for(int ibin = 1; ibin < h_all->GetNbinsX()+1; ibin++){
                double origin_stats = h_all->GetBinError(ibin);
                double content = h_all->GetBinContent(ibin);
                double sys_unc = sys*content;
                double stats = only_sys?sys_unc*sys_unc:origin_stats*origin_stats + sys_unc*sys_unc;
                h_all->SetBinError(ibin, stats);
            }
        } else {
            for(int ibin = 1; ibin < hist->GetNbinsX()+1; ibin++){
                double content = hist->GetBinContent(ibin);
                double stats = hist->GetBinError(ibin);
                double origin_stats = h_all->GetBinError(ibin);
                double sys_unc = sys*content;
                double uncertainty = only_sys?origin_stats+sys_unc*sys_unc:origin_stats + stats*stats + sys_unc*sys_unc;
                double origin_content = h_all->GetBinContent(ibin);
                h_all->SetBinError(ibin, (float) uncertainty);
                h_all->SetBinContent(ibin, (float) origin_content+content);
            }
        }
        icount ++;
    }
    for (int ibin = 1; ibin < h_all->GetNbinsX()+1; ibin++){
        double uncertainty = h_all->GetBinError(ibin);
        h_all->SetBinError(ibin, sqrt(uncertainty)); 
    }
    return h_all;
}

void add_hist(TList* objarray, const string& hist_name)
{
    TH1F* h1 = (TH1F*) gDirectory->Get(hist_name.c_str());
    if(!h1){
        cout <<"Cannot find: " << hist_name << endl;
        return;
    }
    add_hist(objarray, h1);
}

double get_significance(double s, double b)
{
    if(s < 0 || b < 0){ 
        return -1;
    }
    if(fabs(b) < 1e-6){ 
        return -1;
    }
    return TMath::Sqrt(2*((s+b)*TMath::Log(1+s/b)-s));
}

double get_significance_with_sysB(double s, double b, double sigmaB){
    if(s < 0 || b < 0 || sigmaB < 0){ 
        return -1;
    }
    if(fabs(b) < 1e-6){ 
        return -1;
    }
    double sqB = sigmaB*sigmaB;
    return TMath::Sqrt(2*((s+b)*TMath::Log((s+b)*(b+sqB)/(b*b+(s+b)*sqB))
                - b*b/sqB*TMath::Log(1 + sqB*s/b/(b+sqB)))); 
}

int get_roc(TH1F* sig, TH1F* bkg, bool reverse = false)
{
    int n_bins = 100;
    if(sig->GetNbinsX() != bkg->GetNbinsX())
    {
        cout << "signal and background have different bins" << endl;
        return -1 ;
    }
    n_bins = sig->GetNbinsX();
    TH1F* roc = (TH1F*) sig->Clone("roc");
    TH1F* eff_sig = (TH1F*) sig->Clone("eff_sig");
    TH1F* eff_bkg = (TH1F*) bkg->Clone("eff_bkg");
    eff_sig->GetYaxis()->SetRangeUser(0, 1.5);
    double max_significance = -99;
    int max_index = -1;
    double sig_total = sig->Integral();
    double bkg_total = bkg->Integral();
    for(int i = 1; i <= n_bins; i++)
        // for(int i = n_bins; i >= 1; i--)
    {
        double s = sig->Integral(i, n_bins);
        double b = bkg->Integral(i, n_bins);
        // double s = sig->Integral(1, i);
        // double b = bkg->Integral(1, i);
        double s_over_b = get_significance(s, b);
        if(s_over_b+1 != s_over_b 
                && !TMath::IsNaN(s_over_b)
                && s_over_b > max_significance)
        { 
            max_significance = s_over_b;
            max_index = i;
        }
        if(s_over_b+1 != s_over_b && !TMath::IsNaN(s_over_b)){
            roc->SetBinContent(i, s_over_b);
        }else{ 
            roc->SetBinContent(i, 0); 
        }
        eff_sig->SetBinContent(i, s/sig_total);
        eff_bkg->SetBinContent(i, b/bkg_total);
    }
    roc->GetYaxis()->SetRangeUser(0, max_significance*1.4);
    cout << "Max: " << max_significance << " at " << sig->GetBinCenter(max_index) << " with value: " << sig->GetBinLowEdge(max_index) << endl;
    TCanvas* canvas_eff = new TCanvas("canvas_eff", "canvas", 600, 600);
    SetAtlasStyleCanvas(canvas_eff, true);
    eff_sig->SetFillColor(0);
    eff_sig->Draw("HIST");
    eff_bkg->SetFillColor(0);
    eff_bkg->Draw("HISTsame");
    // canvas_eff->Update();
    int roc_color = kGreen + 1;
    TAxis* xaxis = roc->GetXaxis();
    TAxis* yaxis = roc->GetYaxis();

    roc->SetLineColor(roc_color);
    roc->SetFillColor(0);
    Float_t rightmax = roc->GetMaximum();
    // Float_t scale = gPad->GetUymax()/rightmax;
    Float_t scale = yaxis->GetXmax()/rightmax;
    cout << "scaling factor: "<< scale << endl;
    roc->Scale(scale);
    roc->Draw("HISTsame");
    cout << gPad->GetUxmax() << endl;
    /*
       TGaxis* axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
       gPad->GetUxmax(), gPad->GetUymax(),
       0, rightmax, eff_sig->GetNdivisions("Y"), "+L");
       */
    TGaxis* axis = new TGaxis(xaxis->GetXmax(), yaxis->GetXmin(),
            xaxis->GetXmax(), yaxis->GetXmax(),
            0, rightmax, eff_sig->GetNdivisions("Y"), "+L");
    printf("kGreen + 1: %d\n", roc_color);
    axis->SetLineColor(roc_color);
    axis->SetLabelColor(roc_color);
    axis->Draw();
    float x_offset = 0.6;
    float y_offset = 0.9, t_size = 0.04;
    myLineText(x_offset, y_offset, *eff_sig, "signal");
    myLineText(x_offset, y_offset-t_size, *eff_bkg, "background");
    myLineText(x_offset, y_offset-t_size*2, *roc, "significance");
    myText(0.2, 0.85, 1, "#bf{#it{ATLAS}} Internal");
    canvas_eff->SaveAs(Form("%s.pdf", roc->GetName()));
    return max_index;
}

void print_after_cut(const string& name, TH1F* h1, int cutbin)
{
    int tbins = h1->GetNbinsX();
    printf("%s: %.3f (%.1f%%)\n", name.c_str(), 
            h1->Integral(cutbin, tbins), h1->Integral(cutbin, tbins)*100/h1->Integral());
}

void get_ratio_and_error(float a, float b, float& f, float& error)
{
    f = a/b;
    error = f*sqrt((a+b)/(a*b));
    return;
}

TH1F* generate_th(TH1F* h_s1, const char* hist_name, float low_, float hi_ = 0)
{
    TCanvas* canvas = new TCanvas("canvas_new", "canvas", 600, 600);
    h_s1->Draw();
    if(hi_ == 0){
        hi_ = h_s1->GetBinLowEdge(h_s1->GetNbinsX()+1);
    } 
    cout << "Range: " << low_ << " " << hi_ << endl;
    // TF1* fa = new TF1("fa", "[0]*exp([1]*x)", low_, hi_);
    TF1* fa = new TF1("fa", "[0]*exp([1]*x)");
    fa->SetParameter(0, 1.32351e-03);
    fa->SetParameter(1, -7.64122e-02);

    TFitResultPtr fit_res = h_s1->Fit("fa", "SR","", low_, hi_);
    // cout << "[0]: " << fa->GetParameter(0) << endl;
    // cout << "[1]: " << fa->GetParameter(1) << endl;

    int low_bin = h_s1->FindBin(low_);
    int hi_bin = h_s1->FindBin(hi_);
    TAxis* x_axis = h_s1->GetXaxis();
    RooRealVar met("met","met", x_axis->GetXmin(), x_axis->GetXmax());

    TH1F* h1 = nullptr;
    if(h_s1->GetXaxis()->GetXbins()->GetArray())
    {
        RooBinning* binning = new RooBinning(h_s1->GetNbinsX(), h_s1->GetXaxis()->GetXbins()->GetArray(),"newbins");
        met.setBinning(*binning, "newbins");
        RooAbsPdf* roo_exp = RooFit::bindPdf(fa, met);
        h1 = (TH1F*) roo_exp->createHistogram(Form("%s",hist_name), met, RooFit::Binning(*binning));
        delete binning;
        delete roo_exp;
    } else {
        met.setBins(x_axis->GetNbins(), "newbins");
        RooAbsPdf* roo_exp = RooFit::bindPdf(fa, met);
        h1 = (TH1F*) roo_exp->createHistogram(Form("%s",hist_name), met, RooFit::Binning("newbins"));
    }
    h1->SetName(hist_name);
    h1->Scale(h_s1->Integral(low_bin, hi_bin)/h1->Integral(low_bin, hi_bin));

    for(int i = 1; i <= h_s1->GetNbinsX(); i++){
        if( i <= low_bin || i >= hi_bin ){
            h1->SetBinContent(i, h_s1->GetBinContent(i));
        }
    }

    h1->SetMarkerSize(h_s1->GetMarkerSize());
    h1->SetLineColor(h_s1->GetLineColor());
    h1->SetFillColor(h_s1->GetFillColor());
    // h1->Fit("fa", "","", low_, hi_);
    delete fa;
    delete canvas;
    return h1;
}

void print_correlation(TH2* h2)
{
    if(h2 == NULL) return;
    printf("%s: correlation factor: %.2f\n", h2->GetName(), h2->GetCorrelationFactor());
}

double get_mT(double v1_x, double v1_y,  double v2_x, double v2_y)
{
    TVector2 v1_vec(v1_x, v1_y);
    TVector2 v2_vec(v2_x, v2_y);
    return sqrt(2*v1_vec.Mod()*v2_vec.Mod()*(1-TMath::Cos(v1_vec.DeltaPhi(v2_vec))));
}

void compare_two_hists(TH1* h1, TH1* h2,
        const char* x_title, const char* h1_name, const char* h2_name, bool is_log)
{
    if(!h1 || !h2) return ;
    cout << "binning: "<< h1->GetNbinsX() << " " << h2->GetNbinsX() << endl;
    TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
    double max_y = h1->GetMaximum() > h2->GetMaximum()? h1->GetMaximum(): h2->GetMaximum();
    double min_y = h1->GetMinimum() - 0.05;
    if (is_log) {
        min_y = 1e-6;
        // min_y = 1;
        max_y *= 10;
    } else {
        max_y *= 1.1;
    }
    h1->GetYaxis()->SetRangeUser(min_y, max_y);
    h1->SetXTitle(x_title);
    h1->SetLineColor(2);
    h2->SetLineColor(4);
    h2->SetLineStyle(2);
    // h2->Sumw2();
    // h2->Scale(h1->Integral()/h2->Integral());
    TPad* pad = add_ratio_pad(h1, h2);
    pad->cd();
    if (is_log) pad->SetLogy();
    h1->Draw("HIST");
    h2->Draw("same HIST");

    TLegend* legend = myLegend(0.45, 0.8, 0.80, 0.9);
    legend->AddEntry(h1, h1_name, "L");
    legend->AddEntry(h2, h2_name, "L");
    legend->Draw();
    if (is_log){
        canvas->SaveAs(Form("%s_log.pdf", h1->GetName()));
    } else {
        canvas->SaveAs(Form("%s.pdf", h1->GetName()));
    }
    delete canvas;
}

double GetMinOfHist(TH1* h1){
    int total_bins = h1->GetNbinsX();
    double min_value = 9E9;
    for(int ibin = 1; ibin <= total_bins; ibin++){
        double value = (double) h1->GetBinContent(ibin);
        // cout << ibin << " " << value << endl;
        if(fabs(value) < 1E-6) continue;
        min_value = value < min_value? value: min_value;
    }
    return min_value;
}

void GetMaxMin(const TList& hists, double& max, double& min)
{
    max = -999999;
    min = 999999;
    TIter next(hists.MakeIterator());
    TH1* hist = NULL;
    while ((hist = (TH1*) next())) {
        double ff = hist->GetBinContent(hist->GetMaximumBin());
        double jj = GetMinOfHist(hist);
        max = ff>max?ff:max;
        min = jj<min?jj:min;
    }
    return ;
}

void compare_hists(const TList& histograms, const vector<string>& tags,
        const char* x_title, const char* y_title,
        bool is_log, bool color_me=true, bool add_ratio=true, bool is_log_x=false)
{
    int num_hists = histograms.GetSize();
    if(num_hists < 2) return;

    TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
    TH1* h1 = (TH1*) histograms.At(0);
    h1->SetXTitle(x_title);
    h1->SetYTitle(y_title);
    cout <<"total histograms: " << num_hists << endl;
    TH1* hist = NULL;
    TIter next(histograms.MakeIterator());
    if(color_me){
        int icount = 0;
        while ((hist = (TH1*) next())) {
            hist->SetLineColor(my_color_list[icount]);
            hist->SetMarkerColor(my_color_list[icount]);
            hist->SetMarkerSize(0.0);
            icount ++;
        }
    }
    TPad* pad = NULL;
    if (add_ratio) {
        pad = add_ratio_pad(h1, histograms);
        pad->cd();
        if(is_log) pad->SetLogy();
        if(is_log_x) pad->SetLogx();
    } else {
        if (is_log) canvas->SetLogy();
        if (is_log_x) canvas->SetLogx();
    }

    // set range of y-axis
    double max_y, min_y;
    GetMaxMin(histograms, max_y, min_y);
    if (is_log) {
        if(fabs(min_y) < 1E-6) min_y = 1;
        // max_y *= 5;
        // min_y /= 10;
        max_y = 100;
        min_y = 1;
    } else {
        max_y *= 1.1;
    }
    h1->GetYaxis()->SetRangeUser(min_y, max_y);
    float leg_x = 0.43, leg_y = 0.75;
    TLegend* legend = myLegend(leg_x, leg_y, leg_x+0.2, leg_y+0.045*num_hists);

    next.Reset();
    int icount = 0;
    while ((hist = (TH1*) next())) {
        if(icount == 0){
            hist->Draw("HIST");
            // legend->AddEntry(hist, tags.at(icount).c_str(), "LEP");
            legend->AddEntry(hist, tags.at(icount).c_str(), "L");
        } else {
            hist->Draw("same HIST");
            legend->AddEntry(hist, tags.at(icount).c_str(), "L");
        }
        icount += 1;
    }

    legend->Draw("same");
    if (is_log){
        canvas->SaveAs(Form("%s_log.pdf", h1->GetName()));
        canvas->SaveAs(Form("%s_log.eps", h1->GetName()));
    } else {
        canvas->SaveAs(Form("%s.pdf", h1->GetName()));
        canvas->SaveAs(Form("%s.eps", h1->GetName()));
    }
    delete canvas;
}

TH1* merge_hist_list(const TList& hists)
{
    TIter next(hists.MakeIterator());
    TH1* hist = NULL;
    TH1* res = NULL;
    int icount = 0;
    while ((hist = (TH1*) next())){
        if(icount == 0) {
            res = (TH1*) hist->Clone("my_all");
        } else {
            res->Add(hist);
        }
    }
    return res;
}

TH1* merge_stack(const THStack& stack){
    TList* hists = stack.GetHists();
    return merge_hist_list(*hists);
}

TH1F * DrawOverflow(TH1F *h)
{
    // This function paint the histogram h with an extra bin for overflows
    UInt_t nx    = h->GetNbinsX()+1;
    Double_t *xbins= new Double_t[nx+1];
    for (UInt_t i=0;i<nx;i++)
        xbins[i]=h->GetBinLowEdge(i+1);
    xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx)*0.5;
    // Book a temporary histogram having ab extra bin for overflows
    TH1F *htmp = new TH1F(Form("%swtOF",h->GetName()), h->GetTitle(), nx, xbins);
    // Reset the axis labels
    htmp->SetXTitle(h->GetXaxis()->GetTitle());
    htmp->SetYTitle(h->GetYaxis()->GetTitle());

    for (UInt_t i=1; i<=nx; i++)
        htmp->SetBinContent(i, h->GetBinContent(i));

    // FillStyle and color
    htmp->SetFillStyle(h->GetFillStyle());
    htmp->SetFillColor(h->GetFillColor());
    htmp->SetLineStyle(h->GetFillStyle());
    htmp->SetLineColor(h->GetFillColor());
    return htmp;
}

void clean_TList(TList& list){
    TIter next(list.MakeIterator());
    TObject* obj;
    while ((obj = next())){
        delete obj;
    }
}

void compare_hists_from_files(const vector<string>& file_names,
        const vector<string>& tag_name, bool is_log=true, bool add_ratio=false)
{
    if(file_names.size() < 2) return;
    const string& f1_name = file_names.at(0);
    TFile* f1 = TFile::Open(f1_name.c_str(), "read");

    TIter next(f1->GetListOfKeys());
    TH1* h1 = NULL;
    /* do not use vector<const char*> the memeory will be a mess */
    vector<string>* hist_names = new vector<string>();
    while ((h1 = (TH1*) next())) {
        hist_names->push_back(string(h1->GetName()));
    }
    f1->Close();

    vector<TFile*> file_handles;
    for(int i = 0; i < file_names.size(); i++){
        file_handles.push_back(TFile::Open(file_names.at(i).c_str(), "read"));
    }

    bool shape_only = true;
    bool color_me = true;
    for(int ihist = 0; ihist < hist_names->size(); ihist++){
        const char* hist_name = hist_names->at(ihist).c_str();
        if(!hist_name) continue;

        TList* list = new TList();
        for(int i=0; i < file_handles.size(); i++){
            TH1* h_tmp = (TH1*)(file_handles.at(i)->Get(hist_name));
            if(!h_tmp) {
                cout << hist_name << "does not exist in " << i << endl;
                break;
            }
            h_tmp->SetDirectory(0);
            h_tmp->SetName(Form("%s_%s", hist_name, tag_name.at(i).c_str()));
            if (shape_only) h_tmp->Scale(1./h_tmp->Integral());
            h_tmp->SetMarkerSize(0.10);
            list->Add(h_tmp);
        }

        compare_hists(*list, tag_name, hist_name, " ", is_log, color_me, add_ratio);
        clean_TList(*list);
        delete list;
    }
    delete hist_names;
}

void print_graph(TGraph* gr){
    if (!gr) return;
    int n_points = gr->GetN();
    cout << "N points: " << n_points << endl;
    cout << "Loop over: " <<  gr->GetName() << endl;
    double sum = 0;
    double pre_mass = -1;
    for (int i = 0; i <= n_points; i++)
    {
       Double_t x = -1, y = -1;
       gr->GetPoint(i, x, y);
       // cout << i << " " <<  x << " " <<  y << endl;
       if (x >= 200 && x <= 2000) {
            if (pre_mass < 0) {
                pre_mass = x;
                sum = y;
            } else {
                sum += y;
                pre_mass = x;
            }
       }
    }
    cout <<"Summation: " << sum << endl;
    return;
}

double sum_graph_entries(RooCurve* gr, double low, double hi, int nbins)
{
    // cout << "In summation of graph" << endl;
    double sum = 0;
    double step = (hi - low) / nbins;
    for(int i = 0; i < nbins; i++)
    {
        double begin =  low + i*step;
        double end = begin + step;
        double average = gr->average(begin, end);
        sum += average;
        // cout << begin << " " << end << " " << average << " " << endl;
    }
    return sum;
}

#endif
