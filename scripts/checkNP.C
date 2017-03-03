#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <exception>
#include <sys/types.h>
#include <boost/algorithm/string.hpp>
#include <dirent.h>
#include <math.h>
#include <errno.h>

#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TText.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TVectorD.h"

bool readConfig(const char* filename, char delim, map<string, map<string, string> >& all_dic);
void tokenizeString(const string& str, char delim, vector<string>& tokens);


void checkNP(const string& configPath, TString& pathToFiles, int n_NPs){

 // # of NP
 // int n = 27; 
 int n = n_NPs; 
 string* NPlabes = new string[n];

// TString pathToFiles = "/afs/cern.ch/user/x/xju/work/h4l/highmass/systematics/syst_inputs20170303/HighMass_VBF/bkg/";
// string configPath = "../config_qqZZ.ini";

 // read config file
 map<string, map<string, string> > p_dic;
 // string configPath = "../config.ini";
 readConfig(configPath.c_str(),'=', p_dic);

 // get norm files list
 vector<string> normFilesList, sampleList; 
 tokenizeString(p_dic["main"]["samples"].c_str(), ',', normFilesList);
 for(int i=0; i<(int)normFilesList.size(); ++i) {  
    sampleList.push_back(normFilesList[i]);
    normFilesList[i] = pathToFiles + "norm_" + normFilesList[i] + ".txt"; 
 }

// loop over samples
for(int i=0; i<(int)normFilesList.size(); ++i){

  typedef map<string, string>::iterator it_type;
  map<string, map<string, string> > NP_dic;
  readConfig(normFilesList[i].c_str(),'=', NP_dic);

  vector<string> catList; 
  tokenizeString(p_dic["main"]["categories"].c_str(), ',', catList);

  for(int j=0; j<(int)catList.size(); ++j){

      float xl[n], xh[n];
      
      cout << catList[j] <<endl;
      int count(0);
      float max_sys =  0.;
      for(it_type iterator = NP_dic[catList[j].c_str()].begin(); iterator != NP_dic[catList[j].c_str()].end(); iterator++) 
      {
    
         istringstream iss(iterator->second);
         float down, up;
         iss >> down;
         iss >> up;

         std::string t = iterator->first;
         std::string s = "ATLAS_";

         std::string::size_type is = t.find(s);

         if (is != std::string::npos)
             t.erase(is, s.length());
 
         NPlabes[count] = t;

         xl[count] = (down-1.)*100.;
         xh[count] = (up-1.)*100.;
         cout <<  iterator->first << " = " << (down-1.)*100. << "/" << (up-1.)*100. <<endl;
         if (abs(xl[count]) > max_sys){
             max_sys = abs(xl[count]);
         }
         if (abs(xh[count]) > max_sys){
             max_sys = abs(xh[count]);
         }

         count = count + 1;   
       }
   string hist_downName = catList[j] + "down"; 
   string hist_upName   = catList[j] + "up"; 

   gStyle->SetHistMinimumZero();
   cout << "max before ceil: " << max_sys << endl;
   max_sys = ceil(max_sys);
   cout << "max after ceil: " << max_sys << endl;
   float min_sys = -1*max_sys;

   TH1F *hist_down = new TH1F(hist_downName.c_str(),hist_downName.c_str(),n,0,n);
    hist_down->SetFillColor(4);
    hist_down->SetBarWidth(0.5);
    hist_down->SetStats(0);
    hist_down->SetMinimum(min_sys);
    hist_down->SetMaximum(max_sys);
   TH1F *hist_up   = new TH1F(hist_upName.c_str(),hist_upName.c_str(),n,0,n);
    hist_up->SetFillColor(38);
    hist_up->SetBarWidth(0.5);

   for (int t=1; t<=n; t++) {
      hist_down->Fill(NPlabes[t-1].c_str(), xl[t-1]);
      hist_down->GetXaxis()->SetBinLabel(t,NPlabes[t-1].c_str());
      hist_up->Fill(NPlabes[t-1].c_str(), xh[t-1]);
     }

   TString cName = sampleList[i]+"_"+catList[j];
   TCanvas* cv = new TCanvas(cName,cName,450,800);
    gPad->SetLeftMargin(0.5);
    gPad->SetBottomMargin(0.1);
    hist_down->Draw("hbar");
    hist_up->Draw("hbar same");
    hist_down->GetXaxis()->SetLabelSize(0.03);
    hist_down->GetYaxis()->SetTitle("Down/up [%]");
    hist_down->GetYaxis()->SetTitleOffset(1);
    hist_down->GetYaxis()->CenterTitle();
   TLegend *leg=new TLegend(.8214,.9586,.9911,.9871);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->AddEntry(hist_down, "Down", "F");
    leg->AddEntry(hist_up,   "Up", "F");
    leg->Draw(); 
   TLatex * lm = new TLatex(-10., 77., catList[j].c_str()); 
    lm->SetTextSize(0.035);
    lm->Draw();
   
   TString fileN = sampleList[i]+"_"+catList[j]+".eps";

   cv->Print(fileN);

   }// end cat

} //end samples

}

bool readConfig(const char* filename, char delim,
        map<string, map<string, string> >& all_dic)
{
    cout << "Reading: " << filename << endl;
    ifstream file(filename, ifstream::in);
    if(file.fail() || file.bad()){
        cout <<"File: "<< filename << " is bad." << endl; 
        return false;
    }
    string line;
    int lineCount = 0;
    map<string, string> section_dic;
    string section_name;
    while ( getline(file, line) ){ 
        boost::algorithm::trim(line);
        if( line[0] == '[' ){
            if( lineCount < 1 ){
                section_name = string(line.begin()+1, line.end()-1); 
            }else{
                all_dic[section_name] = section_dic;
                section_dic.clear();
                section_name = string(line.begin()+1, line.end()-1);
            }
        }else if( line[0] == '#' ){
            continue;
        }else{
            size_t delim_pos = line.find(delim);
            if (delim_pos != string::npos) {
                string tagName = line.substr(0, delim_pos-1);
                string token = line.substr(delim_pos+1, line.size());
                boost::algorithm::trim(tagName);
                boost::algorithm::trim(token);
                section_dic[tagName] = token;
            } 
            else {
                cerr << line << " Does not have delimeter '" << delim << "', ignored" << endl;
            }
        }
        lineCount ++ ;
    }
    all_dic[section_name] = section_dic;  //pick up the last section
    file.close();
    return true;
}


void tokenizeString(const string& str, char delim, vector<string>& tokens)
{
    tokens.clear();
    istringstream iss(str);
    string token;
    while ( getline(iss, token, delim) ){
        boost::algorithm::trim(token);
        tokens.push_back(token);
    }
}
