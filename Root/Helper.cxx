#include "MyXAODTools/Helper.h"

#include <string>
#include <fstream>
#include <iostream>

#include <TString.h>
#include <TVector2.h>

using namespace MyXAODTools;
using namespace std;

TChain* Helper::loader(const char* inFile_name, const char* chain_name)
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

float MyXAODTools::delta_r(float eta1, float phi1, float eta2, float phi2)
{
    double dphi = fabs(TVector2::Phi_mpi_pi( TVector2::Phi_mpi_pi(phi1) - TVector2::Phi_mpi_pi(phi2) ));
    double deta = eta1 - eta2;
    return sqrt(deta * deta + dphi * dphi);
}
