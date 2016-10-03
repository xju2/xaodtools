/*
  H4l/MultiLepton 2012 MENU
==========================================================================================================
  double Et_cl           =  el_cl_E/cosh(el_etas2);
  double etas2           =  el_etas2;
  double rHad            =  el_Ethad/Et_cl;
  double rHad1           =  el_Ethad1/Et_cl;
  double Reta            =  el_reta;
  double w2              =  el_weta2;
  double f1              =  el_f1;
  double f3              =  el_f3;
  double wstot           =  el_wstot;
  double DEmaxs1         = (el_emaxs1-el_Emax2)/(el_emaxs1+el_Emax2);
  double deltaEta        = el_deltaeta1;
  double deltaPhiRescaled = el_deltaphiRescaled;
  int    nSi              = el_nSiHits;
  int    nSiDeadSensors   = el_nSCTDeadSensors+el_nPixelDeadSensors;
  int    nPix             = el_nPixHits;
  int    nPixDeadSensors  = el_nPixelDeadSensors;
  int    nTRThigh         = el_nTRTHighTHits;
  int    nTRThighOutliers = el_nTRTHighTOutliers;
  int    nTRT             = el_nTRTHits;  
  int    nTRTOutliers     = el_nTRTOutliers; 
  bool   expectBlayer     = el_expectBLayerHit; 
  int    nBlayerHits      = el_nBLHits;

  double rTRT = 0.;
  rTRT = (nTRT+nTRTOutliers) > 0 ?  ((double) (nTRThigh+nTRThighOutliers)/(nTRT+nTRTOutliers) ) : 0.;
  int nTRTTotal = nTRT+nTRTOutliers;

  double dpOverp =0;
  for (unsigned int i = 0; i<el_refittedTrack_LMqoverp.size;++i)
    {
      if((el_refittedTrack_author).at(i)==4)
	{
	  dpOverp= 1-(el_trackqoverp/(el_refittedTrack_LMqoverp.at(i)));
	}
    }
  return passMultiLepton(etas2, Et_cl,
                            rHad, rHad1, Reta, w2,
                            f1,f3, wstot, DEmaxs1, 
			    deltaEta, nSi,nSiDeadSensors, nPix, 
			    nPixDeadSensors,deltaPhiRescaled,dpOverp,
			    rTRT,nTRTTotal,nBlayerHits,expectBlayer);
===============================================================================================
*/

// #include "MultiLeptonDefs_new.h"
#include <iostream>
#include <cmath>

namespace{
  const double GeV = 1000;
}
bool passMultiLepton(double eta, double eT,
		     double rHad, double rHad1, double Reta, double w2, 
		     double f1, double f3, double wstot, double DEmaxs1, 
		     double deltaEta, int nSi, int nSiDeadSensors, int nPix, 
		     int nPixDeadSensors, double deltaPhiRes, double dpOverp, 
		     double TRratio,  int nTRTTotal,  int nBlayerHits, bool expectBlayer,bool debug=false ); 
//Helper Functions 
bool passRHad(double rhad, double rhad1, unsigned int etbin, unsigned int etabin);
bool passF3(double f3, unsigned int etbin, unsigned int etabin);
bool passReta(double reta, unsigned int etbin, unsigned int etabin);
bool passW2(double w2, unsigned int etbin, unsigned int etabin);
bool passWstot(double wstot, unsigned int etbin, unsigned int etabin);
bool passEratio(double demaxs1, unsigned int etbin, unsigned int etabin);
bool passDeltaEta(double deltaEta, unsigned int etbin, unsigned int etabin);
bool passDeltaPhiRes(double deltaPhiRes, bool isBrem, unsigned int etbin, unsigned int etabin);
bool passTR(double TRratio, double eta, unsigned int  nTRTTotal );
bool passConversion(double deltaPhiRes, bool expectBlayer,  int nBlayerHits ,unsigned int eTBin, unsigned int etaBin);
// Helper Fuctions
unsigned int getEtaBinH4l(double eta);
unsigned int getEtBinH4l(double eT);
bool getBremCategoryH4l(double dpOverp, unsigned int etbin, unsigned int etabin);


//----------------------------------------------------------------------------------------
bool passMultiLepton(double eta, double eT,
		     double rHad, double rHad1, double Reta, double w2, 
		     double f1, double f3, double wstot, double DEmaxs1, 
		     double deltaEta, int nSi, int nSiDeadSensors, int nPix, 
		     int nPixDeadSensors, double deltaPhiRes, double dpOverp, 
		     double TRratio, int nTRTTotal,int nBlayerHits, bool expectBlayer,
		     bool debug  
			) {
         		      
  // Sanity Check on eta
  if( fabs( eta ) > 2.47 ) return false;
  
  //Get eta/et bins
  unsigned int eTBin = getEtBinH4l(eT);
  unsigned int etaBin = getEtaBinH4l(eta);
  
  //High Low Brem
  bool isBrem = getBremCategoryH4l(dpOverp, eTBin, etaBin);

  // RHad
  if(!passRHad(rHad,rHad1,eTBin,etaBin)){
    if(debug){
      std::cout << "Failed RHad " << rHad << " " << eT << " " << eta << " " << std::endl;
    }
    return false;
  }
  // f3
  if(eT<90*GeV){
    if(!passF3(f3,eTBin,etaBin)){
      if(debug){
	std::cout << "Failed f3 " << f3 << " " << eT << " " << eta << " " << std::endl;
      }
      return false;
    }
  }
  // Reta 
  if(!passReta(Reta,eTBin,etaBin)){
    if(debug){
      std::cout << "Failed Reta " << Reta << " " << eT << " " << eta << " " << std::endl;
    }
    return false;
  }
  // w2
  if(!passW2(w2,eTBin,etaBin)){
    if(debug){
      std::cout << "Failed w2 " << w2 << " " << eT << " " << eta << " " << std::endl;
    }
    return false;
  }
  // Check the energy in the strips before cutting on it
  if( f1 > 0.005 ){
    if(!passWstot(wstot,eTBin,etaBin)){
      if(debug)
	std::cout << "Failed wstot " << wstot << " " << eT << " " << eta << " " << std::endl;
      return false;
    }
    // Eratio
    if(!passEratio(DEmaxs1,eTBin,etaBin)){
      if(debug){
	std::cout << "Failed DEmaxs1 " << DEmaxs1 << " " << eT << " " << eta << " " << std::endl;
      }
      return false;
    }
  }
  
  // Delta Eta
  if(!passDeltaEta(deltaEta, eTBin, etaBin)){
    if(debug){
      std::cout << "Failed dEta " << deltaEta << " " << eT << " " << eta << " " << std::endl;
    }
    return false;
  }  
  // Rescale deltaPhi
  /***
  if(!passDeltaPhiRes(deltaPhiRes, isBrem, eTBin, etaBin)){
      if(debug){
	std::cout << "Failed dPhiRes " << deltaPhiRes << " " << isBrem << " " << eT << " " << eta << " " << std::endl;
      }
      return false;
  }
  ***/
  
  // Si
  if((nSi + nSiDeadSensors)  < 7){
    if(debug){
      std::cout << "Failed nSi " << nSi + nSiDeadSensors << " " << eT << " " << eta << " " << std::endl;
    }
    return false;
  }
 
   // Pix  & B-layer
   // Pix  
//   if((nPix+nPixDeadSensors) < 1){
//     if(debug){
//       std::cout << "Failed nPix " << nPix << " " << eT << " " << eta << " " << std::endl;
//     }
//     return false;
//   }
 /*** 
   if (expectBlayer) {
 	  if (nBlayerHits < 1 || (nPix+nPixDeadSensors) < 2) {    
 		  if(debug){
 			std::cout << "Failed nPix " << nPix << " " << eT << " " << eta << " " << std::endl;
 		  }
 		  return false;
 	  }
   }
   else {
     if((nPix+nPixDeadSensors) < 1){
 		if(debug){
 		  std::cout << "Failed nPix " << nPix << " " << eT << " " << eta << " " << std::endl;
 		}
 		return false;
 	}
   }
   ***/

  //TRT Ratio in crack
  if (fabs(eta)>1.37 && fabs(eta)<1.52 ){
    if(!passTR(TRratio,eta,nTRTTotal)){
      if(debug) std::cout << "Failed TR " << TRratio << " "<<nTRTTotal<<" "<< eT << " " << eta << " " << std::endl;
      return false;
    }
  }

  //Conversion cut 
  /****
  if (!passConversion(deltaPhiRes,  expectBlayer, nBlayerHits ,eTBin, etaBin)){  
    if(debug) std::cout << "Failed conversions " << expectBlayer << " "<<nBlayerHits<<" "<< deltaPhiRes <<" "<< eT << " " << eta << " " << std::endl;
    return false;
  }
  ***/

  return true;
}

//---------------------------------------------------------------------------------------
// Gets the Eta bin [0-9] given the eta
unsigned int getEtaBinH4l(double eta){

  const unsigned int nEtaBins = 10;
  static const double etaBins[nEtaBins] = {0.1,0.6,0.8,1.15,1.37,1.52,1.81,2.01,2.37,2.47};
  
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin){
    if(fabs(eta) < etaBins[etaBin])
      return etaBin;
  }
  return 9;
}

//---------------------------------------------------------------------------------------
// Gets the Et bin [0-10] given the et (MeV)
unsigned int getEtBinH4l(double eT){
  const unsigned int nEtBins = 10;
  static const double eTBins[nEtBins] = {5*GeV,10*GeV,15*GeV,20*GeV,30*GeV,40*GeV,50*GeV,60*GeV,70*GeV,80*GeV};
  
  for(unsigned int eTBin = 0; eTBin < nEtBins; ++eTBin){
    if(eT < eTBins[eTBin])
      return eTBin;
  }
  
  return 10;
}

//----------------------------------------------------------------------------------------

// Determine whether electron is high- or low-brem using dp/p
bool getBremCategoryH4l(double dpOverp, unsigned int eTBin, unsigned int etaBin){

    //                                      0.0    0.1    0.6    0.8    1.15   1.37   1.52   1.81   2.01   2.37    2.47

 static const double bremThresholds[11][10] = {{0.094, 0.094, 0.142, 0.218, 0.222, 0.358, 0.402, 0.530, 0.298, 0.358}, // 0-5 GeV
					       {0.102, 0.102, 0.162, 0.322, 0.326, 0.370, 0.506, 0.614, 0.314, 0.386}, // 5-10 GeV
					       {0.122, 0.122, 0.198, 0.442, 0.442, 0.486, 0.638, 0.690, 0.338, 0.442}, // 10-15 GeV
					       {0.130, 0.130, 0.214, 0.514, 0.558, 0.530, 0.718, 0.738, 0.334, 0.462}, // 15-20 GeV
					       {0.158, 0.158, 0.278, 0.574, 0.546, 0.602, 0.778, 0.742, 0.350, 0.506}, // 20-30 GeV
					       {0.342, 0.342, 0.474, 0.730, 0.762, 0.798, 0.854, 0.838, 0.386, 0.526}, // 30-40 GeV
					       {0.342, 0.342, 0.474, 0.730, 0.762, 0.798, 0.854, 0.838, 0.386, 0.526}, // 40-50 GeV
					       {0.342, 0.342, 0.474, 0.730, 0.762, 0.798, 0.854, 0.838, 0.386, 0.526}, // 50-60 GeV
					       {0.342, 0.342, 0.474, 0.730, 0.762, 0.798, 0.854, 0.838, 0.386, 0.526}, // 60-70 GeV
					       {0.342, 0.342, 0.474, 0.730, 0.762, 0.798, 0.854, 0.838, 0.386, 0.526}, // 70-80 GeV
					       {0.342, 0.342, 0.474, 0.730, 0.762, 0.798, 0.854, 0.838, 0.386, 0.526}};// >80 GeV
	
  if(dpOverp < bremThresholds[eTBin][etaBin]) return false; // Low-brem electron
  return true; // High-brem electron
}

//----------------------------------------------------------------------------------------
bool passRHad(double rhad, double rhad1,unsigned int etbin,unsigned int etabin){

  //                           0.0     0.1      0.6       0.8        1.15      1.37     1.52      1.81      2.01      2.37     2.47
  static const double cutrhad[11][10] = {{ 0.258,  0.258, 0.208,  0.206,  0.256,  0.188,  0.254,  0.254,  0.226,  0.176 } // 0 - 5 GeV
					,{ 0.155,  0.155, 0.115,  0.125,  0.135,  0.100,  0.140,  0.135,  0.125,  0.105 } // 5 - 10 GeV
					,{ 0.075,  0.075, 0.065,  0.065,  0.065,  0.055,  0.080,  0.080,  0.065,  0.050 } // 10 - 15 GeV
					,{ 0.055,  0.055, 0.045,  0.045,  0.045,  0.040,  0.055,  0.055,  0.050,  0.040 } // 15 - 20 GeV
					,{ 0.038,  0.038, 0.032,  0.032,  0.032,  0.026,  0.040,  0.040,  0.035,  0.030 } // 20 - 30 GeV
					,{ 0.025,  0.025, 0.022,  0.022,  0.022,  0.019,  0.030,  0.030,  0.025,  0.020 } // 30 - 40 GeV
					,{ 0.025,  0.025, 0.021,  0.022,  0.021,  0.019,  0.028,  0.028,  0.023,  0.020 } // 40 - 50 GeV
					,{ 0.025,  0.025, 0.021,  0.022,  0.021,  0.019,  0.028,  0.028,  0.023,  0.020 } // 50 - 60 GeV
					,{ 0.025,  0.025, 0.021,  0.022,  0.021,  0.019,  0.028,  0.028,  0.023,  0.020 } // 60 - 70 GeV
					,{ 0.025,  0.025, 0.021,  0.022,  0.021,  0.019,  0.028,  0.028,  0.023,  0.020 } // 70 - 80 GeV
					,{ 0.025,  0.025, 0.021,  0.022,  0.021,  0.019,  0.028,  0.028,  0.023,  0.020 } }; // 80 - inf GeV

  if(etabin == 3 || etabin == 4){
    if (rhad > cutrhad[etbin][etabin])
      return false;
  } else{
    if(rhad1 > cutrhad[etbin][etabin])
      return false;
  }
  return true;
}

//----------------------------------------------------------------------------------------
bool passF3(double f3,unsigned int etbin,unsigned int etabin){

  //                           0.0     0.1      0.6       0.8        1.15      1.37     1.52      1.81      2.01      2.37     2.47
 static const double cutf3[11][10] = {{  0.0593 , 0.0593 , 0.0540 , 0.0430 , 0.0481 ,  9999 , 0.0363 , 0.0391 , 0.0420 , 9999 }     // 0 - 5 GeV
				      ,{ 0.0377 , 0.0377 , 0.0356 , 0.0327 , 0.0365 ,  9999 , 0.0275 , 0.0315 , 0.0360 , 9999 } // 5 - 10 GeV
				      ,{ 0.0279 , 0.0279 , 0.0261 , 0.0239 , 0.0267 ,  9999 , 0.0217 , 0.0261 , 0.0270 , 9999 } // 10 - 15 GeV
				      ,{ 0.0259 , 0.0259 , 0.0219 , 0.0211 , 0.0239 ,  9999 , 0.0203 , 0.0253 , 0.0270 , 9999 } // 15 - 20 GeV
				      ,{ 0.0252 , 0.0252 , 0.0199 , 0.0196 , 0.0225 ,  9999 , 0.0207 , 0.0261 , 0.0270 , 9999 } // 20 - 30 GeV
				      ,{ 0.0259 , 0.0259 , 0.0197 , 0.0193 , 0.0215 ,  9999 , 0.0223 , 0.0274 , 0.0270 , 9999 } // 30 - 40 GeV
				      ,{ 0.0265 , 0.0265 , 0.0201 , 0.0201 , 0.0222 ,  9999 , 0.0240 , 0.0291 , 0.0290 , 9999 } // 40 - 50 GeV
				      ,{ 0.0265 , 0.0265 , 0.0201 , 0.0201 , 0.0222 ,  9999 , 0.0240 , 0.0291 , 0.0290 , 9999 } // 50 - 60 GeV
				      ,{ 0.0292 , 0.0292 , 0.0219 , 0.0215 , 0.0241 ,  9999 , 0.0264 , 0.0327 , 0.0315 , 9999 } // 60 - 70 GeV
				      ,{ 0.0292 , 0.0292 , 0.0219 , 0.0215 , 0.0241 ,  9999 , 0.0264 , 0.0327 , 0.0315 , 9999 } // 70 - 80 GeV
				      ,{  9999  ,  9999  ,  9999  ,  9999  ,  9999  ,  9999 ,  9999  ,  9999  ,  9999  , 9999 } }; // 80 - inf GeV

  if(f3 > cutf3[etbin][etabin]) return false;
  
  return true;
}

//----------------------------------------------------------------------------------------
bool passReta(double rEta, unsigned int eTBin, unsigned int etaBin){
 
  //				     		0.0     0.1      0.6    0.8        1.15    1.37     1.52      1.81      2.01      2.37     2.47
  static const double cutReta37[11][10] = {{     0.581 ,  0.581 ,  0.526 , 0.600  , 0.650  ,  0.650 ,  0.750 ,  0.639 ,  0.690 ,  0.690 }    // 0 - 5 GeV  
                                            , {  0.748 ,  0.748 ,  0.705 , 0.700  , 0.750  ,  0.700 ,  0.750 ,  0.732 ,  0.775 ,  0.775 }    // 5 - 10 GeV 
                                            , {  0.827 ,  0.827 ,  0.804 , 0.800  , 0.850  ,  0.750 ,  0.800 ,  0.802 ,  0.830 ,  0.830 }    // 10 - 15 GeV
                                            , {  0.863 ,  0.863 ,  0.845 , 0.826  , 0.850  ,  0.750 ,  0.850 ,  0.847 ,  0.853 ,  0.853 }    // 15 - 20 GeV
                                            , {  0.893 ,  0.893 ,  0.878 , 0.864  , 0.850  ,  0.800 ,  0.850 ,  0.873 ,  0.879 ,  0.878 }    // 20 - 30 GeV
                                            , {  0.893 ,  0.893 ,  0.878 , 0.864  , 0.850  ,  0.800 ,  0.900 ,  0.873 ,  0.879 ,  0.878 }    // 30 - 40 GeV
                                            , {  0.917 ,  0.917 ,  0.908 , 0.920  , 0.910  ,  0.800 ,  0.900 ,  0.898 ,  0.898 ,  0.896 }    // 40 - 50 GeV
                                            , {  0.917 ,  0.917 ,  0.908 , 0.920  , 0.910  ,  0.800 ,  0.900 ,  0.898 ,  0.898 ,  0.896 }    // 50 - 60 GeV
                                            , {  0.917 ,  0.917 ,  0.908 , 0.920  , 0.910  ,  0.800 ,  0.900 ,  0.898 ,  0.898 ,  0.896 }    // 60 - 70 GeV
                                            , {  0.917 ,  0.917 ,  0.908 , 0.920  , 0.910  ,  0.800 ,  0.900 ,  0.898 ,  0.898 ,  0.896 }    // 70 - 80 GeV
                                            , {  0.917 ,  0.917 ,  0.908 , 0.920  , 0.910  ,  0.800 ,  0.900 ,  0.898 ,  0.898 ,  0.896 } }; // 80 - inf GeV  
  
  if(rEta < cutReta37[eTBin][etaBin]) return false;

  return true;
}

//----------------------------------------------------------------------------------------

bool passW2(double w2, unsigned int eTBin, unsigned int etaBin){

  //                                  0.0     0.1      0.6       0.8       1.15      1.37      1.52      1.81      2.01      2.37     2.47
  static const double cutWeta2[11][10] = { {   0.0166 , 0.0166 , 0.0172 , 0.0167 , 0.0170 , 0.0385 , 0.0164 , 0.0152 , 0.0156 , 0.0157 }    // 0 - 5 GeV
					   , { 0.0145 , 0.0145 , 0.0152 , 0.0154 , 0.0158 , 0.0347 , 0.0159 , 0.0140 , 0.0150 , 0.0150 }    // 5 - 10 GeV
					   , { 0.0129 , 0.0129 , 0.0137 , 0.0141 , 0.0146 , 0.0311 , 0.0151 , 0.0133 , 0.0140 , 0.0140 }    // 10 - 15 GeV
					   , { 0.0122 , 0.0122 , 0.0129 , 0.0133 , 0.0139 , 0.0278 , 0.0145 , 0.0128 , 0.0140 , 0.0140 }    // 15 - 20 GeV
					   , { 0.0117 , 0.0117 , 0.0123 , 0.0126 , 0.0131 , 0.0257 , 0.0139 , 0.0124 , 0.0135 , 0.0135 }    // 20 - 30 GeV
					   , { 0.0117 , 0.0117 , 0.0123 , 0.0126 , 0.0131 , 0.0257 , 0.0139 , 0.0124 , 0.0135 , 0.0135 }    // 30 - 40 GeV
					   , { 0.0112 , 0.0112 , 0.0118 , 0.0121 , 0.0125 , 0.0247 , 0.0132 , 0.0120 , 0.0130 , 0.0135 }    // 40 - 50 GeV
					   , { 0.0112 , 0.0112 , 0.0118 , 0.0121 , 0.0125 , 0.0247 , 0.0132 , 0.0120 , 0.0130 , 0.0135 }    // 50 - 60 GeV
					   , { 0.0112 , 0.0112 , 0.0118 , 0.0121 , 0.0125 , 0.0247 , 0.0132 , 0.0120 , 0.0130 , 0.0135 }    // 60 - 70 GeV
					   , { 0.0112 , 0.0112 , 0.0118 , 0.0121 , 0.0125 , 0.0247 , 0.0132 , 0.0120 , 0.0130 , 0.0135 }    // 70 - 80 GeV
					   , { 0.0112 , 0.0112 , 0.0118 , 0.0121 , 0.0125 , 0.0247 , 0.0132 , 0.0120 , 0.0130 , 0.0135 } }; // 80 - inf GeV
  

  if(w2 > cutWeta2[eTBin][etaBin]) return false;
    
  return  true;
}

//----------------------------------------------------------------------------------------
bool passWstot(double wstot, unsigned int eTBin, unsigned int etaBin){  

  //                                  0.0     0.1      0.6       0.8       1.15      1.37      1.52      1.81      2.01      2.37     2.47
  static const double cutWstot[11][10] = { {  3.926,   3.926,   4.069,   4.501,   4.986,   9999,   4.650,   3.190,   1.966,   9999 }    // 0 - 5 GeV
					  , { 3.296,   3.296,   3.427,   3.936,   4.309,   9999,   4.313,   2.845,   1.818,   9999 }    // 5 - 10 GeV
					  , { 3.095,   3.095,   3.202,   3.708,   4.095,   9999,   3.968,   2.692,   1.754,   9999 }    // 10 - 15 GeV
					  , { 3.035,   3.035,   3.129,   3.553,   3.941,   9999,   3.758,   2.555,   1.714,   9999 }    // 15 - 20 GeV
					  , { 3.035,   3.035,   3.129,   3.508,   3.793,   9999,   3.609,   2.505,   1.703,   9999 }    // 20 - 30 GeV
					  , { 2.881,   2.881,   2.941,   3.319,   3.506,   9999,   3.380,   2.381,   1.644,   9999 }    // 30 - 40 GeV
					  , { 2.881,   2.881,   2.941,   3.319,   3.506,   9999,   3.380,   2.381,   1.644,   9999 }    // 40 - 50 GeV
					  , { 2.881,   2.881,   2.941,   3.319,   3.506,   9999,   3.380,   2.381,   1.644,   9999 }    // 50 - 60 GeV
					  , { 2.881,   2.881,   2.941,   3.319,   3.506,   9999,   3.380,   2.381,   1.644,   9999 }    // 60 - 70 GeV
					  , { 2.881,   2.881,   2.941,   3.319,   3.506,   9999,   3.380,   2.381,   1.644,   9999 }    // 70 - 80 GeV
					  , { 2.881,   2.881,   2.941,   3.319,   3.506,   9999,   3.380,   2.381,   1.644,   9999 } }; // 80 - inf GeV
  
  
  if(wstot > cutWstot[eTBin][etaBin]) return false;
    
  return  true;
}
//----------------------------------------------------------------------------------------
bool passEratio(double DEmaxs1, unsigned int eTBin, unsigned int etaBin){  
  
  //                                       0.0     0.1      0.6      0.8      1.15      1.37    1.52     1.81    2.01    2.37    2.47
  static const double cutDEmaxs1[11][10] = { {    0.278 ,  0.278 ,  0.122 , 0.150 ,  0.250 , -9999 ,  0.100  ,  0.136 ,  0.492 , -9999 }    // 0 - 5 GeV
					     , {  0.506 ,  0.506 ,  0.320 , 0.300 ,  0.250 , -9999 ,  0.150  ,  0.196 ,  0.543 , -9999 }    // 5 - 10 GeV
					     , {  0.587 ,  0.587 ,  0.509 , 0.500 ,  0.350 , -9999 ,  0.250  ,  0.281 ,  0.616 , -9999 }    // 10 - 15 GeV
					     , {  0.591 ,  0.591 ,  0.556 , 0.550 ,  0.350 , -9999 ,  0.350  ,  0.369 ,  0.639 , -9999 }    // 15 - 20 GeV
					     , {  0.627 ,  0.627 ,  0.617 , 0.610 ,  0.600 , -9999 ,  0.600  ,  0.505 ,  0.653 , -9999 }    // 20 - 30 GeV
					     , {  0.627 ,  0.627 ,  0.617 , 0.610 ,  0.600 , -9999 ,  0.600  ,  0.505 ,  0.653 , -9999 }    // 30 - 40 GeV
					     , {  0.627 ,  0.627 ,  0.617 , 0.610 ,  0.600 , -9999 ,  0.600  ,  0.505 ,  0.653 , -9999 }    // 40 - 50 GeV
					     , {  0.627 ,  0.627 ,  0.617 , 0.610 ,  0.600 , -9999 ,  0.600  ,  0.505 ,  0.653 , -9999 }    // 50 - 60 GeV
					     , {  0.627 ,  0.627 ,  0.617 , 0.610 ,  0.600 , -9999 ,  0.600  ,  0.505 ,  0.653 , -9999 }    // 60 - 70 GeV
					     , {  0.627 ,  0.627 ,  0.617 , 0.610 ,  0.600 , -9999 ,  0.600  ,  0.505 ,  0.653 , -9999 }    // 70 - 80 GeV
					     , {  0.627 ,  0.627 ,  0.617 , 0.610 ,  0.600 , -9999 ,  0.600  ,  0.505 ,  0.653 , -9999 } }; // 80  GeV
  
  if(DEmaxs1 < cutDEmaxs1[eTBin][etaBin]) return false;
  
  return  true;
}

//----------------------------------------------------------------------------------------
bool passDeltaEta(double deltaEta, unsigned int eTBin, unsigned int etaBin){
  
  //                                         0.0     0.1      0.6     0.8      1.15     1.37     1.52     1.81     2.01    2.37     2.47
  static const double cutDeltaEta[11][10]  = {{0.020,   0.020,   0.020,   0.020,   0.020,   0.018,   0.020,   0.020,   0.020,   0.020}, // 0
					      {0.012,   0.009,   0.009,   0.011,   0.015,   0.016,   0.014,   0.011,   0.015,   0.018}, // 5
					      {0.012,   0.009,   0.008,   0.008,   0.008,   0.013,    0.01,   0.007,    0.01,   0.014}, // 10
					      {0.011,   0.008,   0.007,   0.007,   0.007,   0.012,   0.007,   0.007,   0.008,   0.012}, // 15
					      {0.011,   0.008,   0.007,   0.007,   0.009,   0.011,   0.009,   0.007,   0.007,    0.01}, // 20
					      {0.011,   0.008,   0.007,   0.007,   0.007,   0.009,   0.007,   0.007,   0.007,   0.008}, // 30
					      {0.011,   0.008,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007}, // 40
					      {0.011,   0.008,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007}, // 50
					      {0.011,   0.008,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007}, // 60
					      {0.011,   0.008,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007}, // 70
					      {0.011,   0.008,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007,   0.007}}; // 80
  
  if (fabs(deltaEta) > cutDeltaEta[eTBin][etaBin]) return false;  
  return true;

}

//----------------------------------------------------------------------------------------
bool passDeltaPhiRes(double deltaPhiRes, bool isBrem, unsigned int eTBin, unsigned int etaBin){
  
            //                               0.0            0.1      0.6       0.8     1.15      1.37     1.52      1.81     2.01    2.37     2.47   
  static const double maxDeltaPhiResHigh[11][10]  = {{0.085,   0.085,   0.073,   0.052,   0.045,   0.050,   0.049,    0.04,   0.034,   0.046}, // 0
						     {0.040,   0.040,   0.034,   0.031,   0.025,   0.030,   0.028,   0.025,   0.028,   0.043}, // 5
						     {0.022,   0.022,   0.022,   0.019,   0.015,   0.025,   0.022,   0.022,   0.028,   0.031}, // 10
						     {0.019,   0.019,   0.019,   0.016,   0.016,   0.019,   0.016,   0.019,   0.025,   0.031}, // 15
						     {0.019,   0.019,   0.019,   0.016,   0.016,   0.015,   0.016,   0.022,   0.025,   0.028}, // 20
						     {0.019,   0.019,   0.019,    0.01,   0.013,   0.013,   0.016,   0.019,   0.025,   0.025}, // 30
						     {0.019,   0.019,   0.019,    0.01,    0.01,   0.013,   0.013,   0.022,   0.022,   0.025}, // 40
						     {0.019,   0.019,   0.019,    0.01,    0.01,   0.013,   0.013,   0.022,   0.022,   0.025}, // 50
						     {0.019,   0.019,   0.019,    0.01,    0.01,   0.013,   0.013,   0.022,   0.022,   0.025}, // 60
						     {0.019,   0.019,   0.019,    0.01,    0.01,   0.013,   0.013,   0.022,   0.022,   0.025}, // 70
						     {0.019,   0.019,   0.019,    0.01,    0.01,   0.013,   0.013,   0.022,   0.022,   0.025}}; // 80
  
  static const double minDeltaPhiResHigh[11][10]  = {{-0.133,  -0.133,  -0.139,  -0.121,  -0.115,  -0.120,  -0.091,  -0.073,  -0.052,  -0.040}, // 0
						     {-0.088,  -0.088,  -0.094,  -0.085,  -0.085,  -0.055,  -0.058,  -0.040,  -0.028,  -0.028}, // 5
						     {-0.049,  -0.049,  -0.052,  -0.046,  -0.049,  -0.040,  -0.034,  -0.022,  -0.019,  -0.016}, // 10
						     {-0.034,  -0.034,  -0.037,  -0.031,  -0.034,  -0.026,  -0.025,  -0.019,  -0.013,  -0.016}, // 15
						     {-0.034,  -0.034,  -0.037,  -0.025,  -0.025,  -0.022,  -0.022,  -0.016,  -0.010,  -0.013}, // 20
						     {-0.034,  -0.034,  -0.037,  -0.019,  -0.019,  -0.018,  -0.016,  -0.013,  -0.010,  -0.010}, // 30
						     {-0.034,  -0.034,  -0.037,  -0.013,  -0.013,  -0.016,  -0.013,  -0.010,  -0.007,  -0.007}, // 40
						     {-0.034,  -0.034,  -0.037,  -0.013,  -0.013,  -0.016,  -0.013,  -0.010,  -0.007,  -0.007}, // 50
						     {-0.034,  -0.034,  -0.037,  -0.013,  -0.013,  -0.016,  -0.013,  -0.010,  -0.007,  -0.007}, // 60
						     {-0.034,  -0.034,  -0.037,  -0.013,  -0.013,  -0.016,  -0.013,  -0.010,  -0.007,  -0.007}, // 70
						     {-0.034,  -0.034,  -0.037,  -0.013,  -0.013,  -0.016,  -0.013,  -0.010,  -0.007,  -0.007}}; // 80
  
  static const double maxDeltaPhiResLow[11][10] = {{0.103,   0.103,   0.103,   0.076,   0.079,   0.081,   0.076,   0.064,   0.049,   0.070}, // 0
						   {0.043,   0.043,   0.043,   0.033,   0.029,   0.040,   0.030,   0.031,   0.037,   0.070}, // 5
						   {0.022,   0.022,   0.022,   0.020,   0.019,   0.028,   0.025,   0.025,   0.034,   0.064}, // 10
						   {0.019,   0.019,   0.019,   0.019,   0.017,   0.022,   0.020,   0.022,   0.031,   0.055}, // 15
						   {0.019,   0.019,   0.019,   0.017,   0.017,   0.020,   0.020,   0.025,   0.028,   0.049}, // 20
						   {0.019,   0.019,   0.019,   0.015,   0.015,   0.017,   0.017,   0.022,   0.022,   0.043}, // 30
						   {0.019,   0.019,   0.019,   0.013,   0.013,   0.016,   0.016,   0.019,   0.019,   0.037}, // 40
						   {0.019,   0.019,   0.019,   0.013,   0.013,   0.016,   0.016,   0.019,   0.019,   0.037}, // 50
						   {0.019,   0.019,   0.019,   0.013,   0.013,   0.016,   0.016,   0.019,   0.019,   0.037}, // 60
						   {0.019,   0.019,   0.019,   0.013,   0.013,   0.016,   0.016,   0.019,   0.019,   0.037}, // 70
						   {0.019,   0.019,   0.019,   0.013,   0.013,   0.016,   0.016,   0.019,   0.019,   0.037}}; // 80

   static const double minDeltaPhiResLow[11][10] = {{-0.100,  -0.100,  -0.103,  -0.100,  -0.094,  -0.082,  -0.070,  -0.061,  -0.049,  -0.043}, // 0
						    {-0.076,  -0.076,  -0.076,  -0.073,  -0.059,  -0.050,  -0.040,  -0.028,  -0.022,  -0.022}, // 5
						    {-0.043,  -0.043,   -0.04,  -0.040,  -0.030,  -0.030,  -0.023,  -0.013,  -0.013,  -0.016}, // 10
						    {-0.025,  -0.025,  -0.028,  -0.025,  -0.020,  -0.022,  -0.017,  -0.010,  -0.010,  -0.013}, // 15
						    {-0.025,  -0.025,  -0.028,  -0.019,  -0.016,  -0.019,  -0.014,  -0.010,  -0.007,  -0.010}, // 20
						    {-0.025,  -0.025,  -0.028,  -0.013,  -0.013,  -0.015,  -0.010,  -0.007,  -0.007,  -0.007}, // 30
						    {-0.025,  -0.025,  -0.028,  -0.010,  -0.010,  -0.013,  -0.007,  -0.007,  -0.007,  -0.007}, // 40
						    {-0.025,  -0.025,  -0.028,  -0.010,  -0.010,  -0.013,  -0.007,  -0.007,  -0.007,  -0.007}, // 50
						    {-0.025,  -0.025,  -0.028,  -0.010,  -0.010,  -0.013,  -0.007,  -0.007,  -0.007,  -0.007}, // 60
						    {-0.025,  -0.025,  -0.028,  -0.010,  -0.010,  -0.013,  -0.007,  -0.007,  -0.007,  -0.007}, // 70
						    {-0.025,  -0.025,  -0.028,  -0.010,  -0.010,  -0.013,  -0.007,  -0.007,  -0.007,  -0.007}}; // 80
 
  if(isBrem){
    if(deltaPhiRes < minDeltaPhiResHigh[eTBin][etaBin] || deltaPhiRes > maxDeltaPhiResHigh[eTBin][etaBin]) return false;
  }
  if(!isBrem){
    if(deltaPhiRes < minDeltaPhiResLow[eTBin][etaBin] || deltaPhiRes > maxDeltaPhiResLow[eTBin][etaBin]) return false;
  }
  return true;
  
}

//----------------------------------------------------------------------------------------

bool passTR(double TRratio, double eta, unsigned int  nTRTTotal ){  
  if (fabs(eta)>1.37 && fabs(eta)<1.52 ){
    return (nTRTTotal >0 && TRratio > 0.08);
  }
  return true; 
}


bool passConversion(double deltaPhiRes, bool expectBlayer,  int nBlayerHits ,unsigned int eTBin, unsigned int etaBin){  

                //                           0.0     0.1      0.6      0.8   1.15      1.37   1.52    1.81    2.01    2.37     2.47 
  static const double maxDeltaPhiResBL[11][10] = {{  0.020 , 0.020,  0.020,  0.020,  0.020,  0.016,  0.020,  0.023,  0.023,  0.032}, // 0  
						  {  0.014,  0.014,  0.014,  0.014,  0.014,  0.015,  0.020,  0.022,  0.022,  0.032}, // 5  
						  {  0.008,  0.008,  0.008,  0.009,  0.009,  0.011,  0.015,  0.015,  0.018,  0.030}, // 10 
						  {  0.006,  0.006,  0.007,  0.008,  0.008,  0.019,  0.013,  0.015,  0.017,  0.025}, // 15 
						  {  0.006,  0.006,  0.006,  0.006,  0.006,  0.008,  0.012,  0.013,  0.015,  0.021}, // 20 
						  {  0.006 , 0.006,  0.006,  0.006,  0.006,  0.007,  0.012,  0.013,  0.014,  0.020}, // 30 
						  {  0.005,  0.005,  0.005,  0.005,  0.005,  0.006,  0.012,  0.013,  0.014,  0.020}, // 40 
						  {  0.005,  0.005,  0.005,  0.005,  0.005,  0.006,  0.012,  0.013,  0.014,  0.020}, // 50 
						  {  0.005,  0.005,  0.005,  0.005,  0.005,  0.006,  0.012,  0.013,  0.014,  0.020}, // 60 
						  {  0.005,  0.005,  0.005,  0.005,  0.005,  0.006,  0.012,  0.013,  0.014,  0.020}, // 70 
						  {  0.005,  0.005,  0.005,  0.005,  0.005,  0.006,  0.012,  0.013,  0.014,  0.020}}; // 80


  if (expectBlayer && nBlayerHits<1) {	
     if(deltaPhiRes > maxDeltaPhiResBL[eTBin][etaBin]){
       return false;
     }
   }
  

  return true;
}
//============================================================================================================================================

