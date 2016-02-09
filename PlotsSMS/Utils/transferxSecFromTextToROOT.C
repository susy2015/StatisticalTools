#include "TRint.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TStyle.h"

#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2F.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TPaveText.h"

#include <string>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cassert>

#include "TObjArray.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TObjString.h"

using namespace std;

void transferxSecFromTextToROOT(std::string inputStr="xSec_T3G.txt"){

   ifstream fin(inputStr.c_str());
   char line[200];
   TFile *xSecProspinoFile =0;
   TH1D *xSecProspino =0; 
   int minMom = 1000000, maxMom = 0;
   int nBins = 0;

   std::vector<int> momVec; std::vector<double> xSecVec, xSecErrVec;
   while( fin.getline(line, 200) ){
      TString lineT(line);
      if( lineT.Contains("Interactions") ) continue;
      TObjArray *vlist = lineT.Tokenize(" ");
      int nEntries = vlist->GetEntries();
      int mMom;
      double xSec =0, xSecRelErr =0;
      for(int ie=0; ie<nEntries; ie++){
         TObjString* perObj = dynamic_cast<TObjString*>(vlist->At(ie));
         TString perStr = perObj->GetString();
         if( ie==0 ){
            mMom = perStr.Atoi();
            if( minMom > mMom ) minMom = mMom;
            if( maxMom < mMom ) maxMom = mMom;
         }
         if( ie==1 ) xSec = perStr.Atof();
         if( ie==2 ) xSecRelErr = perStr.Atof();
      }
      nBins ++;
      momVec.push_back(mMom); xSecVec.push_back(xSec); xSecErrVec.push_back(xSec*xSecRelErr/100.);
//      std::cout<<"mMom : "<<mMom<<"  xSec : "<<xSec<<"  xSecRelErr : "<<xSecRelErr<<std::endl;
   }
   double divBin = 1.0*(maxMom - minMom)/(nBins-1);
   double lowMom = minMom-divBin/2.0, highMom = maxMom+divBin/2.0;
   std::cout<<"nBins : "<<nBins<<"  minMom : "<<minMom<<"  maxMom : "<<maxMom<<"  divBin : "<<divBin<<"  lowMom : "<<lowMom<<"  highMom : "<<highMom<<std::endl;

   TString rootStrT(inputStr);
   rootStrT.ReplaceAll("txt", "root");

   std::cout<<"root file : "<<rootStrT<<std::endl;
   xSecProspinoFile = new TFile(rootStrT, "RECREATE");

   if( rootStrT.Contains("T1") ) xSecProspino = new TH1D("gluino_xsection", "gluino_xsection", nBins, lowMom, highMom);
   if( rootStrT.Contains("T2") ) xSecProspino = new TH1D("squark_xsection", "squark_xsection", nBins, lowMom, highMom);
   if( rootStrT.Contains("T3G") ) xSecProspino = new TH1D("stop_xsection", "stop_xsection", nBins, lowMom, highMom);

   for(int iv=0; iv<(int)momVec.size(); iv++){
      double mMom = (double)momVec[iv];
      int ib = xSecProspino->FindFixBin(mMom);
      xSecProspino->SetBinContent(ib, xSecVec[iv]);
      xSecProspino->SetBinError(ib, xSecErrVec[iv]);
   }

//   xSecProspino->Write();
   xSecProspinoFile->Write(); xSecProspinoFile->Close();

}
