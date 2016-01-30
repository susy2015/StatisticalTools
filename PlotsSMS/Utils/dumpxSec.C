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

using namespace std;

void dumpxSec(std::string topoStr="T1"){ 

   TFile *xSecProspinoFile =0;
   TH1D *xSecProspino =0; 
   if( topoStr == "T1" ){
      xSecProspinoFile = new TFile("reference_gluino_xSec.root"); 
      xSecProspino = (TH1D*) xSecProspinoFile->Get("gluino_xsection"); 
   }

   for(int ib=1; ib<=xSecProspino->GetXaxis()->GetNbins(); ib++){
      double xVal = xSecProspino->GetXaxis()->GetBinCenter(ib);
      double xSec = xSecProspino->GetBinContent(ib);
      double xSecErr = xSecProspino->GetBinError(ib);
//      std::cout<<"Mass "<<xVal<<" xsec "<<xSec<<" err "<<xSecErr<<std::endl;
      printf("Mass : %4.0f  xSec : %8.6e\n", xVal, xSec);
   }

}
