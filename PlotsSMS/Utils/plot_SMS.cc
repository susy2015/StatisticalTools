#include "plot.h"
#include "SusyScan.h"
#include "GeneratorMasses.h"
#include "PlotTools.h"
#include "TheLimits.h"
#include "GlobalFunctions.h"
#include "StyleSettings.h"

#include "interpolate.h"

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

#include "TColor.h"
#include "TPaletteAxis.h"

#include <string>
#include <cmath>
#include <stdio.h>
#include <sys/stat.h>

static const int nRegions = 5;
static const std::string  disptStrs[] = {">=1b + #slash{E}_{T}>200 GeV", ">=1b + #slash{E}_{T}>350 GeV", ">=2b + #slash{E}_{T}>200 GeV", ">=2b + #slash{E}_{T}>350 GeV", "combined"};
static const std::string regionDirs[] = {"baseline", "highMET", "baseline2b", "highMET2b", "combined"};
static const int     inclFlags[] = {    1,          1,          1,            1,           1    };

static const int      colors[] = {kGreen, kMagenta+1, kBlue, kRed, kTeal+4};

TH1D *xSecProspino_T1, *xSecProspino_T2, *xSecProspino_stop, *xSecProspino;
TString topoStr = "T1";
TString keyStr;

// 25 GeV per X bin; 25 GeV per Y bin
int nXbins = 65, nYbins = 73;
//const double loX = 100, hiX = 1225;
//const double loY = 50, hiY = 1225;
double loX = 587.5, hiX = 2212.5;
double loY = -12.5, hiY = 1812.5;

int nXbinsLG = 52, nYbinsLG = 49;
double loXLG = 612.5, hiXLG = 1912.5;
double loYLG = -12.5, hiYLG = 1212.5;

//const double scaleLegendFont = 1.8;
const double scaleLegendFont = 1.8;
const double widthXbin = 25, widthYbin = 25, widthBin = 25;
/*
const int nXbins = 60, nYbins = 60;
const double loX = 0, hiX = 1500;
const double loY = 0, hiY = 1500;
const double scaleLegendFont = 1.3;
*/

TGraphErrors *diagonalGraph =0;
TGraph  *diagonalCoverBand =0;

int plot(int argc, char** argv)
{
   //interactive root session
   //TApplication theApp("App", 0, 0);
   if (gROOT->IsBatch()) {
      fprintf(stderr, "%s: cannot run in batch mode\n", argv[0]);
      return 1;
   }
   if( argc < 3 ){
      std::cout<<"\n################USAGE##########################"<<std::endl;
      std::cout<<  "#./plot_SMS xxx yyy keyStr                    #"<<std::endl;
      std::cout<<  "# xxx: T1 or T2                               #"<<std::endl;
      std::cout<<  "# yyy: T2tt BR                                #"<<std::endl;
      std::cout<<  "# keystr: comb [(comb,stop,sbottom)]          #"<<std::endl;
      std::cout<<  "################USAGE##########################"<<std::endl;
      std::cout<<std::endl;
      return 1;
   }

   topoStr = argv[1];
   std::cout<<"\nWorking on topology "<<topoStr<<std::endl<<std::endl;

   double T2ttBR = atof(argv[2]);
   std::cout<<"\nT2ttBR : "<<T2ttBR<<std::endl<<std::endl;

   if( argc == 4 ) keyStr = argv[3];
   else keyStr = "NA";
   std::cout<<"\nthe keyStr for tuning : "<<keyStr<<std::endl<<std::endl; 

// For T1, it's reference_gluino_xSec.root
// For T2, it's reference_4squark_xSec.root
// For T2tt, it's reference_stop_xSec.root
   TFile *xSecProspinoFile = new TFile("reference_xSec.root");
   if( topoStr == "T1" || topoStr == "T1qqqq" || topoStr == "T5ZZInc" || topoStr == "T1tttt" || topoStr == "T1ttbb" || topoStr == "T1bbbb" || topoStr == "T5ttcc" || topoStr == "T5ttttDM175" || topoStr == "T5tttt_degen" ) xSecProspino_T1 = (TH1D*) xSecProspinoFile->Get("gluino_xsection");
   if( topoStr == "T2" || topoStr == "T2qq" ) xSecProspino_T2 = (TH1D*) xSecProspinoFile->Get("squark_xsection");
   if( topoStr == "T2tt" || topoStr == "T2bb" || topoStr == "T2tb" || topoStr == "T6ttWW" || topoStr == "T2fbd"){
      xSecProspino_stop = (TH1D*) xSecProspinoFile->Get("stop_xsection");
   }
   if( topoStr == "T1" || topoStr == "T1qqqq" || topoStr == "T5ZZInc" || topoStr == "T1tttt" || topoStr == "T1ttbb" || topoStr == "T1bbbb" || topoStr == "T5ttcc" || topoStr == "T5ttttDM175" || topoStr == "T5tttt_degen" ) xSecProspino = (TH1D*)xSecProspino_T1->Clone();
   if( topoStr == "T2" || topoStr == "T2qq" ) xSecProspino = (TH1D*)xSecProspino_T2->Clone();
   if( topoStr == "T2tt" || topoStr == "T2bb" || topoStr == "T2tb" || topoStr == "T6ttWW" || topoStr == "T2fbd") xSecProspino = (TH1D*)xSecProspino_stop->Clone();
   if( topoStr == "TGQ" ) xSecProspino = (TH1D*)xSecProspino_T2->Clone();

   for(int ib=0; ib<xSecProspino->GetNbinsX()+1; ib++){
      xSecProspino->SetBinContent(ib, xSecProspino->GetBinContent(ib));
      xSecProspino->SetBinError(ib, xSecProspino->GetBinError(ib));
   }

// 25 GeV per X bin; 25 GeV per Y bin
   if( topoStr == "T2tt" || topoStr == "T2tb" ){
      nXbins = 80; nYbins = 40;
      loX =  0; hiX = 2000;
      loY = 0; hiY = 1000;
 
      nXbinsLG = 47; nYbinsLG = 34;
      loXLG =  87.5; hiXLG = 1262.5;
      loYLG = -12.5; hiYLG = 837.5;
   }
   if( topoStr == "T2fbd"){
      nXbins = 80; nYbins = 40;
      loX =  0; hiX = 2000;
      loY = 0; hiY = 1000;
 
      nXbinsLG = 47; nYbinsLG = 34;
      loXLG =  87.5; hiXLG = 1262.5;
      loYLG = -12.5; hiYLG = 837.5;
   }
   if( topoStr == "T1tttt"){
      nXbins = 120; nYbins = 80;
      loX =  0; hiX = 3000;
      loY = 0; hiY = 2000;
 
      nXbinsLG = 47; nYbinsLG = 34;
      loXLG =  87.5; hiXLG = 1262.5;
      loYLG = -12.5; hiYLG = 837.5;
   }
   if( topoStr == "T6ttWW" ){
      nXbins = 27; nYbins = 36;
      loX = 287.5; hiX = 962.5;
      loY =  62.5; hiY = 962.5;

      nXbinsLG = 27; nYbinsLG = 36;
      loXLG = 287.5; hiXLG = 962.5;
      loYLG =  62.5; hiYLG = 962.5;
   }

   util::StyleSettings::paperNoTitle();
   gStyle->SetTitleSize(0.055,"XYZ");
   gStyle->SetTitleXOffset(1.02); gStyle->SetTitleYOffset(1.55); gStyle->SetTitleOffset(1.28, "Z");
   gStyle->SetPadLeftMargin(0.170); gStyle->SetPadRightMargin(0.125);
   gStyle->SetPadTopMargin(0.05); gStyle->SetPadBottomMargin(0.125); 

   gStyle->SetLabelSize(0.050, "XYZ");

   gStyle->SetNdivisions(505, "X");
 
   gStyle->SetTitleBorderSize(0);
   gStyle->SetTitleAlign(13);
   gStyle->SetTitleX(0.15);
   gStyle->SetTitleY(0.998);
   gStyle->SetTitleH(0.048);

   int nColors = 999;
   int MyPalette[nColors];
   Double_t r[]    = {0.0, 0.0, 0.0, 1.0, 1.0};
   Double_t g[]    = {0.0, 1.0, 1.0, 1.0, 0.0};
   Double_t b[]    = {1.0, 1.0, 0.0, 0.0, 0.0};
   Double_t stop[] = {0.000, 0.200, 0.400, 0.600, 0.800, 1.000};
   int FI = TColor::CreateGradientColorTable(5, stop, r, g, b, nColors);
   std::cout<<"\nColor FI : "<<FI<<std::endl;
   for(int i=0; i<nColors; i++) MyPalette[i] = FI+i;
   gStyle->SetNumberContours(nColors);

//   int initHighestColorIndex = gStyle->GetNumberOfColors();
//   std::cout<<"\nInitial color highestIndex : "<<initHighestColorIndex<<std::endl;

   gStyle->SetPaintTextFormat("3.1f");

   TCanvas * c1 = new TCanvas("c1","c1",600,600);
   c1->cd();
   
   //Get limits from signal scan ---------------------------------------------------
   char filename[200], dirname[200], tmpStr[200];

// Make a diagonal graph
   double hiDiagRange = -1, loDiagRange = 0;
   if( hiDiagRange < hiX ) hiDiagRange = hiX;
   if( hiDiagRange < hiY ) hiDiagRange = hiY;
   int nDiagBins = (hiDiagRange - loDiagRange)/widthBin +1;
   double xPtsDiag[nDiagBins], xErrsDiag[nDiagBins], yPtsDiag[nDiagBins], yErrsDiag[nDiagBins]; 
   for(int ib=0; ib<nDiagBins; ib++){
      xPtsDiag[ib] = loDiagRange + ib*widthBin; xErrsDiag[ib] = 0.0;
//      yPtsDiag[ib] = loDiagRange + (ib+1)*widthBin; yErrsDiag[ib] = 0.0;
      if( topoStr == "T2tt" ){
         yPtsDiag[ib] = loDiagRange + ib*widthBin - 175; yErrsDiag[ib] = 0.0;
      }else{
         yPtsDiag[ib] = loDiagRange + ib*widthBin; yErrsDiag[ib] = 0.0;
      }
   }
   diagonalGraph = new TGraphErrors(nDiagBins, xPtsDiag, yPtsDiag, xErrsDiag, yErrsDiag);
   diagonalGraph->SetLineWidth(2); diagonalGraph->SetLineStyle(7); diagonalGraph->SetLineColor(15);
   diagonalGraph->SetName("diagonalGraph");

   std::vector<TheLimits *> genpoints;
   std::vector<TString> region, regionBase, dispt;
   int nTotPlots = 0;

   for(int ir=0; ir<nRegions; ir++){
      if( inclFlags[ir] ){
         TheLimits *tmpLimits = new TheLimits();
         genpoints.push_back( tmpLimits );
         sprintf(tmpStr, "%s_", regionDirs[ir].c_str());
         nTotPlots++; region.push_back(tmpStr);
         sprintf(tmpStr, "%s", regionDirs[ir].c_str());
         regionBase.push_back(tmpStr);
         sprintf(tmpStr, "%s", disptStrs[ir].c_str());
         dispt.push_back(tmpStr);
      }
   }
//   nTotPlots++; region.push_back("combined_"); regionBase.push_back("combined"); dispt.push_back("combined");

   int idxCnt =-1;
   for(int ir=0; ir<nRegions; ir++){
      if( inclFlags[ir] ){
         sprintf(dirname, "%s_plots_%s/", regionDirs[ir].c_str(), topoStr.Data());
         sprintf(filename, "%s_plots_%s/filelist.txt", regionDirs[ir].c_str(), topoStr.Data());
         idxCnt++; genpoints[idxCnt]->Fill(filename, dirname);
//         genpoints[idxCnt]->FillEmptyPointsByInterpolationInM0M12();
//Get SUSY masses and x-sections from generator scan ----------------------------
//and match them to the signal scan
//         genpoints[idxCnt]->FillGeneratorMasses("GenScan_tb10.dat");
//         genpoints[idxCnt]->match();
      }
   }

   TString outDir = "results_"+topoStr;

   struct stat stFile;
   char commandline[200];
   if( stat(outDir.Data(), &stFile) == -1 ){ sprintf(commandline, "mkdir %s", outDir.Data()); system(commandline); }

   TGraph * gCLsExpExcl_NLO[nTotPlots];
   TGraph * gCLsObsExcl_NLO[nTotPlots];
   TGraph * gCLsExpExclm1[nTotPlots];
   TGraph * gCLsExpExclp1[nTotPlots];

   TH2* ObsExcl_NLO[nTotPlots];
   TH2* ExpExcl_NLO[nTotPlots];
   TH2* ExpExclm1[nTotPlots];
   TH2* ExpExclp1[nTotPlots];

   PlotTools<SusyScan> * plotTools[nTotPlots];

   std::vector<TH2F*> obslimitVec(nTotPlots), explimitVec(nTotPlots);
   std::vector<TGraph*> obsExclOneTimesProspinoVec(nTotPlots), expExclOneTimesProspinoVec(nTotPlots);
   std::vector<TGraph*> obsExclThreeTimesProspinoVec(nTotPlots), expExclThreeTimesProspinoVec(nTotPlots);
   std::vector<TGraph*> obsExclOneThirdProspinoVec(nTotPlots), expExclOneThirdProspinoVec(nTotPlots);
   std::vector<TGraph*> expExclPlusOneSigmaProspinoVec(nTotPlots), expExclMinusOneSigmaProspinoVec(nTotPlots);
   std::vector<TGraph*> obsExclPlusSysErrProspinoVec(nTotPlots), obsExclMinusSysErrProspinoVec(nTotPlots);

   std::vector<TGraph*> expExclOneTimesProspino_baseVec(nTotPlots), expExclOneTimesProspino_extraVec(nTotPlots);
   std::vector<TGraph*> expExclPlusOneSigmaProspino_baseVec(nTotPlots), expExclPlusOneSigmaProspino_extraVec(nTotPlots);
   std::vector<TGraph*> expExclMinusOneSigmaProspino_baseVec(nTotPlots), expExclMinusOneSigmaProspino_extraVec(nTotPlots);

   std::vector<TGraph*> obsExclOneTimesProspino_baseVec(nTotPlots), obsExclOneTimesProspino_extraVec(nTotPlots);
   std::vector<TGraph*> obsExclPlusSysErrProspino_baseVec(nTotPlots), obsExclPlusSysErrProspino_extraVec(nTotPlots);
   std::vector<TGraph*> obsExclMinusSysErrProspino_baseVec(nTotPlots), obsExclMinusSysErrProspino_extraVec(nTotPlots);

   TH2F *minCatsObsExclNLO, *minCatsExpExclNLO;
   TH2F *minObsExclNLO, *minExpExclNLO, *minObsExclNLOfollowminExp;
 
   TH2F *obsExclOneTimesxSecProspinoTestComb =0, *expExclOneTimesxSecProspinoTestComb =0;

   //Was not originally designed as a loop, so this will give you warnings about possible memory leaks.

//   TLegend* legexp = new TLegend(0.15,0.63,0.48,0.71,NULL,"brNDC"); 
   TLegend* legexp = new TLegend(0.18,0.65,0.48,0.80,NULL,"brNDC"); 
   legexp->SetFillColor(0);legexp->SetShadowColor(0);legexp->SetFillStyle(4000);legexp->SetTextFont(42);legexp->SetTextSize(0.023*scaleLegendFont);legexp->SetBorderSize(0);

   TLatex cmsPreTex; cmsPreTex.SetTextSize(0.025*scaleLegendFont); cmsPreTex.SetTextColor(1); cmsPreTex.SetNDC(); cmsPreTex.SetTextFont(42);
   TLatex cmsLumiTex; cmsLumiTex.SetTextSize(0.025*scaleLegendFont); cmsLumiTex.SetTextColor(1); cmsLumiTex.SetNDC(); cmsLumiTex.SetTextFont(42);
   TLatex disptTex; disptTex.SetTextSize(0.025*scaleLegendFont); disptTex.SetTextColor(1); disptTex.SetNDC(); disptTex.SetTextFont(42);

   TLatex indexTex; indexTex.SetTextSize(0.025*scaleLegendFont); indexTex.SetTextColor(1); indexTex.SetNDC(); indexTex.SetTextFont(42);

   for(int i = 0; i < nTotPlots; i++){
   //the plotting ----------------------------------------------------------------------
   //plotting helper functions
      plotTools[i] = new PlotTools<SusyScan>(genpoints[i]->GetScan());
      PlotTools<GeneratorMasses> * plotMasses = new PlotTools<GeneratorMasses>(genpoints[i]->GetGeneratorMasses());

      if(i < nTotPlots-1){

         bool dointerpolation = true;
         int dorebin = 1;
         if( keyStr == "raw" ){ dointerpolation = false; dorebin =0; }

   //iso mass lines
   //usage: Line( <x-var>, <y-var>, <mass-function>, <mass-value>, <matching-width> )
         TGraph * gl500 = plotMasses->Line(Mzero, Mhalf, MGluino, 500, 0.8 );
         TGraph * sq500 = plotMasses->Line(Mzero, Mhalf, MSquarkL, 500, 1);
         TGraph * chi100 = plotMasses->Line(Mzero, Mhalf, MChi1, 50, 20);
 
   //the histograms ---------------------------------------------------------------------
   // cross-section in M0 - M1/2
         c1->SetLogz(1);
         c1->SetRightMargin ( 0.19 );
         c1->SetTopMargin ( 0.085 );

   // Observed Limit in M0 - M1/2
         TH2F*hobslimit = 0;
         if( topoStr == "T1" || topoStr == "T1qqqq" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}#; m(#tilde{q})>>m(#tilde{g})");
         }else if( topoStr == "T2" || topoStr == "T2qq" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{q}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{q}#tilde{q}, #tilde{q}#rightarrow q#tilde{#chi}^{0}#; m(#tilde{g})>>m(#tilde{q})");
         }else if( topoStr == "T2tt" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{t}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
         }else if( topoStr == "T2fbd" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{t}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
         }else if( topoStr == "T2bb" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{b}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{b}#tilde{b}, #tilde{b}#rightarrow b#tilde{#chi}^{0}#; m(#tilde{b})>m(#tilde{#chi}^{0})");
         } else if( topoStr == "T2tb" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{t}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; #tilde{t}#rightarrow b#tilde{#chi}^{#pm}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
         } else if( topoStr == "TGQ" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{g}) [GeV]; m(#tilde{q}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{q}#tilde{g}, #tilde{q}#rightarrow q#tilde{#chi}^{0}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}");
         }else if( topoStr == "T5ZZInc" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}_{2}, #tilde{#chi}^{0}_{2}#rightarrow Z#tilde{#chi}^{0}");
         }else if( topoStr == "T1tttt" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }else if( topoStr == "T1ttbb" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }else if( topoStr == "T1bbbb" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow bb#tilde{#chi}^{0}#; m(#tilde{b})>>m(#tilde{g})");
         }else if( topoStr == "T5ttcc" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }else if( topoStr == "T5ttttDM175" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }else if( topoStr == "T5tttt_degen" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }else if( topoStr == "T6ttWW" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{b}) [GeV]; m(#tilde{#chi}^{#pm}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }
         plotTools[i]->Area(hobslimit, Mzero, Mhalf, NLOObsxSecCL, false);

         if( dointerpolation ){
            hobslimit = interpolate(hobslimit);
         }
         for(int ir=0; ir<dorebin; ir++){
            hobslimit = rebin(hobslimit);
         }

         hobslimit->SetMinimum(0.0001);
         hobslimit->SetMaximum(20);
         hobslimit->Draw("colz");
         TString obsHistName = region[i]+"obsLimit";
         obslimitVec[i] = (TH2F*) hobslimit->Clone(obsHistName.Data());

         TGraph *obsExclOneTimesxSecProspino = set_ProspinoExcl(xSecProspino, (TH2D*) hobslimit, 1.0);
         TGraph *obsExclThreeTimesxSecProspino = set_ProspinoExcl(xSecProspino, (TH2D*) hobslimit, 3.0);

         obsHistName = region[i]+"obsExclOneTimesProspino"; obsExclOneTimesProspinoVec[i] = (TGraph*) obsExclOneTimesxSecProspino->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclThreeTimesProspino"; obsExclThreeTimesProspinoVec[i] = (TGraph*) obsExclThreeTimesxSecProspino->Clone(obsHistName.Data());

         if( topoStr != "TGQ" ){
            obsExclOneTimesxSecProspino->SetLineWidth(2); obsExclOneTimesxSecProspino->Draw("same");
//            obsExclThreeTimesxSecProspino->SetLineWidth(2); obsExclThreeTimesxSecProspino->SetLineStyle(2); obsExclThreeTimesxSecProspino->Draw("same");
         }
         cmsPreTex.DrawLatex(0.17, 0.84, "CMS Preliminary");
         cmsLumiTex.DrawLatex(0.17, 0.79, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
         disptTex.DrawLatex(0.17, 0.74, dispt[i]);
/*
         cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
         cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
         disptTex.DrawLatex(350, 960, dispt[i]);
*/
         legexp->Clear();
         legexp->AddEntry(obsExclOneTimesxSecProspino, "#sigma^{prod}=#sigma^{NLO-QCD}", "l");
//         legexp->AddEntry(obsExclThreeTimesxSecProspino, "#sigma^{prod}=3#times#sigma^{NLO-QCD}", "l");
         if( topoStr != "TGQ" ) legexp->Draw();
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ObsLimit_mMother_mLSP.pdf");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ObsLimit_mMother_mLSP.png");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ObsLimit_mMother_mLSP.C");
   
   // Expected Limit in M0 - M1/2
//         TH2F*hexplimit = new TH2F("explimit",";m_{mother} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]", 48, 300, 1500, 56, 100, 1500);
         TH2F*hexplimit = 0;
         if( topoStr == "T1" || topoStr == "T1qqqq" ){
            hexplimit = new TH2F("explimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}#; m(#tilde{q})>>m(#tilde{g})");
         }else if( topoStr == "T2" || topoStr == "T2qq" ){
            hexplimit = new TH2F("explimit",";m(#tilde{q}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{q}#tilde{q}, #tilde{q}#rightarrow q#tilde{#chi}^{0}#; m(#tilde{g})>>m(#tilde{q})");
         }else if( topoStr == "T2tt" ){
            hexplimit = new TH2F("explimit",";m(#tilde{t}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
         }else if( topoStr == "T2fbd" ){
            hexplimit = new TH2F("explimit",";m(#tilde{t}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
         }else if( topoStr == "T2bb" ){
            hexplimit = new TH2F("explimit",";m(#tilde{b}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{b}#tilde{b}, #tilde{b}#rightarrow b#tilde{#chi}^{0}#; m(#tilde{b})>m(#tilde{#chi}^{0})");
         }else if( topoStr == "T2tb" ){
            hexplimit = new TH2F("explimit",";m(#tilde{t}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; #tilde{t}#rightarrow b#tilde{#chi}^{#pm}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
         } else if( topoStr == "TGQ" ){
            hexplimit = new TH2F("explimit",";m(#tilde{g}) [GeV]; m(#tilde{q}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{q}#tilde{g}, #tilde{q}#rightarrow q#tilde{#chi}^{0}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}");
         }else if( topoStr == "T5ZZInc"){
            hexplimit = new TH2F("explimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}_{2}, #tilde{#chi}^{0}_{2}#rightarrow Z#tilde{#chi}^{0}");
         }else if( topoStr == "T1tttt" ){
            hexplimit = new TH2F("explimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }else if( topoStr == "T1ttbb" ){
            hexplimit = new TH2F("explimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }else if( topoStr == "T1bbbb" ){
            hexplimit = new TH2F("explimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow bb#tilde{#chi}^{0}#; m(#tilde{b})>>m(#tilde{g})");
         }else if( topoStr == "T5ttcc" ){
            hexplimit = new TH2F("explimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }else if( topoStr == "T5ttttDM175" ){
            hexplimit = new TH2F("explimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }else if( topoStr == "T5tttt_degen" ){
            hexplimit = new TH2F("explimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }else if( topoStr == "T6ttWW" ){
            hexplimit = new TH2F("explimit",";m(#tilde{b}) [GeV]; m(#tilde{#chi}^{#pm}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }
         plotTools[i]->Area(hexplimit, Mzero, Mhalf, NLOExpxSecCL, false);

         if( dointerpolation ){
            hexplimit = interpolate(hexplimit);
         }
         for(int ir=0; ir<dorebin; ir++){
            hexplimit = rebin(hexplimit);
         }

	 hexplimit->SetTitle("");
	 hexplimit->GetXaxis()->SetRangeUser(150,1500);
	 hexplimit->GetYaxis()->SetRangeUser(20,800);
	 hexplimit->GetXaxis()->SetTitle("mStop [GeV]");
	 hexplimit->GetXaxis()->SetTitleSize(0.045);;
	 hexplimit->GetYaxis()->SetTitle("mLSP [GeV]");
	 hexplimit->GetYaxis()->SetTitleSize(0.045);;
	 hexplimit->GetXaxis()->SetLabelSize(0.04);
	 hexplimit->GetXaxis()->SetNdivisions(510);
	 if( topoStr == "T1tttt" )
         {
		hexplimit->GetXaxis()->SetRangeUser(600,2500);
		hexplimit->GetYaxis()->SetRangeUser(20,1800);
		hexplimit->GetXaxis()->SetTitle("mGluino [GeV]");
	 }
	 if( topoStr == "T2fbd" )
         {
		hexplimit->GetXaxis()->SetRangeUser(200,800);
		hexplimit->GetYaxis()->SetRangeUser(100,600);
	 }

         hexplimit->Draw("colz");
         TString expHistName = region[i]+"expLimit";
         explimitVec[i] = (TH2F*) hexplimit->Clone(expHistName.Data());

         TGraph *expExclOneTimesxSecProspino = set_ProspinoExcl(xSecProspino, (TH2D*) hexplimit, 1.0);
         TGraph *expExclThreeTimesxSecProspino = set_ProspinoExcl(xSecProspino, (TH2D*) hexplimit, 3.0);
         expHistName = region[i]+"expExclOneTimesProspino"; expExclOneTimesProspinoVec[i] = (TGraph*) expExclOneTimesxSecProspino->Clone(expHistName.Data());
         expHistName = region[i]+"expExclThreeTimesProspino"; expExclThreeTimesProspinoVec[i] = (TGraph*) expExclThreeTimesxSecProspino->Clone(expHistName.Data());

         if( topoStr != "TGQ" ){
            expExclOneTimesxSecProspino->SetLineWidth(2); expExclOneTimesxSecProspino->Draw("same");
         }

	 TLatex latex;
	 latex.SetTextSize(0.04);
	 latex.SetNDC();
         latex.DrawLatex(0.20, 0.85, "CMS Preliminary");
         latex.DrawLatex(0.20, 0.79, "L_{int} = 36 fb^{-1}, #sqrt{s} = 13 TeV");

         //cmsPreTex.DrawLatex(0.20, 0.85, "CMS Preliminary");
         //cmsLumiTex.DrawLatex(0.20, 0.79, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
//         disptTex.DrawLatex(0.17, 0.74, dispt[i]);
/*
         cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
         cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
         disptTex.DrawLatex(350, 960, dispt[i]);
*/
         //legexp->Clear();
         //legexp->AddEntry(expExclOneTimesxSecProspino, "#sigma^{prod}=#sigma^{NLO-QCD}", "l");
//         legexp->AddEntry(expExclThreeTimesxSecProspino, "#sigma^{prod}=3#times#sigma^{NLO-QCD}", "l");
         //if( topoStr != "TGQ" ) legexp->Draw();
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ExpLimit_mMother_mLSP.pdf");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ExpLimit_mMother_mLSP.png");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ExpLimit_mMother_mLSP.C");

         if( i==0 ){
//            minCatsObsExclNLO = (TH2F*) obslimitVec[i]->Clone("minCatsObsExclNLO"); minCatsExpExclNLO = (TH2F*) explimitVec[i]->Clone("minCatsExpExclNLO"); minCatsObsExclNLO->Reset(); minCatsExpExclNLO->Reset();
            minCatsObsExclNLO = new TH2F("minCatsObsExclNLO",";m(#tilde{t}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbinsLG, loXLG, hiXLG, nYbinsLG, loYLG, hiYLG); minCatsObsExclNLO->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
            minCatsExpExclNLO = new TH2F("minCatsExpExclNLO",";m(#tilde{t}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbinsLG, loXLG, hiXLG, nYbinsLG, loYLG, hiYLG); minCatsExpExclNLO->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0})");

            minObsExclNLO = (TH2F*) obslimitVec[i]->Clone("minObsExclNLO"); minExpExclNLO = (TH2F*) explimitVec[i]->Clone("minExpExclNLO");
            minObsExclNLOfollowminExp = (TH2F*) obslimitVec[i]->Clone("minObsExclNLOfollowminExp");

            TAxis *obsXaxis = (TAxis*) minObsExclNLO->GetXaxis(); TAxis *obsYaxis = (TAxis*) minObsExclNLO->GetYaxis();
            int obsXnbins = obsXaxis->GetNbins(), obsYnbins = obsYaxis->GetNbins();
            for(int ix=0; ix<obsXnbins; ix++){
               for(int iy=0; iy<obsYnbins; iy++){
                  if( obslimitVec[i]->GetBinContent(ix+1, iy+1) == 0 ) continue;
                  int xBinCats = minCatsObsExclNLO->GetXaxis()->FindBin( obslimitVec[i]->GetXaxis()->GetBinCenter(ix+1) );
                  int yBinCats = minCatsObsExclNLO->GetYaxis()->FindBin( obslimitVec[i]->GetYaxis()->GetBinCenter(iy+1) );
//                  minCatsObsExclNLO->SetBinContent(ix+1, iy+1, 1.0*i+1.0);
                  minCatsObsExclNLO->SetBinContent(xBinCats, yBinCats, 1.0*i+1.0);
               }
            }

            TAxis *expXaxis = (TAxis*) minExpExclNLO->GetXaxis(); TAxis *expYaxis = (TAxis*) minExpExclNLO->GetYaxis();
            int expXnbins = expXaxis->GetNbins(), expYnbins = expYaxis->GetNbins();
            for(int ix=0; ix<expXnbins; ix++){
               for(int iy=0; iy<expYnbins; iy++){
                  if( explimitVec[i]->GetBinContent(ix+1, iy+1) == 0 ) continue;
                  int xBinCats = minCatsExpExclNLO->GetXaxis()->FindBin( explimitVec[i]->GetXaxis()->GetBinCenter(ix+1) );
                  int yBinCats = minCatsExpExclNLO->GetYaxis()->FindBin( explimitVec[i]->GetYaxis()->GetBinCenter(iy+1) );
//                  minCatsExpExclNLO->SetBinContent(ix+1, iy+1, 1.0*i+1.0);
                  minCatsExpExclNLO->SetBinContent(xBinCats, yBinCats, 1.0*i+1.0);
               }
            }
         }else{
            TAxis *obsXaxis = (TAxis*) minObsExclNLO->GetXaxis(); TAxis *obsYaxis = (TAxis*) minObsExclNLO->GetYaxis();
            int obsXnbins = obsXaxis->GetNbins(), obsYnbins = obsYaxis->GetNbins();
            for(int ix=0; ix<obsXnbins; ix++){
               for(int iy=0; iy<obsYnbins; iy++){
                  double curMinVal = minObsExclNLO->GetBinContent(ix+1, iy+1);
                  double curVal = obslimitVec[i]->GetBinContent(ix+1, iy+1);
                  int xBinCats = minCatsObsExclNLO->GetXaxis()->FindBin( obslimitVec[i]->GetXaxis()->GetBinCenter(ix+1) );
                  int yBinCats = minCatsObsExclNLO->GetYaxis()->FindBin( obslimitVec[i]->GetYaxis()->GetBinCenter(iy+1) );
                  if( curVal < curMinVal ){ minCatsObsExclNLO->SetBinContent(xBinCats, yBinCats, 1.0*i+1.0); minObsExclNLO->SetBinContent(ix+1, iy+1, curVal); }
               }
            }

            TAxis *expXaxis = (TAxis*) minExpExclNLO->GetXaxis(); TAxis *expYaxis = (TAxis*) minExpExclNLO->GetYaxis();
            int expXnbins = expXaxis->GetNbins(), expYnbins = expYaxis->GetNbins();
            for(int ix=0; ix<expXnbins; ix++){
               for(int iy=0; iy<expYnbins; iy++){
                  double curMinVal = minExpExclNLO->GetBinContent(ix+1, iy+1);
                  double curVal = explimitVec[i]->GetBinContent(ix+1, iy+1);
                  int xBinCats = minCatsExpExclNLO->GetXaxis()->FindBin( explimitVec[i]->GetXaxis()->GetBinCenter(ix+1) );
                  int yBinCats = minCatsExpExclNLO->GetYaxis()->FindBin( explimitVec[i]->GetYaxis()->GetBinCenter(iy+1) );
                  if( curVal < curMinVal ){
                     minCatsExpExclNLO->SetBinContent(xBinCats, yBinCats, 1.0*i+1.0); minExpExclNLO->SetBinContent(ix+1, iy+1, curVal); 
                     minObsExclNLOfollowminExp->SetBinContent(ix+1, iy+1, obslimitVec[i]->GetBinContent(ix+1, iy+1));
                  }
               }
            }
         }
      }
      if(i == nTotPlots-1){
   // Observed Limit in M0 - M1/2
         bool doMyFill = true;
// xBinLargery is a relative one. It depends on starting points of x and y axis...
//         int xBinLargery = 2;
         int xBinLargery = 0;
         if( topoStr == "T5ZZInc" ) xBinLargery = 2;
         if( topoStr == "T1bbbb" ) xBinLargery = 0;
         if( topoStr == "T1tttt" ) xBinLargery = 0;

         const int moreISRremoval = 3;
//         xBinLargery += moreISRremoval;
//         xBinLargery =0;

         TString obsHistName = region[i]+"obsLimit";
         TH2F * minObsExclNLOcloned = (TH2F*) minObsExclNLO->Clone(); minObsExclNLOcloned->Reset();
         for(int ix=0; ix<=minObsExclNLO->GetXaxis()->GetNbins(); ix++){
            if( ix < xBinLargery ) continue;
            for(int iy=0; iy<=minObsExclNLO->GetYaxis()->GetNbins(); iy++){
//               if( ix-iy < xBinLargery ) continue;
               minObsExclNLOcloned->SetBinContent(ix, iy, minObsExclNLO->GetBinContent(ix, iy));
            }
         }

         bool doSmooth = true;

         int radiusToSmooth = 10;
//         int radiusToSmooth = 6;
         int radiusToSmoothTH2 = 4;
         bool dointerpolation = true;
// smoothMethod : 0 --> the "k5a",  1 --> the "k3a",  2 --> C.A.'s method
//                -1 --> no smoothing
         int smoothMethod = -1;
         int dorebin = 2;

         if( keyStr == "raw" ){ doSmooth = false; dointerpolation = false; dorebin = 0; }

         if( topoStr == "T5ZZInc" ) radiusToSmooth = 16;
         if( topoStr == "T1bbbb" ) radiusToSmooth = 16;
         if( topoStr == "T1tttt" ) radiusToSmooth = 8;

         double obsOneTimesCutOffX = -1, obsThreeTimesCutOffX = -1, obsOneThirdCutOffX = -1;
         double obsPlusSysErrCutOffX = -1, obsMinusSysErrCutOffX = -1;

         TH2F * obsExclOneTimesxSecProspinoSmooth = new TH2F("obsExclOneTimesxSecProspinoSmooth", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[i]->Area(obsExclOneTimesxSecProspinoSmooth, Mzero, Mhalf, NLOObsxSecCL, false);
         if( dointerpolation ){
            obsExclOneTimesxSecProspinoSmooth = interpolate(obsExclOneTimesxSecProspinoSmooth); minObsExclNLOcloned = interpolate(minObsExclNLOcloned);
         }
         for(int ir=0; ir<dorebin; ir++){ 
            obsExclOneTimesxSecProspinoSmooth = rebin(obsExclOneTimesxSecProspinoSmooth); //minObsExclNLOcloned = rebin(minObsExclNLOcloned);
         }
         if( smoothMethod == 0 ){
            obsExclOneTimesxSecProspinoSmooth->Smooth(1, "k5a"); minObsExclNLOcloned->Smooth(1, "k5a");
         }else if( smoothMethod == 1 ){
            obsExclOneTimesxSecProspinoSmooth->Smooth(1, "k3a"); minObsExclNLOcloned->Smooth(1, "k3a");
         }else if( smoothMethod == 2 ){
            Smooth(obsExclOneTimesxSecProspinoSmooth, radiusToSmoothTH2); Smooth(minObsExclNLOcloned, radiusToSmoothTH2);
         }

         obslimitVec[i] = (TH2F*) minObsExclNLOcloned->Clone(obsHistName.Data());
         obslimitVec[i]->SetMinimum(0.0001);
         obslimitVec[i]->SetMaximum(20);

         TH2F * obsExclOneTimesxSecProspinoTest = 0;
         if( !radiusToSmoothTH2 && !dorebin ){
            obsExclOneTimesxSecProspinoTest = new TH2F("obsExclOneTimesxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
            plotTools[i]->Area(obsExclOneTimesxSecProspinoTest, Mzero, Mhalf, NLOObsxSecCL, xSecProspino, 1.0, doMyFill, xBinLargery);
         }else{
            obsExclOneTimesxSecProspinoTest = (TH2F*) obsExclOneTimesxSecProspinoSmooth->Clone(); obsExclOneTimesxSecProspinoTest->Reset();
            plotTools[i]->Area(obsExclOneTimesxSecProspinoTest, Mzero, Mhalf, NLOObsxSecCL, xSecProspino, 1.0, doMyFill, xBinLargery, 0, obsExclOneTimesxSecProspinoSmooth);
         }

         int ixb = obsExclOneTimesxSecProspinoTest->GetXaxis()->FindBin(575);
         int iyb = obsExclOneTimesxSecProspinoTest->GetYaxis()->FindBin(0.0);
//         obsExclOneTimesxSecProspinoTest->SetBinContent(ixb, iyb, 1);
         ixb = obsExclOneTimesxSecProspinoTest->GetXaxis()->FindBin(575);
         iyb = obsExclOneTimesxSecProspinoTest->GetYaxis()->FindBin(0.0);
//         obsExclOneTimesxSecProspinoTest->SetBinContent(ixb, iyb, 1);

         TGraph * obsExclOneTimesxSecProspinoTMP_base = plotTools[i]->GetContour(obsExclOneTimesxSecProspinoTest, 3, 1, !doMyFill);
         TGraph * obsExclOneTimesxSecProspinoTMP_extra = plotTools[i]->GetContour(obsExclOneTimesxSecProspinoTest, 3, 2, !doMyFill);
         if( keyStr == "fill" ) obsExclOneTimesxSecProspinoTMP_base = plotTools[i]->GetContour(obsExclOneTimesxSecProspinoTest, 3, 1, doMyFill);
         double tmpxgr=-1, tmpygr=-1;
//         obsExclOneTimesxSecProspinoTMP_base->GetPoint(obsExclOneTimesxSecProspinoTMP_base->GetN()-1, tmpxgr, tmpygr);
//         obsExclOneTimesxSecProspinoTMP_base->SetPoint(obsExclOneTimesxSecProspinoTMP_base->GetN(), tmpxgr, 0);
         if( doSmooth ) Smooth(obsExclOneTimesxSecProspinoTMP_base, radiusToSmooth+2);
         if( doSmooth ) Smooth(obsExclOneTimesxSecProspinoTMP_extra, radiusToSmooth+2);

         std::vector<double> xPts_obsExclOneTimesxSecProspino_base, yPts_obsExclOneTimesxSecProspino_base;
         for(int ip=0; ip<obsExclOneTimesxSecProspinoTMP_base->GetN(); ip++){
            double xgr=-1, ygr=-1;
            obsExclOneTimesxSecProspinoTMP_base->GetPoint(ip, xgr, ygr);
            if( xgr >= obsOneTimesCutOffX ){
               if( xgr != obsOneTimesCutOffX && xPts_obsExclOneTimesxSecProspino_base.empty() ){
                  xPts_obsExclOneTimesxSecProspino_base.push_back(obsOneTimesCutOffX); yPts_obsExclOneTimesxSecProspino_base.push_back(ygr);
               }
               xPts_obsExclOneTimesxSecProspino_base.push_back(xgr); yPts_obsExclOneTimesxSecProspino_base.push_back(ygr);
            }
         }
         int nPts_obsExclOneTimesxSecProspino_base = (int)xPts_obsExclOneTimesxSecProspino_base.size();
         TGraph *obsExclOneTimesxSecProspino_base = new TGraph(nPts_obsExclOneTimesxSecProspino_base, &xPts_obsExclOneTimesxSecProspino_base[0], &yPts_obsExclOneTimesxSecProspino_base[0]);

         std::vector<double> xPts_obsExclOneTimesxSecProspino_extra, yPts_obsExclOneTimesxSecProspino_extra;
         for(int ip=0; ip<obsExclOneTimesxSecProspinoTMP_extra->GetN(); ip++){
            double xgr=-1, ygr=-1;
            obsExclOneTimesxSecProspinoTMP_extra->GetPoint(ip, xgr, ygr);
            if( xgr >= obsOneTimesCutOffX ){
               if( xgr != obsOneTimesCutOffX && xPts_obsExclOneTimesxSecProspino_extra.empty() ){
                  xPts_obsExclOneTimesxSecProspino_extra.push_back(obsOneTimesCutOffX); yPts_obsExclOneTimesxSecProspino_extra.push_back(ygr);
               }
               xPts_obsExclOneTimesxSecProspino_extra.push_back(xgr); yPts_obsExclOneTimesxSecProspino_extra.push_back(ygr);
            }
         }
         int nPts_obsExclOneTimesxSecProspino_extra = (int)xPts_obsExclOneTimesxSecProspino_extra.size();
         TGraph *obsExclOneTimesxSecProspino_extra = new TGraph(nPts_obsExclOneTimesxSecProspino_extra, &xPts_obsExclOneTimesxSecProspino_extra[0], &yPts_obsExclOneTimesxSecProspino_extra[0]);

         TGraph *obsExclOneTimesxSecProspino = MakeBand(obsExclOneTimesxSecProspino_base, obsExclOneTimesxSecProspino_extra);
/*
         if( topoStr == "T2" || topoStr == "T2qq" ){ obsExclOneTimesxSecProspino->RemovePoint(0); obsExclOneTimesxSecProspino->RemovePoint(0); }
         if( topoStr == "T5ZZInc" ){ obsExclOneTimesxSecProspino->RemovePoint(0); }
         if( topoStr == "T1bbbb" ){ obsExclOneTimesxSecProspino->RemovePoint(0); obsExclOneTimesxSecProspino->RemovePoint(0); }
*/

         obsExclOneTimesxSecProspinoTestComb = (TH2F*) obsExclOneTimesxSecProspinoTest->Clone("obsExclOneTimesxSecTestComb");

         obsExclOneTimesxSecProspinoTest->Draw("colz");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"_Obs_TestContour.pdf");

         TH2F * obsExclThreeTimesxSecProspinoTest = new TH2F("obsExclThreeTimesxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[i]->Area(obsExclThreeTimesxSecProspinoTest, Mzero, Mhalf, NLOObsxSecCL, xSecProspino, 3.0, doMyFill, xBinLargery);
         TGraph * obsExclThreeTimesxSecProspinoTMP = plotTools[i]->GetContour(obsExclThreeTimesxSecProspinoTest, 3, 1, !doMyFill);
         if( doSmooth ) Smooth(obsExclThreeTimesxSecProspinoTMP, radiusToSmooth+2);
         std::vector<double> xPts_obsExclThreeTimesxSecProspino, yPts_obsExclThreeTimesxSecProspino;
         for(int ip=0; ip<obsExclThreeTimesxSecProspinoTMP->GetN(); ip++){
            double xgr=-1, ygr=-1;
            obsExclThreeTimesxSecProspinoTMP->GetPoint(ip, xgr, ygr);
            if( xgr >= obsThreeTimesCutOffX ){
               if( xgr != obsThreeTimesCutOffX && xPts_obsExclThreeTimesxSecProspino.empty() ){
                  xPts_obsExclThreeTimesxSecProspino.push_back(obsThreeTimesCutOffX); yPts_obsExclThreeTimesxSecProspino.push_back(ygr);
               }
               xPts_obsExclThreeTimesxSecProspino.push_back(xgr); yPts_obsExclThreeTimesxSecProspino.push_back(ygr); 
            }
         }
         int nPts_obsExclThreeTimesxSecProspino = (int)xPts_obsExclThreeTimesxSecProspino.size();
         TGraph *obsExclThreeTimesxSecProspino = new TGraph(nPts_obsExclThreeTimesxSecProspino, &xPts_obsExclThreeTimesxSecProspino[0], &yPts_obsExclThreeTimesxSecProspino[0]);

         TH2F * obsExclOneThirdxSecProspinoTest = new TH2F("obsExclOneThirdxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[i]->Area(obsExclOneThirdxSecProspinoTest, Mzero, Mhalf, NLOObsxSecCL, xSecProspino, 1./3.0, doMyFill, xBinLargery);
         TGraph * obsExclOneThirdxSecProspinoTMP = plotTools[i]->GetContour(obsExclOneThirdxSecProspinoTest, 3, 1, !doMyFill);
         if( doSmooth ) Smooth(obsExclOneThirdxSecProspinoTMP, radiusToSmooth);
         std::vector<double> xPts_obsExclOneThirdxSecProspino, yPts_obsExclOneThirdxSecProspino;
         for(int ip=0; ip<obsExclOneThirdxSecProspinoTMP->GetN(); ip++){
            double xgr=-1, ygr=-1;
            obsExclOneThirdxSecProspinoTMP->GetPoint(ip, xgr, ygr);
            if( xgr >= obsOneThirdCutOffX ){
               if( xgr != obsOneThirdCutOffX && xPts_obsExclOneThirdxSecProspino.empty() ){
                  xPts_obsExclOneThirdxSecProspino.push_back(obsOneThirdCutOffX); yPts_obsExclOneThirdxSecProspino.push_back(ygr);
               }
               xPts_obsExclOneThirdxSecProspino.push_back(xgr); yPts_obsExclOneThirdxSecProspino.push_back(ygr); 
            }
         }
         int nPts_obsExclOneThirdxSecProspino = (int)xPts_obsExclOneThirdxSecProspino.size();
         TGraph *obsExclOneThirdxSecProspino = new TGraph(nPts_obsExclOneThirdxSecProspino, &xPts_obsExclOneThirdxSecProspino[0], &yPts_obsExclOneThirdxSecProspino[0]);

// Sys Error on obs limit
         TH2F * obsExclPlusSysErrxSecProspinoSmooth = new TH2F("obsExclPlusSysErrxSecProspinoSmooth", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[i]->Area(obsExclPlusSysErrxSecProspinoSmooth, Mzero, Mhalf, NLOObsxSecCL, false);
         if( dointerpolation ){
            obsExclPlusSysErrxSecProspinoSmooth = interpolate(obsExclPlusSysErrxSecProspinoSmooth);
         }
         for(int ir=0; ir<dorebin; ir++){ 
            obsExclPlusSysErrxSecProspinoSmooth = rebin(obsExclPlusSysErrxSecProspinoSmooth);
         }  
         if( smoothMethod == 0 ){
            obsExclPlusSysErrxSecProspinoSmooth->Smooth(1, "k5a");
         }else if( smoothMethod == 1 ){
            obsExclPlusSysErrxSecProspinoSmooth->Smooth(1, "k3a");
         }else if( smoothMethod == 2 ){
            Smooth(obsExclPlusSysErrxSecProspinoSmooth, radiusToSmoothTH2);
         }

         int tmpixb = obsExclPlusSysErrxSecProspinoSmooth->GetXaxis()->FindBin(500); 
         int tmpiyb = obsExclPlusSysErrxSecProspinoSmooth->GetYaxis()->FindBin(200); 
//         obsExclPlusSysErrxSecProspinoSmooth->SetBinContent(tmpixb, tmpiyb, 0.154);

         TH2F * obsExclPlusSysErrxSecProspinoTest = 0;
         if( !radiusToSmoothTH2 && !dorebin ){
            obsExclPlusSysErrxSecProspinoTest = new TH2F("obsExclPlusSysErrxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
            plotTools[i]->Area(obsExclPlusSysErrxSecProspinoTest, Mzero, Mhalf, NLOObsxSecCL, xSecProspino, 1.0, doMyFill, xBinLargery, 1);
         }else{
            obsExclPlusSysErrxSecProspinoTest = (TH2F*) obsExclPlusSysErrxSecProspinoSmooth->Clone(); obsExclPlusSysErrxSecProspinoTest->Reset();       
            plotTools[i]->Area(obsExclPlusSysErrxSecProspinoTest, Mzero, Mhalf, NLOObsxSecCL, xSecProspino, 1.0, doMyFill, xBinLargery, 1, obsExclPlusSysErrxSecProspinoSmooth);
         }

         TGraph * obsExclPlusSysErrxSecProspinoTMP_base = plotTools[i]->GetContour(obsExclPlusSysErrxSecProspinoTest, 3, 1, !doMyFill);
         TGraph * obsExclPlusSysErrxSecProspinoTMP_extra = plotTools[i]->GetContour(obsExclPlusSysErrxSecProspinoTest, 3, 2, !doMyFill);
         if( keyStr == "fill" ) obsExclPlusSysErrxSecProspinoTMP_base = plotTools[i]->GetContour(obsExclPlusSysErrxSecProspinoTest, 3, 1, doMyFill);
//         obsExclPlusSysErrxSecProspinoTMP_base->GetPoint(obsExclPlusSysErrxSecProspinoTMP_base->GetN()-1, tmpxgr, tmpygr);
//         obsExclPlusSysErrxSecProspinoTMP->SetPoint(obsExclPlusSysErrxSecProspinoTMP->GetN(), tmpxgr, 0);
         if( doSmooth ) Smooth(obsExclPlusSysErrxSecProspinoTMP_base, radiusToSmooth+3);
         if( doSmooth ) Smooth(obsExclPlusSysErrxSecProspinoTMP_extra, radiusToSmooth+3);

         std::vector<double> xPts_obsExclPlusSysErrxSecProspino_base, yPts_obsExclPlusSysErrxSecProspino_base;
         for(int ip=0; ip<obsExclPlusSysErrxSecProspinoTMP_base->GetN(); ip++){
            double xgr=-1, ygr=-1;
            obsExclPlusSysErrxSecProspinoTMP_base->GetPoint(ip, xgr, ygr);
            if( xgr >= obsPlusSysErrCutOffX ){
               if( xgr != obsPlusSysErrCutOffX && xPts_obsExclPlusSysErrxSecProspino_base.empty() ){
                  xPts_obsExclPlusSysErrxSecProspino_base.push_back(obsPlusSysErrCutOffX); yPts_obsExclPlusSysErrxSecProspino_base.push_back(ygr);
               }
               xPts_obsExclPlusSysErrxSecProspino_base.push_back(xgr); yPts_obsExclPlusSysErrxSecProspino_base.push_back(ygr);
            }
         }
         int nPts_obsExclPlusSysErrxSecProspino_base = (int)xPts_obsExclPlusSysErrxSecProspino_base.size();
         TGraph *obsExclPlusSysErrxSecProspino_base = new TGraph(nPts_obsExclPlusSysErrxSecProspino_base, &xPts_obsExclPlusSysErrxSecProspino_base[0], &yPts_obsExclPlusSysErrxSecProspino_base[0]);

         std::vector<double> xPts_obsExclPlusSysErrxSecProspino_extra, yPts_obsExclPlusSysErrxSecProspino_extra;
         for(int ip=0; ip<obsExclPlusSysErrxSecProspinoTMP_extra->GetN(); ip++){
            double xgr=-1, ygr=-1;
            obsExclPlusSysErrxSecProspinoTMP_extra->GetPoint(ip, xgr, ygr);
            if( xgr >= obsPlusSysErrCutOffX ){
               if( xgr != obsPlusSysErrCutOffX && xPts_obsExclPlusSysErrxSecProspino_base.empty() ){
                  xPts_obsExclPlusSysErrxSecProspino_extra.push_back(obsPlusSysErrCutOffX); yPts_obsExclPlusSysErrxSecProspino_extra.push_back(ygr);
               }
               xPts_obsExclPlusSysErrxSecProspino_extra.push_back(xgr); yPts_obsExclPlusSysErrxSecProspino_extra.push_back(ygr);
            }
         }
         int nPts_obsExclPlusSysErrxSecProspino_extra = (int)xPts_obsExclPlusSysErrxSecProspino_extra.size();
         TGraph *obsExclPlusSysErrxSecProspino_extra = new TGraph(nPts_obsExclPlusSysErrxSecProspino_extra, &xPts_obsExclPlusSysErrxSecProspino_extra[0], &yPts_obsExclPlusSysErrxSecProspino_extra[0]);

         TGraph *obsExclPlusSysErrxSecProspino = MakeBand(obsExclPlusSysErrxSecProspino_base, obsExclPlusSysErrxSecProspino_extra);
/*
         if( topoStr == "T2" || topoStr == "T2qq" ){ obsExclPlusSysErrxSecProspino->RemovePoint(0); }
         if( topoStr == "T5ZZInc" ){ obsExclPlusSysErrxSecProspino->RemovePoint(0); }
         if( topoStr == "T1bbbb" ){ obsExclPlusSysErrxSecProspino->RemovePoint(0); obsExclPlusSysErrxSecProspino->RemovePoint(0); }
*/


         TH2F * obsExclMinusSysErrxSecProspinoSmooth = new TH2F("obsExclMinusSysErrxSecProspinoSmooth", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[i]->Area(obsExclMinusSysErrxSecProspinoSmooth, Mzero, Mhalf, NLOObsxSecCL, false);
         if( dointerpolation ){
            obsExclMinusSysErrxSecProspinoSmooth = interpolate(obsExclMinusSysErrxSecProspinoSmooth);
         }
         for(int ir=0; ir<dorebin; ir++){ 
            obsExclMinusSysErrxSecProspinoSmooth = rebin(obsExclMinusSysErrxSecProspinoSmooth);
         }
         if( smoothMethod == 0 ){
            obsExclMinusSysErrxSecProspinoSmooth->Smooth(1, "k5a");
         }else if( smoothMethod == 1 ){
            obsExclMinusSysErrxSecProspinoSmooth->Smooth(1, "k3a");
         }else if( smoothMethod == 2 ){
            Smooth(obsExclMinusSysErrxSecProspinoSmooth, radiusToSmoothTH2);
         }

         tmpixb = obsExclMinusSysErrxSecProspinoSmooth->GetXaxis()->FindBin(475); 
         tmpiyb = obsExclMinusSysErrxSecProspinoSmooth->GetYaxis()->FindBin(175); 
//         obsExclMinusSysErrxSecProspinoSmooth->SetBinContent(tmpixb, tmpiyb, 0.130);

         TH2F * obsExclMinusSysErrxSecProspinoTest = 0;
         if( !radiusToSmoothTH2 && !dorebin ){
            obsExclMinusSysErrxSecProspinoTest = new TH2F("obsExclMinusSysErrxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
            plotTools[i]->Area(obsExclMinusSysErrxSecProspinoTest, Mzero, Mhalf, NLOObsxSecCL, xSecProspino, 1.0, doMyFill, xBinLargery, -1);
         }else{
            obsExclMinusSysErrxSecProspinoTest = (TH2F*) obsExclMinusSysErrxSecProspinoSmooth->Clone(); obsExclMinusSysErrxSecProspinoTest->Reset();
            plotTools[i]->Area(obsExclMinusSysErrxSecProspinoTest, Mzero, Mhalf, NLOObsxSecCL, xSecProspino, 1.0, doMyFill, xBinLargery, -1, obsExclMinusSysErrxSecProspinoSmooth);
         }

         TGraph * obsExclMinusSysErrxSecProspinoTMP_base = plotTools[i]->GetContour(obsExclMinusSysErrxSecProspinoTest, 3, 1, !doMyFill);
         TGraph * obsExclMinusSysErrxSecProspinoTMP_extra = plotTools[i]->GetContour(obsExclMinusSysErrxSecProspinoTest, 3, 2, !doMyFill);
         if( keyStr == "fill" ) obsExclMinusSysErrxSecProspinoTMP_base = plotTools[i]->GetContour(obsExclMinusSysErrxSecProspinoTest, 3, 1, doMyFill);
//         obsExclMinusSysErrxSecProspinoTMP_base->GetPoint(obsExclMinusSysErrxSecProspinoTMP_base->GetN()-1, tmpxgr, tmpygr);
//         obsExclMinusSysErrxSecProspinoTMP_base->SetPoint(obsExclMinusSysErrxSecProspinoTMP_base->GetN(), tmpxgr, 0);
         for(int ip=0; ip<obsExclMinusSysErrxSecProspinoTMP_base->GetN(); ip++){
            double xgr=-1, ygr=-1;
            obsExclMinusSysErrxSecProspinoTMP_base->GetPoint(ip, xgr, ygr);
//            if( xgr <= 400 && xgr >=290 ) obsExclMinusSysErrxSecProspinoTMP->SetPoint(ip, xgr, ygr-8);
         }
         if( doSmooth ) Smooth(obsExclMinusSysErrxSecProspinoTMP_base, radiusToSmooth+3);
         if( doSmooth ) Smooth(obsExclMinusSysErrxSecProspinoTMP_extra, radiusToSmooth+3);

         std::vector<double> xPts_obsExclMinusSysErrxSecProspino_base, yPts_obsExclMinusSysErrxSecProspino_base;
         for(int ip=0; ip<obsExclMinusSysErrxSecProspinoTMP_base->GetN(); ip++){
            double xgr=-1, ygr=-1;
            obsExclMinusSysErrxSecProspinoTMP_base->GetPoint(ip, xgr, ygr);
            if( xgr >= obsMinusSysErrCutOffX ){
               if( xgr != obsMinusSysErrCutOffX && xPts_obsExclMinusSysErrxSecProspino_base.empty() ){
                  xPts_obsExclMinusSysErrxSecProspino_base.push_back(obsMinusSysErrCutOffX); yPts_obsExclMinusSysErrxSecProspino_base.push_back(ygr);
               }
               xPts_obsExclMinusSysErrxSecProspino_base.push_back(xgr); yPts_obsExclMinusSysErrxSecProspino_base.push_back(ygr); 
            }
         }
         int nPts_obsExclMinusSysErrxSecProspino_base = (int)xPts_obsExclMinusSysErrxSecProspino_base.size();
         TGraph *obsExclMinusSysErrxSecProspino_base = new TGraph(nPts_obsExclMinusSysErrxSecProspino_base, &xPts_obsExclMinusSysErrxSecProspino_base[0], &yPts_obsExclMinusSysErrxSecProspino_base[0]);

         std::vector<double> xPts_obsExclMinusSysErrxSecProspino_extra, yPts_obsExclMinusSysErrxSecProspino_extra;
         for(int ip=0; ip<obsExclMinusSysErrxSecProspinoTMP_extra->GetN(); ip++){
            double xgr=-1, ygr=-1;
            obsExclMinusSysErrxSecProspinoTMP_extra->GetPoint(ip, xgr, ygr);
            if( xgr >= obsMinusSysErrCutOffX ){
               if( xgr != obsMinusSysErrCutOffX && xPts_obsExclMinusSysErrxSecProspino_extra.empty() ){
                  xPts_obsExclMinusSysErrxSecProspino_extra.push_back(obsMinusSysErrCutOffX); yPts_obsExclMinusSysErrxSecProspino_extra.push_back(ygr);
               }
               xPts_obsExclMinusSysErrxSecProspino_extra.push_back(xgr); yPts_obsExclMinusSysErrxSecProspino_extra.push_back(ygr);
            }
         }
         int nPts_obsExclMinusSysErrxSecProspino_extra = (int)xPts_obsExclMinusSysErrxSecProspino_extra.size();
         TGraph *obsExclMinusSysErrxSecProspino_extra = new TGraph(nPts_obsExclMinusSysErrxSecProspino_extra, &xPts_obsExclMinusSysErrxSecProspino_extra[0], &yPts_obsExclMinusSysErrxSecProspino_extra[0]);

         TGraph *obsExclMinusSysErrxSecProspino = MakeBand(obsExclMinusSysErrxSecProspino_base, obsExclMinusSysErrxSecProspino_extra);
/*
         if( topoStr == "T5ZZInc" ){ obsExclMinusSysErrxSecProspino->RemovePoint(0); obsExclMinusSysErrxSecProspino->RemovePoint(0); }
         if( topoStr == "T1bbbb" ){ obsExclMinusSysErrxSecProspino->RemovePoint(0); obsExclMinusSysErrxSecProspino->RemovePoint(0); }
*/

/*
         if( T2ttBR !=0 ){
            obsExclOneTimesxSecProspino->RemovePoint(0); obsExclPlusSysErrxSecProspino->RemovePoint(0); obsExclMinusSysErrxSecProspino->RemovePoint(0);
         }else if(keyStr == "comb" || keyStr == "sbottom" ){
            obsExclMinusSysErrxSecProspino->SetPoint(0, 100, 5.5); obsExclMinusSysErrxSecProspino->SetPoint(obsExclMinusSysErrxSecProspino->GetN(), 115, 18.0);
            obsExclOneTimesxSecProspino->SetPoint(0, 100, 5.5); obsExclOneTimesxSecProspino->SetPoint(obsExclOneTimesxSecProspino->GetN(), 115, 18.0);
            obsExclPlusSysErrxSecProspino->SetPoint(0, 100, 5.5);
         }
         if( T2ttBR ==0 && keyStr == "stop" ){
            obsExclOneTimesxSecProspino->RemovePoint(0); obsExclMinusSysErrxSecProspino->RemovePoint(0); obsExclPlusSysErrxSecProspino->RemovePoint(0);
//            Smooth(obsExclOneTimesxSecProspino, 4); Smooth(obsExclPlusSysErrxSecProspino, 4); Smooth(obsExclMinusSysErrxSecProspino, 4);
            obsExclOneTimesxSecProspino->GetPoint(0, tmpxgr, tmpygr); obsExclOneTimesxSecProspino->SetPoint(obsExclOneTimesxSecProspino->GetN(), tmpxgr, tmpygr);
            obsExclPlusSysErrxSecProspino->GetPoint(0, tmpxgr, tmpygr); obsExclPlusSysErrxSecProspino->SetPoint(obsExclPlusSysErrxSecProspino->GetN(), tmpxgr, tmpygr);
            obsExclMinusSysErrxSecProspino->GetPoint(0, tmpxgr, tmpygr); obsExclMinusSysErrxSecProspino->SetPoint(obsExclMinusSysErrxSecProspino->GetN(), tmpxgr, tmpygr);
         }
*/
         if( keyStr == "fill" ){
            for(unsigned int ir=0; ir<47; ir++){ obsExclOneTimesxSecProspino->RemovePoint(0); }
            for(unsigned int ir=0; ir<47; ir++){ obsExclPlusSysErrxSecProspino->RemovePoint(0); }
            for(unsigned int ir=0; ir<43; ir++){ obsExclMinusSysErrxSecProspino->RemovePoint(0); }
         }

         if( keyStr == "adjT2tt" ){
            double ptx, pty;
            obsExclOneTimesxSecProspino->GetPoint(0, ptx, pty);
            obsExclOneTimesxSecProspino->SetPoint(0, ptx, -1);
            obsExclOneTimesxSecProspino->GetPoint(obsExclOneTimesxSecProspino->GetN()-1, ptx, pty);
            obsExclOneTimesxSecProspino->SetPoint(obsExclOneTimesxSecProspino->GetN(), ptx, -1);

            obsExclPlusSysErrxSecProspino->GetPoint(0, ptx, pty);
            obsExclPlusSysErrxSecProspino->SetPoint(0, ptx, -1);
            obsExclPlusSysErrxSecProspino->GetPoint(obsExclPlusSysErrxSecProspino->GetN()-1, ptx, pty);
            obsExclPlusSysErrxSecProspino->SetPoint(obsExclPlusSysErrxSecProspino->GetN(), ptx, -1);

            obsExclMinusSysErrxSecProspino->GetPoint(0, ptx, pty);
            obsExclMinusSysErrxSecProspino->SetPoint(0, ptx, -1);
            obsExclMinusSysErrxSecProspino->GetPoint(obsExclMinusSysErrxSecProspino->GetN()-1, ptx, pty);
            obsExclMinusSysErrxSecProspino->SetPoint(obsExclMinusSysErrxSecProspino->GetN(), ptx, -1);

            obsExclMinusSysErrxSecProspino_extra->RemovePoint(0);
            obsExclOneTimesxSecProspino_extra->RemovePoint(0);

            obsExclOneTimesxSecProspino_base->GetPoint(obsExclOneTimesxSecProspino_base->GetN()-1, ptx, pty);
            obsExclOneTimesxSecProspino_base->SetPoint(obsExclOneTimesxSecProspino_base->GetN(), ptx, -1);
//            obsExclOneTimesxSecProspino_extra->GetPoint(obsExclOneTimesxSecProspino_extra->GetN()-1, ptx, pty);
//            obsExclOneTimesxSecProspino_extra->SetPoint(obsExclOneTimesxSecProspino_extra->GetN(), ptx, -1);

            obsExclPlusSysErrxSecProspino_base->GetPoint(obsExclPlusSysErrxSecProspino_base->GetN()-1, ptx, pty);
            obsExclPlusSysErrxSecProspino_base->SetPoint(obsExclPlusSysErrxSecProspino_base->GetN(), ptx, -1);
//            obsExclPlusSysErrxSecProspino_extra->GetPoint(0, ptx, pty);
//            obsExclPlusSysErrxSecProspino_extra->SetPoint(0, ptx, -1);
//            obsExclPlusSysErrxSecProspino_extra->GetPoint(obsExclPlusSysErrxSecProspino_extra->GetN()-1, ptx, pty);
//            obsExclPlusSysErrxSecProspino_extra->SetPoint(obsExclPlusSysErrxSecProspino_extra->GetN(), ptx, -1);

            obsExclMinusSysErrxSecProspino_base->GetPoint(obsExclMinusSysErrxSecProspino_base->GetN()-1, ptx, pty);
            obsExclMinusSysErrxSecProspino_base->SetPoint(obsExclMinusSysErrxSecProspino_base->GetN(), ptx, -1);
//            obsExclMinusSysErrxSecProspino_extra->GetPoint(obsExclMinusSysErrxSecProspino_extra->GetN()-1, ptx, pty);
//            obsExclMinusSysErrxSecProspino_extra->SetPoint(obsExclMinusSysErrxSecProspino_extra->GetN(), ptx, -1);
            obsExclMinusSysErrxSecProspino_extra->GetPoint(0, ptx, pty);
            obsExclMinusSysErrxSecProspino_extra->SetPoint(0, ptx, -1);
         }

         obsHistName = region[i]+"obsExclOneTimesProspino"; obsExclOneTimesProspinoVec[i] = (TGraph*) obsExclOneTimesxSecProspino->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclOneTimesProspino_base"; obsExclOneTimesProspino_baseVec[i] = (TGraph*) obsExclOneTimesxSecProspino_base->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclOneTimesProspino_extra"; obsExclOneTimesProspino_extraVec[i] = (TGraph*) obsExclOneTimesxSecProspino_extra->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclThreeTimesProspino"; obsExclThreeTimesProspinoVec[i] = (TGraph*) obsExclThreeTimesxSecProspino->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclOneThirdProspino"; obsExclOneThirdProspinoVec[i] = (TGraph*) obsExclOneThirdxSecProspino->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclPlusSysErrProspino"; obsExclPlusSysErrProspinoVec[i] = (TGraph*) obsExclPlusSysErrxSecProspino->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclPlusSysErrProspino_base"; obsExclPlusSysErrProspino_baseVec[i] = (TGraph*) obsExclPlusSysErrxSecProspino_base->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclPlusSysErrProspino_extra"; obsExclPlusSysErrProspino_extraVec[i] = (TGraph*) obsExclPlusSysErrxSecProspino_extra->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclMinusSysErrProspino"; obsExclMinusSysErrProspinoVec[i] = (TGraph*) obsExclMinusSysErrxSecProspino->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclMinusSysErrProspino_base"; obsExclMinusSysErrProspino_baseVec[i] = (TGraph*) obsExclMinusSysErrxSecProspino_base->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclMinusSysErrProspino_extra"; obsExclMinusSysErrProspino_extraVec[i] = (TGraph*) obsExclMinusSysErrxSecProspino_extra->Clone(obsHistName.Data());

         minObsExclNLOcloned->SetMinimum(0.0001);
         minObsExclNLOcloned->SetMaximum(20);
         minObsExclNLOcloned->Draw("colz");
         diagonalGraph->Draw("l");

         if( topoStr != "TGQ" ){
            obsExclOneTimesxSecProspino->SetLineWidth(2); obsExclOneTimesxSecProspino->Draw("l");
//            obsExclThreeTimesxSecProspino->SetLineWidth(2); obsExclThreeTimesxSecProspino->SetLineStyle(2); obsExclThreeTimesxSecProspino->Draw("l");
//            obsExclOneThirdxSecProspino->SetLineWidth(2); obsExclOneThirdxSecProspino->SetLineStyle(3); obsExclOneThirdxSecProspino->Draw("l");
            obsExclPlusSysErrxSecProspino->SetLineWidth(1); obsExclPlusSysErrxSecProspino->SetLineColor(kMagenta); obsExclPlusSysErrxSecProspino->SetLineStyle(1); obsExclPlusSysErrxSecProspino->Draw("l");
            obsExclMinusSysErrxSecProspino->SetLineWidth(1); obsExclMinusSysErrxSecProspino->SetLineColor(kMagenta); obsExclMinusSysErrxSecProspino->SetLineStyle(1); obsExclMinusSysErrxSecProspino->Draw("l");
         }
/*
         cmsPreTex.DrawLatex(0.20, 0.84, "CMS Preliminary");
         cmsLumiTex.DrawLatex(0.20, 0.77, "4.98 fb^{-1}, #sqrt{s} = 7 TeV");
         disptTex.DrawLatex(0.20, 0.70, "Jets+#slash{H}_{T}");
*/
         indexTex.DrawLatex(0.25, 0.55, "(b)");

         legexp->Clear();
         legexp->AddEntry(obsExclOneTimesxSecProspino, "#sigma^{prod}=#sigma^{NLO-QCD}", "l");
//         legexp->AddEntry(obsExclThreeTimesxSecProspino, "#sigma^{prod}=3#times#sigma^{NLO-QCD}", "l");
//         legexp->AddEntry(obsExclOneThirdxSecProspino, "#sigma^{prod}=1/3#times#sigma^{NLO-QCD}", "l");
         if( topoStr != "TGQ" ) legexp->Draw();

         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ObsLimit_mMother_mLSP.pdf");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ObsLimit_mMother_mLSP.png");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ObsLimit_mMother_mLSP.C");
   
//         radiusToSmooth = 7;
//         dorebin = 0;
   // Expected Limit in M0 - M1/2
         TString expHistName = region[i]+"expLimit";
         TH2F * minExpExclNLOcloned = (TH2F*) minExpExclNLO->Clone(); minExpExclNLOcloned->Reset();
         for(int ix=0; ix<=minExpExclNLO->GetXaxis()->GetNbins(); ix++){
            if( ix < xBinLargery ) continue;
            for(int iy=0; iy<=minExpExclNLO->GetYaxis()->GetNbins(); iy++){
//               if( ix-iy < xBinLargery ) continue;
               minExpExclNLOcloned->SetBinContent(ix, iy, minExpExclNLO->GetBinContent(ix, iy));
            }
         }

         double expOneTimesCutOffX = -1, expThreeTimesCutOffX = -1, expOneThirdCutOffX = -1;
         double expPlusOneSigmaCutOffX = -1, expMinusOneSigmaCutOffX = -1;

         TH2F * expExclOneTimesxSecProspinoSmooth = new TH2F("expExclOneTimesxSecProspinoSmooth", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[i]->Area(expExclOneTimesxSecProspinoSmooth, Mzero, Mhalf, NLOExpxSecCL, false);
         if( dointerpolation ){
            expExclOneTimesxSecProspinoSmooth = interpolate(expExclOneTimesxSecProspinoSmooth); minExpExclNLOcloned = interpolate(minExpExclNLOcloned);
         }
         for(int ir=0; ir<dorebin; ir++){ 
            expExclOneTimesxSecProspinoSmooth = rebin(expExclOneTimesxSecProspinoSmooth); //minExpExclNLOcloned = rebin(minExpExclNLOcloned);
         }
         if( smoothMethod == 0 ){
            expExclOneTimesxSecProspinoSmooth->Smooth(1, "k5a"); minExpExclNLOcloned->Smooth(1, "k5a");
         }else if( smoothMethod == 1 ){
            expExclOneTimesxSecProspinoSmooth->Smooth(1, "k3a"); minExpExclNLOcloned->Smooth(1, "k3a");
         }else if( smoothMethod == 2 ){
            Smooth(expExclOneTimesxSecProspinoSmooth, radiusToSmoothTH2); Smooth(minExpExclNLOcloned, radiusToSmoothTH2);
         }

         explimitVec[i] = (TH2F*) minExpExclNLOcloned->Clone(expHistName.Data());

         explimitVec[i]->SetMinimum(0.0001);
         explimitVec[i]->SetMaximum(20);

         TH2F * expExclOneTimesxSecProspinoTest = 0;
         if( !radiusToSmoothTH2 && !dorebin ){
            expExclOneTimesxSecProspinoTest = new TH2F("expExclOneTimesxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
            plotTools[i]->Area(expExclOneTimesxSecProspinoTest, Mzero, Mhalf, NLOExpxSecCL, xSecProspino, 1.0, doMyFill, xBinLargery);
         }else{
            expExclOneTimesxSecProspinoTest = (TH2F*) expExclOneTimesxSecProspinoSmooth->Clone(); expExclOneTimesxSecProspinoTest->Reset();
            plotTools[i]->Area(expExclOneTimesxSecProspinoTest, Mzero, Mhalf, NLOExpxSecCL, xSecProspino, 1.0, doMyFill, xBinLargery, 0, expExclOneTimesxSecProspinoSmooth);
         }
/*
         ixb = expExclOneTimesxSecProspinoTest->GetXaxis()->FindBin(500);
         iyb = expExclOneTimesxSecProspinoTest->GetYaxis()->FindBin(150);
         expExclOneTimesxSecProspinoTest->SetBinContent(ixb, iyb, 1);
         if( dorebin ){
            ixb = expExclOneTimesxSecProspinoTest->GetXaxis()->FindBin(500-12.5);
            iyb = expExclOneTimesxSecProspinoTest->GetYaxis()->FindBin(150);
            expExclOneTimesxSecProspinoTest->SetBinContent(ixb, iyb, 1);
         }
*/
         TGraph * expExclOneTimesxSecProspinoTMP_base = plotTools[i]->GetContour(expExclOneTimesxSecProspinoTest, 3, 1, !doMyFill);
         TGraph * expExclOneTimesxSecProspinoTMP_extra = plotTools[i]->GetContour(expExclOneTimesxSecProspinoTest, 3, 2, !doMyFill);
         expExclOneTimesxSecProspinoTMP_base->GetPoint(expExclOneTimesxSecProspinoTMP_base->GetN()-1, tmpxgr, tmpygr);
//         expExclOneTimesxSecProspinoTMP_base->SetPoint(expExclOneTimesxSecProspinoTMP_base->GetN(), tmpxgr, 0);
         if( doSmooth ) Smooth(expExclOneTimesxSecProspinoTMP_base, radiusToSmooth);
         if( doSmooth ) Smooth(expExclOneTimesxSecProspinoTMP_extra, radiusToSmooth);

         std::vector<double> xPts_expExclOneTimesxSecProspino_base, yPts_expExclOneTimesxSecProspino_base;
         for(int ip=0; ip<expExclOneTimesxSecProspinoTMP_base->GetN(); ip++){
            double xgr=-1, ygr=-1;
            expExclOneTimesxSecProspinoTMP_base->GetPoint(ip, xgr, ygr);
            if( xgr >= expOneTimesCutOffX ){
               if( xgr != expOneTimesCutOffX && xPts_expExclOneTimesxSecProspino_base.empty() ){
                  xPts_expExclOneTimesxSecProspino_base.push_back(expOneTimesCutOffX); yPts_expExclOneTimesxSecProspino_base.push_back(ygr);
               }
               xPts_expExclOneTimesxSecProspino_base.push_back(xgr); yPts_expExclOneTimesxSecProspino_base.push_back(ygr); 
            }
         }
         int nPts_expExclOneTimesxSecProspino_base = (int)xPts_expExclOneTimesxSecProspino_base.size();
         TGraph *expExclOneTimesxSecProspino_base = new TGraph(nPts_expExclOneTimesxSecProspino_base, &xPts_expExclOneTimesxSecProspino_base[0], &yPts_expExclOneTimesxSecProspino_base[0]);

         std::vector<double> xPts_expExclOneTimesxSecProspino_extra, yPts_expExclOneTimesxSecProspino_extra;
         for(int ip=0; ip<expExclOneTimesxSecProspinoTMP_extra->GetN(); ip++){
            double xgr=-1, ygr=-1;
            expExclOneTimesxSecProspinoTMP_extra->GetPoint(ip, xgr, ygr);
            if( xgr >= expOneTimesCutOffX ){
               if( xgr != expOneTimesCutOffX && xPts_expExclOneTimesxSecProspino_extra.empty() ){
                  xPts_expExclOneTimesxSecProspino_extra.push_back(expOneTimesCutOffX); yPts_expExclOneTimesxSecProspino_extra.push_back(ygr);
               }
               xPts_expExclOneTimesxSecProspino_extra.push_back(xgr); yPts_expExclOneTimesxSecProspino_extra.push_back(ygr); 
            }
         }
         int nPts_expExclOneTimesxSecProspino_extra = (int)xPts_expExclOneTimesxSecProspino_extra.size();
         TGraph *expExclOneTimesxSecProspino_extra = new TGraph(nPts_expExclOneTimesxSecProspino_extra, &xPts_expExclOneTimesxSecProspino_extra[0], &yPts_expExclOneTimesxSecProspino_extra[0]);

         TGraph *expExclOneTimesxSecProspino = MakeBand(expExclOneTimesxSecProspino_base, expExclOneTimesxSecProspino_extra);
//         TGraph *expExclOneTimesxSecProspino = expExclOneTimesxSecProspino_base;

         expExclOneTimesxSecProspinoTestComb = (TH2F*) expExclOneTimesxSecProspinoTest->Clone("expExclOneTimesxSecTestComb");

         expExclOneTimesxSecProspinoTest->Draw("colz");
         expExclOneTimesxSecProspino->Draw("same");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"_Exp_TestContour.pdf");

// +1 Sigma
         TH2F * expExclPlusOneSigmaProspinoSmooth = new TH2F("expExclPlusOneSigmaProspinoSmooth", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[i]->Area(expExclPlusOneSigmaProspinoSmooth, Mzero, Mhalf, NLOExpxSecp1sigma, false);
         if( dointerpolation ){
            expExclPlusOneSigmaProspinoSmooth = interpolate(expExclPlusOneSigmaProspinoSmooth);
         }
         for(int ir=0; ir<dorebin; ir++){ 
            expExclPlusOneSigmaProspinoSmooth = rebin(expExclPlusOneSigmaProspinoSmooth);
         }
         if( smoothMethod == 0 ){
            expExclPlusOneSigmaProspinoSmooth->Smooth(1, "k5a");
         }else if( smoothMethod == 1 ){
            expExclPlusOneSigmaProspinoSmooth->Smooth(1, "k3a");
         }else if( smoothMethod == 2 ){
            Smooth(expExclPlusOneSigmaProspinoSmooth, radiusToSmoothTH2);
         }
         
         TH2F * expExclPlusOneSigmaProspinoTest = 0;
         if( !radiusToSmoothTH2 && !dorebin ){
            expExclPlusOneSigmaProspinoTest = new TH2F("expExclPlusOneSigmaProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
            plotTools[i]->Area(expExclPlusOneSigmaProspinoTest, Mzero, Mhalf, NLOExpxSecp1sigma, xSecProspino, 1.0, doMyFill, xBinLargery);
         }else{   
            expExclPlusOneSigmaProspinoTest = (TH2F*) expExclPlusOneSigmaProspinoSmooth->Clone(); expExclPlusOneSigmaProspinoTest->Reset();
            plotTools[i]->Area(expExclPlusOneSigmaProspinoTest, Mzero, Mhalf, NLOExpxSecp1sigma, xSecProspino, 1.0, doMyFill, xBinLargery, 0, expExclPlusOneSigmaProspinoSmooth);
         }
         
         TGraph * expExclPlusOneSigmaProspinoTMP_base = plotTools[i]->GetContour(expExclPlusOneSigmaProspinoTest, 3, 1, !doMyFill);
         TGraph * expExclPlusOneSigmaProspinoTMP_extra = plotTools[i]->GetContour(expExclPlusOneSigmaProspinoTest, 2, 2, !doMyFill);
         expExclPlusOneSigmaProspinoTMP_base->GetPoint(expExclPlusOneSigmaProspinoTMP_base->GetN()-1, tmpxgr, tmpygr);
//         expExclPlusOneSigmaProspinoTMP_base->SetPoint(expExclPlusOneSigmaProspinoTMP_base->GetN(), tmpxgr, 0);
         if( doSmooth ) Smooth(expExclPlusOneSigmaProspinoTMP_base, radiusToSmooth+2);
         if( doSmooth ) Smooth(expExclPlusOneSigmaProspinoTMP_extra, radiusToSmooth+2);

         std::vector<double> xPts_expExclPlusOneSigmaProspino_base, yPts_expExclPlusOneSigmaProspino_base;
         for(int ip=0; ip<expExclPlusOneSigmaProspinoTMP_base->GetN(); ip++){
            double xgr=-1, ygr=-1;
            expExclPlusOneSigmaProspinoTMP_base->GetPoint(ip, xgr, ygr);
            if( xgr >= expPlusOneSigmaCutOffX ){
               if( xgr != expPlusOneSigmaCutOffX && xPts_expExclPlusOneSigmaProspino_base.empty() ){
                  xPts_expExclPlusOneSigmaProspino_base.push_back(expPlusOneSigmaCutOffX); yPts_expExclPlusOneSigmaProspino_base.push_back(ygr);
               }
               xPts_expExclPlusOneSigmaProspino_base.push_back(xgr); yPts_expExclPlusOneSigmaProspino_base.push_back(ygr); 
            }
         }
         int nPts_expExclPlusOneSigmaProspino_base = (int)xPts_expExclPlusOneSigmaProspino_base.size();
         TGraph *expExclPlusOneSigmaProspino_base = new TGraph(nPts_expExclPlusOneSigmaProspino_base, &xPts_expExclPlusOneSigmaProspino_base[0], &yPts_expExclPlusOneSigmaProspino_base[0]);

         std::vector<double> xPts_expExclPlusOneSigmaProspino_extra, yPts_expExclPlusOneSigmaProspino_extra;
         for(int ip=0; ip<expExclPlusOneSigmaProspinoTMP_extra->GetN(); ip++){
            double xgr=-1, ygr=-1;
            expExclPlusOneSigmaProspinoTMP_extra->GetPoint(ip, xgr, ygr);
            if( xgr >= expPlusOneSigmaCutOffX ){
               if( xgr != expPlusOneSigmaCutOffX && xPts_expExclPlusOneSigmaProspino_extra.empty() ){
                  xPts_expExclPlusOneSigmaProspino_extra.push_back(expPlusOneSigmaCutOffX); yPts_expExclPlusOneSigmaProspino_extra.push_back(ygr);
               }
               xPts_expExclPlusOneSigmaProspino_extra.push_back(xgr); yPts_expExclPlusOneSigmaProspino_extra.push_back(ygr); 
            }
         }
         int nPts_expExclPlusOneSigmaProspino_extra = (int)xPts_expExclPlusOneSigmaProspino_extra.size();
         TGraph *expExclPlusOneSigmaProspino_extra = new TGraph(nPts_expExclPlusOneSigmaProspino_extra, &xPts_expExclPlusOneSigmaProspino_extra[0], &yPts_expExclPlusOneSigmaProspino_extra[0]);

         TGraph *expExclPlusOneSigmaProspino = MakeBand(expExclPlusOneSigmaProspino_base, expExclPlusOneSigmaProspino_extra);
/*
         if( topoStr == "T1" || topoStr == "T1qqqq" ){ expExclPlusOneSigmaProspino->RemovePoint(0); expExclPlusOneSigmaProspino->RemovePoint(0); }
         if( topoStr == "T1bbbb" ){ expExclPlusOneSigmaProspino->RemovePoint(0); expExclPlusOneSigmaProspino->RemovePoint(0); }
*/
         expExclPlusOneSigmaProspinoTest->Draw("colz");
         expExclPlusOneSigmaProspino->Draw("same");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"_ExpPlusOneSigma_TestContour.pdf");

// -1 Sigma
         TH2F * expExclMinusOneSigmaProspinoSmooth = new TH2F("expExclMinusOneSigmaProspinoSmooth", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[i]->Area(expExclMinusOneSigmaProspinoSmooth, Mzero, Mhalf, NLOExpxSecm1sigma, false);
         if( dointerpolation ){
            expExclMinusOneSigmaProspinoSmooth = interpolate(expExclMinusOneSigmaProspinoSmooth);
         }
         for(int ir=0; ir<dorebin; ir++){ 
            expExclMinusOneSigmaProspinoSmooth = rebin(expExclMinusOneSigmaProspinoSmooth);
         }
         if( smoothMethod == 0 ){
            expExclMinusOneSigmaProspinoSmooth->Smooth(1, "k5a");
         }else if( smoothMethod == 1 ){
            expExclMinusOneSigmaProspinoSmooth->Smooth(1, "k3a");
         }else if( smoothMethod == 2 ){
            Smooth(expExclMinusOneSigmaProspinoSmooth, radiusToSmoothTH2);
         }

         TH2F * expExclMinusOneSigmaProspinoTest = 0;
         if( !radiusToSmoothTH2 && !dorebin ){
            expExclMinusOneSigmaProspinoTest = new TH2F("expExclMinusOneSigmaProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
            plotTools[i]->Area(expExclMinusOneSigmaProspinoTest, Mzero, Mhalf, NLOExpxSecm1sigma, xSecProspino, 1.0, doMyFill, xBinLargery);
         }else{
            expExclMinusOneSigmaProspinoTest = (TH2F*) expExclMinusOneSigmaProspinoSmooth->Clone(); expExclMinusOneSigmaProspinoTest->Reset();
            plotTools[i]->Area(expExclMinusOneSigmaProspinoTest, Mzero, Mhalf, NLOExpxSecm1sigma, xSecProspino, 1.0, doMyFill, xBinLargery, 0, expExclMinusOneSigmaProspinoSmooth);
         }
                  
         TGraph * expExclMinusOneSigmaProspinoTMP_base = plotTools[i]->GetContour(expExclMinusOneSigmaProspinoTest, 3, 1, !doMyFill);
         TGraph * expExclMinusOneSigmaProspinoTMP_extra = plotTools[i]->GetContour(expExclMinusOneSigmaProspinoTest, 3, 2, !doMyFill);
         expExclMinusOneSigmaProspinoTMP_base->GetPoint(expExclMinusOneSigmaProspinoTMP_base->GetN()-1, tmpxgr, tmpygr);
//         expExclMinusOneSigmaProspinoTMP->SetPoint(expExclMinusOneSigmaProspinoTMP->GetN(), tmpxgr, 0);
         if( doSmooth ) Smooth(expExclMinusOneSigmaProspinoTMP_base, radiusToSmooth+10);
         if( doSmooth ) Smooth(expExclMinusOneSigmaProspinoTMP_extra, radiusToSmooth+10);

         std::vector<double> xPts_expExclMinusOneSigmaProspino_base, yPts_expExclMinusOneSigmaProspino_base;
         for(int ip=0; ip<expExclMinusOneSigmaProspinoTMP_base->GetN(); ip++){
            double xgr=-1, ygr=-1;
            expExclMinusOneSigmaProspinoTMP_base->GetPoint(ip, xgr, ygr);
            if( xgr >= expMinusOneSigmaCutOffX ){
               if( xgr != expMinusOneSigmaCutOffX && xPts_expExclMinusOneSigmaProspino_base.empty() ){
                  xPts_expExclMinusOneSigmaProspino_base.push_back(expMinusOneSigmaCutOffX); yPts_expExclMinusOneSigmaProspino_base.push_back(ygr);
               }
               xPts_expExclMinusOneSigmaProspino_base.push_back(xgr); yPts_expExclMinusOneSigmaProspino_base.push_back(ygr); 
            }
         }
         int nPts_expExclMinusOneSigmaProspino_base = (int)xPts_expExclMinusOneSigmaProspino_base.size();
         TGraph *expExclMinusOneSigmaProspino_base = new TGraph(nPts_expExclMinusOneSigmaProspino_base, &xPts_expExclMinusOneSigmaProspino_base[0], &yPts_expExclMinusOneSigmaProspino_base[0]);

         std::vector<double> xPts_expExclMinusOneSigmaProspino_extra, yPts_expExclMinusOneSigmaProspino_extra;
         for(int ip=0; ip<expExclMinusOneSigmaProspinoTMP_extra->GetN(); ip++){
            double xgr=-1, ygr=-1;
            expExclMinusOneSigmaProspinoTMP_extra->GetPoint(ip, xgr, ygr);
            if( xgr >= expMinusOneSigmaCutOffX ){
               if( xgr != expMinusOneSigmaCutOffX && xPts_expExclMinusOneSigmaProspino_extra.empty() ){
                  xPts_expExclMinusOneSigmaProspino_extra.push_back(expMinusOneSigmaCutOffX); yPts_expExclMinusOneSigmaProspino_extra.push_back(ygr);
               }
               xPts_expExclMinusOneSigmaProspino_extra.push_back(xgr); yPts_expExclMinusOneSigmaProspino_extra.push_back(ygr);
            }
         }
         int nPts_expExclMinusOneSigmaProspino_extra = (int)xPts_expExclMinusOneSigmaProspino_extra.size();
         TGraph *expExclMinusOneSigmaProspino_extra = new TGraph(nPts_expExclMinusOneSigmaProspino_extra, &xPts_expExclMinusOneSigmaProspino_extra[0], &yPts_expExclMinusOneSigmaProspino_extra[0]);

         TGraph *expExclMinusOneSigmaProspino = MakeBand(expExclMinusOneSigmaProspino_base, expExclMinusOneSigmaProspino_extra);
/*
         if( topoStr == "T2" || topoStr == "T2qq" ){ expExclMinusOneSigmaProspino->RemovePoint(0); expExclMinusOneSigmaProspino->RemovePoint(0); }
         if( topoStr == "T5ZZInc" ){ expExclMinusOneSigmaProspino->RemovePoint(0); }
         if( topoStr == "T1bbbb" ){ expExclMinusOneSigmaProspino->RemovePoint(0); expExclMinusOneSigmaProspino->RemovePoint(0); }
*/

         expExclMinusOneSigmaProspinoTest->Draw("colz");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"_ExpMinusOneSigma_TestContour.pdf");


         TH2F * expExclThreeTimesxSecProspinoTest = new TH2F("expExclThreeTimesxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[i]->Area(expExclThreeTimesxSecProspinoTest, Mzero, Mhalf, NLOExpxSecCL, xSecProspino, 3.0, doMyFill, xBinLargery);
         TGraph * expExclThreeTimesxSecProspinoTMP = plotTools[i]->GetContour(expExclThreeTimesxSecProspinoTest, 3, 1, !doMyFill);
         if( doSmooth ) Smooth(expExclThreeTimesxSecProspinoTMP, radiusToSmooth);

         std::vector<double> xPts_expExclThreeTimesxSecProspino, yPts_expExclThreeTimesxSecProspino;
         for(int ip=0; ip<expExclThreeTimesxSecProspinoTMP->GetN(); ip++){
            double xgr=-1, ygr=-1;
            expExclThreeTimesxSecProspinoTMP->GetPoint(ip, xgr, ygr);
            if( xgr >= expThreeTimesCutOffX ){
               if( xgr != expThreeTimesCutOffX && xPts_expExclThreeTimesxSecProspino.empty() ){
                  xPts_expExclThreeTimesxSecProspino.push_back(expThreeTimesCutOffX); yPts_expExclThreeTimesxSecProspino.push_back(ygr);
               }
               xPts_expExclThreeTimesxSecProspino.push_back(xgr); yPts_expExclThreeTimesxSecProspino.push_back(ygr); 
            }
         }
         int nPts_expExclThreeTimesxSecProspino = (int)xPts_expExclThreeTimesxSecProspino.size();
         TGraph *expExclThreeTimesxSecProspino = new TGraph(nPts_expExclThreeTimesxSecProspino, &xPts_expExclThreeTimesxSecProspino[0], &yPts_expExclThreeTimesxSecProspino[0]);
//         if( topoStr == "T2" || topoStr == "T2qq" ){ expExclThreeTimesxSecProspino->RemovePoint(0); }

         TH2F * expExclOneThirdxSecProspinoTest = new TH2F("expExclOneThirdxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[i]->Area(expExclOneThirdxSecProspinoTest, Mzero, Mhalf, NLOExpxSecCL, xSecProspino, 1./3.0, doMyFill, xBinLargery);
         TGraph * expExclOneThirdxSecProspinoTMP = plotTools[i]->GetContour(expExclOneThirdxSecProspinoTest, 3, 1, !doMyFill);
         if( doSmooth ) Smooth(expExclOneThirdxSecProspinoTMP, radiusToSmooth);

         std::vector<double> xPts_expExclOneThirdxSecProspino, yPts_expExclOneThirdxSecProspino;
         for(int ip=0; ip<expExclOneThirdxSecProspinoTMP->GetN(); ip++){
            double xgr=-1, ygr=-1;
            expExclOneThirdxSecProspinoTMP->GetPoint(ip, xgr, ygr);
            if( xgr >= expOneThirdCutOffX ){
               if( xgr != expOneThirdCutOffX && xPts_expExclOneThirdxSecProspino.empty() ){
                  xPts_expExclOneThirdxSecProspino.push_back(expOneThirdCutOffX); yPts_expExclOneThirdxSecProspino.push_back(ygr);
               }
               xPts_expExclOneThirdxSecProspino.push_back(xgr); yPts_expExclOneThirdxSecProspino.push_back(ygr); 
            }
         }
         int nPts_expExclOneThirdxSecProspino = (int)xPts_expExclOneThirdxSecProspino.size();
         TGraph *expExclOneThirdxSecProspino = new TGraph(nPts_expExclOneThirdxSecProspino, &xPts_expExclOneThirdxSecProspino[0], &yPts_expExclOneThirdxSecProspino[0]);

         if( keyStr == "adjT2tt" ){
            double ptx, pty;
            expExclPlusOneSigmaProspino_extra->RemovePoint(0);
            expExclPlusOneSigmaProspino_base->RemovePoint(0);

            expExclPlusOneSigmaProspino_base->GetPoint(0, ptx, pty);
            expExclPlusOneSigmaProspino_base->SetPoint(0, ptx, -1);
            expExclPlusOneSigmaProspino_base->GetPoint(expExclPlusOneSigmaProspino_base->GetN()-1, ptx, pty);
            expExclPlusOneSigmaProspino_base->SetPoint(expExclPlusOneSigmaProspino_base->GetN(), ptx, -1);

            expExclMinusOneSigmaProspino->GetPoint(0, ptx, pty);
            expExclMinusOneSigmaProspino->SetPoint(0, ptx, -1);
            expExclMinusOneSigmaProspino->GetPoint(expExclMinusOneSigmaProspino->GetN()-1, ptx, pty);
            expExclMinusOneSigmaProspino->SetPoint(expExclMinusOneSigmaProspino->GetN(), ptx, -1);

            expExclMinusOneSigmaProspino_base->GetPoint(0, ptx, pty);
            expExclMinusOneSigmaProspino_base->SetPoint(0, ptx, -1);
            expExclMinusOneSigmaProspino_base->GetPoint(expExclMinusOneSigmaProspino_base->GetN()-1, ptx, pty);
            expExclMinusOneSigmaProspino_base->SetPoint(expExclMinusOneSigmaProspino_base->GetN(), ptx, -1);

            expExclOneTimesxSecProspino_extra->RemovePoint(0);

            expExclOneTimesxSecProspino_base->GetPoint(expExclOneTimesxSecProspino_base->GetN()-1, ptx, pty);
            expExclOneTimesxSecProspino_base->SetPoint(expExclOneTimesxSecProspino_base->GetN(), ptx, -1); 
            expExclOneTimesxSecProspino_extra->GetPoint(0, ptx, pty);
            expExclOneTimesxSecProspino_extra->SetPoint(0, ptx, -1);
            expExclOneTimesxSecProspino_extra->GetPoint(expExclOneTimesxSecProspino_extra->GetN()-1, ptx, pty);
            expExclOneTimesxSecProspino_extra->SetPoint(expExclOneTimesxSecProspino_extra->GetN(), ptx, -1); 

            expExclOneTimesxSecProspino->GetPoint(0, ptx, pty);
            expExclOneTimesxSecProspino->SetPoint(0, ptx, -1);
            expExclOneTimesxSecProspino->GetPoint(expExclOneTimesxSecProspino->GetN()-1, ptx, pty);
            expExclOneTimesxSecProspino->SetPoint(expExclOneTimesxSecProspino->GetN(), ptx, -1);
         }

         expHistName = region[i]+"expExclOneTimesProspino"; expExclOneTimesProspinoVec[i] = (TGraph*) expExclOneTimesxSecProspino->Clone(expHistName.Data());
         expHistName = region[i]+"expExclOneTimesProspino_base"; expExclOneTimesProspino_baseVec[i] = (TGraph*) expExclOneTimesxSecProspino_base->Clone(expHistName.Data());
         expHistName = region[i]+"expExclOneTimesProspino_extra"; expExclOneTimesProspino_extraVec[i] = (TGraph*) expExclOneTimesxSecProspino_extra->Clone(expHistName.Data());
         expHistName = region[i]+"expExclThreeTimesProspino"; expExclThreeTimesProspinoVec[i] = (TGraph*) expExclThreeTimesxSecProspino->Clone(expHistName.Data());
         expHistName = region[i]+"expExclOneThirdProspino"; expExclOneThirdProspinoVec[i] = (TGraph*) expExclOneThirdxSecProspino->Clone(expHistName.Data());

         expHistName = region[i]+"expExclPlusOneSigmaProspino"; expExclPlusOneSigmaProspinoVec[i] = (TGraph*) expExclPlusOneSigmaProspino->Clone(expHistName.Data());
         expHistName = region[i]+"expExclPlusOneSigmaProspino_base"; expExclPlusOneSigmaProspino_baseVec[i] = (TGraph*) expExclPlusOneSigmaProspino_base->Clone(expHistName.Data());
         expHistName = region[i]+"expExclPlusOneSigmaProspino_extra"; expExclPlusOneSigmaProspino_extraVec[i] = (TGraph*) expExclPlusOneSigmaProspino_extra->Clone(expHistName.Data());
         expHistName = region[i]+"expExclMinusOneSigmaProspino"; expExclMinusOneSigmaProspinoVec[i] = (TGraph*) expExclMinusOneSigmaProspino->Clone(expHistName.Data());
         expHistName = region[i]+"expExclMinusOneSigmaProspino_base"; expExclMinusOneSigmaProspino_baseVec[i] = (TGraph*) expExclMinusOneSigmaProspino_base->Clone(expHistName.Data());
         expHistName = region[i]+"expExclMinusOneSigmaProspino_extra"; expExclMinusOneSigmaProspino_extraVec[i] = (TGraph*) expExclMinusOneSigmaProspino_extra->Clone(expHistName.Data());

         minExpExclNLOcloned->SetMinimum(0.0001);
         minExpExclNLOcloned->SetMaximum(20);
	 if( topoStr == "T1tttt" )
         {
         minExpExclNLOcloned->SetMinimum(0.0001);
         minExpExclNLOcloned->SetMaximum(2);
	 }

         minExpExclNLOcloned->Draw("colz");
         diagonalGraph->Draw("l");

         if( topoStr != "TGQ" ){
            expExclOneTimesxSecProspino->SetLineWidth(2); expExclOneTimesxSecProspino->Draw("l");
//            expExclThreeTimesxSecProspino->SetLineWidth(2); expExclThreeTimesxSecProspino->SetLineStyle(2); expExclThreeTimesxSecProspino->Draw("l");
//            expExclOneThirdxSecProspino->SetLineWidth(2); expExclOneThirdxSecProspino->SetLineStyle(3); expExclOneThirdxSecProspino->Draw("l");

            expExclPlusOneSigmaProspino->SetLineWidth(2); expExclPlusOneSigmaProspino->SetLineStyle(8); expExclPlusOneSigmaProspino->Draw("l");
            expExclMinusOneSigmaProspino->SetLineWidth(2); expExclMinusOneSigmaProspino->SetLineStyle(8); expExclMinusOneSigmaProspino->Draw("l");
         }
/*
         cmsPreTex.DrawLatex(0.17, 0.84, "CMS Preliminary");
         cmsLumiTex.DrawLatex(0.17, 0.79, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
         disptTex.DrawLatex(0.17, 0.74, "Jets+#slash{H}_{T}");
*/
         legexp->Clear();
         legexp->AddEntry(expExclOneTimesxSecProspino, "#sigma^{prod}=#sigma^{NLO-QCD}", "l");
//         legexp->AddEntry(expExclThreeTimesxSecProspino, "#sigma^{prod}=3#times#sigma^{NLO-QCD}", "l");
//         legexp->AddEntry(expExclOneThirdxSecProspino, "#sigma^{prod}=1/3#times#sigma^{NLO-QCD}", "l");
         if( topoStr != "TGQ" ) legexp->Draw();
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ExpLimit_mMother_mLSP.pdf");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ExpLimit_mMother_mLSP.png");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ExpLimit_mMother_mLSP.C");
      }
   }

// Draw the limit curves in one canvas for all regions
   c1->Clear();
   c1->SetLogz(1);
   c1->SetRightMargin ( 0.1 );
   c1->SetTopMargin ( 0.09 );

   TH2F* hexcl = 0;
   if( topoStr == "T1" || topoStr == "T1qqqq" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}#; m(#tilde{q})>>m(#tilde{g})");
   }else if( topoStr == "T2" || topoStr == "T2qq" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{q}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{q}#tilde{q}, #tilde{q}#rightarrow q#tilde{#chi}^{0}#; m(#tilde{g})>>m(#tilde{q})");
   }else if( topoStr == "T2tt" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{t}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
   }else if( topoStr == "T2fbd" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{t}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
   }else if( topoStr == "T2tb" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{t}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; #tilde{t}#rightarrow b#tilde{#chi}^{#pm}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
   }else if( topoStr == "TGQ" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{g}} [GeV]; m_{#tilde{q}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{q}#tilde{g}, #tilde{q}#rightarrow q#tilde{#chi}^{0}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}");
   }else if( topoStr == "T5ZZInc" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}_{2}, #tilde{#chi}^{0}_{2}#rightarrow Z#tilde{#chi}^{0}");
   }else if( topoStr == "T1tttt" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
   }else if( topoStr == "T1ttbb" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
   }else if( topoStr == "T1bbbb" ){
      hexcl = new TH2F("hexcl",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow bb#tilde{#chi}^{0}#; m(#tilde{b})>>m(#tilde{g})");
   }else if( topoStr == "T5ttcc" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
   }else if( topoStr == "T5ttttDM175" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
   }else if( topoStr == "T5tttt_degen" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
   }else if( topoStr == "T6ttWW" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{b}} [GeV]; m_{#tilde{#chi}^{#pm}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
   }

/*
   hexcl->SetTitle("CMS Preliminary");
   hexcl->GetYaxis()->SetTitleOffset(1.3);
   hexcl->GetXaxis()->SetTitleOffset(0.92);
   hexcl->GetYaxis()->SetLabelSize(0.03);
   hexcl->GetXaxis()->SetLabelSize(0.03);
*/
   hexcl->Draw("colz");

   TLegend* leg = new TLegend(0.20,0.60,0.48,0.77,NULL,"brNDC"); 
   leg->SetFillColor(0);leg->SetShadowColor(0);leg->SetFillStyle(4000);leg->SetTextFont(42);leg->SetTextSize(0.020*scaleLegendFont);leg->SetBorderSize(0);

   for(int i = 0; i < nTotPlots; i++){
      if(i < nTotPlots-1){
         TGraph *obsCloned = (TGraph*) obsExclOneTimesProspinoVec[i]->Clone();
         TGraph *obs3Cloned = (TGraph*) obsExclThreeTimesProspinoVec[i]->Clone();
         obsCloned->SetLineWidth(2); obsCloned->SetLineColor(colors[i]); 
         obsCloned->Draw("l same");
         obs3Cloned->SetLineWidth(2); obs3Cloned->SetLineColor(colors[i]); obs3Cloned->SetLineStyle(2);
         obs3Cloned->Draw("l same");
         sprintf(tmpStr, "%s (%s)", dispt[i].Data(), "#sigma^{prod}=#sigma^{NLO-QCD}");        
         leg->AddEntry(obsCloned,tmpStr,"l");
         sprintf(tmpStr, "%s (%s)", dispt[i].Data(), "#sigma^{prod}=3#times#sigma^{NLO-QCD}");
         leg->AddEntry(obs3Cloned,tmpStr,"l");
      }
   }
   cmsPreTex.DrawLatex(0.22, 0.82, "CMS Preliminary");
   cmsLumiTex.DrawLatex(0.22, 0.78, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
   disptTex.DrawLatex(0.22, 0.74, "Jets+#slash{H}_{T} Observed Exclusion");
/*
   cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
   cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
   disptTex.DrawLatex(350, 960, "Jets+#slash{H}_{T} Observed Exclusion");
*/
   leg->Draw();
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "comp_Obs_Exclusion_mMother_mLSP.pdf");
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "comp_Obs_Exclusion_mMother_mLSP.png");
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "comp_Obs_Exclusion_mMother_mLSP.C");

   if( topoStr == "T1" || topoStr == "T1qqqq" ){
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}#; m(#tilde{q})>>m(#tilde{g});m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T2" || topoStr == "T2qq" ){
      hexcl->SetTitle("pp#rightarrow#tilde{q}#tilde{q}, #tilde{q}#rightarrow q#tilde{#chi}^{0}#; m(#tilde{g})>>m(#tilde{q});m_{#tilde{q}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T2tt" ){
      hexcl->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0});m_{#tilde{t}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T2fbd" ){
      hexcl->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0});m_{#tilde{t}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T2tb" ){
      hexcl->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}#; #tilde{t}#rightarrow b#tilde{#chi}^{#pm}#; m(#tilde{t})>m(#tilde{#chi}^{0});m_{#tilde{t}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "TGQ" ){
      hexcl->SetTitle("pp#rightarrow#tilde{q}#tilde{g}, #tilde{q}#rightarrow q#tilde{#chi}^{0}, #tilde{g}#rightarrow qq#tilde{#chi}^{0};m_{#tilde{g}} [GeV]; m_{#tilde{q}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T5ZZInc" ){
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}_{2}, #tilde{#chi}^{0}_{2}#rightarrow Z#tilde{#chi}^{0}; m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T1tttt" ){
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g}); m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T1ttbb" ){
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g}); m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T1bbbb" ){
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow bb#tilde{#chi}^{0}#; m(#tilde{b})>>m(#tilde{g})");
   }else if( topoStr == "T5ttcc" ){
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g}); m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T5ttttDM175" ){
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g}); m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T5tttt_degen" ){
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g}); m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T6ttWW" ){
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g}); m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }


   hexcl->Draw("colz");
   leg->Clear();

   for(int i = 0; i < nTotPlots; i++){
      if(i < nTotPlots-1){
         TGraph *expCloned = (TGraph*) expExclOneTimesProspinoVec[i]->Clone();
         TGraph *exp3Cloned = (TGraph*) expExclThreeTimesProspinoVec[i]->Clone();
         expCloned->SetLineWidth(2); expCloned->SetLineColor(colors[i]); 
         expCloned->Draw("l same");
         exp3Cloned->SetLineWidth(2); exp3Cloned->SetLineColor(colors[i]); exp3Cloned->SetLineStyle(2);
         exp3Cloned->Draw("l same");
         sprintf(tmpStr, "%s (%s)", dispt[i].Data(), "#sigma^{prod}=#sigma^{NLO-QCD}");        
         leg->AddEntry(expCloned,tmpStr,"l");
         sprintf(tmpStr, "%s (%s)", dispt[i].Data(), "#sigma^{prod}=3#times#sigma^{NLO-QCD}");
         leg->AddEntry(exp3Cloned,tmpStr,"l");
      }
   }
   cmsPreTex.DrawLatex(0.22, 0.85, "CMS Preliminary");
   cmsLumiTex.DrawLatex(0.22, 0.78, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
   disptTex.DrawLatex(0.22, 0.74, "Jets+#slash{H}_{T} Expected Exclusion");
/*
   cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
   cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
   disptTex.DrawLatex(350, 960, "Jets+#slash{H}_{T} Expected Exclusion");
*/
   leg->Draw();
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "comp_Exp_Exclusion_mMother_mLSP.pdf");
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "comp_Exp_Exclusion_mMother_mLSP.png");
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "comp_Exp_Exclusion_mMother_mLSP.C");

   c1->SetLogz(0);
   c1->SetRightMargin ( 0.125 );
   TLegendEntry *entry =0;

   minCatsObsExclNLO->SetMinimum(0);
   minCatsObsExclNLO->SetMaximum(nTotPlots);
   minCatsObsExclNLO->Draw("colz");
   leg->Clear();
// FIXME: have to manually set fill color!?
   entry = leg->AddEntry(region[0], dispt[0], "f"); entry->SetFillColor(7); entry->SetFillStyle(1001);
   entry = leg->AddEntry(region[1], dispt[1], "f"); entry->SetFillColor(3); entry->SetFillStyle(1001);
   cmsPreTex.DrawLatex(0.22, 0.82, "CMS Preliminary");
   cmsLumiTex.DrawLatex(0.22, 0.78, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
   disptTex.DrawLatex(0.22, 0.74, "Jets+#slash{H}_{T}");
/*
   cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
   cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
//   disptTex.DrawLatex(350, 960, "Jets+#slash{H}_{T} 95% CL Observed Limit [pb]");
   disptTex.DrawLatex(350, 960, "Jets+#slash{H}_{T}");
*/
   leg->Draw();
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "min_Obs_Exclusion_mMother_mLSP.pdf");
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "min_Obs_Exclusion_mMother_mLSP.png");
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "min_Obs_Exclusion_mMother_mLSP.C");

   gStyle->SetPaintTextFormat("1.0f");

   minCatsExpExclNLO->SetMinimum(0);
//   minCatsExpExclNLO->SetMaximum(nTotPlots);
   minCatsExpExclNLO->GetYaxis()->SetRangeUser(12.5, hiYLG);
   minCatsExpExclNLO->Draw("text colz");
   leg->Clear();
   int ncolors = gStyle->GetNumberOfColors();
   int ndiv = gStyle->GetNumberContours();
   double zmin = minCatsExpExclNLO->GetMinimum(), zmax = minCatsExpExclNLO->GetMaximum();
   double dz = zmax-zmin;
   double scale = ndiv/dz;
   int theColor1 =int( (int(0.01+(1-zmin)*scale)+0.99)*1.0*double(ncolors)/double(ndiv) ); if( theColor1 > ncolors-1 ) theColor1 = ncolors-1;
   int theColor2 =int( (int(0.01+(2-zmin)*scale)+0.99)*1.0*double(ncolors)/double(ndiv) ); if( theColor2 > ncolors-1 ) theColor2 = ncolors-1;
   int theColor3 =int( (int(0.01+(3-zmin)*scale)+0.99)*1.0*double(ncolors)/double(ndiv) ); if( theColor3 > ncolors-1 ) theColor3 = ncolors-1;
   int theColor4 =int( (int(0.01+(4-zmin)*scale)+0.99)*1.0*double(ncolors)/double(ndiv) ); if( theColor4 > ncolors-1 ) theColor4 = ncolors-1;
   std::cout<<"\nncolors : "<<ncolors<<"  ndiv : "<<ndiv<<"  zmin : "<<zmin<<"  zmax : "<<zmax<<"  dz : "<<dz<<"  scale : "<<scale<<"  theColor1 : "<<theColor1<<"  theColor2 : "<<theColor2<<"  theColor3 : "<<theColor3<<"  theColor4 : "<<theColor4<<std::endl;
// FIXME: have to manually set fill color!?
   entry = leg->AddEntry(region[0], dispt[0], "f"); entry->SetFillColor(gStyle->GetColorPalette(theColor1)); entry->SetFillStyle(1001);
   entry = leg->AddEntry(region[1], dispt[1], "f"); entry->SetFillColor(gStyle->GetColorPalette(theColor2)); entry->SetFillStyle(1001);
   entry = leg->AddEntry(region[2], dispt[2], "f"); entry->SetFillColor(gStyle->GetColorPalette(theColor3)); entry->SetFillStyle(1001);
   entry = leg->AddEntry(region[3], dispt[3], "f"); entry->SetFillColor(gStyle->GetColorPalette(theColor4)); entry->SetFillStyle(1001);
   cmsPreTex.DrawLatex(0.22, 0.85, "CMS Preliminary");
   cmsLumiTex.DrawLatex(0.22, 0.80, "L_{int} = 19.4 fb^{-1}, #sqrt{s} = 8 TeV");

   leg->Draw();
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "min_Exp_Exclusion_mMother_mLSP.pdf");
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "min_Exp_Exclusion_mMother_mLSP.png");
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "min_Exp_Exclusion_mMother_mLSP.C");

   TCanvas * c2 = new TCanvas("c2","c2",600,600);
   c2->cd();
   c2->SetLogz(1);
   c2->SetRightMargin ( 0.15 );
   c2->SetTopMargin ( 0.09 );
   TLatex disptTex4ratio; disptTex4ratio.SetTextSize(0.020*scaleLegendFont); disptTex4ratio.SetTextColor(1);

   c2->Print(outDir + "/" + topoStr + "_" + "ratio_limits_mMother_mLSP.pdf[");
   for(int ip=0; ip<nTotPlots; ip++){
      for(int jp=ip+1; jp<nTotPlots; jp++){
         TH2F *obsRatio2D = (TH2F*) obslimitVec[ip]->Clone(); obsRatio2D->Reset();
         TH2F *expRatio2D = (TH2F*) explimitVec[ip]->Clone(); expRatio2D->Reset();

         TAxis *obsXaxis = (TAxis*) obsRatio2D->GetXaxis(); TAxis *obsYaxis = (TAxis*) obsRatio2D->GetYaxis();
         int obsXnbins = obsXaxis->GetNbins(), obsYnbins = obsYaxis->GetNbins();

         TAxis *expXaxis = (TAxis*) expRatio2D->GetXaxis(); TAxis *expYaxis = (TAxis*) expRatio2D->GetYaxis();
         int expXnbins = expXaxis->GetNbins(), expYnbins = expYaxis->GetNbins();

         for(int ix=0; ix<obsXnbins; ix++){
            for(int iy=0; iy<obsYnbins; iy++){
               double numCont = obslimitVec[jp]->GetBinContent(ix+1, iy+1);
               double denCont = obslimitVec[ip]->GetBinContent(ix+1, iy+1);
               if( denCont !=0 ){
                  obsRatio2D->SetBinContent(ix+1, iy+1, numCont/denCont);
               }
            }
         }
         for(int ix=0; ix<expXnbins; ix++){
            for(int iy=0; iy<expYnbins; iy++){
               double numCont = explimitVec[jp]->GetBinContent(ix+1, iy+1);
               double denCont = explimitVec[ip]->GetBinContent(ix+1, iy+1);
               if( denCont !=0 ){
                  expRatio2D->SetBinContent(ix+1, iy+1, numCont/denCont);
               }
            }
         }
        
         obsRatio2D->SetMinimum(0.02); obsRatio2D->SetMaximum(5.0);
         obsRatio2D->SetMarkerSize(0.5);
         obsRatio2D->Draw("Text colz");
         cmsPreTex.DrawLatex(0.22, 0.82, "CMS Preliminary");
         cmsLumiTex.DrawLatex(0.22, 0.78, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
         disptTex.DrawLatex(0.22, 0.74, "Observed Ratio");
         disptTex4ratio.DrawLatex(0.22, 0.71, "("+dispt[jp]+")"+"/"+"("+dispt[ip]+")");
/*
         cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
         cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
         disptTex.DrawLatex(350, 960, "Observed Ratio");
         disptTex4ratio.DrawLatex(350, 890, "("+dispt[jp]+")"+"/"+"("+dispt[ip]+")");
*/
         c2->Print(outDir + "/" + topoStr + "_" + "ratio_limits_mMother_mLSP.pdf");

         expRatio2D->SetMinimum(0.02); expRatio2D->SetMaximum(5.0);
         expRatio2D->SetMarkerSize(0.5);
         expRatio2D->Draw("Text colz");
         cmsPreTex.DrawLatex(0.22, 0.82, "CMS Preliminary");
         cmsLumiTex.DrawLatex(0.22, 0.78, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
         disptTex.DrawLatex(0.22, 0.74, "Expected Ratio");
         disptTex4ratio.DrawLatex(0.22, 0.71, "("+dispt[jp]+")"+"/"+"("+dispt[ip]+")");
/*
         cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
         cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
         disptTex.DrawLatex(350, 960, "Expected Ratio");
         disptTex4ratio.DrawLatex(350, 890, "("+dispt[jp]+")"+"/"+"("+dispt[ip]+")");
*/
         c2->Print(outDir + "/" + topoStr + "_" + "ratio_limits_mMother_mLSP.pdf");
      }
   }
   c2->Print(outDir + "/" + topoStr + "_" + "ratio_limits_mMother_mLSP.pdf]");

// Define the box and lines

   const Int_t np = 2;

   Double_t xcb[np], ycb[np], ycbmin[np], ycbmax[np];

   xcb[0]=150.; ycbmin[0]=-50.; ycbmax[0]=0.; 

   xcb[1]=262.5; ycbmin[1]=100.; ycbmax[1]=112.5;

   if(topoStr == "T2tb"){
//      ycbmax[0] = 50.; ycbmax[1] = 550.;
      ycbmin[1] = 97.5;
   }

   TGraph *grmin = new TGraph(np,xcb,ycbmin);
   TGraph *grmax = new TGraph(np,xcb,ycbmax);

   diagonalCoverBand = new TGraph(2*np);
   diagonalCoverBand->SetName("diagonalCoverBand");

   for (int ip=0; ip<np; ip++) {
      diagonalCoverBand->SetPoint(ip,xcb[ip],ycbmax[ip]);
      if( ip ==0 )
         diagonalCoverBand->SetPoint(np+ip,xcb[np-ip-1]+37.5,ycbmin[np-ip-1]);
      else 
         diagonalCoverBand->SetPoint(np+ip,xcb[np-ip-1],ycbmin[np-ip-1]);
   }
   diagonalCoverBand->SetFillColor(0);

//   TLatex tl; tl.SetTextFont(43); tl.SetTextSize(18); tl.SetTextAngle(47.);
//   tl.DrawLatex(425.,263.,"m_{#tilde{t}} = m_{t} + m_{#tilde{#chi}_{1}^{0}}");

   sprintf(tmpStr, "CLs_SMS_BR%02dpct.root", int(T2ttBR*100));
   TFile *outRootFile = new TFile(outDir+"/"+tmpStr, "RECREATE");

   diagonalGraph->Write();
   diagonalCoverBand->Write();

   for(int i=0; i<nTotPlots; i++){
      if( i!= nTotPlots -1 ) continue;

      sprintf(tmpStr, "%s_BR%dpct", obslimitVec[i]->GetName(), int(T2ttBR*100)); obslimitVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", explimitVec[i]->GetName(), int(T2ttBR*100)); explimitVec[i]->SetName(tmpStr);
         
      sprintf(tmpStr, "%s_BR%dpct", obsExclOneTimesProspinoVec[i]->GetName(), int(T2ttBR*100)); obsExclOneTimesProspinoVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", obsExclOneTimesProspino_baseVec[i]->GetName(), int(T2ttBR*100)); obsExclOneTimesProspino_baseVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", obsExclOneTimesProspino_extraVec[i]->GetName(), int(T2ttBR*100)); obsExclOneTimesProspino_extraVec[i]->SetName(tmpStr);

      sprintf(tmpStr, "%s_BR%dpct", expExclOneTimesProspinoVec[i]->GetName(), int(T2ttBR*100)); expExclOneTimesProspinoVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", expExclOneTimesProspino_baseVec[i]->GetName(), int(T2ttBR*100)); expExclOneTimesProspino_baseVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", expExclOneTimesProspino_extraVec[i]->GetName(), int(T2ttBR*100)); expExclOneTimesProspino_extraVec[i]->SetName(tmpStr);

      sprintf(tmpStr, "%s_BR%dpct", expExclPlusOneSigmaProspinoVec[i]->GetName(), int(T2ttBR*100)); expExclPlusOneSigmaProspinoVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", expExclPlusOneSigmaProspino_baseVec[i]->GetName(), int(T2ttBR*100)); expExclPlusOneSigmaProspino_baseVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", expExclPlusOneSigmaProspino_extraVec[i]->GetName(), int(T2ttBR*100)); expExclPlusOneSigmaProspino_extraVec[i]->SetName(tmpStr);

      sprintf(tmpStr, "%s_BR%dpct", expExclMinusOneSigmaProspinoVec[i]->GetName(), int(T2ttBR*100)); expExclMinusOneSigmaProspinoVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", expExclMinusOneSigmaProspino_baseVec[i]->GetName(), int(T2ttBR*100)); expExclMinusOneSigmaProspino_baseVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", expExclMinusOneSigmaProspino_extraVec[i]->GetName(), int(T2ttBR*100)); expExclMinusOneSigmaProspino_extraVec[i]->SetName(tmpStr);

      sprintf(tmpStr, "%s_BR%dpct", obsExclPlusSysErrProspinoVec[i]->GetName(), int(T2ttBR*100)); obsExclPlusSysErrProspinoVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", obsExclPlusSysErrProspino_baseVec[i]->GetName(), int(T2ttBR*100)); obsExclPlusSysErrProspino_baseVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", obsExclPlusSysErrProspino_extraVec[i]->GetName(), int(T2ttBR*100)); obsExclPlusSysErrProspino_extraVec[i]->SetName(tmpStr);

      sprintf(tmpStr, "%s_BR%dpct", obsExclMinusSysErrProspinoVec[i]->GetName(), int(T2ttBR*100)); obsExclMinusSysErrProspinoVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", obsExclMinusSysErrProspino_baseVec[i]->GetName(), int(T2ttBR*100)); obsExclMinusSysErrProspino_baseVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", obsExclMinusSysErrProspino_extraVec[i]->GetName(), int(T2ttBR*100)); obsExclMinusSysErrProspino_extraVec[i]->SetName(tmpStr);

      sprintf(tmpStr, "%s_BR%dpct", expExclThreeTimesProspinoVec[i]->GetName(), int(T2ttBR*100)); expExclThreeTimesProspinoVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", expExclOneThirdProspinoVec[i]->GetName(), int(T2ttBR*100)); expExclOneThirdProspinoVec[i]->SetName(tmpStr);

      sprintf(tmpStr, "%s_BR%dpct", obsExclThreeTimesProspinoVec[i]->GetName(), int(T2ttBR*100)); obsExclThreeTimesProspinoVec[i]->SetName(tmpStr);
      sprintf(tmpStr, "%s_BR%dpct", obsExclOneThirdProspinoVec[i]->GetName(), int(T2ttBR*100)); obsExclOneThirdProspinoVec[i]->SetName(tmpStr);

      obslimitVec[i]->Write(); explimitVec[i]->Write();

      if( keyStr != "raw" ){         
         obsExclOneTimesProspinoVec[i]->Write(); obsExclOneTimesProspino_baseVec[i]->Write(); obsExclOneTimesProspino_extraVec[i]->Write();
         //obsExclThreeTimesProspinoVec[i]->Write(); obsExclOneThirdProspinoVec[i]->Write();
         expExclOneTimesProspinoVec[i]->Write(); expExclOneTimesProspino_baseVec[i]->Write(); expExclOneTimesProspino_extraVec[i]->Write();
         //expExclThreeTimesProspinoVec[i]->Write(); expExclOneThirdProspinoVec[i]->Write();

         expExclPlusOneSigmaProspinoVec[i]->Write(); expExclPlusOneSigmaProspino_baseVec[i]->Write(); expExclPlusOneSigmaProspino_extraVec[i]->Write();
         expExclMinusOneSigmaProspinoVec[i]->Write(); expExclMinusOneSigmaProspino_baseVec[i]->Write(); expExclMinusOneSigmaProspino_extraVec[i]->Write();

         obsExclPlusSysErrProspinoVec[i]->Write(); obsExclPlusSysErrProspino_baseVec[i]->Write(); obsExclPlusSysErrProspino_extraVec[i]->Write(); 
         obsExclMinusSysErrProspinoVec[i]->Write(); obsExclMinusSysErrProspino_baseVec[i]->Write(); obsExclMinusSysErrProspino_extraVec[i]->Write();
      }
   }
   outRootFile->Write(); outRootFile->Close();

}

int main(int argc, char** argv)
{
  return plot(argc, argv);
}
