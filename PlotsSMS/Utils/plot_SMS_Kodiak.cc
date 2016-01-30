#include "plot.h"
#include "SusyScan.h"
#include "GeneratorMasses.h"
#include "PlotTools.h"
#include "TheLimits.h"
#include "GlobalFunctions.h"
#include "StyleSettings.h"

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
static const string  disptStrs[] = {"H_{T}>500, #slash{H}_{T}>350 GeV", "H_{T}>800, #slash{H}_{T}>200 GeV", "H_{T}>800, #slash{H}_{T}>500 GeV", "H_{T}>1000, #slash{H}_{T}>400 GeV", "H_{T}>1200, #slash{H}_{T}>200 GeV"};
static const string regionDirs[] = {"mediumHTMHT", "highHT", "highHTMHT", "HT1000MHT400", "HT1200MHT200"};
static const int     inclFlags[] = {     1,           1,          1,           0,               0       };  

static const int      colors[] = {kGreen, kMagenta+1, kBlue, kRed, kTeal+4};

TH1D *xSecProspino_T1, *xSecProspino_T2, *xSecProspino_stop, *xSecProspino;
TString topoStr = "T1";

// 25 GeV per X bin; 25 GeV per Y bin
// last working one for RA2
//int nXbins = 52, nYbins = 40;

// Future SUSY
//int nXbins = 72, nYbins = 52;

//int nXbins = 52, nYbins = 47;
//const double loX = 100, hiX = 1225;
//const double loY = 50, hiY = 1225;

// last working one for RA2
//double loX = 212.5, hiX = 1512.5;
//double loY = 212.5, hiY = 1212.5;


 int nXbins = 52, nYbins = 47;
 double loX = 37.5, hiX = 1512.5;
 double loY = 37.5, hiY = 1212.5;




// Future SUSY
//double loX = 212.5, hiX = 2012.5;
//double loY = 212.5, hiY = 1512.5;

//double loY = 37.5, hiY = 1212.5;
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

int plot(int argc, char** argv)
{
   //interactive root session
   //TApplication theApp("App", 0, 0);
   if (gROOT->IsBatch()) {
      fprintf(stderr, "%s: cannot run in batch mode\n", argv[0]);
      return 1;
   }
   if( argc != 2 ){
      std::cout<<"\n################USAGE##########################"<<std::endl;
      std::cout<<  "#./plot_SMS xxx                               #"<<std::endl;
      std::cout<<  "# xxx: T1 or T2                               #"<<std::endl;
      std::cout<<  "################USAGE##########################"<<std::endl;
      std::cout<<std::endl;
      return 1;
   }

   topoStr = argv[1];
   std::cout<<"\nWorking on topology "<<topoStr<<std::endl<<std::endl;

// For T1, it's reference_gluino_xSec.root
// For T2, it's reference_4squark_xSec.root
// For T2tt, it's reference_stop_xSec.root
   TFile *xSecProspinoFile = new TFile("reference_xSec.root");
   if( topoStr == "T1" || topoStr == "T1qqqq" || topoStr == "T5ZZInc" || topoStr == "T1tttt" || topoStr == "T1bbbb" ) xSecProspino_T1 = (TH1D*) xSecProspinoFile->Get("gluino_xsection");
   if( topoStr == "T2" || topoStr == "T2qq" ) xSecProspino_T2 = (TH1D*) xSecProspinoFile->Get("squark_xsection");
   if( topoStr == "T2tt" ){
      xSecProspino_stop = (TH1D*) xSecProspinoFile->Get("stop_xsection");
   }
   if( topoStr == "T1" || topoStr == "T1qqqq" || topoStr == "T5ZZInc" || topoStr == "T1tttt" || topoStr == "T1bbbb" ) xSecProspino = (TH1D*)xSecProspino_T1->Clone();
   if( topoStr == "T2" || topoStr == "T2qq" ) xSecProspino = (TH1D*)xSecProspino_T2->Clone();
   if( topoStr == "T2tt" ) xSecProspino = (TH1D*)xSecProspino_stop->Clone();
   if( topoStr == "TGQ" ) xSecProspino = (TH1D*)xSecProspino_T2->Clone();

// 25 GeV per X bin; 25 GeV per Y bin
   if( topoStr == "T2tt" ){
      nXbins = 29; nYbins = 17;
      loX = 187.5; hiX = 912.5;
      loY = -12.5; hiY = 412.5;
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
   for(int i=0; i<nColors; i++) MyPalette[i] = FI+i;
   gStyle->SetNumberContours(nColors);

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
      yPtsDiag[ib] = loDiagRange + ib*widthBin; yErrsDiag[ib] = 0.0;
   }
   diagonalGraph = new TGraphErrors(nDiagBins, xPtsDiag, yPtsDiag, xErrsDiag, yErrsDiag);
   diagonalGraph->SetLineWidth(2); diagonalGraph->SetLineStyle(7);
   diagonalGraph->SetName("diagonalGraph");

   vector<TheLimits *> genpoints;
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
   nTotPlots++; region.push_back("combined_"); regionBase.push_back("combined"); dispt.push_back("combined");

   int idxCnt =-1;
   for(int ir=0; ir<nRegions; ir++){
      if( inclFlags[ir] ){
         sprintf(dirname, "%s_plots_%s/", regionDirs[ir].c_str(), topoStr.Data());
         sprintf(filename, "%s_plots_%s/filelist.txt", regionDirs[ir].c_str(), topoStr.Data());
         idxCnt++; genpoints[idxCnt]->Fill(filename, dirname);
//         genpoints[idxCnt]->FillEmptyPointsByInterpolationInM0M12();
//Get SUSY masses and x-sections from generator scan ----------------------------
//and match them to the signal scan
         genpoints[idxCnt]->FillGeneratorMasses("GenScan_tb10.dat");
         genpoints[idxCnt]->match();
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

   TH2F *minCatsObsExclNLO, *minCatsExpExclNLO;
   TH2F *minObsExclNLO, *minExpExclNLO, *minObsExclNLOfollowminExp;

   //Was not originally designed as a loop, so this will give you warnings about possible memory leaks.

//   TLegend* legexp = new TLegend(0.15,0.63,0.48,0.71,NULL,"brNDC"); 
   TLegend* legexp = new TLegend(0.18,0.65,0.48,0.80,NULL,"brNDC"); 
   legexp->SetFillColor(0);legexp->SetShadowColor(0);legexp->SetFillStyle(4000);legexp->SetTextFont(42);legexp->SetTextSize(0.023*scaleLegendFont);legexp->SetBorderSize(0);

   TLatex cmsPreTex; cmsPreTex.SetTextSize(0.025*scaleLegendFont); cmsPreTex.SetTextColor(1); cmsPreTex.SetNDC(); cmsPreTex.SetTextFont(42);
   TLatex cmsLumiTex; cmsLumiTex.SetTextSize(0.025*scaleLegendFont); cmsLumiTex.SetTextColor(1); cmsLumiTex.SetNDC(); cmsLumiTex.SetTextFont(42);
   TLatex disptTex; disptTex.SetTextSize(0.025*scaleLegendFont); disptTex.SetTextColor(1); disptTex.SetNDC(); disptTex.SetTextFont(42);

   TLatex indexTex; indexTex.SetTextSize(0.025*scaleLegendFont); indexTex.SetTextColor(1); indexTex.SetNDC(); indexTex.SetTextFont(42);

   for(int i = 0; i < nTotPlots; i++){
      if(i < nTotPlots-1){
   //the plotting ----------------------------------------------------------------------
   //plotting helper functions
         plotTools[i] = new PlotTools<SusyScan>(genpoints[i]->GetScan());
         PlotTools<GeneratorMasses> * plotMasses = new PlotTools<GeneratorMasses>(genpoints[i]->GetGeneratorMasses());

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
            hobslimit->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow q#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
         }else if( topoStr == "TGQ" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{g}) [GeV]; m(#tilde{q}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{q}#tilde{g}, #tilde{q}#rightarrow q#tilde{#chi}^{0}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}");
         }else if( topoStr == "T5ZZInc" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}_{2}, #tilde{#chi}^{0}_{2}#rightarrow Z#tilde{#chi}^{0}");
         }else if( topoStr == "T1tttt" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }else if( topoStr == "T1bbbb" ){
            hobslimit = new TH2F("obslimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow bb#tilde{#chi}^{0}#; m(#tilde{b})>>m(#tilde{g})");
         }
         plotTools[i]->Area(hobslimit, Mzero, Mhalf, NLOObsxSecCL, false);
         hobslimit->SetMinimum(0.001);
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
            obsExclThreeTimesxSecProspino->SetLineWidth(2); obsExclThreeTimesxSecProspino->SetLineStyle(2); obsExclThreeTimesxSecProspino->Draw("same");
         }
         cmsPreTex.DrawLatex(0.17, 0.84, "CMS Preliminary");
         cmsLumiTex.DrawLatex(0.17, 0.79, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
         disptTex.DrawLatex(0.17, 0.74, dispt[i]);
/*
         cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
         cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
         disptTex.DrawLatex(350, 960, dispt[i]);
*/
         legexp->Clear();
         legexp->AddEntry(obsExclOneTimesxSecProspino, "#sigma^{prod}=#sigma^{NLO-QCD}", "l");
         legexp->AddEntry(obsExclThreeTimesxSecProspino, "#sigma^{prod}=3#times#sigma^{NLO-QCD}", "l");
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
            hexplimit->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow q#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
         }else if( topoStr == "TGQ" ){
            hexplimit = new TH2F("explimit",";m(#tilde{g}) [GeV]; m(#tilde{q}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{q}#tilde{g}, #tilde{q}#rightarrow q#tilde{#chi}^{0}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}");
         }else if( topoStr == "T5ZZInc"){
            hexplimit = new TH2F("explimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}_{2}, #tilde{#chi}^{0}_{2}#rightarrow Z#tilde{#chi}^{0}");
         }else if( topoStr == "T1tttt" ){
            hexplimit = new TH2F("explimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
         }else if( topoStr == "T1bbbb" ){
            hexplimit = new TH2F("explimit",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow bb#tilde{#chi}^{0}#; m(#tilde{b})>>m(#tilde{g})");
         }
         plotTools[i]->Area(hexplimit, Mzero, Mhalf, NLOExpxSecCL, false);
         hexplimit->SetMinimum(0.001);
         hexplimit->SetMaximum(20);
         hexplimit->Draw("colz");
         TString expHistName = region[i]+"expLimit";
         explimitVec[i] = (TH2F*) hexplimit->Clone(expHistName.Data());

         TGraph *expExclOneTimesxSecProspino = set_ProspinoExcl(xSecProspino, (TH2D*) hexplimit, 1.0);
         TGraph *expExclThreeTimesxSecProspino = set_ProspinoExcl(xSecProspino, (TH2D*) hexplimit, 3.0);
         expHistName = region[i]+"expExclOneTimesProspino"; expExclOneTimesProspinoVec[i] = (TGraph*) expExclOneTimesxSecProspino->Clone(expHistName.Data());
         expHistName = region[i]+"expExclThreeTimesProspino"; expExclThreeTimesProspinoVec[i] = (TGraph*) expExclThreeTimesxSecProspino->Clone(expHistName.Data());

         if( topoStr != "TGQ" ){
            expExclOneTimesxSecProspino->SetLineWidth(2); expExclOneTimesxSecProspino->Draw("same");
//            expExclThreeTimesxSecProspino->SetLineWidth(2); expExclThreeTimesxSecProspino->SetLineStyle(2); expExclThreeTimesxSecProspino->Draw("same");
         }
         cmsPreTex.DrawLatex(0.20, 0.85, "CMS Preliminary");
         cmsLumiTex.DrawLatex(0.20, 0.79, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
//         disptTex.DrawLatex(0.17, 0.74, dispt[i]);
/*
         cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
         cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
         disptTex.DrawLatex(350, 960, dispt[i]);
*/
         legexp->Clear();
         legexp->AddEntry(expExclOneTimesxSecProspino, "#sigma^{prod}=#sigma^{NLO-QCD}", "l");
//         legexp->AddEntry(expExclThreeTimesxSecProspino, "#sigma^{prod}=3#times#sigma^{NLO-QCD}", "l");
         if( topoStr != "TGQ" ) legexp->Draw();
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ExpLimit_mMother_mLSP.pdf");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ExpLimit_mMother_mLSP.png");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ExpLimit_mMother_mLSP.C");

         if( i==0 ){
            minCatsObsExclNLO = (TH2F*) obslimitVec[i]->Clone("minCatsObsExclNLO"); minCatsExpExclNLO = (TH2F*) explimitVec[i]->Clone("minCatsExpExclNLO");
            minObsExclNLO = (TH2F*) obslimitVec[i]->Clone("minObsExclNLO"); minExpExclNLO = (TH2F*) explimitVec[i]->Clone("minExpExclNLO");
            minObsExclNLOfollowminExp = (TH2F*) obslimitVec[i]->Clone("minObsExclNLOfollowminExp");
            TAxis *obsXaxis = (TAxis*) minCatsObsExclNLO->GetXaxis(); TAxis *obsYaxis = (TAxis*) minCatsObsExclNLO->GetYaxis();
            int obsXnbins = obsXaxis->GetNbins(), obsYnbins = obsYaxis->GetNbins();
            for(int ix=0; ix<obsXnbins; ix++){
               for(int iy=0; iy<obsYnbins; iy++){
                  if( obsYaxis->GetBinLowEdge(iy+1) > obsXaxis->GetBinLowEdge(ix+1) -50 ) continue;
                  minCatsObsExclNLO->SetBinContent(ix+1, iy+1, 1.0*i+1.0);
               }
            }
            TAxis *expXaxis = (TAxis*) minCatsExpExclNLO->GetXaxis(); TAxis *expYaxis = (TAxis*) minCatsExpExclNLO->GetYaxis();
            int expXnbins = expXaxis->GetNbins(), expYnbins = expYaxis->GetNbins();
            for(int ix=0; ix<expXnbins; ix++){
               for(int iy=0; iy<expYnbins; iy++){
                  if( obsYaxis->GetBinLowEdge(iy+1) > obsXaxis->GetBinLowEdge(ix+1) -50 ) continue;
                  minCatsExpExclNLO->SetBinContent(ix+1, iy+1, 1.0*i+1.0);
               }
            }
         }else{
            TAxis *obsXaxis = (TAxis*) minCatsObsExclNLO->GetXaxis(); TAxis *obsYaxis = (TAxis*) minCatsObsExclNLO->GetYaxis();
            int obsXnbins = obsXaxis->GetNbins(), obsYnbins = obsYaxis->GetNbins();
            for(int ix=0; ix<obsXnbins; ix++){
               for(int iy=0; iy<obsYnbins; iy++){
                  if( obsYaxis->GetBinLowEdge(iy+1) > obsXaxis->GetBinLowEdge(ix+1) -50 ) continue;
                  double curMinVal = minObsExclNLO->GetBinContent(ix+1, iy+1);
                  double curVal = obslimitVec[i]->GetBinContent(ix+1, iy+1);
                  if( curVal < curMinVal ){ minCatsObsExclNLO->SetBinContent(ix+1, iy+1, 1.0*i+1.0); minObsExclNLO->SetBinContent(ix+1, iy+1, curVal); }
               }
            }

            TAxis *expXaxis = (TAxis*) minCatsExpExclNLO->GetXaxis(); TAxis *expYaxis = (TAxis*) minCatsExpExclNLO->GetYaxis();
            int expXnbins = expXaxis->GetNbins(), expYnbins = expYaxis->GetNbins();
            for(int ix=0; ix<expXnbins; ix++){
               for(int iy=0; iy<expYnbins; iy++){
                  if( obsYaxis->GetBinLowEdge(iy+1) > obsXaxis->GetBinLowEdge(ix+1) -50 ) continue;
                  double curMinVal = minExpExclNLO->GetBinContent(ix+1, iy+1);
                  double curVal = explimitVec[i]->GetBinContent(ix+1, iy+1);
                  if( curVal < curMinVal ){
                     minCatsExpExclNLO->SetBinContent(ix+1, iy+1, 1.0*i+1.0); minExpExclNLO->SetBinContent(ix+1, iy+1, curVal); 
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
               if( ix-iy < xBinLargery ) continue;
               minObsExclNLOcloned->SetBinContent(ix, iy, minObsExclNLO->GetBinContent(ix, iy));
            }
         }
//         minObsExclNLOcloned->GetXaxis()->SetRangeUser(300, 1200);
         obslimitVec[i] = (TH2F*) minObsExclNLOcloned->Clone(obsHistName.Data());

//         int radiusToSmooth = 5;
//         int radiusToSmooth = 8;
         int radiusToSmooth = 30;

         if( topoStr == "T5ZZInc" ) radiusToSmooth = 16;
         if( topoStr == "T1bbbb" ) radiusToSmooth = 16;
         if( topoStr == "T1tttt" ) radiusToSmooth = 8;

         double obsOneTimesCutOffX = -1, obsThreeTimesCutOffX = -1, obsOneThirdCutOffX = -1;
         double obsPlusSysErrCutOffX = -1, obsMinusSysErrCutOffX = -1;
/*
         if( topoStr == "T1" || topoStr == "T1qqqq" ){
//            obsOneTimesCutOffX = 537.5, obsThreeTimesCutOffX = 687.5, obsOneThirdCutOffX = 412.5;
            obsOneTimesCutOffX = 535.5, obsThreeTimesCutOffX = 685.5, obsOneThirdCutOffX = 410.5;
//            obsPlusSysErrCutOffX = 562.5, obsMinusSysErrCutOffX = 512.5;
            obsPlusSysErrCutOffX = 560.5, obsMinusSysErrCutOffX = 510.5;
         }
         if( topoStr == "T2" || topoStr == "T2qq" ){
//            obsOneTimesCutOffX = 437.5, obsThreeTimesCutOffX = 562.5, obsOneThirdCutOffX = 312.5;
            obsOneTimesCutOffX = 402.5, obsThreeTimesCutOffX = 560.5, obsOneThirdCutOffX = 310.5;
//            obsPlusSysErrCutOffX = 437.5, obsMinusSysErrCutOffX = 387.5;
            obsPlusSysErrCutOffX = 430.5, obsMinusSysErrCutOffX = 385.5;
         }
         if( topoStr == "T5ZZInc" ){
            obsOneTimesCutOffX = 452.5, obsThreeTimesCutOffX = 562.5, obsOneThirdCutOffX = 362.5;
            obsPlusSysErrCutOffX = 452.5, obsMinusSysErrCutOffX = 422.5;
         }
         if( topoStr == "T1bbbb" ){
            obsOneTimesCutOffX = 472.5, obsThreeTimesCutOffX = 562.5, obsOneThirdCutOffX = 362.5;
            obsPlusSysErrCutOffX = 522.5, obsMinusSysErrCutOffX = 472.5;
         }
*/
         TH2F * obsExclOneTimesxSecProspinoTest = new TH2F("obsExclOneTimesxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[0]->Area(obsExclOneTimesxSecProspinoTest, Mzero, Mhalf, NLOObsxSecCL, xSecProspino, 1.0, doMyFill, xBinLargery);
         TGraph * obsExclOneTimesxSecProspinoTMP = plotTools[i]->GetContour(obsExclOneTimesxSecProspinoTest, 3, 1, !doMyFill);
         Smooth(obsExclOneTimesxSecProspinoTMP, radiusToSmooth);
         std::vector<double> xPts_obsExclOneTimesxSecProspino, yPts_obsExclOneTimesxSecProspino;
         for(int ip=0; ip<obsExclOneTimesxSecProspinoTMP->GetN(); ip++){
            double xgr=-1, ygr=-1;
            obsExclOneTimesxSecProspinoTMP->GetPoint(ip, xgr, ygr);
            if( xgr >= obsOneTimesCutOffX ){
               if( xgr != obsOneTimesCutOffX && xPts_obsExclOneTimesxSecProspino.empty() ){
                  xPts_obsExclOneTimesxSecProspino.push_back(obsOneTimesCutOffX); yPts_obsExclOneTimesxSecProspino.push_back(ygr);
               }
               xPts_obsExclOneTimesxSecProspino.push_back(xgr); yPts_obsExclOneTimesxSecProspino.push_back(ygr);
            }
         }
         int nPts_obsExclOneTimesxSecProspino = (int)xPts_obsExclOneTimesxSecProspino.size();
         TGraph *obsExclOneTimesxSecProspino = new TGraph(nPts_obsExclOneTimesxSecProspino, &xPts_obsExclOneTimesxSecProspino[0], &yPts_obsExclOneTimesxSecProspino[0]);
         if( topoStr == "T2" || topoStr == "T2qq" ){ obsExclOneTimesxSecProspino->RemovePoint(0); obsExclOneTimesxSecProspino->RemovePoint(0); }
         if( topoStr == "T5ZZInc" ){ obsExclOneTimesxSecProspino->RemovePoint(0); }
         if( topoStr == "T1bbbb" ){ obsExclOneTimesxSecProspino->RemovePoint(0); obsExclOneTimesxSecProspino->RemovePoint(0); }

         obsExclOneTimesxSecProspinoTest->Draw("colz");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"_TestContour.pdf");

         TH2F * obsExclThreeTimesxSecProspinoTest = new TH2F("obsExclThreeTimesxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[0]->Area(obsExclThreeTimesxSecProspinoTest, Mzero, Mhalf, NLOObsxSecCL, xSecProspino, 3.0, doMyFill, xBinLargery);
         TGraph * obsExclThreeTimesxSecProspinoTMP = plotTools[i]->GetContour(obsExclThreeTimesxSecProspinoTest, 3, 1, !doMyFill);
         Smooth(obsExclThreeTimesxSecProspinoTMP, radiusToSmooth);
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
         plotTools[0]->Area(obsExclOneThirdxSecProspinoTest, Mzero, Mhalf, NLOObsxSecCL, xSecProspino, 1./3.0, doMyFill, xBinLargery);
         TGraph * obsExclOneThirdxSecProspinoTMP = plotTools[i]->GetContour(obsExclOneThirdxSecProspinoTest, 3, 1, !doMyFill);
         Smooth(obsExclOneThirdxSecProspinoTMP, radiusToSmooth);
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
         TH2F * obsExclPlusSysErrxSecProspinoTest = new TH2F("obsExclPlusSysErrxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[0]->Area(obsExclPlusSysErrxSecProspinoTest, Mzero, Mhalf, NLOObsxSecCL, xSecProspino, 1.0, doMyFill, xBinLargery, 1);
         TGraph * obsExclPlusSysErrxSecProspinoTMP = plotTools[i]->GetContour(obsExclPlusSysErrxSecProspinoTest, 3, 1, !doMyFill);
         Smooth(obsExclPlusSysErrxSecProspinoTMP, radiusToSmooth);
         std::vector<double> xPts_obsExclPlusSysErrxSecProspino, yPts_obsExclPlusSysErrxSecProspino;
         for(int ip=0; ip<obsExclPlusSysErrxSecProspinoTMP->GetN(); ip++){
            double xgr=-1, ygr=-1;
            obsExclPlusSysErrxSecProspinoTMP->GetPoint(ip, xgr, ygr);
            if( xgr >= obsPlusSysErrCutOffX ){
               if( xgr != obsPlusSysErrCutOffX && xPts_obsExclPlusSysErrxSecProspino.empty() ){
                  xPts_obsExclPlusSysErrxSecProspino.push_back(obsPlusSysErrCutOffX); yPts_obsExclPlusSysErrxSecProspino.push_back(ygr);
               }
               xPts_obsExclPlusSysErrxSecProspino.push_back(xgr); yPts_obsExclPlusSysErrxSecProspino.push_back(ygr); 
            }
         }
         int nPts_obsExclPlusSysErrxSecProspino = (int)xPts_obsExclPlusSysErrxSecProspino.size();
         TGraph *obsExclPlusSysErrxSecProspino = new TGraph(nPts_obsExclPlusSysErrxSecProspino, &xPts_obsExclPlusSysErrxSecProspino[0], &yPts_obsExclPlusSysErrxSecProspino[0]);
         if( topoStr == "T2" || topoStr == "T2qq" ){ obsExclPlusSysErrxSecProspino->RemovePoint(0); }
         if( topoStr == "T5ZZInc" ){ obsExclPlusSysErrxSecProspino->RemovePoint(0); }
         if( topoStr == "T1bbbb" ){ obsExclPlusSysErrxSecProspino->RemovePoint(0); obsExclPlusSysErrxSecProspino->RemovePoint(0); }

         TH2F * obsExclMinusSysErrxSecProspinoTest = new TH2F("obsExclMinusSysErrxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[0]->Area(obsExclMinusSysErrxSecProspinoTest, Mzero, Mhalf, NLOObsxSecCL, xSecProspino, 1.0, doMyFill, xBinLargery, -1);
         TGraph * obsExclMinusSysErrxSecProspinoTMP = plotTools[i]->GetContour(obsExclMinusSysErrxSecProspinoTest, 3, 1, !doMyFill);
         Smooth(obsExclMinusSysErrxSecProspinoTMP, radiusToSmooth);
         std::vector<double> xPts_obsExclMinusSysErrxSecProspino, yPts_obsExclMinusSysErrxSecProspino;
         for(int ip=0; ip<obsExclMinusSysErrxSecProspinoTMP->GetN(); ip++){
            double xgr=-1, ygr=-1;
            obsExclMinusSysErrxSecProspinoTMP->GetPoint(ip, xgr, ygr);
            if( xgr >= obsMinusSysErrCutOffX ){
               if( xgr != obsMinusSysErrCutOffX && xPts_obsExclMinusSysErrxSecProspino.empty() ){
                  xPts_obsExclMinusSysErrxSecProspino.push_back(obsMinusSysErrCutOffX); yPts_obsExclMinusSysErrxSecProspino.push_back(ygr);
               }
               xPts_obsExclMinusSysErrxSecProspino.push_back(xgr); yPts_obsExclMinusSysErrxSecProspino.push_back(ygr); 
            }
         }
         int nPts_obsExclMinusSysErrxSecProspino = (int)xPts_obsExclMinusSysErrxSecProspino.size();
         TGraph *obsExclMinusSysErrxSecProspino = new TGraph(nPts_obsExclMinusSysErrxSecProspino, &xPts_obsExclMinusSysErrxSecProspino[0], &yPts_obsExclMinusSysErrxSecProspino[0]);

         if( topoStr == "T5ZZInc" ){ obsExclMinusSysErrxSecProspino->RemovePoint(0); obsExclMinusSysErrxSecProspino->RemovePoint(0); }
         if( topoStr == "T1bbbb" ){ obsExclMinusSysErrxSecProspino->RemovePoint(0); obsExclMinusSysErrxSecProspino->RemovePoint(0); }

         obsHistName = region[i]+"obsExclOneTimesProspino"; obsExclOneTimesProspinoVec[i] = (TGraph*) obsExclOneTimesxSecProspino->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclThreeTimesProspino"; obsExclThreeTimesProspinoVec[i] = (TGraph*) obsExclThreeTimesxSecProspino->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclOneThirdProspino"; obsExclOneThirdProspinoVec[i] = (TGraph*) obsExclOneThirdxSecProspino->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclPlusSysErrProspino"; obsExclPlusSysErrProspinoVec[i] = (TGraph*) obsExclPlusSysErrxSecProspino->Clone(obsHistName.Data());
         obsHistName = region[i]+"obsExclMinusSysErrProspino"; obsExclMinusSysErrProspinoVec[i] = (TGraph*) obsExclMinusSysErrxSecProspino->Clone(obsHistName.Data());

         minObsExclNLOcloned->SetMinimum(0.001);
         minObsExclNLOcloned->SetMaximum(20);
         minObsExclNLOcloned->Draw("colz");
         diagonalGraph->Draw("l");

         if( topoStr != "TGQ" ){
            obsExclOneTimesxSecProspino->SetLineWidth(2); obsExclOneTimesxSecProspino->Draw("l");
            obsExclThreeTimesxSecProspino->SetLineWidth(2); obsExclThreeTimesxSecProspino->SetLineStyle(2); obsExclThreeTimesxSecProspino->Draw("l");
            obsExclOneThirdxSecProspino->SetLineWidth(2); obsExclOneThirdxSecProspino->SetLineStyle(3); obsExclOneThirdxSecProspino->Draw("l");
            obsExclPlusSysErrxSecProspino->SetLineWidth(1); obsExclPlusSysErrxSecProspino->SetLineColor(kMagenta); obsExclPlusSysErrxSecProspino->SetLineStyle(1); obsExclPlusSysErrxSecProspino->Draw("l");
            obsExclMinusSysErrxSecProspino->SetLineWidth(1); obsExclMinusSysErrxSecProspino->SetLineColor(kMagenta); obsExclMinusSysErrxSecProspino->SetLineStyle(1); obsExclMinusSysErrxSecProspino->Draw("l");
         }
/*
         cmsPreTex.DrawLatex(0.20, 0.84, "CMS Preliminary");
         cmsLumiTex.DrawLatex(0.20, 0.77, "19.47 fb^{-1}, #sqrt{s} = 8 TeV");
         disptTex.DrawLatex(0.20, 0.70, "Jets+#slash{H}_{T}");
*/
         indexTex.DrawLatex(0.25, 0.55, "(b)");

         legexp->Clear();
         legexp->AddEntry(obsExclOneTimesxSecProspino, "#sigma^{prod}=#sigma^{NLO-QCD}", "l");
         legexp->AddEntry(obsExclThreeTimesxSecProspino, "#sigma^{prod}=3#times#sigma^{NLO-QCD}", "l");
         legexp->AddEntry(obsExclOneThirdxSecProspino, "#sigma^{prod}=1/3#times#sigma^{NLO-QCD}", "l");
         if( topoStr != "TGQ" ) legexp->Draw();

         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ObsLimit_mMother_mLSP.pdf");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ObsLimit_mMother_mLSP.png");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ObsLimit_mMother_mLSP.C");
   
   // Expected Limit in M0 - M1/2
         TString expHistName = region[i]+"expLimit";
         TH2F * minExpExclNLOcloned = (TH2F*) minExpExclNLO->Clone(); minExpExclNLOcloned->Reset();
         for(int ix=0; ix<=minExpExclNLO->GetXaxis()->GetNbins(); ix++){
            if( ix < xBinLargery ) continue;
            for(int iy=0; iy<=minExpExclNLO->GetYaxis()->GetNbins(); iy++){
               if( ix-iy < xBinLargery ) continue;
               minExpExclNLOcloned->SetBinContent(ix, iy, minExpExclNLO->GetBinContent(ix, iy));
            }
         }
//         minExpExclNLOcloned->GetXaxis()->SetRangeUser(300, 1200);
         explimitVec[i] = (TH2F*) minExpExclNLOcloned->Clone(expHistName.Data());

         double expOneTimesCutOffX = -1, expThreeTimesCutOffX = -1, expOneThirdCutOffX = -1;
         double expPlusOneSigmaCutOffX = -1, expMinusOneSigmaCutOffX = -1;
/*
         if( topoStr == "T1" || topoStr == "T1qqqq" ){
//            expOneTimesCutOffX = 537.5, expThreeTimesCutOffX = 662.5, expOneThirdCutOffX = 437.5;
            expOneTimesCutOffX = 506.5, expThreeTimesCutOffX = 660.5, expOneThirdCutOffX = 435.5;
//            expPlusOneSigmaCutOffX = 512.5, expMinusOneSigmaCutOffX = 562.5;
            expPlusOneSigmaCutOffX = 478.5, expMinusOneSigmaCutOffX = 560.5;
         }
         if( topoStr == "T2" || topoStr == "T2qq" ){
//            expOneTimesCutOffX = 437.5, expThreeTimesCutOffX = 537.5, expOneThirdCutOffX = 337.5;
            expOneTimesCutOffX = 400.5, expThreeTimesCutOffX = 531.5, expOneThirdCutOffX = 335.5;
//            expPlusOneSigmaCutOffX = 387.5, expMinusOneSigmaCutOffX = 487.5;
            expPlusOneSigmaCutOffX = 385.5, expMinusOneSigmaCutOffX = 450.5;
         }
         if( topoStr == "T5ZZInc" ){
            expOneTimesCutOffX = 478.5, expThreeTimesCutOffX = 587.5, expOneThirdCutOffX = 362.5;
            expPlusOneSigmaCutOffX = 445.5, expMinusOneSigmaCutOffX = 500.5;
         }
         if( topoStr == "T1bbbb" ){
            expOneTimesCutOffX = 472.5, expThreeTimesCutOffX = 587.5, expOneThirdCutOffX = 362.5;
            expPlusOneSigmaCutOffX = 470.5, expMinusOneSigmaCutOffX = 520.5;
         }
*/

         TH2F * expExclOneTimesxSecProspinoTest = new TH2F("expExclOneTimesxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[0]->Area(expExclOneTimesxSecProspinoTest, Mzero, Mhalf, NLOExpxSecCL, xSecProspino, 1.0, doMyFill, xBinLargery);
         TGraph * expExclOneTimesxSecProspinoTMP = plotTools[i]->GetContour(expExclOneTimesxSecProspinoTest, 3, 1, !doMyFill);
         Smooth(expExclOneTimesxSecProspinoTMP, radiusToSmooth);
         std::vector<double> xPts_expExclOneTimesxSecProspino, yPts_expExclOneTimesxSecProspino;
         for(int ip=0; ip<expExclOneTimesxSecProspinoTMP->GetN(); ip++){
            double xgr=-1, ygr=-1;
            expExclOneTimesxSecProspinoTMP->GetPoint(ip, xgr, ygr);
            if( xgr >= expOneTimesCutOffX ){
               if( xgr != expOneTimesCutOffX && xPts_expExclOneTimesxSecProspino.empty() ){
                  xPts_expExclOneTimesxSecProspino.push_back(expOneTimesCutOffX); yPts_expExclOneTimesxSecProspino.push_back(ygr);
               }
               xPts_expExclOneTimesxSecProspino.push_back(xgr); yPts_expExclOneTimesxSecProspino.push_back(ygr); 
            }
         }
         int nPts_expExclOneTimesxSecProspino = (int)xPts_expExclOneTimesxSecProspino.size();
         TGraph *expExclOneTimesxSecProspino = new TGraph(nPts_expExclOneTimesxSecProspino, &xPts_expExclOneTimesxSecProspino[0], &yPts_expExclOneTimesxSecProspino[0]);
         if( topoStr == "T1" || topoStr == "T1qqqq" ){ expExclOneTimesxSecProspino->RemovePoint(0); }
         if( topoStr == "T2" || topoStr == "T2qq" ){ expExclOneTimesxSecProspino->RemovePoint(0); expExclOneTimesxSecProspino->RemovePoint(0); }
         if( topoStr == "T1bbbb" ){ expExclOneTimesxSecProspino->RemovePoint(0); expExclOneTimesxSecProspino->RemovePoint(0); }

         expExclOneTimesxSecProspinoTest->Draw("colz");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"_TestContour.pdf");

// +1 Sigma
         TH2F * expExclPlusOneSigmaProspinoTest = new TH2F("expExclPlusOneSigmaProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[0]->Area(expExclPlusOneSigmaProspinoTest, Mzero, Mhalf, NLOExpxSecp1sigma, xSecProspino, 1.0, doMyFill, xBinLargery);
         TGraph * expExclPlusOneSigmaProspinoTMP = plotTools[i]->GetContour(expExclPlusOneSigmaProspinoTest, 3, 1, !doMyFill);
         Smooth(expExclPlusOneSigmaProspinoTMP, radiusToSmooth);
         std::vector<double> xPts_expExclPlusOneSigmaProspino, yPts_expExclPlusOneSigmaProspino;
         for(int ip=0; ip<expExclPlusOneSigmaProspinoTMP->GetN(); ip++){
            double xgr=-1, ygr=-1;
            expExclPlusOneSigmaProspinoTMP->GetPoint(ip, xgr, ygr);
            if( xgr >= expPlusOneSigmaCutOffX ){
               if( xgr != expPlusOneSigmaCutOffX && xPts_expExclPlusOneSigmaProspino.empty() ){
                  xPts_expExclPlusOneSigmaProspino.push_back(expPlusOneSigmaCutOffX); yPts_expExclPlusOneSigmaProspino.push_back(ygr);
               }
               xPts_expExclPlusOneSigmaProspino.push_back(xgr); yPts_expExclPlusOneSigmaProspino.push_back(ygr); 
            }
         }
         int nPts_expExclPlusOneSigmaProspino = (int)xPts_expExclPlusOneSigmaProspino.size();
         TGraph *expExclPlusOneSigmaProspino = new TGraph(nPts_expExclPlusOneSigmaProspino, &xPts_expExclPlusOneSigmaProspino[0], &yPts_expExclPlusOneSigmaProspino[0]);
         if( topoStr == "T1" || topoStr == "T1qqqq" ){ expExclPlusOneSigmaProspino->RemovePoint(0); expExclPlusOneSigmaProspino->RemovePoint(0); }
         if( topoStr == "T1bbbb" ){ expExclPlusOneSigmaProspino->RemovePoint(0); expExclPlusOneSigmaProspino->RemovePoint(0); }

// -1 Sigma
         TH2F * expExclMinusOneSigmaProspinoTest = new TH2F("expExclMinusOneSigmaProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[0]->Area(expExclMinusOneSigmaProspinoTest, Mzero, Mhalf, NLOExpxSecm1sigma, xSecProspino, 1.0, doMyFill, xBinLargery);
         TGraph * expExclMinusOneSigmaProspinoTMP = plotTools[i]->GetContour(expExclMinusOneSigmaProspinoTest, 3, 1, !doMyFill);
         Smooth(expExclMinusOneSigmaProspinoTMP, radiusToSmooth);
         std::vector<double> xPts_expExclMinusOneSigmaProspino, yPts_expExclMinusOneSigmaProspino;
         for(int ip=0; ip<expExclMinusOneSigmaProspinoTMP->GetN(); ip++){
            double xgr=-1, ygr=-1;
            expExclMinusOneSigmaProspinoTMP->GetPoint(ip, xgr, ygr);
            if( xgr >= expMinusOneSigmaCutOffX ){
               if( xgr != expMinusOneSigmaCutOffX && xPts_expExclMinusOneSigmaProspino.empty() ){
                  xPts_expExclMinusOneSigmaProspino.push_back(expMinusOneSigmaCutOffX); yPts_expExclMinusOneSigmaProspino.push_back(ygr);
               }
               xPts_expExclMinusOneSigmaProspino.push_back(xgr); yPts_expExclMinusOneSigmaProspino.push_back(ygr); 
            }
         }
         int nPts_expExclMinusOneSigmaProspino = (int)xPts_expExclMinusOneSigmaProspino.size();
         TGraph *expExclMinusOneSigmaProspino = new TGraph(nPts_expExclMinusOneSigmaProspino, &xPts_expExclMinusOneSigmaProspino[0], &yPts_expExclMinusOneSigmaProspino[0]);
         if( topoStr == "T2" || topoStr == "T2qq" ){ expExclMinusOneSigmaProspino->RemovePoint(0); expExclMinusOneSigmaProspino->RemovePoint(0); }
         if( topoStr == "T5ZZInc" ){ expExclMinusOneSigmaProspino->RemovePoint(0); }
         if( topoStr == "T1bbbb" ){ expExclMinusOneSigmaProspino->RemovePoint(0); expExclMinusOneSigmaProspino->RemovePoint(0); }


         TH2F * expExclThreeTimesxSecProspinoTest = new TH2F("expExclThreeTimesxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[0]->Area(expExclThreeTimesxSecProspinoTest, Mzero, Mhalf, NLOExpxSecCL, xSecProspino, 3.0, doMyFill, xBinLargery);
         TGraph * expExclThreeTimesxSecProspinoTMP = plotTools[i]->GetContour(expExclThreeTimesxSecProspinoTest, 3, 1, !doMyFill);
         Smooth(expExclThreeTimesxSecProspinoTMP, radiusToSmooth);
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
         if( topoStr == "T2" || topoStr == "T2qq" ){ expExclThreeTimesxSecProspino->RemovePoint(0); }

         TH2F * expExclOneThirdxSecProspinoTest = new TH2F("expExclOneThirdxSecProspinoTest", "", nXbins, loX, hiX, nYbins, loY, hiY);
         plotTools[0]->Area(expExclOneThirdxSecProspinoTest, Mzero, Mhalf, NLOExpxSecCL, xSecProspino, 1./3.0, doMyFill, xBinLargery);
         TGraph * expExclOneThirdxSecProspinoTMP = plotTools[i]->GetContour(expExclOneThirdxSecProspinoTest, 3, 1, !doMyFill);
         Smooth(expExclOneThirdxSecProspinoTMP, radiusToSmooth);
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


         expHistName = region[i]+"expExclOneTimesProspino"; expExclOneTimesProspinoVec[i] = (TGraph*) expExclOneTimesxSecProspino->Clone(expHistName.Data());
         expHistName = region[i]+"expExclThreeTimesProspino"; expExclThreeTimesProspinoVec[i] = (TGraph*) expExclThreeTimesxSecProspino->Clone(expHistName.Data());
         expHistName = region[i]+"expExclOneThirdProspino"; expExclOneThirdProspinoVec[i] = (TGraph*) expExclOneThirdxSecProspino->Clone(expHistName.Data());

         expHistName = region[i]+"expExclPlusOneSigmaProspino"; expExclPlusOneSigmaProspinoVec[i] = (TGraph*) expExclPlusOneSigmaProspino->Clone(expHistName.Data());
         expHistName = region[i]+"expExclMinusOneSigmaProspino"; expExclMinusOneSigmaProspinoVec[i] = (TGraph*) expExclMinusOneSigmaProspino->Clone(expHistName.Data());

         minExpExclNLOcloned->SetMinimum(0.001);
         minExpExclNLOcloned->SetMaximum(20);
         minExpExclNLOcloned->Draw("colz");
         diagonalGraph->Draw("l");

         if( topoStr != "TGQ" ){
            expExclOneTimesxSecProspino->SetLineWidth(2); expExclOneTimesxSecProspino->Draw("l");
            expExclThreeTimesxSecProspino->SetLineWidth(2); expExclThreeTimesxSecProspino->SetLineStyle(2); expExclThreeTimesxSecProspino->Draw("l");
            expExclOneThirdxSecProspino->SetLineWidth(2); expExclOneThirdxSecProspino->SetLineStyle(3); expExclOneThirdxSecProspino->Draw("l");

            expExclPlusOneSigmaProspino->SetLineWidth(2); expExclPlusOneSigmaProspino->SetLineStyle(8); expExclPlusOneSigmaProspino->Draw("l");
            expExclMinusOneSigmaProspino->SetLineWidth(2); expExclMinusOneSigmaProspino->SetLineStyle(8); expExclMinusOneSigmaProspino->Draw("l");
         }
/*
         cmsPreTex.DrawLatex(0.17, 0.84, "CMS Preliminary");
         cmsLumiTex.DrawLatex(0.17, 0.79, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
         disptTex.DrawLatex(0.17, 0.74, "Jets+#slash{H}_{T}");
*/
         legexp->Clear();
         legexp->AddEntry(expExclOneTimesxSecProspino, "#sigma^{prod}=#sigma^{NLO-QCD}", "l");
         legexp->AddEntry(expExclThreeTimesxSecProspino, "#sigma^{prod}=3#times#sigma^{NLO-QCD}", "l");
         legexp->AddEntry(expExclOneThirdxSecProspino, "#sigma^{prod}=1/3#times#sigma^{NLO-QCD}", "l");
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
      hexcl->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow q#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
   }else if( topoStr == "TGQ" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{g}} [GeV]; m_{#tilde{q}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{q}#tilde{g}, #tilde{q}#rightarrow q#tilde{#chi}^{0}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}");
   }else if( topoStr == "T5ZZInc" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}_{2}, #tilde{#chi}^{0}_{2}#rightarrow Z#tilde{#chi}^{0}");
   }else if( topoStr == "T1tttt" ){
      hexcl = new TH2F("hexcl",";m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g})");
   }else if( topoStr == "T1bbbb" ){
      hexcl = new TH2F("hexcl",";m(#tilde{g}) [GeV]; m(#tilde{#chi}^{0}) [GeV]; 95% CL Upper Limit on #sigma [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow bb#tilde{#chi}^{0}#; m(#tilde{b})>>m(#tilde{g})");
   }

/*
   hexcl->SetTitle("CMS Preliminary");
   hexcl->GetYaxis()->SetTitleOffset(1.3);
   hexcl->GetXaxis()->SetTitleOffset(0.92);
   hexcl->GetYaxis()->SetLabelSize(0.03);
   hexcl->GetXaxis()->SetLabelSize(0.03);
*/
   hexcl->Draw("colz");

   TLegend* leg = new TLegend(0.20,0.55,0.48,0.72,NULL,"brNDC"); 
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
   cmsLumiTex.DrawLatex(0.22, 0.78, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
   disptTex.DrawLatex(0.22, 0.74, "Jets+#slash{H}_{T} Observed Exclusion");
/*
   cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
   cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
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
      hexcl->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow q#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0});m_{#tilde{q}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "TGQ" ){
      hexcl->SetTitle("pp#rightarrow#tilde{q}#tilde{g}, #tilde{q}#rightarrow q#tilde{#chi}^{0}, #tilde{g}#rightarrow qq#tilde{#chi}^{0};m_{#tilde{g}} [GeV]; m_{#tilde{q}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T5ZZInc" ){
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}_{2}, #tilde{#chi}^{0}_{2}#rightarrow Z#tilde{#chi}^{0}; m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T1tttt" ){
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow tt#tilde{#chi}^{0}#; m(#tilde{t})>>m(#tilde{g}); m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]");
   }else if( topoStr == "T1bbbb" ){
      hexcl->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow bb#tilde{#chi}^{0}#; m(#tilde{b})>>m(#tilde{g})");
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
   cmsPreTex.DrawLatex(0.22, 0.82, "CMS Preliminary");
   cmsLumiTex.DrawLatex(0.22, 0.78, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
   disptTex.DrawLatex(0.22, 0.74, "Jets+#slash{H}_{T} Expected Exclusion");
/*
   cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
   cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
   disptTex.DrawLatex(350, 960, "Jets+#slash{H}_{T} Expected Exclusion");
*/
   leg->Draw();
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "comp_Exp_Exclusion_mMother_mLSP.pdf");
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "comp_Exp_Exclusion_mMother_mLSP.png");
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "comp_Exp_Exclusion_mMother_mLSP.C");

   c1->SetLogz(0);
   c1->SetRightMargin ( 0.1 );
   TLegendEntry *entry =0;

   minCatsObsExclNLO->SetMinimum(0);
   minCatsObsExclNLO->SetMaximum(nTotPlots);
   minCatsObsExclNLO->Draw("col");
   leg->Clear();
// FIXME: have to manually set fill color!?
   entry = leg->AddEntry(region[0], dispt[0], "f"); entry->SetFillColor(7); entry->SetFillStyle(1001);
   entry = leg->AddEntry(region[1], dispt[1], "f"); entry->SetFillColor(3); entry->SetFillStyle(1001);
   entry = leg->AddEntry(region[2], dispt[2], "f"); entry->SetFillColor(5); entry->SetFillStyle(1001);
   cmsPreTex.DrawLatex(0.22, 0.82, "CMS Preliminary");
   cmsLumiTex.DrawLatex(0.22, 0.78, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
   disptTex.DrawLatex(0.22, 0.74, "Jets+#slash{H}_{T}");
/*
   cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
   cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
//   disptTex.DrawLatex(350, 960, "Jets+#slash{H}_{T} 95% CL Observed Limit [pb]");
   disptTex.DrawLatex(350, 960, "Jets+#slash{H}_{T}");
*/
   leg->Draw();
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "min_Obs_Exclusion_mMother_mLSP.pdf");
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "min_Obs_Exclusion_mMother_mLSP.png");
   c1->SaveAs(outDir + "/"  +topoStr +"_" + "min_Obs_Exclusion_mMother_mLSP.C");

   minCatsExpExclNLO->SetMinimum(0);
   minCatsExpExclNLO->SetMaximum(nTotPlots);
   minCatsExpExclNLO->Draw("col");
   leg->Clear();
// FIXME: have to manually set fill color!?
   entry = leg->AddEntry(region[0], dispt[0], "f"); entry->SetFillColor(7); entry->SetFillStyle(1001);
   entry = leg->AddEntry(region[1], dispt[1], "f"); entry->SetFillColor(3); entry->SetFillStyle(1001);
   entry = leg->AddEntry(region[2], dispt[2], "f"); entry->SetFillColor(5); entry->SetFillStyle(1001);
   cmsPreTex.DrawLatex(0.22, 0.82, "CMS Preliminary");
   cmsLumiTex.DrawLatex(0.22, 0.78, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
   disptTex.DrawLatex(0.22, 0.74, "Jets+#slash{H}_{T}");
/*
   cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
   cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
//   disptTex.DrawLatex(350, 960, "Jets+#slash{H}_{T} 95% CL Expected Limit [pb]");
   disptTex.DrawLatex(350, 960, "Jets+#slash{H}_{T}");
*/
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
         cmsLumiTex.DrawLatex(0.22, 0.78, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
         disptTex.DrawLatex(0.22, 0.74, "Observed Ratio");
         disptTex4ratio.DrawLatex(0.22, 0.71, "("+dispt[jp]+")"+"/"+"("+dispt[ip]+")");
/*
         cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
         cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
         disptTex.DrawLatex(350, 960, "Observed Ratio");
         disptTex4ratio.DrawLatex(350, 890, "("+dispt[jp]+")"+"/"+"("+dispt[ip]+")");
*/
         c2->Print(outDir + "/" + topoStr + "_" + "ratio_limits_mMother_mLSP.pdf");

         expRatio2D->SetMinimum(0.02); expRatio2D->SetMaximum(5.0);
         expRatio2D->SetMarkerSize(0.5);
         expRatio2D->Draw("Text colz");
         cmsPreTex.DrawLatex(0.22, 0.82, "CMS Preliminary");
         cmsLumiTex.DrawLatex(0.22, 0.78, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
         disptTex.DrawLatex(0.22, 0.74, "Expected Ratio");
         disptTex4ratio.DrawLatex(0.22, 0.71, "("+dispt[jp]+")"+"/"+"("+dispt[ip]+")");
/*
         cmsPreTex.DrawLatex(350, 1100, "CMS Preliminary");
         cmsLumiTex.DrawLatex(350, 1030, "L_{int} = 19.47 fb^{-1}, #sqrt{s} = 8 TeV");
         disptTex.DrawLatex(350, 960, "Expected Ratio");
         disptTex4ratio.DrawLatex(350, 890, "("+dispt[jp]+")"+"/"+"("+dispt[ip]+")");
*/
         c2->Print(outDir + "/" + topoStr + "_" + "ratio_limits_mMother_mLSP.pdf");
      }
   }
   c2->Print(outDir + "/" + topoStr + "_" + "ratio_limits_mMother_mLSP.pdf]");

   TFile *outRootFile = new TFile(outDir+"/"+"CLs_SMS.root", "RECREATE");
   TDirectory *outdir = outRootFile->mkdir(topoStr, topoStr);
   outdir->cd();
   diagonalGraph->Write();
   for(int i=0; i<nTotPlots; i++){
      if( i!= nTotPlots -1 ) continue;
      obslimitVec[i]->Write(); explimitVec[i]->Write();
      obsExclOneTimesProspinoVec[i]->Write(); obsExclThreeTimesProspinoVec[i]->Write(); obsExclOneThirdProspinoVec[i]->Write();
      expExclOneTimesProspinoVec[i]->Write(); expExclThreeTimesProspinoVec[i]->Write(); expExclOneThirdProspinoVec[i]->Write();
      expExclPlusOneSigmaProspinoVec[i]->Write(); expExclMinusOneSigmaProspinoVec[i]->Write();
      obsExclPlusSysErrProspinoVec[i]->Write(); obsExclMinusSysErrProspinoVec[i]->Write();
   }
   outRootFile->Write(); outRootFile->Close();

}

int main(int argc, char** argv)
{
  return plot(argc, argv);
}
