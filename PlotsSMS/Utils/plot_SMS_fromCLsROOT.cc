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
static const std::string  disptStrs[] = {"H_{T}>500, #slash{H}_{T}>350 GeV", "H_{T}>800, #slash{H}_{T}>200 GeV", "H_{T}>800, #slash{H}_{T}>500 GeV", "H_{T}>1000, #slash{H}_{T}>400 GeV", "H_{T}>1200, #slash{H}_{T}>200 GeV"};
static const std::string regionDirs[] = {"mediumHTMHT", "highHT", "highHTMHT", "HT1000MHT400", "HT1200MHT200"};
static const int     inclFlags[] = {     1,           0,          0,           0,               0       };  

static const int      colors[] = {kGreen, kMagenta+1, kBlue, kRed, kTeal+4};

TH1D *xSecProspino_T1, *xSecProspino_T2, *xSecProspino;
TString topoStr = "T1";

// 25 GeV per X bin; 25 GeV per Y bin
const int nXbins = 45, nYbins = 47;
const double loX = 100, hiX = 1225;
const double loY = 50, hiY = 1225;
//const double scaleLegendFont = 1.8;
const double scaleLegendFont = 2.3;
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
   if( topoStr == "T1" || topoStr == "T5ZZInc" ) xSecProspino_T1 = (TH1D*) xSecProspinoFile->Get("gluino_xsection");
   if( topoStr == "T2" ) xSecProspino_T2 = (TH1D*) xSecProspinoFile->Get("squark_xsection");
   if( topoStr == "T2tt" ){
      TFile *xSecProspinoTGQFile = new TFile("reference_xSec_stop.root");
      xSecProspino_T2 = (TH1D*) xSecProspinoTGQFile->Get("stop_xsection");
   }
   if( topoStr == "T1" || topoStr == "T5ZZInc" ) xSecProspino = (TH1D*)xSecProspino_T1->Clone();
   if( topoStr == "T2" ) xSecProspino = (TH1D*)xSecProspino_T2->Clone();
   if( topoStr == "T2tt" ) xSecProspino = (TH1D*)xSecProspino_T2->Clone();
   if( topoStr == "TGQ" ) xSecProspino = (TH1D*)xSecProspino_T2->Clone();

   util::StyleSettings::paperNoTitle();
   gStyle->SetTitleSize(0.055,"XYZ");
   gStyle->SetTitleXOffset(1.02); gStyle->SetTitleYOffset(1.55); gStyle->SetTitleOffset(1.28, "Z");
   gStyle->SetPadLeftMargin(0.170); gStyle->SetPadRightMargin(0.125);
   gStyle->SetPadTopMargin(0.05); gStyle->SetPadBottomMargin(0.125); 

   gStyle->SetLabelSize(0.050, "XYZ");
 
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
/*
   gStyle->SetPadGridX(true);
   gStyle->SetPadGridY(true);
   gStyle->SetGridColor(0);
   gStyle->SetGridStyle(3);
   gStyle->SetGridWidth(1);
*/
   TCanvas * c1 = new TCanvas("c1","c1",800,800);
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
      yPtsDiag[ib] = loDiagRange + (ib+1)*widthBin; yErrsDiag[ib] = 0.0;
   }
   diagonalGraph = new TGraphErrors(nDiagBins, xPtsDiag, yPtsDiag, xErrsDiag, yErrsDiag);
   diagonalGraph->SetLineWidth(2); diagonalGraph->SetLineStyle(7);

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
   nTotPlots++; region.push_back("combined_"); regionBase.push_back("combined"); dispt.push_back("combined");
   std::cout<<"\nnTotPlots : "<<nTotPlots<<std::endl;

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

   TString outDir = "fromCLsROOT_plots";

   struct stat stFile;
   char commandline[200];
   if( stat(outDir.Data(), &stFile) == -1 ){ sprintf(commandline, "mkdir %s", outDir.Data()); system(commandline); }

   TH2* ObsExcl_NLO[nTotPlots];
   TH2* ExpExcl_NLO[nTotPlots];

   PlotTools<SusyScan> * plotTools[nTotPlots];

   std::vector<TH2F*> obslimitVec(nTotPlots), explimitVec(nTotPlots);
   std::vector<TH2F*> obslimitOriVec(nTotPlots), explimitOriVec(nTotPlots);
   std::vector<TGraph*> obsExclOneTimesProspinoVec(nTotPlots), expExclOneTimesProspinoVec(nTotPlots);
   std::vector<TGraph*> obsExclThreeTimesProspinoVec(nTotPlots), expExclThreeTimesProspinoVec(nTotPlots);
   std::vector<TGraph*> expExclPlusOneSigmaProspinoVec(nTotPlots), expExclMinusOneSigmaProspinoVec(nTotPlots);

   TH2F *minCatsObsExclNLO, *minCatsExpExclNLO;
   TH2F *minObsExclNLO, *minExpExclNLO, *minObsExclNLOfollowminExp;

   //Was not originally designed as a loop, so this will give you warnings about possible memory leaks.

//   TLegend* legexp = new TLegend(0.15,0.63,0.48,0.71,NULL,"brNDC"); 
   TLegend* legexp = new TLegend(0.18,0.70,0.48,0.88,NULL,"brNDC"); 
   legexp->SetFillColor(0);legexp->SetShadowColor(0);legexp->SetFillStyle(4000);legexp->SetTextFont(42);legexp->SetTextSize(0.023*scaleLegendFont);legexp->SetBorderSize(0);

   TLatex cmsPreTex; cmsPreTex.SetTextSize(0.025*scaleLegendFont); cmsPreTex.SetTextColor(1); cmsPreTex.SetNDC(); cmsPreTex.SetTextFont(42);
   TLatex cmsLumiTex; cmsLumiTex.SetTextSize(0.025*scaleLegendFont); cmsLumiTex.SetTextColor(1); cmsLumiTex.SetNDC(); cmsLumiTex.SetTextFont(42);
   TLatex disptTex; disptTex.SetTextSize(0.025*scaleLegendFont); disptTex.SetTextColor(1); disptTex.SetNDC(); disptTex.SetTextFont(42);

   TLatex indexTex; indexTex.SetTextSize(0.025*scaleLegendFont); indexTex.SetTextColor(1); indexTex.SetNDC(); indexTex.SetTextFont(42);

   TLatex cmsPreTexOnePad; cmsPreTexOnePad.SetTextSize(0.0200*scaleLegendFont); cmsPreTexOnePad.SetTextColor(1); cmsPreTexOnePad.SetNDC(); cmsPreTexOnePad.SetTextFont(42);
   TLatex cmsLumiTexOnePad; cmsLumiTexOnePad.SetTextSize(0.0200*scaleLegendFont); cmsLumiTexOnePad.SetTextColor(1); cmsLumiTexOnePad.SetNDC(); cmsLumiTexOnePad.SetTextFont(42);
   TLatex disptTexOnePad; disptTexOnePad.SetTextSize(0.0200*scaleLegendFont); disptTexOnePad.SetTextColor(1); disptTexOnePad.SetNDC(); disptTexOnePad.SetTextFont(42);

   std::vector<std::string> topoVec;
   std::vector<std::vector<TGraph*> > combProspinoVec, directCombProspinoVec;
   std::vector<std::vector<TGraphErrors*> > combDiagonalVec, directCombDiagonalVec;
   std::vector<std::vector<TH2F *> > combLimitsVec, directCombLimitsVec;
   std::vector<std::vector<TGraph*> > combExpBandsVec, directCombExpBandsVec;
   std::vector<std::vector<TGraph*> > combSysErrBandsVec, directCombSysErrBandsVec;

   char rtContourNames[200], rtLimitNames[200];

   std::map<std::string, int> mapIndexT1T2;

   TObjArray *vlist = topoStr.Tokenize(" ");
   int vlistEntries = vlist->GetEntries();

   int validIdxCnt =-1;
   for(int is=0; is<vlistEntries; is++){

      TObjString* perObj = dynamic_cast<TObjString*>(vlist->At(is)); 
      TString perStr = perObj->GetString();

      sprintf(rtContourNames, "CLs_SMS_contour_%s.root", perStr.Data());
      sprintf(rtLimitNames, "CLs_SMS_limit_%s.root", perStr.Data());

      if( stat(rtContourNames, &stFile) == -1 || stat(rtLimitNames, &stFile) == -1 ){
         std::cout<<"Not all root files related to topo : "<<perStr<<"  can be found!?"<<std::endl;
         continue;
      } 

      validIdxCnt ++;

      if( perStr == "T1" ) mapIndexT1T2["T1"] = validIdxCnt;
      if( perStr == "T2" ) mapIndexT1T2["T2"] = validIdxCnt;

      topoVec.push_back(perStr.Data());
      std::vector<TGraph*> tmpCombProspinoVec, tmpDirectCombProspinoVec;
      std::vector<TH2F*> tmpCombLimitsVec, tmpDirectCombLimitsVec;
      std::vector<TGraphErrors*> tmpCombDiagonalVec, tmpDirectCombDiagonalVec;
      std::vector<TGraph*> tmpCombExpBandsVec, tmpDirectCombExpBandsVec;
      std::vector<TGraph*> tmpCombSysErrBandsVec, tmpDirectCombSysErrBandsVec;

      TFile *contourFile = new TFile(rtContourNames);
      tmpCombProspinoVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_obsExclOneThirdProspino")->Clone() );
      tmpCombProspinoVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_obsExclOneTimesProspino")->Clone() );
      tmpCombProspinoVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_obsExclThreeTimesProspino")->Clone() );
      tmpCombProspinoVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_expExclOneThirdProspino")->Clone() );
      tmpCombProspinoVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_expExclOneTimesProspino")->Clone() );
      tmpCombProspinoVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_expExclThreeTimesProspino")->Clone() );
      tmpCombDiagonalVec.push_back( (TGraphErrors*) contourFile->Get(perStr + "/" + "diagonalGraph")->Clone() );
      tmpCombLimitsVec.push_back( (TH2F*) contourFile->Get(perStr + "/" + "combined_obsLimit")->Clone() );
      tmpCombLimitsVec.push_back( (TH2F*) contourFile->Get(perStr + "/" + "combined_expLimit")->Clone() );
      tmpCombExpBandsVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_expExclPlusOneSigmaProspino")->Clone() );
      tmpCombExpBandsVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_expExclMinusOneSigmaProspino")->Clone() );
      tmpCombSysErrBandsVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_obsExclPlusSysErrProspino")->Clone() );
      tmpCombSysErrBandsVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_obsExclMinusSysErrProspino")->Clone() );

      TFile *limitFile = new TFile(rtLimitNames);
      tmpDirectCombProspinoVec.push_back( (TGraph*) limitFile->Get(perStr + "/" + "combined_obsExclOneThirdProspino")->Clone() );
      tmpDirectCombProspinoVec.push_back( (TGraph*) limitFile->Get(perStr + "/" + "combined_obsExclOneTimesProspino")->Clone() );
      tmpDirectCombProspinoVec.push_back( (TGraph*) limitFile->Get(perStr + "/" + "combined_obsExclThreeTimesProspino")->Clone() );
      tmpDirectCombProspinoVec.push_back( (TGraph*) limitFile->Get(perStr + "/" + "combined_expExclOneThirdProspino")->Clone() );
      tmpDirectCombProspinoVec.push_back( (TGraph*) limitFile->Get(perStr + "/" + "combined_expExclOneTimesProspino")->Clone() );
      tmpDirectCombProspinoVec.push_back( (TGraph*) limitFile->Get(perStr + "/" + "combined_expExclThreeTimesProspino")->Clone() );
      tmpDirectCombDiagonalVec.push_back( (TGraphErrors*) limitFile->Get(perStr + "/" + "diagonalGraph")->Clone() );
      tmpDirectCombLimitsVec.push_back( (TH2F*) limitFile->Get(perStr + "/" + "combined_obsLimit")->Clone() );
      tmpDirectCombLimitsVec.push_back( (TH2F*) limitFile->Get(perStr + "/" + "combined_expLimit")->Clone() );
      tmpDirectCombExpBandsVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_expExclPlusOneSigmaProspino")->Clone() );
      tmpDirectCombExpBandsVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_expExclMinusOneSigmaProspino")->Clone() );
      tmpDirectCombSysErrBandsVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_obsExclPlusSysErrProspino")->Clone() );
      tmpDirectCombSysErrBandsVec.push_back( (TGraph*) contourFile->Get(perStr + "/" + "combined_obsExclMinusSysErrProspino")->Clone() );
      
      combProspinoVec.push_back(tmpCombProspinoVec); directCombProspinoVec.push_back(tmpDirectCombProspinoVec);
      combLimitsVec.push_back(tmpCombLimitsVec); directCombLimitsVec.push_back(tmpDirectCombLimitsVec);
      combDiagonalVec.push_back(tmpCombDiagonalVec); directCombDiagonalVec.push_back(tmpDirectCombDiagonalVec);
      combExpBandsVec.push_back(tmpCombExpBandsVec); directCombExpBandsVec.push_back(tmpDirectCombExpBandsVec); 
      combSysErrBandsVec.push_back(tmpCombSysErrBandsVec); directCombSysErrBandsVec.push_back(tmpDirectCombSysErrBandsVec); 
   }

   c1->SetLogz(1);
/*
//   c1->SetRightMargin ( 0.19 );
   c1->SetRightMargin ( 0.195 );
   c1->SetTopMargin ( 0.085 );
*/
   c1->SetRightMargin ( 0.192 );
   c1->SetTopMargin ( 0.095 );

   int validTopoSize = (int) topoVec.size();
   for(int it=0; it<validTopoSize; it++){

      std::string perTopo = topoVec[it]; 

// Observed Limit in M0 - M1/2
      directCombLimitsVec[it][0]->SetMinimum(0.001);
      directCombLimitsVec[it][0]->SetMaximum(200);
      directCombLimitsVec[it][0]->Draw("colz");
      directCombDiagonalVec[it][0]->Draw("l");
      
      directCombProspinoVec[it][0]->SetLineWidth(2); directCombProspinoVec[it][0]->SetLineStyle(3); directCombProspinoVec[it][0]->Draw("l"); 
      directCombProspinoVec[it][1]->SetLineWidth(2); directCombProspinoVec[it][1]->Draw("l"); 
      directCombProspinoVec[it][2]->SetLineWidth(2); directCombProspinoVec[it][2]->SetLineStyle(2); directCombProspinoVec[it][2]->Draw("l");
      directCombSysErrBandsVec[it][0]->SetLineWidth(1); directCombSysErrBandsVec[it][0]->SetLineColor(kMagenta); directCombSysErrBandsVec[it][0]->Draw("l");
      directCombSysErrBandsVec[it][1]->SetLineWidth(1); directCombSysErrBandsVec[it][1]->SetLineColor(kMagenta); directCombSysErrBandsVec[it][1]->Draw("l");

      legexp->Clear();
      legexp->AddEntry(directCombProspinoVec[it][1], "#sigma^{prod}=#sigma^{NLO-QCD}", "l");
      legexp->AddEntry(directCombProspinoVec[it][2], "#sigma^{prod}=3#times#sigma^{NLO-QCD}", "l");
      legexp->AddEntry(directCombProspinoVec[it][0], "#sigma^{prod}=1/3#times#sigma^{NLO-QCD}", "l");

      legexp->Draw();

      c1->SaveAs(outDir + "/" +perTopo +"_" + region[nTotPlots-1]+"ObsLimit_mMother_mLSP.pdf");
      c1->SaveAs(outDir + "/" +perTopo +"_" + region[nTotPlots-1]+"ObsLimit_mMother_mLSP.png");
      c1->SaveAs(outDir + "/" +perTopo +"_" + region[nTotPlots-1]+"ObsLimit_mMother_mLSP.C");

// Expected Limit in M0 - M1/2
      directCombLimitsVec[it][1]->SetMinimum(0.001);
      directCombLimitsVec[it][1]->SetMaximum(200);
      directCombLimitsVec[it][1]->Draw("colz");
      directCombDiagonalVec[it][0]->Draw("l");
      
      directCombProspinoVec[it][3]->SetLineWidth(2); directCombProspinoVec[it][3]->SetLineStyle(3); directCombProspinoVec[it][3]->Draw("l"); 
      directCombProspinoVec[it][4]->SetLineWidth(2); directCombProspinoVec[it][4]->Draw("l"); 
      directCombProspinoVec[it][5]->SetLineWidth(2); directCombProspinoVec[it][5]->SetLineStyle(2); directCombProspinoVec[it][5]->Draw("l"); 

      if( perTopo != "T5ZZInc" ){
         directCombExpBandsVec[it][0]->SetLineWidth(2); directCombExpBandsVec[it][0]->SetLineStyle(7); directCombExpBandsVec[it][0]->Draw("l");
         directCombExpBandsVec[it][1]->SetLineWidth(2); directCombExpBandsVec[it][1]->SetLineStyle(7); directCombExpBandsVec[it][1]->Draw("l");
      }

      legexp->Clear();
      legexp->AddEntry(directCombProspinoVec[it][4], "#sigma^{prod}=#sigma^{NLO-QCD}", "l");
      legexp->AddEntry(directCombProspinoVec[it][5], "#sigma^{prod}=3#times#sigma^{NLO-QCD}", "l");
      legexp->AddEntry(directCombProspinoVec[it][3], "#sigma^{prod}=1/3#times#sigma^{NLO-QCD}", "l");

      legexp->Draw();

      c1->SaveAs(outDir + "/" +perTopo +"_" + region[nTotPlots-1]+"ExpLimit_mMother_mLSP.pdf");
      c1->SaveAs(outDir + "/" +perTopo +"_" + region[nTotPlots-1]+"ExpLimit_mMother_mLSP.png");
      c1->SaveAs(outDir + "/" +perTopo +"_" + region[nTotPlots-1]+"ExpLimit_mMother_mLSP.C");

// Both obs and exp limits are in the same plots...
      directCombLimitsVec[it][0]->SetMinimum(0.001);
      directCombLimitsVec[it][0]->SetMaximum(10);
      directCombLimitsVec[it][0]->Draw("colz");
      directCombDiagonalVec[it][0]->Draw("l");

      directCombProspinoVec[it][1]->SetLineWidth(2); directCombProspinoVec[it][1]->Draw("l");
      directCombSysErrBandsVec[it][0]->SetLineWidth(1); directCombSysErrBandsVec[it][0]->SetLineStyle(1); directCombSysErrBandsVec[it][0]->SetLineColor(kBlack); directCombSysErrBandsVec[it][0]->Draw("l");
      directCombSysErrBandsVec[it][1]->SetLineWidth(1); directCombSysErrBandsVec[it][1]->SetLineStyle(1); directCombSysErrBandsVec[it][1]->SetLineColor(kBlack); directCombSysErrBandsVec[it][1]->Draw("l");
      
      cmsPreTexOnePad.DrawLatex(0.19, 0.84, "CMS, 4.98 fb^{-1}, #sqrt{s} = 7 TeV");
      
      directCombProspinoVec[it][4]->SetLineWidth(2); directCombProspinoVec[it][4]->SetLineColor(kMagenta+2); directCombProspinoVec[it][4]->SetLineStyle(2); directCombProspinoVec[it][4]->Draw("l");

      directCombExpBandsVec[it][0]->SetLineWidth(1); directCombExpBandsVec[it][0]->SetLineColor(kMagenta+2); directCombExpBandsVec[it][0]->SetLineStyle(3); directCombExpBandsVec[it][0]->Draw("l");
      directCombExpBandsVec[it][1]->SetLineWidth(1); directCombExpBandsVec[it][1]->SetLineColor(kMagenta+2); directCombExpBandsVec[it][1]->SetLineStyle(3); directCombExpBandsVec[it][1]->Draw("l");

      TLegend* legComb = new TLegend(0.19, 0.62, 0.39, 0.82, NULL, "brNDC");
      legComb->SetFillColor(0);legComb->SetShadowColor(0);legComb->SetFillStyle(4000);legComb->SetTextFont(42);legComb->SetTextSize(0.0200*scaleLegendFont);legComb->SetBorderSize(0);
      legComb->Clear();
      legComb->AddEntry(directCombProspinoVec[it][1], "#sigma^{NLO+NLL} #pm 1#sigma sig. theory", "l");
      legComb->AddEntry(directCombProspinoVec[it][4], "Exp. limit #pm 1#sigma exp.", "l");
      legComb->AddEntry(directCombProspinoVec[it][0], " ", "");
      legComb->AddEntry(directCombProspinoVec[it][2], " ", "");
      legComb->Draw();

      TLegend* legObsSysErr1 = new TLegend(0.19, 0.75+0.04, 0.39, 0.75+0.07, NULL,"brNDC");
      legObsSysErr1->SetFillColor(0);legObsSysErr1->SetShadowColor(0);legObsSysErr1->SetFillStyle(4000);legObsSysErr1->SetTextFont(42);legObsSysErr1->SetTextSize(0.0200*scaleLegendFont);legObsSysErr1->SetBorderSize(0);
      legObsSysErr1->AddEntry(directCombSysErrBandsVec[it][0], " ", "l");
      legObsSysErr1->Draw();

      TLegend* legObsSysErr2 = new TLegend(0.19, 0.80-0.03, 0.39, 0.80, NULL,"brNDC");
      legObsSysErr2->SetFillColor(0);legObsSysErr2->SetShadowColor(0);legObsSysErr2->SetFillStyle(4000);legObsSysErr2->SetTextFont(42);legObsSysErr2->SetTextSize(0.0200*scaleLegendFont);legObsSysErr2->SetBorderSize(0);
      legObsSysErr2->AddEntry(directCombSysErrBandsVec[it][1], " ", "l");
      legObsSysErr2->Draw();

      TLegend* legExpBand1 = new TLegend(0.19, 0.70+0.04, 0.39, 0.70+0.07, NULL,"brNDC");
      legExpBand1->SetFillColor(0);legExpBand1->SetShadowColor(0);legExpBand1->SetFillStyle(4000);legExpBand1->SetTextFont(42);legExpBand1->SetTextSize(0.0200*scaleLegendFont);legExpBand1->SetBorderSize(0);
      legExpBand1->AddEntry(directCombExpBandsVec[it][0], " ", "l");
      legExpBand1->Draw();

      TLegend* legExpBand2 = new TLegend(0.19, 0.75-0.03, 0.39, 0.75, NULL,"brNDC");
      legExpBand2->SetFillColor(0);legExpBand2->SetShadowColor(0);legExpBand2->SetFillStyle(4000);legExpBand2->SetTextFont(42);legExpBand2->SetTextSize(0.0200*scaleLegendFont);legExpBand2->SetBorderSize(0);
      legExpBand2->AddEntry(directCombExpBandsVec[it][1], " ", "l");
      legExpBand2->Draw();

      disptTexOnePad.DrawLatex(0.19, 0.67, "Jets+#slash{H}_{T}");

      c1->SaveAs(outDir + "/" +perTopo +"_" + region[nTotPlots-1]+"allINone.pdf");
      c1->SaveAs(outDir + "/" +perTopo +"_" + region[nTotPlots-1]+"allINone.png");
      c1->SaveAs(outDir + "/" +perTopo +"_" + region[nTotPlots-1]+"allINone.C");

   }
/*
   const double t_ = 0.02;
   const double b_ = 0.15;
   const double y_ = (1.-t_-b_);
   const double l_ = 0.07;
   const double r_ = 0.01;
   const double x_ = (1.-l_-r_)/2.;
   const double fracHeightTopPad = 0.8;
   const double ratioMin = 0.1;
   const double ratioMax = 2.4;
*/
   const double t_ = 0.08;
   const double b_ = 0.125;
   const double y_ = (1.-t_-b_);
//   const double l_ = 0.085;
   const double l_ = 0.088;
//   const double r_ = 0.10;
   const double r_ = 0.096;
   const double x_ = (1.-l_-r_)/2.;
   const double fracHeightTopPad = 0.8;
   const double ratioMin = 0.1;
   const double ratioMax = 2.4;

   if( mapIndexT1T2.size() == 2 ){

      int idxT1 = mapIndexT1T2.find("T1")->second;
      int idxT2 = mapIndexT1T2.find("T2")->second;

      TCanvas* canComb = new TCanvas("canComb","T1 + T2",800,800*x_/y_);

      TPad* T1Pad = new TPad("T1Pad","",0.,0.,l_+x_,1.);
      T1Pad->SetTopMargin(t_);
      T1Pad->SetLeftMargin(2.*l_);
      T1Pad->SetBottomMargin(b_);
      T1Pad->SetRightMargin(0.);
      T1Pad->SetFillStyle(0);
      T1Pad->SetFrameFillColor(10);
      T1Pad->SetFrameBorderMode(0);

// Observed limits
      canComb->cd(); T1Pad->Draw(); T1Pad->cd(); T1Pad->SetLogz(1);

      directCombLimitsVec[idxT1][0]->SetMinimum(0.001);
      directCombLimitsVec[idxT1][0]->SetMaximum(10);
      directCombLimitsVec[idxT1][0]->SetTitleOffset(1.58, "Y");
      directCombLimitsVec[idxT1][0]->Draw("col");
      directCombDiagonalVec[idxT1][0]->Draw("l");

      directCombProspinoVec[idxT1][0]->SetLineWidth(2); directCombProspinoVec[idxT1][0]->SetLineStyle(3); directCombProspinoVec[idxT1][0]->Draw("l");
      directCombProspinoVec[idxT1][1]->SetLineWidth(2); directCombProspinoVec[idxT1][1]->Draw("l");
      directCombProspinoVec[idxT1][2]->SetLineWidth(2); directCombProspinoVec[idxT1][2]->SetLineStyle(2); directCombProspinoVec[idxT1][2]->Draw("l");

      cmsPreTex.DrawLatex(0.20, 0.84, "CMS Preliminary");
//      cmsPreTex.DrawLatex(0.20, 0.84, "CMS");
      cmsLumiTex.DrawLatex(0.20, 0.77, "4.98 fb^{-1}, #sqrt{s} = 7 TeV");
      disptTex.DrawLatex(0.20, 0.70, "Jets+#slash{H}_{T}");

      indexTex.DrawLatex(0.25, 0.50, "(a)");

      TPad* T2Pad = new TPad("T2Pad","",l_+x_,0.,1.,1.);
      T2Pad->SetTopMargin(t_);
      T2Pad->SetLeftMargin(0.);
      T2Pad->SetBottomMargin(b_);
      T2Pad->SetRightMargin(2.*r_);
      T2Pad->SetFillStyle(0);
      T2Pad->SetFrameFillColor(10);
      T2Pad->SetFrameBorderMode(0);

      canComb->cd(); T2Pad->Draw(); T2Pad->cd(); T2Pad->SetLogz(1);

      directCombLimitsVec[idxT2][0]->SetMinimum(0.001);
      directCombLimitsVec[idxT2][0]->SetMaximum(10);
      directCombLimitsVec[idxT2][0]->Draw("colz");
      T2Pad->Update();
      TPaletteAxis *obsT2PaletteAxis = (TPaletteAxis*) directCombLimitsVec[idxT2][0]->GetListOfFunctions()->FindObject("palette");
      obsT2PaletteAxis->SetTitleOffset(1.32);
      obsT2PaletteAxis->SetLabelOffset(0.002);
      T2Pad->Draw();

      TPaveText *titleT2Pave = (TPaveText*) T2Pad->GetPrimitive("title")->Clone();

      TPaveText* coverT2title = new TPaveText(titleT2Pave->GetX1NDC(), titleT2Pave->GetY1NDC(), titleT2Pave->GetX2NDC(), titleT2Pave->GetY2NDC(), "NDC");
      coverT2title->SetBorderSize(0);
      coverT2title->SetFillColor(0);
      coverT2title->Draw();

      double changeXNDCT2 = titleT2Pave->GetX1NDC()-titleT2Pave->GetX1NDC()*(l_+x_)+l_;
      titleT2Pave->SetX1NDC(titleT2Pave->GetX1NDC()-changeXNDCT2);
      titleT2Pave->SetX2NDC(titleT2Pave->GetX2NDC()-changeXNDCT2);
      titleT2Pave->Draw(); 

      directCombDiagonalVec[idxT2][0]->Draw("l");

      directCombProspinoVec[idxT2][0]->SetLineWidth(2); directCombProspinoVec[idxT2][0]->SetLineStyle(3); directCombProspinoVec[idxT2][0]->Draw("l");
      directCombProspinoVec[idxT2][1]->SetLineWidth(2); directCombProspinoVec[idxT2][1]->Draw("l");
      directCombProspinoVec[idxT2][2]->SetLineWidth(2); directCombProspinoVec[idxT2][2]->SetLineStyle(2); directCombProspinoVec[idxT2][2]->Draw("l");

      double cacheX1legexp = legexp->GetX1NDC(), cacheX2legexp = legexp->GetX2NDC();
      double changeXNDClegexp = cacheX1legexp-cacheX1legexp*(l_+x_)+l_;
      TLegend* legCombT2 = new TLegend(cacheX1legexp-changeXNDClegexp, legexp->GetY1NDC(), cacheX2legexp, legexp->GetY2NDC(),NULL,"brNDC");
      legCombT2->SetFillColor(0);legCombT2->SetShadowColor(0);legCombT2->SetFillStyle(4000);legCombT2->SetTextFont(42);legCombT2->SetTextSize(0.025*scaleLegendFont);legCombT2->SetBorderSize(0);

      legCombT2->Clear();
      legCombT2->AddEntry(directCombProspinoVec[idxT2][0], "#sigma^{prod}=1/3#times#sigma^{NLO-QCD}", "l");
      legCombT2->AddEntry(directCombProspinoVec[idxT2][1], "#sigma^{prod}=#sigma^{NLO-QCD}", "l");
      legCombT2->AddEntry(directCombProspinoVec[idxT2][2], "#sigma^{prod}=3#times#sigma^{NLO-QCD}", "l");
      legCombT2->Draw();

      indexTex.DrawLatex(0.25*(l_+x_), 0.50, "(b)");

      TAxis *axisXT1 = (TAxis*)directCombLimitsVec[idxT1][0]->GetXaxis();
      TPaveText* axisXmodifT1 = new TPaveText(0.925, gStyle->GetPadBottomMargin()-0.05, 1.0, gStyle->GetPadBottomMargin()-0.01, "NDC");
      axisXmodifT1->SetTextSize(axisXT1->GetLabelSize()); axisXmodifT1->SetTextColor(axisXT1->GetLabelColor()); axisXmodifT1->SetTextFont(axisXT1->GetLabelFont());
      axisXmodifT1->SetBorderSize(0);
      axisXmodifT1->SetFillColor(0);

      canComb->cd(); T1Pad->cd(); axisXmodifT1->Draw();
  
      canComb->SaveAs(outDir + "/" + "T1T2" +"_" + region[nTotPlots-1]+"ObsLimit_mMother_mLSP.pdf");
      canComb->SaveAs(outDir + "/" + "T1T2" +"_" + region[nTotPlots-1]+"ObsLimit_mMother_mLSP.png");
      canComb->SaveAs(outDir + "/" + "T1T2" +"_" + region[nTotPlots-1]+"ObsLimit_mMother_mLSP.C");

// Expected limits
      canComb->cd(); T1Pad->Draw(); T1Pad->cd(); T1Pad->SetLogz(1);

      directCombLimitsVec[idxT1][1]->SetMinimum(0.001);
      directCombLimitsVec[idxT1][1]->SetMaximum(10);
      directCombLimitsVec[idxT1][1]->SetTitleOffset(1.58, "Y");
      directCombLimitsVec[idxT1][1]->Draw("col");
      directCombDiagonalVec[idxT1][0]->Draw("l");

      directCombProspinoVec[idxT1][3]->SetLineWidth(2); directCombProspinoVec[idxT1][3]->SetLineStyle(3); directCombProspinoVec[idxT1][3]->Draw("l");
      directCombProspinoVec[idxT1][4]->SetLineWidth(2); directCombProspinoVec[idxT1][4]->Draw("l");
      directCombProspinoVec[idxT1][5]->SetLineWidth(2); directCombProspinoVec[idxT1][5]->SetLineStyle(2); directCombProspinoVec[idxT1][5]->Draw("l");

      directCombExpBandsVec[idxT1][0]->SetLineWidth(2); directCombExpBandsVec[idxT1][0]->SetLineStyle(7); directCombExpBandsVec[idxT1][0]->Draw("l");
      directCombExpBandsVec[idxT1][1]->SetLineWidth(2); directCombExpBandsVec[idxT1][1]->SetLineStyle(7); directCombExpBandsVec[idxT1][1]->Draw("l");

      cmsPreTex.DrawLatex(0.20, 0.84, "CMS Preliminary");
      cmsPreTex.DrawLatex(0.20, 0.84, "CMS");
      cmsLumiTex.DrawLatex(0.20, 0.77, "4.98 fb^{-1}, #sqrt{s} = 7 TeV");
      disptTex.DrawLatex(0.20, 0.70, "Jets+#slash{H}_{T}");

      indexTex.DrawLatex(0.25, 0.50, "(a)");

      canComb->cd(); T2Pad->Draw(); T2Pad->cd(); T2Pad->SetLogz(1);

      directCombLimitsVec[idxT2][1]->SetMinimum(0.001);
      directCombLimitsVec[idxT2][1]->SetMaximum(10);
      directCombLimitsVec[idxT2][1]->Draw("colz");
      T2Pad->Update();
      TPaletteAxis *expT2PaletteAxis = (TPaletteAxis*) directCombLimitsVec[idxT2][1]->GetListOfFunctions()->FindObject("palette");
      expT2PaletteAxis->SetTitleOffset(1.32);
      expT2PaletteAxis->SetLabelOffset(0.002);
      T2Pad->Draw();

      titleT2Pave = (TPaveText*) T2Pad->GetPrimitive("title")->Clone();

      coverT2title = new TPaveText(titleT2Pave->GetX1NDC(), titleT2Pave->GetY1NDC(), titleT2Pave->GetX2NDC(), titleT2Pave->GetY2NDC(), "NDC");
      coverT2title->SetBorderSize(0);
      coverT2title->SetFillColor(0);
      coverT2title->Draw();

      changeXNDCT2 = titleT2Pave->GetX1NDC()-titleT2Pave->GetX1NDC()*(l_+x_)+l_;
      titleT2Pave->SetX1NDC(titleT2Pave->GetX1NDC()-changeXNDCT2);
      titleT2Pave->SetX2NDC(titleT2Pave->GetX2NDC()-changeXNDCT2);
      titleT2Pave->Draw(); 

      directCombDiagonalVec[idxT2][0]->Draw("l");

      directCombProspinoVec[idxT2][3]->SetLineWidth(2); directCombProspinoVec[idxT2][3]->SetLineStyle(3); directCombProspinoVec[idxT2][3]->Draw("l");
      directCombProspinoVec[idxT2][4]->SetLineWidth(2); directCombProspinoVec[idxT2][4]->Draw("l");
      directCombProspinoVec[idxT2][5]->SetLineWidth(2); directCombProspinoVec[idxT2][5]->SetLineStyle(2); directCombProspinoVec[idxT2][5]->Draw("l");

      directCombExpBandsVec[idxT2][0]->SetLineWidth(2); directCombExpBandsVec[idxT2][0]->SetLineStyle(7); directCombExpBandsVec[idxT2][0]->Draw("l");
      directCombExpBandsVec[idxT2][1]->SetLineWidth(2); directCombExpBandsVec[idxT2][1]->SetLineStyle(7); directCombExpBandsVec[idxT2][1]->Draw("l");

      cacheX1legexp = legexp->GetX1NDC(), cacheX2legexp = legexp->GetX2NDC();
      changeXNDClegexp = cacheX1legexp-cacheX1legexp*(l_+x_)+l_ - 0.01;
      legCombT2 = new TLegend(cacheX1legexp-changeXNDClegexp, legexp->GetY1NDC()-0.08, cacheX2legexp, legexp->GetY2NDC(),NULL,"brNDC");
      legCombT2->SetFillColor(0);legCombT2->SetShadowColor(0);legCombT2->SetFillStyle(4000);legCombT2->SetTextFont(42);legCombT2->SetTextSize(0.025*scaleLegendFont);legCombT2->SetBorderSize(0);

      legCombT2->Clear();
      legCombT2->AddEntry(directCombProspinoVec[idxT2][3], "#sigma^{prod}=1/3#times#sigma^{NLO-QCD}", "l");
      legCombT2->AddEntry(directCombProspinoVec[idxT2][4], "#sigma^{prod}=#sigma^{NLO-QCD}", "l");
      legCombT2->AddEntry(directCombProspinoVec[idxT2][5], "#sigma^{prod}=3#times#sigma^{NLO-QCD}", "l");
      legCombT2->Draw();

      indexTex.DrawLatex(0.25*(l_+x_), 0.50, "(b)");

      canComb->cd(); T1Pad->cd(); axisXmodifT1->Draw();
  
      canComb->SaveAs(outDir + "/" + "T1T2" +"_" + region[nTotPlots-1]+"ExpLimit_mMother_mLSP.pdf");
      canComb->SaveAs(outDir + "/" + "T1T2" +"_" + region[nTotPlots-1]+"ExpLimit_mMother_mLSP.png");
      canComb->SaveAs(outDir + "/" + "T1T2" +"_" + region[nTotPlots-1]+"ExpLimit_mMother_mLSP.C");

// Both obs and exp limits are in the same plots...
      T1Pad->cd();
      directCombLimitsVec[idxT1][0]->SetMinimum(0.001);
      directCombLimitsVec[idxT1][0]->SetMaximum(10);
      directCombLimitsVec[idxT1][0]->Draw("col");
      directCombDiagonalVec[idxT1][0]->Draw("l");

//      directCombProspinoVec[idxT1][0]->SetLineWidth(2); directCombProspinoVec[idxT1][0]->SetLineStyle(1); directCombProspinoVec[idxT1][0]->SetLineColor(kMagenta-4); directCombProspinoVec[idxT1][0]->Draw("l");
      directCombProspinoVec[idxT1][1]->SetLineWidth(2); directCombProspinoVec[idxT1][1]->Draw("l");
//      directCombProspinoVec[idxT1][2]->SetLineWidth(2); directCombProspinoVec[idxT1][2]->SetLineStyle(1); directCombProspinoVec[idxT1][2]->SetLineColor(kMagenta+3); directCombProspinoVec[idxT1][2]->Draw("l");
      directCombSysErrBandsVec[idxT1][0]->SetLineWidth(1); directCombSysErrBandsVec[idxT1][0]->SetLineStyle(1); directCombSysErrBandsVec[idxT1][0]->SetLineColor(kBlack); directCombSysErrBandsVec[idxT1][0]->Draw("l");
      directCombSysErrBandsVec[idxT1][1]->SetLineWidth(1); directCombSysErrBandsVec[idxT1][1]->SetLineStyle(1); directCombSysErrBandsVec[idxT1][1]->SetLineColor(kBlack); directCombSysErrBandsVec[idxT1][1]->Draw("l");

/*
//      cmsPreTex.DrawLatex(0.20, 0.84, "CMS Preliminary");
      cmsPreTex.DrawLatex(0.20, 0.84, "CMS");
      cmsLumiTex.DrawLatex(0.20, 0.77, "4.98 fb^{-1}, #sqrt{s} = 7 TeV");
      disptTex.DrawLatex(0.20, 0.70, "Jets+#slash{H}_{T}");
*/
//      cmsPreTex.DrawLatex(0.21, 0.84, "CMS Preliminary");
      cmsPreTex.DrawLatex(0.21, 0.84, "CMS");
      cmsLumiTex.DrawLatex(0.21, 0.77, "4.98 fb^{-1}, #sqrt{s} = 7 TeV");
      disptTex.DrawLatex(0.21, 0.70, "Jets + #slash{H}_{T}");

      indexTex.DrawLatex(0.25, 0.50, "(a)");
//      indexTex.DrawLatex(0.29, 0.50, "(a)");

      directCombProspinoVec[idxT1][4]->SetLineWidth(2); directCombProspinoVec[idxT1][4]->SetLineColor(kMagenta+2); directCombProspinoVec[idxT1][4]->SetLineStyle(2); directCombProspinoVec[idxT1][4]->Draw("l");

      directCombExpBandsVec[idxT1][0]->SetLineWidth(1); directCombExpBandsVec[idxT1][0]->SetLineColor(kMagenta+2); directCombExpBandsVec[idxT1][0]->SetLineStyle(3); directCombExpBandsVec[idxT1][0]->Draw("l");
      directCombExpBandsVec[idxT1][1]->SetLineWidth(1); directCombExpBandsVec[idxT1][1]->SetLineColor(kMagenta+2); directCombExpBandsVec[idxT1][1]->SetLineStyle(3); directCombExpBandsVec[idxT1][1]->Draw("l");

      axisXmodifT1->Draw();

      T2Pad->cd();
      directCombLimitsVec[idxT2][0]->SetMinimum(0.001);
      directCombLimitsVec[idxT2][0]->SetMaximum(10);
      directCombLimitsVec[idxT2][0]->Draw("colz");
      directCombDiagonalVec[idxT2][0]->Draw("l");

      coverT2title->Draw(); titleT2Pave->Draw();

//      directCombProspinoVec[idxT2][0]->SetLineWidth(2); directCombProspinoVec[idxT2][0]->SetLineStyle(1); directCombProspinoVec[idxT2][0]->SetLineColor(kMagenta-4); directCombProspinoVec[idxT2][0]->Draw("l");
      directCombProspinoVec[idxT2][1]->SetLineWidth(2); directCombProspinoVec[idxT2][1]->Draw("l");
//      directCombProspinoVec[idxT2][2]->SetLineWidth(2); directCombProspinoVec[idxT2][2]->SetLineStyle(1); directCombProspinoVec[idxT2][2]->SetLineColor(kMagenta+3); directCombProspinoVec[idxT2][2]->Draw("l");
      directCombSysErrBandsVec[idxT2][0]->SetLineWidth(1); directCombSysErrBandsVec[idxT2][0]->SetLineStyle(1); directCombSysErrBandsVec[idxT2][0]->SetLineColor(kBlack); directCombSysErrBandsVec[idxT2][0]->Draw("l");
      directCombSysErrBandsVec[idxT2][1]->SetLineWidth(1); directCombSysErrBandsVec[idxT2][1]->SetLineStyle(1); directCombSysErrBandsVec[idxT2][1]->SetLineColor(kBlack); directCombSysErrBandsVec[idxT2][1]->Draw("l");

      directCombProspinoVec[idxT2][4]->SetLineWidth(2); directCombProspinoVec[idxT2][4]->SetLineColor(kMagenta+2); directCombProspinoVec[idxT2][4]->SetLineStyle(2); directCombProspinoVec[idxT2][4]->Draw("l");

      directCombExpBandsVec[idxT2][0]->SetLineWidth(1); directCombExpBandsVec[idxT2][0]->SetLineColor(kMagenta+2); directCombExpBandsVec[idxT2][0]->SetLineStyle(3); directCombExpBandsVec[idxT2][0]->Draw("l");
      directCombExpBandsVec[idxT2][1]->SetLineWidth(1); directCombExpBandsVec[idxT2][1]->SetLineColor(kMagenta+2); directCombExpBandsVec[idxT2][1]->SetLineStyle(3); directCombExpBandsVec[idxT2][1]->Draw("l");

      legCombT2->Clear();
      legCombT2->AddEntry(directCombProspinoVec[idxT2][1], "#sigma^{NLO+NLL} #pm 1#sigma sig. theory", "l");
      legCombT2->AddEntry(directCombProspinoVec[idxT2][4], "Exp. limit #pm 1#sigma exp.", "l");
////      legCombT2->AddEntry(directCombExpBandsVec[idxT2][0], "Exp. limit #pm 1#sigma exp.", "l");
      legCombT2->AddEntry(directCombProspinoVec[idxT2][0], " ", "");
      legCombT2->AddEntry(directCombProspinoVec[idxT2][2], " ", "");
//      legCombT2->AddEntry(directCombProspinoVec[idxT2][0], "1/3#times#sigma^{NLO+NLL}", "l");
//      legCombT2->AddEntry(directCombProspinoVec[idxT2][2], "3#times#sigma^{NLO+NLL}", "l");
      legCombT2->Draw();

      TLegend* legObsSysErr1 = new TLegend(cacheX1legexp-changeXNDClegexp, legexp->GetY1NDC()+0.10, cacheX2legexp, legexp->GetY2NDC()+0.05, NULL,"brNDC");
      legObsSysErr1->SetFillColor(0);legObsSysErr1->SetShadowColor(0);legObsSysErr1->SetFillStyle(4000);legObsSysErr1->SetTextFont(42);legObsSysErr1->SetTextSize(0.025*scaleLegendFont);legObsSysErr1->SetBorderSize(0);
      legObsSysErr1->AddEntry(directCombSysErrBandsVec[idxT2][0], " ", "l");
      legObsSysErr1->Draw();

      TLegend* legObsSysErr2 = new TLegend(cacheX1legexp-changeXNDClegexp, legexp->GetY1NDC()+0.10, cacheX2legexp, legexp->GetY2NDC()-0.02, NULL,"brNDC");
      legObsSysErr2->SetFillColor(0);legObsSysErr2->SetShadowColor(0);legObsSysErr2->SetFillStyle(4000);legObsSysErr2->SetTextFont(42);legObsSysErr2->SetTextSize(0.025*scaleLegendFont);legObsSysErr2->SetBorderSize(0);
      legObsSysErr2->AddEntry(directCombSysErrBandsVec[idxT2][1], " ", "l");
      legObsSysErr2->Draw();

      TLegend* legExpBand1 = new TLegend(cacheX1legexp-changeXNDClegexp, legexp->GetY1NDC()+0.05, cacheX2legexp, legexp->GetY2NDC()-0.025, NULL,"brNDC");
      legExpBand1->SetFillColor(0);legExpBand1->SetShadowColor(0);legExpBand1->SetFillStyle(4000);legExpBand1->SetTextFont(42);legExpBand1->SetTextSize(0.025*scaleLegendFont);legExpBand1->SetBorderSize(0);
      legExpBand1->AddEntry(directCombExpBandsVec[idxT2][0], " ", "l");
      legExpBand1->Draw();

      TLegend* legExpBand2 = new TLegend(cacheX1legexp-changeXNDClegexp, legexp->GetY1NDC()+0.05, cacheX2legexp, legexp->GetY2NDC()-0.105, NULL,"brNDC");
      legExpBand2->SetFillColor(0);legExpBand2->SetShadowColor(0);legExpBand2->SetFillStyle(4000);legExpBand2->SetTextFont(42);legExpBand2->SetTextSize(0.025*scaleLegendFont);legExpBand2->SetBorderSize(0);
      legExpBand2->AddEntry(directCombExpBandsVec[idxT2][1], " ", "l");
      legExpBand2->Draw();

//      indexTex.DrawLatex(0.25*(l_+x_), 0.50, "(b)");
      indexTex.DrawLatex(0.25*(l_+x_)-0.04, 0.50, "(b)");

      canComb->SaveAs(outDir + "/" + "T1T2" +"_" + region[nTotPlots-1]+"allINone.pdf");
      canComb->SaveAs(outDir + "/" + "T1T2" +"_" + region[nTotPlots-1]+"allINone.png");
      canComb->SaveAs(outDir + "/" + "T1T2" +"_" + region[nTotPlots-1]+"allINone.C");

   }
  
}

int main(int argc, char** argv)
{
  return plot(argc, argv);
}
