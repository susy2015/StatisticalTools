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

#include <string>
#include <cmath>
#include <stdio.h>
#include <sys/stat.h>

static const int nRegions = 5;
static const string  disptStrs[] = {"H_{T}>500, #slash{H}_{T}>350 GeV", "H_{T}>800, #slash{H}_{T}>200 GeV", "H_{T}>800, #slash{H}_{T}>500 GeV", "H_{T}>1000, #slash{H}_{T}>400 GeV", "H_{T}>1200, #slash{H}_{T}>200 GeV"};
static const string regionDirs[] = {"mediumHTMHT", "highHT", "highHTMHT", "HT1000MHT400", "HT1200MHT200"};
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
      yPtsDiag[ib] = loDiagRange + (ib+1)*widthBin; yErrsDiag[ib] = 0.0;
   }
   diagonalGraph = new TGraphErrors(nDiagBins, xPtsDiag, yPtsDiag, xErrsDiag, yErrsDiag);
   diagonalGraph->SetLineWidth(2); diagonalGraph->SetLineStyle(7);

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

   TString outDir = "smoothed_plots_"+topoStr;

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

   TH2F *minCatsObsExclNLO, *minCatsExpExclNLO;
   TH2F *minObsExclNLO, *minExpExclNLO, *minObsExclNLOfollowminExp;

   //Was not originally designed as a loop, so this will give you warnings about possible memory leaks.

//   TLegend* legexp = new TLegend(0.15,0.63,0.48,0.71,NULL,"brNDC"); 
   TLegend* legexp = new TLegend(0.18,0.70,0.48,0.85,NULL,"brNDC"); 
   legexp->SetFillColor(0);legexp->SetShadowColor(0);legexp->SetFillStyle(4000);legexp->SetTextFont(42);legexp->SetTextSize(0.025*scaleLegendFont);legexp->SetBorderSize(0);

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

   //the histograms ---------------------------------------------------------------------
   // cross-section in M0 - M1/2
         c1->SetLogz(1);
         c1->SetRightMargin ( 0.19 );
         c1->SetTopMargin ( 0.085 );

   // Observed Limit in M0 - M1/2
         TH2F*hobslimit = 0;
         if( topoStr == "T1" ){
            hobslimit = new TH2F("obslimit",";m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}#; m(#tilde{q})>>m(#tilde{g})");
         }else if( topoStr == "T2" ){
            hobslimit = new TH2F("obslimit",";m_{#tilde{q}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{q}#tilde{q}, #tilde{q}#rightarrow q#tilde{#chi}^{0}#; m(#tilde{g})>>m(#tilde{q})");
         }else if( topoStr == "T2tt" ){
            hobslimit = new TH2F("obslimit",";m_{#tilde{t}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow q#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
         }else if( topoStr == "TGQ" ){
            hobslimit = new TH2F("obslimit",";m_{#tilde{g}} [GeV]; m_{#tilde{q}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{q}#tilde{g}, #tilde{q}#rightarrow q#tilde{#chi}^{0}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}");
         }else if( topoStr == "T5ZZInc" ){
            hobslimit = new TH2F("obslimit",";m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Observed Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hobslimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}_{2}, #tilde{#chi}^{0}_{2}#rightarrow Z#tilde{#chi}^{0}");
         }
         plotTools[i]->Area(hobslimit, Mzero, Mhalf, NLOObsxSecCL, false);
         obslimitOriVec[i] = (TH2F*) hobslimit->Clone();
         Smooth(hobslimit, 2);
         hobslimit->SetMinimum(0.001);
         hobslimit->SetMaximum(200);
         hobslimit->Draw("colz");
         TString obsHistName = region[i]+"obsLimit";
         obslimitVec[i] = (TH2F*) hobslimit->Clone(obsHistName.Data());

         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ObsLimit_mMother_mLSP.pdf");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ObsLimit_mMother_mLSP.png");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ObsLimit_mMother_mLSP.C");
   
   // Expected Limit in M0 - M1/2
         TH2F*hexplimit = 0;
         if( topoStr == "T1" ){
            hexplimit = new TH2F("explimit",";m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}#; m(#tilde{q})>>m(#tilde{g})");
         }else if( topoStr == "T2" ){
            hexplimit = new TH2F("explimit",";m_{#tilde{q}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{q}#tilde{q}, #tilde{q}#rightarrow q#tilde{#chi}^{0}#; m(#tilde{g})>>m(#tilde{q})");
         }else if( topoStr == "T2tt" ){
            hexplimit = new TH2F("explimit",";m_{#tilde{t}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow q#tilde{#chi}^{0}#; m(#tilde{t})>m(#tilde{#chi}^{0})");
         }else if( topoStr == "TGQ" ){
            hexplimit = new TH2F("explimit",";m_{#tilde{g}} [GeV]; m_{#tilde{q}} [GeV]; 95% CL Expected Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{q}#tilde{g}, #tilde{q}#rightarrow q#tilde{#chi}^{0}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}");
         }else if( topoStr == "T5ZZInc"){
            hexplimit = new TH2F("explimit",";m_{#tilde{g}} [GeV]; m_{#tilde{#chi}^{0}} [GeV]; 95% CL Expected Limit [pb]", nXbins, loX, hiX, nYbins, loY, hiY);
            hexplimit->SetTitle("pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow qq#tilde{#chi}^{0}_{2}, #tilde{#chi}^{0}_{2}#rightarrow Z#tilde{#chi}^{0}");
         }
         plotTools[i]->Area(hexplimit, Mzero, Mhalf, NLOExpxSecCL, false);
         explimitOriVec[i] = (TH2F*) hexplimit->Clone();
         Smooth(hexplimit, 2);
         hexplimit->SetMinimum(0.001);
         hexplimit->SetMaximum(200);
         hexplimit->Draw("colz");
         TString expHistName = region[i]+"expLimit";
         explimitVec[i] = (TH2F*) hexplimit->Clone(expHistName.Data());

         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ExpLimit_mMother_mLSP.pdf");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ExpLimit_mMother_mLSP.png");
         c1->SaveAs(outDir + "/" +topoStr +"_" + region[i]+"ExpLimit_mMother_mLSP.C");

      }
   }

   TFile *outRootFile = new TFile(outDir+"/"+"CLs_SMS.root", "RECREATE");
   for(int i=0; i<nTotPlots; i++){
      if( i= nTotPlots -1 ) continue;
      obslimitOriVec[i]->Write(); explimitOriVec[i]->Write();
      obslimitVec[i]->Write(); explimitVec[i]->Write();
   }
   outRootFile->Write(); outRootFile->Close();

// Make new data cards using the smoothed plots, take the first one
   ifstream inlist; char line[200], resline[200];
   inlist.open(filename, ifstream::in);
   std::cout<<"\nfilename : "<<filename<<std::endl<<std::endl;

   int obsRelDiffSmtCnt =0, expRelDiffSmtCnt =0;
   while(inlist.getline(line, 200)){
      ifstream resfile;
      char newline[200];
      sprintf(newline, "%s%s", dirname, line);
      resfile.open(newline, ifstream::in);

      int Mzero = -1, Mhalf = -1, Xsection = -1;
      double obs = -1, exp = -1, obserror = -1;
      while(resfile.getline(resline, 200)){
         TString reslineT(resline);
         TObjArray *vlist = reslineT.Tokenize("= ");
         int vlistEntries = vlist->GetEntries();

         TObjString* firstObj = dynamic_cast<TObjString*>(vlist->At(0)); TString firstStr = firstObj->GetString();
         TObjString* lastObj = dynamic_cast<TObjString*>(vlist->At(vlistEntries-1)); TString lastStr = lastObj->GetString();
         if(firstStr.Contains("Mzero")){ Mzero = lastStr.Atoi(); }
         if(firstStr.Contains("Mhalf")){ Mhalf = lastStr.Atoi(); }
         if(firstStr.Contains("Xsection")){ Xsection = lastStr.Atoi(); }

         if(!firstStr.CompareTo("limit.cls.observed")){ obs = lastStr.Atof(); }
         if(!firstStr.CompareTo("limit.cls.observed.error")){ obserror = lastStr.Atof(); }
         if(!firstStr.CompareTo("limit.cls.expected")){ exp = lastStr.Atof(); }

      }
      int Azero = 0, Mu = 1;

      int obsIdxOri = obslimitOriVec[0]->FindFixBin(Mzero, Mhalf);
      int obsIdxSmt = obslimitVec[0]->FindFixBin(Mzero, Mhalf);

      double obsValOri = obslimitOriVec[0]->GetBinContent(obsIdxOri);
      double obsValSmt = obslimitVec[0]->GetBinContent(obsIdxSmt);

      double obsRelDiff = 0, obsRelVery =0;
      if( obsValOri !=0 ) obsRelDiff = fabs(obsValOri-obsValSmt)/obsValOri;
      if( obsValOri !=0 ) obsRelVery = fabs(obs-obsValOri)/obsValOri;
      if( obsRelVery > 1e-6 ){
         std::cout<<"===> Diff between reading from histogram and datacards ?! for newline : "<<newline<<std::endl;
         printf("Mzero : %4d  Mhalf : %4d  obs : %10.8e  obsValOri : %10.8e  obsValSmt : %10.8e\n\n", Mzero, Mhalf, obs, obsValOri, obsValSmt);
      }
      if( obsRelDiff > 0.20 ){
         obsRelDiffSmtCnt ++;
//         std::cout<<"===> INFO obs change between original and smoothed > 0.20. for newline : "<<newline<<std::endl;
//         printf("Mzero : %4d  Mhalf : %4d  obs : %10.8e  obsValOri : %10.8e  obsValSmt : %10.8e\n\n", Mzero, Mhalf, obs, obsValOri, obsValSmt);
      }

      int expIdxOri = explimitOriVec[0]->FindFixBin(Mzero, Mhalf);
      int expIdxSmt = explimitVec[0]->FindFixBin(Mzero, Mhalf);

      double expValOri = explimitOriVec[0]->GetBinContent(expIdxOri);
      double expValSmt = explimitVec[0]->GetBinContent(expIdxSmt);

      double expRelDiff = 0, expRelVery =0;
      if( expValOri !=0 ) expRelDiff = fabs(expValOri-expValSmt)/expValOri;
      if( expValOri !=0 ) expRelVery = fabs(exp-expValOri)/expValOri;
      if( expRelVery > 1e-6 ){
         std::cout<<"===> Diff between reading from histogram and datacards ?! for newline : "<<newline<<std::endl;
         printf("Mzero : %4d  Mhalf : %4d  exp : %10.8e  expValOri : %10.8e  expValSmt : %10.8e\n\n", Mzero, Mhalf, exp, expValOri, expValSmt);
      }
      if( expRelDiff > 0.20 ){
         expRelDiffSmtCnt ++;
//         std::cout<<"===> INFO exp change between original and smoothed > 0.20. for newline : "<<newline<<std::endl;
//         printf("Mzero : %4d  Mhalf : %4d  exp : %10.8e  expValOri : %10.8e  expValSmt : %10.8e\n\n", Mzero, Mhalf, exp, expValOri, expValSmt);
      }

      char outnames[200];
      sprintf(outnames, "%s/%s", outDir.Data(), line);
      ofstream outfile;
      outfile.open(outnames);
      outfile<<"Azero = "<<Azero<<std::endl; 
      outfile<<"Mu = "<<Mu<<std::endl;
      outfile<<"Mzero = "<<Mzero<<std::endl;
      outfile<<"Mhalf = "<<Mhalf<<std::endl;
      outfile<<"Xsection = "<<Xsection<<std::endl;
      outfile<<"limit.cls.observed = "<<obsValSmt<<std::endl;
      outfile<<"limit.cls.observed.error = "<<obserror<<std::endl;
      outfile<<"limit.cls.expected = "<<expValSmt<<std::endl;
      outfile<<"limit.cls.expected.1sigmaCoverage = -1"<<std::endl;
      outfile<<"limit.cls.expected.2sigmaCoverage = -1"<<std::endl;
      outfile<<"limit.cls.expected.m1sigma = -1"<<std::endl;
      outfile<<"limit.cls.expected.m2sigma = -1"<<std::endl;
      outfile<<"limit.cls.expected.p1sigma = -1"<<std::endl;
      outfile<<"limit.cls.expected.p2sigma = -1"<<std::endl;
      outfile.close();
    } 
    std::cout<<"obsRelDiffSmtCnt : "<<obsRelDiffSmtCnt<<std::endl;
    std::cout<<"expRelDiffSmtCnt : "<<expRelDiffSmtCnt<<std::endl<<std::endl;

    char command[200];
    sprintf(command, "cp %s %s", filename, outDir.Data());
    system(command);
}

int main(int argc, char** argv)
{
  return plot(argc, argv);
}
