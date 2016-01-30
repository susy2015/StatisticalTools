#include "plot.h"
#include "SusyScan.h"
#include "GeneratorMasses.h"
#include "simpleGenMasses.h"
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
static const int     inclFlags[] = {     1,           1,          1,           0,               0       };

static const int      colors[] = {kGreen, kMagenta+1, kBlue, kRed, kTeal+4};

// For m0-m1/2 plane
TGraph *LEP_ch_m0m12, *LEP_sl_m0m12, *TEV_sg_cdf_m0m12, *TEV_sg_d0_m0m12, *stau_m0m12;
TGraph *Atlas_m0m12, *Observed2010_m0m12;
TGraphErrors *RA1_m0m12;
TGraph * LEP_ch_m0m12_cloned;

// For guino-squark plane
TGraph *LEP_glsq, *TEV_glsq, *CDF_glsq, *DEZ_glsq, *NoSol_glsq, *NoSol_glsq_aux, *NoSol_glsq_Cloned;

int glb_xnbins, glb_ynbins;
double glb_x_lo, glb_x_hi, glb_y_lo, glb_y_hi;
char plotKeyStr[200];
TString plotKeyStrT;

int tanbeta_input = 10;
int dom0m12 = 1;
bool debug = false;

void loadAUXm0m12(int tanbeta=10){

   std::cout<<"\nloadAUXm0m12 ... "<<std::endl;

   if( tanbeta == 10 ){
      RA1_m0m12 = RA1_NLO();
      RA1_m0m12->SetLineColor(kRed+1);
      RA1_m0m12->SetLineStyle(4);
      RA1_m0m12->SetLineWidth(3);

      Observed2010_m0m12 = Obs2010_NLO();
      Observed2010_m0m12->SetLineColor(kBlue);
      Observed2010_m0m12->SetLineWidth(1);

      Atlas_m0m12 = Atlas_m0_m12_tb3_obs();
      Atlas_m0m12->SetLineColor( 28 );

      LEP_ch_m0m12 = set_lep_ch(tanbeta);
      LEP_sl_m0m12 = set_lep_sl(tanbeta);//slepton curve
      TEV_sg_cdf_m0m12 = set_tev_sg_cdf(tanbeta);//squark gluino cdf
      TEV_sg_d0_m0m12 = set_tev_sg_d0(tanbeta);//squark gluino d0
      stau_m0m12 = set_tev_stau(tanbeta);//stau 
   }
   if( tanbeta == 40 ){
      std::cout<<"WARNING ... no setting for RA1, Observed2010 and Atlas ..."<<std::endl;
      RA1_m0m12 = RA1_NLO();
      RA1_m0m12->SetLineColor(kRed+1);
      RA1_m0m12->SetLineStyle(4);
      RA1_m0m12->SetLineWidth(3);

      Observed2010_m0m12 = Obs2010_NLO();
      Observed2010_m0m12->SetLineColor(kBlue);
      Observed2010_m0m12->SetLineWidth(1);

      Atlas_m0m12 = Atlas_m0_m12_tb3_obs();
      Atlas_m0m12->SetLineColor( 28 );

      LEP_ch_m0m12 = set_lep_ch(10);
      LEP_sl_m0m12 = set_lep_sl(10);//slepton curve
      TEV_sg_cdf_m0m12 = set_tev_sg_cdf(10);//squark gluino cdf
      TEV_sg_d0_m0m12 = set_tev_sg_d0(10);//squark gluino d0
      stau_m0m12 = set_tev_stau(10);//stau 
   }

   LEP_ch_m0m12_cloned = (TGraph*) LEP_ch_m0m12->Clone(); 
   LEP_ch_m0m12_cloned->SetFillStyle(1001); LEP_ch_m0m12_cloned->SetFillColor(kWhite);
}

void loadAUXglsq(int tanbeta=10){

   std::cout<<"\nloadAUXglsq ... "<<std::endl;

   LEP_glsq = gl_LEP();
   TEV_glsq = sq_TEV();
   CDF_glsq = sq_CDF();
   DEZ_glsq = sq_DEZ();
   NoSol_glsq = glsq_NoSol();
   NoSol_glsq_aux = glsq_NoSol_aux();
   NoSol_glsq_Cloned = (TGraph*) NoSol_glsq->Clone();
   NoSol_glsq_Cloned->SetFillStyle(1001); NoSol_glsq_Cloned->SetFillColor(kWhite);
}

int plot(int argc, char** argv)
{
   //interactive root session
   //TApplication theApp("App", 0, 0);
   if (gROOT->IsBatch()) {
      fprintf(stderr, "%s: cannot run in batch mode\n", argv[0]);
      return 1;
   }

   if( argc !=3 ){
      std::cout<<"\n################USAGE##########################"<<std::endl;
      std::cout<<  "#./plot_mSUGRA dom0m12 tanbeta                #"<<std::endl;
      std::cout<<  "# dom0m12  1: m0-m12   0 : gluino-squark      #"<<std::endl;
      std::cout<<  "# tanbeta  10 or 40                           #"<<std::endl;
      std::cout<<  "################USAGE##########################"<<std::endl;
      std::cout<<std::endl;
      return 1;
   }

   dom0m12 = atoi( argv[1] );
   tanbeta_input = atoi( argv[2] );
   std::cout<<"\nWorking on dom0m12 : "<<dom0m12<<"  tanbeta : "<<tanbeta_input<<std::endl<<std::endl;

   util::StyleSettings::paperNoTitle();
   gStyle->SetPadBottomMargin(0.18);
 
   TCanvas * c1 = new TCanvas("c1","c1",600,600);
   c1->cd();
   
   //Get limits from signal scan ---------------------------------------------------
   char filename[200], dirname[200], tmpStr[200];

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
         sprintf(dirname, "%s_plots_tanbeta%d/", regionDirs[ir].c_str(), tanbeta_input);
         sprintf(filename, "%s_plots_tanbeta%d/filelist.txt", regionDirs[ir].c_str(), tanbeta_input);
         idxCnt++; genpoints[idxCnt]->Fill(filename, dirname);
//         if( tanbeta_input == 10 ){ if( ir==0 || ir==1 || ir==2 ) genpoints[idxCnt]->FillEmptyPointsByInterpolationInM0M12();}
//Get SUSY masses and x-sections from generator scan ----------------------------
//and match them to the signal scan
         genpoints[idxCnt]->FillGeneratorMasses("GenScan_tb10.dat");
         genpoints[idxCnt]->match();
//Read in SUSY masses in simpler format (produced by Steve Mrenna)
//and match them to the signal scan
         genpoints[idxCnt]->FillSimpleGenMasses("steveMrenna_scanInputs.txt");
         genpoints[idxCnt]->matchSimple();
      }
   }

   TString outDir = "results";

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

   loadAUXm0m12(tanbeta_input); loadAUXglsq(tanbeta_input);

   //Was not originally designed as a loop, so this will give you warnings about possible memory leaks.
   for(int i = 0; i < nTotPlots; i++){
      if(i < nTotPlots-1){
   //the plotting ----------------------------------------------------------------------
   //plotting helper functions
         plotTools[i] = new PlotTools<SusyScan>(genpoints[i]->GetScan());
         PlotTools<GeneratorMasses> * plotMasses = new PlotTools<GeneratorMasses>(genpoints[i]->GetGeneratorMasses());

   //the histograms ---------------------------------------------------------------------
         c1->SetLogz(0);
         c1->SetRightMargin ( 0.1 );
         c1->SetTopMargin ( 0.09 );

   // Histograms for combination plot later ==================================================================
         if( dom0m12 ){
            std::cout<<"Get contour of "<<regionBase[i]<<" for m0-m1/2"<<std::endl;
            int nm0 = 100, nm12 = 38;
            double m0_lo= 0, m0_hi = 1850.9, m12_lo = 0, m12_hi = 770.;
            if( tanbeta_input == 40 ){
               m0_lo= 0, m0_hi = 1827.0, m12_lo = 0, m12_hi = 770.; 
            }
/*// For old tanbeta10 results
            int nm0 = 100, nm12 = 38;
            double m0_lo= 0, m0_hi = 1870.9, m12_lo = 0, m12_hi = 770.;
            if( tanbeta_input == 40 ){
               m0_lo= 0, m0_hi = 1827.0, m12_lo = 0, m12_hi = 770.; 
            }
*/
            ObsExcl_NLO[i] = new TH2F(region[i]+"obs",";m_{0} [GeV]; m_{1/2} [GeV]; 95% CL Expected Exclusion", nm0, m0_lo, m0_hi, nm12, m12_lo, m12_hi);
            plotTools[i]->Area(ObsExcl_NLO[i],Mzero,Mhalf,NLOObsExclCL);
            ExpExcl_NLO[i] = new TH2F(region[i]+"exp",";m_{0} [GeV]; m_{1/2} [GeV]; 95% CL Expected Exclusion", nm0, m0_lo, m0_hi, nm12, m12_lo, m12_hi);
            plotTools[i]->Area(ExpExcl_NLO[i],Mzero,Mhalf,NLOExpExclCL);
            ExpExclm1[i] = new TH2F(region[i]+"expm1",";m_{0} [GeV]; m_{1/2} [GeV]; 95% CL Expected Exclusion", nm0, m0_lo, m0_hi, nm12, m12_lo, m12_hi);
            plotTools[i]->Area(ExpExclm1[i],Mzero,Mhalf,NLOExpExclCLm1sigma);
            ExpExclp1[i] = new TH2F(region[i]+"expp1",";m_{0} [GeV]; m_{1/2} [GeV]; 95% CL Expected Exclusion", nm0, m0_lo, m0_hi, nm12, m12_lo, m12_hi);
            plotTools[i]->Area(ExpExclp1[i],Mzero,Mhalf,NLOExpExclCLp1sigma);

   // Exclusion in M0 - M1/2 =================================================================================
   // hs defined range in which the contrours are calculated, the binning matters! binning should be equal to signal-scan 
// For old tanbeta10 results
//            nm0 = 95; m0_lo = 100;
            TH2F*hs = new TH2F(region[i]+"hs","",nm0, m0_lo, m0_hi, nm12, m12_lo, m12_hi);
   // usage: GetContour( <range-hist>, <x-var>, <y-var>, <limit-function>, <contour-levels>, <contour index>, <color>, <style> )
            gCLsExpExcl_NLO[i]     = plotTools[i]->GetContour(hs,Mzero,Mhalf,NLOExpExclCL, 3,0, 1,2); 
            gCLsObsExcl_NLO[i]     = plotTools[i]->GetContour(hs,Mzero,Mhalf,NLOObsExclCL, 3,0, 1,1);
            gCLsExpExclm1[i]       = plotTools[i]->GetContour(hs,Mzero,Mhalf,NLOExpExclCLm1sigma,3,0, 5,2); 
            gCLsExpExclp1[i]       = plotTools[i]->GetContour(hs,Mzero,Mhalf,NLOExpExclCLp1sigma,3,0, 5,2);
         }else{
            std::cout<<"Get contour of "<<regionBase[i]<<" for gluino-squark"<<std::endl;
            int ngl = 100, nsq = 100;
            double gl_lo= 0, gl_hi = 1855., sq_lo = 0, sq_hi = 1855.;
            if( tanbeta_input == 40 ){
               gl_lo= 0, gl_hi = 1870., sq_lo = 0, sq_hi = 1870.;
            }
/* // For the old tanbeta=10 results
            int ngl = 100, nsq = 100;
            double gl_lo= 0, gl_hi = 1870., sq_lo = 0, sq_hi = 1870.;
            if( tanbeta_input == 40 ){
               gl_lo= 0, gl_hi = 1870., sq_lo = 0, sq_hi = 1870.;
            }
*/
            ObsExcl_NLO[i] = new TH2F(region[i]+"obs",";m_{#tilde{g}} (GeV); m_{#tilde{q}} (GeV); 95% CL Expected Exclusion", ngl, gl_lo, gl_hi, nsq, sq_lo, sq_hi);
            plotTools[i]->Area(ObsExcl_NLO[i],MGluino,MSquark,NLOObsExclCL);
            ExpExcl_NLO[i] = new TH2F(region[i]+"exp",";m_{#tilde{g}} (GeV); m_{#tilde{q}} (GeV); 95% CL Expected Exclusion", ngl, gl_lo, gl_hi, nsq, sq_lo, sq_hi);
            plotTools[i]->Area(ExpExcl_NLO[i],MGluino,MSquark,NLOExpExclCL);
            ExpExclm1[i] = new TH2F(region[i]+"expm1",";m_{#tilde{g}} (GeV); m_{#tilde{q}} (GeV); 95% CL Expected Exclusion", ngl, gl_lo, gl_hi, nsq, sq_lo, sq_hi);
            plotTools[i]->Area(ExpExclm1[i],MGluino,MSquark,NLOExpExclCLm1sigma);
            ExpExclp1[i] = new TH2F(region[i]+"expp1",";m_{#tilde{g}} (GeV); m_{#tilde{q}} (GeV); 95% CL Expected Exclusion", ngl, gl_lo, gl_hi, nsq, sq_lo, sq_hi);
            plotTools[i]->Area(ExpExclp1[i],MGluino,MSquark,NLOExpExclCLp1sigma);

            TH2F*hs = new TH2F(region[i]+"hs","",ngl, gl_lo, gl_hi, nsq, sq_lo, sq_hi);
            gCLsExpExcl_NLO[i]     = plotTools[i]->GetContour(hs,MGluino,MSquark,NLOExpExclCL, 3,0, 1,2); 
            gCLsObsExcl_NLO[i]     = plotTools[i]->GetContour(hs,MGluino,MSquark,NLOObsExclCL, 3,0, 1,1);
            gCLsExpExclm1[i]       = plotTools[i]->GetContour(hs,MGluino,MSquark,NLOExpExclCLm1sigma,3,0, 5,2); 
            gCLsExpExclp1[i]       = plotTools[i]->GetContour(hs,MGluino,MSquark,NLOExpExclCLp1sigma,3,0, 5,2);
         }
      }

      TH2F*TestCont = new TH2F(region[i]+"texcl",";m_{0} [GeV]; m_{1/2} [GeV]; 95% CL Expected Exclusion", 100,0,1800,100,0,1800);
      if( i == nTotPlots-1 ){
         TestCont = (TH2F*)ObsExcl_NLO[i];
         plotTools[i]->GetContour(TestCont,3,0);
      }
      if( i < nTotPlots-1 ) plotTools[i]->Area(TestCont,MGluino,MSquark,NLOObsExclCL,false);
      std::vector<TGraph*> contours = plotTools[i]->GetContours(TestCont,3);
      TestCont->Draw("colz");
      if( dom0m12 ) sprintf(plotKeyStr, "ExclusionTestContours_m0_m12_tb%d", tanbeta_input);
      else sprintf(plotKeyStr, "ExclusionTestContours_gl_sq_tb%d", tanbeta_input);
      if( debug ){
        c1->SaveAs("results/"+region[i]+plotKeyStr+".pdf");
        c1->SaveAs("results/"+region[i]+plotKeyStr+".png");
        c1->SaveAs("results/"+region[i]+plotKeyStr+".C");
      }

     // hexcl is drawn to define axis range and labels:
      TH2F *hexcl;
      if( dom0m12 ){
         glb_xnbins = 100; glb_ynbins = 45;
         glb_x_lo = 0; glb_x_hi = 1800; glb_y_lo = 110; glb_y_hi = 700;
         hexcl = new TH2F(region[i]+"hexcl",";m_{0} (GeV); m_{1/2} (GeV); 95% CL Exclusion", glb_xnbins, glb_x_lo, glb_x_hi, glb_ynbins, glb_y_lo, glb_y_hi);
      }else{
         glb_xnbins = 100; glb_ynbins = 100;
         glb_x_lo = 0; glb_x_hi = 1700; glb_y_lo = 0; glb_y_hi = 1700;
         hexcl = new TH2F(region[i]+"hexcl",";m_{#tilde{g}} (GeV); m_{#tilde{q}} (GeV); 95% CL Exclusion", glb_xnbins, glb_x_lo, glb_x_hi, glb_ynbins, glb_y_lo, glb_y_hi);
      }

      if( dom0m12 ){
         hexcl->SetTitle("CMS Preliminary");
   
         hexcl->GetYaxis()->SetTitleOffset(1.3);
         hexcl->GetXaxis()->SetTitleOffset(0.92);
         hexcl->GetYaxis()->SetLabelSize(0.03);
         hexcl->GetXaxis()->SetLabelSize(0.03);
         hexcl->Draw("colz");
    
         TLatex ms; ms.SetTextSize(0.025); ms.SetTextFont(42);//ms.SetTextColor(12);
         sprintf(tmpStr, "tan#beta=%d, #mu>0, A_{0}=0", tanbeta_input); ms.DrawLatex(450,540,tmpStr); 
    
         gCLsObsExcl_NLO[i]->SetLineWidth(2);
         gCLsExpExcl_NLO[i]->SetLineWidth(2);
       
       //smooth contours (2D Gaussian smoothing)
         Smooth( gCLsExpExcl_NLO[i], 27 );
         Smooth( gCLsObsExcl_NLO[i], 27 );
         Smooth( gCLsExpExclm1[i], 27 );
         Smooth( gCLsExpExclp1[i], 27 );
    
       //"Combine" HT and  limit by choosing the better one
         TGraph * gCLsExp1Sigma = MakeBand(gCLsExpExclp1[i], gCLsExpExclm1[i]);gCLsExp1Sigma->SetFillStyle(4010);
         gCLsExp1Sigma->Draw("F");
    
         LEP_ch_m0m12->Draw("fsame");
         LEP_sl_m0m12->Draw("fsame");
         TEV_sg_cdf_m0m12->Draw("fsame");
         TEV_sg_cdf_m0m12->Draw("lsame");
         TEV_sg_d0_m0m12->Draw("fsame");
         TEV_sg_d0_m0m12->Draw("lsame");
    
         TLegend* legexp = new TLegend(0.61,0.64,0.93,0.89,NULL,"brNDC");    legexp->SetFillColor(0);legexp->SetShadowColor(0);legexp->SetFillStyle(4000);legexp->SetTextFont(42);legexp->SetTextSize(0.025);legexp->SetBorderSize(0);
         legexp->AddEntry(TEV_sg_cdf_m0m12,"CDF  #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=5, #mu<0}","f"); 
         legexp->AddEntry(TEV_sg_d0_m0m12,"D0   #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=3, #mu<0}","f");  
         legexp->AddEntry(LEP_ch_m0m12,"LEP2   #tilde{#chi}_{1}^{#pm}","f");   
         legexp->AddEntry(LEP_sl_m0m12,"LEP2   #tilde{#font[12]{l}}^{#pm}","f"); //NOT FOR tb=50!
         legexp->AddEntry(RA1_m0m12, "#alpha_{T} CMS 1.1 fb^{-1}","l");   
         legexp->AddEntry(Observed2010_m0m12,"Observed 2010","l");
         legexp->Draw();
    
         RA1_m0m12->Draw("C");
         //Old RA2 contour:
         Observed2010_m0m12->Draw("l");
    
         //draw RA2 limit contours:
         gCLsObsExcl_NLO[i]->Draw("l");
         // Error band needs to move to top so it doesn't cover other things.
         gCLsExpExcl_NLO[i]->Draw("l");
    
         stau_m0m12->Draw("fsame");
         TLatex b; b.SetTextSize(0.02); b.SetTextColor(1); b.SetTextAngle(+82);
         b.DrawLatex( 80,570,"#tilde{#tau} LSP"); 
    
         //constant ssqquark and gluino lines
         TF1* lnsq[4];
         TF1* lngl[4];
         TLatex sqt; sqt.SetTextSize(0.02); sqt.SetTextColor(kGray+2); 
         sqt.SetTextAngle(-28); sqt.DrawLatex(148,230,"#font[92]{#tilde{q}(500)GeV}");
         sqt.SetTextAngle(-22); sqt.DrawLatex(148,360,"#font[92]{#tilde{q}(750)GeV}");
         sqt.SetTextAngle(-16); sqt.DrawLatex(148,485,"#font[92]{#tilde{q}(1000)GeV}");
         TLatex glt; glt.SetTextSize(0.02); sqt.SetTextAngle(-4); glt.SetTextColor(kGray+2);
         glt.DrawLatex(1200,170,"#font[92]{#tilde{g}(500)GeV}");
         glt.DrawLatex(1200,275,"#font[92]{#tilde{g}(750)GeV}");
         glt.DrawLatex(1200,385,"#font[92]{#tilde{g}(1000)GeV}");
         for(int j = 0; j < 4; j++){
            lnsq[j] = constant_squark(tanbeta_input,j);
            lngl[j] = constant_gluino(tanbeta_input,j);
            lngl[j]->Draw("same");   
            lnsq[j]->Draw("same");
         }
    
         TLegend * leg = new TLegend(0.25,0.77,0.6,0.89);
         leg->SetBorderSize(0);leg->SetFillColor(0);leg->SetFillStyle(4000);leg->SetTextFont(42);
         TGraph * expLeg = (TGraph*)gCLsExpExcl_NLO[i]->Clone();expLeg->SetFillStyle(gCLsExp1Sigma->GetFillStyle());expLeg->SetFillColor(gCLsExp1Sigma->GetFillColor());
         leg->SetHeader("L_{int} = 1.1 fb^{-1}, #sqrt{s} = 7 TeV");
         leg->AddEntry(gCLsObsExcl_NLO[i],"Observed","l");
         leg->AddEntry(expLeg,   "Expected #pm 1#sigma","lf");
         leg->Draw();
    
         gPad->RedrawAxis();
  
         sprintf(plotKeyStr, "Exclusion_m0_m12_tb%d", tanbeta_input);
      }else{
         gStyle->SetTitleBorderSize(0);
         gStyle->SetTitleAlign(13);
         gStyle->SetTitleX(0.25);
         gStyle->SetTitleY(0.96);
         gStyle->SetTitleH(0.025);
   
         sprintf(tmpStr, "L_{int} = 1.1 fb^{-1}, #sqrt{s} = 7 TeV              tan#beta=%d, #mu>0, A_{0}=0", tanbeta_input); 
         hexcl->SetTitle(tmpStr);
         hexcl->GetYaxis()->SetTitleOffset(1.3);
         hexcl->GetXaxis()->SetTitleOffset(0.92);
         hexcl->GetYaxis()->SetLabelSize(0.03);
         hexcl->GetXaxis()->SetLabelSize(0.03);
         hexcl->Draw("colz");
    
         gCLsObsExcl_NLO[i]->SetLineWidth(2);
         gCLsExpExcl_NLO[i]->SetLineWidth(2);
    
       //smooth contours (2D Gaussian smoothing)
         Smooth( gCLsExpExcl_NLO[i], 27 );
         Smooth( gCLsObsExcl_NLO[i], 27 );
         Smooth( gCLsExpExclm1[i], 27 );
         Smooth( gCLsExpExclp1[i], 27 );
       //"Combine" HT and  limit by choosing the better one
         TGraph * gCLsExp1Sigma = MakeBand(gCLsExpExclp1[i], gCLsExpExclm1[i]);gCLsExp1Sigma->SetFillStyle(4010);
         gCLsExp1Sigma->Draw("F");
    
         TGraph * expLeg = (TGraph*)gCLsExpExcl_NLO[i]->Clone();expLeg->SetFillStyle(gCLsExp1Sigma->GetFillStyle());expLeg->SetFillColor(gCLsExp1Sigma->GetFillColor());
    
         TLegend* legexp = new TLegend(0.19,0.65,0.39,0.90,NULL,"brNDC");
         legexp->SetFillColor(0);legexp->SetShadowColor(0); //legexp->SetFillStyle(4000);
         legexp->SetTextFont(42);legexp->SetTextSize(0.025);legexp->SetBorderSize(0);
         legexp->SetHeader("CMS preliminary");
         legexp->AddEntry(gCLsObsExcl_NLO[i],"Observed","l");
         legexp->AddEntry(expLeg,   "Expected #pm 1#sigma","lf");
         legexp->AddEntry(TEV_glsq, "Tevatron RunI", "f");
         legexp->AddEntry(CDF_glsq, "CDF RunII","f");
         legexp->AddEntry(DEZ_glsq, "D0 RunII","f");
         legexp->AddEntry(LEP_glsq, "LEP2","f");
    
         gCLsObsExcl_NLO[i]->Draw("l");
         // Error band needs to move to top so it doesn't cover other things.
         gCLsExpExcl_NLO[i]->Draw("l");
    
         NoSol_glsq_Cloned->Draw("fsame");
    
         DEZ_glsq->Draw("fsame");
         CDF_glsq->Draw("fsame");
         TEV_glsq->Draw("fsame");
         LEP_glsq->Draw("fsame");
         NoSol_glsq->Draw("fsame"); NoSol_glsq_aux->Draw("l");
    
         legexp->Draw();
    
         TLatex noSolTex; noSolTex.SetTextSize(0.03); noSolTex.SetTextColor(1);
         noSolTex.DrawLatex( 800, 600, "no CMSSM solution");
    
         gPad->RedrawAxis();
         sprintf(plotKeyStr, "Exclusion_gl_sq_tb%d", tanbeta_input);
      }
      if( i == nTotPlots-1 || debug ){
         c1->SaveAs("results/"+region[i]+plotKeyStr+".pdf");
         c1->SaveAs("results/"+region[i]+plotKeyStr+".png");
         c1->SaveAs("results/"+region[i]+plotKeyStr+".C");
      }
 
      // ===============================================================================================================
      if( i == nTotPlots-2 ){
        //combined exclusion for all five seach regions.
         ObsExcl_NLO[nTotPlots-1] = plotTools[i]->BinWiseOr(ObsExcl_NLO[0],ObsExcl_NLO[1]);
         for(int ic=2; ic<nTotPlots-1; ic++){
            ObsExcl_NLO[nTotPlots-1] = plotTools[i]->BinWiseOr(ObsExcl_NLO[nTotPlots-1], ObsExcl_NLO[ic]);
         }
         gCLsObsExcl_NLO[nTotPlots-1]     = plotTools[i]->GetContour(ObsExcl_NLO[nTotPlots-1],3,0);
         gCLsObsExcl_NLO[nTotPlots-1]->SetLineColor(1);
         gCLsObsExcl_NLO[nTotPlots-1]->SetLineStyle(1);
   
         ExpExcl_NLO[nTotPlots-1] = ExpExcl_NLO[0];
         ExpExclm1[nTotPlots-1] = ExpExclm1[0];
         ExpExclp1[nTotPlots-1] = ExpExclp1[0];
         for(int ic=1; ic<nTotPlots-1; ic++){
            plotTools[i]->BinWiseOr(ExpExcl_NLO[nTotPlots-1],ExpExcl_NLO[ic],ExpExclm1[nTotPlots-1],ExpExclm1[ic],ExpExclp1[nTotPlots-1],ExpExclp1[ic]);
         }
         gCLsExpExcl_NLO[nTotPlots-1]     = plotTools[i]->GetContour(ExpExcl_NLO[nTotPlots-1],3,0);
         gCLsExpExcl_NLO[nTotPlots-1]->SetLineColor(1);
         gCLsExpExcl_NLO[nTotPlots-1]->SetLineStyle(2);
        
         gCLsExpExclm1[nTotPlots-1] = plotTools[i]->GetContour(ExpExclm1[nTotPlots-1],3,0);
         gCLsExpExclm1[nTotPlots-1]->SetLineColor(5);
         gCLsExpExclm1[nTotPlots-1]->SetLineStyle(2);
        
         gCLsExpExclp1[nTotPlots-1] = plotTools[i]->GetContour(ExpExclp1[nTotPlots-1],3,0);
         gCLsExpExclp1[nTotPlots-1]->SetLineColor(5);
         gCLsExpExclp1[nTotPlots-1]->SetLineStyle(2);
      }
   }

// Draw the limit curves in one canvas for all regions
   c1->Clear();
   c1->SetRightMargin ( 0.1 );
   c1->SetTopMargin ( 0.09 );

   TH2F * hexcl = 0;
   if( dom0m12 ){
      hexcl = new TH2F("hexcl",";m_{0} (GeV); m_{1/2} (GeV); 95% CL Exclusion", glb_xnbins, glb_x_lo, glb_x_hi, glb_ynbins, glb_y_lo, glb_y_hi);
      hexcl->SetTitle("CMS Preliminary");
      hexcl->GetYaxis()->SetTitleOffset(1.3);
      hexcl->GetXaxis()->SetTitleOffset(0.92);
      hexcl->GetYaxis()->SetLabelSize(0.03);
      hexcl->GetXaxis()->SetLabelSize(0.03);
      hexcl->Draw("colz");
   }else{
      hexcl = new TH2F("hexcl",";m_{#tilde{g}} (GeV); m_{#tilde{q}} (GeV); 95% CL Exclusion", glb_xnbins, glb_x_lo, glb_x_hi, glb_ynbins, glb_y_lo, glb_y_hi);
      sprintf(tmpStr, "L_{int} = 1.1 fb^{-1}, #sqrt{s} = 7 TeV              tan#beta=%d, #mu>0, A_{0}=0", tanbeta_input); 
      hexcl->SetTitle(tmpStr);
      hexcl->GetYaxis()->SetTitleOffset(1.3);
      hexcl->GetXaxis()->SetTitleOffset(0.92);
      hexcl->GetYaxis()->SetLabelSize(0.03);
      hexcl->GetXaxis()->SetLabelSize(0.03);
      hexcl->Draw("colz");
   }

   TLegend * leg = 0; 

   if( dom0m12 ){
      leg = new TLegend(0.25,0.76,0.6,0.88);
      leg->SetBorderSize(0);leg->SetFillColor(0);leg->SetFillStyle(4000);leg->SetTextFont(42);leg->SetTextSize(0.020);
      leg->SetHeader("L_{int} = 1.1 fb^{-1}, #sqrt{s} = 7 TeV");
   }else{
      leg = new TLegend(0.50,0.78,0.82,0.88);
      leg->SetBorderSize(0);leg->SetFillColor(0);leg->SetFillStyle(4000);leg->SetTextFont(42);leg->SetTextSize(0.020);
   }

   for(int i = 0; i < nTotPlots; i++){
      if(i < nTotPlots-1){
         TGraph *obsCloned = (TGraph*) gCLsObsExcl_NLO[i]->Clone();
         obsCloned->SetLineColor(colors[i]);
         obsCloned->Draw("l same");
//         sprintf(tmpStr, "%s%s", region[i].Data(), "Obs");        
//         sprintf(tmpStr, "%s", regionBase[i].Data());
         sprintf(tmpStr, "Observed (%s)", dispt[i].Data());
         leg->AddEntry(obsCloned,tmpStr,"l");
      }
   }
   if( dom0m12 ){
      LEP_ch_m0m12_cloned->Draw("fsame");
      LEP_ch_m0m12->Draw("fsame");
      LEP_sl_m0m12->Draw("fsame");
      TEV_sg_cdf_m0m12->Draw("fsame");
      TEV_sg_cdf_m0m12->Draw("lsame");
      TEV_sg_d0_m0m12->Draw("fsame");
      TEV_sg_d0_m0m12->Draw("lsame");
      stau_m0m12->Draw("fsame");

      TLegend* legexp = new TLegend(0.61,0.70,0.93,0.88,NULL,"brNDC");    legexp->SetFillColor(0);legexp->SetShadowColor(0);legexp->SetFillStyle(4000);legexp->SetTextFont(42);legexp->SetTextSize(0.025);legexp->SetBorderSize(0);
      legexp->AddEntry(TEV_sg_cdf_m0m12,"CDF  #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=5, #mu<0}","f");
      legexp->AddEntry(TEV_sg_d0_m0m12,"D0   #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=3, #mu<0}","f");
      legexp->AddEntry(LEP_ch_m0m12,"LEP2   #tilde{#chi}_{1}^{#pm}","f");
      legexp->AddEntry(LEP_sl_m0m12,"LEP2   #tilde{#font[12]{l}}^{#pm}","f"); //NOT FOR tb=50!
      legexp->Draw();

      TLatex b; b.SetTextSize(0.02); b.SetTextColor(1); b.SetTextAngle(+82);
      b.DrawLatex( 80,570,"#tilde{#tau} LSP"); 

      TLatex ms; ms.SetTextSize(0.025); ms.SetTextFont(42);//ms.SetTextColor(12);
      sprintf(tmpStr, "tan#beta=%d, #mu>0, A_{0}=0", tanbeta_input); ms.DrawLatex(450,540,tmpStr); 
    
      sprintf(plotKeyStr, "Exclusion_m0_m12_tb%d", tanbeta_input); plotKeyStrT = TString(plotKeyStr);
   }else{
      NoSol_glsq_Cloned->Draw("fsame");

      DEZ_glsq->Draw("fsame");
      CDF_glsq->Draw("fsame");
      TEV_glsq->Draw("fsame");
      LEP_glsq->Draw("fsame");
      NoSol_glsq->Draw("fsame"); NoSol_glsq_aux->Draw("l");

      TLegend* legexp = new TLegend(0.19,0.48,0.35,0.73,NULL,"brNDC");
      legexp->SetFillColor(0);legexp->SetShadowColor(0); //legexp->SetFillStyle(4000);
      legexp->SetTextFont(42);legexp->SetTextSize(0.025);legexp->SetBorderSize(0);
      legexp->AddEntry(TEV_glsq, "Tevatron RunI", "f");
      legexp->AddEntry(CDF_glsq, "CDF RunII","f");
      legexp->AddEntry(DEZ_glsq, "D0 RunII","f");
      legexp->AddEntry(LEP_glsq, "LEP2","f");
      legexp->Draw();

      TLatex noSolTex; noSolTex.SetTextSize(0.03); noSolTex.SetTextColor(1);
      noSolTex.DrawLatex( 800, 600, "no CMSSM solution");
    
      sprintf(plotKeyStr, "Exclusion_gl_sq_tb%d", tanbeta_input); plotKeyStrT = TString(plotKeyStr);
   }
   leg->Draw();
   gPad->RedrawAxis();
   c1->SaveAs("results/comp_Obs_"+plotKeyStrT+".pdf");
   c1->SaveAs("results/comp_Obs_"+plotKeyStrT+".png");
   c1->SaveAs("results/comp_Obs_"+plotKeyStrT+".C");

   hexcl->Draw("colz");

   leg->Clear();
   if( dom0m12 ){
      leg->SetBorderSize(0);leg->SetFillColor(0);leg->SetFillStyle(4000);leg->SetTextFont(42);leg->SetTextSize(0.020);
      leg->SetHeader("L_{int} = 1.1 fb^{-1}, #sqrt{s} = 7 TeV");
   }else{
      leg->SetBorderSize(0);leg->SetFillColor(0);leg->SetFillStyle(4000);leg->SetTextFont(42);leg->SetTextSize(0.020);
   }

   for(int i = 0; i < nTotPlots; i++){
      if(i < nTotPlots-1){
         TGraph *expCloned = (TGraph*) gCLsExpExcl_NLO[i]->Clone();
         expCloned->SetLineColor(colors[i]);
         expCloned->Draw("l same");
//         sprintf(tmpStr, "%s%s", region[i].Data(), "Exp");
//         sprintf(tmpStr, "%s", regionBase[i].Data());
         sprintf(tmpStr, "Expected (%s)", dispt[i].Data());
         leg->AddEntry(expCloned,tmpStr,"l");
      }
   }
   if( dom0m12 ){
      LEP_ch_m0m12_cloned->Draw("fsame");
      LEP_ch_m0m12->Draw("fsame");
      LEP_sl_m0m12->Draw("fsame");
      TEV_sg_cdf_m0m12->Draw("fsame");
      TEV_sg_cdf_m0m12->Draw("lsame");
      TEV_sg_d0_m0m12->Draw("fsame");
      TEV_sg_d0_m0m12->Draw("lsame");
      stau_m0m12->Draw("fsame");

      TLegend* legexp = new TLegend(0.61,0.70,0.93,0.88,NULL,"brNDC");    legexp->SetFillColor(0);legexp->SetShadowColor(0);legexp->SetFillStyle(4000);legexp->SetTextFont(42);legexp->SetTextSize(0.025);legexp->SetBorderSize(0);
      legexp->AddEntry(TEV_sg_cdf_m0m12,"CDF  #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=5, #mu<0}","f");
      legexp->AddEntry(TEV_sg_d0_m0m12,"D0   #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=3, #mu<0}","f");
      legexp->AddEntry(LEP_ch_m0m12,"LEP2   #tilde{#chi}_{1}^{#pm}","f");
      legexp->AddEntry(LEP_sl_m0m12,"LEP2   #tilde{#font[12]{l}}^{#pm}","f"); //NOT FOR tb=50!
      legexp->Draw();

      TLatex b; b.SetTextSize(0.02); b.SetTextColor(1); b.SetTextAngle(+82);
      b.DrawLatex( 80,570,"#tilde{#tau} LSP");

      TLatex ms; ms.SetTextSize(0.025); ms.SetTextFont(42);//ms.SetTextColor(12);
      sprintf(tmpStr, "tan#beta=%d, #mu>0, A_{0}=0", tanbeta_input); ms.DrawLatex(450,540,tmpStr); 
    
      sprintf(plotKeyStr, "Exclusion_m0_m12_tb%d", tanbeta_input); plotKeyStrT = TString(plotKeyStr);
   }else{
      NoSol_glsq_Cloned->Draw("fsame");

      DEZ_glsq->Draw("fsame");
      CDF_glsq->Draw("fsame");
      TEV_glsq->Draw("fsame");
      LEP_glsq->Draw("fsame");
      NoSol_glsq->Draw("fsame"); NoSol_glsq_aux->Draw("l");

      TLegend* legexp = new TLegend(0.19,0.48,0.35,0.73,NULL,"brNDC");
      legexp->SetFillColor(0);legexp->SetShadowColor(0); //legexp->SetFillStyle(4000);
      legexp->SetTextFont(42);legexp->SetTextSize(0.025);legexp->SetBorderSize(0);
      legexp->AddEntry(TEV_glsq, "Tevatron RunI", "f");
      legexp->AddEntry(CDF_glsq, "CDF RunII","f");
      legexp->AddEntry(DEZ_glsq, "D0 RunII","f");
      legexp->AddEntry(LEP_glsq, "LEP2","f");
      legexp->Draw();

      TLatex noSolTex; noSolTex.SetTextSize(0.03); noSolTex.SetTextColor(1);
      noSolTex.DrawLatex( 800, 600, "no CMSSM solution");

      sprintf(plotKeyStr, "Exclusion_gl_sq_tb%d", tanbeta_input); plotKeyStrT = TString(plotKeyStr);
   }
   leg->Draw();
   gPad->RedrawAxis();
   c1->SaveAs("results/comp_Exp_"+plotKeyStrT+".pdf");
   c1->SaveAs("results/comp_Exp_"+plotKeyStrT+".png");
   c1->SaveAs("results/comp_Exp_"+plotKeyStrT+".C");

}

int main(int argc, char** argv)
{
  return plot(argc, argv);
}
