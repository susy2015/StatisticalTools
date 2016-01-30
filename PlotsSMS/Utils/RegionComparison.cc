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
#include "TLatex.h"
#include "TMarker.h"
#include "TPaveText.h"

#include <string>
#include <cmath>
#include <stdio.h>

int RegionComparison(int argc, char** argv)
{
   //interactive root session
   //TApplication theApp("App", 0, 0);
   if (gROOT->IsBatch()) {
      fprintf(stderr, "%s: cannot run in batch mode\n", argv[0]);
      return 1;
   }

   util::StyleSettings::paperNoTitle();
   gStyle->SetPadBottomMargin(0.18);
 
   TCanvas * c1 = new TCanvas("c1","c1",600,600);
   c1->cd();
   
   //Get limits from signal scan ---------------------------------------------------
   TheLimits * genpoints[5];
   genpoints[0] = new TheLimits();
   genpoints[1] = new TheLimits();
   genpoints[2] = new TheLimits();
   genpoints[3] = new TheLimits();
   genpoints[4] = new TheLimits();
   //genpoints->Fill("limits_HT500_MHT350/filelist.txt"); 
   //genpoints->Fill("mediumHTMHT_v3/filelist.txt"); 
   //genpoints->Fill("highHT_v1/filelist.txt"); 

   //genpoints->Fill("limits_HT500_MHT350_output_v1/filelist.txt"); 
   //genpoints->Fill("filelist.txt"); 
   //genpoints->Fill("limits_HT800_MHT500_output_v1/filelist.txt"); 
   
   //genpoints[0]->Fill("filelist.txt");
   genpoints[0]->Fill("mediumHTMHT_plots/filelist.txt","mediumHTMHT_plots/");
   genpoints[1]->Fill("highHT_plots/filelist.txt","highHT_plots/");
   genpoints[2]->Fill("highHTMHT_plots/filelist.txt","highHTMHT_plots/");
   genpoints[3]->Fill("HT1000MHT400_plots/filelist.txt","HT1000MHT400_plots/");
   genpoints[4]->Fill("HT1200MHT200_plots/filelist.txt","HT1200MHT200_plots/");

   genpoints[0]->FillEmptyPointsByInterpolationInM0M12();
   genpoints[1]->FillEmptyPointsByInterpolationInM0M12();
   genpoints[2]->FillEmptyPointsByInterpolationInM0M12();
   //genpoints[3]->FillEmptyPointsByInterpolationInM0M12();
   //genpoints[4]->FillEmptyPointsByInterpolationInM0M12();

   //Get SUSY masses and x-sections from generator scan ----------------------------
   //and match them to the signal scan
   genpoints[0]->FillGeneratorMasses("GenScan_tb10.dat");
   genpoints[0]->match();
   genpoints[1]->FillGeneratorMasses("GenScan_tb10.dat");
   genpoints[1]->match();
   genpoints[2]->FillGeneratorMasses("GenScan_tb10.dat");
   genpoints[2]->match();
   genpoints[3]->FillGeneratorMasses("GenScan_tb10.dat");
   genpoints[3]->match();
   genpoints[4]->FillGeneratorMasses("GenScan_tb10.dat");
   genpoints[4]->match();

   TString region[6];
   region[0] = "mediumHTMHT_";
   region[1] = "highHT_";
   region[2] = "highHTMHT_";
   region[3] = "HT1000MHT400_";
   region[4] = "HT1200MHT200_";
   region[5] = "combined_";

   TGraph * gCLsExpExcl_NLO[6];
   TGraph * gCLsObsExcl_NLO[6];
   //TGraph * gCLsObsExcl_LO[4];
   TGraph * gCLsExpExclm1[6];
   TGraph * gCLsExpExclp1[6];

   TH2* ObsExcl_NLO[6];
   TH2* ExpExcl_NLO[6];
   TH2* ExpExclm1[6];
   TH2* ExpExclp1[6];

   PlotTools<SusyScan> * plotTools[6];

   //Was not originally designed as a loop, so this will give you warnings about possible memory leaks.
   for(int i = 0; i < 6; i++){
     if(i < 5){
   //the plotting ----------------------------------------------------------------------
   //plotting helper functions
   plotTools[i] = new PlotTools<SusyScan>(genpoints[i]->GetScan());
   PlotTools<GeneratorMasses> * plotMasses = new PlotTools<GeneratorMasses>(genpoints[i]->GetGeneratorMasses());

   //iso mass lines
   //usage: Line( <x-var>, <y-var>, <mass-function>, <mass-value>, <matching-width> )
   TGraph * gl500 = plotMasses->Line(Mzero, Mhalf, MGluino, 500, 0.8 );
   TGraph * sq500 = plotMasses->Line(Mzero, Mhalf, MSquarkL, 500, 1);
   TGraph * chi100 = plotMasses->Line(Mzero, Mhalf, MChi1, 50, 20); 

   // Histograms for combination plot later ==================================================================
   ObsExcl_NLO[i] = new TH2F("obs",";m_{0} [GeV]; m_{1/2} [GeV]; 95% CL Expected Exclusion",
			    100,0,2009.9,38,0,770);
   plotTools[i]->Area(ObsExcl_NLO[i], Mzero, Mhalf,  NLOObsExclCL);
   ExpExcl_NLO[i] = new TH2F("exp",";m_{0} [GeV]; m_{1/2} [GeV]; 95% CL Expected Exclusion",
			     100,0,2009.9,38,0,770);
   plotTools[i]->Area(ExpExcl_NLO[i], Mzero, Mhalf,  NLOExpExclCL);
   ExpExclm1[i] = new TH2F("expm1",";m_{0} [GeV]; m_{1/2} [GeV]; 95% CL Expected Exclusion",
			     100,0,2009.9,38,0,770);
   plotTools[i]->Area(ExpExclm1[i], Mzero, Mhalf,  NLOExpExclCLm1sigma);
   ExpExclp1[i] = new TH2F("expp1",";m_{0} [GeV]; m_{1/2} [GeV]; 95% CL Expected Exclusion",
			     100,0,2009.9,38,0,770);
   plotTools[i]->Area(ExpExclp1[i], Mzero, Mhalf,  NLOExpExclCLp1sigma);

   // Exclusion in M0 - M1/2 =================================================================================
   // hs defined range in which the contrours are calculated, the binning matters! binning should be equal to signal-scan 
   TH2F*hs = new TH2F("hs","",95,100,2009.9,38,0,770);
   // usage: GetContour( <range-hist>, <x-var>, <y-var>, <limit-function>, <contour-levels>, <contour index>, <color>, <style> )
   gCLsExpExcl_NLO[i]     = plotTools[i]->GetContour(hs,Mzero,Mhalf,NLOExpExclCL, 3,0, 1,2); 
   gCLsObsExcl_NLO[i]     = plotTools[i]->GetContour(hs,Mzero,Mhalf,NLOObsExclCL, 3,0, 1,1);
   //gCLsObsExcl_LO[i]      = plotTools->GetContour(hs,Mzero,Mhalf,ObsExclCL,    3,0, 2,1);
   gCLsExpExclm1[i]       = plotTools[i]->GetContour(hs,Mzero,Mhalf,NLOExpExclCLm1sigma,3,0, 5,2); 
   gCLsExpExclp1[i]       = plotTools[i]->GetContour(hs,Mzero,Mhalf,NLOExpExclCLp1sigma,3,0, 5,2); 
     }

     // hexcl is drawn to define axis range and labels:
     TH2F*hexcl = new TH2F("hexcl",";m_{0} (GeV); m_{1/2} (GeV); 95% CL Exclusion", 100,0,1509.9,45,110,700);
     hexcl->SetTitle("CMS Preliminary");


   gCLsObsExcl_NLO[i]->SetLineWidth(2);
   gCLsExpExcl_NLO[i]->SetLineWidth(2);
   
   //smooth contours (2D Gaussian smoothing)
   Smooth( gCLsExpExcl_NLO[i], 27 );
   //Smooth( gCLsObsExcl_LO, 27 );
   Smooth( gCLsObsExcl_NLO[i], 27 );
   Smooth( gCLsExpExclm1[i], 27 );
   Smooth( gCLsExpExclp1[i], 27 );
   //"Combine" HT and  limit by choosing the better one
   TGraph * gCLsExp1Sigma = MakeBand(gCLsExpExclp1[i], gCLsExpExclm1[i]);gCLsExp1Sigma->SetFillStyle(4010);
   //gCLsExp1Sigma->Print("all");
   //draw old exclusion Limits
   hexcl->GetYaxis()->SetTitleOffset(1.3);
   hexcl->GetXaxis()->SetTitleOffset(0.92);
   hexcl->Draw("colz");
   TGraphErrors* RA1 = RA1_NLO();
   RA1->SetLineColor(kRed+1);
   RA1->SetLineStyle(1);
   RA1->SetLineWidth(2);
   TGraph* LEP_ch = set_lep_ch(10);
   TGraph* LEP_sl = set_lep_sl(10);//slepton curve
   TGraph* TEV_sg_cdf = set_tev_sg_cdf(10);//squark gluino cdf
   TGraph* TEV_sg_d0 = set_tev_sg_d0(10);//squark gluino d0
   //TGraph* TEV_tlp_cdf = set_tev_tlp_cdf(10);//trilepton cdf
   //TGraph* TEV_tlp_d0 = set_tev_tlp_d0(10);//trilepton d0
   TGraph* Atlas = Atlas_m0_m12_tb3_obs();
   Atlas->SetLineColor( 28 );
   TGraph* stau = set_tev_stau(10);//stau 
   TGraph* TEV_sn_d0_1 = set_sneutrino_d0_1(10);
   TGraph* TEV_sn_d0_2 = set_sneutrino_d0_2(10);
   TGraphErrors* First  = getObserved_NLO_tanBeta10();
   TGraphErrors* Second = getExpected_NLO_tanBeta10();//getLO_jetMultis();
   TGraphErrors* Third  = getLO_tanBeta10();
   First->GetXaxis()->SetRangeUser(0,505);
   First->GetYaxis()->SetRangeUser(80,500);
   First->GetXaxis()->SetTitle("m_{0} (GeV)");
   First->GetYaxis()->SetTitle("m_{1/2} (GeV)");
   TSpline3 *sFirst = new TSpline3("sFirst",First); sFirst->SetLineColor(kBlue);sFirst->SetLineWidth(1);
   TSpline3 *sSecond = new TSpline3("sSecond",Second);
   sSecond->SetLineColor(kBlue); sSecond->SetLineStyle(2); sSecond->SetLineWidth(1);
   TSpline3 *sThird = new TSpline3("sThird",Third);
   sThird->SetLineColor(kBlue); sThird->SetLineStyle(4); sThird->SetLineWidth(1);
   //TEV_sn_d0_1->Draw("fsame");
   //TEV_sn_d0_2->Draw("fsame"); //only for tb=3
   //gCLsExp1Sigma->Draw("F");
   LEP_ch->Draw("fsame");
   LEP_sl->Draw("fsame");
   TEV_sg_cdf->Draw("fsame");
   TEV_sg_cdf->Draw("lsame");
   TEV_sg_d0->Draw("fsame");
   TEV_sg_d0->Draw("lsame");
   TLatex b; b.SetTextSize(0.02); b.SetTextColor(1); b.SetTextAngle(+82);
   //sFirst->Draw("same");
   //sSecond->Draw("same");
   //sThird->Draw("same");
   //Atlas->Draw("c,same");
   TLatex ms; ms.SetTextSize(0.025); ms.SetTextFont(42);//ms.SetTextColor(12);
   ms.DrawLatex(450,540,"tan#beta=10, #mu>0, A_{0}=0"); 
   //draw LM points
   TMarker* LM0 = new TMarker(200.,160.,20);
   TMarker* LM1 = new TMarker(60.,250.,20);
   LM0->SetMarkerSize(0.7);
   LM1->SetMarkerSize(0.7);
   TLatex* tLM0 = new TLatex(206.,162.,"LM0");
   tLM0->SetTextSize(0.02);
   TLatex* tLM1 = new TLatex(66.,256.,"LM1");
   tLM1->SetTextSize(0.02);
   //LM0->Draw("same");   
   //tLM0->Draw("same");
   //LM1->Draw("same");   
   //tLM1->Draw("same");
   TLegend* legexp = new TLegend(0.61,0.61,0.93,0.86,NULL,"brNDC");    legexp->SetFillColor(0);legexp->SetShadowColor(0);legexp->SetFillStyle(4000);legexp->SetTextFont(42);legexp->SetTextSize(0.025);legexp->SetBorderSize(0);
   //TEV_sg_cdf.SetLineColor(1);  
   //legexp->SetHeader("CMS");
   legexp->AddEntry(TEV_sg_cdf,"CDF  #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=5, #mu<0}","f"); 
   legexp->AddEntry(TEV_sg_d0,"D0   #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=3, #mu<0}","f");  
   //ch_gr.SetLineColor(1); 
   legexp->AddEntry(LEP_ch,"LEP2   #tilde{#chi}_{1}^{#pm}","f");   
   //sl_gr.SetLineColor(1); 
   legexp->AddEntry(LEP_sl,"LEP2   #tilde{#font[12]{l}}^{#pm}","f"); //NOT FOR tb=50!
   //if(tanbeta == 3) 
   //legexp->AddEntry(TEV_sn_d0_1,"D0  #chi^{#pm}_{1}, #chi^{0}_{2}","f");  
   //legexp->AddEntry(sFirst, "CMS #alpha_{T}");
   //legexp->AddEntry(Atlas,  "Atlas, #scale[0.8]{tan#beta=3}","l");
   //legexp->Draw();
   //RA1 limit contour:
   RA1->Draw("C");
   //draw RA2 limit contours:
   //gCLsObsExcl_NLO[i]->Draw("l");
   //CLsObsLO->Draw("l");
   // Error band needs to move to top so it doesn't cover other things.
   //gCLsExp1Sigma->Draw("F");
   //gCLsExpExcl_NLO[i]->Draw("l");
   //gCLsObsExclMHT_NLO->Draw("l");
   stau->Draw("fsame");
   b.DrawLatex( 80,570,"#tilde{#tau} LSP"); 
   //constant ssqquark and gluino lines
   TF1* lnsq[4];
   TF1* lngl[4];
   TLatex sqt; sqt.SetTextSize(0.02); sqt.SetTextAngle(-28);sqt.SetTextColor(kGray+2);
   sqt.DrawLatex(148,230,"#font[92]{#tilde{q}(500)GeV}");
   sqt.SetTextAngle(-22);
   sqt.DrawLatex(148,360,"#font[92]{#tilde{q}(750)GeV}");
   sqt.SetTextAngle(-16);
   sqt.DrawLatex(148,485,"#font[92]{#tilde{q}(1000)GeV}");
   TLatex glt; glt.SetTextSize(0.02); sqt.SetTextAngle(-4); glt.SetTextColor(kGray+2);
   glt.DrawLatex(1200,170,"#font[92]{#tilde{g}(500)GeV}");
   glt.DrawLatex(1200,275,"#font[92]{#tilde{g}(750)GeV}");
   glt.DrawLatex(1200,385,"#font[92]{#tilde{g}(1000)GeV}");
   int tanBeta_=10;
   for(int j = 0; j < 4; j++){
    lnsq[j] = constant_squark(tanBeta_,j);
    lngl[j] = constant_gluino(tanBeta_,j);
    lngl[j]->Draw("same");   
    lnsq[j]->Draw("same");
   }
   TLegend * leg = new TLegend(0.25,0.74,0.6,0.86);
   leg->SetBorderSize(0);leg->SetFillColor(0);leg->SetFillStyle(4000);leg->SetTextFont(42);legexp->SetTextSize(0.025);
   TGraph * expLeg = (TGraph*)gCLsExpExcl_NLO[i]->Clone();expLeg->SetFillStyle(gCLsExp1Sigma->GetFillStyle());expLeg->SetFillColor(gCLsExp1Sigma->GetFillColor());
   leg->SetHeader("L_{int} = 1.1 fb^{-1}, #sqrt{s} = 7 TeV");
   leg->AddEntry(gCLsObsExcl_NLO[i],"Observed (NLO), CL_{s}","l");
   //leg->AddEntry(CLsObsLO, "Observed, LO","l");
   //leg->AddEntry(gCLsExpExcl_NLO, "Expected, NLO","l");
   leg->AddEntry(expLeg,   "Expected #pm 1#sigma (NLO), CL_{s}","lf");
   leg->AddEntry(RA1, "RA1 2011 (NLO), CL_{s}","l");
   //leg->AddEntry(CLsExpNoSNLO,   "Expected, no-signal hyp., NLO","l");
   //leg->Draw();
   gPad->RedrawAxis();
   //if( i == 5 ){c1->SaveAs("results/ExclusionComparison_m0_m12_tb10.pdf");}
   // ===============================================================================================================


   //theApp.Run();
   
   if( i == 4 ){
   //combined exclusion for all five seach regions.
   ObsExcl_NLO[5] = plotTools[i]->BinWiseOr(ObsExcl_NLO[0],ObsExcl_NLO[1]);
   ObsExcl_NLO[5] = plotTools[i]->BinWiseOr(ObsExcl_NLO[5],ObsExcl_NLO[2]);
   ObsExcl_NLO[5] = plotTools[i]->BinWiseOr(ObsExcl_NLO[5],ObsExcl_NLO[3]);   
   ObsExcl_NLO[5] = plotTools[i]->BinWiseOr(ObsExcl_NLO[5],ObsExcl_NLO[4]);
   gCLsObsExcl_NLO[5]     = plotTools[i]->GetContour(ObsExcl_NLO[5],3,0);
   gCLsObsExcl_NLO[5]->SetLineColor(1);
   gCLsObsExcl_NLO[5]->SetLineStyle(1);

   ExpExcl_NLO[5] = plotTools[i]->BinWiseOr(ExpExcl_NLO[0],ExpExcl_NLO[1]);
   ExpExcl_NLO[5] = plotTools[i]->BinWiseOr(ExpExcl_NLO[5],ExpExcl_NLO[2]);
   ExpExcl_NLO[5] = plotTools[i]->BinWiseOr(ExpExcl_NLO[5],ExpExcl_NLO[3]);
   ExpExcl_NLO[5] = plotTools[i]->BinWiseOr(ExpExcl_NLO[5],ExpExcl_NLO[4]);
   gCLsExpExcl_NLO[5]     = plotTools[i]->GetContour(ExpExcl_NLO[5],3,0);
   gCLsExpExcl_NLO[5]->SetLineColor(1);
   gCLsExpExcl_NLO[5]->SetLineStyle(2);

   ExpExclm1[5] = plotTools[i]->BinWiseOr(ExpExclm1[0],ExpExclm1[1]);
   ExpExclm1[5] = plotTools[i]->BinWiseOr(ExpExclm1[5],ExpExclm1[2]);
   ExpExclm1[5] = plotTools[i]->BinWiseOr(ExpExclm1[5],ExpExclm1[3]);
   ExpExclm1[5] = plotTools[i]->BinWiseOr(ExpExclm1[5],ExpExclm1[4]);
   gCLsExpExclm1[5]       = plotTools[i]->GetContour(ExpExclm1[5],3,0);
   gCLsExpExclm1[5]->SetLineColor(5);
   gCLsExpExclm1[5]->SetLineStyle(2);

   ExpExclp1[5] = plotTools[i]->BinWiseOr(ExpExclp1[0],ExpExclp1[1]);
   ExpExclp1[5] = plotTools[i]->BinWiseOr(ExpExclp1[5],ExpExclp1[2]);
   ExpExclp1[5] = plotTools[i]->BinWiseOr(ExpExclp1[5],ExpExclp1[3]);
   ExpExclp1[5] = plotTools[i]->BinWiseOr(ExpExclp1[5],ExpExclp1[4]);
   gCLsExpExclp1[5]       = plotTools[i]->GetContour(ExpExclp1[5],3,0);
   gCLsExpExclp1[5]->SetLineColor(5);
   gCLsExpExclp1[5]->SetLineStyle(2);
   }
   
   }

   for( int i = 0; i < 6 ; i++ ){
   gCLsExpExcl_NLO[i]     = plotTools[i]->GetContour(hs,Mzero,Mhalf,NLOExpExclCL, 3,0, 1,2);
   gCLsObsExcl_NLO[i]     = plotTools[i]->GetContour(hs,Mzero,Mhalf,NLOObsExclCL, 3,0, 1,1);

   gCLsObsExcl_NLO[i]->Draw("l");
   }

   c1->SaveAs("results/ExclusionComparison.pdf");

}


int main(int argc, char** argv)
{
  return RegionComparison(argc, argv);
}
