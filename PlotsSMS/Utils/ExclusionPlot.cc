#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>

#include "ExclusionPlot.hh"
 
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "vector.h"
#include "TMath.h"
#include "TList.h"
#include "TGraph.h"
#include "TObjArray.h"

void ExclusionPlot(){
  gStyle->SetPalette(1);

  //get yield plot
  //  hev = yieldPlot(mSuGraFile,mSuGraDir, mSuGraHist);
  //setPlottingStyle(*hev);

  // setPlottingStyle(*exclusionPlot,lineStyle);

  Int_t tanBeta = 3;
   
  /* TH1F* First = getHisto_1d("./","ExclusionLimit","Significance_NLO_expected_tanBeta50.root");
  setPlottingStyle(*First);
  First->SetLineStyle(1);
  First->SetLineWidth(1);

  TH1F* Second = getHisto_1d("./","ExclusionLimit","Significance_NLO_observed_tanBeta50.root"); 
  Second->SetLineStyle(2);
  Second->SetLineWidth(1);
  Second->SetLineColor(kRed);
 
  
   TH1F* Third = getHisto_1d("./","ExclusionLimit","Significance_LO_observed_tanBeta50.root");
  Third->SetLineStyle(1);
  Third->SetLineWidth(3);
  Third->SetLineColor(kGreen+2);
 
  
  exclusionPlots.push_back(First);
    exclusionPlots.push_back(Second);
    exclusionPlots.push_back(Third);*/
  
  
  CommandMSUGRA("35pb_expected_11.root",tanBeta);


}



void CommandMSUGRA(TString plotName_,Int_t tanBeta_, Int_t sigma){
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1); 
  gStyle->SetTextFont(42);
  
  //convert tanb value to string
  std::stringstream tmp;
  tmp << tanBeta_;
  TString tanb( tmp.str() );
  
  
  // Output file
  cout << " create " << plotName_ << endl;
  TFile* output = new TFile( plotName_, "RECREATE" );
  if ( !output || output->IsZombie() ) { std::cout << " zombie alarm output is a zombie " << std::endl; }
  

  //set old exclusion Limits
  TGraph* LEP_ch = set_lep_ch(tanBeta_);
  TGraph* LEP_sl = set_lep_sl(tanBeta_);//slepton curve
  TGraph* TEV_sg_cdf = set_tev_sg_cdf(tanBeta_);//squark gluino cdf
  TGraph* TEV_sg_d0 = set_tev_sg_d0(tanBeta_);//squark gluino d0
  //  TGraph* TEV_tlp_cdf = set_tev_tlp_cdf(tanBeta_);//trilepton cdf
  //  TGraph* TEV_tlp_d0 = set_tev_tlp_d0(tanBeta_);//trilepton d0
  TGraph* stau = set_tev_stau(tanBeta_);//stau 

  TGraph* TEV_sn_d0_1 = set_sneutrino_d0_1(tanBeta_);
  TGraph* TEV_sn_d0_2 = set_sneutrino_d0_2(tanBeta_);

  //constant ssqquark and gluino lines
  TF1* lnsq[4];
  TF1* lngl[4];
  
  TLatex* sq_text[4];
  TLatex* gl_text[4];

  for(int i = 0; i < 4; i++){
    lnsq[i] = constant_squark(tanBeta_,i);
    sq_text[i] = constant_squark_text(i,*lnsq[i],tanBeta_);
    lngl[i] = constant_gluino(tanBeta_,i);
    gl_text[i] = constant_gluino_text(i,*lngl[i]);
  }
  
  //Legends
  TLegend* legst = makeStauLegend(0.05,tanBeta_);
  TLegend* legexp = makeExpLegend( *TEV_sg_cdf,*TEV_sg_d0,*LEP_ch,*LEP_sl,*TEV_sn_d0_1,0.03,tanBeta_);
  
 
  //make Canvas
  TCanvas* cvsSys = new TCanvas("cvsnm","cvsnm",0,0,800,600);
  gStyle->SetOptTitle(0);
  cvsSys->SetFillColor(0);
  cvsSys->GetPad(0)->SetRightMargin(0.07);
  cvsSys->Range(-120.5298,26.16437,736.0927,500);
  //  cvsSys->Range(-50.5298,26.16437,736.0927,500);
  cvsSys->SetFillColor(0);
  cvsSys->SetBorderMode(0);
  cvsSys->GetPad(0)->SetBorderSize(2);
  cvsSys->GetPad(0)->SetLeftMargin(0.1407035);
  cvsSys->GetPad(0)->SetTopMargin(0.08);
  cvsSys->GetPad(0)->SetBottomMargin(0.13);

  cvsSys->SetTitle("tan#beta="+tanb);
 
  output->cd();
  
  //and now
  //the exclusion limits
  TGraphErrors* First ;
  TGraphErrors* Second;
  TGraphErrors* Third;
  TGraphErrors* Second_up;
  TGraphErrors* Second_low;

  if(tanBeta_ == 3){
    First = getObserved_NLO_tanBeta3_2();
    Second = getExpected_NLO_tanBeta3_2();
    Second_up = getExpected_NLO_tanBeta3_up();
    Second_low = getExpected_NLO_tanBeta3_low();
    Third = getLO_tanBeta3();
  }
  if(tanBeta_ == 10){
    First = getObserved_NLO_tanBeta10_2();
    Second = getExpected_NLO_tanBeta10_2();
    Second_up = getExpected_NLO_tanBeta10_up();
    Second_low = getExpected_NLO_tanBeta10_low();
    Third = getLO_tanBeta10();
  }
  if(tanBeta_ == 50){
    First = getObserved_NLO_tanBeta50_2();
    Second = getExpected_NLO_tanBeta50_2();
    Second_up = getExpected_NLO_tanBeta50_up();
    Second_low = getExpected_NLO_tanBeta50_low();
    Third = getLO_tanBeta50();
  }

//   First->SetMarkerColor(kWhite);
//   First->GetXaxis()->SetRangeUser(2.,500.);
//   First->GetYaxis()->SetRangeUser(80,500);
//   if(tanBeta_ == 50) First->GetXaxis()->SetRangeUser(200,500);
//   First->GetXaxis()->SetTitle("m_{0} (GeV)");
//   First->GetYaxis()->SetTitle("m_{1/2} (GeV)");
//   First->GetYaxis()->SetTitleOffset(0.8);

  double m0min = 0;
  if (tanBeta_ == 50) m0min=200;

  //  TH2D* hist = new TH2D("h","h",100,m0min,500,100,80,500);
 
  //for efficiency plot
  //  TH2F* hist = getHisto("./","m0_m12_0","TanBeta3_efficiency.root");
  //CHF DEFAULT RESULT
  //  TH2F* hist = getHisto("/home/hep/elaird1/public_html/57_stat_plots/05_cms_sm_fc/","ExclusionLimit","profileLikelihood_tanBeta"+tanb+"_nlo_1HtBin_expR.root");

   TH2F* hist = getHisto("/home/hep/elaird1/public_html/57_stat_plots/05_cmssm_fc/","ExclusionLimit","feldmanCousins_tanBeta"+tanb+"_nlo_1HtBin_expR.root");

  
  if (sigma) {
    //orig histos
    //    TH2F* exp = getHisto("./TedsHistos/","yield_signal_median","SigmaBand_tanBeta"+tanb+"_2HTBins_expR.root");
    //    TH2F* up = getHisto("./TedsHistos/","yield_signal_upper","SigmaBand_tanBeta"+tanb+"_2HTBins_expR.root");
    //    TH2F* down = getHisto("./TedsHistos/","yield_signal_lower","SigmaBand_tanBeta"+tanb+"_2HTBins_expR.root");

    //latest sigma bands for new result, 2 HT bins, Rconst
    TH2F* exp = getHisto("./TedsHistos/","dsMedian","profileLikelihood_tanBeta"+tanb+"_nlo_1HtBin_expR_computeExpectedLimit_results.root");
    TH2F* up = getHisto("./TedsHistos/","dsMedianPlusOneSigma","profileLikelihood_tanBeta"+tanb+"_nlo_1HtBin_expR_computeExpectedLimit_results.root");
    TH2F* down = getHisto("./TedsHistos/","dsMedianMinusOneSigma","profileLikelihood_tanBeta"+tanb+"_nlo_1HtBin_expR_computeExpectedLimit_results.root");

  


    //    TH2F* exp = getHisto("./TedsHistos/","yield_signal_median","SigmaBand_tanBeta3_1HTBin_exponentialR.root");
    //   TH2F* up = getHisto("./TedsHistos/","yield_signal_upper","SigmaBand_tanBeta3_1HTBin_exponentialR.root");
    // TH2F* down = getHisto("./TedsHistos/","yield_signal_lower","SigmaBand_tanBeta3_1HTBin_exponentialR.root");


    //CHF hack for limit evolution 
    //HT1 exp
    //    TH2F* hist = getHisto("./TedsHistos/story/","ExclusionLimit","profileLikelihood_tanBeta3_nlo_1HtBin_expR.root");
    //  TH2F* hist = getHisto("./TedsHistos/paper/","ExclusionLimit","Significance_tanBeta3_NLO_observed.root");
    //HT1 exp RA2BG
    //    TH2F* exp = getHisto("./TedsHistos/story/","ExclusionLimit","profileLikelihood_tanBeta3_nlo_1HtBin_expR_Ra2SyncHack.root");
  //HT1 const
    //    TH2F* up = getHisto("./TedsHistos/story/","ExclusionLimit","profileLikelihood_tanBeta3_nlo_1HtBin_constantR.root");
  //HT2 exp
    //    TH2F* down = getHisto("./TedsHistos/story/","ExclusionLimit","profileLikelihood_tanBeta3_nlo_1HtBin_expR.root");

}

  gStyle->SetPalette(1);

  hist->GetXaxis()->SetRangeUser(0,500);
  if (tanBeta_ == 10)  hist->GetXaxis()->SetRangeUser(0,600);
  if (tanBeta_ == 50)  hist->GetXaxis()->SetRangeUser(200,600);
  hist->GetYaxis()->SetRangeUser(0,500);
  // hist->Draw("colz");  
  
  hist->SetLineWidth(3);
  hist->SetLineColor(9);
  hist->SetLineStyle(2);
  //  hist->SetContourLevel(1,0.05);
  //  hist->SetContourLevel(7,0.32);

  TH2F* clone = fillHoles(hist);
  
  if (sigma) {
    fillHoles(exp);
    fillHoles(up);
    fillHoles(down);
  }
  //  TH2F* clone = (TH2F*)hist->Clone("clone");

  clone->GetXaxis()->SetTitle("");
  clone->GetYaxis()->SetTitle("");

  output.cd();
  //  cvsSys->Update();
  //  hist->Reset();
  //  hist->Draw();

  //CHF change contour level for tanb3!!!

  TGraph* myGraph = getGraph( clone , -0.5 );
  TGraph* smg = (TGraph*)myGraph->Clone("smg");
  //  smg->Print("all");
  if (sigma) {
    TGraph* myexp = getGraph( exp , -0.5);
    TGraph* smexp = (TGraph*)myexp->Clone("smexp");
    // smexp->Print("all");
    TGraph* myup = getGraph( up , -0.5 );
    TGraph* smup = (TGraph*)myup->Clone("smup");
    //smup->Print("all");
    TGraph* mydown = getGraph( down , -0.5 );
    TGraph* smdown = (TGraph*)mydown->Clone("smdown");
    // smdown->Print("all");
  }

  output.cd();
  hist->Reset();
  hist->Draw();

  cvsSys->Update();

  hist->GetXaxis()->SetTitle("m_{0} (GeV)");
  hist->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hist->GetYaxis()->SetTitleOffset(1.);
  hist->GetXaxis()->SetNdivisions(506);
  //  if (tanBeta_ == 50)  hist->GetXaxis()->SetNdivisions(504);
  hist->GetYaxis()->SetNdivisions(506);



  int col[]={2,3,4};

  TSpline3 *sFirst = new TSpline3("sFirst",First);
  sFirst->SetLineColor(kRed);
  sFirst->SetLineWidth(3);

  TSpline3 *sSecond = new TSpline3("sSecond",Second);
  sSecond->SetLineColor(kBlue+3);
  sSecond->SetLineStyle(4);
  sSecond->SetLineWidth(3);

  TSpline3 *sSecond_up = new TSpline3("sSecond_up",Second_up);
  sSecond_up->SetLineColor(kCyan);
  sSecond_up->SetLineStyle(1);
  sSecond_up->SetLineWidth(3);

  TSpline3 *sSecond_low = new TSpline3("sSecond_low",Second_low);
  sSecond_low->SetLineColor(kCyan);
  sSecond_low->SetLineStyle(1);
  sSecond_low->SetLineWidth(3);
  
  TSpline3 *sThird = new TSpline3("sThird",Third);
  sThird->SetLineColor(kGreen+2);
  sThird->SetLineStyle(4);
  sThird->SetLineWidth(3);


    sSecond_up->Draw("h same");
   sSecond_low->Draw("h same");
      
 
  //constant squark and gluino mass contours
  for (int it=1;it<4;it++) {   
    lngl[it]->Draw("same");   
    lnsq[it]->Draw("same");
    sq_text[it]->Draw();
    gl_text[it]->Draw();
  }
   
    
  //exclusion limits previous experiments
  if(tanBeta_ == 3){
    TEV_sn_d0_1->Draw("fsame");
    TEV_sn_d0_2->Draw("fsame");
  }
  LEP_ch->Draw("fsame");
  if (tanBeta_ != 50) LEP_sl->Draw("fsame");

  TEV_sg_cdf->Draw("fsame");
  TEV_sg_d0->Draw("same");  
  TEV_sg_d0->Draw("fsame");


  //other labels
  Double_t xpos = 0;
  Double_t xposi = 0;
  Double_t ypos = 0;
  if(tanBeta_ == 50) xposi = 100;
  if(tanBeta_ == 50) xpos = 200;
  if(tanBeta_ == 50) ypos = -10;
  
  TLatex* lumilabel = new TLatex(105.+xposi,510.,"CMS preliminary    L_{int} = 35 pb^{-1}, #sqrt{s} = 7 TeV");
  lumilabel->SetTextSize(0.05);

  lumilabel->Draw("same");

  TString text_tanBeta;
  text_tanBeta =  "tan#beta = "+tanb+", A_{0} = 0, #mu > 0";
  TLatex* cmssmpars = new TLatex(70.+xpos,340.+ypos,text_tanBeta);
  cmssmpars->SetTextSize(0.04);

  cmssmpars->Draw("same");
 
  //LM points
  TMarker* LM0 = new TMarker(200.,160.,20);
  TMarker* LM1 = new TMarker(60.,250.,20);
    
  LM0->SetMarkerSize(1.2);
  LM1->SetMarkerSize(1.2);
    
  TLatex* tLM0 = new TLatex(205.,160.,"LM0");
  tLM0->SetTextSize(0.035);
    
  TLatex* tLM1 = new TLatex(65.,243.,"LM1");
  tLM1->SetTextSize(0.035);
  
  if (tanBeta_ != 50){
    LM0->Draw("same");   
    tLM0->Draw("same");
    LM1->Draw("same");   
    tLM1->Draw("same");
  }
  


  //expected and observed (LO & NLO) contours
  
  TLegend* myleg = new TLegend(0.25,0.75,0.5,0.9,NULL,"brNDC");
  myleg->SetFillColor(0); 
  myleg->SetShadowColor(0);
  myleg->SetTextSize(0.03);
  myleg->SetBorderSize(0);
  
  //  myleg->AddEntry(sSecond,"NLO Expected Limit","L");
  //  myleg->AddEntry(sFirst,"NLO Observed Limit","L"); 
  //  myleg->AddEntry(sThird,"LO Observed Limit","L");
  
  myleg->AddEntry(sFirst,"Observed Limit, NLO","L");
  //  myleg->AddEntry(smg,"Obs. Limit, NLO, 1 HT bins, R_{exp}","L");
  //  myleg->AddEntry(sFirst,"Obs. Limit, NLO, 1 HT bin, R_{exp}","L"); 
  if (sigma) {
    myleg->AddEntry(sSecond,"Median Expected Limit","L");
    //    myleg->AddEntry(smexp,"Median Exp. Limit, NLO, 2 HT bins, R_{const}","L");
    myleg->AddEntry(sSecond_up,"Expected Limit #pm 1 #sigma","F");
    //    myleg->AddEntry(smexp,"Obs. Limit NLO, 1HT bin, R_{exp}, RA2 EWK","L");
    //myleg->AddEntry(smup,"Obs. Limit NLO, 1 HT bin, R_{const}","L");
    //myleg->AddEntry(smdown,"Obs. Limit NLO, 2 HT bins, R_{exp}","L");
  }
  //  myleg->AddEntry(sSecond,"Old Exp. Limit #pm 1 #sigma, NLO, 1 HT bin, R_{exp}","L");

  sSecond_up->SetFillStyle(4010);
  sSecond_up->SetFillColor(kCyan-10);

  // sSecond_low->Draw("Csame");

  TGraph* smSec = (TGraph*)sSecond_low->Clone("smSecondlow");
  smSec->SetFillColor(0);

  smSec->Draw("Csame");

  sSecond_low->SetFillStyle(1001);
  sSecond_low->SetFillColor(10);

   sFirst->Draw("same");    
   sSecond->Draw("same");  

  //sThird->Draw("same");

  smg->SetLineStyle(1);
  smg->SetLineColor(kRed);
  // smg->Draw("same");
  if (sigma) {
    smexp->SetLineColor(9);
    smexp->SetLineWidth(3);
    smexp->SetLineStyle(2);
    //smexp->Draw("same");
    smup->SetLineWidth(3);
    smup->SetFillStyle(3001);
    smup->SetFillColor(kRed);
    smup->SetLineStyle(3);
    smup->SetLineColor(9);
    //smup->Draw("same");
    smdown->SetLineWidth(3);
    smdown->SetLineStyle(3);
    smdown->SetLineStyle(3);
    smdown->SetLineColor(9);
    //smdown->Draw("same");
  }


  //stau=LSP contour
  stau->Draw("fsame");
  
  //legends
  legexp->Draw();
  legst->Draw();
  myleg->Draw();

  hist->Draw("sameaxis");
  cvsSys->RedrawAxis();

  cvsSys->Update();

  cvsSys->Write();
  
  
  cvsSys->SaveAs("RA1_ExclusionLimit_tanb"+tanb+".pdf");
  cvsSys->SaveAs("RA1_ExclusionLimit_tanb"+tanb+".png");
  
  
  output->Write();
  output->Close();
  delete output; 
  
}


void setPlottingStyle(TH1F& hsig){
  
  hsig.SetStats(kFALSE);
  
  hsig.SetAxisRange(80,500,"Y");
  hsig.SetAxisRange(0,520,"X");
  hsig.SetAxisRange(200,520,"X");

  hsig.GetXaxis()->SetTitle("m_{0} (GeV)");
  hsig.GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hsig.GetYaxis()->SetTitleOffset(0.8);
  hsig.GetYaxis()->SetTitleSize(0.06);
  hsig.GetYaxis()->SetLabelSize(0.06);
  hsig.GetXaxis()->SetTitleOffset(0.9);
  hsig.GetXaxis()->SetTitleSize(0.06);
  hsig.GetXaxis()->SetLabelSize(0.06);

  hsig.SetLineWidth(1);  
  hsig.SetLineColor(kBlue);  
  
}




TGraph* set_sneutrino_d0_1(Int_t tanBeta){
  double sn_m0[14]= {0,  0, 48, 55, 80, 90,100,105,109,105,100, 72, 55,0};
  double sn_m12[14]={0,140,210,220,237,241,242,241,230,220,210,170,150,0};

  TGraph* sn_d0_gr = new TGraph(14,sn_m0,sn_m12);

  sn_d0_gr->SetFillColor(kGreen+3);
  sn_d0_gr->SetFillStyle(3001);

  return sn_d0_gr;
}

TGraph* set_sneutrino_d0_2(Int_t tanBeta){
  double sn_m0[9]= {0, 45, 75,115,130,150,163,185,0};
  double sn_m12[9]={0,140,170,213,202,183,168,140,0};

  TGraph* sn_d0_gr_2 = new TGraph(9,sn_m0,sn_m12);

  sn_d0_gr_2->SetFillColor(kGreen+3);
  sn_d0_gr_2->SetFillStyle(3001);

  return sn_d0_gr_2;
}

TGraph* set_lep_ch(Int_t tanBeta){
  if(tanBeta == 3) return set_lep_ch_tanBeta3();
  if(tanBeta == 10) return set_lep_ch_tanBeta10();
  if(tanBeta == 50) return set_lep_ch_tanBeta50();
}

TGraph* set_lep_ch_tanBeta10(){

  double ch_m0[11];
  double ch_m12[11];

  ch_m0[0] = 0;
  ch_m0[1] = 100;
  ch_m0[2] = 200;
  ch_m0[3] = 300;
  ch_m0[4] = 400;
  ch_m0[5] = 500;
  ch_m0[6] = 600;
  ch_m0[7] = 700;
  ch_m0[8] = 800; 
  ch_m0[9] = 800;
  ch_m0[10] = 0;

  ch_m12[0] = 163;
  ch_m12[1] = 162;
  ch_m12[2] = 161;
  ch_m12[3] = 160;
  ch_m12[4] = 159;
  ch_m12[5] = 158;
  ch_m12[6] = 157;
  ch_m12[7] = 156;
  ch_m12[8] = 155.4;
  ch_m12[9] = 0;
  ch_m12[10] = 0;
  
  
  TGraph* ch_gr = new TGraph(11,ch_m0,ch_m12);

  ch_gr->SetFillColor(3);
  ch_gr->SetLineColor(3);
  //  ch_gr->SetLineWidth(3);
  ch_gr->SetFillStyle(3001);

  return ch_gr;

}



TGraph* set_lep_ch_tanBeta3(){

  double ch_m0[17];
  double ch_m12[17];

  ch_m0[0] = 0;
  ch_m0[1] = 100;
  ch_m0[2] = 150;
  ch_m0[3] = 200;
  ch_m0[4] = 250;
  ch_m0[5] = 300;
  ch_m0[6] = 350;
  ch_m0[7] = 400;
  ch_m0[8] = 450;
  ch_m0[9] = 500;
  ch_m0[10] = 550;
  ch_m0[11] = 600;
  ch_m0[12] = 650;
  ch_m0[13] = 700;
  ch_m0[14] = 750;
  ch_m0[15] = 750;
  ch_m0[16] = 0;
  
  ch_m12[0] = 170;
  ch_m12[1] = 168;
  ch_m12[2] = 167;
  ch_m12[3] = 165;
  ch_m12[4] = 163;
  ch_m12[5] = 161;
  ch_m12[6] = 158;
  ch_m12[7] = 156;
  ch_m12[8] = 154;
  ch_m12[9] = 152;
  ch_m12[10] = 150;
  ch_m12[11] = 148;
  ch_m12[12] = 147;
  ch_m12[13] = 145;
  ch_m12[14] = 144;
  ch_m12[15] = 0;
  ch_m12[16] = 0;
  
  TGraph* ch_gr = new TGraph(17,ch_m0,ch_m12);

  ch_gr->SetFillColor(3);
  ch_gr->SetLineColor(3);
  // ch_gr->SetLineWidth(3);
  ch_gr->SetFillStyle(3001);

  return ch_gr;

}


TGraph* set_lep_ch_tanBeta50(){

  double ch_m0[21];
  double ch_m12[21];

  ch_m0[0] = 200;
  ch_m0[1] = 250;
  ch_m0[2] = 300;
  ch_m0[3] = 350;
  ch_m0[4] = 400;
  ch_m0[5] = 450;
  ch_m0[6] = 500;
  ch_m0[7] = 550;
  ch_m0[8] = 600;
  ch_m0[9] = 650;
  ch_m0[10] = 700;
  ch_m0[11] = 750;
  ch_m0[12] = 800;
  ch_m0[13] =850;
  ch_m0[14] = 900;
  ch_m0[15] = 950;
  ch_m0[16] = 1000;
  ch_m0[17] = 1050;
  ch_m0[18] = 1100;
  ch_m0[19] = 1100;
  ch_m0[20] = 200;
 
  ch_m12[0] = 157;
  ch_m12[1] = 156;
  ch_m12[2] = 156;
  ch_m12[3] = 155;
  ch_m12[4] = 155;
  ch_m12[5] = 154;
  ch_m12[6] = 154;
  ch_m12[7] = 153;
  ch_m12[8] = 153;
  ch_m12[9] = 152;
  ch_m12[10] = 152;
  ch_m12[11] = 152;
  ch_m12[12] = 152;
  ch_m12[13] = 152;
  ch_m12[14] = 152;
  ch_m12[15] = 153;
  ch_m12[16] = 153;
  ch_m12[17] = 153;
  ch_m12[18] = 154;
  ch_m12[19] = 0;
  ch_m12[20] = 0;
  
  
  TGraph* ch_gr = new TGraph(21,ch_m0,ch_m12);

  ch_gr->SetFillColor(3);
  ch_gr->SetLineColor(3);
  ch_gr->SetFillStyle(3001);

  return ch_gr;

}




TGraph* set_lep_sl(Int_t tanBeta){

  // CMS SUSY Summer2010 implementation
  //  double sl_m0[] =  {0,  0, 30, 50, 60, 75, 80,90,100};
  //  double sl_m12[] = {0,245,240,220,200,150,100,50,0}; 
  
  //contour from D0 trilepton paper (PLB 680 (2009) 34-43)
  if (tanBeta==3){
    double sl_m0[] ={0,  0, 10, 20, 30, 40, 50, 60, 70, 77,88,95};
    double sl_m12[]={0,245,242,239,232,222,209,189,165,140,60,0};
    int n = 12;
  }
  //CMS PTDR-II
  //* Selectron_R line mass=99, ISASUGRA7.69, A0=0, m_top=175, tan(beta]=10
  if (tanBeta==10 || tanBeta==50){
    double sl_m0[]={ 0,  0, 11, 20, 24, 49, 70, 82,88,90};
    double sl_m12[]={0,240,237,233,230,200,150,100,50,0};
    int n = 10;
  }

  TGraph* lep_sl = new TGraph(n,sl_m0,sl_m12);

  lep_sl->SetFillColor(5);
  lep_sl->SetLineColor(5);
  lep_sl->SetFillStyle(3001);
  
  return lep_sl;
}


TGraph* set_tev_sg_cdf(Int_t tanBeta){

  //  double sg_m0[] =  {0,  0, 20, 50,100,150,200,250,300,350,400,450,500,550,600,600};
  //  double sg_m12[] = {0,160,169,170,160,155,150,122,116,112,110,106,105,100, 98,  0};
  //  int np=16;
  //New CHF from CDF plot in ICHEP2010 talk (E. Halkiadakis)
  double sg_m0[]= {0,  0, 30, 75,150,185,225,310,360,400,430,500,600,600};
  double sg_m12[]={0,162,168,170,160,150,130,120,109,108,100, 96, 95,  0};
  int np=14;

  TGraph* sg_gr = new TGraph(np,sg_m0,sg_m12);

  //  gStyle->SetHatchesLineWidth(3);

  sg_gr->SetFillColor(2);
  sg_gr->SetLineColor(2);
  //  sg_gr->SetLineWidth(3);
  sg_gr->SetFillStyle(3001); 

  return sg_gr;

}

TGraph* set_tev_sg_d0(Int_t tanBeta){

  //official D0 contour from P. Verdier
  double sgd_m0[]= {0,  0., 25., 80.,100.,150.,192.,250.,300. ,350.,400.,450. ,500.,600.,600.,0.};
  double sgd_m12[]={0,167.,167.,163.,162.,157.,149.,136.,125.5,116.,109.,106.5,105.,105.,  0.,0.};
  int npd=16;

  TGraph* sgd_gr = new TGraph(npd,sgd_m0,sgd_m12);

  gStyle->SetHatchesLineWidth(3);

  sgd_gr->SetFillColor(kMagenta+1);
  sgd_gr->SetLineColor(kMagenta+1);
  sgd_gr->SetLineWidth(3);
  sgd_gr->SetFillStyle(3335);

  return sgd_gr;

}

// TGraph* set_tev_tlp_cdf(Int_t tanBeta){
//   double tlp1_m0[] = {   0, 20, 40, 60, 70, 80, 90, 80, 70, 60};
//   double tlp1_m12[] = {170,185,200,215,220,215,210,190,175,160};
//   TGraph* tlp1_gr = new TGraph(10,tlp1_m0,tlp1_m12);

//   tlp1_gr->SetFillColor(4);
//   tlp1_gr->SetLineColor(4);
//   tlp1_gr->SetFillStyle(1001);

//   return tlp1_gr;
// }

// TGraph* set_tev_tlp_d0(Int_t tanBeta){
//   double tlp2_m0[] = {  70, 80, 90,100,105,110,120,130,140};
//   double tlp2_m12[] = {160,172,184,196,205,195,185,173,160};
//   TGraph* tlp2_gr = new TGraph(9,tlp2_m0,tlp2_m12);

//   tlp2_gr->SetFillColor(4);
//   tlp2_gr->SetFillStyle(1001); 

//   return tlp2_gr;

// }





TGraph* set_tev_stau(Int_t tanBeta){

    double st_m0_tanBeta3[] = {0,10,20,30,40,50,60,70,80,90,100,0};
    double st_m12_tanBeta3[] = {337,341,356,378,406,439,473,510,548,587,626,626};   

    double st_m0_tanBeta10[] = {0,10,20,30,40,50,60,70,80,90,100,0};
    double st_m12_tanBeta10[] = {213,220,240,275,312,351,393,435,476,518,559,559};

    double st_m0_tanBeta50[] = {200,210,220,230,240,250,260,270,280,290,310,325,200,200};
    double st_m12_tanBeta50[] = {206,226,246,267,288,310,332,354,376,399,450,500,500,206};


    TGraph* st_gr_tanBeta3 = new TGraph(12,st_m0_tanBeta3,st_m12_tanBeta3);
    TGraph* st_gr_tanBeta10 = new TGraph(12,st_m0_tanBeta10,st_m12_tanBeta10);
    TGraph* st_gr_tanBeta50 = new TGraph(14,st_m0_tanBeta50,st_m12_tanBeta50);

    st_gr_tanBeta3->SetFillColor(40);
    st_gr_tanBeta3->SetFillStyle(1001);

    st_gr_tanBeta50->SetFillColor(40);
    st_gr_tanBeta50->SetFillStyle(1001);
    
    st_gr_tanBeta10->SetFillColor(40);
    st_gr_tanBeta10->SetFillStyle(1001);


    if(tanBeta == 3)return st_gr_tanBeta3;
    if(tanBeta == 10)return st_gr_tanBeta10;
    if(tanBeta == 50)return st_gr_tanBeta50;

}




TF1* constant_squark(int tanBeta,int i){
//---lines of constant gluino/squark
  double coef1 = 0.35;
  double coef2[] = {5,5,4.6,4.1};

  char hname[200];

  sprintf(hname,"lnsq_%i",i); 

  
  TF1* lnsq = new TF1(hname,"sqrt([0]-x*x*[1]-[2])",0,1000);
  lnsq->SetParameter(0,(500+150*(i-1))*(500+150*(i-1))/coef2[i]);
  lnsq->SetParameter(1,1./coef2[i]);
  lnsq->SetParameter(2,-coef1*91*91*(2*TMath::Cos(TMath::ATan(tanBeta)))/coef2[i]);//--tanbeta=10 --> cos2beta = -99/101
  lnsq->SetLineWidth(1);


  lnsq->SetLineColor(kGray);

  return lnsq;
}


TF1* constant_gluino(int tanBeta,int i){
//---lines of constant gluino/squark
  double coef1 = 0.35;
  double coef2[] = {5,5,4.6,4.1};

  char hname[200];

  sprintf(hname,"lngl_%i",i); 
    
  TF1* lngl = new TF1(hname,"[0]+x*[1]",0,1000);
  lngl->SetParameter(0,(500+150.*(i-1))/2.4);
  lngl->SetParameter(1,-40./1400);
  lngl->SetLineWidth(1);
  lngl->SetLineColor(kGray);

  return lngl;
}


TLatex* constant_squark_text(Int_t it,TF1& lnsq,Int_t tanBeta_){
  char legnm[200];

  sprintf(legnm,"#font[92]{#tilde{q}(%i)GeV}",500+150*(it-1));
  Double_t place_x = 170;
  if(tanBeta_ == 50)place_x = 290;
  TLatex* t3 = new TLatex(place_x+10*(it-1),lnsq.Eval(place_x+10*(it-1))+5,legnm);
  t3->SetTextSize(0.03);
  t3->SetTextAngle(-8);
  t3->SetTextColor(kGray+2);


  
  return t3;
}

TLatex* constant_gluino_text(Int_t it,TF1& lngl){
  char legnm[200];

  sprintf(legnm,"#font[12]{#tilde{g}}#font[92]{(%i)GeV}",500+150*(it-1));
  TLatex* t4 = new TLatex(423,18+lngl.Eval(480),legnm);
  t4->SetTextSize(0.03);
  t4->SetTextAlign(13);
  t4->SetTextColor(kGray+2);

  return t4;
}



TLegend* makeStauLegend(Double_t txtsz,Int_t tanBeta_){
  Double_t ypos_1 = 0.86;
  Double_t ypos_2 = 0.88;
  Double_t xpos_1 = 0.16;
  Double_t xpos_2 = 0.17;
  if(tanBeta_ == 50){
    xpos_1 = 0.17;
    xpos_2 = 0.18;
    ypos_1 = 0.76;
    ypos_2 = 0.78;

  }
  TLegend* legst = new TLegend(xpos_1,ypos_1,xpos_2,ypos_2);
  legst->SetHeader("#tilde{#tau} = LSP");
  legst->SetFillStyle(0);
  legst->SetBorderSize(0);
  legst->SetTextSize(0.03);

  return legst;
}


TLegend* makeExpLegend(TGraph& sg_gr, TGraph& sgd_gr,TGraph& ch_gr,TGraph& sl_gr,TGraph& tev_sn,Double_t txtsz,Int_t tanbeta){
  TLegend* legexp = new TLegend(0.64,0.65,0.99,0.9,NULL,"brNDC");
  legexp->SetFillColor(0);
  legexp->SetShadowColor(0);
  legexp->SetTextSize(txtsz);
  legexp->SetBorderSize(0);

  sg_gr.SetLineColor(1);
  legexp->AddEntry(&sg_gr,"CDF  #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=5, #mu<0}","f"); 
  //  sgd_gr.SetLineColor(1);
  //  sgd_gr.SetLineWidth(1);

  legexp->AddEntry(&sgd_gr,"D0   #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=3, #mu<0}","f");  
  ch_gr.SetLineColor(1);
  legexp->AddEntry(&ch_gr,"LEP2   #tilde{#chi}_{1}^{#pm}","f");  
  
  sl_gr.SetLineColor(1);
  if(tanbeta != 50) legexp->AddEntry(&sl_gr,"LEP2   #tilde{#font[12]{l}}^{#pm}","f"); 
  if(tanbeta == 3) legexp->AddEntry(&tev_sn,"D0  #chi^{#pm}_{1}, #chi^{0}_{2}","f");  
 

  return legexp;

}


TGraph* getGraph(TH2F* h1, double level){

  h1->SetContour(1);
  h1->SetContourLevel(0,level);
  h1->Draw("CONT LIST");
  gPad->Update();
  
  TObjArray* contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");

  cout << "contours: " << contours << endl;
  // Draw contours
  TList* graphList = (TList*)(contours->At(0));
  cout << "number of graphs: " << graphList->GetSize() << endl;
  for (int igraph = 0; igraph<graphList->GetSize();++igraph) {
    TGraph* myGraph = (TGraph*)graphList->At(igraph);

    std::cout << " - graph " << igraph << " has " << myGraph->GetN() << " points" << std::endl;
    if (myGraph->GetN() > 50){
      cout << "Drawing " << myGraph->GetN() <<" points" << endl;
      //      myGraph->Print("all");
      //      myGraph->SetLineColor(9);
      //      myGraph->SetLineColor(46);
      //myGraph->SetLineWidth(3);
      //      myGraph->SetLineStyle(2);
      //      myGraph->SetMarkerStyle(20);
      //      myGraph->SetMarkerColor(9);

      //      myGraph->Draw("C");
      //      cvsSys->Update();

    //    TString graphName("graph"+name+"_");
    //    graphName += igraph;
    //    myGraph->SetName(graphName);
      break;
    }
  }
  return myGraph;
}

TH2F* fillHoles(TH2F* h2){

  int nx = h2->GetNbinsX();
  int ny = h2->GetNbinsY();
  cout << "Nbins: " << nx << " " << ny << endl;

  for (int i=1;i<nx+1;i++){
    for (int j=ny;j>0;j--){
      int pos = j;
      int count =0;
      if (h2->GetBinContent(i,j)==0 && h2->GetBinContent(i-1,j)==-1 && h2->GetBinContent(i+1,j)==-1) h2->SetBinContent(i,j,-1.);

      if (h2->GetBinContent(i,j)==1 
	  && h2->GetBinContent(i-1,j)==-1 
	  && h2->GetBinContent(i+1,j)==-1 
	  && h2->GetBinContent(i,j+1)==-1 
	  && h2->GetBinContent(i,j-1)==-1) 
	h2->SetBinContent(i,j,-1.);

      //      if (h2->GetBinContent(i,j)<0){
      if (h2->GetBinContent(i,j)<0 && h2->GetBinContent(i,j-1)<0){
       	//      cout <<" bin content " << i <<" " <<  j<< " " << pos << " " << h2->GetBinContent(i,j) <<endl;
       	for (int k=pos;k>0;k--){
       	  h2->SetBinContent(i,k,-1.);
	}
	break;
      }
      
      
    }
  }
  return h2;
}
