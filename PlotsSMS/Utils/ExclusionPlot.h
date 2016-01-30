#ifndef MAKE_EXCLUSIONPLOT_HH
#define MAKE_EXCLUSIONPLOT_HH

#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TSpline.h"
#include "TGraphErrors.h"

#include <stdlib.h>
#include <vector>





void ExclusionPlot();
void CommandMSUGRA(TString plotName,Int_t tanBeta);
void setPlottingStyle(TH1F& hsig);

TGraph* set_lep_ch(Int_t tanBeta);
TGraph* set_lep_ch_tanBeta3();
TGraph* set_lep_ch_tanBeta10();
TGraph* set_lep_ch_tanBeta50();
TGraph* set_lep_sl(Int_t tanBeta);//slepton curve
TGraph* set_tev_sg_cdf(Int_t tanBeta);//squark gluino cdf
TGraph* set_tev_sg_d0(Int_t tanBeta);//squark gluino d0
TGraph* set_tev_stau(Int_t tanBeta);//stau 
TGraph* set_sneutrino_d0_1(Int_t tanBeta);
TGraph* set_sneutrino_d0_2(Int_t tanBeta);
TF1* constant_squark(int tanBeta,int i);
TF1* constant_gluino(int tanBeta,int i);
TLatex* constant_squark_text(Int_t it,TF1& lnsq,Int_t tanBeta);
TLatex* constant_gluino_text(Int_t it,TF1& lngl);
TLegend* makeStauLegend(Double_t txtsz,Int_t tanBeta_);
TLegend* makeExpLegend(TGraph& sg_gr, TGraph& sgd_gr,TGraph& ch_gr,TGraph& sl_gr,TGraph&,Double_t txtsz,Int_t tanbeta);







TGraphErrors* getLO_tanBeta3(){



  Int_t nl = 9;
  Double_t xl[9];
  Double_t yl[9];
  Double_t exl[9];
  Double_t eyl[9];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
 
   xl[0] = 0;
  yl[0] = 265;
  xl[1] = 100;
  yl[1] = 258;
  xl[2] = 200;
  yl[2] = 250;
  xl[3] = 250;
  yl[3] = 240;
  xl[4] = 300;
  yl[4] = 210;
  xl[5] = 340;
  yl[5] = 177;
  xl[6] = 400;
  yl[6] = 140;
  xl[7] = 450;
  yl[7] = 120;
  xl[8] = 520;
  yl[8] =100;

  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kGreen+2);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kGreen+2);
  s->SetLineStyle(4);
  s->SetLineWidth(3);
  

  return gr1;



}

TGraphErrors* getLO_tanBeta10(){



  Int_t nl = 10;
  Double_t xl[10];
  Double_t yl[10];
  Double_t exl[10];
  Double_t eyl[10];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
  xl[0] = 0;
  yl[0] = 270;
  xl[1] = 100;
  yl[1] = 260;
  xl[2] = 200;
  yl[2] = 250;
  xl[3] = 250;
  yl[3] = 240;
  xl[4] = 300;
  yl[4] = 210;
  xl[5] = 350;
  yl[5] = 174;
  xl[6] = 400;
  yl[6] = 147;
  xl[7] = 450;
  yl[7] = 127;
  xl[8] = 500;
  yl[8] =115;
  xl[9] = 520;
  yl[9] = 112;
  
  
  
  
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kGreen+2);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kGreen+2);
  s->SetLineStyle(4);
  s->SetLineWidth(3);
  

  return gr1;



}

TGraphErrors* getLO_tanBeta50(){



  Int_t nl = 10;
  Double_t xl[10];
  Double_t yl[10];
  Double_t exl[10];
  Double_t eyl[10];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  

  xl[0] = 200;
  yl[0] = 239;
  xl[1] = 210;
  yl[1] = 249;
  xl[2] = 229;
  yl[2] = 260;
  xl[3] = 250;
  yl[3] = 245;
  xl[4] = 300;
  yl[4] = 210;
  xl[5] = 350;
  yl[5] = 180;
  xl[6] = 400;
  yl[6] = 160;
  xl[7] = 450;
  yl[7] = 150;
  xl[8] = 500;
  yl[8] =140;
  xl[9] = 520;
  yl[9] = 137;
  
  
  
  
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kGreen+2);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kGreen+2);
  s->SetLineStyle(4);
  s->SetLineWidth(3);
  

  return gr1;



}



TGraphErrors* getExpected_NLO_tanBeta3(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
    xl[0] = 0;
  yl[0] = 283;
  xl[1] = 100;
  yl[1] = 280;
  xl[2] = 150;
  yl[2] = 279;
  xl[3] = 200;
  yl[3] = 275;
  xl[4] = 250;
  yl[4] = 270;
  xl[5] = 300;
  yl[5] = 255;
  xl[6] = 350;
  yl[6] = 225;
  xl[7] = 400;
  yl[7] = 195;
  xl[8] = 450;
  yl[8] = 175;
  xl[9] = 500;
  yl[9] = 155;
  xl[10] = 550;
  yl[10] = 150;

 
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}

TGraphErrors* getExpected_NLO_tanBeta10(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 0;
  yl[0] = 283;
  xl[1] = 100;
  yl[1] = 280;
  xl[2] = 150;
  yl[2] = 279;
  xl[3] = 200;
  yl[3] = 275;
  xl[4] = 250;
  yl[4] = 270;
  xl[5] = 300;
  yl[5] = 255;
  xl[6] = 350;
  yl[6] = 225;
  xl[7] = 400;
  yl[7] = 195;
  xl[8] = 450;
  yl[8] = 175;
  xl[9] = 500;
  yl[9] = 165;
  xl[10] = 550;
  yl[10] = 150;



  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}


TGraphErrors* getExpected_NLO_tanBeta10_2(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 0;
  yl[0] = 300;
  xl[1] = 100;
  yl[1] = 288;
  xl[2] = 150;
  yl[2] = 283;
  xl[3] = 200;
  yl[3] = 278;
  xl[4] = 250;
  yl[4] = 268;
  xl[5] = 300;
  yl[5] = 255;
  xl[6] = 350;
  yl[6] = 225;
  xl[7] = 400;
  yl[7] = 197;
  xl[8] = 450;
  yl[8] = 175;
  xl[9] = 500;
  yl[9] = 160;
  xl[10] = 650;
  yl[10] = 145;


  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}


TGraphErrors* getExpected_NLO_tanBeta3_2(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 0;
  yl[0] = 295;
  xl[1] = 100;
  yl[1] = 290;
  xl[2] = 150;
  yl[2] = 283;
  xl[3] = 200;
  yl[3] = 275;
  xl[4] = 250;
  yl[4] = 270;
  xl[5] = 300;
  yl[5] = 255;
  xl[6] = 350;
  yl[6] = 225;
  xl[7] = 400;
  yl[7] = 195;
  xl[8] = 450;
  yl[8] = 175;
  xl[9] = 500;
  yl[9] = 160;
  xl[10] = 600;
  yl[10] = 140;


  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}


TGraphErrors* getExpected_NLO_tanBeta10_up(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 0;
  yl[0] = 318;
  xl[1] = 100;
  yl[1] = 318;
  xl[2] = 150;
  yl[2] = 312;
  xl[3] = 200;
  yl[3] = 308;
  xl[4] = 250;
  yl[4] = 300;
  xl[5] = 300;
  yl[5] = 290;
  xl[6] = 350;
  yl[6] = 265;
  xl[7] = 400;
  yl[7] = 239;
  xl[8] = 450;
  yl[8] = 208;
  xl[9] = 500;
  yl[9] = 190;
  xl[10] = 620;
  yl[10] = 180;


  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}



TGraphErrors* getExpected_NLO_tanBeta3_up(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 0;
  yl[0] = 320;
  xl[1] = 100;
  yl[1] = 313;
  xl[2] = 150;
  yl[2] = 310;
  xl[3] = 200;
  yl[3] = 307;
  xl[4] = 250;
  yl[4] = 300;
  xl[5] = 300;
  yl[5] = 290;
  xl[6] = 350;
  yl[6] = 258;
  xl[7] = 400;
  yl[7] = 230;
  xl[8] = 450;
  yl[8] = 200;
  xl[9] = 500;
  yl[9] = 188;
  xl[10] = 650;
  yl[10] = 175;


  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}



TGraphErrors* getExpected_NLO_tanBeta10_low(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 0;
  yl[0] = 275;
  xl[1] = 100;
  yl[1] = 268;
  xl[2] = 150;
  yl[2] = 262;
  xl[3] = 200;
  yl[3] = 255;
  xl[4] = 250;
  yl[4] = 250;
  xl[5] = 300;
  yl[5] = 230;
  xl[6] = 350;
  yl[6] = 190;
  xl[7] = 400;
  yl[7] = 160;
  xl[8] = 450;
  yl[8] = 148;
  xl[9] = 500;
  yl[9] = 140;
  xl[10] = 650;
  yl[10] = 120;


  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}

TGraphErrors* getExpected_NLO_tanBeta3_low(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 0;
  yl[0] = 265;
  xl[1] = 100;
  yl[1] = 260;
  xl[2] = 150;
  yl[2] = 258;
  xl[3] = 200;
  yl[3] = 255;
  xl[4] = 250;
  yl[4] = 250;
  xl[5] = 300;
  yl[5] = 230;
  xl[6] = 350;
  yl[6] = 185;
  xl[7] = 400;
  yl[7] = 155;
  xl[8] = 450;
  yl[8] = 138;
  xl[9] = 480;
  yl[9] = 128;
  xl[10] = 530;
  yl[10] = 118;

  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  
  return gr1;
}




TGraphErrors* getExpected_NLO_tanBeta50(){

 Int_t nl = 10;
  Double_t xl[10];
  Double_t yl[10];
  Double_t exl[10];
  Double_t eyl[10];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 200;
  yl[0] = 287;
  xl[1] = 220;
  yl[1] = 287;
  xl[2] = 245;
  yl[2] = 287;
  xl[3] = 270;
  yl[3] = 265;
  xl[4] = 300;
  yl[4] = 245;
  xl[5] = 350;
  yl[5] = 222;
  xl[6] = 400;
  yl[6] = 197;
  xl[7] = 450;
  yl[7] = 180;
  xl[8] = 500;
  yl[8] = 168;
  xl[9] = 550;
  yl[9] = 145;



  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}




TGraphErrors* getExpected_NLO_tanBeta50_2(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 200;
  yl[0] = 287;
  xl[1] = 220;
  yl[1] = 287;
  xl[2] = 245;
  yl[2] = 287;
  xl[3] = 270;
  yl[3] = 270;
  xl[4] = 300;
  yl[4] = 255;
  xl[5] = 350;
  yl[5] = 225;
  xl[6] = 400;
  yl[6] = 203;
  xl[7] = 450;
  yl[7] = 180;
  xl[8] = 500;
  yl[8] = 158;
  xl[9] = 550;
  yl[9] = 150;
  xl[10] = 650;
  yl[10] = 140;



  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}


TGraphErrors* getExpected_NLO_tanBeta50_up(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 200;
  yl[0] = 340;
  xl[1] = 220;
  yl[1] = 320;
  xl[2] = 245;
  yl[2] = 320;
  xl[3] = 270;
  yl[3] = 310;
  xl[4] = 300;
  yl[4] = 290;
  xl[5] = 350;
  yl[5] = 270;
  xl[6] = 400;
  yl[6] = 245;
  xl[7] = 450;
  yl[7] = 220;
  xl[8] = 500;
  yl[8] = 200;
  xl[9] = 550;
  yl[9] = 180;
  xl[10] = 650;
  yl[10] = 170;


  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}

TGraphErrors* getExpected_NLO_tanBeta50_low(){

 Int_t nl = 8;
  Double_t xl[8];
  Double_t yl[8];
  Double_t exl[8];
  Double_t eyl[8];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 200;
  yl[0] = 280;
  xl[1] = 220;
  yl[1] = 270;
  xl[2] = 245;
  yl[2] = 262;
  xl[3] = 270;
  yl[3] = 242;
  xl[4] = 300;
  yl[4] = 227;
  xl[5] = 350;
  yl[5] = 197;
  xl[6] = 400;
  yl[6] = 180;
  xl[7] = 450;
  yl[7] = 150;
  //  xl[8] = 500;
  //yl[8] = 168;
  //xl[9] = 550;
  //yl[9] = 145;



  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}






TGraphErrors* getObserved_NLO_tanBeta3(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  

  
  xl[0] = 0;
  yl[0] = 274;
  xl[1] = 100;
  yl[1] = 270;
  xl[2] = 150;
  yl[2] = 268;
  xl[3] = 200;
  yl[3] = 265;
  xl[4] = 250;
  yl[4] = 255;
  xl[5] = 300;
  yl[5] = 230;
  xl[6] = 350;
  yl[6] = 195;
  xl[7] = 400;
  yl[7] = 160;
  xl[8] = 450;
  yl[8] = 140;
  xl[9] = 480;
  yl[9] = 130;
  xl[10] = 530;
  yl[10] = 120;
 
  
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kRed);
  // s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}

TGraphErrors* getObserved_NLO_tanBeta3_2(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  

  
  xl[0] = 0;
  yl[0] = 277;
  xl[1] = 100;
  yl[1] = 272;
  xl[2] = 150;
  yl[2] = 268;
  xl[3] = 200;
  yl[3] = 263;
  xl[4] = 250;
  yl[4] = 258;
  xl[5] = 300;
  yl[5] = 237;
  xl[6] = 350;
  yl[6] = 198;
  xl[7] = 400;
  yl[7] = 168;
  xl[8] = 450;
  yl[8] = 152;
  xl[9] = 480;
  yl[9] = 145;
  xl[10] = 530;
  yl[10] = 133;
 
  
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kRed);
  // s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}



TGraphErrors* getObserved_NLO_tanBeta10(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 0;
  yl[0] = 278;
  xl[1] = 100;
  yl[1] = 270;
  xl[2] = 150;
  yl[2] = 267;
  xl[3] = 200;
  yl[3] = 262;
  xl[4] = 250;
  yl[4] = 250;
  xl[5] = 300;
  yl[5] = 225;
  xl[6] = 350;
  yl[6] = 192;
  xl[7] = 400;
  yl[7] = 163;
  xl[8] = 450;
  yl[8] = 148;
  xl[9] = 500;
  yl[9] = 140;
  xl[10] = 520;
  yl[10] = 137;
  
 
 
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kRed);
  //  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}


TGraphErrors* getObserved_NLO_tanBeta10_2(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 0;
  yl[0] = 288;
  xl[1] = 100;
  yl[1] = 280;
  xl[2] = 150;
  yl[2] = 275;
  xl[3] = 200;
  yl[3] = 268;
  xl[4] = 250;
  yl[4] = 260;
  xl[5] = 300;
  yl[5] = 237;
  xl[6] = 350;
  yl[6] = 203;
  xl[7] = 400;
  yl[7] = 172;
  xl[8] = 450;
  yl[8] = 156;
  xl[9] = 500;
  yl[9] = 147;
  xl[10] = 650;
  yl[10] = 131;
 
 
 
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kRed);
  //  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}


TGraphErrors* getObserved_NLO_tanBeta50(){

 Int_t nl = 10;
  Double_t xl[10];
  Double_t yl[10];
  Double_t exl[10];
  Double_t eyl[10];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  

  xl[0] = 200;
  yl[0] = 243;
  xl[1] = 220;
  yl[1] = 266;
  xl[2] = 235;
  yl[2] = 278;
  xl[3] = 250;
  yl[3] = 267;
  xl[4] = 300;
  yl[4] = 230;
  xl[5] = 350;
  yl[5] = 205;
  xl[6] = 400;
  yl[6] = 184;
  xl[7] = 450;
  yl[7] = 168;
  xl[8] = 500;
  yl[8] = 156;
  xl[9] = 520;
  yl[9] = 148;
  
 
 
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kRed);
  //  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}


TGraphErrors* getObserved_NLO_tanBeta50_2(){

 Int_t nl = 9;
  Double_t xl[9];
  Double_t yl[9];
  Double_t exl[9];
  Double_t eyl[9];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  

  xl[0] = 200;
  yl[0] = 240;
  xl[1] = 220;
  yl[1] = 266;
  xl[2] = 235;
  yl[2] = 275;
  xl[3] = 250;
  yl[3] = 270;
  xl[4] = 300;
  yl[4] = 240;
  xl[5] = 350;
  yl[5] = 210;
  xl[6] = 400;
  yl[6] = 190;
  xl[7] = 450;
  yl[7] = 168;
  xl[8] = 500;
  yl[8] = 150;

  
 
 
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kRed);
  //  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}




TH1F* getHisto_1d( TString path,
	       TString nameHist,
	       TString nameFile ) {
  TString name = path + nameFile;
  TFile* file =  new TFile(name);
  //  TDirectory* dir = (TDirectory*)file->Get(Dirname);
  TH2F* hist = (TH2F*)file->Get(nameHist);
  if (!hist) {
    std::cout << " name: " << nameHist
	      << " file: " << nameFile
	      << std::endl;
    abort();

  }

  

  TH1F* Onedhist = new TH1F(nameHist,nameHist,int(hist->GetXaxis()->GetNbins()+0.5),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

  for(int x = 0;x < hist->GetXaxis()->GetNbins();x++){

    bool firsthit = false;
    for(int y = hist->GetYaxis()->GetNbins(); y>0; y--){

      double y_height = hist->GetYaxis()->GetXmin()+hist->GetYaxis()->GetBinWidth(0)*y;

      if(firsthit == false && hist->GetBinContent(x+1,y) > 0){
	Onedhist->SetBinContent(x+1,y_height);
	firsthit = true;
	}


    }

  }


  Onedhist->SetLineWidth(1);

  Onedhist->GetXaxis()->SetTitleSize(0.055);
  Onedhist->GetYaxis()->SetTitleSize(0.055);
  Onedhist->GetXaxis()->SetLabelSize(0.05);
  Onedhist->GetYaxis()->SetLabelSize(0.05);
  Onedhist->SetStats(kFALSE);

  
  

  return Onedhist;
}


#endif
