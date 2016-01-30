#include "SusyScan.h"
#include "GeneratorMasses.h"
#include "simpleGenMasses.h"
#include "ExclusionPlot.h"

#include "TGraph.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRint.h"
#include "TROOT.h"
#include "TPad.h"

#include <cmath>
#include <iostream>


double Luminosity = 1143; //[pb^-1]
double Mzero(const SusyScan* p){ return p->Mzero; }
double Mhalf(const SusyScan* p){ return p->Mhalf; }
double MGluino(const SusyScan* p){ return p->Mgluino; }
double MSquark(const SusyScan* p){ return p->Msquark; }
double MLSP(const SusyScan* p){ return p->Mlsp; }
//double MGluino(const SusyScan* p){ return p->MGL; }
double MSquarkL(const SusyScan* p){ return p->MUL; }
double MSquarkR(const SusyScan* p){ return p->MUR; }
double MChi1(const SusyScan* p){ return p->MZ1; }
double MChi2(const SusyScan* p){ return p->MZ2; }
double MChi3(const SusyScan* p){ return p->MZ3; }
double MChi4(const SusyScan* p){ return p->MZ4; }
double MCha1(const SusyScan* p){ return p->MW1; }
double MCha2(const SusyScan* p){ return p->MW2; }
double MTAU1(const SusyScan* p){ return p->MTAU1; }
double ChargedLSP(const SusyScan* p){ return (fabs(p->MTAU1) < fabs(p->MZ1) ? 0.01 : 1); }
double SignalKfactor(const SusyScan* p){return p->signal_kfactor; }

double Xsection(const SusyScan* p){ return p->Xsection; }
double ExpXsecLimit(const SusyScan* p){ return p->ExpXsecLimit; }
double ExpXsecBayes(const SusyScan* p){ return p->ExpXsecBayes; }
double ObsXsecLimit(const SusyScan* p){ return p->ObsXsecLimit; }
double ObsXsecBayes(const SusyScan* p){ return p->ObsXsecBayes; }
double ExpXsecLimit_n1(const SusyScan* p){ return p->ExpXsecLimit_n1; }
double ExpXsecLimit_p1(const SusyScan* p){ return p->ExpXsecLimit_p1; }
double ExpXsecLimit_n2(const SusyScan* p){ return p->ExpXsecLimit_n2; }
double ExpXsecLimit_p2(const SusyScan* p){ return p->ExpXsecLimit_p2; }
double Signal(const SusyScan* p){ return p->signal; }
double SignalUncertainty(const SusyScan* p){ return p->signal_uncertainty; }
double SignalRelUncertainty(const SusyScan* p){ return p->signal_uncertainty/p->signal; }
double ExpExclusion(const SusyScan* p){ return (ExpXsecLimit(p)<Xsection(p)&&ExpXsecLimit(p)>0.01?1:0.01); }
double CLsObsExclusion(const SusyScan* p){ return (ObsXsecLimit(p)<Xsection(p)&&ObsXsecLimit(p)>0.01?1:0.01); }
double BayesObsExclusion(const SusyScan* p){ return (ObsXsecBayes(p)<Xsection(p)&&ObsXsecBayes(p)>0.01?1:0.01); }
double ExpExclCL(const SusyScan* p){ return (p->CLs_b_xsec<=0.05 ? 1:0.01); }
double ExpExclCLm2sigma(const SusyScan* p){ return (p->CLs_b_n2_xsec<=0.05 ? 1:0.01); }
double ExpExclCLm1sigma(const SusyScan* p){ return (p->CLs_b_n1_xsec<=0.05 ? 1:0.01); }
double ExpExclCLp1sigma(const SusyScan* p){ return (p->CLs_b_p1_xsec<=0.05 ? 1:0.01); }
double ExpExclCLp2sigma(const SusyScan* p){ return (p->CLs_b_p2_xsec<=0.05 ? 1:0.01); }
double ObsExclCL(const SusyScan* p){ return (p->CLs_xsec<=0.05 ? 1:0.01); }
double SoverSqrtB(const SusyScan* p){ return p->signal/(sqrt(p->background)+p->background_uncertainty+p->signal_uncertainty); }
double SignalAcceptance(const SusyScan* p){ return  p->signal / (Luminosity*Xsection(p)); }
double SignalContamination(const SusyScan* p){return p->signal_contamination; }

double NLOXsection(const SusyScan* p){ return p->Xsection * p->signal_kfactor; }
double NLOExpXsecLimit(const SusyScan* p){ return p->NLO_ExpXsecLimit; }
double NLOExpXsecBayes(const SusyScan* p){ return p->NLO_ExpXsecBayes; }
double NLOObsXsecLimit(const SusyScan* p){ return p->NLO_ObsXsecLimit; }
double NLOObsXsecBayes(const SusyScan* p){ return p->NLO_ObsXsecBayes; }
double NLOExpXsecLimit_n1(const SusyScan* p){ return p->NLO_ExpXsecLimit_n1; }
double NLOExpXsecLimit_p1(const SusyScan* p){ return p->NLO_ExpXsecLimit_p1; }
double NLOExpXsecLimit_n2(const SusyScan* p){ return p->NLO_ExpXsecLimit_n2; }
double NLOExpXsecLimit_p2(const SusyScan* p){ return p->NLO_ExpXsecLimit_p2; }
double NLOSignal(const SusyScan* p){ return p->NLO_signal; }
double NLOSignalUncertainty(const SusyScan* p){ return p->NLO_signal_uncertainty; }
double NLOSignalRelUncertainty(const SusyScan* p){ return p->NLO_signal_uncertainty/p->NLO_signal; }
double NLOExpExclusion(const SusyScan* p){ return (NLOExpXsecLimit(p)<NLOXsection(p)&&NLOExpXsecLimit(p)>0.01?1:0.01); }
double NLOCLsObsExclusion(const SusyScan* p){ return (NLOObsXsecLimit(p)<NLOXsection(p)&&NLOObsXsecLimit(p)>0.01?1:0.01); }
double NLOBayesObsExclusion(const SusyScan* p){ 
	return (NLOObsXsecBayes(p)<NLOXsection(p)&&NLOObsXsecBayes(p)>0.001?1:0.01); }
double NLOExpCL(const SusyScan* p){ return (p->NLO_CLs_b_xsec); }
double NLOExpExclCL(const SusyScan* p) { return ((p->NLO_CLs_b_xsec <= p->Xsection * p->signal_kfactor) ? 1:0.01); } //{ return (p->NLO_CLs_b_xsec<=0.05 ? 1:0.01); }
double NLOExpExclCLm2sigma(const SusyScan* p) { return ((p->NLO_CLs_b_n2_xsec <= p->Xsection * p->signal_kfactor) ? 1:0.01); }//{ return (p->NLO_CLs_b_n2_xsec<=0.05 ? 1:0.01); }
double NLOExpExclCLm1sigma(const SusyScan* p) { return ((p->NLO_CLs_b_n1_xsec <= p->Xsection * p->signal_kfactor) ? 1:0.01); }//{ return (p->NLO_CLs_b_n1_xsec<=0.05 ? 1:0.01); }
double NLOExpExclCLp1sigma(const SusyScan* p) { return ((p->NLO_CLs_b_p1_xsec <= p->Xsection * p->signal_kfactor) ? 1:0.01); }//{ return (p->NLO_CLs_b_p1_xsec<=0.05 ? 1:0.01); }
double NLOExpExclCLp2sigma(const SusyScan* p) { return ((p->NLO_CLs_b_p2_xsec <= p->Xsection * p->signal_kfactor) ? 1:0.01); }//{ return (p->NLO_CLs_b_p2_xsec<=0.05 ? 1:0.01); }
double NLOObsExclCL(const SusyScan* p){ return ((p->NLO_CLs_xsec <= p->Xsection * p->signal_kfactor) ? 1:0.01); }//return (p->NLO_CLs_xsec<=0.05 ? 1:0.01); }
double NLOSoverSqrtB(const SusyScan* p){ return p->NLO_signal/(sqrt(p->background)+p->background_uncertainty+p->NLO_signal_uncertainty); }
double NLOSignalAcceptance(const SusyScan* p){ return  p->NLO_signal / (Luminosity*NLOXsection(p)); }
double NLOSignalContamination(const SusyScan* p){return p->NLO_signal_contamination; }

double NLOObsxSecCL(const SusyScan* p){ return p->NLO_CLs_xsec; }
double NLOExpxSecCL(const SusyScan* p){ return p->NLO_CLs_b_xsec; }

double NLOExpxSecp1sigma(const SusyScan* p){ return p->NLO_CLs_b_p1_xsec; }
double NLOExpxSecm1sigma(const SusyScan* p){ return p->NLO_CLs_b_n1_xsec; }

double Mzero(const GeneratorMasses* p){ return p->Mzero; }
double Mhalf(const GeneratorMasses* p){ return p->Mhalf; }
double MGluino(const GeneratorMasses* p){ return p->MGL; }
double MSquarkL(const GeneratorMasses* p){ return p->MUL; }
double MSquarkR(const GeneratorMasses* p){ return p->MUR; }
double MChi1(const GeneratorMasses* p){ return p->MZ1; }
double MChi2(const GeneratorMasses* p){ return p->MZ2; }
double MChi3(const GeneratorMasses* p){ return p->MZ3; }
double MChi4(const GeneratorMasses* p){ return p->MZ4; }
double MCha1(const GeneratorMasses* p){ return p->MW1; }
double MCha2(const GeneratorMasses* p){ return p->MW2; }



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


TGraph* set_ProspinoExcl(const TH1D *ProspinoXsec, const TH2D *limit, double multiplier =1.0){

   TAxis *xLimit = (TAxis*) limit->GetXaxis();
   TAxis *yLimit = (TAxis*) limit->GetYaxis();
   int nxBinsLimit = xLimit->GetNbins();
   int nyBinsLimit = yLimit->GetNbins();

   std::vector<double> exclxPointsProspinoVec, exclyPointsProspinoVec;

   TAxis *xProspino = (TAxis*) ProspinoXsec->GetXaxis();
   for(int ib=0; ib<nxBinsLimit; ib++){
      double mother = xLimit->GetBinLowEdge(ib+1);
      double nextMother = xLimit->GetBinUpEdge(ib+1);
//      int iBinProspino = xProspino->FindBin(mother);
      int iBinProspino = xProspino->FindBin((mother+nextMother)*0.5);
      double xSec = ProspinoXsec->GetBinContent( iBinProspino );
      double xSecToExcl = xSec*multiplier;
      int findExclPoint = 0;
      for(int jb=0; jb<nyBinsLimit; jb++){
         double xSecUL = limit->GetBinContent(ib+1, nyBinsLimit-jb);
//         if( xSecUL >= xSecToExcl && !findExclPoint ){
         if( xSecUL == 0 ) continue;
         if( xSecUL < xSecToExcl && !findExclPoint ){
            findExclPoint =1;
            double lsp = yLimit->GetBinLowEdge(nyBinsLimit-jb);
            double nextLsp = yLimit->GetBinUpEdge(nyBinsLimit-jb);
//            exclxPointsProspinoVec.push_back(0.5*(mother+nextMother));
//            exclyPointsProspinoVec.push_back(0.5*(lsp+nextLsp));
            exclxPointsProspinoVec.push_back(nextMother);
            exclyPointsProspinoVec.push_back(nextLsp);
         }
      }
   }

   int size = (int) exclxPointsProspinoVec.size();
   double exclxPointsProspinoArr[size+1], exclyPointsProspinoArr[size+1];
   for(int is=0; is<size; is++){
      exclxPointsProspinoArr[is] = exclxPointsProspinoVec[is];
      exclyPointsProspinoArr[is] = exclyPointsProspinoVec[is];
      if( is== size-1 ){
         exclxPointsProspinoArr[is+1] = exclxPointsProspinoVec[is];
         exclyPointsProspinoArr[is+1] = 0.0;
      }
   }

   TGraph * excl_gr = new TGraph(size+1, exclxPointsProspinoArr, exclyPointsProspinoArr);

   return excl_gr;
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
  double ch_m0[] ={0  ,100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,2000,0.};  
  double ch_m12[]={162,162, 161, 160, 160, 159, 158, 157, 156, 155, 154, 154, 153, 152, 150, 149, 148,  147 , 146 , 146 , 146 , 147 , 148 , 149 , 151 , 154 ,159, 0., 0.}; 

  TGraph* ch_gr = new TGraph(29,ch_m0,ch_m12);

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
  TGraph* lep_sl = 0;
  if (tanBeta==3){
    double sl_m0[] ={0,  0, 10, 20, 30, 40, 50, 60, 70, 77,88,95};
    double sl_m12[]={0,245,242,239,232,222,209,189,165,140,60,0};
    int n = 12;
    lep_sl = new TGraph(n,sl_m0,sl_m12);
  }
  //CMS PTDR-II
  //* Selectron_R line mass=99, ISASUGRA7.69, A0=0, m_top=175, tan(beta]=10
  if (tanBeta==10 || tanBeta==50){
    double sl_m0[]={ 0,  0, 11, 20, 24, 49, 70, 82,88,90};
    double sl_m12[]={0,240,237,233,230,200,150,100,50,0};
    int n = 10;
    lep_sl = new TGraph(n,sl_m0,sl_m12);
  }


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
  //sgd_gr->SetFillStyle(3001);

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
  // taken from Frederic Ronga's calculation

  double st_m0_tanBeta10[] = {0,   10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110 , 130,  147, 0, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590, 610, 630, 6500};
  double st_m12_tanBeta10[] = {213,220,240,275,312,351,393,435,476,518, 559, 600., 682., 750.,750, 842., 921., 999., 1076, 1152, 1228, 1304, 1378, 1453, 1527, 1600, 1673, 1746, 1818, 1890, 1962, 2034, 2105, 2175, 2246, 2316, 2386, 2456, 2526, 2595}; 


  double st_m0_tanBeta40[] = {240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880};
  double st_m12_tanBeta40[] = {186, 256, 329, 400, 470, 537, 603, 666, 727, 787, 845, 902, 958, 1013, 1067, 1121, 1174, 1226, 1278, 1330, 1381, 1431, 1481, 1531, 1581, 1630, 1679, 1728, 1779, 1825, 1874, 1920, 1971};

  TGraph* st_gr_tanBeta10 = new TGraph(15,st_m0_tanBeta10,st_m12_tanBeta10);
  TGraph* st_gr_tanBeta40 = new TGraph(10,st_m0_tanBeta40,st_m12_tanBeta40);

    
  st_gr_tanBeta40->SetFillColor(40);
  st_gr_tanBeta40->SetFillStyle(1001);
    
  st_gr_tanBeta10->SetFillColor(40);
  st_gr_tanBeta10->SetFillStyle(1001);


  if(tanBeta == 10)return st_gr_tanBeta10;
  if(tanBeta == 40)return st_gr_tanBeta40;
}




//From Sanjay
TF1* constant_squark(int tanBeta,int i){
  //---lines of constant gluino/squark.
  // Min squark mass from 1st and 2nd generations using fit for tanbeta = 10.

  double coef1[] = {2.67058e+04, 6.39642e+04, 1.16565e+05, 1.95737e+05, 2.86190e+05};
  double coef2[] = {1.98772e-01, 2.11242e-01, 2.17734e-01, 2.39535e-01, 2.39768e-01};
  double coef3[] = {2.67058e+04, 6.39641e+04, 1.16565e+05, 1.95736e+05, 2.86189e+05};
 
  char hname[200];

  sprintf(hname,"lnsq_%i",i);
  TF1* lnsq = new TF1(hname,"sqrt([0]-x*x*[1]+[2])",0,2000);
  lnsq->SetParameter(0,coef1[i-1]);
  lnsq->SetParameter(1,coef2[i-1]);
  lnsq->SetParameter(2,coef3[i-1]);
  lnsq->SetLineWidth(1);
  lnsq->SetLineColor(kGray);

  return lnsq;
}

TF1* constant_gluino(int tanBeta,int i){
  //---lines of constant gluino/squark
  char hname[200];
  sprintf(hname,"lngl_%i",i);

  double coef1[] = {201.77, 311.027, 431.582, 553.895, 676.137};
  double coef2[] = {-0.0146608, -0.01677, -0.022244, -0.0271851, -0.0292212};
   
  TF1* lngl = new TF1(hname,"[0]+x*[1]",0,2000);
  lngl->SetParameter(0,coef1[i-1]);
  lngl->SetParameter(1,coef2[i-1]);
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
  t3->SetTextSize(0.02);
  t3->SetTextAngle(-8);
  t3->SetTextColor(kGray+2);


  
  return t3;
}

TLatex* constant_gluino_text(Int_t it,TF1& lngl){
  char legnm[200];

  sprintf(legnm,"#font[12]{#tilde{g}}#font[92]{(%i)GeV}",500+150*(it-1));
  TLatex* t4 = new TLatex(423,18+lngl.Eval(480),legnm);
  t4->SetTextSize(0.02);
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

  std::cout << "contours: " << contours << std::endl;
  // Draw contours
  TList* graphList = (TList*)(contours->At(0));
  std::cout << "number of graphs: " << graphList->GetSize() << std::endl;
  TGraph* myGraph = 0;
  for (int igraph = 0; igraph<graphList->GetSize();++igraph) {
    myGraph = (TGraph*)graphList->At(igraph);

    std::cout << " - graph " << igraph << " has " << myGraph->GetN() << " points" << std::endl;
    if (myGraph->GetN() > 50){
      std::cout << "Drawing " << myGraph->GetN() <<" points" << std::endl;
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
  std::cout << "Nbins: " << nx << " " << ny << std::endl;

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












//old------------------------------------------------------------------------

TGraph* sq_LEP(){//sq-gl
    double sq[] = {0,0,100,100};
    double gl[] = {0,2000,2000,0};   
    TGraph* res = new TGraph(4,sq,gl);
    res->SetFillColor(kBlue);
    return res;
}
TGraph* gl_TEV(){//sq-gl
    double sq[] = {0,2000,2000,0};
    double gl[] = {0,0,190,190};   
    TGraph* res = new TGraph(4,sq,gl);
    res->SetFillColor(kGreen+2);
    return res;
}

TGraph* gl_CDF(){//sq-gl
    double sq[] = {0,2000,2000,0};
    double gl[] = {190,190,230,230};   
    TGraph* res = new TGraph(4,sq,gl);
    res->SetFillColor(kOrange+5);
    return res;
}

TGraph* gl_DEZ(){//sq-gl
    double sq[] = {0,2000,2000,0};
    double gl[] = {230,230,255,255};   
    TGraph* res = new TGraph(4,sq,gl);
    res->SetFillColor(kYellow-5);
    return res;
}

TGraph* gl_WHT(){//sq-gl
    double sq[] = {101,2000,2000,101};
    double gl[] = {256,256,400,400};   
    TGraph* res = new TGraph(4,sq,gl);
    res->SetFillColor(kWhite);
    return res;
}


////////////////////////////////////////////
TGraph* gl_LEP(){//gl-sq
    double sq[] = {0,0,100,100};
    double gl[] = {0,2000,2000,0};
    TGraph* res = new TGraph(4,gl,sq);
    res->SetFillColor(kBlue);
    return res;
}

TGraph* sq_TEV(){//gl-sq
    double sq[] = {0,2000,2000,330,250,300,200,150,100,0};
    double gl[] = {0,0,190,190,260,300,500,560,500,500};

    TGraph* res = new TGraph(10,gl,sq);
    res->SetFillColor(kGreen+2);
    return res;
}

TGraph* sq_CDF(){//gl-sq
    double sq[] = {0,2000,2000,480,460,420,410,380,390,290,0};
    double gl[] = {0,0,280,280,300,310,330,340,440,320,320};
    TGraph* res = new TGraph(11,gl,sq);
    res->SetFillColor(kOrange+5);
    return res;
}
TGraph* sq_DEZ(){//gl-sq
    double sq[] = {0,2000,2000,460,430,400,390,290,0};
    double gl[] = {0,0,305,305,320,350,440,320,320};
    TGraph* res = new TGraph(9,gl,sq);
    res->SetFillColor(kYellow-5);
    return res;
}

TGraph* glsq_NoSol(){//sq-gl
    gStyle->SetHatchesSpacing(2.0);
    gStyle->SetHatchesLineWidth(1);
//    double sq[] = {83,83,110,1297.6,0,   0};
//    double gl[] = { 0,63,120,1466,  1466,0};   
    double sq[] = {83,83,110,1602.3,0,   0};
    double gl[] = { 0,63,120,1800,  1800,0};
    TGraph* res = new TGraph(6,gl,sq);
    res->SetLineColor(1);
    res->SetLineWidth(2);
    res->SetFillStyle(3354);
    return res;
}


TGraph* glsq_NoSol_aux(){//sq-gl : aux for drawing
    gStyle->SetHatchesSpacing(2.0);
    gStyle->SetHatchesLineWidth(1);
    double sq[] = {83,83,110,1602.3};
    double gl[] = { 0,63,120,1800  };
    TGraph* res = new TGraph(4,gl,sq);
    res->SetLineColor(1);
    res->SetLineWidth(2);
    return res;
}


TGraph * Atlas_m0_m12_tb3_obs()
{
   TGraph *graph = new TGraph(71);
   graph->SetLineWidth(2);
   graph->SetPoint(0,54,357.341);
   graph->SetPoint(1,82,358.4171);
   graph->SetPoint(2,110,359.8179);
   graph->SetPoint(3,110.934,359.875);
   graph->SetPoint(4,138,361.5302);
   graph->SetPoint(5,166,362.5696);
   graph->SetPoint(6,189.0676,359.875);
   graph->SetPoint(7,194,359.1415);
   graph->SetPoint(8,222,354.5211);
   graph->SetPoint(9,248.1714,351.625);
   graph->SetPoint(10,250,351.4583);
   graph->SetPoint(11,278,349.4991);
   graph->SetPoint(12,306,343.4051);
   graph->SetPoint(13,306.1313,343.375);
   graph->SetPoint(14,334,335.3855);
   graph->SetPoint(15,334.7966,335.125);
   graph->SetPoint(16,350.8793,326.875);
   graph->SetPoint(17,362,319.7823);
   graph->SetPoint(18,364.6568,318.625);
   graph->SetPoint(19,383.5969,310.375);
   graph->SetPoint(20,390,308.6978);
   graph->SetPoint(21,415.0931,302.125);
   graph->SetPoint(22,418,300.2992);
   graph->SetPoint(23,423.9359,293.875);
   graph->SetPoint(24,431.0348,285.625);
   graph->SetPoint(25,441.6066,277.375);
   graph->SetPoint(26,446,275.3629);
   graph->SetPoint(27,468.0485,269.125);
   graph->SetPoint(28,474,267.2393);
   graph->SetPoint(29,486.3632,260.875);
   graph->SetPoint(30,495.8534,252.625);
   graph->SetPoint(31,502,244.7684);
   graph->SetPoint(32,502.3699,244.375);
   graph->SetPoint(33,508.424,236.125);
   graph->SetPoint(34,514.4781,227.875);
   graph->SetPoint(35,518.8059,219.625);
   graph->SetPoint(36,521.5056,211.375);
   graph->SetPoint(37,524.324,203.125);
   graph->SetPoint(38,527.1424,194.875);
   graph->SetPoint(39,530,193.5052);
   graph->SetPoint(40,552.5782,186.625);
   graph->SetPoint(41,558,185.8767);
   graph->SetPoint(42,586,180.1476);
   graph->SetPoint(43,592.2836,178.375);
   graph->SetPoint(44,614,171.422);
   graph->SetPoint(45,617.5927,170.125);
   graph->SetPoint(46,627.5396,161.875);
   graph->SetPoint(47,642,157.9876);
   graph->SetPoint(48,658.2277,153.625);
   graph->SetPoint(49,670,151.1225);
   graph->SetPoint(50,698,146.9886);
   graph->SetPoint(51,711.5425,145.375);
   graph->SetPoint(52,726,143.6935);
   graph->SetPoint(53,754,140.0648);
   graph->SetPoint(54,782,137.4563);
   graph->SetPoint(55,783.0426,137.125);
   graph->SetPoint(56,796.0326,128.875);
   graph->SetPoint(57,810,125.0328);
   graph->SetPoint(58,826.0235,120.625);
   graph->SetPoint(59,838,117.3305);
   graph->SetPoint(60,866,116.7462);
   graph->SetPoint(61,894,116.1044);
   graph->SetPoint(62,922,114.6645);
   graph->SetPoint(63,950,116.937);
   graph->SetPoint(64,978,117.6956);
   graph->SetPoint(65,1006,115.198);
   graph->SetPoint(66,1034,113.1166);
   graph->SetPoint(67,1038.967,112.375);
   graph->SetPoint(68,1062,108.0997);
   graph->SetPoint(69,1090,104.373);
   graph->SetPoint(70,1096.429,104.125);
   graph->SetFillColor(0);
   return graph;
}
TGraph * Atlas_mGl_mSq_obs()
{
   //x:gluino mass, y:squark mass
   TGraph *graph = new TGraph(64);
   graph->SetLineWidth(2);
   graph->SetPoint(0,511.2569,1976.25);
   graph->SetPoint(1,513.125,1936.181);
   graph->SetPoint(2,513.4714,1928.75);
   graph->SetPoint(3,515.6859,1881.25);
   graph->SetPoint(4,517.9004,1833.75);
   graph->SetPoint(5,520.1148,1786.25);
   graph->SetPoint(6,522.3293,1738.75);
   graph->SetPoint(7,524.5438,1691.25);
   graph->SetPoint(8,526.1095,1643.75);
   graph->SetPoint(9,526.5267,1596.25);
   graph->SetPoint(10,526.9924,1548.75);
   graph->SetPoint(11,527.2745,1501.25);
   graph->SetPoint(12,526.9487,1453.75);
   graph->SetPoint(13,526.6228,1406.25);
   graph->SetPoint(14,532.902,1358.75);
   graph->SetPoint(15,540.182,1311.25);
   graph->SetPoint(16,547.4619,1263.75);
   graph->SetPoint(17,558.0438,1216.25);
   graph->SetPoint(18,561.875,1203.002);
   graph->SetPoint(19,579.4475,1168.75);
   graph->SetPoint(20,603.8166,1121.25);
   graph->SetPoint(21,610.625,1107.979);
   graph->SetPoint(22,630.8311,1073.75);
   graph->SetPoint(23,658.8712,1026.25);
   graph->SetPoint(24,659.375,1025.691);
   graph->SetPoint(25,703.8609,978.75);
   graph->SetPoint(26,708.125,975.5966);
   graph->SetPoint(27,756.875,945.0332);
   graph->SetPoint(28,778.5082,931.25);
   graph->SetPoint(29,805.625,916.3566);
   graph->SetPoint(30,847.5791,883.75);
   graph->SetPoint(31,854.375,879.4212);
   graph->SetPoint(32,903.125,848.0476);
   graph->SetPoint(33,927.7959,836.25);
   graph->SetPoint(34,951.875,823.2943);
   graph->SetPoint(35,1000.625,790.1523);
   graph->SetPoint(36,1004.488,788.75);
   graph->SetPoint(37,1049.375,774.2142);
   graph->SetPoint(38,1098.125,755.1404);
   graph->SetPoint(39,1130.904,741.25);
   graph->SetPoint(40,1146.875,736.2226);
   graph->SetPoint(41,1195.625,720.3208);
   graph->SetPoint(42,1244.375,701.1234);
   graph->SetPoint(43,1266.694,693.75);
   graph->SetPoint(44,1293.125,685.9941);
   graph->SetPoint(45,1341.875,668.9126);
   graph->SetPoint(46,1390.625,652.8072);
   graph->SetPoint(47,1439.375,647.4204);
   graph->SetPoint(48,1448.107,646.25);
   graph->SetPoint(49,1488.125,638.5908);
   graph->SetPoint(50,1536.875,620.8084);
   graph->SetPoint(51,1585.625,603.026);
   graph->SetPoint(52,1597.347,598.75);
   graph->SetPoint(53,1634.375,581.9133);
   graph->SetPoint(54,1683.125,559.7463);
   graph->SetPoint(55,1701.81,551.25);
   graph->SetPoint(56,1731.875,537.5793);
   graph->SetPoint(57,1780.625,515.4123);
   graph->SetPoint(58,1806.273,503.75);
   graph->SetPoint(59,1829.375,493.2101);
   graph->SetPoint(60,1878.125,471.0784);
   graph->SetPoint(61,1910.736,456.25);
   graph->SetPoint(62,1926.875,448.241);
   graph->SetPoint(63,1975.625,426.7444);
   return graph;
}


//- Jim Lungu ----------------------------------------------

TGraph* Jim_mht_tb3(int mode){
  Int_t n = 10;
  Double_t x[10] = {15,65,115,165,215,265,315,365,415,465};
  Double_t yobs[10] = {330,330,340,340,330,320,310,300,280,260};
  Double_t yexp[10] = {330,320,320,320,320,310,300,290,270,250};
  Double_t yexpp[10] = {350,350,350,340,340,330,330,310,300,280};
  Double_t yexpn[10] = {310,310,310,300,300,300,280,270,240,240};
  Double_t y[10];
  for (int i=0;i<n;i++) {
    if (mode==0) y[i] = yobs[i];
    if (mode==1) y[i] = yexp[i];
    if (mode==2) y[i] = yexpp[i];
    if (mode==3) y[i] = yexpn[i];
  }
  TGraph* gr = new TGraph(n,x,y);
  gr->SetLineWidth(2);
  if (mode==1) gr->SetLineStyle(2);
  if (mode>2) {gr->SetFillColor(8); gr->SetLineColor(8);}
  else 
    gr->SetLineColor(2);
  return gr;
}

TGraph* Jim_ht_tb3(int mode){
  Int_t n = 10;
  Double_t x[10] = {15,65,115,165,215,265,315,365,415,465};
  Double_t yobs[10] = {300,300,290,290,300,290,280,270,260,250};
  Double_t yexp[10] = {300,300,290,300,300,290,280,270,260,250};
  Double_t yexpp[10] = {320,320,320,320,310,310,300,290,280,270};
  Double_t yexpn[10] = {290,290,280,290,280,280,270,250,250,240};
  Double_t y[10];
  for (int i=0;i<n;i++) {
    if (mode==0) y[i] = yobs[i];
    if (mode==1) y[i] = yexp[i];
    if (mode==2) y[i] = yexpp[i];
    if (mode==3) y[i] = yexpn[i];
  }
  TGraph* gr = new TGraph(n,x,y);
  gr->SetLineWidth(2);
  if (mode==1) gr->SetLineStyle(2);
  if (mode>2) {gr->SetFillColor(8); gr->SetLineColor(8);}
  else 
    gr->SetLineColor(2);
  return gr;
}

TGraph* Jim_mht_tb50(int mode){
  Int_t n = 9;
  Double_t x[9] = {205,255,305,355,405,455,505,555,605};
  Double_t yobs[9] = {330,330,310,300,300,270,250,240,210};
  Double_t yexp[9] = {320,320,300,290,280,250,240,230,210};
  Double_t yexpp[9] = {350,340,330,320,310,290,260,250,230};
  Double_t yexpn[9] = {300,300,290,280,260,240,230,200,190};
  Double_t y[9];
  for (int i=0;i<n;i++) {
    if (mode==0) y[i] = yobs[i];
    if (mode==1) y[i] = yexp[i];
    if (mode==2) y[i] = yexpp[i];
    if (mode==3) y[i] = yexpn[i];
  }
  TGraph* gr = new TGraph(n,x,y);
  gr->SetLineWidth(2);
  if (mode==1) gr->SetLineStyle(2);
  if (mode>2) {gr->SetFillColor(8); gr->SetLineColor(8);}
  else 
    gr->SetLineColor(2);
  return gr;
}

TGraph* Jim_ht_tb50(int mode){
  Int_t n = 9;
  Double_t x[9] = {205,255,305,355,405,455,505,555,605};
  Double_t yobs[9] = {300,290,280,280,260,250,230,220,220};
  Double_t yexp[9] = {300,290,280,280,260,250,240,230,220};
  Double_t yexpp[9] = {320,310,300,290,290,270,260,250,240};
  Double_t yexpn[9] = {280,270,270,260,250,240,230,220,210};
  Double_t y[10];
  for (int i=0;i<n;i++) {
    if (mode==0) y[i] = yobs[i];
    if (mode==1) y[i] = yexp[i];
    if (mode==2) y[i] = yexpp[i];
    if (mode==3) y[i] = yexpn[i];
  }
  TGraph* gr = new TGraph(n,x,y);
  gr->SetLineWidth(2);
  if (mode==1) gr->SetLineStyle(2);
  if (mode>2) {gr->SetFillColor(8); gr->SetLineColor(8);}
  else 
    gr->SetLineColor(2);
  return gr;
}

TGraph* Jim_mht_tb10(int mode){
  Int_t n = 20;
  Double_t x[20] = {15,65,115,165,215,265,315,365,415,465,515,565,615,665,715,765,815,865,915,965};
  Double_t yobs[20] = {330,340,330,330,330,320,310,280,280,260,240,230,220,220,210,190,190,180,180,180};
  Double_t yexp[20] = {330,320,320,320,320,310,300,280,270,250,240,220,220,200,190,190,180,160,170,160};
  Double_t yexpp[20] = {330,340,350,340,330,330,330,310,300,280,270,240,240,230,220,210,200,180,190,180};
  Double_t yexpn[20] = {290,280,310,300,300,290,280,270,250,240,220,210,200,190,160,170,160,150,150,140};
  Double_t y[20];
  for (int i=0;i<n;i++) {
    if (mode==0) y[i] = yobs[i];
    if (mode==1) y[i] = yexp[i];
    if (mode==2) y[i] = yexpp[i];
    if (mode==3) y[i] = yexpn[i];
  }
  TGraph* gr = new TGraph(n,x,y);
  gr->SetLineWidth(2);
  if (mode==1) gr->SetLineStyle(2);
  if (mode>2) {gr->SetFillColor(8); gr->SetLineColor(8);}
  else 
    gr->SetLineColor(2);
  return gr;
}

TGraph* Jim_ht_tb10(int mode){
  Int_t n = 20;
  Double_t x[20] = {15,65,115,165,215,265,315,365,415,465,515,565,615,665,715,765,815,865,915,965};
  Double_t yobs[20] = {290,280,300,290,300,290,280,270,250,250,240,230,220,220,210,190,200,180,190,180};
  Double_t yexp[20] = {290,280,290,300,300,290,280,270,250,250,240,240,220,220,210,190,200,180,190,180};
  Double_t yexpp[20] = {330,320,310,310,310,310,300,280,280,260,260,250,240,230,220,220,210,210,200,190};
  Double_t yexpn[20] = {290,280,280,280,280,280,270,260,240,240,230,220,200,200,190,190,180,180,170,170};
  Double_t y[20];
  for (int i=0;i<n;i++) {
    if (mode==0) y[i] = yobs[i];
    if (mode==1) y[i] = yexp[i];
    if (mode==2) y[i] = yexpp[i];
    if (mode==3) y[i] = yexpn[i];
  }
  TGraph* gr = new TGraph(n,x,y);
  gr->SetLineWidth(2);
  if (mode==1) gr->SetLineStyle(2);
  if (mode>2) {gr->SetFillColor(8); gr->SetLineColor(8);}
  else 
    gr->SetLineColor(2);
  return gr;
}

TGraphErrors* RA1_NLO( ){

  const unsigned int nNLO=50;

  Double_t xNLO[nNLO],yNLO[nNLO],xerr[nNLO],yerr[nNLO];
  for( unsigned int ierr = 0 ; ierr < nNLO ; ++ierr ){
    xerr[ierr] = 0.;
    yerr[ierr] = 0.;
    xNLO[ierr] = 0.;
    yNLO[ierr] = 0.;
  }

  Double_t xl[12]={  0,100,300,500,700,900,1100,1300,1500,1700,1900,2000};
  Double_t yl[12]={540,535,520,490,400,280, 240, 220, 210, 200, 195, 190};

  TGraphErrors* grtb10  = new TGraphErrors(12,xl, yl,xerr,yerr);
  grtb10->SetMarkerColor(kWhite);
  //grtb10->SetMarkerStyle(21);
  return grtb10;


}


TGraph* Obs2010_NLO(){
  TGraph *graph = new TGraph(130);
  graph->SetName("CLsNLO_ObservedTb10");
  graph->SetTitle("Graph");
  graph->SetFillColor(100);
  graph->SetLineWidth(3);
  graph->SetPoint(0,5.0495,313.8242);
  graph->SetPoint(1,15.1485,314.3105);
  graph->SetPoint(2,25.2475,314.8432);
  graph->SetPoint(3,35.3465,315.4141);
  graph->SetPoint(4,45.4455,315.9472);
  graph->SetPoint(5,55.5445,316.9393);
  graph->SetPoint(6,62.24316,317.2839);
  graph->SetPoint(7,65.6435,317.4589);
  graph->SetPoint(8,75.7425,317.9413);
  graph->SetPoint(9,85.8415,318.3767);
  graph->SetPoint(10,95.9405,318.7587);
  graph->SetPoint(11,106.0395,319.0845);
  graph->SetPoint(12,116.1385,319.2899);
  graph->SetPoint(13,124.4685,319.4044);
  graph->SetPoint(14,133.3826,319.4634);
  graph->SetPoint(15,142.2486,319.4935);
  graph->SetPoint(16,151.0729,319.4177);
  graph->SetPoint(17,159.7999,319.2646);
  graph->SetPoint(18,168.4544,319.0443);
  graph->SetPoint(19,177.0062,318.7368);
  graph->SetPoint(20,185.439,318.3313);
  graph->SetPoint(21,193.7495,317.781);
  graph->SetPoint(22,201.917,317.0872);
  graph->SetPoint(23,209.9348,316.3023);
  graph->SetPoint(24,217.7284,315.4762);
  graph->SetPoint(25,225.267,314.574);
  graph->SetPoint(26,232.5481,313.6103);
  graph->SetPoint(27,239.5453,312.5919);
  graph->SetPoint(28,246.2657,311.4896);
  graph->SetPoint(29,252.7315,310.3486);
  graph->SetPoint(30,258.9505,309.2365);
  graph->SetPoint(31,264.9422,308.0927);
  graph->SetPoint(32,270.7631,306.8801);
  graph->SetPoint(33,276.4743,305.6246);
  graph->SetPoint(34,282.0051,304.389);
  graph->SetPoint(35,287.4176,303.1749);
  graph->SetPoint(36,292.7761,301.9636);
  graph->SetPoint(37,298.1372,300.7121);
  graph->SetPoint(38,303.562,299.3596);
  graph->SetPoint(39,309.0409,297.9735);
  graph->SetPoint(40,314.6325,296.5453);
  graph->SetPoint(41,320.3938,295.0364);
  graph->SetPoint(42,326.3474,293.3707);
  graph->SetPoint(43,332.5077,291.4995);
  graph->SetPoint(44,338.9301,289.4635);
  graph->SetPoint(45,345.5467,287.2513);
  graph->SetPoint(46,352.3343,284.8557);
  graph->SetPoint(47,359.3478,282.2877);
  graph->SetPoint(48,366.5652,279.5413);
  graph->SetPoint(49,373.8549,276.734);
  graph->SetPoint(50,381.2286,273.942);
  graph->SetPoint(51,388.7544,271.0892);
  graph->SetPoint(52,396.4155,268.0935);
  graph->SetPoint(53,404.1744,265.0521);
  graph->SetPoint(54,412.0284,262.1176);
  graph->SetPoint(55,419.9529,259.2966);
  graph->SetPoint(56,427.8657,256.6903);
  graph->SetPoint(57,435.7926,254.343);
  graph->SetPoint(58,443.732,252.1696);
  graph->SetPoint(59,451.6475,250.1392);
  graph->SetPoint(60,459.522,248.209);
  graph->SetPoint(61,467.4147,246.3386);
  graph->SetPoint(62,475.3556,244.5038);
  graph->SetPoint(63,483.2743,242.6911);
  graph->SetPoint(64,491.1437,240.9127);
  graph->SetPoint(65,498.9059,239.1429);
  graph->SetPoint(66,506.6303,237.878);
  graph->SetPoint(67,514.3106,237.0588);
  graph->SetPoint(68,521.9283,236.1634);
  graph->SetPoint(69,529.4356,235.1685);
  graph->SetPoint(70,536.8647,234.0413);
  graph->SetPoint(71,544.2577,232.834);
  graph->SetPoint(72,551.6591,231.503);
  graph->SetPoint(73,559.1076,230.0904);
  graph->SetPoint(74,566.5661,228.6313);
  graph->SetPoint(75,574.0618,227.1308);
  graph->SetPoint(76,581.6695,225.6218);
  graph->SetPoint(77,589.3114,224.1005);
  graph->SetPoint(78,596.9926,222.4802);
  graph->SetPoint(79,604.7216,220.8893);
  graph->SetPoint(80,612.4689,219.333);
  graph->SetPoint(81,620.2864,217.8096);
  graph->SetPoint(82,628.1088,216.3082);
  graph->SetPoint(83,635.9936,214.8319);
  graph->SetPoint(84,643.9187,213.4581);
  graph->SetPoint(85,651.8702,212.0884);
  graph->SetPoint(86,659.8609,210.7696);
  graph->SetPoint(87,667.8767,209.5601);
  graph->SetPoint(88,675.8856,208.4753);
  graph->SetPoint(89,683.8681,207.4519);
  graph->SetPoint(90,691.8409,206.4889);
  graph->SetPoint(91,699.8237,205.5841);
  graph->SetPoint(92,707.8874,204.7591);
  graph->SetPoint(93,715.9544,203.9953);
  graph->SetPoint(94,724.0324,203.1892);
  graph->SetPoint(95,732.168,202.3859);
  graph->SetPoint(96,740.4279,201.6071);
  graph->SetPoint(97,748.7894,200.8673);
  graph->SetPoint(98,757.2338,200.1053);
  graph->SetPoint(99,765.7054,199.2916);
  graph->SetPoint(100,774.2199,198.3878);
  graph->SetPoint(101,782.7635,197.5044);
  graph->SetPoint(102,791.2852,196.6289);
  graph->SetPoint(103,799.8181,195.7428);
  graph->SetPoint(104,808.3957,194.8782);
  graph->SetPoint(105,816.9027,194.0288);
  graph->SetPoint(106,825.3378,193.0564);
  graph->SetPoint(107,833.736,192.0561);
  graph->SetPoint(108,842.0355,191.0663);
  graph->SetPoint(109,850.2923,190.1167);
  graph->SetPoint(110,858.4818,189.2294);
  graph->SetPoint(111,866.5929,188.3491);
  graph->SetPoint(112,874.6509,187.4837);
  graph->SetPoint(113,882.6188,186.6064);
  graph->SetPoint(114,890.538,185.68);
  graph->SetPoint(115,898.4241,184.8146);
  graph->SetPoint(116,906.258,183.9957);
  graph->SetPoint(117,913.9595,183.2018);
  graph->SetPoint(118,917.3598,182.4365);
  graph->SetPoint(119,924.0585,181.7045);
  graph->SetPoint(120,934.1575,181.0356);
  graph->SetPoint(121,944.2565,180.4533);
  graph->SetPoint(122,954.3555,179.9049);
  graph->SetPoint(123,964.4545,179.4012);
  graph->SetPoint(124,974.5535,178.9513);
  graph->SetPoint(125,984.6525,178.562);
  graph->SetPoint(126,988.0528,178.2617);
  graph->SetPoint(127,994.7515,178.0687);
  graph->SetPoint(128,1001.45,177.874);
  graph->SetPoint(129,1004.851,177.7141);

  return graph;
}

TGraph* Exp2010_NLO(){
  TGraph *graph = new TGraph(137);
  graph->SetName("CLsNLO_ExpectedTb10");
  graph->SetTitle("Graph");
  graph->SetFillColor(100);
  graph->SetLineStyle(2);
  graph->SetLineWidth(3);
  graph->SetPoint(0,5.0495,300.5129);
  graph->SetPoint(1,15.1485,300.9677);
  graph->SetPoint(2,21.84716,301.2948);
  graph->SetPoint(3,25.2475,301.4608);
  graph->SetPoint(4,35.3465,301.9853);
  graph->SetPoint(5,45.4455,303.0626);
  graph->SetPoint(6,55.5445,304.0056);
  graph->SetPoint(7,62.24316,304.4887);
  graph->SetPoint(8,65.6435,304.9372);
  graph->SetPoint(9,75.7425,305.3447);
  graph->SetPoint(10,85.8415,305.7067);
  graph->SetPoint(11,89.24184,305.8128);
  graph->SetPoint(12,95.9405,306.0217);
  graph->SetPoint(13,97.59791,306.0627);
  graph->SetPoint(14,103.888,306.2183);
  graph->SetPoint(15,110.1274,306.3697);
  graph->SetPoint(16,116.3795,306.4989);
  graph->SetPoint(17,122.6683,306.562);
  graph->SetPoint(18,128.9313,306.5791);
  graph->SetPoint(19,135.1892,306.5627);
  graph->SetPoint(20,141.487,306.5009);
  graph->SetPoint(21,147.8971,306.3795);
  graph->SetPoint(22,154.4551,306.2034);
  graph->SetPoint(23,161.1149,306.0154);
  graph->SetPoint(24,167.8912,305.8059);
  graph->SetPoint(25,174.7781,305.5289);
  graph->SetPoint(26,181.7589,305.133);
  graph->SetPoint(27,188.7904,304.663);
  graph->SetPoint(28,195.8412,304.1339);
  graph->SetPoint(29,202.9348,303.551);
  graph->SetPoint(30,210.0235,302.9191);
  graph->SetPoint(31,217.0261,302.2364);
  graph->SetPoint(32,223.9634,301.4868);
  graph->SetPoint(33,230.8665,300.6491);
  graph->SetPoint(34,237.6992,299.7569);
  graph->SetPoint(35,244.462,298.8047);
  graph->SetPoint(36,251.1199,297.7871);
  graph->SetPoint(37,257.6422,296.7234);
  graph->SetPoint(38,264.0752,295.6414);
  graph->SetPoint(39,270.4454,294.5133);
  graph->SetPoint(40,276.7675,293.2863);
  graph->SetPoint(41,283.0517,291.963);
  graph->SetPoint(42,289.3272,290.5521);
  graph->SetPoint(43,295.6325,288.9925);
  graph->SetPoint(44,301.9649,287.3513);
  graph->SetPoint(45,308.3892,285.5887);
  graph->SetPoint(46,314.9451,283.6627);
  graph->SetPoint(47,321.6225,281.5525);
  graph->SetPoint(48,328.4148,279.2871);
  graph->SetPoint(49,335.3427,276.8879);
  graph->SetPoint(50,342.3888,274.3826);
  graph->SetPoint(51,349.5217,271.7413);
  graph->SetPoint(52,356.8404,268.9502);
  graph->SetPoint(53,364.242,266.086);
  graph->SetPoint(54,371.7563,263.1322);
  graph->SetPoint(55,379.4145,260.0843);
  graph->SetPoint(56,387.146,257.1272);
  graph->SetPoint(57,394.9464,254.3058);
  graph->SetPoint(58,402.8615,251.5706);
  graph->SetPoint(59,410.7879,248.9807);
  graph->SetPoint(60,418.6862,246.5155);
  graph->SetPoint(61,426.5895,244.1677);
  graph->SetPoint(62,434.5064,242.8891);
  graph->SetPoint(63,442.4945,241.9938);
  graph->SetPoint(64,450.5271,241.1735);
  graph->SetPoint(65,458.5219,240.3202);
  graph->SetPoint(66,466.5067,239.3827);
  graph->SetPoint(67,474.5647,238.4213);
  graph->SetPoint(68,482.6686,237.4448);
  graph->SetPoint(69,490.7993,236.3902);
  graph->SetPoint(70,498.9696,235.2629);
  graph->SetPoint(71,507.1744,234.0642);
  graph->SetPoint(72,515.3984,232.7821);
  graph->SetPoint(73,523.6416,231.4153);
  graph->SetPoint(74,531.9107,229.9799);
  graph->SetPoint(75,540.2113,228.4931);
  graph->SetPoint(76,548.5486,226.9718);
  graph->SetPoint(77,556.9526,225.3851);
  graph->SetPoint(78,565.446,223.7861);
  graph->SetPoint(79,573.9101,222.1862);
  graph->SetPoint(80,582.4165,220.579);
  graph->SetPoint(81,590.958,219.035);
  graph->SetPoint(82,599.5117,217.4998);
  graph->SetPoint(83,608.0864,215.9827);
  graph->SetPoint(84,616.6249,214.4591);
  graph->SetPoint(85,625.1458,212.8941);
  graph->SetPoint(86,633.7226,211.3607);
  graph->SetPoint(87,642.2853,209.9407);
  graph->SetPoint(88,650.7635,208.563);
  graph->SetPoint(89,659.1396,207.1816);
  graph->SetPoint(90,667.4248,205.7657);
  graph->SetPoint(91,675.6389,204.3693);
  graph->SetPoint(92,683.8275,203.0231);
  graph->SetPoint(93,691.9539,201.7508);
  graph->SetPoint(94,699.9947,200.4999);
  graph->SetPoint(95,707.9632,199.2818);
  graph->SetPoint(96,715.8771,198.1068);
  graph->SetPoint(97,723.7582,196.9841);
  graph->SetPoint(98,731.6149,195.9127);
  graph->SetPoint(99,739.4497,194.8737);
  graph->SetPoint(100,747.2823,193.8541);
  graph->SetPoint(101,755.1398,192.8549);
  graph->SetPoint(102,763.0472,191.8751);
  graph->SetPoint(103,771.0268,190.9126);
  graph->SetPoint(104,779.0968,189.9647);
  graph->SetPoint(105,787.2718,189.0291);
  graph->SetPoint(106,795.6106,188.1287);
  graph->SetPoint(107,804.0175,187.248);
  graph->SetPoint(108,812.4795,186.2878);
  graph->SetPoint(109,820.9869,185.3628);
  graph->SetPoint(110,829.4936,184.4651);
  graph->SetPoint(111,838.0852,183.6049);
  graph->SetPoint(112,846.6475,182.7776);
  graph->SetPoint(113,855.1522,181.893);
  graph->SetPoint(114,863.6039,181.0256);
  graph->SetPoint(115,871.9675,180.2085);
  graph->SetPoint(116,880.2025,179.4492);
  graph->SetPoint(117,888.2795,178.7289);
  graph->SetPoint(118,896.1251,178.0026);
  graph->SetPoint(119,903.7555,177.23);
  graph->SetPoint(120,911.1671,176.5176);
  graph->SetPoint(121,918.3289,175.8459);
  graph->SetPoint(122,925.2691,175.252);
  graph->SetPoint(123,931.9649,174.7216);
  graph->SetPoint(124,944.2565,174.2344);
  graph->SetPoint(125,947.6568,173.8205);
  graph->SetPoint(126,954.3555,173.5037);
  graph->SetPoint(127,961.0542,173.233);
  graph->SetPoint(128,964.4545,173.0213);
  graph->SetPoint(129,967.8548,172.8802);
  graph->SetPoint(130,974.5535,172.8193);
  graph->SetPoint(131,977.9538,172.8445);
  graph->SetPoint(132,984.6525,172.9579);
  graph->SetPoint(133,991.3512,173.1569);
  graph->SetPoint(134,994.7515,173.4587);
  graph->SetPoint(135,1001.45,173.8689);
  graph->SetPoint(136,1004.851,174.2656);

  return graph;
}
