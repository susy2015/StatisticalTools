#ifndef SUSYSCAN_H
#define SUSYSCAN_H

#include <string>
#include <vector>

class SusyScan{
 public:
  SusyScan();
  SusyScan(const SusyScan&);
  SusyScan(std::string file);
  
 private:	
  std::vector<double*> p;
 public:
  const SusyScan operator*(const double f) const;
  const SusyScan operator+(const SusyScan& f) const;
  
  double Mzero;
  double Mhalf;
  double Mu;
  double TanBeta;
  double Azero;
  double Run;
  double data;
  double Xsection;
  double background;
  double background_uncertainty;
  double signal;
  double signal_uncertainty;
  double signal_contamination;
  double signal_kfactor;
  double ExpXsecLimit;
  double ObsXsecLimit;
  double ObsXsecBayes;
  double ExpXsecBayes;
  double ExpXsecLimit_n1;
  double ExpXsecLimit_p1;
  double ExpXsecLimit_n2;
  double ExpXsecLimit_p2;
  double ExpXsecBayes_n1;
  double ExpXsecBayes_p1;
  double ExpXsecBayes_n2;
  double ExpXsecBayes_p2;
  double CLb_b_xsec;
  double CLs_xsec;
  double CLs_b_xsec;
  double CLs_b_n1_xsec;
  double CLs_b_n2_xsec;
  double CLs_b_p1_xsec;
  double CLs_b_p2_xsec;
  double CLsb_b_xsec;
  
  double NLO_signal;
  double NLO_signal_uncertainty;
  double NLO_signal_contamination;
  double NLO_ExpXsecLimit;
  double NLO_ObsXsecLimit;
  double NLO_ObsXsecBayes;
  double NLO_ExpXsecBayes;
  double NLO_ExpXsecLimit_n1;
  double NLO_ExpXsecLimit_p1;
  double NLO_ExpXsecLimit_n2;
  double NLO_ExpXsecLimit_p2;
  double NLO_ExpXsecBayes_n1;
  double NLO_ExpXsecBayes_p1;
  double NLO_ExpXsecBayes_n2;
  double NLO_ExpXsecBayes_p2;
  double NLO_CLb_b_xsec;
  double NLO_CLs_xsec;
  double NLO_CLs_b_xsec;
  double NLO_CLs_error_exp;
  double NLO_CLs_b_n1_xsec;
  double NLO_CLs_b_n2_xsec;
  double NLO_CLs_b_p1_xsec;
  double NLO_CLs_b_p2_xsec;
  double NLO_CLsb_b_xsec;

  double M1;
  double M2;
  double M3;
  double MGL;
  double MUL;
  double MB1;
  double MSN;
  double MNTAU;
  double MZ1;
  double MW1;
  double MHL;
  double MUR;
  double MB2;
  double MEL;
  double MTAU1;
  double MZ2;
  double MW2;
  double MHH;
  double MDL;
  double MT1;
  double MER;
  double MTAU2;
  double MZ3;
  double MHA;
  double MDR;
  double MT2;
  double MZ4;
  double MHp;

  double Mgluino;
  double Msquark;
  double Mlsp;

 private:
  void SetPtr();

};


#endif
