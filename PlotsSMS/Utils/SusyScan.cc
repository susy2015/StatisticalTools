#include "SusyScan.h"
#include "ConfigFile.h"


SusyScan::SusyScan()
{
  SetPtr();
  for (std::vector<double*>::iterator it=p.begin(); it!=p.end(); ++it)
    **it = 0.;
}
SusyScan::SusyScan(const SusyScan& c)
{
  SetPtr();
  //std::cout<<Mzero<<"=m0, "<<Mhalf <<"=m12, " <<p.size()<<"<-this, c:"<<c.p.size()<<std::endl; 

  std::vector<double*>::iterator it=p.begin();
  std::vector<double*>::const_iterator ct=c.p.begin();
  for (; ct!=c.p.end(); ++it, ++ct)
    **it = **ct;
}


SusyScan::SusyScan(std::string filename)
{
  SetPtr();
  ConfigFile config(filename);
  Mzero 	= config.read<double>("Mzero", 0); 
  Mhalf	 	= config.read<double>("Mhalf", 0); 
  Mu 		= config.read<double>("Mu", 0);
  TanBeta 	= config.read<double>("TanBeta", 0);
  Azero 	= config.read<double>("Azero", 0);
  Run 		= config.read<double>("Run", 0);
  background 	= config.read<double>("background", 0);
  background_uncertainty = config.read<double>("background.uncertainty", 0);
  data 		= config.read<double>("data", 0);
  signal_contamination 	= config.read<double>("signal.LO.signalregion.IsoMuon", 0) + 
                          config.read<double>("signal.LO.signalregion.Tau", 0);
  signal 	= config.read<double>("signal.LO", 0);
  signal_kfactor 	= config.read<double>("signal.kFactor", 1);
  signal_uncertainty = config.read<double>("signal.LO.uncertainty", 0);
  Xsection 	= config.read<double>("Xsection", 0);
  ExpXsecLimit 	= config.read<double>("LO.CLs.expected", 0);
  ObsXsecLimit 	= config.read<double>("LO.CLs.observed", 0);
  ObsXsecBayes 	= config.read<double>("LO.Bayes.observed", 0);
  ExpXsecBayes 	= config.read<double>("LO.Bayes.expected", 0);
  ExpXsecLimit_n1   = config.read<double>("LO.CLs.expected.n1", 0);
  ExpXsecLimit_p1   = config.read<double>("LO.CLs.expected.p1", 0);
  ExpXsecLimit_n2   = config.read<double>("LO.CLs.expected.n2", 0);
  ExpXsecLimit_p2   = config.read<double>("LO.CLs.expected.p2", 0);
  ExpXsecBayes_n1   = config.read<double>("LO.Bayes.expected.n1", 0);
  ExpXsecBayes_p1   = config.read<double>("LO.Bayes.expected.p1", 0);
  ExpXsecBayes_n2   = config.read<double>("LO.Bayes.expected.n2", 0);
  ExpXsecBayes_p2   = config.read<double>("LO.Bayes.expected.p2", 0);
  CLb_b_xsec 	= config.read<double>("LO.CLb_b@xsec", 0);
  CLs_xsec 	= config.read<double>("LO.CLs@xsec", 1);
  CLs_b_xsec 	= config.read<double>("LO.CLsExp", 1);
  CLs_b_n1_xsec = config.read<double>("LO.CLs_b_n1@xsec", 0);
  CLs_b_n2_xsec	= config.read<double>("LO.CLs_b_n2@xsec", 0);
  CLs_b_p1_xsec	= config.read<double>("LO.CLs_b_p1@xsec", 0);
  CLs_b_p2_xsec	= config.read<double>("LO.CLs_b_p2@xsec", 0);
  CLsb_b_xsec	= config.read<double>("LO.CLsb_b@xsec", 0);

  NLO_signal_contamination 	= config.read<double>("signal.NLO.signalregion.IsoMuon", 0) +
                                  config.read<double>("signal.NLO.signalregion.Tau", 0);
  NLO_signal 	        = config.read<double>("signal.NLO", 0);
  NLO_signal_uncertainty= config.read<double>("signal.NLO.uncertainty", 0);
  NLO_ExpXsecLimit 	= config.read<double>("NLO.CLs.expected", 0);
  NLO_ObsXsecLimit 	= config.read<double>("NLO.CLs.observed", 0);
  NLO_ObsXsecBayes 	= config.read<double>("NLO.Bayes.observed", 0);
  NLO_ExpXsecBayes 	= config.read<double>("NLO.Bayes.expected", 0);
  NLO_ExpXsecLimit_n1	= config.read<double>("NLO.CLs.expected.n1", 0);
  NLO_ExpXsecLimit_p1	= config.read<double>("NLO.CLs.expected.p1", 0);
  NLO_ExpXsecLimit_n2	= config.read<double>("NLO.CLs.expected.n2", 0);
  NLO_ExpXsecLimit_p2	= config.read<double>("NLO.CLs.expected.p2", 0);
  NLO_ExpXsecBayes_n1	= config.read<double>("NLO.Bayes.expected.n1", 0);
  NLO_ExpXsecBayes_p1	= config.read<double>("NLO.Bayes.expected.p1", 0);
  NLO_ExpXsecBayes_n2	= config.read<double>("NLO.Bayes.expected.n2", 0);
  NLO_ExpXsecBayes_p2	= config.read<double>("NLO.Bayes.expected.p2", 0);  
  NLO_CLb_b_xsec 	= config.read<double>("NLO.CLb_b@xsec", 1);
  NLO_CLs_xsec 	        = config.read<double>("limit.cls.observed", 1);
  NLO_CLs_b_xsec 	= config.read<double>("limit.cls.expected", 0);
  
  //CLs values need to call different names for RA3 script.
  //NLO_CLs_error_exp     = config.read<double>("NLO.CLs_errorExp",0);
  NLO_CLs_b_n1_xsec     = config.read<double>("limit.cls.expected.m1sigma");
  NLO_CLs_b_n2_xsec	= config.read<double>("limit.cls.expected.m2sigma", 0);
  NLO_CLs_b_p1_xsec	= config.read<double>("limit.cls.expected.p1sigma");
  NLO_CLs_b_p2_xsec	= config.read<double>("limit.cls.expected.p2sigma", 0);
  //printf("%8.3e %8.3e %8.3e %8.3e\n",NLO_CLs_b_xsec,NLO_CLs_error_exp,NLO_CLs_b_n1_xsec,NLO_CLs_b_p1_xsec);
 
  NLO_CLsb_b_xsec	= config.read<double>("NLO.CLsb_b@xsec", 0);
}

void SusyScan::SetPtr()
{
  if (p.size()!=0) return; 
  p.push_back( &Mzero);
  p.push_back( &Mhalf);
  p.push_back( &Mu);
  p.push_back( &TanBeta);
  p.push_back( &Azero);
  p.push_back( &Run);
  p.push_back( &data);
  p.push_back( &Xsection);
  p.push_back( &background);
  p.push_back( &background_uncertainty);
  p.push_back( &signal_kfactor);
  p.push_back( &signal);
  p.push_back( &signal_uncertainty);
  p.push_back( &signal_contamination);
  p.push_back( &ExpXsecLimit);
  p.push_back( &ExpXsecBayes);
  p.push_back( &ObsXsecLimit);
  p.push_back( &ObsXsecBayes);
  p.push_back( &ExpXsecLimit_n1);
  p.push_back( &ExpXsecLimit_p1);
  p.push_back( &ExpXsecLimit_n2);
  p.push_back( &ExpXsecLimit_p2);
  p.push_back( &ExpXsecBayes_n1);
  p.push_back( &ExpXsecBayes_p1);
  p.push_back( &ExpXsecBayes_n2);
  p.push_back( &ExpXsecBayes_p2);
  p.push_back( &CLb_b_xsec);
  p.push_back( &CLs_xsec);
  p.push_back( &CLs_b_xsec);
  p.push_back( &CLs_b_n1_xsec);
  p.push_back( &CLs_b_n2_xsec);
  p.push_back( &CLs_b_p1_xsec);
  p.push_back( &CLs_b_p2_xsec);
  p.push_back( &CLsb_b_xsec);
  p.push_back( &NLO_signal);
  p.push_back( &NLO_signal_uncertainty);
  p.push_back( &NLO_signal_contamination);
  p.push_back( &NLO_ExpXsecLimit);
  p.push_back( &NLO_ExpXsecBayes);
  p.push_back( &NLO_ObsXsecLimit);
  p.push_back( &NLO_ObsXsecBayes);
  p.push_back( &NLO_ExpXsecLimit_n1);
  p.push_back( &NLO_ExpXsecLimit_p1);
  p.push_back( &NLO_ExpXsecLimit_n2);
  p.push_back( &NLO_ExpXsecLimit_p2);
  p.push_back( &NLO_ExpXsecBayes_n1);
  p.push_back( &NLO_ExpXsecBayes_p1);
  p.push_back( &NLO_ExpXsecBayes_n2);
  p.push_back( &NLO_ExpXsecBayes_p2);
  p.push_back( &NLO_CLb_b_xsec);
  p.push_back( &NLO_CLs_xsec);
  p.push_back( &NLO_CLs_b_xsec);
  p.push_back( &NLO_CLs_b_n1_xsec);
  p.push_back( &NLO_CLs_b_n2_xsec);
  p.push_back( &NLO_CLs_b_p1_xsec);
  p.push_back( &NLO_CLs_b_p2_xsec);
  p.push_back( &NLO_CLsb_b_xsec);
  p.push_back( &M1);
  p.push_back( &M2);
  p.push_back( &M3);
  p.push_back( &MGL);
  p.push_back( &MUL);
  p.push_back( &MB1);
  p.push_back( &MSN);
  p.push_back( &MNTAU);
  p.push_back( &MZ1);
  p.push_back( &MW1);
  p.push_back( &MHL);
  p.push_back( &MUR);
  p.push_back( &MB2);
  p.push_back( &MEL);
  p.push_back( &MTAU1);
  p.push_back( &MZ2);
  p.push_back( &MW2);
  p.push_back( &MHH);
  p.push_back( &MDL);
  p.push_back( &MT1);
  p.push_back( &MER);
  p.push_back( &MTAU2);
  p.push_back( &MZ3);
  p.push_back( &MHA);
  p.push_back( &MDR);
  p.push_back( &MT2);
  p.push_back( &MZ4);
  p.push_back( &MHp);

}

const SusyScan SusyScan::operator*(const double f) const
{
  SusyScan res(*this);
  for (std::vector<double*>::iterator it=res.p.begin(); it!=res.p.end(); ++it)
    **it *= f;
  return res;
}

const SusyScan SusyScan::operator+(const SusyScan& f) const
{
  SusyScan res(*this);
  std::vector<double*>::iterator it=res.p.begin(); 
  std::vector<double*>::const_iterator fi=  f.p.begin(); 
  for (; it!=res.p.end(); ++it, ++fi)
    **it += **fi;
  return res;
}
