#!/usr/bin/python

#from os import listdir
#from os import sys
#from os import system
#from os import path
import os
import math
from optparse import OptionParser # Command line parsing

import ROOT as rt

glb_data_rate = []
glb_lostle_rate = []
glb_hadtau_rate = []
glb_zinv_rate = []
glb_qcd_rate = []
glb_ttz_rate = []
glb_rare_rate = []
glb_tot_rate = []

glb_lostle_stat_up = []
glb_hadtau_stat_up = []
glb_zinv_stat_up = []
glb_qcd_stat_up = []
glb_ttz_stat_up = []
glb_rare_stat_up = []
glb_lostle_stat_dn = []
glb_hadtau_stat_dn = []
glb_zinv_stat_dn = []
glb_qcd_stat_dn = []
glb_ttz_stat_dn = []
glb_rare_stat_dn = []

glb_tot_stat_up = []
glb_tot_stat_dn = []

glb_lostle_syst_up = []
glb_hadtau_syst_up = []
glb_zinv_syst_up = []
glb_qcd_syst_up = []
glb_ttz_syst_up = []
glb_rare_syst_up = []
glb_lostle_syst_dn = []
glb_hadtau_syst_dn = []
glb_zinv_syst_dn = []
glb_qcd_syst_dn = []
glb_ttz_syst_dn = []
glb_rare_syst_dn = []

glb_tot_syst_up = []
glb_tot_syst_dn = []

glb_signal_rate = []
glb_signal_stat_up = []
glb_signal_stat_dn = []
glb_signal_syst_up = []
glb_signal_syst_dn = []

# Remove the "=" part and anything after "#"
def procline(inputline):
   outline = []
   tmpline = inputline.split()
   isComments = False
   for each in tmpline:
      if each == "#":
         isComments = True
      if each == "=": 
         continue
      if isComments == False :
         outline.append(each)
   return outline

def prodCardPerChn(signal_key, outputdir="", lostle_file ="lostle.txt", hadtau_file ="hadtau.txt", zinv_file="zinv.txt", qcd_file="qcd.txt", ttz_file="qcd.txt", rare_file="rare.txt", data_file="data.txt", signal_file ="signal.txt"):

   del glb_data_rate[:]
   del glb_lostle_rate[:]
   del glb_hadtau_rate[:]
   del glb_zinv_rate[:]
   del glb_qcd_rate[:]
   del glb_ttz_rate[:]
   del glb_rare_rate[:]
   del glb_tot_rate[:]
#   
   del glb_lostle_stat_up[:]
   del glb_hadtau_stat_up[:]
   del glb_zinv_stat_up[:]
   del glb_qcd_stat_up[:]
   del glb_ttz_stat_up[:]
   del glb_rare_stat_up[:]
   del glb_lostle_stat_dn[:]
   del glb_hadtau_stat_dn[:]
   del glb_zinv_stat_dn[:]
   del glb_qcd_stat_dn[:]
   del glb_ttz_stat_dn[:]
   del glb_rare_stat_dn[:]
#   
   del glb_tot_stat_up[:]
   del glb_tot_stat_dn[:]
#   
   del glb_lostle_syst_up[:]
   del glb_hadtau_syst_up[:]
   del glb_zinv_syst_up[:]
   del glb_qcd_syst_up[:]
   del glb_ttz_syst_up[:]
   del glb_rare_syst_up[:]
   del glb_lostle_syst_dn[:]
   del glb_hadtau_syst_dn[:]
   del glb_zinv_syst_dn[:]
   del glb_qcd_syst_dn[:]
   del glb_ttz_syst_dn[:]
   del glb_rare_syst_dn[:]
#   
   del glb_tot_syst_up[:]
   del glb_tot_syst_dn[:]
#   
   del glb_signal_rate[:]
   del glb_signal_stat_up[:]
   del glb_signal_stat_dn[:]
   del glb_signal_syst_up[:]
   del glb_signal_syst_dn[:]

   doExpLimitOnly = False

   outbase_filename = "comb_" + signal_key
   print '\noutbase_filename : %s\n' % (outbase_filename)

   prt_table_file = open("prt_table.txt", "w")

   # parsing data file first to get basic information
   global channels
   for line in data_file:
      splitline = procline(line)
      if splitline: 
#         print splitline
         if splitline[0] == "luminosity" : lumi = splitline[1]
         if splitline[0] == "channels" : channels = splitline[1]

   val_channels = int(channels)

   rtfile = rt.TFile("searchBins.root", "RECREATE") 
   h1_data = rt.TH1D("data", "data", val_channels, 0, val_channels)
   h1_data.Sumw2()
#   h1_data.SetBinErrorOption(rt.TH1.kPoisson)
   h1_lostle = rt.TH1D("lostle", "lostle", val_channels, 0, val_channels)
   h1_lostle.Sumw2()
   h1_hadtau = rt.TH1D("hadtau", "hadtau", val_channels, 0, val_channels)
   h1_hadtau.Sumw2()
   h1_zinv = rt.TH1D("zinv", "zinv", val_channels, 0, val_channels)
   h1_zinv.Sumw2()
   h1_qcd = rt.TH1D("qcd", "qcd", val_channels, 0, val_channels)
   h1_qcd.Sumw2()
   h1_ttz = rt.TH1D("ttz", "ttz", val_channels, 0, val_channels)
   h1_ttz.Sumw2()
   h1_rare = rt.TH1D("rare", "rare", val_channels, 0, val_channels)
   h1_rare.Sumw2()

   h1_lostle_syst = rt.TH1D("lostle_syst", "lostle_syst", val_channels, 0, val_channels)
   h1_lostle_syst.Sumw2()
   h1_hadtau_syst = rt.TH1D("hadtau_syst", "hadtau_syst", val_channels, 0, val_channels)
   h1_hadtau_syst.Sumw2()
   h1_zinv_syst = rt.TH1D("zinv_syst", "zinv_syst", val_channels, 0, val_channels)
   h1_zinv_syst.Sumw2()
   h1_qcd_syst = rt.TH1D("qcd_syst", "qcd_syst", val_channels, 0, val_channels)
   h1_qcd_syst.Sumw2()
   h1_ttz_syst = rt.TH1D("ttz_syst", "ttz_syst", val_channels, 0, val_channels)
   h1_ttz_syst.Sumw2()
   h1_rare_syst = rt.TH1D("rare_syst", "rare_syst", val_channels, 0, val_channels)
   h1_rare_syst.Sumw2()

   qcd_tfactor_cnt = 1
   qcd_tfactor_dict = {}

   for chn in range(1, int(channels)+1):
      if len(outputdir) !=0:
         outfile_perchn = open(outputdir + "/" + outbase_filename + "_ch" + str(chn) + ".txt", "w")
      else:
         outfile_perchn = open(outbase_filename + "_ch" + str(chn) + ".txt", "w")
      outfile_perchn.write("imax 1 # number of channels\n")
      outfile_perchn.write("jmax 8 # number of backgrounds\n")
      outfile_perchn.write("kmax * nuissance\n")
      outfile_perchn.write("shapes * * FAKE\n")
      outfile_perchn.write("----------------\n")
      outfile_perchn.write("bin bin%d\n" % (chn))

      data_file.seek(0, 0)
      for data_line in data_file:
         data_splitline = procline(data_line)
         if data_splitline and data_splitline[0] == "rate": 
            data_rate = data_splitline[chn]
      outfile_perchn.write("observation %s\n" % (data_rate))

      h1_data.SetBinContent(chn, float(data_rate))
      h1_data.SetBinError(chn, math.sqrt(float(data_rate)))

      outfile_perchn.write("bin       ")
      for ibkg in range(9):
         outfile_perchn.write("bin%d "% (chn))
      outfile_perchn.write("\n")

      outfile_perchn.write("process   Sig  LostLep  LostLepHighW  HadTau  HadTauHighW  Zinv  QCD  TTZ  RARE\n")
      outfile_perchn.write("process     0        1            2       3           4     5    6    7    8\n")
      outfile_perchn.write("rate      ")

      signal_file.seek(0, 0)
      for signal_line in signal_file:
         signal_splitline = procline(signal_line)
         if signal_splitline and signal_splitline[0] == "rate":
            signal_rate = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "cs_event":
            signal_cs_event = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "avg_weight":
            signal_avg_weight = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "contam":
            signal_contam = float(signal_splitline[chn])
      outfile_perchn.write("%.4f  " % ((signal_cs_event-signal_contam)*signal_avg_weight) )
      signal_corr_rate = (signal_cs_event-signal_contam)*signal_avg_weight
      glb_signal_rate.append(signal_corr_rate)

      pred_tot_rate = 0
      pred_tot_stat = 0
 
      lostle_file.seek(0, 0)
      for lostle_line in lostle_file:
         lostle_splitline = procline(lostle_line)
         if lostle_splitline and lostle_splitline[0] == "rate":
            lostle_rate = float(lostle_splitline[chn])
            if lostle_rate == 0: lostle_rate = 0.0001
         if lostle_splitline and lostle_splitline[0] == "avg_weight":
            lostle_avg_weight = float(lostle_splitline[chn])
      outfile_perchn.write("%.4f  " % (lostle_rate))
#      outfile_perchn.write("0.0001  ")
      outfile_perchn.write("%.4f  " % (lostle_avg_weight))
      pred_tot_rate += lostle_rate

      h1_lostle.SetBinContent(chn, lostle_rate)
      h1_lostle_syst.SetBinContent(chn, lostle_rate)
 
      hadtau_file.seek(0, 0)
      for hadtau_line in hadtau_file:
         hadtau_splitline = procline(hadtau_line)
         if hadtau_splitline and hadtau_splitline[0] == "rate":
            hadtau_rate = float(hadtau_splitline[chn])
            if hadtau_rate ==0: hadtau_rate = 0.0001
      outfile_perchn.write("%.4f  " % (hadtau_rate))
      outfile_perchn.write("0.3  ")
      pred_tot_rate += hadtau_rate
 
      h1_hadtau.SetBinContent(chn, hadtau_rate)
      h1_hadtau_syst.SetBinContent(chn, hadtau_rate)
 
      zinv_file.seek(0, 0)
      for zinv_line in zinv_file:
         zinv_splitline = procline(zinv_line)
         if zinv_splitline and zinv_splitline[0] == "rate":
            zinv_rate = float(zinv_splitline[chn])
            if zinv_rate ==0: zinv_rate = 0.0001
      outfile_perchn.write("%.4f  " % (zinv_rate))
      pred_tot_rate += zinv_rate
 
      h1_zinv.SetBinContent(chn, zinv_rate)
      h1_zinv_syst.SetBinContent(chn, zinv_rate)
 
      qcd_file.seek(0, 0)
      for qcd_line in qcd_file:
         qcd_splitline = procline(qcd_line)
         if qcd_splitline and qcd_splitline[0] == "QCD_Data_CS":
            qcd_cs = float(qcd_splitline[chn])
         if qcd_splitline and qcd_splitline[0] == "QCD_TFactor":
            qcd_tfactor = float(qcd_splitline[chn])
            qcd_tfactor_str = qcd_splitline[chn]
         if qcd_splitline and qcd_splitline[0] == "QCD_otherBG_CS":
            tmp_qcd_contam_pred = float(qcd_splitline[chn])
      if qcd_cs != 0:
         outfile_perchn.write("%.4f  " % (qcd_cs*qcd_tfactor))
      else:
         outfile_perchn.write("%.4f  " % qcd_tfactor)
      if (qcd_cs-tmp_qcd_contam_pred)*qcd_tfactor >= 0: qcd_corr_rate = (qcd_cs-tmp_qcd_contam_pred)*qcd_tfactor
      else: qcd_corr_rate = 0
      pred_tot_rate += qcd_corr_rate

      h1_qcd.SetBinContent(chn, qcd_corr_rate)
      h1_qcd_syst.SetBinContent(chn, qcd_corr_rate)
 
      qcd_stat_unc_avg = math.sqrt(qcd_cs) * qcd_tfactor 
      prt_qcd_stat_unc_up = math.sqrt(qcd_cs) * qcd_tfactor
      prt_qcd_stat_unc_dn = math.sqrt(qcd_cs) * qcd_tfactor
      if qcd_cs ==0 : 
         qcd_stat_unc_avg = 1.84 * qcd_tfactor
         prt_qcd_stat_unc_up = 1.84 * qcd_tfactor
         prt_qcd_stat_unc_dn = 0.0
      pred_tot_stat += qcd_stat_unc_avg*qcd_stat_unc_avg 
      h1_qcd.SetBinError(chn, qcd_stat_unc_avg)
 
      ttz_file.seek(0, 0)
      for ttz_line in ttz_file:
         ttz_splitline = procline(ttz_line)
         if ttz_splitline and ttz_splitline[0] == "rate":
            ttz_rate = float(ttz_splitline[chn])
      outfile_perchn.write("%.4f  " % (ttz_rate))
#      outfile_perchn.write("\n")
      pred_tot_rate += ttz_rate

      h1_ttz.SetBinContent(chn, ttz_rate)
      h1_ttz_syst.SetBinContent(chn, ttz_rate)
 
      rare_file.seek(0, 0)
      for rare_line in rare_file:
         rare_splitline = procline(rare_line)
         if rare_splitline and rare_splitline[0] == "rate":
            rare_rate = float(rare_splitline[chn])
      outfile_perchn.write("%.5f  " % (rare_rate))
      outfile_perchn.write("\n")
      pred_tot_rate += rare_rate

      h1_rare.SetBinContent(chn, rare_rate)
      h1_rare_syst.SetBinContent(chn, rare_rate)
 
      outfile_perchn.write("---------------------------------------\n")
# signal stat. unc.
      outfile_perchn.write("signal_stat_unc_chn%d   gmN %0.0f %0.8f - - - - - - - - \n" % (chn, (signal_cs_event-signal_contam), signal_avg_weight))
      if signal_cs_event-signal_contam <=0 : signal_stat_unc_approx_rel = 0
      else:
         signal_stat_unc_approx_rel = math.sqrt(signal_cs_event-signal_contam)/(signal_cs_event-signal_contam)
      glb_signal_stat_up.append(signal_corr_rate * signal_stat_unc_approx_rel)
      glb_signal_stat_dn.append(signal_corr_rate * signal_stat_unc_approx_rel)

# muonCS stat. unc.
      lostle_file.seek(0, 0)
      for lostle_line in lostle_file:
         lostle_splitline = procline(lostle_line)
         if lostle_splitline and lostle_splitline[0] == "stat_unc_up":
            lostle_stat_unc_up_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "stat_unc_dn":
            lostle_stat_unc_dn_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "cs_event":
            lostle_cs_event = float(lostle_splitline[chn])

      lostle_stat_unc_avg = (lostle_stat_unc_up_rel + lostle_stat_unc_dn_rel)*0.5 * lostle_rate
      prt_lostle_stat_unc_up = lostle_stat_unc_up_rel * lostle_rate
      prt_lostle_stat_unc_dn = lostle_stat_unc_dn_rel * lostle_rate

      prt_lostle_stat_unc_up = math.sqrt(prt_lostle_stat_unc_up*prt_lostle_stat_unc_up + 1.84*lostle_avg_weight*1.84*lostle_avg_weight)
#      prt_lostle_stat_unc_dn = math.sqrt(prt_lostle_stat_unc_dn*prt_lostle_stat_unc_dn + 1.84*lostle_avg_weight*1.84*lostle_avg_weight)

      if prt_lostle_stat_unc_dn > lostle_rate: prt_lostle_stat_unc_dn = lostle_rate
      if lostle_cs_event == 0: 
         lostle_stat_unc_avg = 1.84 * lostle_avg_weight
         prt_lostle_stat_unc_up = 1.84 * lostle_avg_weight
         prt_lostle_stat_unc_dn = 0.0
      pred_tot_stat += lostle_stat_unc_avg*lostle_stat_unc_avg
 
      h1_lostle.SetBinError(chn, lostle_stat_unc_avg)
 
      hadtau_file.seek(0, 0)
      for hadtau_line in hadtau_file:
         hadtau_splitline = procline(hadtau_line)
         if hadtau_splitline and hadtau_splitline[0] == "stat_unc_up":
            hadtau_stat_unc_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "stat_unc_dn":
            hadtau_stat_unc_dn_rel = float(hadtau_splitline[chn])

      hadtau_stat_unc_avg = (hadtau_stat_unc_up_rel + hadtau_stat_unc_dn_rel)*0.5 * hadtau_rate
      pred_tot_stat += (hadtau_stat_unc_avg*hadtau_stat_unc_avg + (1.84*0.30)*(1.84*0.30))

      prt_hadtau_stat_unc_up = math.sqrt(hadtau_stat_unc_up_rel*hadtau_stat_unc_up_rel*hadtau_rate*hadtau_rate + (1.84*0.3)*(1.84*0.3))
      prt_hadtau_stat_unc_dn = hadtau_stat_unc_dn_rel*hadtau_rate
 
      h1_hadtau.SetBinError(chn, math.sqrt(hadtau_stat_unc_avg*hadtau_stat_unc_avg + (1.84*0.30)*(1.84*0.30)))
 
      if lostle_stat_unc_dn_rel ==1: lostle_stat_unc_dn_rel -= 0.001
      if hadtau_stat_unc_dn_rel ==1: hadtau_stat_unc_dn_rel -= 0.001
 
      outfile_perchn.write("muonCS_stat_unc_chn%d    lnN    - %.4f/%.4f - %.4f/%.4f - - - - -\n" % (chn, 1-lostle_stat_unc_dn_rel, 1+lostle_stat_unc_up_rel, 1-hadtau_stat_unc_dn_rel, 1+hadtau_stat_unc_up_rel))
#      if lostle_cs_event !=0: 
#         outfile_perchn.write("stat_unc_HighW_chn%d     gmN  0 - - %.4f - %.4f - - - -\n" % (chn, lostle_avg_weight, 0.30))
#      else:
#         outfile_perchn.write("stat_unc_HighW_chn%d     gmN  1 - - %.4f - %.4f - - - -\n" % (chn, lostle_avg_weight*0.55, 0.30*0.55))
      outfile_perchn.write("stat_unc_HighW_chn%d    lnU    - - 10 - 10 - - - -\n" % (chn))

# zinv stat. unc.
      zinv_file.seek(0, 0)
      for zinv_line in zinv_file:
         zinv_splitline = procline(zinv_line)
         if zinv_splitline and zinv_splitline[0] == "cs_event":
            zinv_cs_event = float(zinv_splitline[chn])
         if zinv_splitline and zinv_splitline[0] == "avg_weight":
            zinv_avg_weight = float(zinv_splitline[chn])

      zinv_stat_unc_avg = math.sqrt(zinv_cs_event) * zinv_avg_weight
      prt_zinv_stat_unc_up = math.sqrt(zinv_cs_event) * zinv_avg_weight
      prt_zinv_stat_unc_dn = math.sqrt(zinv_cs_event) * zinv_avg_weight
      if zinv_cs_event == 0: 
         zinv_stat_unc_avg = 1.84 * zinv_avg_weight
         prt_zinv_stat_unc_up = 1.84 * zinv_avg_weight
         prt_zinv_stat_unc_dn = 0.0
      pred_tot_stat += zinv_stat_unc_avg*zinv_stat_unc_avg
 
      h1_zinv.SetBinError(chn, zinv_stat_unc_avg)

#      if int(channels) != 45: 
#         outfile_perchn.write("zinv_stat_unc_chn%d  gmN %0.0f - - - - - %0.4f - - -\n" % (chn, zinv_cs_event, zinv_avg_weight))

# ttz stat. unc.
      ttz_file.seek(0, 0)
      for ttz_line in ttz_file:
         ttz_splitline = procline(ttz_line)
         if ttz_splitline and ttz_splitline[0] == "cs_event":
            ttz_cs_event = float(ttz_splitline[chn])
         if ttz_splitline and ttz_splitline[0] == "avg_weight":
            ttz_avg_weight = float(ttz_splitline[chn])

      ttz_stat_unc_avg = math.sqrt(ttz_cs_event) * ttz_avg_weight
      prt_ttz_stat_unc_up = math.sqrt(ttz_cs_event) * ttz_avg_weight
      prt_ttz_stat_unc_dn = math.sqrt(ttz_cs_event) * ttz_avg_weight
      if ttz_cs_event == 0: 
         ttz_stat_unc_avg = 1.84 * ttz_avg_weight
         prt_ttz_stat_unc_up = 1.84 * ttz_avg_weight
         prt_ttz_stat_unc_up = 0.0
      pred_tot_stat += ttz_stat_unc_avg*ttz_stat_unc_avg
 
      h1_ttz.SetBinError(chn, ttz_stat_unc_avg)
 
      outfile_perchn.write("ttz_stat_unc_chn%d  gmN %0.0f - - - - - - - %0.4f -\n" % (chn, ttz_cs_event, ttz_avg_weight))
#      ttz_stat_approx_square = ttz_cs_event * 2 # 291495 - 106505 gets roughly 3:1 for pos:neg
#      ttz_stat_approx_rel = math.sqrt(ttz_stat_approx_square)/ttz_cs_event
#      outfile_perchn.write("ttz_stat_unc_chn%d  lnN - - - - - - - %0.4f -\n" % (chn, 1+ttz_stat_approx_rel))

# rare stat. unc.
      rare_file.seek(0, 0)
      for rare_line in rare_file:
         rare_splitline = procline(rare_line)
         if rare_splitline and rare_splitline[0] == "cs_event":
            rare_cs_event = float(rare_splitline[chn])
         if rare_splitline and rare_splitline[0] == "avg_weight":
            rare_avg_weight = float(rare_splitline[chn])

      rare_stat_unc_avg = math.sqrt(rare_cs_event) * rare_avg_weight
      prt_rare_stat_unc_up = math.sqrt(rare_cs_event) * rare_avg_weight
      prt_rare_stat_unc_dn = math.sqrt(rare_cs_event) * rare_avg_weight
      if rare_cs_event == 0: 
         rare_stat_unc_avg = 1.84 * rare_avg_weight
         prt_rare_stat_unc_up = 1.84 * rare_avg_weight
         prt_rare_stat_unc_up = 0.0
      pred_tot_stat += rare_stat_unc_avg*rare_stat_unc_avg
 
      h1_rare.SetBinError(chn, rare_stat_unc_avg)

#      if int(channels) != 45: 
#         outfile_perchn.write("rare_stat_unc_chn%d  gmN %0.0f - - - - - - - - %0.5f\n" % (chn, rare_cs_event, rare_avg_weight))

      pred_tot_stat = math.sqrt(pred_tot_stat)

# Systematic uncertainties
      outfile_perchn.write("-------------------------------------------\n")
# signal
      signal_file.seek(0, 0)
      for signal_line in signal_file:
         signal_splitline = procline(signal_line)
         if signal_splitline and signal_splitline[0] == "lumi_unc_up":
            signal_syst_lumi_up_rel = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "lumi_unc_dn":
            signal_syst_lumi_dn_rel = float(signal_splitline[chn])

         if signal_splitline and signal_splitline[0] == "bTagSF_up":
            signal_syst_bTagSF_up_rel = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "bTagSF_dn":
            signal_syst_bTagSF_dn_rel = float(signal_splitline[chn])
         
         if signal_splitline and signal_splitline[0] == "mistagSF_up":
            signal_syst_mistagSF_up_rel = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "mistagSF_dn":
            signal_syst_mistagSF_dn_rel = float(signal_splitline[chn])
        
         if signal_splitline and signal_splitline[0] == "scaleUnc_up":
            signal_syst_scaleUnc_up_rel = float(signal_splitline[chn]) 
         if signal_splitline and signal_splitline[0] == "scaleUnc_dn":
            signal_syst_scaleUnc_dn_rel = float(signal_splitline[chn]) 

         if signal_splitline and signal_splitline[0] == "isrUnc_up":
            signal_syst_isrUnc_up_rel = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "isrUnc_dn":
            signal_syst_isrUnc_dn_rel = float(signal_splitline[chn])

         if signal_splitline and signal_splitline[0] == "metMag_up":
            signal_syst_metMag_up_rel = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "metMag_dn":
            signal_syst_metMag_dn_rel = float(signal_splitline[chn])

         if signal_splitline and signal_splitline[0] == "jetJEC_up":
            signal_syst_jetJEC_up_rel = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "jetJEC_dn":
            signal_syst_jetJEC_dn_rel = float(signal_splitline[chn])

         if signal_splitline and signal_splitline[0] == "trigUnc_up":
            signal_syst_trigUnc_up_rel = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "trigUnc_dn":
            signal_syst_trigUnc_dn_rel = float(signal_splitline[chn])

         if signal_splitline and signal_splitline[0] == "lepVetoUnc_up":
            signal_syst_lepVetoUnc_up_rel = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "lepVetoUnc_dn":
            signal_syst_lepVetoUnc_dn_rel = float(signal_splitline[chn])

         if signal_splitline and signal_splitline[0] == "genTopSF_up":
            signal_syst_genTopSFUnc_up_rel = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "genTopSF_dn":
            signal_syst_genTopSFUnc_dn_rel = float(signal_splitline[chn])

         if signal_splitline and signal_splitline[0] == "mistaggenTopSF_up":
            signal_syst_mistaggenTopSFUnc_up_rel = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "genTopSF_dn":
            signal_syst_mistaggenTopSFUnc_dn_rel = float(signal_splitline[chn])

         if signal_splitline and signal_splitline[0] == "data_vs_MC_recoTop_unc_up":
            signal_syst_data_vs_MC_recoTop_unc_up_rel = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "data_vs_MC_recoTop_unc_dn":
            signal_syst_data_vs_MC_recoTop_unc_dn_rel = float(signal_splitline[chn])

# Ignore metPhi since it's so small
#         if signal_splitline and signal_splitline[0] == "metPhi_up":
#            signal_syst_metPhi_up_rel = float(signal_splitline[chn])
#         if signal_splitline and signal_splitline[0] == "metPhi_dn":
#            signal_syst_metPhi_dn_rel = float(signal_splitline[chn])

#         if signal_splitline and signal_splitline[0] == "syst_unc_up":
#            signal_syst_up_rel = float(signal_splitline[chn])
#         if signal_splitline and signal_splitline[0] == "syst_unc_dn":
#            signal_syst_dn_rel = float(signal_splitline[chn])

# Symmetrize the scale unc using the envelope the absolute values of both up and dn 
      if abs(signal_syst_scaleUnc_up_rel) > abs(signal_syst_scaleUnc_dn_rel):
         signal_syst_scaleUnc_up_rel = abs(signal_syst_scaleUnc_up_rel)
         signal_syst_scaleUnc_dn_rel = abs(signal_syst_scaleUnc_up_rel) 
      else:
         signal_syst_scaleUnc_up_rel = abs(signal_syst_scaleUnc_dn_rel)
         signal_syst_scaleUnc_dn_rel = abs(signal_syst_scaleUnc_dn_rel) 
         
# Symmetrize the metMag unc using the envelope the absolute values of both up and dn 
      if abs(signal_syst_metMag_up_rel) > abs(signal_syst_metMag_dn_rel):
         signal_syst_metMag_up_rel = abs(signal_syst_metMag_up_rel)
         signal_syst_metMag_dn_rel = abs(signal_syst_metMag_up_rel) 
      else:
         signal_syst_metMag_up_rel = abs(signal_syst_metMag_dn_rel)
         signal_syst_metMag_dn_rel = abs(signal_syst_metMag_dn_rel) 
         
# Symmetrize the jetJEC unc using the envelope the absolute values of both up and dn 
      if abs(signal_syst_jetJEC_up_rel) > abs(signal_syst_jetJEC_dn_rel):
         signal_syst_jetJEC_up_rel = abs(signal_syst_jetJEC_up_rel)
         signal_syst_jetJEC_dn_rel = abs(signal_syst_jetJEC_up_rel) 
      else:
         signal_syst_jetJEC_up_rel = abs(signal_syst_jetJEC_dn_rel)
         signal_syst_jetJEC_dn_rel = abs(signal_syst_jetJEC_dn_rel) 

# If number of signal events less than 5, the uncertainties are subject to stat. flucturation.
# Therefore, for signal events less than 5, their uncertainties are capped at the maximum uncertainties
# measured for number of signal events >= 100
#      if signal_cs_event < 5:
      if True:
#         signal_syst_scaleUnc_up_rel = 0.05
#         signal_syst_scaleUnc_dn_rel = 0.05
#         signal_syst_jetJEC_up_rel = 0.30
#         signal_syst_jetJEC_dn_rel = 0.30
#         signal_syst_metMag_up_rel = 0.20
#         signal_syst_metMag_dn_rel = 0.20
         if signal_syst_scaleUnc_up_rel > 0.035 : signal_syst_scaleUnc_up_rel = 0.035
         if signal_syst_scaleUnc_dn_rel > 0.035 : signal_syst_scaleUnc_dn_rel = 0.035
         if signal_syst_jetJEC_up_rel > 0.31 : signal_syst_jetJEC_up_rel = 0.31
         if signal_syst_jetJEC_dn_rel > 0.31 : signal_syst_jetJEC_dn_rel = 0.31
         if signal_syst_metMag_up_rel > 0.21 : signal_syst_metMag_up_rel = 0.21
         if signal_syst_metMag_dn_rel > 0.21 : signal_syst_metMag_dn_rel = 0.21

      signal_syst_genTopSFUnc_mistaggenTopSFUnc_dn_rel = math.sqrt( signal_syst_genTopSFUnc_dn_rel*signal_syst_genTopSFUnc_dn_rel + signal_syst_mistaggenTopSFUnc_dn_rel*signal_syst_mistaggenTopSFUnc_dn_rel)
      signal_syst_genTopSFUnc_mistaggenTopSFUnc_up_rel = math.sqrt( signal_syst_genTopSFUnc_up_rel*signal_syst_genTopSFUnc_up_rel + signal_syst_mistaggenTopSFUnc_up_rel*signal_syst_mistaggenTopSFUnc_up_rel)

      signal_syst_lumi_trigUnc_dn_rel = math.sqrt( signal_syst_lumi_dn_rel*signal_syst_lumi_dn_rel + signal_syst_trigUnc_dn_rel*signal_syst_trigUnc_dn_rel )
      signal_syst_lumi_trigUnc_up_rel = math.sqrt( signal_syst_lumi_up_rel*signal_syst_lumi_up_rel + signal_syst_trigUnc_up_rel*signal_syst_trigUnc_up_rel )

      signal_syst_lumi_trigUnc_dataVSMCrecoTop_dn_rel = math.sqrt( signal_syst_lumi_trigUnc_dn_rel*signal_syst_lumi_trigUnc_dn_rel + signal_syst_data_vs_MC_recoTop_unc_dn_rel*signal_syst_data_vs_MC_recoTop_unc_dn_rel )
      signal_syst_lumi_trigUnc_dataVSMCrecoTop_up_rel = math.sqrt( signal_syst_lumi_trigUnc_up_rel*signal_syst_lumi_trigUnc_up_rel + signal_syst_data_vs_MC_recoTop_unc_up_rel*signal_syst_data_vs_MC_recoTop_unc_up_rel )

      signal_syst_metMag_jetJEC_dn_rel = math.sqrt( signal_syst_metMag_dn_rel*signal_syst_metMag_dn_rel + signal_syst_jetJEC_dn_rel*signal_syst_jetJEC_dn_rel )
      signal_syst_metMag_jetJEC_up_rel = math.sqrt( signal_syst_metMag_up_rel*signal_syst_metMag_up_rel + signal_syst_jetJEC_up_rel*signal_syst_jetJEC_up_rel )

      signal_tot_syst_dn_rel = math.sqrt( signal_syst_lumi_dn_rel*signal_syst_lumi_dn_rel + signal_syst_scaleUnc_dn_rel*signal_syst_scaleUnc_dn_rel + signal_syst_isrUnc_dn_rel*signal_syst_isrUnc_dn_rel + signal_syst_bTagSF_dn_rel*signal_syst_bTagSF_dn_rel + signal_syst_mistagSF_dn_rel*signal_syst_mistagSF_dn_rel + signal_syst_metMag_dn_rel*signal_syst_metMag_dn_rel + signal_syst_jetJEC_dn_rel*signal_syst_jetJEC_dn_rel + signal_syst_trigUnc_dn_rel*signal_syst_trigUnc_dn_rel + signal_syst_lepVetoUnc_dn_rel*signal_syst_lepVetoUnc_dn_rel + signal_syst_genTopSFUnc_dn_rel*signal_syst_genTopSFUnc_dn_rel + signal_syst_mistaggenTopSFUnc_dn_rel*signal_syst_mistaggenTopSFUnc_dn_rel + signal_syst_data_vs_MC_recoTop_unc_dn_rel*signal_syst_data_vs_MC_recoTop_unc_dn_rel )
      signal_tot_syst_up_rel = math.sqrt( signal_syst_lumi_up_rel*signal_syst_lumi_up_rel + signal_syst_scaleUnc_up_rel*signal_syst_scaleUnc_up_rel + signal_syst_isrUnc_up_rel*signal_syst_isrUnc_up_rel + signal_syst_bTagSF_up_rel*signal_syst_bTagSF_up_rel + signal_syst_mistagSF_up_rel*signal_syst_mistagSF_up_rel + signal_syst_metMag_up_rel*signal_syst_metMag_up_rel + signal_syst_jetJEC_up_rel*signal_syst_jetJEC_up_rel + signal_syst_trigUnc_up_rel*signal_syst_trigUnc_up_rel + signal_syst_lepVetoUnc_up_rel*signal_syst_lepVetoUnc_up_rel + signal_syst_genTopSFUnc_up_rel*signal_syst_genTopSFUnc_up_rel + signal_syst_mistaggenTopSFUnc_up_rel*signal_syst_mistaggenTopSFUnc_up_rel + signal_syst_data_vs_MC_recoTop_unc_up_rel*signal_syst_data_vs_MC_recoTop_unc_up_rel )

      glb_signal_syst_dn.append(signal_corr_rate * signal_tot_syst_dn_rel)
      glb_signal_syst_up.append(signal_corr_rate * signal_tot_syst_up_rel)

      if signal_syst_lumi_dn_rel ==1: signal_syst_lumi_dn_rel -= 0.001
      if signal_syst_scaleUnc_dn_rel ==1: signal_syst_scaleUnc_dn_rel -= 0.001
      if signal_syst_isrUnc_dn_rel ==1: signal_syst_isrUnc_dn_rel -= 0.001
      if signal_syst_bTagSF_dn_rel ==1: signal_syst_bTagSF_dn_rel -= 0.001
      if signal_syst_mistagSF_dn_rel ==1: signal_syst_mistagSF_dn_rel -= 0.001
      if signal_syst_metMag_dn_rel ==1: signal_syst_metMag_dn_rel -= 0.001
      if signal_syst_jetJEC_dn_rel ==1: signal_syst_jetJEC_dn_rel -= 0.001
      if signal_syst_trigUnc_dn_rel ==1: signal_syst_trigUnc_dn_rel -= 0.001
      if signal_syst_lepVetoUnc_dn_rel ==1: signal_syst_lepVetoUnc_dn_rel -= 0.001
      if signal_syst_genTopSFUnc_dn_rel ==1: signal_syst_genTopSFUnc_dn_rel -= 0.001
      if signal_syst_data_vs_MC_recoTop_unc_dn_rel ==1: signal_syst_data_vs_MC_recoTop_unc_dn_rel -= 0.001

      if signal_syst_genTopSFUnc_mistaggenTopSFUnc_dn_rel ==1: signal_syst_genTopSFUnc_mistaggenTopSFUnc_dn_rel -= 0.001
      if signal_syst_lumi_trigUnc_dn_rel ==1: signal_syst_lumi_trigUnc_dn_rel -= 0.001
      if signal_syst_metMag_jetJEC_dn_rel ==1: signal_syst_metMag_jetJEC_dn_rel -= 0.001
      if signal_syst_lumi_trigUnc_dataVSMCrecoTop_dn_rel ==1: signal_syst_lumi_trigUnc_dataVSMCrecoTop_dn_rel -= 0.001

#      outfile_perchn.write("signal_lumi  lnN %.4f/%.4f - - - - - - - -\n" % (1-signal_syst_lumi_dn_rel, 1+signal_syst_lumi_up_rel) )
      outfile_perchn.write("signal_scaleUnc  lnN %.4f/%.4f - - - - - - - -\n" % (1-signal_syst_scaleUnc_dn_rel, 1+signal_syst_scaleUnc_up_rel) )
      outfile_perchn.write("signal_isrUnc  lnN %.4f/%.4f - - - - - - - -\n" % (1-signal_syst_isrUnc_dn_rel, 1+signal_syst_isrUnc_up_rel) )
      outfile_perchn.write("signal_bTagSF  lnN %.4f/%.4f - - - - - - - -\n" % (1-signal_syst_bTagSF_dn_rel, 1+signal_syst_bTagSF_up_rel) )
      outfile_perchn.write("signal_mistagSF  lnN %.4f/%.4f - - - - - - - -\n" % (1-signal_syst_mistagSF_dn_rel, 1+signal_syst_mistagSF_up_rel) )
#      outfile_perchn.write("signal_metMag  lnN %.4f/%.4f - - - - - - - -\n" % (1-signal_syst_metMag_dn_rel, 1+signal_syst_metMag_up_rel) )
#      outfile_perchn.write("signal_jetJEC  lnN %.4f/%.4f - - - - - - - -\n" % (1-signal_syst_jetJEC_dn_rel, 1+signal_syst_jetJEC_up_rel) )
      outfile_perchn.write("signal_metMag_jetJEC  lnN %.4f/%.4f - - - - - - - -\n" % (1-signal_syst_metMag_jetJEC_dn_rel, 1+signal_syst_metMag_jetJEC_up_rel) )
#      outfile_perchn.write("signal_trigUnc  lnN %.4f/%.4f - - - - - - - -\n" % (1-signal_syst_trigUnc_dn_rel, 1+signal_syst_trigUnc_up_rel) )
      outfile_perchn.write("signal_lepVetoUnc  lnN %.4f/%.4f - - - - - - - -\n" % (1-signal_syst_lepVetoUnc_dn_rel, 1+signal_syst_lepVetoUnc_up_rel) )
#      outfile_perchn.write("signal_genTopSFUnc  lnN %.4f/%.4f - - - - - - - -\n" %(1-signal_syst_genTopSFUnc_dn_rel, 1+signal_syst_genTopSFUnc_up_rel) )
#      outfile_perchn.write("signal_mistaggenTopSFUnc  lnN %.4f/%.4f - - - - - - - -\n" %(1-signal_syst_mistaggenTopSFUnc_dn_rel, 1+signal_syst_mistaggenTopSFUnc_up_rel) )
      outfile_perchn.write("signal_lumi_trigUnc_dataVSMCrecoTop  lnN %.4f/%.4f - - - - - - - -\n" % (1-signal_syst_lumi_trigUnc_dataVSMCrecoTop_dn_rel, 1+signal_syst_lumi_trigUnc_dataVSMCrecoTop_up_rel) )
      outfile_perchn.write("signal_genTopSFUnc_mistaggenTopSFUnc  lnN %.4f/%.4f - - - - - - - -\n" %(1-signal_syst_genTopSFUnc_mistaggenTopSFUnc_dn_rel, 1+signal_syst_genTopSFUnc_mistaggenTopSFUnc_up_rel) )
         
# lostle 
      lostle_file.seek(0, 0)
      for lostle_line in lostle_file:
         lostle_splitline = procline(lostle_line)
         if lostle_splitline and lostle_splitline[0] == "syst_unc_closure_up":
            lostle_syst_closure_up_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "syst_unc_closure_dn":
            lostle_syst_closure_dn_rel = float(lostle_splitline[chn])

         if lostle_splitline and lostle_splitline[0] == "syst_unc_purity_up":
            lostle_syst_purity_up_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "syst_unc_purity_dn":
            lostle_syst_purity_dn_rel = float(lostle_splitline[chn])

         if lostle_splitline and lostle_splitline[0] == "syst_unc_dimuon_up":
            lostle_syst_dimuon_up_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "syst_unc_dimuon_dn":
            lostle_syst_dimuon_dn_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "syst_unc_diele_up":
            lostle_syst_diele_up_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "syst_unc_diele_dn":
            lostle_syst_diele_dn_rel = float(lostle_splitline[chn])

         if lostle_splitline and lostle_splitline[0] == "syst_unc_mt_up":
            lostle_syst_mt_up_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "syst_unc_mt_dn":
            lostle_syst_mt_dn_rel = float(lostle_splitline[chn])

         if lostle_splitline and lostle_splitline[0] == "syst_unc_acc_up":
            lostle_syst_acc_up_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "syst_unc_acc_dn":
            lostle_syst_acc_dn_rel = float(lostle_splitline[chn])

         if lostle_splitline and lostle_splitline[0] == "syst_unc_muiso_up":
            lostle_syst_muiso_up_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "syst_unc_muiso_dn":
            lostle_syst_muiso_dn_rel = float(lostle_splitline[chn])

         if lostle_splitline and lostle_splitline[0] == "syst_unc_eiso_up":
            lostle_syst_eiso_up_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "syst_unc_eiso_dn":
            lostle_syst_eiso_dn_rel = float(lostle_splitline[chn])

         if lostle_splitline and lostle_splitline[0] == "syst_unc_mureco_up":
            lostle_syst_mureco_up_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "syst_unc_mureco_dn":
            lostle_syst_mureco_dn_rel = float(lostle_splitline[chn])

         if lostle_splitline and lostle_splitline[0] == "syst_unc_ereco_up":
            lostle_syst_ereco_up_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "syst_unc_ereco_dn":
            lostle_syst_ereco_dn_rel = float(lostle_splitline[chn])

         if lostle_splitline and lostle_splitline[0] == "syst_unc_isotrk_up":
            lostle_syst_isotrk_up_rel = float(lostle_splitline[chn])
         if lostle_splitline and lostle_splitline[0] == "syst_unc_isotrk_dn":
            lostle_syst_isotrk_dn_rel = float(lostle_splitline[chn])

      lostle_syst_dilep_up_rel = math.sqrt(lostle_syst_dimuon_up_rel*lostle_syst_dimuon_up_rel + lostle_syst_diele_up_rel*lostle_syst_diele_up_rel) 
      lostle_syst_dilep_dn_rel = math.sqrt(lostle_syst_dimuon_dn_rel*lostle_syst_dimuon_dn_rel + lostle_syst_diele_dn_rel*lostle_syst_diele_dn_rel) 

      lostle_syst_closure_avg = (lostle_syst_closure_up_rel + lostle_syst_closure_dn_rel)*0.5
      lostle_syst_purity_avg = (lostle_syst_purity_up_rel + lostle_syst_purity_dn_rel)*0.5
      lostle_syst_dilep_avg = (lostle_syst_dilep_up_rel + lostle_syst_dilep_dn_rel)*0.5
      lostle_syst_mt_avg = (lostle_syst_mt_up_rel + lostle_syst_mt_dn_rel)*0.5
      lostle_syst_acc_avg = (lostle_syst_acc_up_rel + lostle_syst_acc_dn_rel)*0.5
      lostle_syst_muiso_avg = (lostle_syst_muiso_up_rel + lostle_syst_muiso_dn_rel)*0.5
      lostle_syst_mureco_avg = (lostle_syst_mureco_up_rel + lostle_syst_mureco_dn_rel)*0.5
      lostle_syst_eiso_avg = (lostle_syst_eiso_up_rel + lostle_syst_eiso_dn_rel)*0.5
      lostle_syst_ereco_avg = (lostle_syst_ereco_up_rel + lostle_syst_ereco_dn_rel)*0.5
      lostle_syst_isotrk_avg = (lostle_syst_isotrk_up_rel + lostle_syst_isotrk_dn_rel)*0.5

      lostle_syst_unc_avg = math.sqrt( lostle_syst_closure_avg * lostle_syst_closure_avg + lostle_syst_purity_avg * lostle_syst_purity_avg + lostle_syst_dilep_avg * lostle_syst_dilep_avg + lostle_syst_mt_avg * lostle_syst_mt_avg + lostle_syst_acc_avg * lostle_syst_acc_avg + lostle_syst_muiso_avg * lostle_syst_muiso_avg + lostle_syst_mureco_avg * lostle_syst_mureco_avg + lostle_syst_eiso_avg * lostle_syst_eiso_avg + lostle_syst_ereco_avg * lostle_syst_ereco_avg + lostle_syst_isotrk_avg * lostle_syst_isotrk_avg)

      prt_lostle_syst_unc_up = lostle_rate * math.sqrt( lostle_syst_closure_up_rel * lostle_syst_closure_up_rel + lostle_syst_purity_up_rel * lostle_syst_purity_up_rel + lostle_syst_dilep_up_rel * lostle_syst_dilep_up_rel + lostle_syst_mt_up_rel * lostle_syst_mt_up_rel + lostle_syst_acc_up_rel * lostle_syst_acc_up_rel + lostle_syst_muiso_up_rel * lostle_syst_muiso_up_rel + lostle_syst_mureco_up_rel * lostle_syst_mureco_up_rel + lostle_syst_eiso_up_rel * lostle_syst_eiso_up_rel + lostle_syst_ereco_up_rel * lostle_syst_ereco_up_rel + lostle_syst_isotrk_up_rel * lostle_syst_isotrk_up_rel )
      prt_lostle_syst_unc_dn = lostle_rate * math.sqrt( lostle_syst_closure_dn_rel * lostle_syst_closure_dn_rel + lostle_syst_purity_dn_rel * lostle_syst_purity_dn_rel + lostle_syst_dilep_dn_rel * lostle_syst_dilep_dn_rel + lostle_syst_mt_dn_rel * lostle_syst_mt_dn_rel + lostle_syst_acc_dn_rel * lostle_syst_acc_dn_rel + lostle_syst_muiso_dn_rel * lostle_syst_muiso_dn_rel + lostle_syst_mureco_dn_rel * lostle_syst_mureco_dn_rel + lostle_syst_eiso_dn_rel * lostle_syst_eiso_dn_rel + lostle_syst_ereco_dn_rel * lostle_syst_ereco_dn_rel + lostle_syst_isotrk_dn_rel * lostle_syst_isotrk_dn_rel )

      if prt_lostle_syst_unc_dn > lostle_rate : prt_lostle_syst_unc_dn = lostle_rate

      h1_lostle_syst.SetBinError(chn, lostle_rate*lostle_syst_unc_avg)

      if lostle_syst_closure_dn_rel >=1: lostle_syst_closure_dn_rel = 1-0.001
      if lostle_syst_purity_dn_rel >=1: lostle_syst_purity_dn_rel = 1-0.001
      if lostle_syst_eiso_dn_rel >=1: lostle_syst_eiso_dn_rel = 1-0.001
      if lostle_syst_ereco_dn_rel >=1: lostle_syst_ereco_dn_rel = 1-0.001
      if lostle_syst_isotrk_dn_rel >=1: lostle_syst_isotrk_dn_rel = 1-0.001
      if lostle_syst_acc_dn_rel >=1: lostle_syst_acc_dn_rel = 1-0.001

      outfile_perchn.write("lostle_closure_chn%d   lnN - %.4f/%.4f - - - - - - -\n" % (chn, 1-lostle_syst_closure_dn_rel, 1+lostle_syst_closure_up_rel))
#      outfile_perchn.write("lostle_purity_chn%d   lnN - %.4f/%.4f - - - - - - -\n" % (chn, 1-lostle_syst_purity_dn_rel, 1+lostle_syst_purity_up_rel))
      outfile_perchn.write("lostle_purity   lnN - %.4f/%.4f - - - - - - -\n" % (1-lostle_syst_purity_dn_rel, 1+lostle_syst_purity_up_rel))
      outfile_perchn.write("lostle_eiso  lnN - %.4f/%.4f - - - - - - -\n" % (1-lostle_syst_eiso_dn_rel, 1+lostle_syst_eiso_up_rel))
      outfile_perchn.write("lostle_ereco  lnN - %.4f/%.4f - - - - - - -\n" % (1-lostle_syst_ereco_dn_rel, 1+lostle_syst_ereco_up_rel))
      outfile_perchn.write("lostle_isotrk  lnN - %.4f/%.4f - - - - - - -\n" % (1-lostle_syst_isotrk_dn_rel, 1+lostle_syst_isotrk_up_rel))
      outfile_perchn.write("lostle_acc  lnN - %.4f/%.4f - - - - - - -\n" % (1-lostle_syst_acc_dn_rel, 1+lostle_syst_acc_up_rel))

# hadtau
      hadtau_file.seek(0, 0)
      for hadtau_line in hadtau_file:
         hadtau_splitline = procline(hadtau_line)
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_pdf_up":
            hadtau_syst_pdf_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_pdf_dn":
            hadtau_syst_pdf_dn_rel = float(hadtau_splitline[chn])

         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_MCScale_up":
            hadtau_syst_MCScale_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_MCScale_dn":
            hadtau_syst_MCScale_dn_rel = float(hadtau_splitline[chn])

         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_Mt_up":
            hadtau_syst_Mt_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_Mt_dn":
            hadtau_syst_Mt_dn_rel = float(hadtau_splitline[chn])

         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_taumu_up":
            hadtau_syst_taumu_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_taumu_dn":
            hadtau_syst_taumu_dn_rel = float(hadtau_splitline[chn])

         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_isotrk_up":
            hadtau_syst_isotrk_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_isotrk_dn":
            hadtau_syst_isotrk_dn_rel = float(hadtau_splitline[chn])

         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_mureco_up":
            hadtau_syst_mureco_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_mureco_dn":
            hadtau_syst_mureco_dn_rel = float(hadtau_splitline[chn])

         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_muiso_up":
            hadtau_syst_muiso_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_muiso_dn":
            hadtau_syst_muiso_dn_rel = float(hadtau_splitline[chn])

         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_mistag_up":
            hadtau_syst_mistag_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_mistag_dn":
            hadtau_syst_mistag_dn_rel = float(hadtau_splitline[chn])

         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_temp_up":
            hadtau_syst_temp_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_temp_dn":
            hadtau_syst_temp_dn_rel = float(hadtau_splitline[chn])

         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_closure_up":
            hadtau_syst_closure_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_closure_dn":
            hadtau_syst_closure_dn_rel = float(hadtau_splitline[chn])

         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_llovr_up":
            hadtau_syst_llovr_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_llovr_dn":
            hadtau_syst_llovr_dn_rel = float(hadtau_splitline[chn])

         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_trg_up":
            hadtau_syst_trg_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_trg_dn":
            hadtau_syst_trg_dn_rel = float(hadtau_splitline[chn])

      hadtau_syst_closure_avg = (hadtau_syst_closure_up_rel + hadtau_syst_closure_dn_rel)*0.5
      hadtau_syst_llovr_avg = (hadtau_syst_llovr_up_rel + hadtau_syst_llovr_dn_rel)*0.5
      hadtau_syst_Mt_avg = (hadtau_syst_Mt_up_rel + hadtau_syst_Mt_dn_rel)*0.5
      hadtau_syst_pdf_avg = (hadtau_syst_pdf_up_rel + hadtau_syst_pdf_dn_rel)*0.5
      hadtau_syst_MCScale_avg = (hadtau_syst_MCScale_up_rel + hadtau_syst_MCScale_dn_rel)*0.5
      hadtau_syst_muiso_avg = (hadtau_syst_muiso_up_rel + hadtau_syst_muiso_dn_rel)*0.5
      hadtau_syst_mureco_avg = (hadtau_syst_mureco_up_rel + hadtau_syst_mureco_dn_rel)*0.5
      hadtau_syst_isotrk_avg = (hadtau_syst_isotrk_up_rel + hadtau_syst_isotrk_dn_rel)*0.5
      hadtau_syst_taumu_avg = (hadtau_syst_taumu_up_rel + hadtau_syst_taumu_dn_rel)*0.5
      hadtau_syst_mistag_avg = (hadtau_syst_mistag_up_rel + hadtau_syst_mistag_dn_rel)*0.5
      hadtau_syst_temp_avg = (hadtau_syst_temp_up_rel + hadtau_syst_temp_dn_rel)*0.5
      hadtau_syst_trg_avg = (hadtau_syst_trg_up_rel + hadtau_syst_trg_dn_rel)*0.5

      hadtau_syst_unc_avg = math.sqrt( hadtau_syst_closure_avg * hadtau_syst_closure_avg + hadtau_syst_llovr_avg * hadtau_syst_llovr_avg + hadtau_syst_Mt_avg * hadtau_syst_Mt_avg + hadtau_syst_pdf_avg * hadtau_syst_pdf_avg + hadtau_syst_MCScale_avg * hadtau_syst_MCScale_avg + hadtau_syst_muiso_avg * hadtau_syst_muiso_avg + hadtau_syst_mureco_avg * hadtau_syst_mureco_avg + hadtau_syst_isotrk_avg * hadtau_syst_isotrk_avg + hadtau_syst_taumu_avg * hadtau_syst_taumu_avg + hadtau_syst_mistag_avg * hadtau_syst_mistag_avg + hadtau_syst_temp_avg * hadtau_syst_temp_avg + hadtau_syst_trg_avg * hadtau_syst_trg_avg )

      prt_hadtau_syst_unc_up = hadtau_rate * math.sqrt( hadtau_syst_closure_up_rel * hadtau_syst_closure_up_rel + hadtau_syst_llovr_up_rel * hadtau_syst_llovr_up_rel + hadtau_syst_Mt_up_rel * hadtau_syst_Mt_up_rel + hadtau_syst_pdf_up_rel * hadtau_syst_pdf_up_rel + hadtau_syst_MCScale_up_rel * hadtau_syst_MCScale_up_rel + hadtau_syst_muiso_up_rel * hadtau_syst_muiso_up_rel + hadtau_syst_mureco_up_rel * hadtau_syst_mureco_up_rel + hadtau_syst_isotrk_up_rel * hadtau_syst_isotrk_up_rel + hadtau_syst_taumu_up_rel * hadtau_syst_taumu_up_rel + hadtau_syst_mistag_up_rel * hadtau_syst_mistag_up_rel + hadtau_syst_temp_up_rel * hadtau_syst_temp_up_rel + hadtau_syst_trg_up_rel * hadtau_syst_trg_up_rel )
      prt_hadtau_syst_unc_dn = hadtau_rate * math.sqrt( hadtau_syst_closure_dn_rel * hadtau_syst_closure_dn_rel + hadtau_syst_llovr_dn_rel * hadtau_syst_llovr_dn_rel + hadtau_syst_Mt_dn_rel * hadtau_syst_Mt_dn_rel + hadtau_syst_pdf_dn_rel * hadtau_syst_pdf_dn_rel + hadtau_syst_MCScale_dn_rel * hadtau_syst_MCScale_dn_rel + hadtau_syst_muiso_dn_rel * hadtau_syst_muiso_dn_rel + hadtau_syst_mureco_dn_rel * hadtau_syst_mureco_dn_rel + hadtau_syst_isotrk_dn_rel * hadtau_syst_isotrk_dn_rel + hadtau_syst_taumu_dn_rel * hadtau_syst_taumu_dn_rel + hadtau_syst_mistag_dn_rel * hadtau_syst_mistag_dn_rel + hadtau_syst_temp_dn_rel * hadtau_syst_temp_dn_rel + hadtau_syst_trg_dn_rel * hadtau_syst_trg_dn_rel )

      if prt_hadtau_syst_unc_dn > hadtau_rate : prt_hadtau_syst_unc_dn = hadtau_rate

      h1_hadtau_syst.SetBinError(chn, hadtau_rate*hadtau_syst_unc_avg)

      hadtau_syst_MCScale_pdf_trig_dn_rel = math.sqrt( hadtau_syst_pdf_dn_rel*hadtau_syst_pdf_dn_rel + hadtau_syst_MCScale_dn_rel*hadtau_syst_MCScale_dn_rel + hadtau_syst_trg_dn_rel*hadtau_syst_trg_dn_rel )
      hadtau_syst_MCScale_pdf_trig_up_rel = math.sqrt( hadtau_syst_pdf_up_rel*hadtau_syst_pdf_up_rel + hadtau_syst_MCScale_up_rel*hadtau_syst_MCScale_up_rel + hadtau_syst_trg_up_rel*hadtau_syst_trg_up_rel )

      hadtau_syst_closure_taumu_dn_rel = math.sqrt( hadtau_syst_closure_dn_rel*hadtau_syst_closure_dn_rel + hadtau_syst_taumu_dn_rel*hadtau_syst_taumu_dn_rel )
      hadtau_syst_closure_taumu_up_rel = math.sqrt( hadtau_syst_closure_up_rel*hadtau_syst_closure_up_rel + hadtau_syst_taumu_up_rel*hadtau_syst_taumu_up_rel )

      if hadtau_syst_closure_dn_rel >= 1: hadtau_syst_closure_dn_rel = 1-0.001
      if hadtau_syst_taumu_dn_rel == 1: hadtau_syst_taumu_dn_rel -= 0.001
      if hadtau_syst_isotrk_dn_rel == 1: hadtau_syst_isotrk_dn_rel -= 0.001
      if hadtau_syst_pdf_dn_rel >= 1: hadtau_syst_pdf_dn_rel = 1-0.001
      if hadtau_syst_MCScale_dn_rel == 1: hadtau_syst_MCScale_dn_rel -= 0.001
      if hadtau_syst_temp_dn_rel == 1: hadtau_syst_temp_dn_rel -= 0.001
      if hadtau_syst_mistag_dn_rel == 1: hadtau_syst_mistag_dn_rel -= 0.001

      if hadtau_syst_MCScale_pdf_trig_dn_rel >= 1: hadtau_syst_MCScale_pdf_trig_dn_rel = 1-0.001
      if hadtau_syst_closure_taumu_dn_rel >= 1: hadtau_syst_closure_taumu_dn_rel = 1-0.001

      if hadtau_syst_Mt_dn_rel == 1: hadtau_syst_Mt_dn_rel -= 0.001
      if hadtau_syst_muiso_dn_rel == 1: hadtau_syst_muiso_dn_rel -= 0.001
      if hadtau_syst_mureco_dn_rel == 1: hadtau_syst_mureco_dn_rel -= 0.001
      if hadtau_syst_llovr_dn_rel == 1: hadtau_syst_llovr_dn_rel -= 0.001
      if hadtau_syst_trg_dn_rel == 1: hadtau_syst_trg_dn_rel -= 0.001

      if lostle_syst_mt_dn_rel == 1: lostle_syst_mt_dn_rel -= 0.001
      if lostle_syst_muiso_dn_rel == 1: lostle_syst_muiso_dn_rel -= 0.001
      if lostle_syst_mureco_dn_rel == 1: lostle_syst_mureco_dn_rel -= 0.001
      if lostle_syst_dilep_dn_rel == 1: lostle_syst_dilep_dn_rel -= 0.001

      outfile_perchn.write("hadtau_closure_taumu_chn%d   lnN - - - %.4f/%.4f - - - - -\n" % (chn, 1-hadtau_syst_closure_taumu_dn_rel, 1+hadtau_syst_closure_taumu_up_rel))
#      outfile_perchn.write("hadtau_closure_chn%d   lnN - - - %.4f/%.4f - - - - -\n" % (chn, 1-hadtau_syst_closure_dn_rel, 1+hadtau_syst_closure_up_rel))
#      outfile_perchn.write("hadtau_taumu_chn%d   lnN - - - %.4f/%.4f - - - - -\n" % (chn, 1-hadtau_syst_taumu_dn_rel, 1+hadtau_syst_taumu_up_rel))
      outfile_perchn.write("hadtau_taumu   lnN - - - %.4f/%.4f - - - - -\n" % (1-hadtau_syst_taumu_dn_rel, 1+hadtau_syst_taumu_up_rel))
      outfile_perchn.write("hadtau_isotrk   lnN - - - %.4f/%.4f - - - - -\n" % (1-hadtau_syst_isotrk_dn_rel, 1+hadtau_syst_isotrk_up_rel))
#      outfile_perchn.write("hadtau_pdf   lnN - - - %.4f/%.4f - - - - -\n" % (1-hadtau_syst_pdf_dn_rel, 1+hadtau_syst_pdf_up_rel))
#      outfile_perchn.write("hadtau_MCScale   lnN - - - %.4f/%.4f - - - - -\n" % (1-hadtau_syst_MCScale_dn_rel, 1+hadtau_syst_MCScale_up_rel))
      outfile_perchn.write("hadtau_temp   lnN - - - %.4f/%.4f - - - - -\n" % (1-hadtau_syst_temp_dn_rel, 1+hadtau_syst_temp_up_rel))
      outfile_perchn.write("hadtau_mistag   lnN - - - %.4f/%.4f - - - - -\n" % (1-hadtau_syst_mistag_dn_rel, 1+hadtau_syst_mistag_up_rel))
#      outfile_perchn.write("hadtau_trg   lnN - - - %.4f/%.4f - - - - -\n" % (1-hadtau_syst_trg_dn_rel, 1+hadtau_syst_trg_up_rel))
      outfile_perchn.write("hadtau_MCScale_pdf_trig   lnN - - - %.4f/%.4f - - - - -\n" % (1-hadtau_syst_MCScale_pdf_trig_dn_rel, 1+hadtau_syst_MCScale_pdf_trig_up_rel))

      outfile_perchn.write("LLHadTau_mt   lnN - %.4f/%.4f - %.4f/%.4f - - - - -\n" % (1-lostle_syst_mt_dn_rel, 1+lostle_syst_mt_up_rel, 1-hadtau_syst_Mt_dn_rel, 1+hadtau_syst_Mt_up_rel))
      outfile_perchn.write("LLHadTau_muiso   lnN - %.4f/%.4f - %.4f/%.4f - - - - -\n" % (1-lostle_syst_muiso_dn_rel, 1+lostle_syst_muiso_up_rel, 1-hadtau_syst_muiso_dn_rel, 1+hadtau_syst_muiso_up_rel))
      outfile_perchn.write("LLHadTau_mureco   lnN - %.4f/%.4f - %.4f/%.4f - - - - -\n" % (1-lostle_syst_mureco_dn_rel, 1+lostle_syst_mureco_up_rel, 1-hadtau_syst_mureco_dn_rel, 1+hadtau_syst_mureco_up_rel))
      outfile_perchn.write("LLHadTau_dilep   lnN - %.4f/%.4f - %.4f/%.4f - - - - -\n" % (1-lostle_syst_dilep_dn_rel, 1+lostle_syst_dilep_up_rel, 1-hadtau_syst_llovr_dn_rel, 1+hadtau_syst_llovr_up_rel))

# zinv
      zinv_file.seek(0, 0)
      for zinv_line in zinv_file:
         zinv_splitline = procline(zinv_line)
         if zinv_splitline and zinv_splitline[0] == "syst_unc_shape_central_up":
            zinv_syst_shape_central_up_rel = float(zinv_splitline[chn])
         if zinv_splitline and zinv_splitline[0] == "syst_unc_shape_central_dn":
            zinv_syst_shape_central_dn_rel = float(zinv_splitline[chn])
      
         if zinv_splitline and zinv_splitline[0] == "syst_unc_shape_stat_up":
            zinv_syst_shape_stat_up_rel = float(zinv_splitline[chn])
         if zinv_splitline and zinv_splitline[0] == "syst_unc_shape_stat_dn":
            zinv_syst_shape_stat_dn_rel = float(zinv_splitline[chn])
      
         if zinv_splitline and zinv_splitline[0] == "syst_unc_jeu_up":
            zinv_syst_jeu_up_rel = float(zinv_splitline[chn])
         if zinv_splitline and zinv_splitline[0] == "syst_unc_jeu_dn":
            zinv_syst_jeu_dn_rel = float(zinv_splitline[chn])
      
         if zinv_splitline and zinv_splitline[0] == "syst_unc_meu_up":
            zinv_syst_meu_up_rel = float(zinv_splitline[chn])
         if zinv_splitline and zinv_splitline[0] == "syst_unc_meu_dn":
            zinv_syst_meu_dn_rel = float(zinv_splitline[chn])
      
         if zinv_splitline and zinv_splitline[0] == "syst_unc_scale_up":
            zinv_syst_scale_up_rel = float(zinv_splitline[chn])
         if zinv_splitline and zinv_splitline[0] == "syst_unc_scale_dn":
            zinv_syst_scale_dn_rel = float(zinv_splitline[chn])
      
         if zinv_splitline and zinv_splitline[0] == "syst_unc_pdf_up":
            zinv_syst_pdf_up_rel = float(zinv_splitline[chn])
         if zinv_splitline and zinv_splitline[0] == "syst_unc_pdf_dn":
            zinv_syst_pdf_dn_rel = float(zinv_splitline[chn])
      
         if zinv_splitline and zinv_splitline[0] == "syst_unc_trig_up":
            zinv_syst_trig_up_rel = float(zinv_splitline[chn])
         if zinv_splitline and zinv_splitline[0] == "syst_unc_trig_dn":
            zinv_syst_trig_dn_rel = float(zinv_splitline[chn])

         if zinv_splitline and zinv_splitline[0] == "syst_unc_norm_up":
            zinv_syst_norm_up_rel = float(zinv_splitline[chn])
         if zinv_splitline and zinv_splitline[0] == "syst_unc_norm_dn":
            zinv_syst_norm_dn_rel = float(zinv_splitline[chn])

         if zinv_splitline and zinv_splitline[0] == "syst_unc_btag_up":
            zinv_syst_btag_up_rel = float(zinv_splitline[chn])
         if zinv_splitline and zinv_splitline[0] == "syst_unc_btag_dn":
            zinv_syst_btag_dn_rel = float(zinv_splitline[chn])

         if zinv_splitline and zinv_splitline[0] == "syst_unc_bmistag_up":
            zinv_syst_bmistag_up_rel = float(zinv_splitline[chn])
         if zinv_splitline and zinv_splitline[0] == "syst_unc_bmistag_dn":
            zinv_syst_bmistag_dn_rel = float(zinv_splitline[chn])

      zinv_syst_shape_central_avg = (zinv_syst_shape_central_up_rel + zinv_syst_shape_central_dn_rel)*0.5
      zinv_syst_shape_stat_avg = (zinv_syst_shape_stat_up_rel + zinv_syst_shape_stat_dn_rel)*0.5
      zinv_syst_jeu_avg = (zinv_syst_jeu_up_rel + zinv_syst_jeu_dn_rel)*0.5
      zinv_syst_meu_avg = (zinv_syst_meu_up_rel + zinv_syst_meu_dn_rel)*0.5
      zinv_syst_scale_avg = (zinv_syst_scale_up_rel + zinv_syst_scale_dn_rel)*0.5
      zinv_syst_pdf_avg = (zinv_syst_pdf_up_rel + zinv_syst_pdf_dn_rel)*0.5
      zinv_syst_trig_avg = (zinv_syst_trig_up_rel + zinv_syst_trig_dn_rel)*0.5
      zinv_syst_norm_avg = (zinv_syst_norm_up_rel + zinv_syst_norm_dn_rel)*0.5

      zinv_syst_unc_avg = math.sqrt( zinv_syst_shape_central_avg*zinv_syst_shape_central_avg + zinv_syst_shape_stat_avg*zinv_syst_shape_stat_avg + zinv_syst_jeu_avg*zinv_syst_jeu_avg + zinv_syst_meu_avg*zinv_syst_meu_avg + zinv_syst_scale_avg*zinv_syst_scale_avg + zinv_syst_pdf_avg*zinv_syst_pdf_avg + zinv_syst_trig_avg*zinv_syst_trig_avg + zinv_syst_norm_avg * zinv_syst_norm_avg )

      prt_zinv_syst_unc_up = zinv_rate * math.sqrt( zinv_syst_shape_central_up_rel*zinv_syst_shape_central_up_rel + zinv_syst_shape_stat_up_rel*zinv_syst_shape_stat_up_rel + zinv_syst_jeu_up_rel*zinv_syst_jeu_up_rel + zinv_syst_meu_up_rel*zinv_syst_meu_up_rel + zinv_syst_scale_up_rel*zinv_syst_scale_up_rel + zinv_syst_pdf_up_rel*zinv_syst_pdf_up_rel + zinv_syst_trig_up_rel*zinv_syst_trig_up_rel + zinv_syst_norm_up_rel * zinv_syst_norm_up_rel + zinv_syst_btag_up_rel*zinv_syst_btag_up_rel + zinv_syst_bmistag_up_rel*zinv_syst_bmistag_up_rel )
      prt_zinv_syst_unc_dn = zinv_rate * math.sqrt( zinv_syst_shape_central_dn_rel*zinv_syst_shape_central_dn_rel + zinv_syst_shape_stat_dn_rel*zinv_syst_shape_stat_dn_rel + zinv_syst_jeu_dn_rel*zinv_syst_jeu_dn_rel + zinv_syst_meu_dn_rel*zinv_syst_meu_dn_rel + zinv_syst_scale_dn_rel*zinv_syst_scale_dn_rel + zinv_syst_pdf_dn_rel*zinv_syst_pdf_dn_rel + zinv_syst_trig_dn_rel*zinv_syst_trig_dn_rel + zinv_syst_norm_dn_rel * zinv_syst_norm_dn_rel + zinv_syst_btag_dn_rel*zinv_syst_btag_dn_rel + zinv_syst_bmistag_dn_rel*zinv_syst_bmistag_dn_rel )

      if prt_zinv_syst_unc_dn > zinv_rate : prt_zinv_syst_unc_dn = zinv_rate

      h1_zinv_syst.SetBinError(chn, zinv_rate*zinv_syst_unc_avg)

      zinv_syst_scale_pdf_trig_dn_rel = math.sqrt( zinv_syst_scale_dn_rel*zinv_syst_scale_dn_rel + zinv_syst_pdf_dn_rel*zinv_syst_pdf_dn_rel + zinv_syst_trig_dn_rel*zinv_syst_trig_dn_rel )
      zinv_syst_scale_pdf_trig_up_rel = math.sqrt( zinv_syst_scale_up_rel*zinv_syst_scale_up_rel + zinv_syst_pdf_up_rel*zinv_syst_pdf_up_rel + zinv_syst_trig_up_rel*zinv_syst_trig_up_rel )

      zinv_syst_norm_scale_pdf_trig_dn_rel = math.sqrt( zinv_syst_norm_dn_rel*zinv_syst_norm_dn_rel + zinv_syst_scale_pdf_trig_dn_rel*zinv_syst_scale_pdf_trig_dn_rel )
      zinv_syst_norm_scale_pdf_trig_up_rel = math.sqrt( zinv_syst_norm_up_rel*zinv_syst_norm_up_rel + zinv_syst_scale_pdf_trig_up_rel*zinv_syst_scale_pdf_trig_up_rel )

      zinv_syst_norm_scale_pdf_trig_btag_bmistag_dn_rel = math.sqrt( zinv_syst_norm_scale_pdf_trig_dn_rel*zinv_syst_norm_scale_pdf_trig_dn_rel + zinv_syst_btag_dn_rel*zinv_syst_btag_dn_rel + zinv_syst_bmistag_dn_rel*zinv_syst_bmistag_dn_rel )
      zinv_syst_norm_scale_pdf_trig_btag_bmistag_up_rel = math.sqrt( zinv_syst_norm_scale_pdf_trig_up_rel*zinv_syst_norm_scale_pdf_trig_up_rel + zinv_syst_btag_up_rel*zinv_syst_btag_up_rel + zinv_syst_bmistag_up_rel*zinv_syst_bmistag_up_rel )

      zinv_syst_jeu_meu_dn_rel = math.sqrt( zinv_syst_jeu_dn_rel*zinv_syst_jeu_dn_rel + zinv_syst_meu_dn_rel*zinv_syst_meu_dn_rel )
      zinv_syst_jeu_meu_up_rel = math.sqrt( zinv_syst_jeu_up_rel*zinv_syst_jeu_up_rel + zinv_syst_meu_up_rel*zinv_syst_meu_up_rel )
      
      if zinv_syst_shape_central_dn_rel == 1: zinv_syst_shape_central_dn_rel -= 0.001
      if zinv_syst_shape_stat_dn_rel == 1: zinv_syst_shape_stat_dn_rel -= 0.001
      if zinv_syst_jeu_dn_rel >= 1: zinv_syst_jeu_dn_rel = 1-0.001
      if zinv_syst_meu_dn_rel == 1: zinv_syst_meu_dn_rel -= 0.001
      if zinv_syst_scale_dn_rel == 1: zinv_syst_scale_dn_rel -= 0.001
      if zinv_syst_pdf_dn_rel >= 1: zinv_syst_pdf_dn_rel = 1-0.001
      if zinv_syst_trig_dn_rel == 1: zinv_syst_trig_dn_rel -= 0.001

      if zinv_syst_scale_pdf_trig_dn_rel >= 1: zinv_syst_scale_pdf_trig_dn_rel = 1-0.001
      if zinv_syst_norm_scale_pdf_trig_dn_rel >= 1: zinv_norm_syst_scale_pdf_trig_dn_rel = 1-0.001
      if zinv_syst_jeu_meu_dn_rel >= 1: zinv_syst_jeu_meu_dn_rel = 1-0.001

      if zinv_syst_norm_scale_pdf_trig_btag_bmistag_dn_rel >=1: zinv_syst_norm_scale_pdf_trig_btag_bmistag_dn_rel = 1-0.001

#      outfile_perchn.write("zinv_norm   lnN - - - - - %.4f/%.4f - - -\n" % (1-zinv_syst_norm_dn_rel, 1+zinv_syst_norm_up_rel))
      outfile_perchn.write("zinv_shape_central   lnN - - - - - %.4f/%.4f - - -\n" % (1-zinv_syst_shape_central_dn_rel, 1+zinv_syst_shape_central_up_rel))
#      outfile_perchn.write("zinv_shape_stat_chn%d   lnN - - - - - %.4f/%.4f - - -\n" % (chn, 1-zinv_syst_shape_stat_dn_rel, 1+zinv_syst_shape_stat_up_rel))
      outfile_perchn.write("zinv_shape_stat   lnN - - - - - %.4f/%.4f - - -\n" % (1-zinv_syst_shape_stat_dn_rel, 1+zinv_syst_shape_stat_up_rel))
#      outfile_perchn.write("zinv_jeu   lnN - - - - - %.4f/%.4f - - -\n" % (1-zinv_syst_jeu_dn_rel, 1+zinv_syst_jeu_up_rel))
#      outfile_perchn.write("zinv_meu   lnN - - - - - %.4f/%.4f - - -\n" % (1-zinv_syst_meu_dn_rel, 1+zinv_syst_meu_up_rel))
      outfile_perchn.write("zinv_jeu_meu   lnN - - - - - %.4f/%.4f - - -\n" % (1-zinv_syst_jeu_meu_dn_rel, 1+zinv_syst_jeu_meu_up_rel))
#      outfile_perchn.write("zinv_scale   lnN - - - - - %.4f/%.4f - - -\n" % (1-zinv_syst_scale_dn_rel, 1+zinv_syst_scale_up_rel))
#      outfile_perchn.write("zinv_pdf   lnN - - - - - %.4f/%.4f - - -\n" % (1-zinv_syst_pdf_dn_rel, 1+zinv_syst_pdf_up_rel))
#      outfile_perchn.write("zinv_trig   lnN - - - - - %.4f/%.4f - - -\n" % (1-zinv_syst_trig_dn_rel, 1+zinv_syst_trig_up_rel))
      outfile_perchn.write("zinv_syst_norm_scale_pdf_trig_btag_bmistag   lnN - - - - - %.4f/%.4f - - -\n" % (1-zinv_syst_norm_scale_pdf_trig_btag_bmistag_dn_rel, 1+zinv_syst_norm_scale_pdf_trig_btag_bmistag_up_rel))

# qcd
      qcd_file.seek(0, 0)
      for qcd_line in qcd_file:
         qcd_splitline = procline(qcd_line)
         if qcd_splitline and qcd_splitline[0] == "QCD_TFactor_relative_err":
            qcd_tfactor_err_rel = float(qcd_splitline[chn])
         if qcd_splitline and qcd_splitline[0] == "QCD_NonClosure_relative_err":
            qcd_nonclosure_err_rel = float(qcd_splitline[chn])

      qcd_syst_unc_avg = math.sqrt( qcd_tfactor_err_rel*qcd_tfactor_err_rel + qcd_nonclosure_err_rel*qcd_nonclosure_err_rel )

      prt_qcd_syst_unc_up = qcd_corr_rate * math.sqrt( qcd_tfactor_err_rel*qcd_tfactor_err_rel + qcd_nonclosure_err_rel*qcd_nonclosure_err_rel )
      prt_qcd_syst_unc_dn = qcd_corr_rate * math.sqrt( qcd_tfactor_err_rel*qcd_tfactor_err_rel + qcd_nonclosure_err_rel*qcd_nonclosure_err_rel )

      h1_qcd_syst.SetBinError(chn, qcd_corr_rate*qcd_syst_unc_avg)

      if qcd_tfactor_str in qcd_tfactor_dict:
         qcd_tfactor_idx = qcd_tfactor_dict[qcd_tfactor_str]
      else:
         qcd_tfactor_dict[qcd_tfactor_str] = qcd_tfactor_cnt
         qcd_tfactor_idx = qcd_tfactor_cnt
         qcd_tfactor_cnt += 1 

      outfile_perchn.write("qcd_tfactor_chn%d   lnN - - - - - - %.4f -  -\n" % (qcd_tfactor_idx, 1+qcd_tfactor_err_rel)) 
      outfile_perchn.write("qcd_nonclosure_chn%d   lnN - - - - - - %.4f -  -\n" % (chn, 1+qcd_nonclosure_err_rel)) 
#      outfile_perchn.write("qcd_nonclosure   lnN - - - - - - %.4f -  -\n" % (1+qcd_nonclosure_err_rel)) 
      outfile_perchn.write("ivtDphiCR_chn%d   lnU - - - - - - 10 -  -\n" % (chn))

# ttz
      ttz_file.seek(0, 0)
      for ttz_line in ttz_file:
         ttz_splitline = procline(ttz_line)
         if ttz_splitline and ttz_splitline[0] == "syst_unc_pdf_up":
            ttz_syst_pdf_up = float(ttz_splitline[chn])
            if ttz_rate !=0: ttz_syst_pdf_up_rel = ttz_syst_pdf_up/ttz_rate
         if ttz_splitline and ttz_splitline[0] == "syst_unc_pdf_down":
            ttz_syst_pdf_dn = float(ttz_splitline[chn])
            if ttz_rate !=0: ttz_syst_pdf_dn_rel = ttz_syst_pdf_dn/ttz_rate

         if ttz_splitline and ttz_splitline[0] == "syst_unc_scale_up":
            ttz_syst_scale_up = float(ttz_splitline[chn])
            if ttz_rate !=0: ttz_syst_scale_up_rel = ttz_syst_scale_up/ttz_rate
         if ttz_splitline and ttz_splitline[0] == "syst_unc_scale_down":
            ttz_syst_scale_dn = float(ttz_splitline[chn])
            if ttz_rate !=0: ttz_syst_scale_dn_rel = ttz_syst_scale_dn/ttz_rate

         if ttz_splitline and ttz_splitline[0] == "syst_unc_rate_up":
            ttz_syst_rate_up = float(ttz_splitline[chn])
            if ttz_rate !=0: ttz_syst_rate_up_rel = ttz_syst_rate_up/ttz_rate
         if ttz_splitline and ttz_splitline[0] == "syst_unc_rate_down":
            ttz_syst_rate_dn = float(ttz_splitline[chn])
            if ttz_rate !=0: ttz_syst_rate_dn_rel = ttz_syst_rate_dn/ttz_rate

      ttz_syst_pdf_avg = (ttz_syst_pdf_up + ttz_syst_pdf_dn)*0.5
      ttz_syst_scale_avg = (ttz_syst_scale_up + ttz_syst_scale_dn)*0.5
      ttz_syst_rate_avg = (ttz_syst_rate_up + ttz_syst_rate_dn)*0.5

      ttz_syst_unc_avg = math.sqrt( ttz_syst_pdf_avg*ttz_syst_pdf_avg + ttz_syst_scale_avg*ttz_syst_scale_avg + ttz_syst_rate_avg*ttz_syst_rate_avg )

      prt_ttz_syst_unc_up = math.sqrt( ttz_syst_pdf_up*ttz_syst_pdf_up + ttz_syst_scale_up*ttz_syst_scale_up + ttz_syst_rate_up*ttz_syst_rate_up )
      prt_ttz_syst_unc_dn = math.sqrt( ttz_syst_pdf_dn*ttz_syst_pdf_dn + ttz_syst_scale_dn*ttz_syst_scale_dn + ttz_syst_rate_dn*ttz_syst_rate_dn )
      if prt_ttz_syst_unc_dn > ttz_rate : prt_ttz_syst_unc_dn = ttz_rate

      h1_ttz_syst.SetBinError(chn, ttz_syst_unc_avg)

      ttz_syst_pdf_scale_dn_rel = math.sqrt( ttz_syst_pdf_dn_rel*ttz_syst_pdf_dn_rel + ttz_syst_scale_dn_rel*ttz_syst_scale_dn_rel )
      ttz_syst_pdf_scale_up_rel = math.sqrt( ttz_syst_pdf_up_rel*ttz_syst_pdf_up_rel + ttz_syst_scale_up_rel*ttz_syst_scale_up_rel )

      ttz_syst_pdf_scale_avg_rel = (ttz_syst_pdf_scale_dn_rel + ttz_syst_pdf_scale_up_rel)*0.5

      if ttz_syst_pdf_dn_rel == 1: ttz_syst_pdf_dn_rel -= 0.001
      if ttz_syst_scale_dn_rel == 1: ttz_syst_scale_dn_rel -= 0.001
      if ttz_syst_rate_dn_rel == 1: ttz_syst_rate_dn_rel -= 0.001

      if ttz_syst_pdf_scale_dn_rel == 1: ttz_syst_pdf_scale_dn_rel -= 0.001

#      outfile_perchn.write("ttz_pdf   lnN - - - - - - - %.4f/%.4f -\n" % (1-ttz_syst_pdf_dn_rel, 1+ttz_syst_pdf_up_rel)) 
#      outfile_perchn.write("ttz_scale   lnN - - - - - - - %.4f/%.4f -\n" % (1-ttz_syst_scale_dn_rel, 1+ttz_syst_scale_up_rel)) 
#      outfile_perchn.write("ttz_pdf_scale   lnN - - - - - - - %.4f/%.4f -\n" % (1-ttz_syst_pdf_scale_dn_rel, 1+ttz_syst_pdf_scale_up_rel)) 
      outfile_perchn.write("ttz_pdf_scale   lnN - - - - - - - %.4f -\n" % (1+ttz_syst_pdf_scale_avg_rel)) 
      outfile_perchn.write("ttz_rate   lnN - - - - - - - %.4f/%.4f -\n" % (1-ttz_syst_rate_dn_rel, 1+ttz_syst_rate_up_rel)) 

# rare
      rare_file.seek(0, 0)
      for rare_line in rare_file:
         rare_splitline = procline(rare_line)
         if rare_splitline and rare_splitline[0] == "syst_unc_pdf_up":
            rare_syst_pdf_up = float(rare_splitline[chn])
            if rare_rate !=0: rare_syst_pdf_up_rel = rare_syst_pdf_up/rare_rate
         if rare_splitline and rare_splitline[0] == "syst_unc_pdf_down":
            rare_syst_pdf_dn = float(rare_splitline[chn])
            if rare_rate !=0: rare_syst_pdf_dn_rel = rare_syst_pdf_dn/rare_rate

         if rare_splitline and rare_splitline[0] == "syst_unc_scale_up":
            rare_syst_scale_up = float(rare_splitline[chn])
            if rare_rate !=0: rare_syst_scale_up_rel = rare_syst_scale_up/rare_rate
         if rare_splitline and rare_splitline[0] == "syst_unc_scale_down":
            rare_syst_scale_dn = float(rare_splitline[chn])
            if rare_rate !=0: rare_syst_scale_dn_rel = rare_syst_scale_dn/rare_rate

      rare_syst_pdf_avg = (rare_syst_pdf_up + rare_syst_pdf_dn)*0.5
      rare_syst_scale_avg = (rare_syst_scale_up + rare_syst_scale_dn)*0.5

      rare_syst_unc_avg = math.sqrt( rare_syst_pdf_avg*rare_syst_pdf_avg + rare_syst_scale_avg*rare_syst_scale_avg )

      prt_rare_syst_unc_up = math.sqrt( rare_syst_pdf_up*rare_syst_pdf_up + rare_syst_scale_up*rare_syst_scale_up )
      prt_rare_syst_unc_dn = math.sqrt( rare_syst_pdf_dn*rare_syst_pdf_dn + rare_syst_scale_dn*rare_syst_scale_dn )
      if prt_rare_syst_unc_dn > rare_rate : prt_rare_syst_unc_dn = rare_rate

      h1_rare_syst.SetBinError(chn, rare_syst_unc_avg)

      rare_syst_pdf_scale_dn_rel = math.sqrt( rare_syst_pdf_dn_rel*rare_syst_pdf_dn_rel + rare_syst_scale_dn_rel*rare_syst_scale_dn_rel )
      rare_syst_pdf_scale_up_rel = math.sqrt( rare_syst_pdf_up_rel*rare_syst_pdf_up_rel + rare_syst_scale_up_rel*rare_syst_scale_up_rel )

      rare_syst_pdf_scale_avg_rel = (rare_syst_pdf_scale_dn_rel + rare_syst_pdf_scale_up_rel)*0.5

      if rare_syst_pdf_dn_rel == 1: rare_syst_pdf_dn_rel -= 0.001
      if rare_syst_scale_dn_rel == 1: rare_syst_scale_dn_rel -= 0.001

      if rare_syst_pdf_scale_dn_rel == 1: rare_syst_pdf_scale_dn_rel -= 0.001

      outfile_perchn.write("rare_pdf_scale   lnN - - - - - - - - %.5f\n" % (1+rare_syst_pdf_scale_avg_rel)) 

      outfile_perchn.close()
# Modify the data card to make only expected limits
      if doExpLimitOnly:
         if len(outputdir) !=0:
            outfile_perchn_ori = open(outputdir + "/" + outbase_filename + "_ch" + str(chn) + ".txt", "r")
            outfile_perchn_modif = open(outputdir + "/" + outbase_filename + "_ch" + str(chn) + ".tmp.txt", "w")
         else:
            outfile_perchn_ori = open(outbase_filename + "_ch" + str(chn) + ".txt", "r")
            outfile_perchn_modif = open(outbase_filename + "_ch" + str(chn) + ".tmp.txt", "w")
   
         for outfile_line in outfile_perchn_ori:
            if "observation" in outfile_line:
               outfile_perchn_modif.write("observation %.4f\n" % (pred_tot_rate))
            else:
               outfile_perchn_modif.write(outfile_line)
   
         outfile_perchn_ori.close()
         outfile_perchn_modif.close()
         if len(outputdir) !=0:
            command_line = "mv "+outputdir + "/" + outbase_filename + "_ch" + str(chn) + ".tmp.txt "+outputdir + "/" + outbase_filename + "_ch" + str(chn) + ".txt"
         else:
            command_line = "mv "+outbase_filename + "_ch" + str(chn) + ".tmp.txt "+outbase_filename + "_ch" + str(chn) + ".txt"
         os.system(command_line)

# Additional output file for lostle and hadtau
      if len(outputdir) !=0:
         LLHadTau_stat_HighW_file_perchn = open(outputdir + "/comb_SLControl_"+ signal_key + "_ch" + str(chn) + ".txt", "w")
      else:
         LLHadTau_stat_HighW_file_perchn = open("comb_SLControl_"+ signal_key + "_ch" + str(chn) + ".txt", "w")
      LLHadTau_stat_HighW_file_perchn.write("imax 1 # number of channels\n")
      LLHadTau_stat_HighW_file_perchn.write("jmax 4 # number of backgrounds\n")
      LLHadTau_stat_HighW_file_perchn.write("kmax * nuissance\n")
      LLHadTau_stat_HighW_file_perchn.write("shapes * * FAKE\n")
      LLHadTau_stat_HighW_file_perchn.write("----------------\n")
      LLHadTau_stat_HighW_file_perchn.write("bin binSLControl%d\n" % (chn))
      LLHadTau_stat_HighW_file_perchn.write("observation 0.0\n")

      LLHadTau_stat_HighW_file_perchn.write("bin       ")
      for ibkg in range(5):
         LLHadTau_stat_HighW_file_perchn.write("binSLControl%d "% (chn))
      LLHadTau_stat_HighW_file_perchn.write("\n")

      LLHadTau_stat_HighW_file_perchn.write("process   Sig   LostLep  HadTau LostLepHighW  HadTauHighW\n")
      LLHadTau_stat_HighW_file_perchn.write("process     0         1       2            3            4\n")
#      LLHadTau_stat_HighW_file_perchn.write("rate   0.0001    0.0001  0.0001      0.50000      0.50000\n")
      LLHadTau_stat_HighW_file_perchn.write("rate   0.0001    0.0001  0.0001      1.00000      1.00000\n")
      LLHadTau_stat_HighW_file_perchn.write("---------------------------------\n")
      LLHadTau_stat_HighW_file_perchn.write("stat_unc_HighW_chn%d    lnU - - - 10 10\n" % (chn))

# Additional output file for QCD
      if len(outputdir) !=0:
         qcd_outfile_perchn = open(outputdir + "/comb_invertDphi_"+ signal_key + "_ch" + str(chn) + ".txt", "w")
      else:
         qcd_outfile_perchn = open("comb_invertDphi_"+ signal_key + "_ch" + str(chn) + ".txt", "w")

      qcd_file.seek(0, 0)
      for qcd_line in qcd_file:
         qcd_splitline = procline(qcd_line)
         if qcd_splitline and qcd_splitline[0] == "QCD_TFactor_err":
            qcd_tfactor_err = float(qcd_splitline[chn])
         if qcd_splitline:
            if qcd_splitline[0] == "QCD_otherBG_CS":
               qcd_contam_pred = float(qcd_splitline[chn])

            if qcd_splitline[0] == "QCD_otherBG_CS_relative_errup":
               qcd_contam_unc_up_rel = float(qcd_splitline[chn])
            if qcd_splitline[0] == "QCD_otherBG_CS_relative_errdown":
               qcd_contam_unc_dn_rel = float(qcd_splitline[chn])

      qcd_outfile_perchn.write("imax 1 # number of channels\n")
      qcd_outfile_perchn.write("jmax 2 # number of backgrounds\n")
      qcd_outfile_perchn.write("kmax * nuissance\n")
      qcd_outfile_perchn.write("shapes * * FAKE\n")
      qcd_outfile_perchn.write("----------------\n")
      qcd_outfile_perchn.write("bin binInvertDphi%d\n" % (chn))
      qcd_outfile_perchn.write("observation %0.0f\n" % (qcd_cs))

      qcd_outfile_perchn.write("bin       ")
      for ibkg in range(3):
         qcd_outfile_perchn.write("binInvertDphi%d "% (chn))
      qcd_outfile_perchn.write("\n")

      qcd_outfile_perchn.write("process   Sig  QCD  Contam\n")
      qcd_outfile_perchn.write("process     0   1      2\n")
      qcd_outfile_perchn.write("rate      0.0001")
      if qcd_cs != 0:
         qcd_outfile_perchn.write("  %.4f  %.4f\n" % (qcd_cs, qcd_contam_pred))
      else:
         qcd_outfile_perchn.write("  %.4f  %.4f\n" % (1.0, qcd_contam_pred))
      qcd_outfile_perchn.write("---------------------------------\n")
      
      if qcd_contam_unc_dn_rel >= 1: qcd_contam_unc_dn_rel = 1-0.001

      qcd_outfile_perchn.write("ivtDphiCR_chn%d  lnU - 10 -\n" % (chn))
      qcd_outfile_perchn.write("contamUnc_chn%d  lnN - - %.4f/%.4f\n" % (chn, 1-qcd_contam_unc_dn_rel, 1+qcd_contam_unc_up_rel))

      qcd_outfile_perchn.close() 

      prt_qcd_syst_unc_up = math.sqrt( prt_qcd_syst_unc_up * prt_qcd_syst_unc_up + qcd_contam_unc_up_rel*qcd_contam_pred*qcd_tfactor*qcd_contam_unc_up_rel*qcd_contam_pred*qcd_tfactor)
      prt_qcd_syst_unc_dn = math.sqrt( prt_qcd_syst_unc_dn * prt_qcd_syst_unc_dn + qcd_contam_unc_dn_rel*qcd_contam_pred*qcd_tfactor*qcd_contam_unc_dn_rel*qcd_contam_pred*qcd_tfactor)
      if prt_qcd_syst_unc_dn > qcd_corr_rate : prt_qcd_syst_unc_dn = qcd_corr_rate

#      print "chn%d  data_rate : %s +- %.4f <-->  sumPred : %.4f +- %.4f" % (chn-1, data_rate, math.sqrt(float(data_rate)), pred_tot_rate, pred_tot_stat)
# Print out results related to event yields
      prt_pred_tot_rate = pred_tot_rate

      prt_pred_tot_stat_up = math.sqrt( prt_qcd_stat_unc_up * prt_qcd_stat_unc_up + prt_lostle_stat_unc_up * prt_lostle_stat_unc_up + prt_hadtau_stat_unc_up * prt_hadtau_stat_unc_up + prt_zinv_stat_unc_up * prt_zinv_stat_unc_up + prt_ttz_stat_unc_up * prt_ttz_stat_unc_up + prt_rare_stat_unc_up * prt_rare_stat_unc_up )
      prt_pred_tot_stat_dn = math.sqrt( prt_qcd_stat_unc_dn * prt_qcd_stat_unc_dn + prt_lostle_stat_unc_dn * prt_lostle_stat_unc_dn + prt_hadtau_stat_unc_dn * prt_hadtau_stat_unc_dn + prt_zinv_stat_unc_dn * prt_zinv_stat_unc_dn + prt_ttz_stat_unc_dn * prt_ttz_stat_unc_dn + prt_rare_stat_unc_dn * prt_rare_stat_unc_dn )

      prt_pred_tot_syst_up = math.sqrt( prt_qcd_syst_unc_up * prt_qcd_syst_unc_up + prt_lostle_syst_unc_up * prt_lostle_syst_unc_up + prt_hadtau_syst_unc_up * prt_hadtau_syst_unc_up + prt_zinv_syst_unc_up * prt_zinv_syst_unc_up + prt_ttz_syst_unc_up * prt_ttz_syst_unc_up + prt_rare_syst_unc_up * prt_rare_syst_unc_up )
      prt_pred_tot_syst_dn = math.sqrt( prt_qcd_syst_unc_dn * prt_qcd_syst_unc_dn + prt_lostle_syst_unc_dn * prt_lostle_syst_unc_dn + prt_hadtau_syst_unc_dn * prt_hadtau_syst_unc_dn + prt_zinv_syst_unc_dn * prt_zinv_syst_unc_dn + prt_ttz_syst_unc_dn * prt_ttz_syst_unc_dn + prt_rare_syst_unc_dn * prt_rare_syst_unc_dn )

      glb_data_rate.append(float(data_rate))
      glb_lostle_rate.append(lostle_rate)
      glb_hadtau_rate.append(hadtau_rate)
      glb_zinv_rate.append(zinv_rate)
      glb_qcd_rate.append(qcd_corr_rate)
      glb_ttz_rate.append(ttz_rate)
      glb_rare_rate.append(rare_rate)

      glb_tot_rate.append(pred_tot_rate)

      glb_lostle_stat_up.append(prt_lostle_stat_unc_up)
      glb_hadtau_stat_up.append(prt_hadtau_stat_unc_up)
      glb_zinv_stat_up.append(prt_zinv_stat_unc_up)
      glb_qcd_stat_up.append(prt_qcd_stat_unc_up)
      glb_ttz_stat_up.append(prt_ttz_stat_unc_up)
      glb_rare_stat_up.append(prt_rare_stat_unc_up)

      glb_lostle_stat_dn.append(prt_lostle_stat_unc_dn)
      glb_hadtau_stat_dn.append(prt_hadtau_stat_unc_dn)
      glb_zinv_stat_dn.append(prt_zinv_stat_unc_dn)
      glb_qcd_stat_dn.append(prt_qcd_stat_unc_dn)
      glb_ttz_stat_dn.append(prt_ttz_stat_unc_dn)
      glb_rare_stat_dn.append(prt_rare_stat_unc_dn)

      glb_tot_stat_up.append(prt_pred_tot_stat_up)
      glb_tot_stat_dn.append(prt_pred_tot_stat_dn)

      glb_lostle_syst_up.append(prt_lostle_syst_unc_up)
      glb_hadtau_syst_up.append(prt_hadtau_syst_unc_up)
      glb_zinv_syst_up.append(prt_zinv_syst_unc_up)
      glb_qcd_syst_up.append(prt_qcd_syst_unc_up)
      glb_ttz_syst_up.append(prt_ttz_syst_unc_up)
      glb_rare_syst_up.append(prt_rare_syst_unc_up)

      glb_lostle_syst_dn.append(prt_lostle_syst_unc_dn)
      glb_hadtau_syst_dn.append(prt_hadtau_syst_unc_dn)
      glb_zinv_syst_dn.append(prt_zinv_syst_unc_dn)
      glb_qcd_syst_dn.append(prt_qcd_syst_unc_dn)
      glb_ttz_syst_dn.append(prt_ttz_syst_unc_dn)
      glb_rare_syst_dn.append(prt_rare_syst_unc_dn)

      glb_tot_syst_up.append(prt_pred_tot_syst_up)
      glb_tot_syst_dn.append(prt_pred_tot_syst_dn)

      print "chn%d   data_rate : %s +- %.4f   sumPred : %.4f + %.4f - %.4f (+ %.4f - %.4f)" % (chn-1, data_rate, math.sqrt(float(data_rate)), prt_pred_tot_rate, prt_pred_tot_stat_up, prt_pred_tot_stat_dn, prt_pred_tot_syst_up, prt_pred_tot_syst_dn)
# Easy parsing format
      prt_table_file.write("%d  %s %.4f  %.4f %.4f %.4f %.4f %.4f    %.4f %.4f %.4f %.4f %.4f    %.4f %.4f %.4f %.4f %.4f    %.4f %.4f %.4f %.4f %.4f    %.4f %.4f %.4f %.4f %.4f    %.4f %.4f %.4f %.4f %.4f     %.5f %.5f %.5f %.5f %.5f\n" % (chn-1, data_rate, math.sqrt(float(data_rate)), prt_pred_tot_rate, prt_pred_tot_stat_up, prt_pred_tot_stat_dn, prt_pred_tot_syst_up, prt_pred_tot_syst_dn,     lostle_rate, prt_lostle_stat_unc_up, prt_lostle_stat_unc_dn, prt_lostle_syst_unc_up, prt_lostle_syst_unc_dn,     hadtau_rate, prt_hadtau_stat_unc_up, prt_hadtau_stat_unc_dn, prt_hadtau_syst_unc_up, prt_hadtau_syst_unc_dn,     zinv_rate, prt_zinv_stat_unc_up, prt_zinv_stat_unc_dn, prt_zinv_syst_unc_up, prt_zinv_syst_unc_dn,     qcd_corr_rate, prt_qcd_stat_unc_up, prt_qcd_stat_unc_dn, prt_qcd_syst_unc_up, prt_qcd_syst_unc_dn,     ttz_rate, prt_ttz_stat_unc_up, prt_ttz_stat_unc_dn, prt_ttz_syst_unc_up, prt_ttz_syst_unc_dn,     rare_rate, prt_rare_stat_unc_up, prt_rare_stat_unc_dn, prt_rare_syst_unc_up, prt_rare_syst_unc_dn ))
#      prt_table_file.write("%d  %s %.4f  %.4f %.4f %.4f %.4f %.4f\n" % (chn-1, data_rate, math.sqrt(float(data_rate)), prt_pred_tot_rate, prt_pred_tot_stat_up, prt_pred_tot_stat_dn, prt_pred_tot_syst_up, prt_pred_tot_syst_dn))

   prt_table_file.close()

   rtfile.cd()
   h1_data.Write()
   h1_lostle.Write()
   h1_hadtau.Write()
   h1_zinv.Write()
   h1_qcd.Write()
   h1_ttz.Write()
   h1_rare.Write()
   h1_lostle_syst.Write()
   h1_hadtau_syst.Write()
   h1_zinv_syst.Write()
   h1_qcd_syst.Write()
   h1_ttz_syst.Write()
   h1_rare_syst.Write()
   rtfile.Close()

def main():
   usage = "usage: %prog options"
   version = "%prog."
   parser = OptionParser(usage=usage,version=version)
   parser.add_option("-l", "--lostle", action="store", dest="lostle", type="string", default="lostle.txt", help="set lostle data card name")
   parser.add_option("-t", "--hadtau", action="store", dest="hadtau", type="string", default="hadtau.txt", help="set hadtau data card name")
   parser.add_option("-z", "--zinv", action="store", dest="zinv", type="string", default="zinv.txt", help="set zinv data card name")
   parser.add_option("-q", "--qcd", action="store", dest="qcd", type="string", default="qcd.txt", help="set qcd data card name")
   parser.add_option("-r", "--ttz", action="store", dest="ttz", type="string", default="ttz.txt", help="set ttz data card name")
   parser.add_option("-b", "--rare", action="store", dest="rare", type="string", default="rare.txt", help="set rare data card name")
   parser.add_option("-d", "--data", action="store", dest="data", type="string", default="data.txt", help="set data data card name")
   parser.add_option("-s", "--signal", action="store", dest="signal", type="string", default="signal.txt", help="set signal data card name")
   parser.add_option("-o", "--outputdir", action="store", dest="outputdir", type="string", default="", help="set combined card output directory")
   
   (options, args) = parser.parse_args()
   
   print 'lostle :', options.lostle
   print 'hadtau :', options.hadtau
   print 'zinv :', options.zinv
   print 'qcd :', options.qcd
   print 'ttz :', options.ttz
   print 'rare :', options.rare
   print 'data :', options.data
   print 'signal : ', options.signal
   print 'outputdir : ', options.outputdir

   if len(options.outputdir) !=0:
      if not os.path.exists(options.outputdir): os.mkdir(options.outputdir) 

   splitsignalinput = options.signal.split("/")
   stripDirSignalInput = splitsignalinput[-1]

   if "signal_" in stripDirSignalInput: tmp_signal_key = stripDirSignalInput.replace("signal_", "")
   elif "signal" in stripDirSignalInput: tmp_signal_key = stripDirSignalInput.replace("signal", "")
   tmp_signal_key = tmp_signal_key.replace(".txt", "")
   
   lostle_file = open(options.lostle)
   hadtau_file = open(options.hadtau)
   zinv_file = open(options.zinv)
   qcd_file = open(options.qcd)
   ttz_file = open(options.ttz)
   rare_file = open(options.rare)
   data_file = open(options.data)
   signal_file = open(options.signal)
   
   prodCardPerChn(tmp_signal_key, options.outputdir, lostle_file, hadtau_file, zinv_file, qcd_file, ttz_file, rare_file, data_file, signal_file)

if __name__ == "__main__":
   main()
