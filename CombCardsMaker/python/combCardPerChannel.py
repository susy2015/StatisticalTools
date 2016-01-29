#!/usr/bin/python

#from os import listdir
#from os import sys
#from os import system
#from os import path
import os
import math
from optparse import OptionParser # Command line parsing

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

def main():
   usage = "usage: %prog options"
   version = "%prog."
   parser = OptionParser(usage=usage,version=version)
   parser.add_option("-l", "--lostle", action="store", dest="lostle", type="string", default="lostle.txt", help="set lostle data card name")
   parser.add_option("-t", "--hadtau", action="store", dest="hadtau", type="string", default="hadtau.txt", help="set hadtau data card name")
   parser.add_option("-z", "--zinv", action="store", dest="zinv", type="string", default="zinv.txt", help="set zinv data card name")
   parser.add_option("-q", "--qcd", action="store", dest="qcd", type="string", default="qcd.txt", help="set qcd data card name")
   parser.add_option("-r", "--ttz", action="store", dest="ttz", type="string", default="ttz.txt", help="set ttz data card name")
   parser.add_option("-d", "--data", action="store", dest="data", type="string", default="data.txt", help="set data data card name")
   parser.add_option("-s", "--signal", action="store", dest="signal", type="string", default="signal.txt", help="set signal data card name")
   parser.add_option("-o", "--outputdir", action="store", dest="outputdir", type="string", default="", help="set combined card output directory")
   
   (options, args) = parser.parse_args()
   
   print 'lostle :', options.lostle
   print 'hadtau :', options.hadtau
   print 'zinv :', options.zinv
   print 'qcd :', options.qcd
   print 'ttz :', options.ttz
   print 'data :', options.data
   print 'signal : ', options.signal
   print 'outputdir : ', options.outputdir

   if len(options.outputdir) !=0:
      if not os.path.exists(options.outputdir): os.mkdir(options.outputdir) 

   if "signal_" in options.signal: tmp_signal_key = options.signal.replace("signal_", "")
   elif "signal" in options.signal: tmp_signal_key = options.signal.replace("signal", "")
   tmp_signal_key = tmp_signal_key.replace(".txt", "")
   outbase_filename = "comb_" + tmp_signal_key
   print '\noutbase_filename : %s\n' % (outbase_filename)
   
   lostle_file = open(options.lostle)
   hadtau_file = open(options.hadtau)
   zinv_file = open(options.zinv)
   qcd_file = open(options.qcd)
   ttz_file = open(options.ttz)
   data_file = open(options.data)
   signal_file = open(options.signal)
   
   # parsing data file first to get basic information
   for line in data_file:
      splitline = procline(line)
      if splitline: 
#         print splitline
         if splitline[0] == "luminosity" : lumi = splitline[1]
         if splitline[0] == "channels" : channels = splitline[1]

   for chn in range(1, int(channels)+1):
      if len(options.outputdir) !=0:
         outfile_perchn = open(options.outputdir + "/" + outbase_filename + "_ch" + str(chn) + ".txt", "w")
      else:
         outfile_perchn = open(outbase_filename + "_ch" + str(chn) + ".txt", "w")
      outfile_perchn.write("imax 1 # number of channels\n")
      outfile_perchn.write("jmax 7 # number of backgrounds\n")
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

      outfile_perchn.write("bin       ")
      for ibkg in range(8):
         outfile_perchn.write("bin%d "% (chn))
      outfile_perchn.write("\n")

      outfile_perchn.write("process   Sig  LostLep  LostLepHighW  HadTau  HadTauHighW  Zinv  QCD  TTZ\n")
      outfile_perchn.write("process     0        1            2       3           4     5    6    7\n")
      outfile_perchn.write("rate      ")

      signal_file.seek(0, 0)
      for signal_line in signal_file:
         signal_splitline = procline(signal_line)
         if signal_splitline and signal_splitline[0] == "rate":
            signal_rate = float(signal_splitline[chn])
      outfile_perchn.write("%.4f  " % (signal_rate))

      pred_tot_rate = 0
 
      lostle_file.seek(0, 0)
      for lostle_line in lostle_file:
         lostle_splitline = procline(lostle_line)
         if lostle_splitline and lostle_splitline[0] == "rate":
            lostle_rate = float(lostle_splitline[chn])
            if lostle_rate == 0: lostle_rate = 0.0001
      outfile_perchn.write("%.4f  " % (lostle_rate))
      outfile_perchn.write("0.0001  ")
      pred_tot_rate += lostle_rate
 
      hadtau_file.seek(0, 0)
      for hadtau_line in hadtau_file:
         hadtau_splitline = procline(hadtau_line)
         if hadtau_splitline and hadtau_splitline[0] == "rate":
            hadtau_rate = float(hadtau_splitline[chn])
            if hadtau_rate ==0: hadtau_rate = 0.0001
      outfile_perchn.write("%.4f  " % (hadtau_rate))
      outfile_perchn.write("0.0001  ")
      pred_tot_rate += hadtau_rate
 
      zinv_file.seek(0, 0)
      for zinv_line in zinv_file:
         zinv_splitline = procline(zinv_line)
         if zinv_splitline and zinv_splitline[0] == "rate":
            zinv_rate = float(zinv_splitline[chn])
            if zinv_rate ==0: zinv_rate = 0.0001
      outfile_perchn.write("%.4f  " % (zinv_rate))
      pred_tot_rate += zinv_rate
 
      qcd_file.seek(0, 0)
      for qcd_line in qcd_file:
         qcd_splitline = procline(qcd_line)
         if qcd_splitline and qcd_splitline[0] == "QCD_Data_CS":
            qcd_cs = float(qcd_splitline[chn])
         if qcd_splitline and qcd_splitline[0] == "QCD_TFactor":
            qcd_tfactor = float(qcd_splitline[chn])
         if qcd_splitline and qcd_splitline[0] == "QCD_otherBG_CS":
            tmp_qcd_contam_pred = float(qcd_splitline[chn])
      if qcd_cs != 0:
         outfile_perchn.write("%.4f  " % (qcd_cs*qcd_tfactor))
      else:
         outfile_perchn.write("%.4f  " % qcd_tfactor)
      pred_tot_rate += (qcd_cs-tmp_qcd_contam_pred)*qcd_tfactor
 
      ttz_file.seek(0, 0)
      for ttz_line in ttz_file:
         ttz_splitline = procline(ttz_line)
         if ttz_splitline and ttz_splitline[0] == "rate":
            ttz_rate = float(ttz_splitline[chn])
      outfile_perchn.write("%.4f  " % (ttz_rate))
      outfile_perchn.write("\n")
      pred_tot_rate += ttz_rate

      outfile_perchn.write("---------------------------------------\n")
# signal stat. unc.
      signal_file.seek(0, 0)
      for signal_line in signal_file:
         signal_splitline = procline(signal_line)
         if signal_splitline and signal_splitline[0] == "cs_event":
            signal_cs_event = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "avg_weight":
            signal_avg_weight = float(signal_splitline[chn])
      outfile_perchn.write("signal_stat_unc_chn%d   gmN %0.0f %0.8f - - - - - - - \n" % (chn, signal_cs_event, signal_avg_weight))

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
         if lostle_splitline and lostle_splitline[0] == "avg_weight":
            lostle_avg_weight = float(lostle_splitline[chn])
 
      hadtau_file.seek(0, 0)
      for hadtau_line in hadtau_file:
         hadtau_splitline = procline(hadtau_line)
         if hadtau_splitline and hadtau_splitline[0] == "stat_unc_up":
            hadtau_stat_unc_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "stat_unc_dn":
            hadtau_stat_unc_dn_rel = float(hadtau_splitline[chn])

      if lostle_stat_unc_dn_rel ==1: lostle_stat_unc_dn_rel -= 0.001
      if hadtau_stat_unc_dn_rel ==1: hadtau_stat_unc_dn_rel -= 0.001
 
      outfile_perchn.write("muonCS_stat_unc_chn%d    lnN    - %.4f/%.4f - %.4f/%.4f - - - -\n" % (chn, 1-lostle_stat_unc_dn_rel, 1+lostle_stat_unc_up_rel, 1-hadtau_stat_unc_dn_rel, 1+hadtau_stat_unc_up_rel))
      outfile_perchn.write("stat_unc_HighW_chn%d     gmN  0 - - %.4f - %.4f - - -\n" % (chn, lostle_avg_weight, 0.25))

# zinv stat. unc.
      zinv_file.seek(0, 0)
      for zinv_line in zinv_file:
         zinv_splitline = procline(zinv_line)
         if zinv_splitline and zinv_splitline[0] == "cs_event":
            zinv_cs_event = float(zinv_splitline[chn])
         if zinv_splitline and zinv_splitline[0] == "avg_weight":
            zinv_avg_weight = float(zinv_splitline[chn])

      outfile_perchn.write("zinv_stat_unc_chn%d  gmN %0.0f - - - - - %0.4f - - \n" % (chn, zinv_cs_event, zinv_avg_weight))

# ttz stat. unc.
      ttz_file.seek(0, 0)
      for ttz_line in ttz_file:
         ttz_splitline = procline(ttz_line)
         if ttz_splitline and ttz_splitline[0] == "cs_event":
            ttz_cs_event = float(ttz_splitline[chn])
         if ttz_splitline and ttz_splitline[0] == "avg_weight":
            ttz_avg_weight = float(ttz_splitline[chn])

      outfile_perchn.write("ttz_stat_unc_chn%d  gmN %0.0f - - - - - - - %0.4f\n" % (chn, ttz_cs_event, ttz_avg_weight))

# Systematic uncertainties
      outfile_perchn.write("-------------------------------------------\n")
# signal
      signal_file.seek(0, 0)
      for signal_line in signal_file:
         signal_splitline = procline(signal_line)
         if signal_splitline and signal_splitline[0] == "syst_unc_up":
            signal_syst_up_rel = float(signal_splitline[chn])
         if signal_splitline and signal_splitline[0] == "syst_unc_dn":
            signal_syst_dn_rel = float(signal_splitline[chn])
      if signal_syst_dn_rel ==1: signal_syst_dn_rel -= 0.001
      outfile_perchn.write("signal_syst  lnN %.4f/%.4f - - - - - - -\n" % (1-signal_syst_dn_rel, 1+signal_syst_up_rel))

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

      if lostle_syst_closure_dn_rel ==1: lostle_syst_closure_dn_rel -= 0.001
      if lostle_syst_purity_dn_rel ==1: lostle_syst_purity_dn_rel -= 0.001
      if lostle_syst_eiso_dn_rel ==1: lostle_syst_eiso_dn_rel -= 0.001
      if lostle_syst_ereco_dn_rel ==1: lostle_syst_ereco_dn_rel -= 0.001
      if lostle_syst_isotrk_dn_rel ==1: lostle_syst_isotrk_dn_rel -= 0.001
      if lostle_syst_acc_dn_rel ==1: lostle_syst_acc_dn_rel -= 0.001

      outfile_perchn.write("lostle_closure_chn%d   lnN - %.4f/%.4f - - - - - -\n" % (chn, 1-lostle_syst_closure_dn_rel, 1+lostle_syst_closure_up_rel))
      outfile_perchn.write("lostle_purity_chn%d   lnN - %.4f/%.4f - - - - - -\n" % (chn, 1-lostle_syst_purity_dn_rel, 1+lostle_syst_purity_up_rel))
      outfile_perchn.write("lostle_eiso  lnN - %.4f/%.4f - - - - - -\n" % (1-lostle_syst_eiso_dn_rel, 1+lostle_syst_eiso_up_rel))
      outfile_perchn.write("lostle_ereco  lnN - %.4f/%.4f - - - - - -\n" % (1-lostle_syst_ereco_dn_rel, 1+lostle_syst_ereco_up_rel))
      outfile_perchn.write("lostle_isotrk  lnN - %.4f/%.4f - - - - - -\n" % (1-lostle_syst_isotrk_dn_rel, 1+lostle_syst_isotrk_up_rel))
      outfile_perchn.write("lostle_acc  lnN - %.4f/%.4f - - - - - -\n" % (1-lostle_syst_acc_dn_rel, 1+lostle_syst_acc_up_rel))

# hadtau
      hadtau_file.seek(0, 0)
      for hadtau_line in hadtau_file:
         hadtau_splitline = procline(hadtau_line)
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_pdf_up":
            hadtau_syst_pdf_up_rel = float(hadtau_splitline[chn])
         if hadtau_splitline and hadtau_splitline[0] == "syst_unc_pdf_dn":
            hadtau_syst_pdf_dn_rel = float(hadtau_splitline[chn])

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

      if hadtau_syst_closure_dn_rel == 1: hadtau_syst_closure_dn_rel -= 0.001
      if hadtau_syst_taumu_dn_rel == 1: hadtau_syst_taumu_dn_rel -= 0.001
      if hadtau_syst_isotrk_dn_rel == 1: hadtau_syst_isotrk_dn_rel -= 0.001
      if hadtau_syst_pdf_dn_rel == 1: hadtau_syst_pdf_dn_rel -= 0.001
      if hadtau_syst_temp_dn_rel == 1: hadtau_syst_temp_dn_rel -= 0.001
      if hadtau_syst_mistag_dn_rel == 1: hadtau_syst_mistag_dn_rel -= 0.001

      if hadtau_syst_Mt_dn_rel == 1: hadtau_syst_Mt_dn_rel -= 0.001
      if hadtau_syst_muiso_dn_rel == 1: hadtau_syst_muiso_dn_rel -= 0.001
      if hadtau_syst_mureco_dn_rel == 1: hadtau_syst_mureco_dn_rel -= 0.001
      if hadtau_syst_llovr_dn_rel == 1: hadtau_syst_llovr_dn_rel -= 0.001

      if lostle_syst_mt_dn_rel == 1: lostle_syst_mt_dn_rel -= 0.001
      if lostle_syst_muiso_dn_rel == 1: lostle_syst_muiso_dn_rel -= 0.001
      if lostle_syst_mureco_dn_rel == 1: lostle_syst_mureco_dn_rel -= 0.001
      if lostle_syst_dilep_dn_rel == 1: lostle_syst_dilep_dn_rel -= 0.001

      outfile_perchn.write("hadtau_closure_chn%d   lnN - - - %.4f/%.4f - - - -\n" % (chn, 1-hadtau_syst_closure_dn_rel, 1+hadtau_syst_closure_up_rel))
      outfile_perchn.write("hadtau_taumu_chn%d   lnN - - - %.4f/%.4f - - - -\n" % (chn, 1-hadtau_syst_taumu_dn_rel, 1+hadtau_syst_taumu_up_rel))
      outfile_perchn.write("hadtau_isotrk   lnN - - - %.4f/%.4f - - - -\n" % (1-hadtau_syst_isotrk_dn_rel, 1+hadtau_syst_isotrk_up_rel))
      outfile_perchn.write("hadtau_pdf   lnN - - - %.4f/%.4f - - - -\n" % (1-hadtau_syst_pdf_dn_rel, 1+hadtau_syst_pdf_up_rel))
      outfile_perchn.write("hadtau_temp   lnN - - - %.4f/%.4f - - - -\n" % (1-hadtau_syst_temp_dn_rel, 1+hadtau_syst_temp_up_rel))
      outfile_perchn.write("hadtau_mistag   lnN - - - %.4f/%.4f - - - -\n" % (1-hadtau_syst_mistag_dn_rel, 1+hadtau_syst_mistag_up_rel))

      outfile_perchn.write("LLHadTau_mt   lnN - %.4f/%.4f - %.4f/%.4f - - - -\n" % (1-lostle_syst_mt_dn_rel, 1+lostle_syst_mt_up_rel, 1-hadtau_syst_Mt_dn_rel, 1+hadtau_syst_Mt_up_rel))
      outfile_perchn.write("LLHadTau_muiso   lnN - %.4f/%.4f - %.4f/%.4f - - - -\n" % (1-lostle_syst_muiso_dn_rel, 1+lostle_syst_muiso_up_rel, 1-hadtau_syst_muiso_dn_rel, 1+hadtau_syst_muiso_up_rel))
      outfile_perchn.write("LLHadTau_mureco   lnN - %.4f/%.4f - %.4f/%.4f - - - -\n" % (1-lostle_syst_mureco_dn_rel, 1+lostle_syst_mureco_up_rel, 1-hadtau_syst_mureco_dn_rel, 1+hadtau_syst_mureco_up_rel))
      outfile_perchn.write("LLHadTau_dilep   lnN - %.4f/%.4f - %.4f/%.4f - - - -\n" % (1-lostle_syst_dilep_dn_rel, 1+lostle_syst_dilep_up_rel, 1-hadtau_syst_llovr_dn_rel, 1+hadtau_syst_llovr_up_rel))

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
      
      if zinv_syst_shape_central_dn_rel == 1: zinv_syst_shape_central_dn_rel -= 0.001
      if zinv_syst_shape_stat_dn_rel == 1: zinv_syst_shape_stat_dn_rel -= 0.001
      if zinv_syst_jeu_dn_rel == 1: zinv_syst_jeu_dn_rel -= 0.001
      if zinv_syst_meu_dn_rel == 1: zinv_syst_meu_dn_rel -= 0.001
      if zinv_syst_scale_dn_rel == 1: zinv_syst_scale_dn_rel -= 0.001
      if zinv_syst_pdf_dn_rel == 1: zinv_syst_pdf_dn_rel -= 0.001
      if zinv_syst_trig_dn_rel == 1: zinv_syst_trig_dn_rel -= 0.001

      outfile_perchn.write("zinv_shape_central   lnN - - - - - %.4f/%.4f - -\n" % (1-zinv_syst_shape_central_dn_rel, 1+zinv_syst_shape_central_up_rel))
      outfile_perchn.write("zinv_shape_stat_chn%d   lnN - - - - - %.4f/%.4f - -\n" % (chn, 1-zinv_syst_shape_stat_dn_rel, 1+zinv_syst_shape_stat_up_rel))
      outfile_perchn.write("zinv_jeu   lnN - - - - - %.4f/%.4f - -\n" % (1-zinv_syst_jeu_dn_rel, 1+zinv_syst_jeu_up_rel))
      outfile_perchn.write("zinv_meu   lnN - - - - - %.4f/%.4f - -\n" % (1-zinv_syst_meu_dn_rel, 1+zinv_syst_meu_up_rel))
      outfile_perchn.write("zinv_scale   lnN - - - - - %.4f/%.4f - -\n" % (1-zinv_syst_scale_dn_rel, 1+zinv_syst_scale_up_rel))
      outfile_perchn.write("zinv_pdf   lnN - - - - - %.4f/%.4f - -\n" % (1-zinv_syst_pdf_dn_rel, 1+zinv_syst_pdf_up_rel))
      outfile_perchn.write("zinv_trig   lnN - - - - - %.4f/%.4f - -\n" % (1-zinv_syst_trig_dn_rel, 1+zinv_syst_trig_up_rel))

# qcd
      qcd_file.seek(0, 0)
      for qcd_line in qcd_file:
         qcd_splitline = procline(qcd_line)
         if qcd_splitline and qcd_splitline[0] == "QCD_TFactor_relative_err":
            qcd_tfactor_err_rel = float(qcd_splitline[chn])
         if qcd_splitline and qcd_splitline[0] == "QCD_NonClosure_relative_err":
            qcd_nonclosure_err_rel = float(qcd_splitline[chn])

      outfile_perchn.write("qcd_tfactor   lnN - - - - - - %.4f - \n" % (1+qcd_tfactor_err_rel)) 
      outfile_perchn.write("qcd_nonclosure_chn%d   lnN - - - - - - %.4f - \n" % (chn, 1+qcd_nonclosure_err_rel)) 
      outfile_perchn.write("ivtDphiCR_chn%d   lnU - - - - - - 10000 - \n" % (chn))

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

      if ttz_syst_pdf_dn_rel == 1: ttz_syst_pdf_dn_rel -= 0.001
      if ttz_syst_scale_dn_rel == 1: ttz_syst_scale_dn_rel -= 0.001
      if ttz_syst_rate_dn_rel == 1: ttz_syst_rate_dn_rel -= 0.001

      outfile_perchn.write("ttz_pdf   lnN - - - - - - - %.4f/%.4f\n" % (1-ttz_syst_pdf_dn_rel, 1+ttz_syst_pdf_up_rel)) 
      outfile_perchn.write("ttz_scale   lnN - - - - - - - %.4f/%.4f\n" % (1-ttz_syst_scale_dn_rel, 1+ttz_syst_scale_up_rel)) 
      outfile_perchn.write("ttz_rate   lnN - - - - - - - %.4f/%.4f\n" % (1-ttz_syst_rate_dn_rel, 1+ttz_syst_rate_up_rel)) 

      outfile_perchn.close()
      if len(options.outputdir) !=0:
         outfile_perchn_ori = open(options.outputdir + "/" + outbase_filename + "_ch" + str(chn) + ".txt", "r")
         outfile_perchn_modif = open(options.outputdir + "/" + outbase_filename + "_ch" + str(chn) + ".tmp.txt", "w")
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
      if len(options.outputdir) !=0:
         command_line = "mv "+options.outputdir + "/" + outbase_filename + "_ch" + str(chn) + ".tmp.txt "+options.outputdir + "/" + outbase_filename + "_ch" + str(chn) + ".txt"
      else:
         command_line = "mv "+outbase_filename + "_ch" + str(chn) + ".tmp.txt "+outbase_filename + "_ch" + str(chn) + ".txt"
      os.system(command_line)

# Additional output file for QCD
      if len(options.outputdir) !=0:
         qcd_outfile_perchn = open(options.outputdir + "/comb_invertDphi_"+ tmp_signal_key + "_ch" + str(chn) + ".txt", "w")
      else:
         qcd_outfile_perchn = open("comb_invertDphi_"+ tmp_signal_key + "_ch" + str(chn) + ".txt", "w")

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
      
      if qcd_contam_unc_dn_rel == 1: qcd_contam_unc_dn_rel -= 0.001

      qcd_outfile_perchn.write("ivtDphiCR_chn%d  lnU - 10000 -\n" % (chn))
      qcd_outfile_perchn.write("contamUnc_chn%d  lnN - - %.4f/%.4f\n" % (chn, 1-qcd_contam_unc_dn_rel, 1+qcd_contam_unc_up_rel))

      qcd_outfile_perchn.close() 

if __name__ == "__main__":
   main()
