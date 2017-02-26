#!/usr/bin/python

import os
import time
import math
import ROOT as rt
from optparse import OptionParser # Command line parsing

def main():
   usage = "usage: %prog options"
   version = "%prog."
   parser = OptionParser(usage=usage,version=version)
   parser.add_option("-s", "--signal", action="store", dest="signaldir", type="string", default="", help="set signal card directory for batch processing")
   parser.add_option("-o", "--outputdir", action="store", dest="outputdir", type="string", default="", help="set combined card output directory")
   parser.add_option("-n", "--runlimit", action="store", dest="runlimit", type="string", default="yes", help="set run limit or not")
   parser.add_option("-m", "--model", action="store", dest="model", type="string", default="T2tt", help="set SMS model type")
   parser.add_option("-u", "--hist", action="store", dest="hist", type="string", default="strength", help="produce 2D hist for obs and exp strength. xSec limits are always produced for exclusion curves")
   
   (options, args) = parser.parse_args()
   
   print 'signaldir : ', options.signaldir
   print 'runlimit : ', options.runlimit


   nXbins = 35
   nYbins = 22
   loX =  87.5
   hiX = 962.5
   loY = -12.5
   hiY = 537.5

   if "T1tt" in options.model or "T5tt" in options.model:
      nXbins = 52
      nYbins = 49
      loX = 612.5
      hiX = 1912.5
      loY = -12.5
      hiY = 1212.5

   cards_for_plotting_dir_name = "cards_for_plotting"
   if not os.path.exists(cards_for_plotting_dir_name): os.mkdir(cards_for_plotting_dir_name)

   if options.hist == "strength":
      signal_strength_file = rt.TFile(cards_for_plotting_dir_name + "/signal_strength.root", "RECREATE")
      h2_exp_strength = rt.TH2D("exp_strength", "exp_strength", nXbins, loX, hiX, nYbins, loY, hiY)
      h2_obs_strength = rt.TH2D("obs_strength", "obs_strength", nXbins, loX, hiX, nYbins, loY, hiY)
   if options.hist == "xsec":
      signal_xsec_UL_file = rt.TFile(cards_for_plotting_dir_name + "/signal_xsec_UL.root", "RECREATE")
      h2_exp_xsec_UL = rt.TH2D("exp_xsec_UL", "exp_xsec_UL", nXbins, loX, hiX, nYbins, loY, hiY)
      h2_expPlusOneSigma_xsec_UL = rt.TH2D("expPlusOneSigma_xsec_UL", "expPlusOneSigma_xsec_UL", nXbins, loX, hiX, nYbins, loY, hiY)
      h2_expMinusOneSigma_xsec_UL = rt.TH2D("expMinusOneSigma_xsec_UL", "expMinusOneSigma_xsec_UL", nXbins, loX, hiX, nYbins, loY, hiY)
      h2_obs_xsec_UL = rt.TH2D("obs_xsec_UL", "obs_xsec_UL", nXbins, loX, hiX, nYbins, loY, hiY)
#      h2_obsPlusSysErr_xsec_UL = rt.TH2D("obsPlusSysErr_xsec_UL", "obsPlusSysErr_xsec_UL", nXbins, loX, hiX, nYbins, loY, hiY)
#      h2_obsMinusSysErr_xsec_UL = rt.TH2D("obsMinusSysErr_xsec_UL", "obsMinusSysErr_xsec_UL", nXbins, loX, hiX, nYbins, loY, hiY)

   signal_key_list = []

# Batch making the cards for individual channels and a combined one for the combination tool   
   for signal_name in os.listdir(os.getcwd()):
      if not ("log_allComb_" in signal_name): continue
      if not (".txt" in signal_name): continue
      if signal_name.startswith('.'): continue

      splitsignalinput = signal_name.split("/")
      stripDirSignalInput = splitsignalinput[-1]

      if "log_allComb_signal_" in stripDirSignalInput: tmp_signal_key = stripDirSignalInput.replace("log_allComb_signal_", "")
      elif "log_allComb_" in stripDirSignalInput: tmp_signal_key = stripDirSignalInput.replace("log_allComb_", "")
      elif "signal" in stripDirSignalInput: tmp_signal_key = stripDirSignalInput.replace("signal", "")
      tmp_signal_key = tmp_signal_key.replace(".txt.lg", "")

      signal_key_list.append(tmp_signal_key)

      tmp_signal_name = signal_name
      full_signal_name = os.path.join(options.signaldir, signal_name)

      outputdir = options.outputdir + "/" + tmp_signal_name.replace(".txt", "")

# Parsing the output files to make new files with correct format for drawing steps
   print "\nsignal_key_list : ", signal_key_list
   filelist_for_plotting = open(cards_for_plotting_dir_name+"/filelist.txt", "w")

   xSec_root_file = rt.TFile("xSec.root")
   if "T2tt" in options.model or "T2bb" in options.model or "T2tb" in options.model or "T6tt" in options.model:
      h1_xSec = xSec_root_file.Get("stop_xsection")
   elif "T1tt" in options.model or "T5tt" in options.model:
      h1_xSec = xSec_root_file.Get("gluino_xsection")

   for log_file_name in os.listdir(os.getcwd()):
      if not ("log_allComb" in log_file_name) and not (".lg" in log_file_name) : continue
      if log_file_name.startswith('.'): continue
      picked_key_name = ''
      for key_name in signal_key_list:
         more_key_name = "signal_"+key_name+".txt"
         if more_key_name in log_file_name: 
            picked_key_name = key_name
            break

      if len(picked_key_name) ==0: continue

      split_key = picked_key_name.split("_")
      momMass = split_key[0]
      dauMass = split_key[1]

      binIdx = h1_xSec.FindBin(float(momMass))
      xSec = h1_xSec.GetBinContent(binIdx)
     
      log_file = open(log_file_name)
      for log_line in log_file:
         log_line_split = log_line.split()
         if "Observed Limit:" in log_line: 
            obsLimit = float(log_line_split[-1])*xSec
            if options.hist == "strength": h2_obs_strength.Fill(float(momMass), float(dauMass), float(log_line_split[-1]))
            if options.hist == "xsec": h2_obs_xsec_UL.Fill(float(momMass), float(dauMass), float(log_line_split[-1])*xSec)
         if "Expected 16.0%:" in log_line:
            expLimit_m1sigma = float(log_line_split[-1])*xSec
            if options.hist == "xsec": h2_expMinusOneSigma_xsec_UL.Fill(float(momMass), float(dauMass), float(log_line_split[-1])*xSec)
         if "Expected 50.0%:" in log_line: 
            expLimit = float(log_line_split[-1])*xSec
            if options.hist == "strength": h2_exp_strength.Fill(float(momMass), float(dauMass), float(log_line_split[-1]))
            if options.hist == "xsec": h2_exp_xsec_UL.Fill(float(momMass), float(dauMass), float(log_line_split[-1])*xSec)
         if "Expected 84.0%:" in log_line: 
            expLimit_p1sigma = float(log_line_split[-1])*xSec
            if options.hist == "xsec": h2_expPlusOneSigma_xsec_UL.Fill(float(momMass), float(dauMass), float(log_line_split[-1])*xSec)

      format_file_name = "mSUGRA_" + picked_key_name + "_10_0_1.dat"
      format_file = open(cards_for_plotting_dir_name+"/"+format_file_name, "w")

      formatLine = """Azero = 0
Mu = 1
Mzero = """ + momMass + """
Mhalf = """ + dauMass + """
Xsection = """ + str(xSec) + """
limit.cls.observed = """ + str(obsLimit) + """
limit.cls.observed.error = 0
limit.cls.expected.m1sigma = """ + str(expLimit_m1sigma) + """
limit.cls.expected = """ + str(expLimit) + """
limit.cls.expected.p1sigma = """ + str(expLimit_p1sigma) + """
limit.cls.expected.1sigmaCoverage = -1
limit.cls.expected.2sigmaCoverage = -1
limit.cls.expected.m2sigma = -1
limit.cls.expected.p2sigma = -1

"""

      format_file.write(formatLine)
      format_file.close()

      filelist_for_plotting.write(format_file_name + "\n")

   filelist_for_plotting.close()      

   if options.hist == "strength":
#      signal_strength_file.cd()
#      h2_obs_strength.Write()
#      h2_exp_strength.Write()
      signal_strength_file.Write()
      signal_strength_file.Close()
   if options.hist == "xsec":
      signal_xsec_UL_file.Write()
      signal_xsec_UL_file.Close()

if __name__ == "__main__":
   main()
