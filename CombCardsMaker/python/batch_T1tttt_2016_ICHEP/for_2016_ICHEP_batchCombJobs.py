#!/usr/bin/python

import os
import time
import math
import ROOT
from optparse import OptionParser # Command line parsing

import for_2016_ICHEP_combCardPerChannel as combCard

def main():
   usage = "usage: %prog options"
   version = "%prog."
   parser = OptionParser(usage=usage,version=version)
   parser.add_option("-s", "--signal", action="store", dest="signaldir", type="string", default="", help="set signal card directory for batch processing")
   parser.add_option("-l", "--lostle", action="store", dest="lostle", type="string", default="lostle.txt", help="set lostle data card name")
   parser.add_option("-t", "--hadtau", action="store", dest="hadtau", type="string", default="hadtau.txt", help="set hadtau data card name")
   parser.add_option("-z", "--zinv", action="store", dest="zinv", type="string", default="zinv.txt", help="set zinv data card name")
   parser.add_option("-q", "--qcd", action="store", dest="qcd", type="string", default="qcd.txt", help="set qcd data card name")
   parser.add_option("-r", "--ttz", action="store", dest="ttz", type="string", default="ttz.txt", help="set ttz data card name")
   parser.add_option("-b", "--rare", action="store", dest="rare", type="string", default="rare.txt", help="set rare data card name")
   parser.add_option("-d", "--data", action="store", dest="data", type="string", default="data.txt", help="set data data card name")
   parser.add_option("-o", "--outputdir", action="store", dest="outputdir", type="string", default="", help="set combined card output directory")
   parser.add_option("-n", "--runlimit", action="store", dest="runlimit", type="string", default="yes", help="set run limit or not")
   parser.add_option("-m", "--model", action="store", dest="model", type="string", default="T2tt", help="set SMS model type")
   
   (options, args) = parser.parse_args()
   
   print 'lostle :', options.lostle
   print 'hadtau :', options.hadtau
   print 'zinv :', options.zinv
   print 'qcd :', options.qcd
   print 'ttz :', options.ttz
   print 'rare :', options.rare
   print 'data :', options.data
   print 'signaldir : ', options.signaldir
   print 'runlimit : ', options.runlimit

   if len(options.outputdir) !=0:
      if not os.path.exists(options.outputdir): os.mkdir(options.outputdir)
   else:
      print "output directory cannot be empty!"
      return 

   lostle_file = open(options.lostle)
   hadtau_file = open(options.hadtau)
   zinv_file = open(options.zinv)
   qcd_file = open(options.qcd)
   ttz_file = open(options.ttz)
   rare_file = open(options.rare)
   data_file = open(options.data)
   
   tarfilename = options.outputdir + ".tgz"

   submitLine = """universe = vanilla
Executable = """ + os.getcwd() + """/goScan.sh
Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = """ + os.getcwd() + """/""" + tarfilename + """
Output = logs/basicCheck_$(Process).stdout
Error = logs/basicCheck_$(Process).stderr
Log = logs/basicCheck_$(Process).log
notify_user = ${LOGNAME}@FNAL.GOV

"""
   fileParts = [submitLine]

   signal_key_list = []

# Batch making the cards for individual channels and a combined one for the combination tool   
   for signal_name in os.listdir(options.signaldir):
      if not ("signal" in signal_name): continue
      if not (".txt" in signal_name): continue
      if signal_name.startswith('.'): continue

      splitsignalinput = signal_name.split("/")
      stripDirSignalInput = splitsignalinput[-1]

      if "signal_" in stripDirSignalInput: tmp_signal_key = stripDirSignalInput.replace("signal_", "")
      elif "signal" in stripDirSignalInput: tmp_signal_key = stripDirSignalInput.replace("signal", "")
      tmp_signal_key = tmp_signal_key.replace(".txt", "")

      signal_key_list.append(tmp_signal_key)

      tmp_signal_name = signal_name
      full_signal_name = os.path.join(options.signaldir, signal_name)

      outputdir = options.outputdir + "/" + tmp_signal_name.replace(".txt", "")
   
      if options.runlimit == "yes":
         if len(outputdir) !=0:
            if not os.path.exists(outputdir): os.mkdir(outputdir)

         print "processing the signal file : %s in output dir : %s" % (signal_name, outputdir)

         signal_file = open(full_signal_name)
# Core function to produce the cards
         combCard.prodCardPerChn(tmp_signal_key, outputdir, lostle_file, hadtau_file, zinv_file, qcd_file, ttz_file, rare_file, data_file, signal_file)
         signal_file.close()

# Making a combined card
         rm_comb_cards_command = "rm "
         make_allComb_command = "combineCards.py "
         for card_name in os.listdir(outputdir):
            if not ("comb" in card_name): continue
            make_allComb_command += (outputdir + "/" + card_name + " ")
            rm_comb_cards_command += (outputdir + "/" + card_name + " ")
         make_allComb_command += (" > " + outputdir + "/allComb_"+tmp_signal_name)
         os.system(make_allComb_command)
         os.system(rm_comb_cards_command)

      allComb_file_dir = outputdir + "/" + "allComb_" + tmp_signal_name
      allComb_run_output_filename = "log_allComb_" + tmp_signal_name + ".lg"
      fileParts.append("Arguments = %s $ENV(CMSSW_BASE) %s\nQueue\n\n"%(allComb_file_dir, allComb_run_output_filename))

   lostle_file.close()
   hadtau_file.close()
   zinv_file.close()
   qcd_file.close()
   ttz_file.close()
   rare_file.close()
   data_file.close()

# Preparing tarred cards, executable and condor configuration files for limit jobs
   if options.runlimit == "yes":
      os.system("tar -czvf "+tarfilename+" "+options.outputdir + "/*/allComb_*.txt")
      os.system("rm -fr "+options.outputdir)

   exeLine = """#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

cd $2/src
eval `scramv1 runtime -sh`

cd ${_CONDOR_SCRATCH_DIR}

tar -zxvf """ + tarfilename + """

combine -M Asymptotic $1 > $3

"""
   exeFile = open("goScan.sh", "w")
   exeFile.write(exeLine)
   exeFile.close()

   condorFile = open("condor_submit.txt", "w")
   condorFile.write(''.join(fileParts))
   condorFile.close()

# Submit the jobs
   if not os.path.exists("logs"): os.mkdir("logs")
   if options.runlimit == "yes":
      os.system("echo 'condor_submit condor_submit.txt'")
      os.system('condor_submit condor_submit.txt')

# Checking if the jobs are done
   while True:
      condor_check_command = "condor_q -submitter " + os.environ["USER"] + " -wide > tmp_condor_check_output.txt" 
      os.system(condor_check_command)
      tmp_condor_file = open("tmp_condor_check_output.txt")
      cnt_runjobs =0
      for condor_line in tmp_condor_file:
         if "goScan.sh" in condor_line and "allComb" in condor_line : cnt_runjobs += 1
      print "current jobs in queue: %d" % (cnt_runjobs)
      if cnt_runjobs ==0: break;
      time.sleep(60)

# Removing intermediate files
   os.system("rm tmp_condor_check_output.txt")
   os.system("rm higgsCombineTest.Asymptotic.mH120.root")
   os.system("rm roostats*.root")

# Parsing the output files to make new files with correct format for drawing steps
   cards_for_plotting_dir_name = "cards_for_plotting"
   if not os.path.exists(cards_for_plotting_dir_name): os.mkdir(cards_for_plotting_dir_name)

   print "\nsignal_key_list : ", signal_key_list
   filelist_for_plotting = open(cards_for_plotting_dir_name+"/filelist.txt", "w")

   xSec_root_file = ROOT.TFile("xSec.root")
   if "T2tt" in options.model or "T2bb" in options.model or "T2tb" in options.model or "T6tt" in options.model:
      h1_xSec = xSec_root_file.Get("stop_xsection")
   elif "T1tt" in options.model or "T5tt" in options.model:
      h1_xSec = xSec_root_file.Get("gluino_xsection")

   for log_file_name in os.listdir(os.getcwd()):
      if not ("log_allComb" in log_file_name) and not (".lg" in log_file_name) : continue
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
         if "Observed Limit:" in log_line: obsLimit = float(log_line_split[-1])*xSec
         if "Expected 16.0%:" in log_line: expLimit_m1sigma = float(log_line_split[-1])*xSec
         if "Expected 50.0%:" in log_line: expLimit = float(log_line_split[-1])*xSec
         if "Expected 84.0%:" in log_line: expLimit_p1sigma = float(log_line_split[-1])*xSec

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

if __name__ == "__main__":
   main()
