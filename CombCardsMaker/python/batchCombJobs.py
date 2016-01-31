#!/usr/bin/python

#from os import listdir
#from os import sys
#from os import system
#from os import path
import os
import math
from optparse import OptionParser # Command line parsing

import combCardPerChannel as combCard

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
   parser.add_option("-d", "--data", action="store", dest="data", type="string", default="data.txt", help="set data data card name")
   parser.add_option("-o", "--outputdir", action="store", dest="outputdir", type="string", default="", help="set combined card output directory")
   
   (options, args) = parser.parse_args()
   
   print 'lostle :', options.lostle
   print 'hadtau :', options.hadtau
   print 'zinv :', options.zinv
   print 'qcd :', options.qcd
   print 'ttz :', options.ttz
   print 'data :', options.data
   print 'signaldir : ', options.signaldir

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
   data_file = open(options.data)
   
   for signal_name in os.listdir(options.signaldir):
      if not ("signal" in signal_name): continue

      splitsignalinput = signal_name.split("/")
      stripDirSignalInput = splitsignalinput[-1]

      if "signal_" in stripDirSignalInput: tmp_signal_key = stripDirSignalInput.replace("signal_", "")
      elif "signal" in stripDirSignalInput: tmp_signal_key = stripDirSignalInput.replace("signal", "")
      tmp_signal_key = tmp_signal_key.replace(".txt", "")

      tmp_signal_name = signal_name
      full_signal_name = os.path.join(options.signaldir, signal_name)

      outputdir = options.outputdir + "/" + tmp_signal_name.replace(".txt", "")
   
      if len(outputdir) !=0:
         if not os.path.exists(outputdir): os.mkdir(outputdir)

      print "processing the signal file : %s in output dir : %s" % (signal_name, outputdir)

      signal_file = open(full_signal_name)

      combCard.prodCardPerChn(tmp_signal_key, outputdir, lostle_file, hadtau_file, zinv_file, qcd_file, ttz_file, data_file, signal_file)

      signal_file.close()

   lostle_file.close()
   hadtau_file.close()
   zinv_file.close()
   qcd_file.close()
   ttz_file.close()
   data_file.close()

   tarfilename = options.outputdir + ".tgz" 
   os.system("tar -czvf "+tarfilename+" "+options.outputdir)

   exeFile = """#!/bin/bash

   export PATH=${PATH}:/cvmfs/cms.cern.ch/common
   export CMS_PATH=/cvmfs/cms.cern.ch

   cd $2/src
   eval `scramv1 runtime -sh`

   cd ${_CONDOR_SCRATCH_DIR}

   tar -zxvf batch.tgz

   combine -M Asymptotic $1 > $3

"""

   submitFile = """universe = vanilla
   Executable = $ENV(CMSSW_BASE)/src/SusyAnaTools/Tools/condor/goMakePlots.sh
   Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
   Should_Transfer_Files = YES
   WhenToTransferOutput = ON_EXIT
   Transfer_Input_Files = $ENV(CMSSW_BASE)/src/SusyAnaTools/Tools/basicCheck, $ENV(CMSSW_BASE)/src/SusyAnaTools/Tools/condor/goMakePlots.sh, $ENV(CMSSW_BASE)/lib/$ENV(SCRAM_ARCH)/librecipeAUXOxbridgeMT2.so
   Output = logs/basicCheck_$(Process).stdout
   Error = logs/basicCheck_$(Process).stderr
   Log = logs/basicCheck_$(Process).log
   notify_user = ${LOGNAME}@FNAL.GOV
   x509userproxy = $ENV(X509_USER_PROXY)

"""
   

if __name__ == "__main__":
   main()
