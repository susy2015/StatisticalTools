#!/usr/bin/python

import os
import time
import math
import ROOT as rt
from optparse import OptionParser # Command line parsing

nTotBins = 37

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
   parser.add_option("-n", "--runlimit", action="store", dest="runlimit", type="string", default="yes", help="set run limit or not")
   parser.add_option("-m", "--model", action="store", dest="model", type="string", default="T2tt", help="set SMS model type")
   parser.add_option("-u", "--hist", action="store", dest="hist", type="string", default="strength", help="produce 2D hist for obs and exp strength. xSec limits are always produced for exclusion curves")
   
   (options, args) = parser.parse_args()
   
   print 'lostle :', options.lostle
   print 'hadtau :', options.hadtau
   print 'zinv :', options.zinv
   print 'qcd :', options.qcd
   print 'ttz :', options.ttz
   print 'data :', options.data
   print 'signaldir : ', options.signaldir
   print 'runlimit : ', options.runlimit


   nXbins = 32
   nYbins = 22
   loX = 112.5
   hiX = 912.5
   loY = -12.5
   hiY = 537.5

   sig_file = rt.TFile("signalScan_SMS-T2tt_37Bins.root")

   output_file = rt.TFile("output.root", "RECREATE")
   h2_acc_baseline = rt.TH2D("acc_baseline", "acceptance after baseline cuts", nXbins, loX, hiX, nYbins, loY, hiY)

   signal_key_list = []

# Batch making the cards for individual channels and a combined one for the combination tool   
   for signal_name in os.listdir(options.signaldir):
      if not ("signal" in signal_name): continue
      if not (".txt" in signal_name): continue
      if signal_name.startswith('.'): continue

      splitsignalinput = signal_name.split("/")
      stripDirSignalInput = splitsignalinput[-1]

      if "signal_" in stripDirSignalInput: signal_key = stripDirSignalInput.replace("signal_", "")
      elif "signal" in stripDirSignalInput: signal_key = stripDirSignalInput.replace("signal", "")
      signal_key = signal_key.replace(".txt", "")

      splitsMass = signal_key.split("_")
      momMass = int(splitsMass[0])
      dauMass = int(splitsMass[1])

      entry_str = "totEntries_"+signal_key
      nSearchBin_str = "baseline_nSearchBin_"+signal_key

      h1_entry = sig_file.Get(entry_str)
      entries = float(h1_entry.GetBinContent(1))

      h1_nSearchBin = sig_file.Get(nSearchBin_str)
      nSearchBin = h1_nSearchBin.GetXaxis().GetNbins()

      print "\nsignal_key : %s   entry_str : %s   nSearchBin_str : %s   entries : %f\n" %(signal_key, entry_str, nSearchBin_str, entries)

      sumEvt_baseline = 0

#      for ib in range(0, nSearchBin):
      for ib in range(0, nTotBins):
         perBin_cont = h1_nSearchBin.GetBinContent(ib+1)

         sumEvt_baseline += perBin_cont

         perBin_err = h1_nSearchBin.GetBinError(ib+1)
         perBin_acc = float(perBin_cont/entries)

         if perBin_cont !=0:
            perBin_acc_err = perBin_acc*math.sqrt(1/perBin_cont - 1/entries)
         else:
            perBin_acc_err = 0

         print "%2dth bin  %9.4f +- %7.4f  ->  acc : %10.8f +- %10.8f\n" %(ib, perBin_cont, perBin_err, perBin_acc, perBin_acc_err)

      signal_key_list.append(signal_key)

      acc_baseline = sumEvt_baseline/entries
      acc_err_baseline = acc_baseline*math.sqrt(1/sumEvt_baseline - 1/entries)

      binXidx = h2_acc_baseline.GetXaxis().FindBin(float(momMass))
      binYidx = h2_acc_baseline.GetYaxis().FindBin(float(dauMass))

      h2_acc_baseline.SetBinContent(binXidx, binYidx, acc_baseline)
      h2_acc_baseline.SetBinError(binXidx, binYidx, acc_err_baseline)

   output_file.Write()
   output_file.Close()

if __name__ == "__main__":
   main()
