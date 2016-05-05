#!/usr/bin/python
import sys
import os
import time
import math
import ROOT as rt
from optparse import OptionParser # Command line parsing

def main():
   usage = "usage: %prog options"
   version = "%prog."
   parser = OptionParser(usage=usage,version=version)
   parser.add_option("-f", "--file", action="store", dest="rtfile", type="string", default="results_T2tt/CLs_SMS_BR100pct.root", help="set root file name")
   parser.add_option("-e", "--extra", action="store", dest="extrafile", type="string", default="flatNtp_v6/cards_for_plotting/signal_strength.root", help="set root file name")
   parser.add_option("-o", "--outputdir", action="store", dest="outputdir", type="string", default="", help="set combined card output directory")
   parser.add_option("-m", "--model", action="store", dest="model", type="string", default="T2tt", help="set SMS model type")
   parser.add_option("-u", "--hist", action="store", dest="hist", type="string", default="strength", help="produce 2D hist for obs and exp strength. xSec limits are always produced for exclusion curves")
   parser.add_option("-t", "--type", action="store", dest="opttype", type="string", default="all", help="options are: all, limit2D, strength2D, curves")
   parser.add_option("-c", "--cover", action="store", dest="cover", type="string", default="yes", help="yes to cover the top corridor")
   
   (options, args) = parser.parse_args()
   
   print 'rtfile : ', options.rtfile
   print 'extrafile : ', options.extrafile
   print 'outputdir : ', options.outputdir
   print 'opttype : ', options.opttype
   print 'cover : ', options.cover

   splitrtinput = options.rtfile.split("/")
   stripDirrtInput = splitrtinput[-1]

   sig_file = rt.TFile.Open(options.rtfile)
   extra_file = rt.TFile.Open(options.extrafile)

   if options.outputdir:
      if "T2tt" in options.rtfile:
         output_filename = options.outputdir + "/T2tt_" + stripDirrtInput
      elif "T2tb" in options.rtfile:
         output_filename = options.outputdir + "/T2tb_" + stripDirrtInput
   else:
      if "T2tt" in options.rtfile:
         output_filename = "T2tt_" + stripDirrtInput
      elif "T2tb" in options.rtfile:
         output_filename = "T2tb_" + stripDirrtInput

   print 'output_filename : ', output_filename

   if output_filename == options.rtfile:
      print "\noutput_filename : %s  SAME AS  input rtfile : %s\n" %(output_filename, options.rtfile)
      sys.exit()

   output_rtfile = rt.TFile(output_filename, "RECREATE")

   if "T2tt" in options.rtfile:
      output_rtfile.cd()

      if options.opttype == "all" or options.opttype == "limit2D":
         combined_obsLimit_BR100pct = sig_file.Get("combined_obsLimit_BR100pct")
         if options.cover == "yes":
            nXBins = combined_obsLimit_BR100pct.GetXaxis().GetNbins()
            nYBins = combined_obsLimit_BR100pct.GetYaxis().GetNbins()
            for iy in xrange(nYBins):
               yCen = combined_obsLimit_BR100pct.GetYaxis().GetBinCenter(iy+1)
               for ix in xrange(nXBins):
                  xCen = combined_obsLimit_BR100pct.GetXaxis().GetBinCenter(ix+1)
                  if xCen <= yCen+175+25 and xCen >= yCen+175-25:
                     combined_obsLimit_BR100pct.SetBinContent(ix+1, iy+1, 0)
                     combined_obsLimit_BR100pct.SetBinError(ix+1, iy+1, 0)
         combined_obsLimit_BR100pct.Write()

         combined_expLimit_BR100pct = sig_file.Get("combined_expLimit_BR100pct")
         if options.cover == "yes":
            nXBins = combined_expLimit_BR100pct.GetXaxis().GetNbins()
            nYBins = combined_expLimit_BR100pct.GetYaxis().GetNbins()
            for iy in xrange(nYBins):
               yCen = combined_expLimit_BR100pct.GetYaxis().GetBinCenter(iy+1)
               for ix in xrange(nXBins):
                  xCen = combined_expLimit_BR100pct.GetXaxis().GetBinCenter(ix+1)
                  if xCen <= yCen+175+25 and xCen >= yCen+175-25:
                     combined_expLimit_BR100pct.SetBinContent(ix+1, iy+1, 0)
                     combined_expLimit_BR100pct.SetBinError(ix+1, iy+1, 0)
         combined_expLimit_BR100pct.Write()

      if options.opttype == "all" or options.opttype == "curves":
         combined_obsExclOneTimesProspino_base_BR100pct = sig_file.Get("combined_obsExclOneTimesProspino_base_BR100pct")
         combined_obsExclOneTimesProspino_base_BR100pct.Write()
         combined_obsExclOneTimesProspino_extra_BR100pct = sig_file.Get("combined_obsExclOneTimesProspino_extra_BR100pct")
         combined_obsExclOneTimesProspino_extra_BR100pct.Write()
   
         combined_expExclOneTimesProspino_base_BR100pct = sig_file.Get("combined_expExclOneTimesProspino_base_BR100pct")
         combined_expExclOneTimesProspino_base_BR100pct.Write()
         combined_expExclOneTimesProspino_extra_BR100pct = sig_file.Get("combined_expExclOneTimesProspino_extra_BR100pct")
         combined_expExclOneTimesProspino_extra_BR100pct.Write()
   
         combined_expExclPlusOneSigmaProspino_base_BR100pct = sig_file.Get("combined_expExclPlusOneSigmaProspino_base_BR100pct")
         combined_expExclPlusOneSigmaProspino_base_BR100pct.Write()
         combined_expExclPlusOneSigmaProspino_extra_BR100pct = sig_file.Get("combined_expExclPlusOneSigmaProspino_extra_BR100pct")
         combined_expExclPlusOneSigmaProspino_extra_BR100pct.Write()
   
         combined_expExclMinusOneSigmaProspino_base_BR100pct = sig_file.Get("combined_expExclMinusOneSigmaProspino_base_BR100pct")
         combined_expExclMinusOneSigmaProspino_base_BR100pct.Write()
   
         combined_obsExclPlusSysErrProspino_base_BR100pct = sig_file.Get("combined_obsExclPlusSysErrProspino_base_BR100pct")
         combined_obsExclPlusSysErrProspino_base_BR100pct.Write()
         combined_obsExclPlusSysErrProspino_extra_BR100pct = sig_file.Get("combined_obsExclPlusSysErrProspino_extra_BR100pct")
         combined_obsExclPlusSysErrProspino_extra_BR100pct.Write()
   
         combined_obsExclMinusSysErrProspino_base_BR100pct = sig_file.Get("combined_obsExclMinusSysErrProspino_base_BR100pct")
         combined_obsExclMinusSysErrProspino_base_BR100pct.Write()
         combined_obsExclMinusSysErrProspino_extra_BR100pct = sig_file.Get("combined_obsExclMinusSysErrProspino_extra_BR100pct")
         combined_obsExclMinusSysErrProspino_extra_BR100pct.Write()

         diagonalCoverBand = sig_file.Get("diagonalCoverBand")
         diagonalCoverBand.Write()
  
      if options.opttype == "all" or options.opttype == "strength2D": 
         exp_strength = extra_file.Get("exp_strength")
         exp_strength.Write()
         obs_strength = extra_file.Get("obs_strength")
         obs_strength.Write()

   elif "T2tb" in options.rtfile:
      output_rtfile.cd()

      if options.opttype == "all" or options.opttype == "limit2D":
         combined_obsLimit_BR50pct = sig_file.Get("combined_obsLimit_BR50pct")
         if options.cover == "yes":
            nXBins = combined_obsLimit_BR50pct.GetXaxis().GetNbins()
            nYBins = combined_obsLimit_BR50pct.GetYaxis().GetNbins()
            for iy in xrange(nYBins):
               yCen = combined_obsLimit_BR50pct.GetYaxis().GetBinCenter(iy+1)
               for ix in xrange(nXBins):
                  xCen = combined_obsLimit_BR50pct.GetXaxis().GetBinCenter(ix+1)
                  if xCen <= yCen+175+25:
                     combined_obsLimit_BR50pct.SetBinContent(ix+1, iy+1, 0)
                     combined_obsLimit_BR50pct.SetBinError(ix+1, iy+1, 0)
         combined_obsLimit_BR50pct.Write()

         combined_expLimit_BR50pct = sig_file.Get("combined_expLimit_BR50pct")
         if options.cover == "yes":
            nXBins = combined_expLimit_BR50pct.GetXaxis().GetNbins()
            nYBins = combined_expLimit_BR50pct.GetYaxis().GetNbins()
            for iy in xrange(nYBins):
               yCen = combined_expLimit_BR50pct.GetYaxis().GetBinCenter(iy+1)
               for ix in xrange(nXBins):
                  xCen = combined_expLimit_BR50pct.GetXaxis().GetBinCenter(ix+1)
                  if xCen <= yCen+175+25:
                     combined_expLimit_BR50pct.SetBinContent(ix+1, iy+1, 0)
                     combined_expLimit_BR50pct.SetBinError(ix+1, iy+1, 0)
         combined_expLimit_BR50pct.Write()

      if options.opttype == "all" or options.opttype == "curves":
         combined_obsExclOneTimesProspino_base_BR50pct = sig_file.Get("combined_obsExclOneTimesProspino_base_BR50pct")
         combined_obsExclOneTimesProspino_base_BR50pct.Write()

         combined_expExclOneTimesProspino_base_BR50pct = sig_file.Get("combined_expExclOneTimesProspino_base_BR50pct")
         combined_expExclOneTimesProspino_base_BR50pct.Write()

         combined_expExclPlusOneSigmaProspino_base_BR50pct = sig_file.Get("combined_expExclPlusOneSigmaProspino_base_BR50pct")
         combined_expExclPlusOneSigmaProspino_base_BR50pct.Write()

         combined_expExclMinusOneSigmaProspino_base_BR50pct = sig_file.Get("combined_expExclMinusOneSigmaProspino_base_BR50pct")
         combined_expExclMinusOneSigmaProspino_base_BR50pct.Write()

         combined_obsExclPlusSysErrProspino_base_BR50pct = sig_file.Get("combined_obsExclPlusSysErrProspino_base_BR50pct")
         combined_obsExclPlusSysErrProspino_base_BR50pct.Write()

         combined_obsExclMinusSysErrProspino_base_BR50pct = sig_file.Get("combined_obsExclMinusSysErrProspino_base_BR50pct")
         combined_obsExclMinusSysErrProspino_base_BR50pct.Write()

         diagonalCoverBand = sig_file.Get("diagonalCoverBand")
         diagonalCoverBand.Write()
  
      if options.opttype == "all" or options.opttype == "strength2D":
         exp_strength = extra_file.Get("exp_strength")
         exp_strength.Write()
         obs_strength = extra_file.Get("obs_strength")
         obs_strength.Write()

   else:
      print "\nNot supported signals!\n"

   output_rtfile.Write()
   output_rtfile.Close()

if __name__ == "__main__":
   main()
