#!/usr/bin/python

import os
import time
import math
import ROOT
from optparse import OptionParser # Command line parsing

import combCardPerChannel as combCard

def main():
   usage = "usage: %prog options"
   version = "%prog."
   parser = OptionParser(usage=usage,version=version)
   parser.add_option("-s", "--signal", action="store", dest="signal", type="string", default="", help="set signal card for batch processing")
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
   print 'signal : ', options.signal
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

#   for data_line in data_file:
#      splitline = combCard.procline(data_line)
#      if splitline:
#         if splitline[0] == "luminosity" : lumi = splitline[1]
#         if splitline[0] == "channels" : channels = splitline[1]
       
   if not ("signal" in options.signal): return
   if not (".txt" in options.signal): return 
   if options.signal.startswith('.'): return

   splitsignalinput = options.signal.split("/")
   stripDirSignalInput = splitsignalinput[-1]

   if "signal_" in stripDirSignalInput: tmp_signal_key = stripDirSignalInput.replace("signal_", "")
   elif "signal" in stripDirSignalInput: tmp_signal_key = stripDirSignalInput.replace("signal", "")
   tmp_signal_key = tmp_signal_key.replace(".txt", "")

   tmp_signal_name = options.signal
   full_signal_name = options.signal

   outputdir = options.outputdir + "/" + tmp_signal_key
   
   if len(outputdir) !=0:
      if not os.path.exists(outputdir): os.mkdir(outputdir)

   print "processing the signal file : %s in output dir : %s" % (options.signal, outputdir)

   signal_file = open(full_signal_name)
# Core function to produce the cards
   combCard.prodCardPerChn(tmp_signal_key, outputdir, lostle_file, hadtau_file, zinv_file, qcd_file, ttz_file, rare_file, data_file, signal_file)
   signal_file.close()

# Making a combined card
   make_allComb_command = "combineCards.py "
   for card_name in os.listdir(outputdir):
      if not ("comb" in card_name): continue
      if card_name.startswith('.'): continue
      make_allComb_command += (outputdir + "/" + card_name + " ")
   make_allComb_command += (" > " + outputdir + "/allComb_"+tmp_signal_key+".txt")
   print make_allComb_command
   os.system(make_allComb_command)

   bestChnDict = {}

   for ch in range(1, int(combCard.channels)+1):
      expLimit = 999.0
      make_perCh_comb_command = "combineCards.py " + outputdir+"/comb_*ch"+str(ch)+".txt > " + outputdir+"/allComb_ch"+str(ch)+".txt"
      print make_perCh_comb_command
      os.system(make_perCh_comb_command)

      limit_out_file = outputdir+"/log_allComb_ch"+str(ch)+".lg"
      if options.runlimit == "yes":
         runlimit_perCh_command = "combine -M Asymptotic " + outputdir+"/allComb_ch"+str(ch)+".txt > " + limit_out_file
         print runlimit_perCh_command
         os.system(runlimit_perCh_command)

      if os.path.isfile(limit_out_file):
         log_file = open(limit_out_file)
         for log_line in log_file:
            log_line_split = log_line.split()
            if "Expected 50.0%:" in log_line: expLimit = float(log_line_split[-1])
         log_file.close()
         bestChnDict[ch] = expLimit

   output_res_file_name = outputdir + "/output_results_"+tmp_signal_key+".txt"
   output_res_latex_file_name = outputdir + "/output_results_"+tmp_signal_key+".tex"
   output_res_file = open(output_res_file_name, "w")
   output_res_latex_file = open(output_res_latex_file_name, "w")
   latex_begin_line = """\\begin{table}[htb]
  \\begin{center}
    \\begin{tabular}{|c|c||c||c|c|c||c|c|c|}
      \\hline
      Bin & limit (r) & Obs. & Pred. & Stat. & Syst. & Signal & Stat. & Syst. \\\\
      \\hline
"""
   output_res_latex_file.write(latex_begin_line)
   for (key, value) in sorted(bestChnDict.items(), key=lambda (k,v): (v,k)):
      data_rate = combCard.glb_data_rate[key-1]
      tot_rate = combCard.glb_tot_rate[key-1]
      tot_stat_up = combCard.glb_tot_stat_up[key-1]
      tot_stat_dn = combCard.glb_tot_stat_dn[key-1]
      tot_syst_up = combCard.glb_tot_syst_up[key-1]
      tot_syst_dn = combCard.glb_tot_syst_dn[key-1]

      signal_rate = combCard.glb_signal_rate[key-1]
      signal_stat_up = combCard.glb_signal_stat_up[key-1]
      signal_stat_dn = combCard.glb_signal_stat_dn[key-1]
      signal_syst_up = combCard.glb_signal_syst_up[key-1]
      signal_syst_dn = combCard.glb_signal_syst_dn[key-1]

      print "key : %2d   value : %7.4f  -->  data : %2.0f  sumPred : %6.3f +%5.3f -%5.3f (+%5.3f -%5.3f)  <-->  signal : %6.3f +%5.3f -%5.3f (+%5.3f -%5.3f)" % (key-1, value, data_rate, tot_rate, tot_stat_up, tot_stat_dn, tot_syst_up, tot_syst_dn, signal_rate, signal_stat_up, signal_stat_dn, signal_syst_up, signal_syst_dn)
      output_res_file.write("key : %2d   value : %7.4f  -->  data : %2.0f  sumPred : %6.3f +%5.3f -%5.3f (+%5.3f -%5.3f)  <-->  signal : %6.3f +%5.3f -%5.3f (+%5.3f -%5.3f)\n" % (key-1, value, data_rate, tot_rate, tot_stat_up, tot_stat_dn, tot_syst_up, tot_syst_dn, signal_rate, signal_stat_up, signal_stat_dn, signal_syst_up, signal_syst_dn) )
      output_res_latex_file.write("      %2d &  %7.4f & %2.0f & %6.3f & $^{+%5.3f}_{-%5.3f}$ & $^{+%5.3f}_{-%5.3f}$ & %6.3f & $^{+%5.3f}_{-%5.3f}$ & $^{+%5.3f}_{-%5.3f}$ \\\\ \n" % (key-1, value, data_rate, tot_rate, tot_stat_up, tot_stat_dn, tot_syst_up, tot_syst_dn, signal_rate, signal_stat_up, signal_stat_dn, signal_syst_up, signal_syst_dn) )
      output_res_latex_file.write("\\hline\n")

   latex_end_line = """     \\end{tabular}
  \\end{center}
\\end{table}
"""
   output_res_latex_file.write(latex_end_line)
 
   output_res_latex_file.close() 
   output_res_file.close() 

#   for ch in range(1, int(combCard.channels)+1):
      
   lostle_file.close()
   hadtau_file.close()
   zinv_file.close()
   qcd_file.close()
   ttz_file.close()
   data_file.close()

if __name__ == "__main__":
   main()
