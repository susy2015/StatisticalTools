import sys
from inputFile import *
from plotBinaryLimits import *
import ROOT as rt

if __name__ == '__main__':
    # read input arguments
    filename = sys.argv[1]
    modelname = sys.argv[1].split("/")[-1].split("_")[0]
    analysisLabel = sys.argv[1].split("/")[-1].split("_")[1]
    outputname = sys.argv[2]

    # read the config file
    fileIN = inputFile(filename)

    # classic temperature histogra
    xsecPlot = plotBinaryLimits(modelname, fileIN.HISTOGRAM, fileIN.OBSERVED, fileIN.EXPECTED, fileIN.ENERGY, fileIN.LUMI, fileIN.PRELIMINARY, "", fileIN.SPECIAL)
    xsecPlot.Draw()

    inputfile = open(filename)
    for line in inputfile:
       tmpLINE =  line[:-1].split(" ")
       if tmpLINE[0] != "HISTOGRAM": continue
       if "exp" in tmpLINE[2] or "EXP" in tmpLINE[2]:
          xsecPlot.Save("EXP_%s_binary_XSEC" %outputname)
       else:
          xsecPlot.Save("OBS_%s_binary_XSEC" %outputname)
    inputfile.close()

    # only lines
#    contPlot = smsPlotCONT(modelname, fileIN.HISTOGRAM, fileIN.OBSERVED, fileIN.EXPECTED, fileIN.ENERGY, fileIN.LUMI, fileIN.PRELIMINARY, "", fileIN.SPECIAL)
#    contPlot.Draw()
#    contPlot.Save("%sCONT" %outputname)

    # brazilian flag (show only 1 sigma)
#    brazilPlot = smsPlotBrazil(modelname, fileIN.HISTOGRAM, fileIN.OBSERVED, fileIN.EXPECTED, fileIN.ENERGY, fileIN.LUMI, fileIN.PRELIMINARY, "", fileIN.SPECIAL)
#    brazilPlot.Draw()
#    brazilPlot.Save("%sBAND" %outputname)
