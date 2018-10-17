import sys
from inputFile import *
from smsPlotXSEC import *
from smsPlotCONT import *
from smsPlotBrazil import *
from smsPlotXSECaddRealExp import *

if __name__ == '__main__':
    # read input arguments
    filename = sys.argv[1]
    modelname = sys.argv[1].split("/")[-1].split("_")[0]
    analysisLabel = sys.argv[1].split("/")[-1].split("_")[1]
    outputname = sys.argv[2]

    # read the config file
    fileIN = inputFile(filename)

    # classic temperature histogra
    xsecPlot = smsPlotXSEC(modelname, fileIN.HISTOGRAM, fileIN.OBSERVED, fileIN.EXPECTED, fileIN.ENERGY, fileIN.LUMI, fileIN.PRELIMINARY, "", fileIN.SPECIAL)
    xsecPlot.Draw()

    inputfile = open(filename)
    for line in inputfile:
       tmpLINE =  line[:-1].split(" ")
       if tmpLINE[0] != "SPECIAL": continue
       foundCovered = False
       for objName in tmpLINE[2:]:
          if "diagonalCoverBand" in objName: foundCovered = True
    inputfile.close()
       
    if foundCovered:
       xsecPlot.Save("Covered_%sXSEC" %outputname)
    else:
       xsecPlot.Save("%sXSEC" %outputname)

    # only lines
#    contPlot = smsPlotCONT(modelname, fileIN.HISTOGRAM, fileIN.OBSERVED, fileIN.EXPECTED, fileIN.ENERGY, fileIN.LUMI, fileIN.PRELIMINARY, "", fileIN.SPECIAL)
#    contPlot.Draw()
#    contPlot.Save("%sCONT" %outputname)

    # brazilian flag (show only 1 sigma)
#    brazilPlot = smsPlotBrazil(modelname, fileIN.HISTOGRAM, fileIN.OBSERVED, fileIN.EXPECTED, fileIN.ENERGY, fileIN.LUMI, fileIN.PRELIMINARY, "", fileIN.SPECIAL)
#    brazilPlot.Draw()
#    brazilPlot.Save("%sBAND" %outputname)
