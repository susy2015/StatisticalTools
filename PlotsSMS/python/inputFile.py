import sys
import ROOT as rt

class inputFile():

    def __init__(self, fileName):
        self.HISTOGRAM = self.findHISTOGRAM(fileName)
        self.EXPECTED = self.findEXPECTED(fileName)
        self.OBSERVED = self.findOBSERVED(fileName)
        self.LUMI = self.findATTRIBUTE(fileName, "LUMI")
        self.ENERGY = self.findATTRIBUTE(fileName, "ENERGY")
        self.PRELIMINARY = self.findATTRIBUTE(fileName, "PRELIMINARY")
        self.REALEXPECTED = self.findREALEXPECTED(fileName)

    def findATTRIBUTE(self, fileName, attribute):
        fileIN = open(fileName)        
        for line in fileIN:
            tmpLINE =  line[:-1].split(" ")
            if tmpLINE[0] != attribute: continue
            fileIN.close()
            return tmpLINE[1]

    def findHISTOGRAM(self, fileName):
        fileIN = open(fileName)        
        for line in fileIN:
            tmpLINE =  line[:-1].split(" ")
            if tmpLINE[0] != "HISTOGRAM": continue
            fileIN.close()
            rootFileIn = rt.TFile.Open(tmpLINE[1])
            return {'histogram': rootFileIn.Get(tmpLINE[2])}
            
    def findEXPECTED(self, fileName):
        fileIN = open(fileName)        
        for line in fileIN:
            tmpLINE =  line[:-1].split(" ")
            if tmpLINE[0] != "EXPECTED": continue
            fileIN.close()
            rootFileIn = rt.TFile.Open(tmpLINE[1])
            return {'nominal': rootFileIn.Get(tmpLINE[2]),
                    'plus': rootFileIn.Get(tmpLINE[3]),
                    'minus': rootFileIn.Get(tmpLINE[4]),
                    'colorLine': tmpLINE[5],
                    'colorArea': tmpLINE[6]}

    def findOBSERVED(self, fileName):
        fileIN = open(fileName)        
        for line in fileIN:
            tmpLINE =  line[:-1].split(" ")
            if tmpLINE[0] != "OBSERVED": continue
            fileIN.close()
            rootFileIn = rt.TFile.Open(tmpLINE[1])
            return {'nominal': rootFileIn.Get(tmpLINE[2]),
                    'plus': rootFileIn.Get(tmpLINE[3]),
                    'minus': rootFileIn.Get(tmpLINE[4]),
                    'colorLine': tmpLINE[5],
                    'colorArea': tmpLINE[6]}

    def findREALEXPECTED(self, fileName):
        fileIN = open(fileName)        
        for line in fileIN:
            tmpLINE =  line[:-1].split(" ")
            if tmpLINE[0] != "REALEXPECTED": continue
            fileIN.close()
            rootFileIn = rt.TFile.Open(tmpLINE[1])
            return {'nominal': rootFileIn.Get(tmpLINE[2]),
                    'plus': rootFileIn.Get(tmpLINE[3]),
                    'minus': rootFileIn.Get(tmpLINE[4]),
                    'colorLine': tmpLINE[5],
                    'colorArea': tmpLINE[6]}
