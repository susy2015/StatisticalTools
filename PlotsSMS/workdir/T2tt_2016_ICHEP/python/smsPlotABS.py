import ROOT as rt
from array import *
from sms import *
from color import *
import CMS_lumi

class smsPlotABS(object):
    # modelname is the sms name (see sms.py)
    # histo is the 2D xsec map
    # obsLimits is a list of opbserved limits [NOMINAL, +1SIGMA, -1SIGMA]
    # expLimits is a list of expected limits [NOMINAL, +1SIGMA, -1SIGMA]
    # label is a label referring to the analysis (e.g. RA1, RA2, RA2b, etc)

    def __init__(self, modelname, histo, obsLimits, expLimits, energy, lumi, preliminary, label, special):
        self.LABEL = label
        self.standardDef(modelname, histo, obsLimits, expLimits, energy, lumi, preliminary, special)
        self.c = rt.TCanvas("cABS_%s" %label,"cABS_%s" %label,300,300)
        self.histo = histo

    def standardDef(self, modelname, histo, obsLimits, expLimits, energy, lumi, preliminary, special):
        # which SMS?
        self.model = sms(modelname)
        self.OBS = obsLimits
        self.EXP = expLimits
        self.lumi = lumi
        self.energy = energy
        self.preliminary = preliminary
        # create the reference empty histo
        self.emptyhisto = self.emptyHistogramFromModel()
        self.SPECIAL = special

    def emptyHistogramFromModel(self):
        self.emptyHisto = rt.TH2D("emptyHisto"+self.LABEL, "", 1, self.model.Xmin, self.model.Xmax, 1, self.model.Ymin, self.model.Ymax)
        
    # define the plot canvas
    def setStyle(self):
        # canvas style
        rt.gStyle.SetOptStat(0)
        rt.gStyle.SetOptTitle(0)        

        self.c.SetLogz()
        self.c.SetTickx(1)
        self.c.SetTicky(1)

        self.c.SetRightMargin(0.19)
        self.c.SetTopMargin(0.08)
        self.c.SetLeftMargin(0.14)
        self.c.SetBottomMargin(0.14)

        # set x axis
        self.emptyHisto.GetXaxis().SetLabelFont(42)
        self.emptyHisto.GetXaxis().SetLabelSize(0.035)
        self.emptyHisto.GetXaxis().SetTitleFont(42)
        self.emptyHisto.GetXaxis().SetTitleSize(0.05)
        self.emptyHisto.GetXaxis().SetTitleOffset(1.2)
        self.emptyHisto.GetXaxis().SetTitle(self.model.sParticle)
        #self.emptyHisto.GetXaxis().CenterTitle(True)

        # set y axis
        self.emptyHisto.GetYaxis().SetLabelFont(42)
        self.emptyHisto.GetYaxis().SetLabelSize(0.035)
        self.emptyHisto.GetYaxis().SetTitleFont(42)
        self.emptyHisto.GetYaxis().SetTitleSize(0.05)
        self.emptyHisto.GetYaxis().SetTitleOffset(1.3)
        self.emptyHisto.GetYaxis().SetTitle(self.model.LSP)
        #self.emptyHisto.GetYaxis().CenterTitle(True)
                
    def DrawText(self, pctYforLegend=0.75):
        #redraw axes
        self.c.RedrawAxis()
        # white background
        graphWhite = rt.TGraph(5)
        graphWhite.SetName("white")
        graphWhite.SetTitle("white")
        graphWhite.SetFillColor(rt.kWhite)
        graphWhite.SetFillStyle(1001)
        graphWhite.SetLineColor(rt.kBlack)
        graphWhite.SetLineStyle(1)
        graphWhite.SetLineWidth(3)
        graphWhite.SetPoint(0,self.model.Xmin, self.model.Ymax)
        graphWhite.SetPoint(1,self.model.Xmax, self.model.Ymax)
        if self.model.label2 == "":
           graphWhite.SetPoint(2,self.model.Xmax, self.model.Ymax*pctYforLegend)
           graphWhite.SetPoint(3,self.model.Xmin, self.model.Ymax*pctYforLegend)
        else:
           graphWhite.SetPoint(2,self.model.Xmax, self.model.Ymax*0.69)
           graphWhite.SetPoint(3,self.model.Xmin, self.model.Ymax*0.69)
        graphWhite.SetPoint(4,self.model.Xmin, self.model.Ymax)
        graphWhite.Draw("FSAME")
        graphWhite.Draw("LSAME")
        self.c.graphWhite = graphWhite
        CMS_lumi.writeExtraText = 0
        CMS_lumi.extraText = self.preliminary
        CMS_lumi.lumi_13TeV = self.lumi+" fb^{-1}"

        CMS_lumi.lumi_sqrtS = self.energy+" TeV"
        iPos=0
        CMS_lumi.CMS_lumi(self.c,4, iPos)
        # CMS LABEL
        textCMS = rt.TLatex(0.25,0.96,"  %s " %(self.preliminary))
        textCMS.SetNDC()
        textCMS.SetTextAlign(13)
        textCMS.SetTextFont(52)
        textCMS.SetTextSize(0.038)
        textCMS.Draw()
        self.c.textCMS = textCMS
        # MODEL LABEL
        if self.model.label2 == "":
           textModelLabel= rt.TLatex(0.15,0.90,"%s  NLO+NLL exclusion" %self.model.label)
           textModelLabel.SetNDC()
           textModelLabel.SetTextAlign(13)
           textModelLabel.SetTextFont(42)
           textModelLabel.SetTextSize(0.035)
           textModelLabel.Draw()
           self.c.textModelLabel = textModelLabel
        else:
           textModelLabel= rt.TLatex(0.15,0.91,"%s" %self.model.label)
           textModelLabel.SetNDC()
           textModelLabel.SetTextAlign(13)
           textModelLabel.SetTextFont(42)
           textModelLabel.SetTextSize(0.035)
           textModelLabel.Draw()
           self.c.textModelLabel = textModelLabel
#           if self.model.label2 == "NLO+NLL exclusion":
#              textModelLabel2= rt.TLatex(0.15,0.845,"%s" %self.model.label2)
           if self.model.label2 == "T2tb":
              textModelLabel2= rt.TLatex(0.15,0.845,"NLO+NLL exclusion")
           else:
              textModelLabel2= rt.TLatex(0.15,0.845,"%s    NLO+NLL exclusion" %self.model.label2)
           textModelLabel2.SetNDC()
           textModelLabel2.SetTextAlign(13)
           textModelLabel2.SetTextFont(42)
           textModelLabel2.SetTextSize(0.035)
           textModelLabel2.Draw()
           self.c.textModelLabel2 = textModelLabel2
        # NLO NLL XSEC
        textNLONLL= rt.TLatex(0.16,0.32,"NLO-NLL exclusion")
        textNLONLL.SetNDC()
        textNLONLL.SetTextAlign(13)
        textNLONLL.SetTextFont(42)
        textNLONLL.SetTextSize(0.04)
        textNLONLL.Draw()
        #self.c.textNLONLL = textNLONLL

    def Save(self,label):
        # save the output
        self.c.SaveAs("%s.pdf" %label)
        self.c.SaveAs("%s.png" %label)
        
    def DrawLegend(self):
        if(self.model.label2 == ""):
            offset = 0
        else:
            offset = -25
        xRange = self.model.Xmax-self.model.Xmin
        yRange = self.model.Ymax-self.model.Ymin
        
        LObs = rt.TGraph(2)
        LObs.SetName("LObs")
        LObs.SetTitle("LObs")
        if self.OBS: LObs.SetLineColor(color(self.OBS['colorLine']))
        LObs.SetLineStyle(1)
        LObs.SetLineWidth(4)
        LObs.SetMarkerStyle(20)
        LObs.SetPoint(0,self.model.Xmin+3*xRange/100, self.model.Ymax-1.35*yRange/100*10+offset)
        LObs.SetPoint(1,self.model.Xmin+10*xRange/100, self.model.Ymax-1.35*yRange/100*10+offset)

        LObsP = rt.TGraph(2)
        LObsP.SetName("LObsP")
        LObsP.SetTitle("LObsP")
        if self.OBS: LObsP.SetLineColor(color(self.OBS['colorLine']))
        LObsP.SetLineStyle(1)
        LObsP.SetLineWidth(2)
        LObsP.SetMarkerStyle(20)
        LObsP.SetPoint(0,self.model.Xmin+3*xRange/100, self.model.Ymax-1.20*yRange/100*10+offset)
        LObsP.SetPoint(1,self.model.Xmin+10*xRange/100, self.model.Ymax-1.20*yRange/100*10+offset)

        LObsM = rt.TGraph(2)
        LObsM.SetName("LObsM")
        LObsM.SetTitle("LObsM")
        if self.OBS: LObsM.SetLineColor(color(self.OBS['colorLine']))
        LObsM.SetLineStyle(1)
        LObsM.SetLineWidth(2)
        LObsM.SetMarkerStyle(20)
        LObsM.SetPoint(0,self.model.Xmin+3*xRange/100, self.model.Ymax-1.50*yRange/100*10+offset)
        LObsM.SetPoint(1,self.model.Xmin+10*xRange/100, self.model.Ymax-1.50*yRange/100*10+offset)

        textObs = rt.TLatex(self.model.Xmin+11*xRange/100, self.model.Ymax-1.50*yRange/100*10+offset, "Observed #pm 1 #sigma_{theory}")
        textObs.SetTextFont(42)
        textObs.SetTextSize(0.040)
        textObs.Draw()
        self.c.textObs = textObs

        if self.model.label2 == "T2tb":
           extraModelLabel1 = rt.TLatex(self.model.Xmin+65*xRange/100, self.model.Ymax-1.50*yRange/100*10+offset, "m_{#tilde{#chi}^{#pm}_{1}} - m_{#tilde{#chi}^{0}_{1}} = 5 GeV")
           extraModelLabel1.SetTextFont(42)
           extraModelLabel1.SetTextSize(0.030)
           extraModelLabel1.Draw()
           self.c.extraModelLabel1 = extraModelLabel1

        LExpP = rt.TGraph(2)
        LExpP.SetName("LExpP")
        LExpP.SetTitle("LExpP")
        if self.EXP: LExpP.SetLineColor(color(self.EXP['colorLine']))
        LExpP.SetLineStyle(7)
        LExpP.SetLineWidth(2)  
        LExpP.SetPoint(0,self.model.Xmin+3*xRange/100, self.model.Ymax-1.85*yRange/100*10+offset)
        LExpP.SetPoint(1,self.model.Xmin+10*xRange/100, self.model.Ymax-1.85*yRange/100*10+offset)

        LExp = rt.TGraph(2)
        LExp.SetName("LExp")
        LExp.SetTitle("LExp")
        if self.EXP: LExp.SetLineColor(color(self.EXP['colorLine']))
        LExp.SetLineStyle(7)
        LExp.SetLineWidth(4)
        LExp.SetPoint(0,self.model.Xmin+3*xRange/100, self.model.Ymax-2.00*yRange/100*10+offset)
        LExp.SetPoint(1,self.model.Xmin+10*xRange/100, self.model.Ymax-2.00*yRange/100*10+offset)
        
        LExpM = rt.TGraph(2)
        LExpM.SetName("LExpM")
        LExpM.SetTitle("LExpM")
        if self.EXP: LExpM.SetLineColor(color(self.EXP['colorLine']))
        LExpM.SetLineStyle(7)
        LExpM.SetLineWidth(2)  
        LExpM.SetPoint(0,self.model.Xmin+3*xRange/100, self.model.Ymax-2.15*yRange/100*10+offset)
        LExpM.SetPoint(1,self.model.Xmin+10*xRange/100, self.model.Ymax-2.15*yRange/100*10+offset)

        textExp = rt.TLatex(self.model.Xmin+11*xRange/100, self.model.Ymax-2.15*yRange/100*10+offset, "Expected #pm 1 #sigma_{experiment}")
        textExp.SetTextFont(42)
        textExp.SetTextSize(0.040)
        textExp.Draw()
        self.c.textExp = textExp

        if self.model.label2 == "T2tb":
           extraModelLabel2 = rt.TLatex(self.model.Xmin+65*xRange/100, self.model.Ymax-2.15*yRange/100*10+offset, "BR(#tilde{t} #rightarrow t #tilde{#chi}^{0}_{1}) = 50%")
           extraModelLabel2.SetTextFont(42)
           extraModelLabel2.SetTextSize(0.030)
           extraModelLabel2.Draw()
           self.c.extraModelLabel2 = extraModelLabel2

        LObs.Draw("LSAME")
        LObsM.Draw("LSAME")
        LObsP.Draw("LSAME")
        LExp.Draw("LSAME")
        LExpM.Draw("LSAME")
        LExpP.Draw("LSAME")
        
        self.c.LObs = LObs
        self.c.LObsM = LObsM
        self.c.LObsP = LObsP
        self.c.LExp = LExp
        self.c.LExpM = LExpM
        self.c.LExpP = LExpP

    def DrawDiagonal(self):
        diagonal = rt.TGraph(2, self.model.diagX, self.model.diagY)
        diagonal.SetName("diagonal")
        diagonal.SetFillColor(rt.kWhite)
        diagonal.SetLineColor(rt.kGray)
        diagonal.SetLineStyle(2)
#        diagonal.Draw("FSAME")
        diagonal.Draw("LSAME")
        self.c.diagonal = diagonal

        if self.model.modelname == "T2tt" or self.model.modelname == "T2tb":
           text_diag = rt.TLatex(self.model.diagX[0] + 500, self.model.diagX[0] + 500 - 170, "m_{#tilde{t}} = m_{t} + m_{#tilde{#chi}^{0}_{1}}")
           text_diag.SetTextFont(42)
           text_diag.SetTextSize(0.023)
           text_diag.SetTextAngle(60.00)
           text_diag.Draw()
           self.c.text_diag = text_diag
        
    def DrawLines(self):
        # observed
        if self.OBS:
           doOBS = True
           if 'nominal_base' in self.OBS:
              doOBS = False
              self.OBS['nominal_base'].SetLineColor(color(self.OBS['colorLine']))
              self.OBS['nominal_base'].SetLineStyle(1)
              self.OBS['nominal_base'].SetLineWidth(4)
           if 'nominal_extra' in self.OBS:
              doOBS = False
              self.OBS['nominal_extra'].SetLineColor(color(self.OBS['colorLine']))
              self.OBS['nominal_extra'].SetLineStyle(1)
              self.OBS['nominal_extra'].SetLineWidth(4)
           if doOBS and 'nominal' in self.OBS:
              self.OBS['nominal'].SetLineColor(color(self.OBS['colorLine']))
              self.OBS['nominal'].SetLineStyle(1)
              self.OBS['nominal'].SetLineWidth(4)
           # observed + 1sigma
           doOBSplus = True
           if 'plus_base' in self.OBS:
              doOBSplus = False
              self.OBS['plus_base'].SetLineColor(color(self.OBS['colorLine']))
              self.OBS['plus_base'].SetLineStyle(1)
              self.OBS['plus_base'].SetLineWidth(2)        
           if 'plus_extra' in self.OBS:
              doOBSplus = False
              self.OBS['plus_extra'].SetLineColor(color(self.OBS['colorLine']))
              self.OBS['plus_extra'].SetLineStyle(1)
              self.OBS['plus_extra'].SetLineWidth(2)        
           if doOBSplus and 'plus' in self.OBS:
              self.OBS['plus'].SetLineColor(color(self.OBS['colorLine']))
              self.OBS['plus'].SetLineStyle(1)
              self.OBS['plus'].SetLineWidth(2)        
           # observed - 1sigma
           doOBSminus = True
           if 'minus_base' in self.OBS:
              doOBSminus = False
              self.OBS['minus_base'].SetLineColor(color(self.OBS['colorLine']))
              self.OBS['minus_base'].SetLineStyle(1)
              self.OBS['minus_base'].SetLineWidth(2)        
           if 'minus_extra' in self.OBS:
              doOBSminus = False
              self.OBS['minus_extra'].SetLineColor(color(self.OBS['colorLine']))
              self.OBS['minus_extra'].SetLineStyle(1)
              self.OBS['minus_extra'].SetLineWidth(2)        
           if doOBSminus and 'minus' in self.OBS:
              self.OBS['minus'].SetLineColor(color(self.OBS['colorLine']))
              self.OBS['minus'].SetLineStyle(1)
              self.OBS['minus'].SetLineWidth(2)        
        if self.EXP:
           # expected + 1sigma
           doEXPplus = True
           if 'plus_base' in self.EXP:
              doEXPplus = False
              self.EXP['plus_base'].SetLineColor(color(self.EXP['colorLine']))
              self.EXP['plus_base'].SetLineStyle(7)
              self.EXP['plus_base'].SetLineWidth(2)
           if 'plus_extra' in self.EXP:
              doEXPplus = False
              self.EXP['plus_extra'].SetLineColor(color(self.EXP['colorLine']))
              self.EXP['plus_extra'].SetLineStyle(7)
              self.EXP['plus_extra'].SetLineWidth(2)
           if doEXPplus and 'plus' in self.EXP:  
              self.EXP['plus'].SetLineColor(color(self.EXP['colorLine']))
              self.EXP['plus'].SetLineStyle(7)
              self.EXP['plus'].SetLineWidth(2)                
           # expected
           doEXP = True
           if 'nominal_base' in self.EXP:
              doEXP = False
              self.EXP['nominal_base'].SetLineColor(color(self.EXP['colorLine']))
              self.EXP['nominal_base'].SetLineStyle(7)
              self.EXP['nominal_base'].SetLineWidth(4)
           if 'nominal_extra' in self.EXP:
              doEXP = False
              self.EXP['nominal_extra'].SetLineColor(color(self.EXP['colorLine']))
              self.EXP['nominal_extra'].SetLineStyle(7)
              self.EXP['nominal_extra'].SetLineWidth(4)
           if doEXP and 'nominal' in self.EXP:
              self.EXP['nominal'].SetLineColor(color(self.EXP['colorLine']))
              self.EXP['nominal'].SetLineStyle(7)
              self.EXP['nominal'].SetLineWidth(4)
           # expected - 1sigma
           doEXPminus = True
           if 'minus_base' in self.EXP:
              doEXPminus = False
              self.EXP['minus_base'].SetLineColor(color(self.EXP['colorLine']))
              self.EXP['minus_base'].SetLineStyle(7)
              self.EXP['minus_base'].SetLineWidth(2)
           if 'minus_extra' in self.EXP:
              doEXPminus = False
              self.EXP['minus_extra'].SetLineColor(color(self.EXP['colorLine']))
              self.EXP['minus_extra'].SetLineStyle(7)
              self.EXP['minus_extra'].SetLineWidth(2)
           if doEXPminus and 'minus' in self.EXP:  
              self.EXP['minus'].SetLineColor(color(self.EXP['colorLine']))
              self.EXP['minus'].SetLineStyle(7)
              self.EXP['minus'].SetLineWidth(2)                
        # DRAW LINES
        if self.OBS:
           if 'nominal_base' in self.OBS: self.OBS['nominal_base'].Draw("LSAME")
           if 'nominal_extra' in self.OBS: self.OBS['nominal_extra'].Draw("LSAME")
           if doOBS and 'nominal' in self.OBS: self.OBS['nominal'].Draw("LSAME")
           
           if 'plus_base' in self.OBS: self.OBS['plus_base'].Draw("LSAME")
           if 'plus_extra' in self.OBS: self.OBS['plus_extra'].Draw("LSAME")
           if doOBSplus and 'plus' in self.OBS: self.OBS['plus'].Draw("LSAME")
           
           if 'minus_base' in self.OBS: self.OBS['minus_base'].Draw("LSAME")
           if 'minus_extra' in self.OBS: self.OBS['minus_extra'].Draw("LSAME")
           if doOBSminus and 'minus' in self.OBS: self.OBS['minus'].Draw("LSAME")
        
        if self.EXP:   
           if 'nominal_base' in self.EXP: self.EXP['nominal_base'].Draw("LSAME")
           if 'nominal_extra' in self.EXP: self.EXP['nominal_extra'].Draw("LSAME")
           if doEXP and 'nominal' in self.EXP: self.EXP['nominal'].Draw("LSAME")
           
           if 'plus_base' in self.EXP: self.EXP['plus_base'].Draw("LSAME")
           if 'plus_extra' in self.EXP: self.EXP['plus_extra'].Draw("LSAME")
           if doEXPplus and 'plus' in self.EXP: self.EXP['plus'].Draw("LSAME")
           
           if 'minus_base' in self.EXP: self.EXP['minus_base'].Draw("LSAME")
           if 'minus_extra' in self.EXP: self.EXP['minus_extra'].Draw("LSAME")
           if doEXPminus and 'minus' in self.EXP: self.EXP['minus'].Draw("LSAME")
           
    def DrawSpecial(self):
        if "diagonalCoverBand" in self.SPECIAL:
           diagonalCoverBand = rt.TGraph(self.SPECIAL["diagonalCoverBand"])
           diagonalCoverBand.SetLineColor(rt.kWhite)
           diagonalCoverBand.SetFillColor(rt.kWhite)
           diagonalCoverBand.Draw("FSAME")
           self.c.diagonalCoverBand = diagonalCoverBand
           self.c.diagonal.Draw("LSAME")
           if self.model.modelname == "T2tt" or self.model.modelname == "T2tb":
              self.c.text_diag.Draw()
