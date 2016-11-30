import ROOT as rt
from array import *
from sms import *
from smsPlotABS import *

# class producing the 2D plot with xsec colors
class smsPlotXSECaddRealExp(smsPlotABS):

    def __init__(self, modelname, histo, obsLimits, expLimits, realExpLimits, energy, lumi, preliminary, label):
        self.standardDef(modelname, histo, obsLimits, expLimits, energy, lumi, preliminary)
        self.REALEXP = realExpLimits 
        self.LABEL = label
        # canvas for the plot
        self.c = rt.TCanvas("cCONT_%s" %label,"cCONT_%s" %label,600,600)
        self.histo = histo['histogram']
        # canvas style
        self.setStyle()
        self.setStyleCOLZ()

    # define the plot canvas
    def setStyleCOLZ(self):        
        # set z axis
        self.histo.GetZaxis().SetLabelFont(42)
        self.histo.GetZaxis().SetTitleFont(42)
        self.histo.GetZaxis().SetLabelSize(0.035)
        self.histo.GetZaxis().SetTitleSize(0.035)

        # define the palette for z axis
        NRGBs = 5
        NCont = 255
        stops = array("d",[0.00, 0.34, 0.61, 0.84, 1.00])
        red= array("d",[0.50, 0.50, 1.00, 1.00, 1.00])
        green = array("d",[ 0.50, 1.00, 1.00, 0.60, 0.50])
        blue = array("d",[1.00, 1.00, 0.50, 0.40, 0.50])
        rt.TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)
        rt.gStyle.SetNumberContours(NCont)
        
        self.c.cd()
        self.histo.Draw("colz")
        
        rt.gPad.Update()
        palette = self.histo.GetListOfFunctions().FindObject("palette")
        palette.SetX1NDC(1.-0.18)
        palette.SetY1NDC(0.14)
        palette.SetX2NDC(1.-0.13)
        palette.SetY2NDC(1.-0.08)
        palette.SetLabelFont(42)
        palette.SetLabelSize(0.035)

    def DrawPaletteLabel(self):
        textCOLZ = rt.TLatex(0.98,0.15,"95% C.L. upper limit on cross section (pb)")
        textCOLZ.SetNDC()
        #textCOLZ.SetTextAlign(13)
        textCOLZ.SetTextFont(42)
        textCOLZ.SetTextSize(0.045)
        textCOLZ.SetTextAngle(90)
        textCOLZ.Draw()
        self.c.textCOLZ = textCOLZ

    def DrawRealExpLegend(self):
        xRange = self.model.Xmax-self.model.Xmin
        yRange = self.model.Ymax-self.model.Ymin

        prefitLExpP = rt.TGraph(2)
        prefitLExpP.SetName("prefitLExpP")
        prefitLExpP.SetTitle("prefitLExpP")
        prefitLExpP.SetLineColor(color(self.REALEXP['colorLine']))
        prefitLExpP.SetLineStyle(5)
        prefitLExpP.SetLineWidth(2)
        prefitLExpP.SetPoint(0,self.model.Xmin+3*xRange/100, self.model.Ymax-2.50*yRange/100*10)
        prefitLExpP.SetPoint(1,self.model.Xmin+10*xRange/100, self.model.Ymax-2.50*yRange/100*10)

        prefitLExp = rt.TGraph(2)
        prefitLExp.SetName("prefitLExp")
        prefitLExp.SetTitle("prefitLExp")
        prefitLExp.SetLineColor(color(self.REALEXP['colorLine']))
        prefitLExp.SetLineStyle(5)
        prefitLExp.SetLineWidth(4)
        prefitLExp.SetPoint(0,self.model.Xmin+3*xRange/100, self.model.Ymax-2.65*yRange/100*10)
        prefitLExp.SetPoint(1,self.model.Xmin+10*xRange/100, self.model.Ymax-2.65*yRange/100*10)

        prefitLExpM = rt.TGraph(2)
        prefitLExpM.SetName("prefitLExpM")
        prefitLExpM.SetTitle("prefitLExpM")
        prefitLExpM.SetLineColor(color(self.REALEXP['colorLine']))
        prefitLExpM.SetLineStyle(5)
        prefitLExpM.SetLineWidth(2)
        prefitLExpM.SetPoint(0,self.model.Xmin+3*xRange/100, self.model.Ymax-2.80*yRange/100*10)
        prefitLExpM.SetPoint(1,self.model.Xmin+10*xRange/100, self.model.Ymax-2.80*yRange/100*10)

        prefittextExp = rt.TLatex(self.model.Xmin+11*xRange/100, self.model.Ymax-2.80*yRange/100*10, "Expected w/o post-fit #pm 1 #sigma_{experiment}")
        prefittextExp.SetTextFont(42)
        prefittextExp.SetTextSize(0.040)
        prefittextExp.Draw()
        self.c.prefittextExp = prefittextExp

        prefitLExp.Draw("LSAME")
        prefitLExpM.Draw("LSAME")
        prefitLExpP.Draw("LSAME")

        self.c.prefitLExp = prefitLExp
        self.c.prefitLExpM = prefitLExpM
        self.c.prefitLExpP = prefitLExpP

    def DrawRealExp(self):
        # expected + 1sigma
        self.REALEXP['plus'].SetLineColor(color(self.REALEXP['colorLine']))
        self.REALEXP['plus'].SetLineStyle(5)
        self.REALEXP['plus'].SetLineWidth(2)
        # expected
        self.REALEXP['nominal'].SetLineColor(color(self.REALEXP['colorLine']))
        self.REALEXP['nominal'].SetLineStyle(5)
        self.REALEXP['nominal'].SetLineWidth(4)
        # expected - 1sigma
        self.REALEXP['minus'].SetLineColor(color(self.REALEXP['colorLine']))
        self.REALEXP['minus'].SetLineStyle(5)
        self.REALEXP['minus'].SetLineWidth(2)
        # DRAW LINES
        self.REALEXP['nominal'].Draw("LSAME")
        self.REALEXP['plus'].Draw("LSAME")
        self.REALEXP['minus'].Draw("LSAME")
            
    def Draw(self):
        self.emptyHisto.GetXaxis().SetRangeUser(self.model.Xmin, self.model.Xmax)
        self.emptyHisto.GetYaxis().SetRangeUser(self.model.Ymin, self.model.Ymax)
        self.emptyHisto.Draw()
        self.histo.Draw("COLZSAME")
        self.DrawDiagonal()
        self.DrawLines()
        self.DrawRealExp()
        self.DrawText(pctYforLegend=0.70)
        self.DrawLegend()
        self.DrawRealExpLegend()
        self.DrawPaletteLabel()
        self.c.RedrawAxis()
