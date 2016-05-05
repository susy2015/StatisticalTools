import ROOT as rt
from array import *
from sms import *
from smsPlotABS import *

# class producing the 2D plot with xsec colors
class plotBinaryLimits(smsPlotABS):

    def __init__(self, modelname, histo, obsLimits, expLimits, energy, lumi, preliminary, label, special):
        self.LABEL = label
        self.standardDef(modelname, histo, obsLimits, expLimits, energy, lumi, preliminary, special)
        # canvas for the plot
        self.c = rt.TCanvas("cCONT_%s" %label,"cCONT_%s" %label,1200,800)
        self.histo = histo['histogram']
        self.histo.GetZaxis().SetTitle("")
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
        NRGBs = 4
        NCont = 255
        stops = array("d",[0.00, 0.495, 0.505, 1.00])
        red= array("d",[0.20, 0.20, 0.50, 1.00])
        green = array("d",[ 0.10, 0.60, 0.15, 0.00])
        blue = array("d",[0.20, 0.70, 0.15, 0.00])
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
            
    def Draw(self):
        rt.gPad.SetLogz(0)
        self.emptyHisto.GetXaxis().SetRangeUser(self.model.Xmin, self.model.Xmax)
        self.emptyHisto.GetYaxis().SetRangeUser(self.model.Ymin, self.model.Ymax)
        self.emptyHisto.Draw()
        rt.gStyle.SetPaintTextFormat("4.2f")
        self.histo.SetMarkerSize(0.75)
        self.histo.GetZaxis().SetRangeUser(0, 2)
        self.histo.Draw("COLZ SAME TEXT")
#        self.DrawDiagonal()
        self.DrawLines()
        if self.EXP and 'nominal' in self.EXP: self.EXP['nominal'].SetLineWidth(2) 
        if self.EXP and 'nominal_base' in self.EXP: self.EXP['nominal_base'].SetLineWidth(2) 
        if self.EXP and 'nominal_extra' in self.EXP: self.EXP['nominal_extra'].SetLineWidth(2) 
        if self.OBS and 'nominal' in self.OBS: self.OBS['nominal'].SetLineWidth(2) 
        if self.OBS and 'nominal_base' in self.OBS: self.OBS['nominal_base'].SetLineWidth(2) 
        if self.OBS and 'nominal_extra' in self.OBS: self.OBS['nominal_extra'].SetLineWidth(2) 
        self.DrawSpecial()
        self.DrawText()
        self.DrawLegend()
        self.DrawPaletteLabel()
