from array import *

class sms():

    def __init__(self, modelname):
        if modelname.find("T1tttt") != -1: self.T1tttt()
        if modelname.find("T1ttbb") != -1: self.T1ttbb()
        if modelname.find("T5ttcc") != -1: self.T5ttcc()
        if modelname.find("T5ttttDM175") != -1: self.T5ttttDM175()
        if modelname.find("T5ttttdegen") != -1: self.T5ttttdegen()
        if modelname.find("T6ttWW") != -1: self.T6ttWW()
        if modelname.find("T1bbbb") != -1: self.T1bbbb()
        if modelname.find("T1qqqq") != -1: self.T1qqqq()
        if modelname.find("T2qq")   != -1: self.T2qq()
        if modelname.find("T2tt")   != -1: self.T2tt()
        if modelname.find("T2bb")   != -1: self.T2bb()
        if modelname.find("T2tb")   != -1: self.T2tb()


    def T1tttt(self):
        # model name
        self.modelname = "T1tttt"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow t #bar{t} "+lsp_s;
        self.label2= "";
        # scan range to plot
        self.Xmin = 700
        self.Xmax = 1950
        self.Ymin = 0
        self.Ymax = 1800
#        self.Ymax = 2100
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} [GeV]"
        # diagonal position: mLSP = mgluino - 2mtop 
        mW = 80
        self.diagX = array('d',[0,20000])
        self.diagY = array('d',[-mW, 20000-mW])        
        
    def T1ttbb(self):
        # model name
        self.modelname = "T1ttbb"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow tb#tilde{#chi}^{#pm}_{1}";
        self.label2= "m_{#tilde{#chi}^{#pm}_{1}}-m_{"+lsp_s+"}= 5 GeV";
        # scan range to plot
        self.Xmin = 700
        self.Xmax = 1950
        self.Ymin = 0
        self.Ymax = 1800
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} [GeV]"
        # diagonal position: mLSP = mgluino - 2mtop 
        mW = 80
        self.diagX = array('d',[0,20000])
        self.diagY = array('d',[-mW, 20000-mW])        
        
    def T5ttcc(self):
        # model name
        self.modelname = "T5ttcc"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow #tilde{t}t, #tilde{t} #rightarrow c"+lsp_s;
        self.label2= "";
        # scan range to plot
        self.Xmin = 650
        self.Xmax = 1700
        self.Ymin = 0
        self.Ymax = 1700
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} [GeV]"
        # diagonal position: mLSP = mgluino - 2mtop 
        mW = 80
        self.diagX = array('d',[0,20000])
        self.diagY = array('d',[-mW, 20000-mW])        
        
    def T5ttttDM175(self):
        # model name
        self.modelname = "T5ttttDM175"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow t #bar{t} "+lsp_s;
        self.label2= "m_{#tilde{t}}-m_{"+lsp_s+"}= 175 GeV";
        # scan range to plot
        self.Xmin = 650
        self.Xmax = 1750
        self.Ymin = 0
        self.Ymax = 1500
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} [GeV]"
        # diagonal position: mLSP = mgluino - 2mtop 
        mW = 80
        self.diagX = array('d',[0,20000])
        self.diagY = array('d',[-mW, 20000-mW])        
        
    def T5ttttdegen(self):
        # model name
        self.modelname = "T5ttttdegen"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow t #bar{t} "+lsp_s;
        self.label2= "m_{#tilde{t}}-m_{"+lsp_s+"}= 20 GeV";
        # scan range to plot
        self.Xmin = 700
        self.Xmax = 1950
        self.Ymin = 0
        self.Ymax = 1800
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} [GeV]"
        # diagonal position: mLSP = mgluino - 2mtop 
        mW = 80
        self.diagX = array('d',[0,20000])
        self.diagY = array('d',[-mW, 20000-mW])        
        
    def T6ttWW(self):
        # model name
        self.modelname = "T6ttWW"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{b}#tilde{b}, #tilde{b} #rightarrow tW"+lsp_s;
        self.label2= "m_{"+lsp_s+"} = 50 GeV";
        # scan range to plot
        self.Xmin = 300.0
        self.Xmax = 900.0
        self.Ymin = 75.0
        self.Ymax = 900.0
        self.Zmin = 0.02
        self.Zmax = 100
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{b}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{#tilde{#chi}^{#pm}_{1}}} [GeV]"
        # diagonal position: mLSP = mgluino - 2mtop 
        mW = 0
        self.diagX = array('d',[0,20000])
        self.diagY = array('d',[-mW, 20000-mW])        
        
    def T1bbbb(self):
        # model name
        self.modelname = "T1bbbb"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow b #bar{b} "+lsp_s;
        self.label2= "";
        # plot boundary. The top 1/4 of the y axis is taken by the legend
        self.Xmin = 600
        self.Xmax = 1950
        self.Ymin = 0
        self.Ymax = 1800
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} [GeV]"
        # diagonal position: mLSP = mgluino - 2mtop
        self.diagX = array('d',[0,20000])
        self.diagY = array('d',[0, 20000])

    def T1qqqq(self):
        # model name
        self.modelname = "T1qqqq"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow q #bar{q} "+lsp_s;
        self.label2= "";
        # plot boundary. The top 1/4 of the y axis is taken by the legend
        self.Xmin = 600
        self.Xmax = 1950
        self.Ymin = 0 
        self.Ymax = 1600
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} [GeV]"
        # diagonal position: mLSP = mgluino - 2mtop
        self.diagX = array('d',[0,20000])
        self.diagY = array('d',[0, 20000])

    def T2qq(self):
        # model name
        self.modelname = "T2qq"
        # decay chain
        self.label= "pp #rightarrow #tilde{q} #tilde{q}, #tilde{q} #rightarrow q #bar{q} #tilde{#chi}^{0}_{1}";
        self.label2= "";
        # plot boundary. The top 1/4 of the y axis is taken by the legend
        self.Xmin = 300
        self.Xmax = 1200
        self.Ymin = 50
        self.Ymax = 800
        # produce sparticle
        self.sParticle = "m_{gluino} (GeV)"
        # LSP
        self.LSP = "m_{LSP} (GeV)"
        # diagonal position: mLSP = mgluino - 2mtop
        self.diagX = array('d',[0,20000])
        self.diagY = array('d',[0, 20000])

    def T2tt(self):
        # model name
        self.modelname = "T2tt"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{t} #tilde{t}, #tilde{t} #rightarrow t "+lsp_s;
        self.label2= "";
        # plot boundary. The top 1/4 of the y axis is taken by the legend
#        self.Xmin = 112.5
#        self.Xmax = 912.5
#        self.Ymin = 12.5
#        self.Ymax = 537.5
        self.Xmin = 100.0
        self.Xmax = 1200.0
        self.Ymin = 0.00
        self.Ymax = 650.0
        self.Zmin = 0.0005
        self.Zmax = 100
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{t}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} [GeV]"
        # diagonal position: mLSP = mgluino - 2mtop
        mtop = 175
        self.diagX = array('d',[0,20000])
        self.diagY = array('d',[-mtop, 20000-mtop])

    def T2tb(self):
        # model name
        self.modelname = "T2tb"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label= "pp #rightarrow #tilde{t} #tilde{t}*, #tilde{t} #rightarrow b #tilde{#chi}^{#pm}_{1} #rightarrow b W^{#pm} #tilde{#chi}^{0}_{1} or #tilde{t} #rightarrow t #tilde{#chi}^{0}_{1}"
#        self.label2= "NLO+NLL exclusion";
        self.label2= "T2tb"
#        self.label= "pp #rightarrow #tilde{t} #tilde{t}, #tilde{t} #rightarrow t #tilde{#chi}^{0}_{1}, #tilde{t} #rightarrow b #tilde{#chi}_{1}^{#pm}";
        # plot boundary. The top 1/4 of the y axis is taken by the legend
        self.Xmin = 200.0
        self.Xmax = 900.0
        self.Ymin = 0.00
        self.Ymax = 500.0
        self.Zmin = 0.02
        self.Zmax = 100
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{t}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{"+lsp_s+"}} [GeV]"
        # diagonal position: mLSP = mgluino - 2mtop
        mtop = 175
        self.diagX = array('d',[0,20000])
        self.diagY = array('d',[-mtop, 20000-mtop])

    def T2bb(self):
        # model name
        self.modelname = "T2bb"
        # decay chain
        self.label= "pp #rightarrow #tilde{b} #tilde{b}, #tilde{b} #rightarrow b #tilde{#chi}^{0}_{1}";
        self.label2= "";
        # plot boundary. The top 1/4 of the y axis is taken by the legend
#        self.Xmin = 187.5
        self.Xmin = 312.5
        self.Xmax = 912.5
        self.Ymin = 12.5
        self.Ymax = 637.5
        # produce sparticle
        self.sParticle = "m_{sbottom} (GeV)"
        # LSP
        self.LSP = "m_{LSP} (GeV)"
        # diagonal position: mLSP = mgluino - 2mtop
        self.diagX = array('d',[0,20000])
        self.diagY = array('d',[0, 20000])
