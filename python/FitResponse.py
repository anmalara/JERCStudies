from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

class FitResponse(VariablesBase):
    def __init__(self, sample):
        VariablesBase.__init__(self)
        self.outdir = self.Path_ANALYSIS+"python/Closure/"
        os.system("mkdir -p "+self.outdir)
        self.sample = sample
        self.color  = {"EGE": ROOT.kAzure+2,
                       "CB": ROOT.kOrange+1,
                       }
        self.flav_bins = ["uds", "c", "b", "g", "pu"]
        self.flav_bins = ["uds", "g", "all", "light", "heavy"]
        # self.flav_bins = ["all"]
        # self.flav_bins = ["uds"]
        self.eta_bins = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0]
        # self.eta_bins = [0.0, 0.5]
        # self.eta_bins = [2.5, 3.0]
        # self.pt_bins = [10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0]
        self.pt_bins = [30.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0]
        self.eta_bins_all = [-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, +0.000, +0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, +0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, +1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500, +2.650, +2.853, +2.964, +3.139, +3.314, +3.489, +3.664, +3.839, +4.013, +4.191, +4.363, +4.538, +4.716, +4.889, +5.191]
        self.pt_bins = [15.0, 17.0, 20.0, 23.0, 27.0, 30.0, 35.0, 40.0, 45.0, 57.0, 72.0, 90.0, 120.0, 150.0, 200.0, 300.0, 400.0, 550.0, 750.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0]

        self.algos = ["Puppi","CHS"]
        self.algos = ["CHS"]
        self.modes = ["nocuts", "weights", "cleaned"]
        self.modes = ["cleaned"]

        self.etapos = ["Neg_HF","Neg_EC","Neg_BA","Pos_BA","Pos_EC","Pos_HF"]

        for algo in self.algos:
            for mode in self.modes:
                for pos in self.etapos:
                    isPos = "Pos" in pos
                    isBA = "BA" in pos
                    isHF = "HF" in pos
                    self.eta_bins = self.eta_bins_all
                    self.eta_bins = FilterVector(self.eta_bins, 0,                                   "" if isPos else "invert")
                    self.eta_bins = FilterVector(self.eta_bins, (1 if isPos else -1)*ROOT.etaBarrel, "invert" if isPos==isBA else "")
                    self.eta_bins = FilterVector(self.eta_bins, (1 if isPos else -1)*ROOT.etaHF,      "" if isPos==isHF else "invert")
                    self.LoadHistos("Response_"+algo+"_"+pos+"_"+mode)
                    # self.Plot("Response_"+algo+"_"+mode+"_"+self.sample)

    def ExtractInfo(self, name):
        if not hasattr(self, 'quant'):
            self.quant   = array('d',[0.5])
            self.quant_y   = array('d',[0.5])
        mean = self.hists[name].GetMean()
        std = self.hists[name].GetStdDev()
        min = mean-2*std
        max = mean+2*std
        self.funcs[name] = ROOT.TF1(name+"gaus", "gaus(0)", min, max)
        self.funcs[name].SetParameters(mean,std)
        for nfit in range(100):
            if min>1 or max < 1: continue
            fitRes = self.hists[name].Fit(name+"gaus", "QMS", "", min, max)
            if self.funcs[name].GetNDF()== 0: continue
            chi2_red = self.funcs[name].GetChisquare()/self.funcs[name].GetNDF()
            if chi2_red>0.5 and chi2_red<2:
                break
            else:
                min = self.funcs[name].GetParameter(1) - 1.5*self.funcs[name].GetParameter(2)
                max = self.funcs[name].GetParameter(1) + 1.5*self.funcs[name].GetParameter(2)
        self.hists[name].GetQuantiles(1,self.quant_y,self.quant)
        return (self.funcs[name].GetParameter(1), mean, self.quant_y[0])

    def LoadHistos(self, folder):
        f_ = ROOT.TFile(self.Path_STORAGE+"ResponseClosure/"+self.PrefixrootFile+"MC."+self.sample+".root")
        self.hists = {}
        self.graphs = {}
        self.funcs = {}
        self.chi2 = {}
        self.errors = {}
        self.means = -1*np.ones((3,len(self.flav_bins),len(self.eta_bins_all),len(self.pt_bins)))
        for f, flav in enumerate(self.flav_bins):
            flav_name = "flav_"+flav
            for e in range(len(self.eta_bins)-1):
                eta_name = "eta_"+ROOT.GetStringFromFloat(self.eta_bins[e])+"to"+ROOT.GetStringFromFloat(self.eta_bins[e+1])
                pts = []
                means = []
                means_err = []
                for p in range(len(self.pt_bins)-1):
                    if self.pt_bins[p+1]<=30: continue
                    pt_name = "pt_"+ROOT.GetStringFromFloat(self.pt_bins[p])+"to"+ROOT.GetStringFromFloat(self.pt_bins[p+1])
                    name = "Resp_"+flav_name+eta_name+pt_name
                    self.hists[name] = f_.Get(folder+"/"+name)
                    if self.hists[name].GetEntries()==0 or self.hists[name].Integral()==0: continue
                    self.hists[name].SetDirectory(0)
                    self.hists[name].Scale(1./self.hists[name].Integral())
                    if self.hists[name].GetEntries()<100:
                        self.means[0,f,e,p] = 0
                        self.means[1,f,e,p] = 0
                        self.means[2,f,e,p] = 0
                        continue
                    gaus, mean, median = self.ExtractInfo(name)
                    if self.means[0,f,e,p] == -1:
                        self.means[0,f,e,p] = gaus
                    else:
                        print "ERROR"
                    pts.append(self.pt_bins[p])
                    means.append(median)
                    means_err.append(np.abs(median-gaus))
                    self.canv = tdrCanvas(folder+name, 0, 2, 0, 1.2*self.hists[name].GetMaximum(), "Response", "A.U.")
                    tex = rt.TLatex()
                    tex.SetNDC()
                    tex.SetTextFont(42)
                    tex.SetTextSize(0.035)
                    tex.DrawLatex(0.75,0.85, "gauss:  "+str(round(gaus,3)))
                    tex.DrawLatex(0.75,0.80, "mean:   "+str(round(mean,3)))
                    tex.DrawLatex(0.75,0.75, "median: "+str(round(median,3)))
                    tdrDraw(self.hists[name], "P", ROOT.kFullCircle, ROOT.kBlack, fstyle=0)
                    self.funcs[name].SetLineColor(ROOT.kRed+1)
                    self.funcs[name].Draw("same")
                    self.canv.SaveAs(self.outdir+folder+name+".pdf")

                if len(pts) ==0 : continue
                self.graphs[folder+flav_name+eta_name] = ROOT.TGraphErrors(len(pts), array('d',pts), array('d',means), array('d',np.zeros(len(pts))), array('d',means_err))
        f_.Close()
        print np.where(self.means == -1)

    def CreateCanvas(self, fname):
        PlotXMin = 11
        PlotXMax = 5500
        PlotYMin = 0.97
        PlotYMax = 1.03
        self.canv = tdrCanvas(fname, PlotXMin, PlotXMax, PlotYMin, PlotYMax, "p_{T}^{ptlc}", "Response")
        self.canv.SetLogx(True)
        GettdrCanvasHist(self.canv).GetXaxis().SetTitleOffset(1.0)
        GettdrCanvasHist(self.canv).GetXaxis().SetNoExponent()
        tex = rt.TLatex()
        tex.SetTextFont(GettdrCanvasHist(self.canv).GetXaxis().GetLabelFont())
        tex.SetTextSize(GettdrCanvasHist(self.canv).GetXaxis().GetLabelSize())
        tex.DrawLatex(25,PlotYMin-0.0037, "30")
        tex.DrawLatex(235,PlotYMin-0.0037, "300")
        tex.DrawLatex(2100,PlotYMin-0.0037, "3000")
        self.leg = tdrLeg(0.60,0.70,0.89,0.90, textSize=0.04)
        self.lines = {}
        for shift in [+0.01, 0.00, -0.01, +0.001, -0.001]:
            self.lines[shift] = ROOT.TLine(PlotXMin, 1+shift, PlotXMax, 1+shift)
            self.lines[shift].SetLineWidth(1)
            self.lines[shift].SetLineStyle(ROOT.kDotted if shift != 0 else ROOT.kDashed)
            self.lines[shift].SetLineColor(ROOT.kBlack)
            self.lines[shift].Draw("same")

    def Plot(self, fname):
        for flav in self.flav_bins:
            flav_name = "flav_"+flav
            TDR.extraText3 = ["flav = "+flav]
            self.CreateCanvas(fname+flav_name)
            for e in range(len(self.eta_bins)-1):
                eta_name = "eta_"+ROOT.GetStringFromFloat(self.eta_bins[e])+"to"+ROOT.GetStringFromFloat(self.eta_bins[e+1])
                name = flav_name+eta_name
                if self.eta_bins[e] == 0.0:
                    color = ROOT.kOrange+1
                if self.eta_bins[e] == 0.5:
                    color = ROOT.kGreen+2
                if self.eta_bins[e] == 1.0:
                    color = ROOT.kRed+1
                if self.eta_bins[e] == 1.5:
                    color = ROOT.kBlue-7
                if self.eta_bins[e] == 2.0:
                    color = ROOT.kViolet-7
                if self.eta_bins[e] == 2.5:
                    color = ROOT.kAzure+2
                if self.eta_bins[e] == 3.0:
                    color = ROOT.kBlack
                color = ROOT.kBlack
                tdrDraw(self.graphs[name], "PC", ROOT.kFullCircle, color, 1, color, 0, color)
                self.leg.AddEntry(self.graphs[name], str(self.eta_bins[e])+"< #eta < "+str(self.eta_bins[e+1]), "lp")
            self.canv.SaveAs(self.outdir+fname+flav+".pdf")

        for e in range(len(self.eta_bins)-1):
            eta_name = "eta_"+ROOT.GetStringFromFloat(self.eta_bins[e])+"to"+ROOT.GetStringFromFloat(self.eta_bins[e+1])
            TDR.extraText3 = [str(self.eta_bins[e])+"< #eta < "+str(self.eta_bins[e+1])]
            self.CreateCanvas(fname+eta_name)
            for flav in self.flav_bins:
                flav_name = "flav_"+flav
                name = flav_name+eta_name
                if "uds" in name: color = ROOT.kOrange+1
                if "c" in name: color = ROOT.kGreen+2
                if "b" in name: color = ROOT.kRed+1
                if "g" in name: color = ROOT.kBlue-7
                if "heavy" in name: color = ROOT.kViolet-7
                if "light" in name: color = ROOT.kAzure+2
                if "all" in name: color = ROOT.kBlack
                tdrDraw(self.graphs[name], "PC", ROOT.kFullCircle, color, 1, color, 0, color)
                self.leg.AddEntry(self.graphs[name], flav, "lp")
            self.canv.SaveAs(self.outdir+fname+eta_name+".pdf")


if __name__ == '__main__':

    samples = ["MC_DY_UL16APV", "MC_QCD_UL16nonAPV"]
    samples = ["MC_QCD_UL16nonAPV"]

    for sample in samples:
        PlotBkg = FitResponse(sample=sample)
    # PlotBkg.DoFits()
