from Utils import *
from Evaluate_MCJER import Evaluate_MCJER

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText = "Work in progress"

def NSC_Modified(x, par):
    N = par[0]
    S = par[1]
    C = par[2]
    d = par[3]
    k = par[4]
    N_RC_MC = par[5]
    N_RC_Data = par[6]
    C_IC_MC = par[7]
    C_IC_Data = par[8]
    mode = par[9]

    N_term = (N*ROOT.TMath.Abs(N)-N_RC_MC*N_RC_MC)/(x[0]*x[0])
    S_term = S*S*ROOT.TMath.Power(x[0],d)
    C_term = C*C - C_IC_MC*C_IC_MC
    Data = ROOT.TMath.Sqrt(N_RC_Data*N_RC_Data/(x[0]*x[0]) + k*k*( N_term + S_term + C_term) + C_IC_Data*C_IC_Data)
    MC   = ROOT.TMath.Sqrt(N_RC_MC  *N_RC_MC  /(x[0]*x[0]) +     ( N_term + S_term + C_term) + C_IC_MC*C_IC_MC)

    if mode == 0: return MC
    elif mode == 1: return Data
    else: return Data/MC

class JERCombination(VariablesBase):
    def __init__(self, year = "UL18"):
        VariablesBase.__init__(self)
        self.year = year
        TDR.cms_lumi_TeV = TDR.commonScheme["legend"][self.year]+" Legacy, "+commonScheme["lumi"][self.year]+" fb^{-1}"
        self.moduleName = self.__class__.__name__
        self.inpdir = self.Path_ANALYSIS+"StorageArea/"+self.moduleName+"/"
        self.outdir = self.Path_ANALYSIS+"python/"+self.moduleName+"/"
        os.system("mkdir -p "+self.outdir)
        self.colors = {"MC": ROOT.kRed+1, "Data": ROOT.kBlue+1, "free":ROOT.kBlack,
                       "rho0":  ROOT.kOrange+2,
                       "rho10": ROOT.kOrange+1,
                       "rho15": ROOT.kGreen+1,
                       "rho20": ROOT.kAzure+2,
                       "rho30": ROOT.kBlue-7,
                       "rho35": ROOT.kViolet,
                       "rho40": ROOT.kViolet-4,
                       "eta0":  ROOT.kRed,
                       "eta1":  ROOT.kRed+1,
                       "eta2":  ROOT.kRed+2,
                       "eta3":  ROOT.kOrange,
                       "eta4":  ROOT.kOrange+1,
                       "eta5":  ROOT.kOrange+2,
                       "eta6":  ROOT.kGreen,
                       "eta7":  ROOT.kGreen+1,
                       "eta8":  ROOT.kGreen+2,
                       "eta9":  ROOT.kGreen+3,
                       "eta10": ROOT.kBlue,
                       "eta11": ROOT.kBlue+1,
                       "eta12": ROOT.kBlue-7,
                       "eta13": ROOT.kAzure+2,
                       "eta14": ROOT.kAzure+7,
                       "eta15": ROOT.kAzure-7,
                       "eta16": ROOT.kViolet,
                       "eta17": ROOT.kViolet-2,
                       "eta18": ROOT.kViolet-4,
                       }
        self.fit_min, self.fit_max = (8., 5000.)
        self.fit_min2 = 15.
        self.inputs = {
            "N":           {"ID": 0, "vec": np.zeros(len(self.etaBinsCommon))},
            "S":           {"ID": 1, "vec": np.zeros(len(self.etaBinsCommon))},
            "C":           {"ID": 2, "vec": np.zeros(len(self.etaBinsCommon))},
            "d":           {"ID": 3, "vec": np.zeros(len(self.etaBinsCommon))},
            "k":           {"ID": 4, "vec": np.zeros(len(self.etaBinsCommon))},
            "N_{RC,MC}":   {"ID": 5, "vec": np.zeros(len(self.etaBinsCommon))},
            "N_{RC,Data}": {"ID": 6, "vec": np.zeros(len(self.etaBinsCommon))},
            "C_{IC,MC}":   {"ID": 7, "vec": np.zeros(len(self.etaBinsCommon))},
            "C_{IC,Data}": {"ID": 8, "vec": np.zeros(len(self.etaBinsCommon))},
        }
        self.Npar = len(self.inputs)+1
        self.JER = Evaluate_MCJER(GetJERfile(self.JERversions[self.year]))
        os.system("mkdir -p "+self.outdir)
        self.ExtractInfo()
        self.DoFits()
        self.DoFinalPlot()
        # print self.inputs.keys()

    def ExtractInfo(self):
        self.functions = {}
        self.ExtractRC()
        self.ExtractCterm()
        # self.ExtractDijet(fname="dijet_andrea")
        self.ExtractDijet()

    def ExtractRC(self, fname="RC"):
        f_ = ROOT.TFile(self.inpdir+fname+".root")
        for mode in ["MC", "Data"]:
            name = "N_{RC,"+("MC" if "MC" in mode else "Data")+"}"
            graph = f_.Get(mode+"/RMS")
            for p, etaRef in enumerate(self.etaBinsCommon):
                eta, N = (ROOT.Double(0),ROOT.Double(0))
                graph.GetPoint(int(p),eta,N)
                if round(eta,4) != round(etaRef,4): raise Exception("Something is not ok. Fix me!"+str(eta)+"!="+str(etaRef))
                self.inputs[name][p] = N
        f_.Close()

    def ExtractCterm(self, fname="Cterm"):
        f_ = ROOT.TFile(self.inpdir+fname+".root")
        for mode in ["mc", "data"]:
            name = "C_{IC,"+("MC" if "mc" in mode else "Data")+"}"
            h_ = f_.Get("jerc_rms_"+mode)
            for p, etaRef in enumerate(self.etaBinsCommon):
                if round(h_.GetBinCenter(p+1),4) != round(etaRef,4): raise Exception("Something is not ok. Fix me!"+str(h_.GetBinCenter(p+1))+"!="+str(etaRef))
                self.inputs[name][p] = h_.GetBinContent(p+1)/100
        f_.Close()

    def ExtractDijet(self, fname="dijet"):
        f_ = ROOT.TFile(self.inpdir+fname+".root")
        for n in range(len(self.etaBinsCommon)):
            dijetbin = 0
            if n<=2    : dijetbin = 1
            elif n<=5  : dijetbin = n-1
            elif n==6  : dijetbin = 4
            elif n<=9  : dijetbin = n-2
            elif n==10 : dijetbin = 7
            elif n<=15 : dijetbin = n-3
            else: dijetbin = 13
            name = "dijet"+FloatToString(self.etaBinsEdges[n])+"to"+FloatToString(self.etaBinsEdges[n+1])
            for mode in ["MC", "Data"]:
                if "andrea" in fname:
                    self.inputs[name+mode] = f_.Get((mode if "MC" in mode else "data")+"_JER_correlated_FE"+str(dijetbin))
                else:
                    self.inputs[name+mode] = f_.Get((mode if "MC" in mode else "data")+"_JER_standard_FE"+str(dijetbin))
                self.inputs[name+mode].SetDirectory(0)
                self.functions[name+mode+"free"] = self.inputs[name+mode].GetListOfFunctions().FindObject("mcFIT" if "MC" in mode else "dtFIT")
        f_.Close()

    def DoFits(self):
        for n in range(len(self.etaBinsCommon)):
            name = "dijet"+FloatToString(self.etaBinsEdges[n])+"to"+FloatToString(self.etaBinsEdges[n+1])
            TDR.extraText3 = []
            TDR.extraText3.append("AK4, PF+CHS")
            TDR.extraText3.append(str(self.etaBinsEdges[n])+" < |#eta| < "+str(self.etaBinsEdges[n+1]))
            self.CreateCanvas(self.moduleName+name)
            lines = {}

            for mode in ["MC", "Data"]:
                hist = self.inputs[name+mode]
                isMC = "MC" in mode
                tdrDraw(hist, "P", ROOT.kFullSquare if isMC else ROOT.kFullCircle, self.colors[mode], 1, self.colors[mode], 0, self.colors[mode])
                self.leg.AddEntry(hist, mode+" dijet ", "p")
                self.functions[name+mode+"free"].SetLineStyle(ROOT.kDashed)
                self.functions[name+mode+"free"].SetLineColor(self.colors["free"])
                self.functions[name+mode+"fix"]   = ROOT.TF1(name+mode+"fix", NSC_Modified, self.fit_min2, self.fit_max, self.Npar)

                if isMC:
                    self.functions[name+mode+"fix"].SetParameter(0, +5.00)
                    self.functions[name+mode+"fix"].SetParameter(1, +1.00)
                    self.functions[name+mode+"fix"].SetParameter(2, +0.05)
                    self.functions[name+mode+"fix"].SetParameter(3, -0.80)

                    self.functions[name+mode+"fix"].SetParLimits(0, -20.000, 20.0)
                    self.functions[name+mode+"fix"].SetParLimits(1,  +0.000,  2.0)
                    self.functions[name+mode+"fix"].SetParLimits(2,  +0.015,  1.0)
                    self.functions[name+mode+"fix"].SetParLimits(3,  -3.000,  0.0)

                    self.functions[name+mode+"fix"].FixParameter(4, 1.)
                    mode = 0
                else:
                    self.functions[name+mode+"fix"].FixParameter(0, self.functions[name+"MC"+"fix"].GetParameter(0))
                    self.functions[name+mode+"fix"].FixParameter(1, self.functions[name+"MC"+"fix"].GetParameter(1))
                    self.functions[name+mode+"fix"].FixParameter(2, self.functions[name+"MC"+"fix"].GetParameter(2))
                    self.functions[name+mode+"fix"].FixParameter(3, self.functions[name+"MC"+"fix"].GetParameter(3))

                    self.functions[name+mode+"fix"].SetParameter(4, +1.10)
                    self.functions[name+mode+"fix"].SetParLimits(4, +0.900,  2.0)
                    mode = 1

                self.functions[name+mode+"fix"].FixParameter(5, self.inputs["N_{RC,MC}"][n])
                self.functions[name+mode+"fix"].FixParameter(6, self.inputs["N_{RC,"+mode+"}"][n])
                self.functions[name+mode+"fix"].FixParameter(7, self.inputs["C_{IC,MC}"][n])
                self.functions[name+mode+"fix"].FixParameter(8, self.inputs["C_{IC,"+mode+"}"][n])
                self.functions[name+mode+"fix"].FixParameter(9, mode)

                hist.Fit(self.functions[name+mode+"fix"], "RMQ+")
                lines[mode]  = ROOT.TLine(1000.,self.functions[name+"MC"+"fix"].GetParameter(2),self.fit_max,self.inputs["C_{IC,"+mode+"}"][n])
                lines[mode].SetLineStyle(rt.kSolid if isMC else rt.kDashed)
                lines[mode].Draw("same")
                if not isMC:
                    self.inputs["N"][n] = self.functions[name+mode+"fix"].GetParameter(0)
                    self.inputs["S"][n] = self.functions[name+mode+"fix"].GetParameter(1)
                    self.inputs["C"][n] = self.functions[name+mode+"fix"].GetParameter(2)
                    self.inputs["d"][n] = self.functions[name+mode+"fix"].GetParameter(3)
                    self.inputs["k"][n] = self.functions[name+mode+"fix"].GetParameter(4)

                # for rho in [0, 10, 15, 20, 30, 35, 40]:
                for rho in [0, 20, 40]:
                    self.JER.setEta(self.GetEtaBinCenter(self.etaBinsEdges[n]))
                    self.JER.setRho(rho)
                    pars = self.JER.GetParameters()

                    lines[mode+"rho"+str(rho)]  = ROOT.TLine(1000.,pars[2],self.fit_max,pars[2])
                    lines[mode+"rho"+str(rho)].SetLineStyle(rt.kSolid if isMC else rt.kDashed)
                    lines[mode+"rho"+str(rho)].Draw("same")

                    self.functions[name+mode+"ref"+str(rho)] = ROOT.TF1(name+mode+"ref"+str(rho), NSC_Modified, self.fit_min2, self.fit_max, self.Npar)
                    for i, par in enumerate(pars):
                        self.functions[name+mode+"ref"+str(rho)].FixParameter(i, par)
                    self.functions[name+mode+"ref"+str(rho)].FixParameter(5, self.inputs["N_{RC,MC}"][n])
                    self.functions[name+mode+"ref"+str(rho)].FixParameter(6, self.inputs["N_{RC,"+mode+"}"][n])
                    self.functions[name+mode+"ref"+str(rho)].FixParameter(7, self.inputs["C_{IC,MC}"][n])
                    self.functions[name+mode+"ref"+str(rho)].FixParameter(8, self.inputs["C_{IC,"+mode+"}"][n])

                    if isMC:
                        self.functions[name+mode+"ref"+str(rho)].FixParameter(4, 1.)
                        self.functions[name+mode+"ref"+str(rho)].FixParameter(9, 0)
                    else:
                        self.functions[name+mode+"ref"+str(rho)].SetParameter(4, +1.10)
                        self.functions[name+mode+"ref"+str(rho)].SetParLimits(4, +0.900,  2.0)
                        self.functions[name+mode+"ref"+str(rho)].FixParameter(9, 1)
                        self.functions[name+mode+"ref"+str(rho)].SetLineStyle(ROOT.kDashed)
                        hist.Fit(self.functions[name+mode+"ref"+str(rho)], "RMQ+")

                    self.functions[name+mode+"ref"+str(rho)].SetLineColor(self.colors["rho"+str(rho)])
                    self.functions[name+mode+"ref"+str(rho)].Draw("same")

                self.functions[name+mode+"fix"].SetLineColor(self.colors[mode])
                self.functions[name+mode+"fix"].Draw("same")
                self.leg.AddEntry(self.functions[name+mode+"fix"], "NC fix", "l")
                self.functions[name+mode+"free"].Draw("same")
                self.leg.AddEntry(self.functions[name+mode+"free"], mode+" no fix", "l")

            self.canv.SaveAs(self.outdir+name+".pdf")
        # print "N", self.inputs["N"]
        # print "S", self.inputs["S"]
        # print "C", self.inputs["C"]
        # print "d", self.inputs["d"]
        # print "k", self.inputs["k"]

    def CreateCanvas(self, fname):
        PlotXMin = self.fit_min
        PlotXMax = self.fit_max
        PlotYMin = 0.
        PlotYMax = 0.30
        self.canv = tdrCanvas(fname, PlotXMin, PlotXMax, PlotYMin, PlotYMax, "p_{T}^{jet}", "#sigma")
        self.leg = tdrLeg(0.60,0.70,0.89,0.90, textSize=0.04)
        self.canv.SetLogx(True)

    def DoFinalPlot(self):
        PlotXMin = self.fit_min
        PlotXMax = self.fit_max
        # PlotYMin = 0.6
        # PlotYMax = 1.70
        PlotYMin = 0.9
        PlotYMax = 1.30
        TDR.extraText3 = []
        self.canv = tdrCanvas(self.moduleName, PlotXMin, PlotXMax, PlotYMin, PlotYMax, "p_{T}^{jet}", "SF")
        self.canv.SetLogx(True)
        GettdrCanvasHist(self.canv).GetXaxis().SetTitleOffset(1.0)
        GettdrCanvasHist(self.canv).GetXaxis().SetNoExponent()
        tex = rt.TLatex()
        tex.SetTextFont(GettdrCanvasHist(self.canv).GetXaxis().GetLabelFont())
        tex.SetTextSize(GettdrCanvasHist(self.canv).GetXaxis().GetLabelSize())
        # To be generalised
        tex.DrawLatex(25,PlotYMin-0.04, "30")
        tex.DrawLatex(235,PlotYMin-0.04, "300")
        tex.DrawLatex(2100,PlotYMin-0.04, "3000")
        # self.leg = tdrLeg(0.40,0.55,0.90,0.88, textSize=0.03)
        self.leg = tdrLeg(0.40,0.65,0.90,0.88, textSize=0.03)
        self.leg.SetNColumns(3)

        for n in range(len(self.etaBinsCommon)):
            if self.etaBinsEdges[n+1]>2.6: continue

            name = "dijet"+FloatToString(self.etaBinsEdges[n])+"to"+FloatToString(self.etaBinsEdges[n+1])
            mode = "ratio"
            self.functions[name+mode] = ROOT.TF1(name+mode, NSC_Modified, self.fit_min, self.fit_max, self.Npar)
            for par in range(self.functions[name+mode].GetNpar()):
                self.functions[name+mode].FixParameter(par, self.functions[name+"Data"+"fix"].GetParameter(par))
            self.functions[name+mode].FixParameter(9, 2)
            self.functions[name+mode].SetLineColor(self.colors["eta"+str(n)])
            self.functions[name+mode].Draw("same")
            self.leg.AddEntry(self.functions[name+mode], str(round(self.etaBinsEdges[n],1))+" < |#eta| < "+str(round(self.etaBinsEdges[n+1],1)), "l")
        self.canv.SaveAs(self.outdir+self.moduleName+".pdf")


if __name__ == '__main__':

    Comb = JERCombination()
    # PlotBkg.DoFits()
