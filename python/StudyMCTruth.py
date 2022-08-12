from Utils import *
from GraphHistUtils import *

tdr.writeExtraText = True
tdr.extraText = 'Work in progress'

debug = False
postfix = ''

def isMC(name):    return 'MC' in name
def isData(name):  return 'Data' in name
def isRatio(name): return 'ratio' in name

def NSC_2D(x, par):
    pt = x[0]
    mu = x[1]
    N_term = (par[0]*par[0]+par[1]*par[1]*mu)/(pt*pt)
    S_term = par[2]*par[2]*ROOT.TMath.Power(pt,par[3])
    C_term = par[4]*par[4]
    return ROOT.TMath.Sqrt(N_term + S_term + C_term)

def NSC_1D(x, par):
    pt = x[0]
    mu = par[5]
    N_term = (par[0]*par[0]+par[1]*par[1]*mu)/(pt*pt)
    S_term = par[2]*par[2]*ROOT.TMath.Power(pt,par[3])
    C_term = par[4]*par[4]
    return ROOT.TMath.Sqrt(N_term + S_term + C_term)



def NSC_kfactor(x, par):
    # N_term = (par[0]*ROOT.TMath.Abs(par[0]))/(x[0]*x[0])
    N_term = (par[0]*par[0])/(x[0]*x[0])
    S_term = par[1]*par[1]*ROOT.TMath.Power(x[0],par[3])
    C_term = par[2]*par[2]

    N_RC_MC   = par[4]*par[4]/(x[0]*x[0])
    N_RC_Data = par[5]*par[5]/(x[0]*x[0])
    C_IC_MC   = par[6]*par[6]
    C_IC_Data = par[7]*par[7]

    kN = par[8]*par[8]
    kS = par[9]*par[9]
    kC = par[10]*par[10]

    isMC = par[11]

    Data = ROOT.TMath.Sqrt(N_RC_Data + kN*N_term + kS*S_term + kC*C_term + C_IC_Data)
    MC   = ROOT.TMath.Sqrt(N_RC_MC   +    N_term +    S_term +    C_term + C_IC_MC)

    if isMC == int(True): return MC
    elif isMC == int(False): return Data
    else: return Data/MC

class StudyMCTruth(InputBase):
    def __init__(self):
        InputBase.__init__(self, inputdir = os.path.join('Path_JERC','Summer20year','MC_JER'))
        self.LoadFiles()
        self.DoPlots()

    def LoadFiles(self, fname='JER_MCtruth', hname='RelResVsJetPt_JetEtaeta_mintoeta_max'):
        self.inputs = {'Mu_avg':{}, 'Mu_bin': {}, '2DGraphs':{}}
        self.funcs = {}
        eta = 0.261
        for year in self.years['UL']:
            inputdir = self.inputdir.replace('year', year).replace('nonAPV','')
            f_bin = ROOT.TFile(os.path.join(inputdir,fname+'_'+year.replace('nonAPV','')+'.root'))
            f_avg = ROOT.TFile(os.path.join(inputdir+'_avg_mu',fname+'_avg_mu_'+year.replace('nonAPV','')+'.root'))
            for type in ['ak4pfchs', 'ak4puppi']:
                for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                    eta_min, eta_max = (str(JERC_Constants.GetEtaBinEdgeMin(etaRef)),str(JERC_Constants.GetEtaBinEdgeMax(etaRef)))
                    if eta_min =='0.0': eta_min='0'
                    name = year+type+eta_min+eta_max
                    self.inputs['Mu_avg'][name] = f_avg.Get(type+'l1l2l3/'+hname.replace('eta_min',eta_min).replace('eta_max',eta_max)+'_Mu0to60')
                    self.inputs['Mu_avg'][name].GetListOfFunctions().RemoveLast()
                    x, y, z, ex, ey, ez = ([],[],[],[],[],[])
                    mu_bins = [0,10,20,30,40,50]
                    if 'UL16APV' in year: mu_bins = [0,10,20,30,40]
                    for mu in mu_bins:
                        name_mu = str(mu)+'to'+str(mu+10)
                        self.inputs['Mu_bin'][name+name_mu] = f_bin.Get(type+'l1l2l3/'+hname.replace('eta_min',eta_min).replace('eta_max',eta_max)+'_Mu'+name_mu)
                        self.inputs['Mu_bin'][name+name_mu].GetListOfFunctions().RemoveLast()
                        for i in range(self.inputs['Mu_bin'][name+name_mu].GetN()):
                            pt, jer = (ROOT.Double(0),ROOT.Double(0))
                            self.inputs['Mu_bin'][name+name_mu].GetPoint(int(i),pt,jer)
                            if (pt * math.cosh(eta)> 6500.): continue
                            if (pt > 3000): continue
                            x.append(pt)
                            y.append(mu+5)
                            z.append(jer)
                            ex.append(pt/20)
                            ey.append(5)
                            ez.append(self.inputs['Mu_bin'][name+name_mu].GetErrorY(i))
                    self.inputs['2DGraphs'][name] = ROOT.TGraph2DErrors(len(x), array('d',x), array('d',y),array('d',z),array('d',ex),array('d',ey),array('d',ez))
                    self.inputs['2DGraphs'][name].SetDirectory(0)
            f_bin.Close()
            f_avg.Close()

    def Plot2DGraph(self):
        SetAlternative2DColor()
        f_ = ROOT.TFile(os.path.join(self.outdir,'JER_2D.root'), 'recreate')
        # funcs = {}
        tdr.extraText3 = []
        for year in self.years['UL']:
            for type in ['ak4pfchs', 'ak4puppi']:
                for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                    eta_min, eta_max = (str(JERC_Constants.GetEtaBinEdgeMin(etaRef)),str(JERC_Constants.GetEtaBinEdgeMax(etaRef)))
                    if eta_min =='0.0': eta_min='0'
                    name = year+type+eta_min+eta_max
                    graph = self.inputs['2DGraphs'][name]
                    canv_name = 'JER_JetEta'+JERC_Constants.BinToString(float(eta_min),float(eta_max))+'_2D'
                    canv = tdrCanvas(canv_name, 0, 4500, 0, 100, 'p_{T}', 'mu', square=kSquare, is2D=True)
                    canv.SetLogx(True)
                    canv.cd()
                    graph.Draw('lego')
                    func = ROOT.TF2(name+'func2D', NSC_2D, 0, 3000, 0,60,5)
                    func.SetParLimits(0, -10,10)
                    func.SetParLimits(1,   0, 1)
                    func.SetParLimits(2,   0, 1)
                    func.SetParLimits(3,  -2, 0)
                    func.SetParLimits(4, 1e-03, 1e-01)
                    func.SetParameter(4, 1e-02)
                    func.FixParameter(0, 1.47/0.87)
                    # func.FixParameter(1, 0.526/0.87)
                    # graph.Fit(func, 'QS+')
                    # func.Draw('same')
                    # funcs[name+'func2D'] = func
                    f_.cd()
                    graph.Write()
                    canv.SaveAs(os.path.join(self.outdir,canv_name+'.pdf'))
        f_.Close()

    def PlotVsPt(self):
        colors = {'avg': (rt.kBlack,    rt.kFullCircle),
                  '0':   (rt.kRed+1,    rt.kFullTriangleDown),
                  '10':  (rt.kOrange+1, rt.kFullTriangleUp),
                  '20':  (rt.kOrange,   rt.kFullSquare),
                  '30':  (rt.kAzure+2,  rt.kFullDiamond),
                  '40':  (rt.kViolet-3, rt.kFullStar),
                  '50':  (rt.kCyan+2,   rt.kFullCross),
                  }
        funcs = {}
        for year in self.years['UL']:
            for type in ['ak4pfchs', 'ak4puppi']:
                for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                    eta_min, eta_max = (str(JERC_Constants.GetEtaBinEdgeMin(etaRef)),str(JERC_Constants.GetEtaBinEdgeMax(etaRef)))
                    if eta_min =='0.0': eta_min='0'
                    name = year+type+eta_min+eta_max
                    canv_name = 'JER_JetEta'+year+'_'+type+'_'+JERC_Constants.BinToString(float(eta_min),float(eta_max))
                    tdr.extraText3 = []
                    tdr.extraText3.append('AK4, PF+CHS')
                    tdr.extraText3.append(eta_min+' < |#eta| < '+eta_max)
                    # canv = tdrCanvas(canv_name, 8, 4500, 0, 0.5, 'p_{T}', 'JER')
                    canv = tdrDiCanvas(canv_name, 8, 3000, 0.00001, 0.60, 0.8, 1.15, 'p_{T}', 'JER', 'Func/hist')
                    canv.cd(1).SetLogx(True)
                    canv.cd(2).SetLogx(True)
                    canv.cd(1)
                    leg = tdrLeg(0.70,0.50,0.89,0.89, textSize=0.035)
                    graphs = {}
                    jer_avg = self.inputs['Mu_avg'][name]
                    color, marker = colors['avg']
                    tdrDraw(jer_avg, 'P', marker=marker, mcolor=color)
                    func = ROOT.TF1(name+'func1D', NSC_1D, 10, 3000, 6)
                    # func.FixParameter(0, 1.47/0.87)
                    # func.FixParameter(1, 0.526 1.053)
                    func.FixParameter(0, 1.2560)
                    func.FixParameter(1, 0.4959)
                    # func.FixParameter(0, 0.0988)
                    # func.FixParameter(1, 0.4686)
                    func.SetParLimits(2, 0, 1)
                    func.SetParLimits(3,-2, 0)
                    func.SetParameter(4, 1e-02)
                    func.SetParLimits(4, 1e-03, 1e-01)
                    func.FixParameter(5,32)
                    fitres = jer_avg.Fit(func, 'RQMS+')
                    ROOT.gStyle.SetOptFit(0)
                    func.SetLineColor(color)
                    func.Draw('same')
                    temp_name = name+'func1D'
                    funcs[temp_name]= func
                    graphs['FitBands'+temp_name], graphs['FitRatioBands'+temp_name] = ComputeHistWithCL(name+'func1D', func, fitres, jer_avg, cl=0.68)
                    graphs['FitBands'+temp_name].SetMarkerSize(0)
                    graphs['FitBands'+temp_name].SetLineWidth(0)
                    tdrDraw(graphs['FitBands'+temp_name], 'e3',  fcolor = color, alpha = 0.35)
                    canv.cd(2)
                    tdrDraw(graphs['FitRatioBands'+temp_name], 'e3', fcolor=color, alpha=0.35)
                    graphs['RatioFitHist'+temp_name] = TGraphRatio(graphs['FitBands'+temp_name], jer_avg)
                    tdrDraw(graphs['RatioFitHist'+temp_name], 'P', marker=marker, mcolor=color)
                    mu_bins = [0,10,20,30,40,50]
                    if 'UL16APV' in year: mu_bins = [0,10,20,30,40]
                    for mu in mu_bins:
                        name_mu = str(mu)+'to'+str(mu+10)
                        color, marker = colors[str(mu)]
                        if 'UL16APVak4pfchs2.8532.96450to60' == name+name_mu: continue
                        if 'UL16APVak4puppi2.8532.96450to60' == name+name_mu: continue
                        if 'UL16APVak4puppi2.9643.13950to60' == name+name_mu: continue
                        print name+name_mu
                        graph = self.inputs['Mu_bin'][name+name_mu]
                        canv.cd(1)
                        tdrDraw(graph, 'P', marker=marker, mcolor=color)
                        func = ROOT.TF1(name+'func1D'+str(mu), NSC_1D, 10, 3000, 6)
                        func.FixParameter(2, funcs[name+'func1D'].GetParameter(2))
                        func.FixParameter(3, funcs[name+'func1D'].GetParameter(3))
                        func.FixParameter(4, funcs[name+'func1D'].GetParameter(4))
                        func.FixParameter(5,mu+5)
                        fitres = graph.Fit(func, 'RQMS+')
                        ROOT.gStyle.SetOptFit(0)
                        func.SetLineColor(color)
                        func.Draw('same')
                        funcs[name+'func1D'+str(mu)] = func
                        temp_name = name+'func1D'+str(mu)
                        graphs['FitBands'+temp_name], graphs['FitRatioBands'+temp_name] = ComputeHistWithCL(name+'func1D'+str(mu), func, fitres, graph, cl=0.68)
                        tdrDraw(graphs['FitBands'+temp_name], 'e3',  fcolor = color, alpha = 0.35)
                        graphs['FitBands'+temp_name].SetMarkerSize(0)
                        graphs['FitBands'+temp_name].SetLineWidth(0)
                        canv.cd(2)
                        tdrDraw(graphs['FitRatioBands'+temp_name], 'e3', fcolor=color, alpha=0.35)
                        graphs['RatioFitHist'+temp_name] = TGraphRatio(graphs['FitBands'+temp_name], graph)
                        tdrDraw(graphs['RatioFitHist'+temp_name], 'P', marker=marker, mcolor=color)
                    print os.path.join(self.outdir,canv_name+'_partfix.pdf')
                    canv.SaveAs(os.path.join(self.outdir,canv_name+'_partfix.pdf'))

    def DoPlots(self):
        # self.Plot2DGraph()
        self.PlotVsPt()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug',         action='store_true', default=False, dest='debug')
    parser.add_argument('--postfix', '-p', action='store',      default='',    dest='postfix')
    args = parser.parse_args()
    debug = args.debug
    postfix = args.postfix

    # Comb = JERCombination()
    # Comb.ExtractInfo()
    # Comb.DoCombination()
    # Comb.PlotVsEta()
