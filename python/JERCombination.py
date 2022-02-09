from Utils import *
from Evaluate_MCJER import Evaluate_MCJER

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText = 'Work in progress'

def isMC(name):
    return 'MC' in name

def Oplus(x,y):
    return math.sqrt(x*x+y*y)

def Ominus(x,y):
    return math.sqrt(x*x-y*y)

def TGraphRatio(graph1,graph2):
    if graph1.GetN()!=graph2.GetN():
        print 'Error in TGraphRatio: Not same points:', graph1.GetN(), graph2.GetN()
    x_vals = []
    y_vals = []
    x_errs = []
    y_errs = []
    for bin in range(graph1.GetN()):
        x1,y1,x2,y2 = (ROOT.Double(0),ROOT.Double(0),ROOT.Double(0),ROOT.Double(0))
        graph1.GetPoint(int(bin),x1,y1)
        graph2.GetPoint(int(bin),x2,y2)
        if x1!=x2:
            print 'Error in TGraphRatio: Not same x', x1,x2
        x_vals.append(x1)
        y_vals.append(y1/y2)
        x_errs.append(graph1.GetErrorX(bin))
        y_errs.append(graph1.GetErrorY(bin))
    return ROOT.TGraphErrors(len(x_vals), array('d',x_vals), array('d',y_vals), array('d',x_errs), array('d',y_errs))

def HistToGraph(hist):
    x_vals = []
    y_vals = []
    x_errs = []
    y_errs = []
    for bin in range(hist.GetNbinsX()+1):
        if hist.GetBinContent(bin)==0:continue
        x_vals.append(hist.GetBinCenter(bin))
        y_vals.append(hist.GetBinContent(bin))
        x_errs.append(hist.GetBinWidth(bin))
        y_errs.append(hist.GetBinError(bin))
    return ROOT.TGraphErrors(len(x_vals), array('d',x_vals), array('d',y_vals), array('d',x_errs), array('d',y_errs))

def HistToGraphUncertainty(hist, uncertainties):
    x_vals = []
    y_vals = []
    x_errs_down = []
    x_errs_up = []
    y_errs_down = []
    y_errs_up = []
    for bin in range(hist.GetNbinsX()+1):
        if hist.GetBinContent(bin)==0:continue
        y_val = hist.GetBinContent(bin)
        x_vals.append(hist.GetBinCenter(bin))
        y_vals.append(y_val)
        x_errs_down.append(hist.GetBinWidth(bin))
        x_errs_up.append(hist.GetBinWidth(bin))
        y_err_down = hist.GetBinError(bin)
        y_err_up = hist.GetBinError(bin)
        for name, unc in uncertainties.items():
            err = unc.GetBinContent(bin)-y_val
            print name, unc.GetBinContent(bin)-y_val
            if err>0: y_err_up = Oplus(y_err_up,err)
            else: y_err_down = Oplus(y_err_down,err)
        y_errs_down.append(y_err_down)
        y_errs_up.append(y_err_up)
    return ROOT.TGraphAsymmErrors(len(x_vals), array('d',x_vals), array('d',y_vals), array('d',x_errs_down), array('d',x_errs_up), array('d',y_errs_down), array('d',y_errs_up))

def ExpandGraph(graph, x_val, y_val, x_err = 0, y_err = 0):
    x_vals = []
    y_vals = []
    x_errs = []
    y_errs = []
    for bin in range(graph.GetN()):
        x,y = (ROOT.Double(0),ROOT.Double(0))
        graph.GetPoint(int(bin),x,y)
        if x>x_val:
            x_vals.append(x_val)
            y_vals.append(y_val)
            x_errs.append(x_err)
            y_errs.append(y_err)
        x_vals.append(x)
        y_vals.append(y)
        x_errs.append(graph.GetErrorX(bin))
        y_errs.append(graph.GetErrorY(bin))
    return ROOT.TGraphErrors(len(x_vals), array('d',x_vals), array('d',y_vals), array('d',x_errs), array('d',y_errs))


def N_Term(x, par):
    return par[0]/x[0]

def S_Term(x, par):
    return par[0]*ROOT.TMath.Power(x[0],par[1]/2)

def C_Term(x, par):
    return par[0]

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
    isMC = par[9]

    # N_term = (N*ROOT.TMath.Abs(N)-N_RC_MC*N_RC_MC)/(x[0]*x[0])
    # S_term = S*S*ROOT.TMath.Power(x[0],d)
    # C_term = C*C - C_IC_MC*C_IC_MC
    # Data = ROOT.TMath.Sqrt(N_RC_Data*N_RC_Data/(x[0]*x[0]) + k*k*( N_term + S_term + C_term) + C_IC_Data*C_IC_Data)
    # MC   = ROOT.TMath.Sqrt(N_RC_MC  *N_RC_MC  /(x[0]*x[0]) +     ( N_term + S_term + C_term) + C_IC_MC*C_IC_MC)

    N_term = (N*ROOT.TMath.Abs(N))/(x[0]*x[0])
    S_term = S*S*ROOT.TMath.Power(x[0],d)
    C_term = C*C
    Data = ROOT.TMath.Sqrt(N_RC_Data*N_RC_Data/(x[0]*x[0]) + k*k*( N_term + S_term + C_term) + C_IC_Data*C_IC_Data)
    MC   = ROOT.TMath.Sqrt(N_RC_MC  *N_RC_MC  /(x[0]*x[0]) +     ( N_term + S_term + C_term) + C_IC_MC*C_IC_MC)

    if isMC == int(True): return MC
    elif isMC == int(False): return Data
    else: return Data/MC


def NSC_kfactor(x, par):
    N_term = (par[0]*ROOT.TMath.Abs(par[0]))/(x[0]*x[0])
    S_term = par[1]*par[1]*ROOT.TMath.Power(x[0],par[3])
    C_term = par[2]*par[2]

    N_RC_MC   = par[4]*par[4]/(x[0]*x[0])
    N_RC_Data = par[5]*par[5]/(x[0]*x[0])
    C_IC_MC   = par[6]*par[6]
    C_IC_Data = par[7]*par[7]

    k_N = par[8]*par[8]
    k_S = par[9]*par[9]
    k_C = par[10]*par[10]

    isMC = par[11]

    Data = ROOT.TMath.Sqrt(N_RC_Data + k_N*N_term + k_S*S_term + k_C*C_term + C_IC_Data)
    MC   = ROOT.TMath.Sqrt(N_RC_MC   +     N_term +     S_term +     C_term + C_IC_MC)

    if isMC == int(True): return MC
    elif isMC == int(False): return Data
    else: return Data/MC

class JERCombiner(Constants):
    def __init__(self, eta):
        Constants.__init__(self)
        self.eta = self.GetEtaBinCenter(eta)
        self.eta_min = self.GetEtaBinEdgeMin(self.eta)
        self.eta_max = self.GetEtaBinEdgeMax(self.eta)
        self.eta_width = self.GetEtaBinWidth(self.eta)
        self.datas = ['Data', 'MC', 'ratio']
        self.dijet = OrderedDict((data, None) for data in self.datas)
        self.JER = OrderedDict((data, None) for data in self.datas)
        self.parameters = OrderedDict((x, None) for x in ['N','S', 'C', 'd', 'k', 'N_{RC,MC}', 'N_{RC,Data}', 'C_{IC,MC}', 'C_{IC,Data}', 'isMC'])
        self.epsilon = 0.1
        self.fit_min, self.fit_max = (8.-self.epsilon, 3500.+self.epsilon)
        self.Npars = len(self.parameters)
        self.NSCs = OrderedDict((data, ROOT.TF1(str(self.eta)+data, NSC_Modified, self.fit_min, self.fit_max, self.Npars)) for data in self.datas)
        for data in self.datas:
            for mode in ['NSC', 'Noise','RC Noise', 'Stochastic','Constant','2D Constant']:
                name = mode+' '+data
                func = N_Term if 'Noise' in mode else (C_Term if 'Constant' in mode else (S_Term if 'Stochastic' in mode else NSC_Modified))
                npar = self.Npars if 'NSC' in mode else 2 if 'Stochastic' in mode else 1
                self.NSCs[name] = ROOT.TF1(str(self.eta)+name,  func, self.fit_min, self.fit_max, npar)

        self.colors = {
            'Data':        ROOT.kBlack,
            'MC':          ROOT.kBlack,
            'dijet':       ROOT.kBlack,
            'NSC':         ROOT.kViolet,
            'Noise':       ROOT.kRed+1,
            'RC Noise':    ROOT.kOrange+1,
            'Stochastic':  ROOT.kGreen+2,
            'Constant':    ROOT.kAzure+7,
            '2D Constant': ROOT.kViolet-4,
            }

    def GetParameter(self,name):
        return self.parameters[name]

    def SetParameter(self,name,var):
        self.parameters[name] = var

    def GetNSC(self, name):
        return self.NSCs[name]

    def CreateCanvas(self):
        TDR.extraText3 = []
        TDR.extraText3.append('AK4, PF+CHS')
        TDR.extraText3.append(str(self.eta_min)+' < |#eta| < '+str(self.eta_max))
        PlotXMin = self.fit_min-self.epsilon
        PlotXMax = self.fit_max+100
        PlotYMin = 0.
        # PlotYMax = 0.30
        PlotYMax = 0.65
        self.canv = tdrDiCanvas(self.__class__.__name__+str(self.eta), PlotXMin, PlotXMax, PlotYMin, PlotYMax, 0.9, 1.3, 'p_{T}^{jet} [GeV]', 'JER', 'Data / MC')
        self.leg = tdrLeg(0.68,0.50,0.89,0.90, textSize=0.04)
        self.leg2 = tdrLeg(0.38,0.50,0.60,0.90, textSize=0.04)
        self.canv.cd(1).SetLogx(True)
        self.canv.cd(2).SetLogx(True)
        self.canv.cd(2)
        self.lines = {}
        self.lines['RefRatio'] = rt.TLine(PlotXMin, 1, PlotXMax, 1)
        self.lines['RefRatio'].SetLineWidth(1)
        self.lines['RefRatio'].SetLineStyle(rt.kDotted)
        self.lines['RefRatio'].SetLineColor(rt.kBlack)
        self.lines['RefRatio'].Draw("same")
        RemoveRootLabels()

    def CreateFinalGraph(self):
        for data in ['Data', 'MC']:
            JER = math.sqrt(self.GetNSC('RC Noise MC').Eval(self.fit_min+1)**2+self.GetNSC('Stochastic MC').Eval(self.fit_min+1)**2)
            JER = 0.53
            pt = self.fit_min
            JER = self.GetNSC('NSC '+data).Eval(pt)
            # print JER
            JER = Ominus(JER,self.GetNSC('RC Noise MC').Eval(pt))
            JER = Oplus(JER,self.GetNSC('RC Noise '+data).Eval(pt))
            # print self.GetParameter('S')
            self.JER['RC Noise '+data] = ROOT.TGraphErrors(len([1]), array('d',[pt]), array('d',[JER]), array('d',[0.1]), array('d',[0.001]))
            self.JER[data] = self.JER['dijet '+data]
            self.JER[data] = ExpandGraph(self.JER[data], pt, JER, 0.1, 0.001)
            pt = self.fit_max
            JER = self.GetNSC('NSC '+data).Eval(pt)
            JER = Ominus(JER,self.GetNSC('2D Constant MC').Eval(pt))
            JER = Oplus(JER,self.GetNSC('2D Constant '+data).Eval(pt))
            self.JER['Constant '+data] = ROOT.TGraphErrors(len([1]), array('d',[pt]), array('d',[JER]), array('d',[0.1]), array('d',[0.001]))
            self.JER[data] = ExpandGraph(self.JER[data], pt, JER, 0.1, 0.001)

        self.JER['ratio'] = TGraphRatio(self.JER['Data'],self.JER['MC'])
        self.JER['RC Noise ratio'] = TGraphRatio(self.JER['RC Noise Data'],self.JER['RC Noise MC'])
        self.JER['Constant ratio'] = TGraphRatio(self.JER['Constant Data'],self.JER['Constant MC'])


    def FitFunctions(self):
        for data in ['MC', 'Data']:
            if isMC(data):
                self.GetNSC(data).SetParameter(0, self.GetParameter('N_{RC,MC}'))
                self.GetNSC(data).SetParameter(1, self.GetParameter('S'))
                self.GetNSC(data).SetParameter(2, self.GetParameter('C'))
                self.GetNSC(data).SetParameter(3, self.GetParameter('d'))
                # self.GetNSC(data).SetParLimits(0, self.GetParameter('N_{RC,MC}'), 20.0)
                self.GetNSC(data).SetParLimits(0, -20, 20.0)
                self.GetNSC(data).SetParLimits(1,  +0.000,  2.0)
                self.GetNSC(data).SetParLimits(2,  0,  0.1)
                # self.GetNSC(data).SetParLimits(2,  self.GetNSC('NSC MC').Eval(self.fit_max)-self.GetParameter('C_{IC,MC}'),  0.1)
                self.GetNSC(data).SetParLimits(3,  -3.000,  0.0)
                self.GetNSC(data).FixParameter(4, 1.)
                self.GetNSC(data).FixParameter(3, -1.)
            else:
                self.GetNSC(data).FixParameter(0, self.GetNSC('MC').GetParameter(0))
                self.GetNSC(data).FixParameter(1, self.GetNSC('MC').GetParameter(1))
                self.GetNSC(data).FixParameter(2, self.GetNSC('MC').GetParameter(2))
                self.GetNSC(data).FixParameter(3, self.GetNSC('MC').GetParameter(3))
                self.GetNSC(data).SetParameter(4, +1.10)
                self.GetNSC(data).SetParLimits(4, +0.900,  2.0)
            self.GetNSC(data).FixParameter(5, self.GetParameter('N_{RC,MC}'))
            self.GetNSC(data).FixParameter(6, self.GetParameter('N_{RC,'+data+'}'))
            self.GetNSC(data).FixParameter(7, self.GetParameter('C_{IC,MC}'))
            self.GetNSC(data).FixParameter(8, self.GetParameter('C_{IC,'+data+'}'))
            self.GetNSC(data).FixParameter(9, int(isMC(data)))
            self.JER[data].Fit(self.GetNSC(data), 'RMQ+')
            self.JER[data].GetListOfFunctions().Remove(self.JER[data].GetListOfFunctions().FindObject(self.GetNSC(data).GetName()))


    def PlotCombination(self, outdir):

        for data in ['Data', 'MC']:
            self.SetParameter('isMC', int(isMC(data)))
            for i, par in enumerate(self.parameters.keys()):
                self.GetNSC('NSC '+data).SetParameter(i, self.GetParameter(par))
            self.GetNSC('RC Noise '+data).SetParameter(0, self.GetParameter('N_{RC,'+data+'}'))
            self.GetNSC('2D Constant '+data).SetParameter(0, self.GetParameter('C_{IC,'+data+'}'))

        self.CreateFinalGraph()
        self.FitFunctions()

        self.GetNSC('Noise MC').SetParameter(0, self.GetNSC('MC').GetParameter(0))
        self.GetNSC('Stochastic MC').SetParameter(0, self.GetNSC('MC').GetParameter(1))
        self.GetNSC('Stochastic MC').SetParameter(1, self.GetNSC('MC').GetParameter(3))
        self.GetNSC('Constant MC').SetParameter(0, self.GetNSC('MC').GetParameter(2))

        self.CreateCanvas()
        self.canv.cd(1)
        self.leg.AddEntry(self.JER['dijet Data'], 'dijet Data', 'P')
        self.leg.AddEntry(self.JER['dijet MC'], 'dijet MC', 'P')
        self.leg2.AddEntry(ROOT.TObject(), 'N^{0} = '+str(round(self.GetNSC('MC').GetParameter(0),2)), 'P')
        self.leg2.AddEntry(ROOT.TObject(), 'N^{0} #oplus N^{RC} = '+str(round(Oplus(self.GetNSC('MC').GetParameter(0),self.GetNSC('MC').GetParameter(5)),2)), 'P')
        self.leg2.AddEntry(ROOT.TObject(), 'S = '+str(round(self.GetNSC('MC').GetParameter(1),2)), 'P')
        self.leg2.AddEntry(ROOT.TObject(), 'C^{0} = '+str(round(self.GetNSC('MC').GetParameter(2),2)), 'P')
        self.leg2.AddEntry(ROOT.TObject(), 'C^{0} #oplus C^{IC} = '+str(round(Oplus(self.GetNSC('MC').GetParameter(2),self.GetNSC('MC').GetParameter(7)),2)), 'P')
        self.leg2.AddEntry(ROOT.TObject(), 'd = '+str(round(self.GetNSC('MC').GetParameter(3),2)), 'P')
        self.leg2.AddEntry(ROOT.TObject(), 'k = '+str(round(self.GetNSC('Data').GetParameter(4),2)), 'P')
        for data in ['Data', 'MC']:
            tdrDraw(self.JER[data], 'P', ROOT.kOpenCircle if isMC(data) else ROOT.kFullCircle, self.colors[data])
            self.GetNSC(data).SetLineColor(self.colors[data])
            self.GetNSC(data).SetLineStyle(ROOT.kDashed if isMC(data) else ROOT.kSolid)
            self.GetNSC(data).Draw('same')
            for mode in ['RC Noise', 'dijet', 'Constant']:
                tdrDraw(self.JER[mode+' '+data], 'P', ROOT.kOpenCircle if isMC(data) else ROOT.kFullCircle, self.colors[mode])
            for mode in ['RC Noise', '2D Constant', 'Noise', 'Stochastic', 'Constant']:
                name = mode+' '+data
                if not ' ' in mode and not isMC(data): continue
                self.GetNSC(name).SetLineColor(self.colors[mode])
                self.GetNSC(name).SetLineStyle(ROOT.kDashed if isMC(data) else ROOT.kSolid)
                self.GetNSC(name).SetLineWidth(2)
                self.GetNSC(name).Draw('same')
                if isMC(data): self.leg.AddEntry(self.GetNSC(name), mode, 'l')

        self.canv.cd(2)
        for par in range(self.Npars):
            self.GetNSC('ratio').FixParameter(par, self.GetNSC('Data').GetParameter(par))
        self.GetNSC('ratio').FixParameter(9, 2)
        self.GetNSC('ratio').SetLineColor(ROOT.kRed+1)
        self.GetNSC('ratio').Draw('same')
        tdrDraw(self.JER['ratio'], 'P', ROOT.kFullCircle, self.colors['Data'])
        tdrDraw(self.JER['RC Noise ratio'], 'P', ROOT.kFullCircle, self.colors['RC Noise'])
        tdrDraw(self.JER['Constant ratio'], 'P', ROOT.kFullCircle, self.colors['Constant'])
        self.canv.SaveAs(outdir+'Combination'+FloatToString(self.eta_min)+'to'+FloatToString(self.eta_max)+'.pdf')









class JERCombination(VariablesBase):
    def __init__(self, year = 'UL18'):
        VariablesBase.__init__(self)
        self.year = year
        TDR.cms_lumi_TeV = TDR.commonScheme['legend'][self.year]+' Legacy, '+commonScheme['lumi'][self.year]+' fb^{-1}'
        self.moduleName = self.__class__.__name__
        self.inpdir = self.Path_ANALYSIS+'StorageArea/'+self.moduleName+'/'
        self.outdir = self.Path_ANALYSIS+'python/'+self.moduleName+'/'
        os.system('mkdir -p '+self.outdir)
        self.JERCombiners = OrderedDict((eta, JERCombiner(eta)) for eta in self.etaBinsCommon)
        # self.JER = Evaluate_MCJER(GetJERfile(self.JERversions[self.year]))
        os.system('mkdir -p '+self.outdir)
        self.variations = ['PUup','PUdown','JECup','JECdown','alpha','JERnominal','gaustails95']

    def ExtractInfo(self):
        self.functions = {}
        self.ExtractRC()
        self.ExtractCterm()
        # self.ExtractDijet(fname='dijet_andrea')
        # self.ExtractDijet(fname='dijet_alex')
        self.ExtractDijet()

    def ExtractRC(self, fname='RC'):
        f_ = ROOT.TFile(self.inpdir+fname+'.root')
        for data in ['MC', 'Data']:
            name = 'N_{RC,'+('MC' if isMC(data) else 'Data')+'}'
            graph = f_.Get(data+'/RMS')
            for n, etaRef in enumerate(self.etaBinsCommon):
                eta, N = (ROOT.Double(0),ROOT.Double(0))
                graph.GetPoint(int(n),eta,N)
                if round(eta,4) != round(etaRef,4): raise Exception('Something is not ok. Fix me!'+str(eta)+'!='+str(etaRef))
                self.JERCombiners[etaRef].SetParameter(name, N)
        f_.Close()

    def ExtractCterm(self, fname='Cterm'):
        f_ = ROOT.TFile(self.inpdir+fname+'.root')
        for data in ['MC', 'Data']:
            name = 'C_{IC,'+data+'}'
            h_ = f_.Get('jerc_rms_'+data.lower())
            for n, etaRef in enumerate(self.etaBinsCommon):
                if round(h_.GetBinCenter(n+1),4) != round(etaRef,4): raise Exception('Something is not ok. Fix me!'+str(h_.GetBinCenter(n+1))+'!='+str(etaRef))
                self.JERCombiners[etaRef].SetParameter(name, h_.GetBinContent(n+1)/100)
        f_.Close()

    def ExtractDijet(self, fname='dijet'):
        f_ = ROOT.TFile(self.inpdir+fname+'.root')
        # print f_.ls()
        for n, etaRef in enumerate(self.etaBinsCommon):
            dijetbin = 0
            if n<=2    : dijetbin = 1
            elif n<=5  : dijetbin = n-1
            elif n==6  : dijetbin = 4
            elif n<=9  : dijetbin = n-2
            elif n==10 : dijetbin = 7
            elif n<=15 : dijetbin = n-3
            else: dijetbin = 13
            name = 'dijet'+FloatToString(self.etaBinsEdges[n])+'to'+FloatToString(self.etaBinsEdges[n+1])
            for data in ['MC', 'Data']:
                if 'andrea' in fname or 'alex' in fname:
                    self.JERCombiners[etaRef].JER['dijet '+data] = HistToGraph(f_.Get((data if isMC(data) else 'data')+'_JER_correlated_FE'+str(dijetbin)))
                    func = f_.Get((data if isMC(data) else 'data')+'_JER_correlated_FE'+str(dijetbin)).GetListOfFunctions().FindObject('mcFIT' if isMC(data) else 'dtFIT')
                    # print f_.Get((data if isMC(data) else 'data')+'_JER_correlated_FE'+str(dijetbin)).GetListOfFunctions().ls()
                else:
                    name = (data if isMC(data) else 'data')+'_JER_standard_FE'+str(dijetbin)
                    self.JERCombiners[etaRef].JER['dijet '+data] = HistToGraph(f_.Get(name))
                    # variations = OrderedDict((var, f_.Get(name.replace('standard',var))) for var in self.variations)
                    # self.JERCombiners[etaRef].JER['dijet '+data] = HistToGraphUncertainty(f_.Get(name), variations)
                    func = f_.Get(name).GetListOfFunctions().FindObject('mcFIT' if isMC(data) else 'dtFIT')
                    # func = f_.Get(name).GetListOfFunctions().FindObject('mcFITpow' if isMC(data) else 'dtFIT')
                    # print f_.Get(name).GetListOfFunctions().ls()
                if isMC(data):
                    self.JERCombiners[etaRef].SetParameter('N', func.GetParameter(0))
                    self.JERCombiners[etaRef].SetParameter('S', func.GetParameter(1))
                    self.JERCombiners[etaRef].SetParameter('C', func.GetParameter(2))
                    if 'andrea' in fname or 'alex' in fname:
                        self.JERCombiners[etaRef].SetParameter('d', -1)
                    else:
                        self.JERCombiners[etaRef].SetParameter('d', -1)
                        # self.JERCombiners[etaRef].SetParameter('d', func.GetParameter(3))
                else:
                    self.JERCombiners[etaRef].SetParameter('k', 1)
                    # self.JERCombiners[etaRef].SetParameter('k', func.GetParameter(0))
                # self.JERCombiners[etaRef].SetParameter('isMC', 1)
                    # self.JERCombiners[etaRef].SetParameter('d', -func.GetParameter(3))
                # self.JERCombiners[etaRef].dijet[data].GetListOfFunctions().Remove(func)
            # prettydic(self.JERCombiners[etaRef].parameters)
        f_.Close()


    def DoCombination(self):
        for n, etaRef in enumerate(self.etaBinsCommon):
            self.JERCombiners[etaRef].PlotCombination(self.outdir)




if __name__ == '__main__':

    Comb = JERCombination()
    Comb.ExtractInfo()
    Comb.DoCombination()
