from Utils import *
from GraphHistUtils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText = 'Work in progress'

debug = False
postfix = ''

def isMC(name):    return 'MC' in name
def isData(name):  return 'Data' in name
def isRatio(name): return 'ratio' in name

def NSC_kfactor(x, par):
    N_term = (par[0]*ROOT.TMath.Abs(par[0]))/(x[0]*x[0])
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


class JERCombiner(Constants):
    def __init__(self, eta, outdir):
        Constants.__init__(self)
        self.outdir = outdir
        self.eta = self.GetEtaBinCenter(eta)
        self.eta_min = self.GetEtaBinEdgeMin(self.eta)
        self.eta_max = self.GetEtaBinEdgeMax(self.eta)
        self.FitFunc = NSC_kfactor
        self.PlotXMin, self.PlotXMax = (7., 4000.)
        self.fit_min, self.fit_max = (8., 3500.)
        fit_range = (7., 3550.)

        self.PlottingStyle = {
            'Parameters Info': True,
            'Prefit': True,
            'Draw NSC Ratio': True,
        }
        # self.doZjet = self.eta<self.GetEtaBinEdgeMax(1.3)
        # self.doZjet = True
        # self.doZjet_Sonly = False
        self.DefaultRange = (None,None)
        self.Parameters = OrderedDict([
            ('N',           {'pos': 0,  'value': None,  'range': (-20., 20.)}),
            ('S',           {'pos': 1,  'value': None,  'range': (+0.0, 2.0)}),
            ('C',           {'pos': 2,  'value': None,  'range': (+0.0, 0.1)}),
            ('d',           {'pos': 3,  'value': -1,    'range': self.DefaultRange}), #Hypothesis
            # ('d',           {'pos': 3,  'value': -1,    'range': (-3,0)}), #Hypothesis
            ('N_{RC,MC}',   {'pos': 4,  'value': None,  'range': self.DefaultRange}), #Take from inputs
            ('N_{RC,Data}', {'pos': 5,  'value': None,  'range': self.DefaultRange}), #Take from inputs
            ('C_{IC,MC}',   {'pos': 6,  'value': None,  'range': self.DefaultRange}), #Take from inputs
            ('C_{IC,Data}', {'pos': 7,  'value': None,  'range': self.DefaultRange}), #Take from inputs
            ('kN',          {'pos': 8,  'value': +1.10, 'range': (+0.0, 5.0)}),
            ('kS',          {'pos': 9,  'value': +1.10, 'range': (+0.0, 5.0)}),
            ('kC',          {'pos': 10, 'value': +1.10, 'range': (+0.0, 5.0)}),
            ('isMC',        {'pos': 11, 'value': None,  'range': self.DefaultRange}), #Fixed given function
        ])

        name_N_PU_0 = 'N_{'+ScaleLeg('PU=0')+'}'
        name_N_PU   = 'N_{'+ScaleLeg('PU!=0')+'}^{'+ScaleLeg('RC')+'}'
        name_C_0    = 'C_{'+ScaleLeg('0')+'}'
        name_C_IC   = 'C_{'+ScaleLeg('IC')+'}^{'+ScaleLeg('2D')+'}'

        self.Contributes = OrderedDict([
            ('Pure Noise',    {'name': 'Noise '+ScaleLeg('(PU=0)'),  'color': ROOT.kRed+1,     'symbol': name_N_PU_0,  'par': 'N',         'parameters': ['N', 'kN']}),
            ('RC Noise',      {'name': 'Noise '+ScaleLeg('(PU!=0)'), 'color': ROOT.kOrange+1,  'symbol': name_N_PU,    'par': 'N_{RC,MC}', 'parameters': ['N_{RC,MC}', 'N_{RC,Data}']}),
            ('Stochastic',    {'name': 'Stochastic',                 'color': ROOT.kGreen+2,   'symbol': 'S',          'par': 'S',         'parameters': ['S','d', 'kS']}),
            ('2D Constant',   {'name': 'Cal. intercal.',             'color': ROOT.kMagenta+2, 'symbol': name_C_IC,    'par': 'C_{IC,MC}', 'parameters': ['C_{IC,MC}', 'C_{IC,Data}']}),
            ('Pure Constant', {'name': 'Constant',                   'color': ROOT.kAzure+7,   'symbol': name_C_0,     'par': 'C',         'parameters': ['C', 'kC']}),
        ])

        self.Constrained = OrderedDict([
            ('RC Noise',    {'pt': self.fit_min, 'pt_err': 0.5, 'err': 0.01,  'color ratio': ROOT.kOrange+1,  'mstyle': ROOT.kFullCircle, 'mstyle MC': ROOT.kOpenCircle}),
            ('2D Constant', {'pt': self.fit_max, 'pt_err': 250, 'err': 0.002, 'color ratio': ROOT.kMagenta+2, 'mstyle': ROOT.kFullCircle, 'mstyle MC': ROOT.kOpenCircle}),
            ])
        for c_ in self.Constrained: self.Constrained[c_].update(self.Contributes[c_])

        self.FittingFunctions = OrderedDict([
            ('Total',   {'name': 'Total',   'color': ROOT.kBlack,  'parameters': self.Parameters.keys()}),
            ('Partial', {'name': 'Partial', 'color': ROOT.kCyan+2, 'parameters': self.Parameters.keys()}),
        ])
        self.Datasets = OrderedDict([
            ('dijet FE', {'name': 'dijet balance', 'color': ROOT.kBlack, 'color ratio': ROOT.kGray+1, 'mstyle': ROOT.kFullCircle,     'mstyle MC': ROOT.kOpenCircle}),
            ('zjet MPF', {'name': 'Z+jet MPF',     'color': ROOT.kBlack, 'color ratio': ROOT.kGray+1, 'mstyle': ROOT.kFullTriangleUp, 'mstyle MC': ROOT.kOpenTriangleUp}),
        ])

        self.AdditionalDatasets = OrderedDict([
            # ('zjet MPF S-only', {'name': 'Z+jet MPF S-only', 'color': ROOT.kGreen+2, 'color ratio': ROOT.kGreen+2, 'mstyle': ROOT.kFullSquare, 'mstyle MC': ROOT.kOpenSquare}),
            # ('dijet SM',        {'name': 'dijet SM',         'color': ROOT.kRed+2, 'color ratio': ROOT.kRed+2, 'mstyle': ROOT.kFullStar,   'mstyle MC': ROOT.kOpenStar}),
            # ('zjet balance',    {'name': 'Z+jet balance', 'color': ROOT.kRed+2, 'color ratio': ROOT.kRed+2, 'mstyle': ROOT.kFullSquare, 'mstyle MC': ROOT.kOpenSquare}),
        ])

        self.Types = OrderedDict([
            ('Data',  {'name': 'Data',    'line': ROOT.TLine(), 'points': ROOT.TGraph(), 'mstyle': ROOT.kFullCircle, 'lstyle': ROOT.kSolid,  'mcolor': ROOT.kGray+1, 'lcolor': ROOT.kGray+1}),
            ('MC',    {'name': 'MC',      'line': ROOT.TLine(), 'points': ROOT.TGraph(), 'mstyle': ROOT.kOpenCircle, 'lstyle': ROOT.kDashed, 'mcolor': ROOT.kGray+1, 'lcolor': ROOT.kGray+1}),
            ('ratio', {'name': 'Data/MC', 'line': ROOT.TLine(), 'points': ROOT.TGraph(), 'mstyle': ROOT.kFullCircle, 'lstyle': ROOT.kDotted, 'mcolor': ROOT.kGray+1, 'lcolor': ROOT.kGray+1}),
            ])

        self.Npars = len(self.Parameters)
        self.JER = {}
        self.fitRes = {}
        self.FitBands = {}
        self.FitRatioBands = {}
        self.NSCs  = dict((mode+' '+type, ROOT.TF1(str(self.eta)+mode+' '+type, self.FitFunc, fit_range[0], fit_range[1], self.Npars)) for type in self.Types for mode in self.Contributes.keys()+self.FittingFunctions.keys())

    def GetParIndex(self,name):                   return self.Parameters[name]['pos']
    def GetParValue(self,name):                   return self.Parameters[name]['value']
    def GetParRange(self,name):                   return self.Parameters[name]['range']
    def SetParValue(self,name,var):               self.Parameters[name]['value'] = var
    def SetParRange(self,name,var):               self.Parameters[name]['range'] = var
    def NSC(self, name):                          return self.NSCs[name]
    def GetNSCParam(self, name, index):           return self.NSC(name).GetParameter(index)
    def SetNSCParam(self, name, index,value):     self.NSC(name).SetParameter(index, value)
    def SetNSCParLimits(self, name, index,range): self.NSC(name).SetParLimits(index, range[0], range[1])
    def FixNSCParam(self, name, index,value):     self.NSC(name).FixParameter(index, value)

    def CreateRatioDataset(self, dataset):
        if debug: debugStr('Create ratio for: '+dataset, color=yellow)
        print 'Start', dataset
        self.JER[dataset+' ratio'] = TGraphRatio(self.JER[dataset+' Data'],self.JER[dataset+' MC'])
        print 'END', dataset

    def CreateRatioInputDatasets(self):
        for ds in self.Datasets.keys()+self.AdditionalDatasets.keys():
            self.CreateRatioDataset(ds)

    def CreateCombinedDatasets(self, func_name):
        for type in self.Types:
            for ds in self.Datasets.keys()+self.Constrained.keys():
                if ds in self.Constrained and func_name == 'Partial': continue
                name  = func_name+' '+type
                if debug: debugStr('Combine '+name+' with '+ds)
                input = self.JER[ds+' '+type]
                self.JER[name] = MergeGraphs(self.JER[name], input) if name in self.JER else input

    def CreateConstraintGraphs(self):
        # print self.NSCs.keys()
        for constr, info in self.Constrained.items():
            for type in reversed(self.Types.keys()):
                if isRatio(type): continue
                NSC_name = constr+' '+type
                if debug: debugStr('Create constraint graph for: '+NSC_name)
                for par in self.Parameters:
                    value = self.GetParValue(par)
                    if not par in self.Contributes[constr]['parameters']: value = 0
                    self.FixNSCParam(NSC_name, self.GetParIndex(par), value)
                self.FixNSCParam(NSC_name, self.GetParIndex('isMC'), int(isMC(type)))
                JER = 0.43 if constr == 'RC Noise' else 0.045
                JER = 0.48 if constr == 'RC Noise' else 0.045
                # JER = 0.53 if constr == 'RC Noise' else 0.045
                JER = Oplus(JER,self.NSC(NSC_name).Eval(info['pt']))
                JER = Ominus(JER,self.NSC(constr+' MC').Eval(info['pt']))
                # print JER
                if isData(type):
                    k_factor = math.sqrt(1.05**2-1)*self.NSC('Partial MC').GetParameter(self.GetParIndex('N'))/info['pt']
                    JER = Oplus(JER,k_factor) if k_factor>0 else Ominus(JER,abs(k_factor))
                    # print JER, 'N', k_factor
                    k_factor = math.sqrt(1.13**2-1)*self.NSC('Partial MC').GetParameter(self.GetParIndex('S'))/math.sqrt(info['pt'])
                    JER = Oplus(JER,k_factor) if k_factor>0 else Ominus(JER,k_factor)
                    # print JER, 'S', k_factor
                self.JER[constr+' '+type] = ROOT.TGraphErrors(len([1]), array('d',[info['pt']]), array('d',[JER]), array('d',[info['pt_err']]), array('d',[info['err']]))
            self.CreateRatioDataset(constr)

    def FitFunction(self, func_name):
        for type in reversed(self.Types.keys()):
            if isRatio(type): continue
            name = func_name+' '+type
            for par in self.Parameters:
                range = self.GetParRange(par)
                index = self.GetParIndex(par)
                value = self.GetParValue(par)
                if isMC(type) and par in ['kN', 'kS', 'kC']:
                    range = self.DefaultRange
                    value = 0
                if isData(type) and par in ['N', 'S', 'C', 'd']:
                    range = self.DefaultRange
                    value = self.GetNSCParam(func_name+' MC', index)
                if par == 'isMC':
                    range = self.DefaultRange
                    value = int(isMC(type))
                if range == self.DefaultRange:
                    self.FixNSCParam(name,index, value)
                else:
                    self.SetNSCParam(name, index, value)
                    self.SetNSCParLimits(name, index, range)
            if debug:
                debugStr('Eta: '+str(self.eta)+' Fitting '+name+' with '+str(self.JER[name].GetN())+' points and '+str(self.NSC(name).GetNumberFreeParameters())+' free parameters')
                for par in self.Parameters:
                    index = self.GetParIndex(par)
                    parmin, parmax = (ROOT.Double(0),ROOT.Double(0))
                    self.NSC(name).GetParLimits(index, parmin,parmax)
                    debugStr('  Prefit: '+par+' = '+str(round(self.GetNSCParam(name, index),3))+('' if parmin!=parmax else '  FIXED'), color=yellow)
            self.fitRes[name] = self.JER[name].Fit(self.NSC(name), 'RMQS+')
            self.FitBands[name], self.FitRatioBands[name] = ComputeHistWithCL(name, self.NSC(name), self.fitRes[name], self.JER[name], cl=0.68)
            # self.JER[name].GetListOfFunctions().Remove(self.JER[name].GetListOfFunctions().FindObject(self.NSC(name).GetName()))
            if debug:
                for par in self.Parameters:
                    index = self.GetParIndex(par)
                    parmin, parmax = (ROOT.Double(0),ROOT.Double(0))
                    self.NSC(name).GetParLimits(index, parmin,parmax)
                    debugStr('  Postfit: '+par+' = '+str(round(self.GetNSCParam(name, index),3))+('' if parmin!=parmax else '  FIXED'), color=yellow)
            if isData(type):
                for par in self.Parameters:
                    index = self.GetParIndex(par)
                    value = self.GetNSCParam(func_name+' Data', index)
                    self.FixNSCParam(name.replace(type,'ratio'), index, value)
                self.FixNSCParam(name.replace(type,'ratio'), self.GetParIndex('isMC'), 2)

        self.FitRatioBands[func_name+' ratio'] = TGraphRatio(self.FitBands[func_name+' Data'], self.FitBands[func_name+' MC'])

    def FixNSCContributions(self):
        for type in self.Types:
            if isData(type): continue
            for mode, info in self.Contributes.items():
                name = mode+' '+type
                if debug: debugStr('Creating NSC contributions for: '+name, color=yellow)
                for par in self.Parameters:
                    index = self.GetParIndex(par)
                    value = self.GetNSCParam('Total '+type, index)
                    if not par in info['parameters']: value = 0
                    self.FixNSCParam(name, index, value)
                self.FixNSCParam(name, self.GetParIndex('isMC'), 2 if isRatio(type) else int(isMC(type)))

    def CreateCanvas(self):
        TDR.extraText3 = []
        TDR.extraText3.append('AK4, PF+CHS')
        TDR.extraText3.append(str(self.eta_min)+' < |#eta| < '+str(self.eta_max))
        self.canv = tdrDiCanvas(self.__class__.__name__+str(self.eta), self.PlotXMin, self.PlotXMax, 0.001, 0.65, 0.89, 1.3, 'p_{T}^{jet} [GeV]', 'JER', 'Data / MC')
        self.leg = tdrLeg(0.63,0.48,0.89,0.93, textSize=0.035)
        self.leg_Par = tdrLeg(0.28,0.60,0.60,0.80, textSize=0.032)
        self.leg_Par.SetNColumns(2)
        self.leg_Types = tdrLeg(0.40,0.8,0.65,0.92, textSize=0.035)
        self.leg_TypesExtra = tdrLeg(0.40,0.81,0.65,0.93, textSize=0.035)
        self.canv.cd(1).SetLogx(True)
        self.canv.cd(2).SetLogx(True)
        self.canv.cd(2)
        self.lines = {}
        self.lines['RefRatio'] = rt.TLine(self.PlotXMin, 1, self.PlotXMax, 1)
        self.lines['RefRatio'].SetLineWidth(1)
        self.lines['RefRatio'].SetLineStyle(rt.kDashed)
        self.lines['RefRatio'].SetLineColor(rt.kBlack)
        self.lines['RefRatio'].Draw("same")
        RemoveRootLabels()

    def DrawLegTypes(self):
        self.canv.cd(1)
        for info in self.Types.values():
            info['line'].SetLineColor(info['lcolor'])
            info['line'].SetLineStyle(info['lstyle'])
            info['line'].SetLineWidth(2)
            info['points'].SetMarkerColor(info['mcolor'])
            info['points'].SetMarkerStyle(info['mstyle'])
            self.leg_Types.AddEntry(info['line'], info['name'], 'L')
            self.leg_TypesExtra.AddEntry(info['points'], '', 'P')

    def DrawLegParams(self):
        if not self.PlottingStyle['Parameters Info']: return
        self.canv.cd(1)
        for mode in ['Pure Noise', 'Stochastic', 'RC Noise', 'kN', 'Pure Constant', 'kS', '2D Constant', 'kC']:
            if 'k' in mode:
                par = mode
                name = 'k_{'+ScaleLeg(mode[1])+'}'
            else:
                par = self.Contributes[mode]['par']
                name = self.Contributes[mode]['symbol']
            name += ' = '+str(round(self.GetNSCParam('Total Data', self.GetParIndex(par)),2))
            self.leg_Par.AddEntry(ROOT.TObject(), name, 'P')


    def DrawGraphs(self):
        for type in self.Types:
            if not self.PlottingStyle['Draw NSC Ratio'] and isRatio(type): continue
            if isRatio(type): self.canv.cd(2)
            else: self.canv.cd(1)
            items = self.Datasets.items()+self.AdditionalDatasets.items()+self.Constrained.items()
            for name, info in items:
                if debug: debugStr('Plotting '+name+' '+type, color=green)
                if isData(type) and not name in self.Constrained:
                    self.leg.AddEntry(self.JER[name+' Data'], info['name'], 'P')
                tdrDraw(self.JER[name+' '+type], 'P', info['mstyle'+(' MC' if isMC(type) else '')], info['color'])
            if isRatio(type): continue
            self.canv.cd(2)
            if not self.PlottingStyle['Draw NSC Ratio']: continue
            self.JER['JER ratio '+type] = TGraphFuncRatio(self.JER['Total '+type],self.NSC('Total '+type))
            tdrDraw(self.JER['JER ratio '+type], 'P', self.Types[type]['mstyle'], self.Types[type]['mcolor'])

    def DrawNSCs(self):
        if debug: debugStr('NCSs:'+str(self.NSCs.keys()))
        if debug: debugStr('FitBands:'+str(self.FitBands.keys()))
        if debug: debugStr('FitBands:'+str(self.FitRatioBands.keys()))
        for type in self.Types:
            if isRatio(type): self.canv.cd(2)
            else: self.canv.cd(1)
            for mode, info in self.FittingFunctions.items()+self.Contributes.items():
                name = mode+' '+type
                lstyle = ROOT.kSolid if isRatio(type) and 'Total' == mode else self.Types[type]['lstyle']
                color = info['color']
                if not self.PlottingStyle['Prefit'] and not mode in self.Constrained: continue
                if isData(type):
                    legInfo = info['name']
                    if 'Total' == mode or 'Partial' == mode:
                        legInfo += ' '+str(round(math.sqrt(self.NSC(name.replace(type, 'MC')).GetChisquare()/self.NSC(name.replace(type, 'MC')).GetNDF()),2))
                        if 'Total' in mode:
                            legInfo += ' '+str(round(math.sqrt(self.NSC(name).GetChisquare()/self.NSC(name).GetNDF()),2))
                    self.leg.AddEntry(self.NSC(name), legInfo, 'l')
                self.NSC(name).SetLineColor(color)
                self.NSC(name).SetLineStyle(lstyle)
                self.NSC(name).SetLineWidth( 2)
                # if isData(type) and mode!='Total': continue
                if not self.PlottingStyle['Draw NSC Ratio'] and isRatio(type): continue
                if debug: debugStr('Draw NSCs '+name, color=green)
                if mode in self.FittingFunctions:
                    if mode != 'Total' and type!="MC": continue
                    if isRatio(type):
                        for type2 in self.Types:
                            tdrDraw(self.FitRatioBands[name.replace(type,type2)], "e3",  fcolor = color, alpha = 0.1 if 'MC'==type2 else 0.35)
                            self.FitRatioBands[name.replace(type,type2)].SetMarkerSize(0)
                            self.FitRatioBands[name.replace(type,type2)].SetLineWidth(0)
                    else:
                        tdrDraw(self.FitBands[name], "e3",  fcolor = color, alpha = 0.35)
                        self.FitBands[name].SetMarkerSize(0)
                        self.FitBands[name].SetLineWidth(0)
                if isRatio(type):
                    if 'Pure Noise' in name:
                        self.FixNSCParam(name, self.GetParIndex('N'), math.fabs(self.GetNSCParam(name, self.GetParIndex('N'))))
                self.NSC(name).Draw('same')

    def PlotCombination(self):
        self.CreateCanvas()
        self.DrawLegTypes()
        self.DrawLegParams()
        self.DrawGraphs()
        self.DrawNSCs()
        self.canv.cd(1)
        TDR.fixOverlay()
        self.canv.cd(2)
        TDR.fixOverlay()
        self.canv.SaveAs(self.outdir+'Combination'+FloatToString(self.eta_min)+'to'+FloatToString(self.eta_max)+('' if postfix=='' else '_'+postfix)+'.pdf')


    def DoCombination(self):

        self.CreateRatioInputDatasets()
        self.CreateCombinedDatasets('Partial')
        self.FitFunction('Partial')
        self.CreateConstraintGraphs()
        self.CreateCombinedDatasets('Total')
        self.FitFunction('Total')
        self.FixNSCContributions()
        self.PlotCombination()
        # Do chi2
        # do errorbands








class JERCombination(VariablesBase):
    def __init__(self, year = 'UL18'):
        VariablesBase.__init__(self)
        self.year = year
        TDR.cms_lumi_TeV = TDR.commonScheme['legend'][self.year]+' Legacy, '+commonScheme['lumi'][self.year]+' fb^{-1}'
        self.moduleName = self.__class__.__name__
        self.inpdir = self.Path_ANALYSIS+'StorageArea/'+self.moduleName+'/'
        self.outdir = self.Path_ANALYSIS+'python/'+self.moduleName+'/'
        os.system('mkdir -p '+self.outdir)
        self.JERCombiners = OrderedDict((eta, JERCombiner(eta, self.outdir)) for eta in self.etaBinsCommon)
        os.system('mkdir -p '+self.outdir)
        self.variations = ['PU', 'PLI', 'JEC', 'alpha', 'gaustails95']
        # self.variations = ['PU', 'PLI', 'JEC', 'alpha']

    def ExtractInfo(self):
        self.functions = {}
        self.ExtractRC()
        self.ExtractCterm()
        self.ExtractDijet()
        self.ExtractZjet()
        self.ExtractZjet2D()

    def ExtractRC(self, fname='RC'):
        f_ = ROOT.TFile(self.inpdir+fname+'_'+self.year+'.root')
        for type in ['MC', 'Data']:
            name = 'N_{RC,'+('MC' if isMC(type) else 'Data')+'}'
            graph = f_.Get(type+'/RMS')
            # graph = f_.Get(type+'/SigmaRC')
            for n, etaRef in enumerate(self.etaBinsCommon):
                eta, N = (ROOT.Double(0),ROOT.Double(0))
                graph.GetPoint(int(n),eta,N)
                if round(eta,4) != round(etaRef,4): raise Exception('Something is not ok. Fix me!'+str(eta)+'!='+str(etaRef))
                self.JERCombiners[etaRef].SetParValue(name, N)
        f_.Close()

    def ExtractCterm(self, fname='Cterm'):
        f_ = ROOT.TFile(self.inpdir+fname+'_'+self.year+'.root')
        for type in ['MC', 'Data']:
            name = 'C_{IC,'+type+'}'
            h_ = f_.Get('jerc_rms_'+type.lower())
            for n, etaRef in enumerate(self.etaBinsCommon):
                if round(h_.GetBinCenter(n+1),4) != round(etaRef,4): raise Exception('Something is not ok. Fix me!'+str(h_.GetBinCenter(n+1))+'!='+str(etaRef))
                self.JERCombiners[etaRef].SetParValue(name, h_.GetBinContent(n+1)/100)
        f_.Close()

    def ExtractDijet(self, fname='dijet'):
        f_ = ROOT.TFile(self.inpdir+fname+'_'+self.year+'.root')
        for n, etaRef in enumerate(self.etaBinsCommon):
            name = 'dijet'+FloatToString(self.etaBinsEdges[n])+'to'+FloatToString(self.etaBinsEdges[n+1])
            for type in ['MC', 'Data']:
                for method in ['FE', 'SM']:
                    name = type+'_jer_dijet_'+method+'_'+str(n+1)+'_nominal'
                    variations = {}
                    for var in self.variations:
                        if var == 'JEC' and isData(type): continue
                        if var == 'PU' and isData(type): continue
                        if var in ['alpha', 'gaustails95']:
                            variations[var] = {
                                'up':   f_.Get(name.replace('nominal',var)),
                                'down': f_.Get(name.replace('nominal',var))
                                }
                        else:
                            variations[var] = {
                                'up':   f_.Get(name.replace('nominal',var+'up')),
                                'down': f_.Get(name.replace('nominal',var+'down'))
                                }
                    self.JERCombiners[etaRef].JER['dijet '+method+' '+type] = HistToGraphUncertainty(f_.Get(name), variations)
                self.JERCombiners[etaRef].SetParValue('N', 0.1)
                self.JERCombiners[etaRef].SetParValue('S', 0.9)
                self.JERCombiners[etaRef].SetParValue('C', 0.004)
                self.JERCombiners[etaRef].SetParValue('d', -1)
                self.JERCombiners[etaRef].SetParValue('kN', 1)
                self.JERCombiners[etaRef].SetParValue('kS', 1)
                self.JERCombiners[etaRef].SetParValue('kC', 1)
        f_.Close()

    def ExtractZjet(self, fname='zjet_MPF'):
        f_ = ROOT.TFile(self.inpdir+fname+'_'+self.year+'.root')
        for n, etaRef in enumerate(self.etaBinsCommon):
            for type in ['MC', 'Data']:
                self.JERCombiners[etaRef].JER['zjet MPF '+type] = HistToGraph(f_.Get('zjet_jer_sn_'+type.lower()), min_val=15, max_val=120 )
                self.JERCombiners[etaRef].JER['zjet MPF S-only '+type] = HistToGraph(f_.Get('zjet_jer_s_'+type.lower()), min_val=15, max_val=120, remove_values=[27.5])
        f_.Close()

    def ExtractZjet2D(self, fname='zjet_balance'):
        for type in ['MC', 'Data']:
            f_ = ROOT.TFile(self.inpdir+fname+'_'+type+'_'+self.year+'.root')
            jer = f_.Get('jer')
            for y_bin in range(1, jer.GetNbinsY()+1):
                if jer.GetYaxis().GetBinCenter(y_bin) > self.etaBinsEdges[-1]: continue
                etaRef = self.GetEtaBinCenter(jer.GetYaxis().GetBinCenter(y_bin))
                x_vals = []
                y_vals = []
                x_errs = []
                y_errs = []
                for x_bin in range(1, jer.GetNbinsX()+1):
                    x_val = jer.GetXaxis().GetBinCenter(x_bin)
                    y_val = jer.GetBinContent(x_bin,y_bin)
                    x_vals.append(x_val)
                    y_vals.append(y_val)
                    x_errs.append(jer.GetXaxis().GetBinWidth(x_bin)/2)
                    y_errs.append(jer.GetBinError(x_bin,y_bin))
                self.JERCombiners[etaRef].JER['zjet balance '+type] = ROOT.TGraphErrors(len(x_vals), array('d',x_vals), array('d',y_vals), array('d',x_errs), array('d',y_errs))
            f_.Close()

    def DoCombination(self):
        for n, etaRef in enumerate(self.etaBinsCommon):
            print 'Combine', etaRef
            self.JERCombiners[etaRef].DoCombination()

    def PlotVsEta(self, canv_name='SF_vs_eta'):
        plotshift = 0.2;
        PlotYMin, PlotYMax = (0.85, 1.65)
        TDR.extraText3 = []
        TDR.extraText3.append('AK4, PF+CHS')
        canv = tdrCanvas(canv_name, self.etaBinsEdges_[0]-plotshift, self.etaBinsEdges_[-1]+plotshift, PlotYMin, PlotYMax, "|#eta|", "JER SF");
        leg = tdrLeg(0.70,0.50,0.89,0.89, textSize=0.035)
        lines = {}
        for eta in [1.31, 2.5, 3.0]:
            lines[eta] = rt.TLine(self.GetEtaBinEdgeMin(eta), PlotYMin, self.GetEtaBinEdgeMin(eta), PlotYMax)
            lines[eta].SetLineWidth(1)
            lines[eta].SetLineStyle(rt.kDashed)
            lines[eta].SetLineColor(rt.kBlack)
            lines[eta].Draw("same")

        objects = OrderedDict([
            ('8',             {'name': 'p_{T} = 8 GeV',    'color': ROOT.kRed+1,    'marker': ROOT.kFullTriangleDown}),
            ('80',            {'name': 'p_{T} = 80 GeV',   'color': ROOT.kAzure+2,  'marker': ROOT.kFullCircle}),
            # ('300',           {'name': 'p_{T} = 300 GeV',  'color': ROOT.kBlack,    'marker': ROOT.kFullTriangleDown}),
            ('3500',          {'name': 'p_{T} = 3500 GeV', 'color': ROOT.kViolet-4, 'marker': ROOT.kFullTriangleUp}),
            ('RC Noise',      {'name': 'RC Noise',         'color': ROOT.kOrange+1, 'marker': ROOT.kOpenDiamond}),
            # ('2D Constant',   {'name': 'Cal. intercal.',   'color': ROOT.kViolet-4, 'marker': ROOT.kOpenSquare}),
            ('dijet',         {'name': 'dijet only',       'color': ROOT.kGreen+2,  'marker': ROOT.kFullSquare}),
            ('dijet_official',{'name': 'JER SF V3',        'color': ROOT.kBlack,    'marker': ROOT.kFullSquare}),
            ])

        dijetSF = {'dijet': [], 'dijet_official': []}
        lines = open('/nfs/dust/cms/user/amalara/WorkingArea/UHH2_106X_v2_UL/CMSSW_10_6_28/src/UHH2/JRDatabase/textFiles/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt').readlines()
        for line in lines[1:]:
            line = line.split()
            if float(line[0])<0: continue
            dijetSF['dijet_official'].append({'min': float(line[0]), 'max': float(line[1]), 'SF': float(line[5])})
        lines = open(self.inpdir+'dijet_SF'+'_'+self.year+'.txt').readlines()
        for line in lines[1:]:
            line = line.split()
            if float(line[0])<0: continue
            dijetSF['dijet'].append({'min': float(line[0]), 'max': float(line[1]), 'SF': float(line[5])})
        graphs = {}
        for name, info in objects.items():
            vals = {x:[] for x in ['x_vals','y_vals','x_errs', 'y_errs']}
            if 'dijet' in name:
                for x in dijetSF[name]:
                    vals['x_vals'].append((x['max']+x['min'])/2)
                    vals['x_errs'].append((x['max']-x['min'])/2)
                    vals['y_vals'].append(x['SF'])
                    vals['y_errs'].append(0.01)
            else:
                for n, etaRef in enumerate(self.etaBinsCommon):
                    SF = self.JERCombiners[etaRef].NSC('Total ratio').Eval(float(name)) if not ' ' in name else self.JERCombiners[etaRef].NSC(name+' ratio').Eval(100)
                    vals['x_vals'].append(self.GetEtaBinCenter(etaRef))
                    vals['x_errs'].append(self.GetEtaBinWidth(etaRef))
                    vals['y_vals'].append(SF)
                    vals['y_errs'].append(0.01)
            graphs[name] = ROOT.TGraphErrors(len(vals['x_vals']), array('d',vals['x_vals']), array('d',vals['y_vals']), array('d',vals['x_errs']), array('d',vals['y_errs']))
            tdrDraw(graphs[name], 'P', info['marker'], info['color'])
            leg.AddEntry(graphs[name], info['name'], 'lp')
        canv.SaveAs(self.outdir+canv_name+'_'+self.year+('' if postfix=='' else '_'+postfix)+'.pdf')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug',         action='store_true', default=False, dest='debug')
    parser.add_argument('--postfix', '-p', action='store',      default='',    dest='postfix')
    args = parser.parse_args()
    debug = args.debug
    postfix = args.postfix

    Comb = JERCombination()
    Comb.ExtractInfo()
    Comb.DoCombination()
    Comb.PlotVsEta()
