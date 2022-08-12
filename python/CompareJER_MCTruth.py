from Utils import *
from GraphHistUtils import *
from Evaluate_MCJER import *


tdr.writeExtraText = True
tdr.extraText = 'Work in progress'
tdr.extraText2 = 'Simulation'

debug = False
postfix = ''



class CompareJER(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.moduleName = self.__class__.__name__
        self.outdir = self.Path_ANALYSIS+'python/'+self.moduleName+'/'
        os.system('mkdir -p '+self.outdir)
        self.PlotXMin, self.PlotXMax = (7., 4000.)
        self.years = ['UL16APV', 'UL16nonAPV', 'UL17', 'UL18']
        self.algos = ['pfchs','puppi']
        self.radii = ['ak4','ak8']
        self.modes = ['dep', 'avg']


    def AddResolutions(self,name,txtfile):
        self.res[name] = rt.JME.JetResolution(txtfile)

    def SetResolutions(self):
        path = os.path.join(os.environ['CMSSW_BASE'],'src/UHH2/JERCProtoLab')
        self.res = OrderedDict()
        for year in self.years:
            for algo in self.algos:
                for radius in self.radii:
                    fpath = 'Summer20'+year.replace('nonAPV','')+'/MC_JER_avg_mu/Summer20'+year+'_PtResolution_'+radius+algo+'l1l2l3.txt'
                    self.AddResolutions(year+algo+radius+'avg', os.path.join(path,fpath))
                    fpath = 'Summer20'+year.replace('nonAPV','')+'/MC_JER/Summer20'+year.replace('nonAPV','')+'_V'+('3' if 'puppi' in algo else '2')+'_PtResolution_'+radius+algo+'l1l2l3.txt'
                    if year == 'UL18' and algo=='pfchs' and radius=='ak4':
                        fpath = fpath.replace('V2','V3')
                        print fpath
                    self.AddResolutions(year+algo+radius+'dep', os.path.join(path,fpath))

    def CreateCanvas(self):
        tdr.extraText3 = []
        tdr.extraText3.append(str(self.eta_min)+' < |#eta| < '+str(self.eta_max))
        # self.canv = tdrCanvas(self.moduleName+str(self.eta), self.PlotXMin, self.PlotXMax, 0.001, 0.65, 'p_{T}^{jet} [GeV]', 'JER')
        self.canv = tdrDiCanvas(self.moduleName+str(self.eta), self.PlotXMin, self.PlotXMax, 0.001, 0.65, 0.9, 1.1, 'p_{T}^{jet} [GeV]', 'JER', 'Ratio')
        self.leg  = tdrLeg(0.63,0.50,0.89,0.89, textSize=0.035)
        self.canv.cd(1).SetLogx(True)
        self.canv.cd(2).SetLogx(True)
        self.lines = {}
        self.lines['RefRatio'] = rt.TLine(self.PlotXMin, 1, self.PlotXMax, 1)
        self.lines['RefRatio'].SetLineWidth(1)
        self.lines['RefRatio'].SetLineStyle(rt.kDashed)
        self.lines['RefRatio'].SetLineColor(rt.kBlack)
        self.lines['RefRatio'].Draw("same")
        for year in self.years:
            self.lines[year] = rt.TLine()
            self.lines[year].SetLineColor(tdr.commonScheme['color'][year])
            self.leg.AddEntry(self.lines[year], tdr.commonScheme['legend'][year], 'l')

    def PlotOverEta(self):
        # for eta in [0,1.3,2.5,3,4.5]:
        for eta in [0]:
            self.eta = self.GetEtaBinCenter(eta)
            self.eta_min = self.GetEtaBinEdgeMin(self.eta)
            self.eta_max = self.GetEtaBinEdgeMax(self.eta)
            self.pt_max = min(self.PlotXMax, 6500/math.cosh(self.eta_min))
            self.Plot()


    def Plot(self):
        self.CreateCanvas()
        formulas = {}
        style = {
            'pfchsak4dep': rt.kSolid,
            'pfchsak4avg': 10,
            'puppiak4dep': 7,
            'puppiak4avg': rt.kDashed,
            'pfchsak8dep': 9,
            'pfchsak8avg': rt.kDotted,
            'puppiak8dep': 8,
            'puppiak8avg': 4,
        }
        for year in self.years:
            for radius in self.radii:
                for algo in self.algos:
                    for mode in self.modes:
                        name = year+algo+radius+mode
                        is_avg = 'avg' in mode
                        is_chs = 'chs' in algo
                        is_ak4 = '4' in radius
                        mu = 30 if '18' in year else (33 if '17' in year else 23)
                        pars_Mu = rt.JME.JetParameters().setJetPt(0).setJetEta(self.eta).setMu(mu)
                        pars_Rho = rt.JME.JetParameters().setJetPt(0).setJetEta(self.eta).setRho(rhoFromMu(mu))
                        if year == 'UL18' and algo=='pfchs' and radius=='ak4' and not is_avg:
                            pars = pars_Mu
                        else:
                            pars = pars_Rho
                        formulas[name] = GetResolutionFormula(self.res[name], pars, name='Resolution'+name, max=self.pt_max)
                        formulas[name].SetLineColor(tdr.commonScheme['color'][year])
                        formulas[name].SetLineStyle(style[algo+radius+mode])
                        self.canv.cd(1)
                        toDraw = True
                        # if '4' in radius: toDraw = False
                        if '8' in radius: toDraw = False
                        # if 'chs' in algo: toDraw = False
                        if 'puppi' in algo: toDraw = False
                        # if is_avg: toDraw = False
                        # if not is_avg: toDraw = False
                        if toDraw:
                            formulas[name].Draw('same')
                        # if is_avg:
                        # if toDraw and is_avg:
                        # if toDraw:
                            r_name = year+algo+radius+'ratio'
                            if year == 'UL18' and algo=='pfchs' and radius=='ak4' and is_avg:
                                formulas[r_name] = GetResolutionRatio(self.res[name], self.res[name.replace('avg','dep')], pars, pars_Mu, name=r_name, max=self.pt_max)
                            else:
                                formulas[r_name] = GetResolutionRatio(self.res[name], self.res[name.replace('avg','dep')], pars, pars, name=r_name, max=self.pt_max)
                            # formulas[r_name] = GetResolutionRatio(self.res[name], self.res[name.replace('pfchs','puppi')], pars, name=r_name, max=self.pt_max)
                            formulas[r_name].SetLineColor(tdr.commonScheme['color'][year])
                            formulas[r_name].SetLineStyle(style[algo+radius+'avg'])
                            self.canv.cd(2)
                            formulas[r_name].Draw('same')
                        if '18' in year:
                            self.lines[name] = rt.TLine()
                            self.lines[name].SetLineStyle(style[algo+radius+mode])
                            if toDraw:
                                self.leg.AddEntry(self.lines[name], radius+' '+algo+' '+mode, 'l')
        extraname = ''
        extraname += '_ak4'
        # extraname += '_ak8'
        # extraname += '_ak4_vs_ak8_'
        # extraname += 'chs_vs_puppi'
        extraname += 'chs'
        # extraname += 'puppi'
        # extraname += '_dep'
        # extraname += '_avg'
        extraname += '_avg_vs_dep'
        print self.outdir+'CompareJER'+FloatToString(self.eta_min)+'to'+FloatToString(self.eta_max)+extraname+'.pdf'
        self.canv.SaveAs(self.outdir+'CompareJER'+FloatToString(self.eta_min)+'to'+FloatToString(self.eta_max)+extraname+'.pdf')



if __name__ == '__main__':
    print 'start'
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug',         action='store_true', default=False, dest='debug')
    parser.add_argument('--postfix', '-p', action='store',      default='',    dest='postfix')
    args = parser.parse_args()
    debug = args.debug
    postfix = args.postfix


    Comp = CompareJER()
    Comp.SetResolutions()
    Comp.PlotOverEta()
