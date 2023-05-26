from Utils import *
from GraphHistUtils import *

tdr.writeExtraText = True
tdr.extraText = 'Work in progress'

class StudyNoiseTerm(InputBase):
    def __init__(self):
        InputBase.__init__(self, inputdir = os.path.join('Path_JERC','Summer20year','JER_noise'))
        self.modes = ['nominal', 'sigrc']
        # self.years['UL'] = ['UL16APV', 'UL17', 'UL18']
        self.LoadFiles()
        self.DoPlots()

    def LoadFiles(self):
        fname='RC_noise'
        self.inputs = {'NoiseVsPu':OrderedDict(), 'NoiseVsEta': OrderedDict()}
        self.funcs = OrderedDict()
        self.FitBands = OrderedDict()
        self.FitRatioBands = OrderedDict()
        self.params = OrderedDict()
        self.graphs = OrderedDict()
        for year in self.years['UL']:
            inputdir = self.inputdir.replace('year', year).replace('nonAPV','')
            f_ = ROOT.TFile(os.path.join(inputdir,fname+'_'+year+'.root'))
            self.params[year] = OrderedDict()
            for type in ['MC', 'Data']:
                self.params[year][type] = OrderedDict()
                for mode in self.modes:
                    name = year+type+mode
                    hname='rc_noiseterm_jer_type_'+mode
                    self.inputs['NoiseVsEta'][name] = f_.Get(hname.replace('type',type))
                    # h2D = f_.Get(hname.replace('type',type).replace('rc_noiseterm', 'rc_noiseterm_vs_rho'))
                    h2D = f_.Get(hname.replace('type',type).replace('rc_noiseterm', 'rc_noiseterm_vs_npu'))
                    h2D_ref = f_.Get(hname.replace('type',type).replace('rc_noiseterm', 'rc_noiseterm_vs_npu').replace(mode,list(set(self.modes)-set([mode]))[0]))
                    pars = { 'eta':[], 'par0':[],'par1':[],'err0':[], 'err1':[]}
                    self.params[year][type][mode] = OrderedDict([('par0',OrderedDict()), ('par1',OrderedDict())])
                    for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                        eta_center = JERC_Constants.GetEtaBinCenter(h2D.GetXaxis().GetBinCenter(eta_bin+1))
                        if etaRef!=eta_center:
                            raise ValueError('This eta binning is not expected: '+str(etaRef)+' -> '+str(eta_center))
                        name_eta = name+GetEtaName(etaRef)
                        hist = h2D.ProjectionY(name_eta, eta_bin+1, eta_bin+1,'e')
                        self.inputs['NoiseVsPu'][name_eta] = hist
                        self.inputs['NoiseVsPu'][name_eta].SetDirectory(0)
                        href = h2D_ref.ProjectionY(name_eta, eta_bin+1, eta_bin+1,'e')
                        href.Add(hist,-1)
                        for bin in range(hist.GetNbinsX()+1):
                            # err = math.fabs(hist.GetBinError(bin)-2*math.fabs(href.GetBinContent(bin)))
                            # hist.SetBinError(bin,err/3)
                            # hist.SetBinError(bin,math.fabs(href.GetBinContent(bin))/2)
                            hist.SetBinContent(bin,math.fabs(hist.GetBinContent(bin))*(1.10 if 'sigrc' in mode else 0.90))
                            hist.SetBinError(bin,math.fabs(hist.GetBinContent(bin))*0.025)
                            if 'UL16' in year and href.GetBinCenter(bin)>40:
                                hist.SetBinContent(bin,0)
                                hist.SetBinError(bin,0)
                        func = rt.TF1(name_eta, 'TMath::Sqrt([0]*[0]+[1]*[1]*x)', 0, 45 if 'UL16' in year else 65)
                        # func.SetParLimits(0,0.5,10)
                        self.funcs[name_eta] = func
                        fitRes = hist.Fit(func,'RQMS', '', 10, 40 if 'UL16' in year else 60)
                        self.FitBands[name_eta], self.FitRatioBands[name_eta] = ComputeHistWithCL(name_eta, func, fitRes, hist, cl=0.68, isGraph=False)
                        hist.GetListOfFunctions().RemoveLast()
                        ROOT.gStyle.SetOptFit(0)
                        par0,err0,par1,err1 = func.GetParameter(0), func.GetParError(0), func.GetParameter(1), func.GetParError(1)
                        par0,err0,par1,err1 = (round(math.fabs(par0),4), round(math.fabs(err0),4), round(math.fabs(par1),4), round(math.fabs(err1),4))
                        self.params[year][type][mode]['par0'][GetEtaName(etaRef)] = (par0,err0)
                        self.params[year][type][mode]['par1'][GetEtaName(etaRef)] = (par1,err1)
                        pars['par0'].append(par0)
                        pars['par1'].append(par1)
                        pars['err0'].append(err0)
                        pars['err1'].append(err1)
                        pars['eta'].append(etaRef)
                    self.graphs[name+'par0'] = ROOT.TGraphErrors(len(pars['eta']), array('d',pars['eta']), array('d',pars['par0']), array('d',[0]*len(pars['eta'])), array('d',pars['err0']))
                    self.graphs[name+'par1'] = ROOT.TGraphErrors(len(pars['eta']), array('d',pars['eta']), array('d',pars['par1']), array('d',[0]*len(pars['eta'])), array('d',pars['err1']))
            f_.Close()
        with open(os.path.join(self.outdir,'NTerm.json'), 'w') as outfile:
            json.dump(self.params, outfile, indent=4)

    def DoPlots(self):
        for eta_bin, etaRef in enumerate(self.etaBinsCommon):
            self.PlotSingle(eta_bin,etaRef)
        self.PlotParametersVsEta()

    def PlotSingle(self, eta_bin, etaRef):
        lines = {}
        l_data = rt.TGraph()
        l_data.SetLineStyle(rt.kDashed)
        l_data.SetMarkerStyle(rt.kFullTriangleUp)
        l_MC = rt.TGraph()
        l_MC.SetLineStyle(rt.kSolid)
        l_MC.SetMarkerStyle(rt.kFullCircle)
        eta_min, eta_max = GetEtaMinMax(etaRef)
        name_eta = GetEtaName(etaRef)
        tdr.extraText3 = []
        tdr.extraText3.append('AK4, PF+CHS')
        tdr.extraText3.append(str(eta_min)+' < |#eta| < '+str(eta_max))
        canv_name = name_eta
        # canv = tdrDiCanvas(canv_name, 0, 70, 0, 10, 0.5, 1.0, 'N_{PU}', 'Noise term', 'RMS/gaus')
        canv = tdrDiCanvas(canv_name, 0, 70, 0.01, 5.5 if eta_bin <=10 else 9, 0.9, 1.3, 'N_{PU}', 'Noise term', 'Data/MC')
        leg = tdrLeg(0.65,0.40,0.90,0.60, textSize=0.035)
        leg2 = tdrLeg(0.45,0.80,0.65,0.90, textSize=0.035)
        leg2.AddEntry(l_MC,'MC','LP')
        leg2.AddEntry(l_data,'Data','LP')
        for name, h in sorted(self.inputs['NoiseVsPu'].items()):
            if not name_eta in name: continue
            if not self.modes[0] in name: continue
            mstyle = rt.kFullCircle if 'MC' in name else rt.kFullTriangleUp
            lstyle = rt.kSolid if 'MC' in name else rt.kDashed
            year = name[:name.find('MC' if 'MC' in name else 'Data')]
            color = tdr.commonScheme['color'][year]
            canv.cd(1)
            tdrDraw(h, 'P', mstyle, color)
            # self.leg_TypesExtra.AddEntry(info['points'], '', 'P')
            tdrDrawLine(self.funcs[name], color, lstyle, 2)
            if 'MC' in name:
                leg.AddEntry(h, year, 'p')
            eta, N = (ROOT.Double(0),ROOT.Double(0))
            self.inputs['NoiseVsEta'][name.strip(name_eta)].GetPoint(int(eta_bin),eta,N)
            eta_center = JERC_Constants.GetEtaBinCenter(eta)
            if etaRef!=eta_center:
                raise ValueError('This eta binning is not expected: '+str(etaRef)+' -> '+str(eta_center))
            lines[name] = rt.TLine(0, N, 5, N)
            rho = self.funcs[name].GetX(N,0,100)
            lines[name+'inter'] = rt.TLine(rho, 0, rho, 1)
            for lname in ['','inter']:
                lines[name+lname].SetLineWidth(2)
                tdrDrawLine(lines[name+lname], color, lstyle, 2)
            if 'MC' in name: continue
            canv.cd(2)
            self.FitRatioBands[name+'ratio'] = TGraphRatio(self.FitBands[name], self.FitBands[name.replace('Data','MC')])
            tdrDraw(self.FitRatioBands[name+'ratio'], 'lc', marker=mstyle, mcolor=color)
            # self.inputs['NoiseVsPu'][name+'ratio'] = self.inputs['NoiseVsPu'][name.replace(self.modes[0],self.modes[1])].Clone(name+'ratio')
            self.inputs['NoiseVsPu'][name+'ratio'] = h.Clone(name+'ratio')
            self.inputs['NoiseVsPu'][name+'ratio'].Divide(self.inputs['NoiseVsPu'][name+'ratio'],self.inputs['NoiseVsPu'][name.replace('Data','MC')], 1,1,'B')
            tdrDraw(self.inputs['NoiseVsPu'][name+'ratio'], 'P', mstyle, color)


        canv.SaveAs(os.path.join(self.outdir,canv_name+'.pdf'))
        del canv

    def PlotParametersVsEta(self):
        tdr.extraText3 = []
        for par in [0,1]:
            canv_name = 'Parameter_'+str(par)+'VsEta'
            PlotYMin, PlotYMax, PlotYMin2, PlotYMax2 = (0.3, 1.4, 0.7, 1.5) if par ==1 else (0.0001,5,0., 2.0)
            canv = tdrDiCanvas(canv_name, 0, 5, PlotYMin, PlotYMax, PlotYMin2, PlotYMax2, '|eta|', 'Parameter '+str(par), 'Data/MC')
            leg = tdrLeg(0.65,0.65,0.90,0.90, textSize=0.035)
            leg1 = tdrLeg(0.45,0.80,0.65,0.90, textSize=0.035)
            extraobjs = {}
            for year in self.years['UL']:
                color = tdr.commonScheme['color'][year]
                for type in ['MC', 'Data']:
                    lstyle = rt.kSolid if 'MC' in type else rt.kDashed
                    name = year+type+self.modes[0]
                    graph = self.graphs[year+type+self.modes[0]+'par'+str(par)]
                    canv.cd(1)
                    tdrDraw(graph, 'l', lcolor = color, lstyle=lstyle)
                    if 'MC' in type:
                        leg.AddEntry(graph, year, 'l')
                    if year == self.years['UL'][0]:
                        extraobjs[type] = ROOT.TLine()
                        extraobjs[type].SetLineStyle(lstyle)
                        leg1.AddEntry(extraobjs[type], type, 'l')
                    if not 'MC' in type:
                        ref_graph = self.graphs[year+'MC'+self.modes[0]+'par'+str(par)]
                        self.graphs[year+type+self.modes[0]+'par'+str(par)+'ratio']= TGraphRatio(graph,ref_graph)
                        canv.cd(2)
                        tdrDraw(self.graphs[year+type+self.modes[0]+'par'+str(par)+'ratio'], 'l', lcolor = color, lstyle=lstyle)
            canv.SaveAs(os.path.join(self.outdir,canv_name+'.pdf'))
        del canv



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug',         action='store_true', default=False, dest='debug')
    parser.add_argument('--postfix', '-p', action='store',      default='',    dest='postfix')
    args = parser.parse_args()
    debug = args.debug
    postfix = args.postfix

    Noise = StudyNoiseTerm()
    # print Noise.params
    # params = OrderedDict(filter(lambda x: 'MC' in x[0] and 'UL17' in x[0], Noise.params.items()))
    # print params
