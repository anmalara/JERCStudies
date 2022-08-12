from Utils import *
from GraphHistUtils import TGraphRatio, ComputeHistWithCL
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gErrorIgnoreLevel = ROOT.kError

tdr.writeExtraText = True
tdr.extraText = 'Work in progress'


def NSC(x, par):
    pt = x[0]
    mu = par[5]
    N_term = (par[0]*par[0]+par[1]*par[1]*mu)/(pt*pt)
    S_term = par[2]*par[2]*ROOT.TMath.Power(pt,par[3])
    C_term = par[4]*par[4]
    return ROOT.TMath.Sqrt(N_term + S_term + C_term)

class CompareJER_DifferentGroups(InputBase):
    def __init__(self):
        InputBase.__init__(self, inputdir = os.path.join('Path_JERC','Summer20year'))
        self.RMS = [0.68, 0.87, 0.95, 0.99]
        self.algos = ['ak4pfchs', 'ak4puppi']
        # self.algos = ['ak4pfchs']
        # self.years['UL'] = ['UL17', 'UL18']
        # self.years['UL'] = ['UL17']
        self.years['UL'] = ['UL18']
        self.etaBinsCommon = self.etaBinsCommon[0:1]
        self.LoadFiles()
        self.DoPlots()

    def LoadFiles(self):
        self.ref_pts = [21.5, 25.0, 28.5, 32.5, 37.5, 42.5, 51.0, 64.5, 81.0, 105.0, 135.0, 175.0, 250.0, 350.0, 475.0, 650.0, 875.0, 1250.0, 1750.0, 2250.0, 2750.0]
        self.graphs = {}
        self.LoadFiles_Ilias()
        self.LoadFiles_Hirak()
        self.LoadFiles_NonGaussTails()
        self.LoadFiles_DijetAnalysis_MCTruth()
        self.ExtractRC()

    def LoadFiles_NonGaussTails(self):
        group = 'NonGaussTails'
        # bins = [2,3,4,5,6,7,8,9,11,12,14,16,19,22,26,30,34,40,45,50,53]
        self.graphs[group] = {}
        for year in self.years['UL']:
            self.graphs[group][year] = {}
            for algo in self.algos:
                self.graphs[group][year][algo] = {}
                # inputdir = self.inputdir.replace('year', year).replace('nonAPV','')
                inputdir = self.inputdir.replace('year', 'UL18').replace('nonAPV','')
                f_ = ROOT.TFile(inputdir+'/JER_inputs/Convert_resHistos_to_resGraph.root')
                for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                    eta_name = GetEtaName(etaRef)
                    _, eta_max = GetEtaMinMax(etaRef)
                    if eta_max<= 0.522: gname = 'JER_vs_pt_EtaBin_0.0_0.5'
                    elif eta_max<= 0.783: gname = 'JER_vs_pt_EtaBin_0.5_0.8'
                    elif eta_max<= 1.044: gname = 'JER_vs_pt_EtaBin_0.8_1.1'
                    elif eta_max<= 1.305: gname = 'JER_vs_pt_EtaBin_1.1_1.3'
                    elif eta_max<= 1.740: gname = 'JER_vs_pt_EtaBin_1.3_1.7'
                    elif eta_max<= 1.930: gname = 'JER_vs_pt_EtaBin_1.7_1.9'
                    elif eta_max<= 2.172: gname = 'JER_vs_pt_EtaBin_1.9_2.1'
                    elif eta_max<= 2.322: gname = 'JER_vs_pt_EtaBin_2.1_2.3'
                    elif eta_max<= 2.500: gname = 'JER_vs_pt_EtaBin_2.3_2.5'
                    elif eta_max<= 2.853: gname = 'JER_vs_pt_EtaBin_2.5_2.8'
                    elif eta_max<= 2.964: gname = 'JER_vs_pt_EtaBin_2.8_3.0'
                    elif eta_max<= 3.139: gname = 'JER_vs_pt_EtaBin_2.8_3.0'
                    elif eta_max<= 3.489: gname = 'JER_vs_pt_EtaBin_2.1_2.3'
                    else: gname = 'JER_vs_pt_EtaBin_0.8_1.1'
                    graph = f_.Get('positiveEtaBins/'+gname)
                    max_ = len(self.ref_pts)
                    if eta_max>= 1.740: max_ = len(self.ref_pts)-1
                    if eta_max>= 1.930: max_ = len(self.ref_pts)-2
                    if eta_max>= 2.043: max_ = len(self.ref_pts)-3
                    if eta_max>= 2.322: max_ = len(self.ref_pts)-4
                    if eta_max>= 2.853: max_ = len(self.ref_pts)-5
                    if eta_max>= 2.964: max_ = len(self.ref_pts)-6
                    if eta_max>= 3.839: max_ = len(self.ref_pts)-8
                    if eta_max>= 5.191: max_ = len(self.ref_pts)-9
                    pts = self.ref_pts[0:max_]
                    # bins_ = bins[0:max_]
                    # self.graphs[group][year][algo][eta_name] = ROOT.TGraphErrors(len(pts), array('d',pts), array('d',[graph.Eval(pt) for pt in pts]), array('d',[0]*len(pts)), array('d',[graph.GetErrorY(bin) for bin in bins_]))
                    self.graphs[group][year][algo][eta_name] = ROOT.TGraphErrors(len(pts), array('d',pts), array('d',[graph.Eval(pt) for pt in pts]), array('d',[0]*len(pts)), array('d',[0]*len(pts)))
                f_.Close()

    def LoadFiles_DijetAnalysis_MCTruth(self):
        group = 'dijet'
        self.graphs[group] = {}
        for year in self.years['UL']:
            self.graphs[group][year] = {}
            for algo in self.algos:
                self.graphs[group][year][algo] = {}
                # f_ = ROOT.TFile(self.outdir+'/JER_MCTruth_Alex_'+algo+'_'+year+'.root')
                f_ = ROOT.TFile(self.outdir+'/JER_MCTruth_Alex_'+algo+'_'+'UL18'+'.root')
                for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                    eta_name = GetEtaName(etaRef)
                    self.graphs[group][year][algo][eta_name] = {}
                    self.graphs[group][year][algo][eta_name]['nominal'] = f_.Get('JER_MCTruth_00sigma_eta'+str(eta_bin+1))
                    self.graphs[group][year][algo][eta_name]['rms_0p68'] = f_.Get('JER_MCTruth_10sigma_eta'+str(eta_bin+1))
                    self.graphs[group][year][algo][eta_name]['rms_0p87'] = f_.Get('JER_MCTruth_15sigma_eta'+str(eta_bin+1))
                    self.graphs[group][year][algo][eta_name]['rms_0p95'] = f_.Get('JER_MCTruth_20sigma_eta'+str(eta_bin+1))
                    self.graphs[group][year][algo][eta_name]['rms_0p99'] = f_.Get('JER_MCTruth_30sigma_eta'+str(eta_bin+1))
                f_.Close()

    def LoadFiles_Ilias(self):
        group = 'Ilias'
        self.graphs[group] = {}
        for year in self.years['UL']:
            self.graphs[group][year] = {}
            f_ = ROOT.TFile('FitResponse_Ilias/MC_JER_mu_avg_'+year+'.root')
            for algo in self.algos:
                self.graphs[group][year][algo] = {}
                for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                    eta_name = GetEtaName(etaRef)
                    self.graphs[group][year][algo][eta_name] = OrderedDict()
                    for mode in ['CB', 'gauss', 'CI']:
                        for index, perc in enumerate(self.RMS):
                            jer_name = 'JER_'+mode+'_{:.2f}'.format(perc).replace('.','p')
                            self.graphs[group][year][algo][eta_name][jer_name] = f_.Get(algo+'/MCTruth_jer_mu_avg_MC_'+eta_name+'_'+jer_name)
            f_.Close()

    def LoadFiles_Hirak(self):
        group = 'Hirak'
        self.graphs[group] = {}
        for year in self.years['UL']:
            self.graphs[group][year] = {}
            inputdir = self.inputdir.replace('year', year).replace('nonAPV','')
            f_ = ROOT.TFile(inputdir+'/MC_JER_avg_mu/JER_MCtruth_avg_mu_'+year.replace('nonAPV','')+'.root')
            for algo in self.algos:
                self.graphs[group][year][algo] = {}
                for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                    eta_name = GetEtaName(etaRef)
                    eta_min, eta_max = GetEtaMinMax(etaRef)
                    graph = f_.Get(algo+'l1l2l3/RelResVsJetPt_JetEta'+(str(eta_min) if eta_min>0.1 else '0')+'to'+str(eta_max)+'_Mu0to60')
                    # self.graphs[group][year][algo][eta_name].GetListOfFunctions().RemoveLast()
                    max_ = len(self.ref_pts)
                    if eta_max>= 1.740: max_ = len(self.ref_pts)-1
                    if eta_max>= 1.930: max_ = len(self.ref_pts)-2
                    if eta_max>= 2.043: max_ = len(self.ref_pts)-3
                    if eta_max>= 2.322: max_ = len(self.ref_pts)-4
                    if eta_max>= 2.853: max_ = len(self.ref_pts)-5
                    if eta_max>= 2.964: max_ = len(self.ref_pts)-6
                    if eta_max>= 3.839: max_ = len(self.ref_pts)-8
                    if eta_max>= 5.191: max_ = len(self.ref_pts)-9
                    pts = self.ref_pts[0:max_]
                    self.graphs[group][year][algo][eta_name] = ROOT.TGraphErrors(len(pts), array('d',pts), array('d',[graph.Eval(pt) for pt in pts]), array('d',[0]*len(pts)), array('d',[0]*len(pts)))
            f_.Close()

    def ExtractRC(self, mode='sigrc'):
        with open('StudyNoiseTerm/NTerm.json') as json_file:
            self.NTerms = json.load(json_file)


    def CreateCanvas(self, name, year,eta_min,eta_max):
        tdr.cms_lumi_TeV = tdr.commonScheme['legend'][year]+' Legacy, '+commonScheme['lumi'][year]+' fb^{-1}'
        tdr.extraText3 = []
        tdr.extraText3.append(eta_min+' < |#eta| < '+eta_max)
        self.canv = tdrDiCanvas(name, 8, 3500, 0.0001, 0.6, 0.9, 1.5, 'p_{T}^{ptcl} [GeV]', 'JER', 'Ratio')
        self.canv.cd(2).SetLogx(True)
        tdrDrawLine(self.ref_line1, lcolor=ROOT.kBlack, lstyle=ROOT.kDashed)
        tdrDrawLine(self.ref_line2, lcolor=ROOT.kGray+1, lstyle=ROOT.kDotted)
        tdrDrawLine(self.ref_line3, lcolor=ROOT.kGray+1, lstyle=ROOT.kDotted)
        self.canv.cd(1).SetLogx(True)
        self.leg = tdrLeg(0.70,0.65,0.89,0.89, textSize=0.035)
        self.leg1 = tdrLeg(0.40,0.65,0.70,0.89, textSize=0.035)
        self.legPars = {}
        self.legPars['CB_95'] = tdrLeg(0.40,0.50,0.55,0.65, textSize=0.035)
        self.legPars['CB_87'] = tdrLeg(0.55,0.50,0.70,0.65, textSize=0.035)
        self.legPars['CB_68'] = tdrLeg(0.70,0.50,0.89,0.65, textSize=0.035)
        self.legPars['CI_95'] = tdrLeg(0.40,0.40,0.55,0.50, textSize=0.035)
        self.legPars['CI_87'] = tdrLeg(0.55,0.40,0.70,0.50, textSize=0.035)
        self.legPars['CI_68'] = tdrLeg(0.70,0.40,0.89,0.50, textSize=0.035)
        self.canv.cd(2)
        self.leg2 = tdrLeg(0.17,0.70,0.30,0.85, textSize=0.04)
        self.canv.cd(1)


    def SaveCanvas(self, year,algo,eta_name):
        self.canv.SaveAs(self.outdir+'/'+'_'.join(['JER_comparison',year,algo,eta_name])+'.pdf')
        del self.canv

    def CreateNSC(self, name, eta_name, year, color, graph, doLeg=False):
        funcs = {}
        fitRes = {}
        parNames = {0: 'N_{#mu=0}',1: 'N',2: 'S',3: 'd',4: 'C',5: '#mu'}
        mu = 32 if '18' in year else (33 if '17' in year else 23)
        modes = ['free', 'const', 'fix']
        modes = ['fix']
        val0, err0 = self.NTerms[year]['MC']['sigrc']['par0'][eta_name]
        val1, err1 = self.NTerms[year]['MC']['sigrc']['par1'][eta_name]
        for mode in modes:
            funcs[name+mode] = ROOT.TF1(name+mode, NSC, 7., 3550., 6)
            funcs[name+mode].FixParameter(0, val0)
            funcs[name+mode].FixParameter(1, val1)
            if mode=='free':
                funcs[name+mode].SetLineStyle(ROOT.kDotted)
            if mode=='const':
                funcs[name+mode].SetParameter(5, mu)
                funcs[name+mode].SetParLimits(5, mu-5, mu+5)
                funcs[name+mode].SetLineStyle(ROOT.kSolid)
            if mode=='fix':
                funcs[name+mode].FixParameter(5, mu)
                funcs[name+mode].SetLineStyle(ROOT.kDashed)
            funcs[name+mode].SetLineColor(color)
            fitRes[mode] = graph.Fit(funcs[name+mode], 'RMQS+')
            ROOT.gStyle.SetOptFit(0)
        if doLeg:
            lname = 'CB_' if 'CB' in name else 'CI_'
            lname += '68' if '68' in name else ('87' if '87' in name else '95')
            for mode in modes:
                self.legPars[lname].AddEntry(funcs[name+mode], lname, 'l')
            for par in range(2,4):
                for mode in modes:
                    val = round(funcs[name+mode].GetParameter(par),2)
                    if par!=3: val = math.fabs(val)
                    self.legPars[lname].AddEntry(ROOT.TObject(), parNames[par]+' = '+str(val), 'P')
        return fitRes

    def DoPlots(self):
        colors ={
            # 'nominal':  ROOT.kBlack,
            '0p68': ROOT.kRed+1,
            '0p87': ROOT.kOrange+1,
            '0p95': ROOT.kGreen+2,
            '0p99': ROOT.kAzure+2,
        }
        self.ref_line1 = ROOT.TLine(8, 1.00, 3500, 1.00)
        self.ref_line2 = ROOT.TLine(8, 1.05, 3500, 1.05)
        self.ref_line3 = ROOT.TLine(8, 0.95, 3500, 0.95)
        ksenia = {}
        for eta_name in ['0p000_0p261', '0p261_0p522', '0p522_0p783', '0p783_1p044', '1p044_1p305']:
            ksenia[eta_name] = ROOT.TGraph(8, array('d',[35.067,104.560,1086.612,24.820,64.741,225.645,455.449,2121.107]), array('d',[0.1615,0.1089,0.0487,0.1910,0.1256,0.0807,0.0641,0.0423]))
        for eta_name in ['1p305_1p566', '1p566_1p740', '1p740_1p930']:
            ksenia[eta_name] = ROOT.TGraph(8, array('d',[25.0990,35.0672,65.4672,104.5600,275.7869,544.3815,1098.7933,1501.3335]), array('d',[0.2474,0.2012,0.1474,0.1153,0.0730,0.0589,0.0474,0.0435]))
        for eta_name in ['1p930_2p043', '2p043_2p172', '2p172_2p322', '2p322_2p500']:
            ksenia[eta_name] = ROOT.TGraph(8, array('d',[25.0990,35.0672,65.4672,104.5600,275.7869,544.3815,1098.7933,1501.3335]), array('d',[0.2474,0.2012,0.1474,0.1153,0.0730,0.0589,0.0474,0.0435]))
        for eta_name in ['2p500_2p650', '2p650_2p853', '2p853_2p964']:
            ksenia[eta_name] = ROOT.TGraph(17, array('d',[24.8208,25.0990,34.6785,34.6785,44.8140,45.3164,54.7722,55.3862,64.7414,65.4672,74.8379,74.8379,84.6014,84.6014,105.7322,275.7869,695.6899]), array('d',[0.3405,0.4392,0.2569,0.3215,0.2088,0.2316,0.1835,0.1962,0.1632,0.1772,0.1430,0.1569,0.1367,0.1493,0.1240,0.0784,0.0531]))
        for eta_name in ['2p964_3p139', '3p139_3p489', '3p489_3p839', '3p839_5p191']:
            ksenia[eta_name] = ROOT.TGraph(17, array('d',[24.8208,25.0990,34.6785,34.6785,44.8140,45.3164,54.7722,55.3862,64.7414,65.4672,74.8379,74.8379,84.6014,84.6014,105.7322,275.7869,695.6899]), array('d',[0.3405,0.4392,0.2569,0.3215,0.2088,0.2316,0.1835,0.1962,0.1632,0.1772,0.1430,0.1569,0.1367,0.1493,0.1240,0.0784,0.0531]))


        for year in self.years['UL']:
            for algo in self.algos:
                for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                    eta_name = GetEtaName(etaRef)
                    eta_min, eta_max = eta_name.replace('p','.').split('_')
                    self.CreateCanvas(eta_name,year,eta_min,eta_max)
                    funcs = {}
                    ref_graph = self.graphs['Hirak'][year][algo][eta_name]
                    # ref_graph = self.graphs['Ilias'][year][algo][eta_name]['rms_0p95']
                    # ref_graph = self.graphs['Ilias'][year][algo][eta_name]['nominal']
                    color = ROOT.kOrange+1
                    # tdrDraw(ksenia[eta_name], 'p', marker=ROOT.kFullTriangleUp, mcolor=color)
                    pts = []
                    jers = []
                    # for bin in range(ksenia[eta_name].GetN()):
                    #     pt, jer = (ROOT.Double(0),ROOT.Double(0))
                    #     ksenia[eta_name].GetPoint(int(bin),pt,jer)
                    #     pts.append(pt)
                    #     jers.append(jer/ref_graph.Eval(pt))
                    # ksenia[eta_name+'ratio'] = ROOT.TGraph(len(pts), array('d',pts), array('d',jers))
                    # self.canv.cd(2)
                    # tdrDraw(ksenia[eta_name+'ratio'], 'p', marker=ROOT.kFullCircle, mcolor=color)
                    group, color = ('Hirak',ROOT.kViolet+1)
                    graph = self.graphs[group][year][algo][eta_name]
                    fitRes = self.CreateNSC(name=group+year+algo+eta_name, eta_name=eta_name, year=year, color=color, graph=graph)
                    self.canv.cd(1)
                    style = ROOT.kFullDiamond
                    tdrDraw(graph, 'p', marker=style, mcolor=color)
                    self.leg1.AddEntry(graph, group, 'p')
                    self.graphs[group][year][algo][eta_name+'ratio'] = TGraphRatio(graph, ref_graph)
                    self.canv.cd(2)
                    tdrDraw(self.graphs[group][year][algo][eta_name+'ratio'], 'p', marker=style, mcolor=color)
                    group = 'Ilias'
                    for mode, graph in self.graphs[group][year][algo][eta_name].items():
                        color = colors[mode.split('_')[-1]]
                        style = ROOT.kFullCircle if 'CB' in mode else (ROOT.kFullTriangleUp if 'CI' in mode else ROOT.kFullTriangleDown)
                        lstyle = ROOT.kSolid if 'CB' in mode else (ROOT.kDashed if 'CI' in mode else ROOT.kDotted)
                        fitRes = self.CreateNSC(name=group+year+algo+eta_name+mode, eta_name=eta_name, year=year, color=color, graph=graph, doLeg=(not 'gaus' in mode and not '99' in mode))
                        self.canv.cd(1)
                        tdrDraw(graph, 'p', marker=style, mcolor=color)
                        if 'CB' in mode:
                            self.leg.AddEntry(graph, mode.split('_')[-1].replace('p','.'), 'p')
                        if '68' in mode:
                            self.leg1.AddEntry(graph, mode.split('_')[1], 'p')
                        self.canv.cd(2)
                        # if mode != 'nominal' and not '95' in mode: continue
                        self.graphs[group][year][algo][eta_name][mode+'ratio'] = TGraphRatio(graph, ref_graph)
                        tdrDraw(self.graphs[group][year][algo][eta_name][mode+'ratio'], 'l', marker=style, mcolor=color, lstyle=lstyle, lcolor=color)
                        if '68' in mode:
                            self.leg2.AddEntry(self.graphs[group][year][algo][eta_name][mode+'ratio'], mode.split('_')[1], 'l')
                    group, color = ('NonGaussTails',ROOT.kGray+1)
                    graph = self.graphs[group][year][algo][eta_name]
                    # fitRes = self.CreateNSC(name=group+year+algo+eta_name, eta_name=eta_name, year=year, color=color, graph=graph)
                    self.canv.cd(1)
                    style = ROOT.kFullSquare
                    tdrDraw(graph, 'p', marker=style, mcolor=color)
                    self.leg1.AddEntry(graph, 'core from CB fit', 'p')
                    self.graphs[group][year][algo][eta_name+'ratio'] = TGraphRatio(graph, ref_graph)
                    self.canv.cd(2)
                    tdrDraw(self.graphs[group][year][algo][eta_name+'ratio'], 'p', marker=style, mcolor=color)
                    # group = 'dijet'
                    # for mode, graph in self.graphs[group][year][algo][eta_name].items():
                    #     color = colors[mode]
                    #     self.canv.cd(1)
                    #     tdrDraw(graph, 'l', lcolor=color, lstyle=ROOT.kDashed)
                    #     if mode == 'nominal':
                    #         self.leg1.AddEntry(graph, 'MCTruth from '+group, 'l')
                    self.SaveCanvas(year=year,algo=algo,eta_name=eta_name)


        algo, algo_ref = ('ak4pfchs', 'ak4puppi')
        for year in self.years['UL']:
            for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                eta_name = GetEtaName(etaRef)
                eta_min, eta_max = eta_name.replace('p','.').split('_')
                self.CreateCanvas(eta_name, year,eta_min,eta_max)
                # tdrDraw(ksenia[eta_name], 'p', marker=ROOT.kFullTriangleUp, mcolor=ROOT.kOrange+1)
                group, color = ('Hirak',ROOT.kViolet+1)
                graph = self.graphs[group][year][algo][eta_name]
                graph_ref = self.graphs[group][year][algo_ref][eta_name]
                graph_ratio = TGraphRatio(graph, graph_ref)
                self.graphs[group][year][algo][eta_name+'ratio_algo'] = graph_ratio
                self.canv.cd(1)
                tdrDraw(graph_ref, 'l', lcolor=color, lstyle=ROOT.kDotted)
                tdrDraw(graph, 'l', lcolor=color, lstyle=ROOT.kSolid)
                self.leg.AddEntry(graph, group, 'l')
                self.canv.cd(2)
                tdrDraw(graph_ratio, 'p', marker=ROOT.kFullCircle, mcolor=color)
                for mode in ['rms_0p95', 'rms_0p68']:
                    group, color = ('Ilias', colors[mode.split('_')[-1]])
                    graph = self.graphs[group][year][algo][eta_name][mode]
                    graph_ref = self.graphs[group][year][algo_ref][eta_name][mode]
                    graph_ratio = TGraphRatio(graph, graph_ref)
                    self.graphs[group][year][algo][eta_name][mode+'ratio_algo'] = graph_ratio
                    self.leg.AddEntry(graph, group+' '+mode, 'l')
                    if mode=='rms_0p68':
                        self.leg1.AddEntry(graph, 'chs', 'l')
                        self.leg1.AddEntry(graph_ref, 'Puppi', 'l')
                    self.canv.cd(1)
                    tdrDraw(graph_ref, 'l', lcolor=color, lstyle=ROOT.kDotted)
                    tdrDraw(graph, 'l', lcolor=color, lstyle=ROOT.kSolid)
                    self.canv.cd(2)
                    tdrDraw(graph_ratio, 'p', marker=ROOT.kFullCircle, mcolor=color)
                self.SaveCanvas(year=year,algo='chs_vs_puppi',eta_name=eta_name)


        year_ref = 'UL18'
        for algo in self.algos:
            for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                eta_name = GetEtaName(etaRef)
                eta_min, eta_max = eta_name.replace('p','.').split('_')
                self.CreateCanvas(eta_name,'Run2',eta_min,eta_max)
                for year in self.years['UL']:
                    color = commonScheme['color'][year]
                    group, style = ('Hirak',ROOT.kDashed)
                    graph = self.graphs[group][year][algo][eta_name]
                    graph_ratio = TGraphRatio(graph, self.graphs[group][year_ref][algo][eta_name])
                    self.graphs[group][year][algo][eta_name+'ratio_years'] = graph_ratio
                    if year =='UL18':
                        self.leg1.AddEntry(graph, group, 'l')
                    self.canv.cd(1)
                    tdrDraw(graph, 'l', lcolor=color, lstyle=style)
                    self.canv.cd(2)
                    tdrDraw(graph_ratio, 'l', lcolor=color, lstyle=style)
                    group = 'Ilias'
                    for mode in ['rms_0p68', 'rms_0p95']:
                        style = ROOT.kSolid if mode=='rms_0p68' else ROOT.kDotted
                        graph = self.graphs[group][year][algo][eta_name][mode]
                        graph_ratio = TGraphRatio(graph, self.graphs[group][year_ref][algo][eta_name][mode])
                        self.graphs[group][year][algo][eta_name][mode+'ratio_years'] = graph_ratio
                        if mode == 'rms_0p68':
                            self.leg.AddEntry(graph, year, 'l')
                        if year =='UL18':
                            self.leg1.AddEntry(graph, group+' '+mode, 'l')
                        self.canv.cd(1)
                        tdrDraw(graph, 'l', lcolor=color, lstyle=style)
                        self.canv.cd(2)
                        tdrDraw(graph_ratio, 'l', lcolor=color, lstyle=style)
                self.SaveCanvas(year='years',algo=algo,eta_name=eta_name)

def main():
    CompareJER_DifferentGroups()

if __name__ == '__main__':
    main()
