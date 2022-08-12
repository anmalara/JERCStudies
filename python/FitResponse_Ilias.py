from Utils import *
from GraphHistUtils import TGraphRatio, ComputeHistWithCL
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gErrorIgnoreLevel = ROOT.kError

tdr.writeExtraText = True
tdr.extraText = 'Simulation'
tdr.extraText2 = 'Work in progress'

colors = {0.68: ROOT.kRed+1,
          0.87: ROOT.kOrange+1,
          0.95: ROOT.kGreen+2,
          0.99: ROOT.kAzure+2,
          }

# def CrystalBall(x, par):
#     return ROOT.Math.crystalball_function(x[0], par[0], par[1], par[2], par[3])*par[4]

def SetParLimits(func, pt_avg, sigma, N):
    func.SetParameter(0,1.)
    func.SetParameter(1,sigma)
    func.SetParameter(4,4)
    func.SetParameter(5,4)
    func.SetParameter(6,N)
    if pt_avg>1500:
        func.SetParameter(2,1.3)
        func.SetParameter(3,2.0)
    elif pt_avg>400:
        func.SetParameter(2,1.8)
        func.SetParameter(3,3.0)
    elif pt_avg>200:
        func.SetParameter(2,1.4)
        func.SetParameter(3,2.0)
    elif pt_avg>30:
        func.SetParameter(2,0.6)
        func.SetParameter(3,0.8)
    else:
        func.SetParameter(2,0.4)
        func.SetParameter(3,1.6)
    func.SetParLimits(0,0.9,1.1)
    func.SetParLimits(1,0.01,1)
    func.SetParLimits(2,0.01,5)
    func.SetParLimits(3,0.01,5)
    func.SetParLimits(4,0.01,20)
    func.SetParLimits(5,0.01,20)

def DoubleSidedCrystalBall(x, par):
    if any([math.isnan(par[i]) for i in range(6)]):
        return 0
    mean = par[0]
    sigma = par[1]
    alpha_l = par[2]
    alpha_h = par[3]
    n_l, n_h = (par[4],par[5])
    N = par[6]
    t = (x[0]-mean)/sigma
    fact1TLessMinosAlphaL = alpha_l/n_l if n_l!=0 else 0
    fact2TLessMinosAlphaL = (n_l/alpha_l) - alpha_l -t
    fact1THihgerAlphaH = alpha_h/n_h if n_h!=0 else 0
    fact2THigherAlphaH = (n_h/alpha_h) - alpha_h +t
    if (-alpha_l <= t and t <= alpha_h):
        result = math.exp(-0.5*t*t)
    elif (t < -alpha_l):
        result = math.exp(-0.5*alpha_l*alpha_l)*math.pow(math.fabs(fact1TLessMinosAlphaL*fact2TLessMinosAlphaL), -n_l)
    elif (t > alpha_h):
        result = math.exp(-0.5*alpha_h*alpha_h)*math.pow(math.fabs(fact1THihgerAlphaH*fact2THigherAlphaH), -n_h)
    return N*result


def Confidence(h1, median, confLevel = 0.95):
    ix = h1.GetXaxis().FindBin(median)
    ixlow = ix
    ixhigh = ix
    nb = h1.GetNbinsX()
    ntot = h1.Integral()
    nsum = h1.GetBinContent(ix)
    width = h1.GetBinWidth(ix)
    if ntot==0: return 0
    while (nsum < confLevel * ntot):
        nlow = h1.GetBinContent(ixlow-1) if ixlow>0 else 0
        nhigh = h1.GetBinContent(ixhigh+1) if ixhigh<nb else 0
        if (nsum+max(nlow,nhigh) < confLevel * ntot):
            if (nlow>=nhigh and ixlow>0):
                nsum += nlow
                ixlow -=1
                width += h1.GetBinWidth(ixlow)
            elif ixhigh<nb:
                nsum += nhigh
                ixhigh+=1
                width += h1.GetBinWidth(ixhigh)
            else: raise ValueError('BOOM')
        else:
            if (nlow>nhigh):
                width += h1.GetBinWidth(ixlow-1) * (confLevel * ntot - nsum) / nlow
            else:
                width += h1.GetBinWidth(ixhigh+1) * (confLevel * ntot - nsum) / nhigh
            nsum = ntot
    return width

def GetBin(h, x_val, doMin):
    shift = 0.001
    x_axis = h.GetXaxis()
    if doMin:
        bin = x_axis.FindBin(x_val+shift)
        return x_axis.FindBin(x_axis.GetBinCenter(bin)-x_axis.GetBinWidth(bin)/2+shift)
    else:
        bin = x_axis.FindBin(x_val-shift)
        return x_axis.FindBin(x_axis.GetBinCenter(bin)+x_axis.GetBinWidth(bin)/2-shift)


class FitResponse_Ilias(InputBase):
    def __init__(self):
        InputBase.__init__(self, inputdir = os.path.join('Path_JERC','Summer20year','L2Relative'))
        self.algos = ['ak4pfchs', 'ak4puppi']
        self.RMS = [0.99, 0.95, 0.87, 0.68]
        self.pars = OrderedDict([(i,OrderedDict()) for i in range(7)])
        self.years['UL'] = ['UL18']
        self.etaBinsCommon = self.etaBinsCommon[0:1]
        self.pt_bins = ['RefPt20to23', 'RefPt23to27', 'RefPt27to30', 'RefPt30to35', 'RefPt35to40', 'RefPt40to45', 'RefPt45to57', 'RefPt57to72', 'RefPt72to90', 'RefPt90to120', 'RefPt120to150', 'RefPt150to200', 'RefPt200to300', 'RefPt300to400', 'RefPt400to550', 'RefPt550to750', 'RefPt750to1000', 'RefPt1000to1500', 'RefPt1500to2000', 'RefPt2000to2500', 'RefPt2500to3000']
        # self.pt_bins = ['RefPt20to23', 'RefPt23to27', 'RefPt27to30', 'RefPt30to35', 'RefPt35to40', 'RefPt40to45']
        # self.pt_bins = ['RefPt45to57', 'RefPt57to72', 'RefPt72to90',]
        for year in self.years['UL']:
            print year
            self.LoadFiles(year)
        print 'check'
        for pi in range(7):
            pdfName = str(pi)
            PlotYMax = 60 if pi>3 else (20 if pi>1 else (2 if pi==0 else 0.5))
            canv = tdrCanvas(pdfName, 0, 3500, 0, PlotYMax,'Response', 'A.U.')
            canv.SetLogx(1)
            graphs = []
            for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                x_val,y_val =([],[])
                for key, data in self.pars[pi].items():
                    if key[1]!=etaRef: continue
                    if len(data)==0: continue
                    data = np.array(list(set([x for x in data if math.isnan(x) == False])))
                    mean, std = (np.mean(data), np.std(data) if len(data)>1 else 0)
                    data = data[abs(data-mean<2*std)]
                    if len(data)==0: continue
                    mean, std = (np.mean(data), np.std(data) if len(data)>1 else 0)
                    pt_min, pt_max = key[0].replace('RefPt','').split('to')
                    pt_avg = (float(pt_min)+float(pt_max))/2
                    for d in data:
                        x_val.append(pt_avg)
                        y_val.append(d)
                if len(y_val)==0: continue
                graph = ROOT.TGraph(len(x_val), array('d',x_val), array('d',y_val))
                tdrDraw(graph, 'lp')
                graphs.append(graph)
            canv.SaveAs(self.outdir+'/'+pdfName+'.pdf')
            del canv

    def LoadFiles(self, year):
        inputdir = self.inputdir.replace('year', year).replace('nonAPV','')
        f_ = {}
        f_['ak4pfchs'] = ROOT.TFile(inputdir+'/CorrectedResponses_'+year+'_AK4CHS.root')
        f_['ak4puppi'] = ROOT.TFile(inputdir+'/CorrectedResponses_'+year+'_AK4PUPPI.root')

        JER  = {'ak4pfchs':OrderedDict(),'ak4puppi':OrderedDict()}

        for pt_ref in self.pt_bins:
            pt_min, pt_max = pt_ref.replace('RefPt','').split('to')
            pt_avg = (float(pt_min)+float(pt_max))/2
            print pt_ref, pt_avg
            for algo in self.algos:
                h = f_[algo].Get(algo+'/RelRspVsJetEta_'+pt_ref)
                for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                    eta_min, eta_max = GetEtaMinMax(etaRef)
                    if float(pt_max)*math.cosh(eta_min)>6500: continue
                    name = eta_name = GetEtaName(etaRef)
                    for pi in range(7):
                        self.pars[pi][(pt_ref,etaRef)] =[]
                    # print pt_avg, algo, name
                    if not name in JER[algo]: JER[algo][name] = {}
                    bin_min = GetBin(h=h, x_val=eta_min, doMin=True)
                    bin_max = GetBin(h=h, x_val=eta_max, doMin=False)
                    h_proj = h.ProjectionY(name+'pos', bin_min, bin_max, 'e')
                    bin_min = GetBin(h=h, x_val=-eta_max, doMin=True)
                    bin_max = GetBin(h=h, x_val=-eta_min, doMin=False)
                    h_proj_neg = h.ProjectionY(name+'neg', bin_min, bin_max, 'e')
                    h_proj.Add(h_proj_neg)
                    if h_proj.GetEntries()<10:
                        print 'skipping', pt_avg, algo, name, pt_max, eta_min, float(pt_max)*math.cosh(eta_min)
                        continue
                    JER[algo][name].setdefault('pt', []).append(pt_avg)
                    xqs,yqs,funcs,widths = ({},{},{},{})
                    skipPlot = True
                    for index, perc in enumerate(self.RMS):
                        xq = array('d', [(1.00-perc)/2,(1.00+perc)/2])
                        yq = array('d', [0,0])
                        h_proj.GetQuantiles(2, yq, xq)
                        xqs[perc] = xq
                        yqs[perc] = yq
                        func_name = 'CI'+str(perc)
                        widths[func_name] = Confidence(h_proj, h_proj.GetMean(), confLevel = perc)/(2*[2.576,1.960,1.514,0.9945][index])
                        JER[algo][name].setdefault('JER_CI_'+'{:.2f}'.format(perc).replace('.','p'), []).append(widths[func_name])
                        JER[algo][name].setdefault('err_CI_'+'{:.2f}'.format(perc).replace('.','p'), []).append(h_proj.GetRMSError())
                        func_name = 'gauss'+str(perc)
                        funcs[func_name] = ROOT.TF1('gauss'+str(perc),'gaus', yq[0], yq[1])
                        h_proj.Fit(funcs[func_name], 'RQMS')
                        widths[func_name] = funcs[func_name].GetParameter(2)
                        sigma = widths[func_name]
                        JER[algo][name].setdefault('JER_gauss_'+'{:.2f}'.format(perc).replace('.','p'), []).append(widths[func_name])
                        JER[algo][name].setdefault('err_gauss_'+'{:.2f}'.format(perc).replace('.','p'), []).append(funcs[func_name].GetParError(2))
                        func_name = 'CB'+str(perc)
                        xq_ref = array('d', [(1.00-0.90)/2,(1.00+0.90)/2])
                        yq_ref = array('d', [0,0])
                        h_proj.GetQuantiles(2, yq_ref,xq_ref)
                        funcs['CB_ref'] = ROOT.TF1('CB_ref',DoubleSidedCrystalBall,yq_ref[0], yq_ref[1],7)
                        funcs[func_name] = ROOT.TF1('CB'+str(perc),DoubleSidedCrystalBall,yq[0], yq[1],7)
                        SetParLimits(funcs['CB_ref'], pt_avg, sigma, h_proj.Integral())
                        SetParLimits(funcs[func_name], pt_avg, sigma, h_proj.Integral())
                        h_proj.Fit(funcs['CB_ref'], 'RQMLS')
                        h_proj.Fit(funcs[func_name], 'RQMLS')
                        chi2_red = funcs[func_name].GetChisquare()/funcs[func_name].GetNDF()
                        if index!=0 and chi2_red>20:
                            for i in range(7):
                                funcs[func_name].SetParameter(i,funcs['CB'+str(self.RMS[0])].GetParameter(i))
                            h_proj.Fit(funcs[func_name], 'RQMLS')
                            chi2_red = funcs[func_name].GetChisquare()/funcs[func_name].GetNDF()
                        if chi2_red>20:
                            for i in range(7):
                                funcs[func_name].SetParameter(i,funcs['CB_ref'].GetParameter(i))
                            h_proj.Fit(funcs[func_name], 'RQMLS')
                            chi2_red = funcs[func_name].GetChisquare()/funcs[func_name].GetNDF()
                        counter = 0
                        while funcs[func_name].GetChisquare()/funcs[func_name].GetNDF()>20 and counter<20:
                            h_proj.Fit(funcs[func_name], 'RQMLS')
                            counter += 1
                        if funcs[func_name].GetChisquare()/funcs[func_name].GetNDF()>50:
                            print pt_ref, algo,name, perc, funcs[func_name].GetChisquare()/funcs[func_name].GetNDF()
                            skipPlot = False
                        for pi in range(7):
                            self.pars[pi][(pt_ref,etaRef)].append(funcs[func_name].GetParameter(pi))
                        #     print round(funcs[func_name].GetParameter(pi),2),
                        # print round(funcs[func_name].GetChisquare()/funcs[func_name].GetNDF(),2)
                        widths[func_name] = funcs[func_name].GetParameter(1)
                        JER[algo][name].setdefault('JER_CB_'+'{:.2f}'.format(perc).replace('.','p'), []).append(widths[func_name])
                        JER[algo][name].setdefault('err_CB_'+'{:.2f}'.format(perc).replace('.','p'), []).append(funcs[func_name].GetParError(1))
                    if skipPlot: continue
                    tdr.cms_lumi_TeV = tdr.commonScheme['legend'][year]+' Legacy, '+commonScheme['lumi'][year]+' fb^{-1}'
                    tdr.extraText3 = []
                    tdr.extraText3.append(name.replace('p','.').split('_')[0]+' < |#eta| < '+name.replace('p','.').split('_')[1])
                    tdr.extraText3.append(pt_min+' < p_{T} < '+pt_max)
                    PlotYMax = h_proj.GetMaximum()*1.2
                    LineYMax = h_proj.GetMaximum()*0.75
                    canv = tdrCanvas(pt_ref+algo+name, 0, 2, 0, PlotYMax,'Response', 'A.U.')
                    leg = tdrLeg(0.80,0.65,0.95,0.89, textSize=0.035)
                    leg_ref = tdrLeg(0.60,0.70,0.80,0.89, textSize=0.035)
                    tdrDraw(h_proj,'')
                    lines = {}
                    for perc in self.RMS:
                        xq = xqs[perc]
                        yq = yqs[perc]
                        color = colors[perc]
                        for par in [0,1]:
                            lines[str(perc)+str(par)] = rt.TLine(yq[par], 0.00001, yq[par], LineYMax)
                            tdrDrawLine(lines[str(perc)+str(par)], lcolor=color, lstyle=ROOT.kDotted)
                        for mode in ['CB', 'gauss', 'CI']:
                            isRef = 'CB'==mode
                            func_name = mode+str(perc)
                            hasFit = func_name in funcs
                            style = ROOT.kSolid if isRef else (ROOT.kDashed if hasFit else ROOT.kDotted)
                            width = widths[func_name]
                            if hasFit:
                                funcs[func_name].SetLineColor(color)
                                funcs[func_name].SetLineStyle(style)
                                funcs[func_name].Draw('same')
                                if isRef:
                                    leg.AddEntry(funcs[func_name], str(perc), 'l')
                            height = LineYMax*(1-(perc if perc!=0.68 else 0.80)*(0.4 if isRef else 0.6 if hasFit else 0.8))
                            lines[func_name] = rt.TLine(1-width, height, 1+width, height)
                            tdrDrawLine(lines[func_name], lcolor=color, lstyle=style)
                            if perc == self.RMS[0]:
                                lines[func_name+'leg']= ROOT.TLine()
                                tdrDrawLine(lines[func_name+'leg'], lcolor=ROOT.kBlack, lstyle=style)
                                leg_ref.AddEntry(lines[func_name+'leg'], str(mode), 'l')
                    ROOT.gStyle.SetOptFit(0)
                    canv.SaveAs(self.outdir+'/'+'_'.join(['JER_fit',year,algo,pt_ref,name])+'.pdf')
                    del h_proj
                    del canv
                    del funcs



        for fname in f_:
            f_[fname].Close()
        f_output = ROOT.TFile(self.outdir+'/MC_JER_mu_avg_'+year+'.root', 'RECREATE')
        for algo in self.algos:
            f_output.mkdir(algo)
            f_output.cd(algo)
            for eta_name in JER[algo].keys():
                x = JER[algo][eta_name]['pt']
                if len(x)==0: continue
                for mode in ['CB', 'gauss', 'CI']:
                    for index, perc in enumerate(self.RMS):
                        jer_name = 'JER_'+mode+'_{:.2f}'.format(perc).replace('.','p')
                        y = JER[algo][eta_name][jer_name]
                        err = JER[algo][eta_name][jer_name.replace('JER_', 'err_')]
                        graph = ROOT.TGraphErrors(len(x), array('d',x), array('d',y), array('d',[0]*len(x)), array('d',err))
                        graph.Write('MCTruth_jer_mu_avg_MC_'+eta_name+'_'+jer_name)

        f_output.Write()
        f_output.Close()

def main():
    FitResponse_Ilias()

if __name__ == '__main__':
    main()
