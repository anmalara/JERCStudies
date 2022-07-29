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
        # self.years['UL'] = ['UL18']
        for year in self.years['UL']:
            print year
            self.LoadFiles(year)

    def LoadFiles(self, year):
        inputdir = self.inputdir.replace('year', year).replace('nonAPV','')
        f_ = {}
        f_['ak4pfchs'] = ROOT.TFile(inputdir+'/CorrectedResponses_'+year+'_AK4CHS.root')
        f_['ak4puppi'] = ROOT.TFile(inputdir+'/CorrectedResponses_'+year+'_AK4PUPPI.root')

        JER  = {'ak4pfchs':OrderedDict(),'ak4puppi':OrderedDict()}

        pt_bins = ['RefPt20to23', 'RefPt23to27', 'RefPt27to30', 'RefPt30to35', 'RefPt35to40', 'RefPt40to45', 'RefPt45to57', 'RefPt57to72', 'RefPt72to90', 'RefPt90to120', 'RefPt120to150', 'RefPt150to200', 'RefPt200to300', 'RefPt300to400', 'RefPt400to550', 'RefPt550to750', 'RefPt750to1000', 'RefPt1000to1500', 'RefPt1500to2000', 'RefPt2000to2500', 'RefPt2500to3000']
        for pt_ref in pt_bins:
            print pt_ref
            pt_min, pt_max = pt_ref.replace('RefPt','').split('to')
            pt_avg = (float(pt_min)+float(pt_max))/2
            for algo in self.algos:
                h = f_[algo].Get(algo+'/RelRspVsJetEta_'+pt_ref)
                for eta_bin, etaRef in enumerate(self.etaBinsCommon):
                    eta_min, eta_max = (JERC_Constants.GetEtaBinEdgeMin(etaRef),JERC_Constants.GetEtaBinEdgeMax(etaRef))
                    if float(pt_max)*math.cosh(eta_min)>6500: continue
                    name = '{:.3f}_{:.3f}'.format(eta_min,eta_max).replace('.','p')
                    if not name in JER[algo]: JER[algo][name] = {}
                    bin_min = GetBin(h=h, x_val=eta_min, doMin=True)
                    bin_max = GetBin(h=h, x_val=eta_max, doMin=False)
                    # x_min = h.GetXaxis().GetBinCenter(bin_min)-h.GetXaxis().GetBinWidth(bin_min)/2
                    # x_max = h.GetXaxis().GetBinCenter(bin_max)+h.GetXaxis().GetBinWidth(bin_max)/2
                    # print round(eta_min-x_min,2), round(eta_max-x_max,2),
                    h_proj = h.ProjectionY(name+'pos', bin_min, bin_max, 'e')
                    bin_min = GetBin(h=h, x_val=-eta_max, doMin=True)
                    bin_max = GetBin(h=h, x_val=-eta_min, doMin=False)
                    # x_min = h.GetXaxis().GetBinCenter(bin_min)-h.GetXaxis().GetBinWidth(bin_min)/2
                    # x_max = h.GetXaxis().GetBinCenter(bin_max)+h.GetXaxis().GetBinWidth(bin_max)/2
                    # print round(eta_max+x_min,2), round(eta_min+x_max,2)
                    h_proj_neg = h.ProjectionY(name+'neg', bin_min, bin_max, 'e')
                    h_proj.Add(h_proj_neg)
                    if h_proj.GetEntries()<10:
                        print 'skipping', pt_avg, algo, name, pt_max, eta_min, float(pt_max)*math.cosh(eta_min)
                        continue
                    x_min, x_max = (h_proj.GetMean()-h_proj.GetRMS(),h_proj.GetMean()+h_proj.GetRMS())
                    gauss = ROOT.TF1('gauss','gaus',x_min,x_max)
                    funcs = {}
                    lines = {}
                    xqs = {}
                    yqs = {}
                    rms = {}
                    for perc in self.RMS:
                        h_proj_trunc = h_proj.Clone('h_proj_trunc_'+str(perc))
                        xq = array('d', [(1.00-perc)/2,(1.00+perc)/2])
                        yq = array('d', [0,0])
                        h_proj_trunc.GetQuantiles(2, yq, xq)
                        h_proj_trunc.GetXaxis().SetRange(h_proj_trunc.FindBin(yq[0]), h_proj_trunc.FindBin(yq[1]));
                        funcs[perc] = ROOT.TF1('gauss'+str(perc),'gaus', yq[0], yq[1])
                        fit_res = h_proj_trunc.Fit(funcs[perc], 'RQMS+', '')
                        JER[algo][name].setdefault('JER_rms_'+'{:.2f}'.format(perc).replace('.','p'), []).append(funcs[perc].GetParameter(2))
                        JER[algo][name].setdefault('err_rms_'+'{:.2f}'.format(perc).replace('.','p'), []).append(funcs[perc].GetParError(2))
                        xqs[perc] = xq
                        yqs[perc] = yq
                        rms[perc] = funcs[perc].GetParameter(2)
                    iterations = 1000
                    last_jer = 100
                    min_diff = 10
                    counter = 0
                    for i_ in range(iterations+1):
                        fit_res = h_proj.Fit(gauss, 'QMS+', '', x_min, x_max)
                        current_jer, chi2 = (gauss.GetParameter(2), gauss.GetChisquare())
                        diff = math.fabs(last_jer-current_jer)
                        if diff<min_diff: min_diff, counter = (diff,0)
                        else: counter += 1
                        if counter>10: break
                        if diff<1e-05: break
                        if i_!=iterations:
                            h_proj.GetListOfFunctions().RemoveLast()
                        last_jer, x_min, x_max = (current_jer, h_proj.GetMean()-current_jer,h_proj.GetMean()+current_jer)
                    JER[algo][name].setdefault('pt', []).append(pt_avg)
                    JER[algo][name].setdefault('JER_nominal', []).append(current_jer)
                    JER[algo][name].setdefault('err_nominal', []).append(gauss.GetParError(2))
                    tdr.cms_lumi_TeV = tdr.commonScheme['legend'][year]+' Legacy, '+commonScheme['lumi'][year]+' fb^{-1}'
                    tdr.extraText3 = []
                    tdr.extraText3.append(name.replace('p','.').split('_')[0]+' < |#eta| < '+name.replace('p','.').split('_')[1])
                    tdr.extraText3.append(pt_min+' < p_{T} < '+pt_max)
                    PlotYMax = h_proj.GetMaximum()*1.2
                    LineYMax = h_proj.GetMaximum()*0.8
                    canv = tdrCanvas(name+pt_ref, 0, 2, 0, PlotYMax,'Response', 'A.U.')
                    leg = tdrLeg(0.80,0.65,0.95,0.89, textSize=0.035)
                    tdrDraw(h_proj,'')
                    lines = {}
                    lines['fit1'] = rt.TLine(1-current_jer, 0.00001, 1-current_jer, LineYMax)
                    lines['fit2'] = rt.TLine(1+current_jer, 0.00001, 1+current_jer, LineYMax)
                    tdrDrawLine(lines['fit1'], lcolor=ROOT.kBlack, lstyle=ROOT.kDashed)
                    tdrDrawLine(lines['fit2'], lcolor=ROOT.kBlack, lstyle=ROOT.kDashed)
                    for perc in self.RMS:
                        funcs[perc].SetLineColor(colors[perc])
                        funcs[perc].Draw('same')
                        leg.AddEntry(funcs[perc], str(perc), 'l')
                        xq = xqs[perc]
                        yq = yqs[perc]
                        width = rms[perc]
                        for par in [0,1]:
                            lines[str(perc)+str(par)] = rt.TLine(yq[par], 0.00001, yq[par], LineYMax)
                            tdrDrawLine(lines[str(perc)+str(par)], lcolor=colors[perc])
                        lines[str(perc)+str(2)] = rt.TLine(1-width, 0.00001, 1-width, LineYMax)
                        lines[str(perc)+str(3)] = rt.TLine(1+width, 0.00001, 1+width, LineYMax)
                        tdrDrawLine(lines[str(perc)+str(2)], lcolor=colors[perc], lstyle=ROOT.kDashed)
                        tdrDrawLine(lines[str(perc)+str(3)], lcolor=colors[perc], lstyle=ROOT.kDashed)
                    ROOT.gStyle.SetOptFit(0)
                    canv.SaveAs(self.outdir+'/'+'_'.join(['JER_fit',year,algo,pt_ref,name])+'.pdf')
                    del canv
                    del h_proj

        for fname in f_:
            f_[fname].Close()
        f_['output'] = ROOT.TFile(self.outdir+'/MC_JER_mu_avg_'+year+'.root', 'RECREATE')
        for algo in self.algos:
            f_['output'].mkdir(algo)
            f_['output'].cd(algo)
            for eta_name in JER[algo].keys():
                x = JER[algo][eta_name]['pt']
                if len(x)==0: continue
                for mode in ['nominal']+['rms_'+'{:.2f}'.format(perc).replace('.','p') for perc in sorted(self.RMS)]:
                    y = JER[algo][eta_name]['JER_'+mode]
                    err = JER[algo][eta_name]['err_'+mode]
                    graph = ROOT.TGraphErrors(len(x), array('d',x), array('d',y), array('d',[0]*len(x)), array('d',err))
                    graph.Write('MCTruth_jer_mu_avg_MC_'+eta_name+'_'+mode)

        f_['output'].Write()
        f_['output'].Close()

def main():
    FitResponse_Ilias()

if __name__ == '__main__':
    main()
