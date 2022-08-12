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


class FitResponse_Segvi(InputBase):
    def __init__(self):
        InputBase.__init__(self, inputdir = '')
        self.inputdir = os.path.join(self.Path_ANALYSIS,'StorageArea','JERCombination')
        self.year = 'UL17'
        self.pt_bins = [20, 23, 27, 27, 30, 35, 40, 45, 57, 72, 90, 120, 150, 200, 300, 400, 550, 750, 1000, 1500, 2000, 2500, 3000]
        self.RMS = [0.99,0.95,0.87,0.68]
        self.LoadFiles()


    def LoadFiles(self):
        f_ = ROOT.TFile(self.inputdir+'/output-MC-mpf_mpfx_2Dand3D.root')
        eta_folders = ['Eta_'+eta for eta in ['0.0-0.5','0.5-1.0','1.0-1.5','1.5-2.0','2.0-2.5','2.5-3.0','3.0-3.5','3.5-4.0','4.0-4.5']]
        histos = {}
        JER = OrderedDict()
        for eta_folder in eta_folders:
            eta_min, eta_max = eta_folder.replace('Eta_','').split('-')
            # histos[eta_folder] = {} 'h2mpf', 'h2mpfx', 'h2mpf_tag', 'h2mpfx_tag',
            for hname in ['h2mpf']:
                h = f_.Get('Standard/'+eta_folder+'/h2mpf')
                for pt_bin in range(len(self.pt_bins)-1):
                    pt_min, pt_max = (self.pt_bins[pt_bin], self.pt_bins[pt_bin+1])
                    if pt_min >3000: continue
                    if float(pt_max)*math.cosh(float(eta_min))>6500: continue
                    pt_avg = (pt_min+pt_max)/2
                    pt_range = str(pt_min)+'_'+str(pt_max)
                    h_proj_name = eta_folder+hname+str(pt_avg)
                    if not h_proj_name in JER: JER[h_proj_name] = {}
                    bin_min = GetBin(h=h, x_val=pt_min, doMin=True)
                    bin_max = GetBin(h=h, x_val=pt_max, doMin=False)
                    x_min = h.GetXaxis().GetBinCenter(bin_min)-h.GetXaxis().GetBinWidth(bin_min)/2
                    x_max = h.GetXaxis().GetBinCenter(bin_max)+h.GetXaxis().GetBinWidth(bin_max)/2
                    h_proj = h.ProjectionY(h_proj_name, bin_min, bin_max, 'e')
                    print h_proj.GetXaxis().GetXmin(), h_proj.GetXaxis().GetXmax()
                    if h_proj.GetEntries()<10:
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
                        JER[h_proj_name].setdefault('JER_rms_'+'{:.2f}'.format(perc).replace('.','p'), []).append(funcs[perc].GetParameter(2))
                        JER[h_proj_name].setdefault('err_rms_'+'{:.2f}'.format(perc).replace('.','p'), []).append(funcs[perc].GetParError(2)+funcs[perc].GetParameter(2)*0.05)
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
                    JER[h_proj_name].setdefault('pt', []).append(pt_avg)
                    JER[h_proj_name].setdefault('JER_nominal', []).append(current_jer)
                    JER[h_proj_name].setdefault('err_nominal', []).append(gauss.GetParError(2)+current_jer*0.05)
                    tdr.cms_lumi_TeV = tdr.commonScheme['legend'][self.year]+' Legacy, '+commonScheme['lumi'][self.year]+' fb^{-1}'
                    tdr.extraText3 = []
                    tdr.extraText3.append(eta_min+' < |#eta| < '+eta_max)
                    tdr.extraText3.append(str(pt_min)+' < p_{T} < '+str(pt_max))
                    PlotYMax = h_proj.GetMaximum()*1.2
                    LineYMax = h_proj.GetMaximum()*0.8
                    canv = tdrCanvas(h_proj_name, 0, 2, 0, PlotYMax,'Response', 'A.U.')
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
                    canv.SaveAs(self.outdir+'/'+'_'.join(['JER_fit',eta_folder,hname,str(pt_range)])+'.pdf')
                    del canv
                    del h_proj
            # for h in histos[eta_folder].values():
            #     h.SetDirectory(0)
        f_.Close()


        # for pt_ref in pt_bins:
        #     print pt_ref
        #     pt_min, pt_max = pt_ref.replace('RefPt','').split('to')
        #     pt_avg = (float(pt_min)+float(pt_max))/2
        #     for algo in self.algos:
        #         h = f_[algo].Get(algo+'/RelRspVsJetEta_'+pt_ref)
        #         for eta_bin, etaRef in enumerate(self.etaBinsCommon):
        #             eta_min, eta_max = GetEtaMinMax(etaRef)
        #             if float(pt_max)*math.cosh(eta_min)>6500: continue
        #             name = eta_name = GetEtaName(etaRef)
        #             if not name in JER[algo]: JER[algo][name] = {}
        #             bin_min = GetBin(h=h, x_val=eta_min, doMin=True)
        #             bin_max = GetBin(h=h, x_val=eta_max, doMin=False)
        #             # x_min = h.GetXaxis().GetBinCenter(bin_min)-h.GetXaxis().GetBinWidth(bin_min)/2
        #             # x_max = h.GetXaxis().GetBinCenter(bin_max)+h.GetXaxis().GetBinWidth(bin_max)/2
        #             # print round(eta_min-x_min,2), round(eta_max-x_max,2),
        #             h_proj = h.ProjectionY(name+'pos', bin_min, bin_max, 'e')
        #             bin_min = GetBin(h=h, x_val=-eta_max, doMin=True)
        #             bin_max = GetBin(h=h, x_val=-eta_min, doMin=False)
        #             # x_min = h.GetXaxis().GetBinCenter(bin_min)-h.GetXaxis().GetBinWidth(bin_min)/2
        #             # x_max = h.GetXaxis().GetBinCenter(bin_max)+h.GetXaxis().GetBinWidth(bin_max)/2
        #             # print round(eta_max+x_min,2), round(eta_min+x_max,2)
        #             h_proj_neg = h.ProjectionY(name+'neg', bin_min, bin_max, 'e')
        #             h_proj.Add(h_proj_neg)
        #             if h_proj.GetEntries()<10:
        #                 print 'skipping', pt_avg, algo, name, pt_max, eta_min, float(pt_max)*math.cosh(eta_min)
        #                 continue
        #             x_min, x_max = (h_proj.GetMean()-h_proj.GetRMS(),h_proj.GetMean()+h_proj.GetRMS())
        #             gauss = ROOT.TF1('gauss','gaus',x_min,x_max)
        #             funcs = {}
        #             lines = {}
        #             xqs = {}
        #             yqs = {}
        #             rms = {}
        #             for perc in self.RMS:
        #                 h_proj_trunc = h_proj.Clone('h_proj_trunc_'+str(perc))
        #                 xq = array('d', [(1.00-perc)/2,(1.00+perc)/2])
        #                 yq = array('d', [0,0])
        #                 h_proj_trunc.GetQuantiles(2, yq, xq)
        #                 h_proj_trunc.GetXaxis().SetRange(h_proj_trunc.FindBin(yq[0]), h_proj_trunc.FindBin(yq[1]));
        #                 funcs[perc] = ROOT.TF1('gauss'+str(perc),'gaus', yq[0], yq[1])
        #                 fit_res = h_proj_trunc.Fit(funcs[perc], 'RQMS+', '')
        #                 JER[algo][name].setdefault('JER_rms_'+'{:.2f}'.format(perc).replace('.','p'), []).append(funcs[perc].GetParameter(2))
        #                 JER[algo][name].setdefault('err_rms_'+'{:.2f}'.format(perc).replace('.','p'), []).append(funcs[perc].GetParError(2)+funcs[perc].GetParameter(2)*0.05)
        #                 xqs[perc] = xq
        #                 yqs[perc] = yq
        #                 rms[perc] = funcs[perc].GetParameter(2)
        #             iterations = 1000
        #             last_jer = 100
        #             min_diff = 10
        #             counter = 0
        #             for i_ in range(iterations+1):
        #                 fit_res = h_proj.Fit(gauss, 'QMS+', '', x_min, x_max)
        #                 current_jer, chi2 = (gauss.GetParameter(2), gauss.GetChisquare())
        #                 diff = math.fabs(last_jer-current_jer)
        #                 if diff<min_diff: min_diff, counter = (diff,0)
        #                 else: counter += 1
        #                 if counter>10: break
        #                 if diff<1e-05: break
        #                 if i_!=iterations:
        #                     h_proj.GetListOfFunctions().RemoveLast()
        #                 last_jer, x_min, x_max = (current_jer, h_proj.GetMean()-current_jer,h_proj.GetMean()+current_jer)
        #             JER[algo][name].setdefault('pt', []).append(pt_avg)
        #             JER[algo][name].setdefault('JER_nominal', []).append(current_jer)
        #             JER[algo][name].setdefault('err_nominal', []).append(gauss.GetParError(2)+current_jer*0.05)
        #             tdr.cms_lumi_TeV = tdr.commonScheme['legend'][year]+' Legacy, '+commonScheme['lumi'][year]+' fb^{-1}'
        #             tdr.extraText3 = []
        #             tdr.extraText3.append(name.replace('p','.').split('_')[0]+' < |#eta| < '+name.replace('p','.').split('_')[1])
        #             tdr.extraText3.append(pt_min+' < p_{T} < '+pt_max)
        #             PlotYMax = h_proj.GetMaximum()*1.2
        #             LineYMax = h_proj.GetMaximum()*0.8
        #             canv = tdrCanvas(name+pt_ref, 0, 2, 0, PlotYMax,'Response', 'A.U.')
        #             leg = tdrLeg(0.80,0.65,0.95,0.89, textSize=0.035)
        #             tdrDraw(h_proj,'')
        #             lines = {}
        #             lines['fit1'] = rt.TLine(1-current_jer, 0.00001, 1-current_jer, LineYMax)
        #             lines['fit2'] = rt.TLine(1+current_jer, 0.00001, 1+current_jer, LineYMax)
        #             tdrDrawLine(lines['fit1'], lcolor=ROOT.kBlack, lstyle=ROOT.kDashed)
        #             tdrDrawLine(lines['fit2'], lcolor=ROOT.kBlack, lstyle=ROOT.kDashed)
        #             for perc in self.RMS:
        #                 funcs[perc].SetLineColor(colors[perc])
        #                 funcs[perc].Draw('same')
        #                 leg.AddEntry(funcs[perc], str(perc), 'l')
        #                 xq = xqs[perc]
        #                 yq = yqs[perc]
        #                 width = rms[perc]
        #                 for par in [0,1]:
        #                     lines[str(perc)+str(par)] = rt.TLine(yq[par], 0.00001, yq[par], LineYMax)
        #                     tdrDrawLine(lines[str(perc)+str(par)], lcolor=colors[perc])
        #                 lines[str(perc)+str(2)] = rt.TLine(1-width, 0.00001, 1-width, LineYMax)
        #                 lines[str(perc)+str(3)] = rt.TLine(1+width, 0.00001, 1+width, LineYMax)
        #                 tdrDrawLine(lines[str(perc)+str(2)], lcolor=colors[perc], lstyle=ROOT.kDashed)
        #                 tdrDrawLine(lines[str(perc)+str(3)], lcolor=colors[perc], lstyle=ROOT.kDashed)
        #             ROOT.gStyle.SetOptFit(0)
        #             canv.SaveAs(self.outdir+'/'+'_'.join(['JER_fit',year,algo,pt_ref,name])+'.pdf')
        #             del canv
        #             del h_proj
        #
        # for fname in f_:
        #     f_[fname].Close()
        # f_['output'] = ROOT.TFile(self.outdir+'/MC_JER_mu_avg_'+year+'.root', 'RECREATE')
        # for algo in self.algos:
        #     f_['output'].mkdir(algo)
        #     f_['output'].cd(algo)
        #     for eta_name in JER[algo].keys():
        #         x = JER[algo][eta_name]['pt']
        #         if len(x)==0: continue
        #         for mode in ['nominal']+['rms_'+'{:.2f}'.format(perc).replace('.','p') for perc in sorted(self.RMS)]:
        #             y = JER[algo][eta_name]['JER_'+mode]
        #             err = JER[algo][eta_name]['err_'+mode]
        #             graph = ROOT.TGraphErrors(len(x), array('d',x), array('d',y), array('d',[0]*len(x)), array('d',err))
        #             graph.Write('MCTruth_jer_mu_avg_MC_'+eta_name+'_'+mode)
        #
        # f_['output'].Write()
        # f_['output'].Close()

def main():
    FitResponse_Segvi()

if __name__ == '__main__':
    main()
