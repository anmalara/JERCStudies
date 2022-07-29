from Utils import *

def TGraphFuncRatio(graph,func):
    x_vals = []
    y_vals = []
    x_errs = []
    y_errs = []
    for bin in range(graph.GetN()):
        x,y = (ROOT.Double(0),ROOT.Double(0))
        graph.GetPoint(int(bin),x,y)
        x_vals.append(x)
        y_vals.append(y/func.Eval(x))
        x_errs.append(graph.GetErrorX(bin))
        y_errs.append(y/func.Eval(x)*Oplus(graph.GetErrorY(bin)/y,graph.GetErrorY(bin)/func.Eval(x)))
    return ROOT.TGraphErrors(len(x_vals), array('d',x_vals), array('d',y_vals), array('d',x_errs), array('d',y_errs))

def TGraphRatio(graph1,graph2, doOffset=False):
    offset1, offset2 = (0,0)
    npoints_1 = graph1.GetN()
    npoints_2 = graph2.GetN()
    if npoints_1!=npoints_2:
        if doOffset:
            if npoints_1-npoints_2 == 6: offset1 = 6
            if npoints_1-npoints_2 == 5: offset1 = 5
            if npoints_1-npoints_2 == 4: offset1 = 3
            if npoints_1-npoints_2 == 3: offset1 = 3
            if npoints_1-npoints_2 == 2: offset1 = 3
            if npoints_1-npoints_2 == 1: offset1 = 3
        else:
            print 'Error in TGraphRatio: Not same points:', npoints_1, npoints_2
    x_vals = []
    y_vals = []
    x_errs = []
    y_errs = []
    for bin in range(npoints_1-offset1):
        x1,y1,x2,y2 = (ROOT.Double(0),ROOT.Double(0),ROOT.Double(0),ROOT.Double(0))
        graph1.GetPoint(int(bin+offset1),x1,y1)
        graph2.GetPoint(int(bin+offset2),x2,y2)
        if y2==0: continue
        if x1>3500: continue
        if x2>3500: continue
        if x1!=x2:
            if (1-x1/x2)*100<2:
                x1 = (x1+x2)/2
                x2 = x1
                y1 = graph1.Eval(x1)
                y2 = graph2.Eval(x2)
            else:
                print 'Error in TGraphRatio: Not same x,', x1,x2, (1-x1/x2)*100, bin+offset1, bin+offset2
        x_vals.append(x1)
        y_vals.append(y1/y2)
        x_errs.append(graph1.GetErrorX(bin+offset1))
        y_errs.append(y1/y2*Oplus(0 if y1==0 else graph1.GetErrorY(bin+offset1)/y1,graph2.GetErrorY(bin+offset2)/y2))
    return ROOT.TGraphErrors(len(x_vals), array('d',x_vals), array('d',y_vals), array('d',x_errs), array('d',y_errs))

def HistToGraph(hist, min_val = None, max_val = None, remove_values=None):
    x_vals = []
    y_vals = []
    x_errs = []
    y_errs = []
    for bin in range(hist.GetNbinsX()+1):
        x_val = hist.GetBinCenter(bin)
        y_val = hist.GetBinContent(bin)
        if y_val==0:continue
        if min_val and x_val<min_val: continue
        if max_val and x_val>max_val: continue
        if remove_values and x_val in remove_values: continue
        x_vals.append(x_val)
        y_vals.append(y_val)
        x_errs.append(hist.GetBinWidth(bin)/2)
        y_errs.append(hist.GetBinError(bin))
    return ROOT.TGraphErrors(len(x_vals), array('d',x_vals), array('d',y_vals), array('d',x_errs), array('d',y_errs))

def HistToGraphUncertainty(hist, variations, min_val = None, max_val = None, remove_values=None):
    x_vals = []
    y_vals = []
    x_errs = []
    y_errs = []
    for bin in range(hist.GetNbinsX()+1):
        x_val = hist.GetBinCenter(bin)
        y_val = hist.GetBinContent(bin)
        if y_val==0:continue
        if min_val and x_val<min_val: continue
        if max_val and x_val>max_val: continue
        if remove_values and x_val in remove_values: continue
        y_err = hist.GetBinError(bin)
        for name, variation in variations.items():
            err_up = variation['up'].GetBinContent(bin)-y_val
            err_down = variation['down'].GetBinContent(bin)-y_val
            if 'gaustails' in name:
                err_up /=2
                err_down /=2
            mean = (err_up+err_down)/2
            diff = (err_up-err_down)/2
            err = Oplus(mean, diff*math.sqrt(2))
            if err/y_val>0.85:continue
            y_val += diff
            y_err = Oplus(y_err, err)
            # if err/y_val<0.2:continue
            # print x_val, err/y_val, round(y_err,3), name, round(err_up,3), round(err_down,3), round(mean,3), round(diff,3), round(err,3)
        x_vals.append(x_val)
        y_vals.append(y_val)
        x_errs.append(10+hist.GetBinWidth(bin)/2)
        y_errs.append(y_err)
    return ROOT.TGraphErrors(len(x_vals), array('d',x_vals), array('d',y_vals), array('d',x_errs), array('d',y_errs))

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


def MergeGraphs(graph1, graph2):
    x_vals = []
    y_vals = []
    x_errs = []
    y_errs = []
    for bin in range(graph1.GetN()):
        x,y = (ROOT.Double(0),ROOT.Double(0))
        graph1.GetPoint(int(bin),x,y)
        x_vals.append(x)
        y_vals.append(y)
        x_errs.append(graph1.GetErrorX(bin))
        y_errs.append(graph1.GetErrorY(bin))
    if sorted(x_vals) != x_vals:
        print x_vals
        raise RuntimeError('Error in ExpandGraph. Graph 1 not sorted')
    for bin in range(graph2.GetN()):
        x,y = (ROOT.Double(0),ROOT.Double(0))
        graph2.GetPoint(int(bin),x,y)
        index = FindIndexVector(x_vals, x)
        x_vals.insert(index, x)
        y_vals.insert(index,y)
        x_errs.insert(index,graph2.GetErrorX(bin))
        y_errs.insert(index,graph2.GetErrorY(bin))
    return ROOT.TGraphErrors(len(x_vals), array('d',x_vals), array('d',y_vals), array('d',x_errs), array('d',y_errs))

def GetConfidenceIntervals(func, fitRes, Nbins, xcenters, ci, cl = 0.68):
    npar = func.GetNpar()
    covmatr = fitRes.GetCovarianceMatrix()
    tStudent = ROOT.TMath.StudentQuantile(0.5 + cl/2, func.GetNDF())
    chindf = ROOT.TMath.Sqrt(func.GetChisquare()/func.GetNDF())
    grad = array('d',[0.]*npar)
    sum_vector = array('d',[0.]*npar)
    for ipoint in range(Nbins):
        ci_=0
        func.GradientPar(array('d',[xcenters[ipoint]]), grad)
        # multiply the covariance matrix by gradient
        for irow in range(npar):
            sum_vector[irow]=0
            for icol in range(npar):
                sum_vector[irow] += covmatr[irow][icol]*grad[icol]
        for i_ in range(npar):
            ci_ += grad[i_]*sum_vector[i_]
        ci_ = ROOT.TMath.Sqrt(ci_)
        ci[ipoint] = ci_*tStudent*chindf

def ComputeHistWithCL(name, func, fitRes, graph, cl=0.68):
    x_vals = array('d',graph.GetX())
    x_errs = []
    y_vals_band = []
    y_vals_ratioband = []
    y_errs_band = []
    y_errs_ratioband = []
    Nbins = len(x_vals)
    ci = array('d',[0.]*Nbins)
    GetConfidenceIntervals(func, fitRes, Nbins, x_vals, ci, cl)
    for bin in range(Nbins):
        y_ = func.Eval(x_vals[bin])
        x_errs.append(graph.GetErrorX(bin))
        y_vals_band.append(y_)
        y_vals_ratioband.append(1)
        y_errs_band.append(ci[bin])
        y_errs_ratioband.append(ci[bin]/y_)
    band = ROOT.TGraphErrors(Nbins, array('d',x_vals), array('d',y_vals_band), array('d',x_errs), array('d',y_errs_band))
    ratioband = ROOT.TGraphErrors(Nbins, array('d',x_vals), array('d',y_vals_ratioband), array('d',x_errs), array('d',y_errs_ratioband))
    return band, ratioband
