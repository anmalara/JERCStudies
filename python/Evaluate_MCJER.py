import ROOT as rt
rt.gInterpreter.ProcessLine('#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"')
rt.gInterpreter.ProcessLine('#include "CondFormats/JetMETObjects/interface/JetResolution.h"')

def rhoFromMu(mu=0):
    # Eta_0.0-1.3, jt320
    return (1.01272 + 0.551183*mu + 0.000362936*mu*mu)

class Evaluate_MCJER():
    def __init__(self, JERfile="", JERSFfile=""):
        self.JERfile   = JERfile
        self.JERSFfile = JERSFfile
        self.eta = 1.0
        self.pt  = 300
        self.mu = 30
        self.rho = rhoFromMu(self.mu)
        self.SetParameters()

    def setPt(self, pt):
        self.pt = pt
        if hasattr(self, 'jp'): self.jp.setJetPt(self.pt)

    def setEta(self, eta):
        self.eta = eta
        if hasattr(self, 'jp'): self.jp.setJetEta(self.eta)

    def setMu(self, mu):
        self.mu = mu
        self.rho = rhoFromMu(self.mu)
        if hasattr(self, 'jp'): self.jp.setRho(self.rho)

    def setRho(self, rho):
        self.rho = rho
        if hasattr(self, 'jp'): self.jp.setRho(self.rho)

    def SetParameters(self):
        self.jp = rt.JME.JetParameters()
        self.jp.setJetPt(self.pt)
        self.jp.setJetEta(self.eta)
        self.jp.setRho(self.rho)

    def getSF(self, variation = "nominal"):
        var = rt.Variation.NOMINAL
        if variation=="down":
            var = rt.Variation.DOWN
        elif variation=="up":
            var = rt.Variation.UP
        sf = 1
        if self.JERSFfile!="":
            sf = rt.JME.JetResolutionScaleFactor(self.JERSFfile).getScaleFactor(self.jp, var)
        return sf

    def getResolution(self, isData=False):
        jer = rt.JME.JetResolution(self.JERfile).getResolution(self.jp)
        if self.JERSFfile!="" and isData:
            jer *= self.getSF()
        return jer

    def GetParameters(self):
        jer = rt.JME.JetResolution(self.JERfile)
        pars_ = jer.getResolutionObject().getRecord(self.jp).getParametersValues()
        pars = []
        for p in range(pars_.size()):
            pars.append(pars_[p])
        return pars


def main():
    #This is an example
    JERfile  = "../../Summer19UL18/MC_JER/Summer19UL18_V2_MC_PtResolution_ak4pfchsl1l2l3.txt"
    JERSFfile = "../../Summer19UL18/JER_SF/Summer19UL18_JRV1_MC_SF_AK4PFchs.txt"
    print Evaluate_MCJER(JERfile).getResolution()
    print Evaluate_MCJER(JERfile,JERSFfile).getResolution()


if __name__ == '__main__':
    main()
