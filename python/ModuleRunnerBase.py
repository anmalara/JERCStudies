import os
import sys
import itertools


class Constants():
    def __init__(self):
        self.etaBinsEdges  = [0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.566, 1.740, 1.930, 2.043, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191]
        self.etaBinsEdges  = [0.000, 0.261]
        self.etaBinsCommon = [round(self.GetEtaBinCenter(x),4) for x in self.etaBinsEdges[:-1]]

    def GetEtaBinEdgeMin(self,val):
        return filter(lambda x: x<=val, self.etaBinsEdges)[-1]

    def GetEtaBinEdgeMax(self,val):
        return filter(lambda x: x>val, self.etaBinsEdges)[0]

    def GetEtaBinCenter(self,val):
        min = self.GetEtaBinEdgeMin(val)
        max = self.GetEtaBinEdgeMax(val)
        return (max+min)/2

    def GetEtaBinWidth(self,val):
        min = self.GetEtaBinEdgeMin(val)
        max = self.GetEtaBinEdgeMax(val)
        return (max-min)/2



class GenericPath:
    ''' Class container for paths '''
    def __init__(self):
        self.user = os.environ["USER"]
        self.cmssw_base = os.environ["CMSSW_BASE"]
        self.Path_UHH2  = self.cmssw_base+"/src/UHH2/"
        self.Path_ANALYSIS  = self.Path_UHH2+"JERCStudies/"
        self.Path_STORAGE   = self.Path_ANALYSIS+"StorageArea/"
        self.JERversions = {
            "UL16APV":    "Summer20UL16APV_JRV3",
            "UL16nonAPV": "Summer20UL16_JRV3",
            "UL17":       "Summer19UL17_JRV2",
            "UL18":       "Summer19UL18_JRV2",
            "EOY16":      "Summer16_25nsV1",
            "EOY17":      "Fall17_V3",
            "EOY18":      "Autumn18_V7b",
            "8TeV":       "Fall15_25nsV2",
            }
    def Get(self, source):
        return getattr(self,source)


class VariablesBase(GenericPath,Constants):
    ''' Class container for list of objects '''
    def __init__(self):
        GenericPath.__init__(self)
        Constants.__init__(self)
        self.PrefixrootFile     = "uhh2.AnalysisModuleRunner."
        self.RunPeriods_Dict    = {"2016":       ["B", "C", "D", "E", "F", "G", "H"],
                                   "2017":       ["B", "C", "D", "E", "F"],
                                   "2018":       ["A", "B", "C", "D"],
                                   "UL16APV":    ["B", "C", "D", "E", "F"],
                                   "UL16nonAPV": ["F", "G", "H"],
                                   "UL17":       ["B", "C", "D", "E", "F"],
                                   "UL18":       ["A", "B", "C", "D"],
                                   }
        self.years = {}
        self.years["all"] = sorted(self.RunPeriods_Dict.keys())
        self.years["EOY"] = list(filter(lambda x: not "UL" in x,self.years["all"]))
        self.years["UL"] = list(filter(lambda x: "UL" in x,self.years["all"]))


class ModuleRunnerBase(VariablesBase):
    ''' Class container for list of objects for particular year '''
    def __init__(self,year="2016", ):
        VariablesBase.__init__(self)
        self.year = year
        self.Path_STORAGE = self.Path_STORAGE+self.year+"/"
        self.ConfigDir    = self.Path_ANALYSIS+"config/"
        self.SubmitDir    = self.ConfigDir+"SubmittedJobs/"+self.year+"/"
