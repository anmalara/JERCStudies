from Utils import *
from fileManipulation import *
from xml.dom.minidom import parseString


YearVars = {}
YearVars["JEC_Version"]         = {"2016":       "Summer16_07Aug2017_V11",
                                   "2017":       "Fall17_17Nov2017_V32",
                                   "2018":       "Autumn18_V19",
                                   "UL16APV":    "Summer19UL16APV_V8",
                                   "UL16nonAPV": "Summer19UL16_V8",
                                   # "UL17":       "Summer19UL16APV_V7",
                                   # "UL18":       "Summer19UL16APV_V7",
                                   }
YearVars["JER_Version"]         = {"2016":       "Summer16_25nsV1",
                                   "2017":       "Fall17_V3",
                                   "2018":       "Autumn18_V7b",
                                   "UL16APV":    "Summer20UL16APV_JRV3",
                                   "UL16nonAPV": "Summer20UL16APV_JRV3",
                                   "UL17":       "Summer20UL16APV_JRV3",
                                   "UL18":       "Summer20UL16APV_JRV3",
                                   }
YearVars["TargetLumi"]          = {"2016":       "36.33",
                                   "2017":       "41.48",
                                   "2018":       "59.83",
                                   "UL16APV":    "19.53",
                                   "UL16nonAPV": "16.80",
                                   "UL17":       "41.48",
                                   "UL18":       "59.83",
                                   }
YearVars["lumi_file"]           = {"2016": os.environ["CMSSW_BASE"]+"/src/UHH2/common/UHH2-data/2016/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.root",
                                   "2017": os.environ["CMSSW_BASE"]+"/src/UHH2/common/UHH2-data/2017/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.root",
                                   "2018": os.environ["CMSSW_BASE"]+"/src/UHH2/common/UHH2-data/2018/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.root",
                                   "UL16APV":    os.environ["CMSSW_BASE"]+"/src/UHH2/common/UHH2-data/2016/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.root",
                                   "UL16nonAPV": os.environ["CMSSW_BASE"]+"/src/UHH2/common/UHH2-data/2016/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.root",
                                   "UL17":       os.environ["CMSSW_BASE"]+"/src/UHH2/common/UHH2-data/2017/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.root",
                                   "UL18":       os.environ["CMSSW_BASE"]+"/src/UHH2/common/UHH2-data/2018/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.root",
                                   }



def newNumber(year,sample,ConfigFile):
    newNumber = 20
    if "DATA" in sample:
        if year=="2016":
            newNumber = 350
            if any(x in sample for x in ["DATA_SingleMuon_RunF", "DATA_SingleMuon_RunG", "DATA_SingleMuon_RunH"]):
                newNumber = 300
    return str(max(1,int(newNumber/(defaulTimePerJob/TimePerJob))))


@timeit
# def CreateConfigFiles(year, samples, all_samples, collections, channels, systematics, controls, SubmitDir, ConfigFile, Path_SFRAME, lumi):
def CreateConfigFiles(year, folder="test",ConfigFile_="JERCStudies"):
    global TimePerJob; global defaulTimePerJob; defaulTimePerJob = 3.
    TimePerJob = defaulTimePerJob
    inputdir = os.environ["CMSSW_BASE"]+"/src/UHH2/JERCStudies/config/"
    outdir = inputdir+"SubmittedJobs/"+folder+"/"
    a = os.system("mkdir -p "+outdir)
    a = os.system("mkdir -p /nfs/dust/cms/user/"+os.environ["USER"]+"/sframe_all/JERCStudies/"+folder)
    ConfigFile = inputdir+ConfigFile_+"Config.xml"
    samples = []
    with open(ConfigFile) as f_:
        for line in f_.readlines():
            if not "<InputData" in line: continue
            samples.append([x for x in line.split() if "Version=" in x][0].replace('Version="',"").replace('"',""))
    print samples
    for sample in samples:
        if not year in sample: continue
        filename = ConfigFile_+"Config_"+sample+".xml"
        a = os.system("cp "+ConfigFile+" "+outdir+filename)
        a = os.system("cp "+inputdir+"JobConfig.dtd "+outdir)
        comments = []
        for el in samples:
            if sample == el: continue
            comments.append(["<InputData", "Type", "MC",   '"'+el+'"'])
        comment_lines(outdir, filename, comments, remove=True)
        changes = []
        # Change anmalara when creating xml. change also email
        changes.append(["Mail=", "USER@mail.desy.de", "USER@mail.desy.de", os.environ["USER"]+"@mail.desy.de"])
        changes.append(["<!ENTITY", "/nfs/dust/cms/user/USER", "USER", os.environ["USER"]])
        changes.append(["<!ENTITY", "CMSSW_BASE", "CMSSW_BASE", os.environ["CMSSW_BASE"]])
        changes.append(["<!ENTITY", "STUDIES", "Studies", folder])
        change_lines(outdir, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
        changes = []
        changes.append(["<ConfigSGE", "Workdir", "workdir_"+ConfigFile_, "workdir_"+filename.replace(".xml","")])
        changes.append(["<ConfigParse", 'FileSplit="20"', 'FileSplit="20"', 'FileSplit="'+newNumber(year,sample,ConfigFile)+'"'])
        changes.append(["<!ENTITY", "OUTDIR", outdir , outdir+"/"+folder])
        changes.append(["<!ENTITY", "YEAR", 'defaultValue', year.replace("16APV", "16preVFP").replace("16nonAPV", "16postVFP")])
        changes.append(["<Cycle", "TargetLumi", 'defaultValue', YearVars["TargetLumi"][year]])
        for var in YearVars:
            changes.append(["<!ENTITY", var, 'defaultValue', YearVars[var][year]])
        change_lines(outdir, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])



if __name__ == '__main__':

    years = ["UL16APV","UL16nonAPV"]

    for year in years:
        CreateConfigFiles(year=year)
