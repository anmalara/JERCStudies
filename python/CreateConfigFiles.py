#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Utils import *
from fileManipulation import *


YearVars = {}
YearVars['JEC_Version']         = {'2016':       'Summer16_07Aug2017_V11',
                                   '2017':       'Fall17_17Nov2017_V32',
                                   '2018':       'Autumn18_V19',
                                   'UL16APV':    'Summer19UL16APV_V9',
                                   'UL16nonAPV': 'Summer19UL16_V9',
                                   'UL17':       'Summer19UL17_V5',
                                   'UL18':       'Summer19UL18_V5',
                                   }
YearVars['JER_Version']         = {'2016':       'Summer16_25nsV1',
                                   '2017':       'Fall17_V3',
                                   '2018':       'Autumn18_V7b',
                                   'UL16APV':    'Summer20UL16APV_JRV3',
                                   'UL16nonAPV': 'Summer20UL16APV_JRV3',
                                   'UL17':       'Summer20UL16APV_JRV2',
                                   'UL18':       'Summer20UL16APV_JRV2',
                                   }
YearVars['TargetLumi']          = {'2016':       '36.33',
                                   '2017':       '41.48',
                                   '2018':       '59.83',
                                   'UL16APV':    str(tdr.commonScheme['lumidec']['UL16APV']),
                                   'UL16nonAPV': str(tdr.commonScheme['lumidec']['UL16nonAPV']),
                                   'UL17':       str(tdr.commonScheme['lumidec']['UL17']),
                                   'UL18':       str(tdr.commonScheme['lumidec']['UL18']),
                                   }
YearVars['lumi_file']           = {'2016': os.environ['CMSSW_BASE']+'/src/UHH2/common/UHH2-data/2016/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.root',
                                   '2017': os.environ['CMSSW_BASE']+'/src/UHH2/common/UHH2-data/2017/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.root',
                                   '2018': os.environ['CMSSW_BASE']+'/src/UHH2/common/UHH2-data/2018/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.root',
                                   'UL16APV':    os.environ['CMSSW_BASE']+'/src/UHH2/common/UHH2-data/2016/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.root',
                                   'UL16nonAPV': os.environ['CMSSW_BASE']+'/src/UHH2/common/UHH2-data/2016/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.root',
                                   'UL17':       os.environ['CMSSW_BASE']+'/src/UHH2/common/UHH2-data/2017/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.root',
                                   'UL18':       os.environ['CMSSW_BASE']+'/src/UHH2/common/UHH2-data/2018/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.root',
                                   }



def newNumber(year,sample,ConfigFile):
    newNumber = 20
    if 'DATA' in sample:
        if year=='2016':
            newNumber = 350
            if any(x in sample for x in ['DATA_SingleMuon_RunF', 'DATA_SingleMuon_RunG', 'DATA_SingleMuon_RunH']):
                newNumber = 300
    return str(max(1,int(newNumber/(defaulTimePerJob/TimePerJob))))

def CreateConfigFiles(year, folder='test', config='JERCStudies', command ='c'):
    global TimePerJob; global defaulTimePerJob; defaulTimePerJob = 3.
    TimePerJob = defaulTimePerJob
    do_create = command == 'c'
    inputdir = os.environ['CMSSW_BASE']+'/src/UHH2/JERCStudies/config/'
    outdir = inputdir+'SubmittedJobs/'+args.folder+'/'
    if do_create:
        a = os.system('mkdir -p '+outdir)
        a = os.system('mkdir -p /nfs/dust/cms/user/'+os.environ['USER']+'/sframe_all/JERCStudies/'+folder)
    ConfigFile = inputdir+config+'Config.xml'
    samples = []
    with open(ConfigFile) as f_:
        for line in f_.readlines():
            if not '<InputData' in line: continue
            samples.append([x for x in line.split() if 'Version=' in x][0].replace('Version=','').replace('"',''))
    commands = []
    for sample in samples:
        if not year in sample: continue
        filename = config+'Config_'+sample+'.xml'
        print(blue('--> Running on: '+filename))
        if do_create:
            a = os.system('cp '+ConfigFile+' '+outdir+filename)
            a = os.system('cp '+inputdir+'JobConfig.dtd '+outdir)
        comments = []
        for el in samples:
            if sample == el: continue
            comments.append(['<InputData', 'Type', 'MC',   '"'+el+'"'])
        if do_create:
            comment_lines(outdir, filename, comments, remove=True)
        changes = []
        # Change anmalara when creating xml. change also email
        changes.append(['Mail=', 'USER@mail.desy.de', 'USER@mail.desy.de', os.environ['USER']+'@mail.desy.de'])
        changes.append(['<!ENTITY', '/nfs/dust/cms/user/USER', 'USER', os.environ['USER']])
        changes.append(['<!ENTITY', 'CMSSW_BASE', 'CMSSW_BASE', os.environ['CMSSW_BASE']])
        changes.append(['<!ENTITY', 'STUDIES', 'Studies', folder])
        if do_create:
            change_lines(outdir, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
        changes = []
        changes.append(['<ConfigSGE', 'Workdir', 'workdir_'+config, 'workdir_'+filename.replace('.xml','')])
        changes.append(['<ConfigParse', 'FileSplit="10"', 'FileSplit="10"', 'FileSplit="'+newNumber(year,sample,ConfigFile)+'"'])
        changes.append(['<!ENTITY', 'OUTDIR', outdir , outdir+'/'+folder])
        changes.append(['<!ENTITY', 'YEAR', 'defaultValue', year.replace('16APV', '16preVFP').replace('16nonAPV', '16postVFP')])
        changes.append(['<Cycle', 'TargetLumi', 'defaultValue', YearVars['TargetLumi'][year]])
        for var in YearVars:
            changes.append(['<!ENTITY', var, 'defaultValue', YearVars[var][year]])
        if do_create:
            change_lines(outdir, filename, [el[0:2] for el in changes ], [el[2:3] for el in changes ], [el[3:4] for el in changes ])
        if not do_create:
            command_list = ['sframe_batch.py', outdir+filename]
            if command!='':
                command_list.append('-'+command)
            process = subprocess.Popen(command_list,  cwd=outdir)
            process.wait()



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--folder',  '-f', action='store', default='JERStudies',  dest='folder',  type=str)
    parser.add_argument('--config',  '-x', action='store', default='JERCStudies', dest='config',  type=str)
    parser.add_argument('--year',    '-y', action='store', default='UL18',        dest='year',    type=str)
    parser.add_argument('--command', '-c', action='store', default='',            dest='command', type=str)

    args = parser.parse_args()

    CreateConfigFiles(year = args.year, folder = args.folder, config = args.config,  command=args.command)
