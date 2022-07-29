import time, sys, os, glob
import functools, argparse

import math
import numpy as np
import pandas as pd
from array import array
from collections import OrderedDict

from ModuleRunnerBase import *

from JERCProtoLab.macros.common_info.common_binning import *
import JERCProtoLab.macros.plotting.tdrstyle_JERC as tdr
from JERCProtoLab.macros.plotting.tdrstyle_JERC import *
import ROOT
ROOT.gInterpreter.ProcessLine('#include "'+os.environ['CMSSW_BASE']+'/src/UHH2/JERCStudies/include/Utils.hpp"')

def RemoveRootLabels():
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gErrorIgnoreLevel = ROOT.kError

RemoveRootLabels()

def modify_printed_string(type,string): return '%s%s\033[0m'%(type,string)
def red(string):     return modify_printed_string('\x1b[0;31m',string)
def green(string):   return modify_printed_string('\x1b[0;32m',string)
def yellow(string):  return modify_printed_string('\x1b[0;33m',string)
def blue(string):    return modify_printed_string('\x1b[0;34m',string)
def magenta(string): return modify_printed_string('\x1b[0;35m',string)
def cyan(string):    return modify_printed_string('\x1b[0;36m',string)
def bold(string):    return modify_printed_string('\033[1m',string)

def debugStr(string, color=blue):
    print(color('  --> '+str(string)))

def prettydict(d, indent=8, color=blue):
    space = max([0]+[len(str(x)) for x in d])+2
    for key, value in d.items():
        print(color(' '*indent + str(key))),
        if isinstance(value, dict):
            print('')
            prettydict(value, len(' '*indent + str(key)+' '*(space+1-len(str(key)))))
        else:
            print(color(' '*(space-len(str(key))) + str(value)))

def timeit(method):
    @functools.wraps(method)
    def timed(*args, **kw):
        print('Start'+method.__name__)
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts))
        else:
            print('%r  %2.2f s' % (method.__name__, (te - ts)))
        return result
    return timed

def PrintFormattedLine(listArgs=[], space=10):
    for x in listArgs: print x+' '*(space-len(str(x)) if space-len(str(x))>0 else 2*space-len(str(x))),
    print('\t')


def FindIndexVector(vec, val):
    return vec.index(FilterVector(vec, val, '')[0]) if len(FilterVector(vec, val, ''))!=0 else len(vec)

def FilterVector(vec, val, invert):
    vec = np.array(vec)
    mask = vec > val if invert == 'invert' else vec < val
    return vec[~mask]

def FloatToString(x):
    return str(int(x))+'p'+str(int((x-int(x))*1000))


def Getfile(ver, mode):
    return GenericPath().Path_UHH2+'JRDatabase/textFiles/{ver}_MC/{ver}_MC_{mode}_AK4PFchs.txt'.format(ver = ver, mode=mode)

def GetJERfile(ver): return Getfile(ver, 'PtResolution')
def GetSFfile(ver): return Getfile(ver, 'SF')
def Oplus(x,y): return math.sqrt(x*x+y*y)
def Ominus(x,y): return math.sqrt(x*x-y*y) if x>=y else x
