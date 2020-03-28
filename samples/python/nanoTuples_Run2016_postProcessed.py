import copy, os, sys
from RootTools.core.Sample import Sample 
import ROOT

# Logging
import logging
logger = logging.getLogger(__name__)

# Data directory
try:
    directory_ = sys.modules['__main__'].directory_
except:
    from tWZ.samples.default_locations import default_locations
    directory_ = default_locations.data_2016

logger.info("Loading data samples from directory %s", directory_)

dirs = {}
for (run, version) in [('B','_ver2'),('C',''),('D',''),('E',''),('F',''),('G',''),('H','')]: # no event that passes json in B_ver1
    runTag = 'Run2016' + run + '_25Oct2019' + version
    dirs["DoubleEG_Run2016"         + run + version ] = ["DoubleEG_"          + runTag ]
    dirs["DoubleMuon_Run2016"       + run + version ] = ["DoubleMuon_"        + runTag ]
    dirs["SingleElectron_Run2016"   + run + version ] = ["SingleElectron_"    + runTag ]
    dirs["SingleMuon_Run2016"       + run + version ] = ["SingleMuon_"        + runTag ]
    dirs["MuonEG_Run2016"           + run + version ] = ["MuonEG_"            + runTag ]

def merge(pd, totalRunName, listOfRuns):
    dirs[pd + '_' + totalRunName] = []
    for run in listOfRuns: dirs[pd + '_' + totalRunName].extend(dirs[pd + '_' + run])

for pd in ['MuonEG', 'DoubleMuon', 'DoubleEG', 'SingleElectron', 'SingleMuon']:
    merge(pd, 'Run2016BCD',    ['Run2016B_ver2', 'Run2016C', 'Run2016D'])
    merge(pd, 'Run2016BCDEFG', ['Run2016BCD', 'Run2016E', 'Run2016F', 'Run2016G'])
    merge(pd, 'Run2016EF', ['Run2016E', 'Run2016F'])
    merge(pd, 'Run2016GH', ['Run2016G', 'Run2016H'])
    merge(pd, 'Run2016',       ['Run2016BCDEFG', 'Run2016H'])

for key in dirs:
    dirs[key] = [ os.path.join( directory_, dir) for dir in dirs[key]]

def getSample(pd, runName, lumi):
    sample      = Sample.fromDirectory(name=(pd + '_' + runName), treeName="Events", texName=(pd + ' (' + runName + ')'), directory=dirs[pd + '_' + runName])
    sample.lumi = lumi
    return sample

l_tot = (5.883+2.646+4.353+4.050+3.124+7.554+8.494+0.217) # old normtag

allSamples_Data25ns = []

DoubleEG_Run2016BCD                = getSample('DoubleEG',         'Run2016BCD',       (5.883+2.646+4.353)*(35.9/l_tot)*1000)
DoubleMuon_Run2016BCD              = getSample('DoubleMuon',       'Run2016BCD',       (5.883+2.646+4.353)*(35.9/l_tot)*1000)
SingleElectron_Run2016BCD          = getSample('SingleElectron',   'Run2016BCD',       (5.883+2.646+4.353)*(35.9/l_tot)*1000)
SingleMuon_Run2016BCD              = getSample('SingleMuon',       'Run2016BCD',       (5.883+2.646+4.353)*(35.9/l_tot)*1000)
MuonEG_Run2016BCD                  = getSample('MuonEG',           'Run2016BCD',       (5.883+2.646+4.353)*(35.9/l_tot)*1000)
allSamples_Data25ns += [MuonEG_Run2016BCD, DoubleEG_Run2016BCD, DoubleMuon_Run2016BCD, SingleElectron_Run2016BCD, SingleMuon_Run2016BCD]

DoubleEG_Run2016EF                = getSample('DoubleEG',         'Run2016EF',         (4.050+3.124)*(35.9/l_tot)*1000)
DoubleMuon_Run2016EF              = getSample('DoubleMuon',       'Run2016EF',         (4.050+3.124)*(35.9/l_tot)*1000)
SingleElectron_Run2016EF          = getSample('SingleElectron',   'Run2016EF',         (4.050+3.124)*(35.9/l_tot)*1000)
SingleMuon_Run2016EF              = getSample('SingleMuon',       'Run2016EF',         (4.050+3.124)*(35.9/l_tot)*1000)
MuonEG_Run2016EF                  = getSample('MuonEG',           'Run2016EF',         (4.050+3.124)*(35.9/l_tot)*1000)
allSamples_Data25ns += [MuonEG_Run2016EF, DoubleEG_Run2016EF, DoubleMuon_Run2016EF, SingleElectron_Run2016EF, SingleMuon_Run2016EF]

DoubleEG_Run2016GH                = getSample('DoubleEG',         'Run2016GH',         (7.554+8.494+0.217)*(35.9/l_tot)*1000)
DoubleMuon_Run2016GH              = getSample('DoubleMuon',       'Run2016GH',         (7.554+8.494+0.217)*(35.9/l_tot)*1000)
SingleElectron_Run2016GH          = getSample('SingleElectron',   'Run2016GH',         (7.554+8.494+0.217)*(35.9/l_tot)*1000)
SingleMuon_Run2016GH              = getSample('SingleMuon',       'Run2016GH',         (7.554+8.494+0.217)*(35.9/l_tot)*1000)
MuonEG_Run2016GH                  = getSample('MuonEG',           'Run2016GH',         (7.554+8.494+0.217)*(35.9/l_tot)*1000)
allSamples_Data25ns += [MuonEG_Run2016GH, DoubleEG_Run2016GH, DoubleMuon_Run2016GH, SingleElectron_Run2016GH, SingleMuon_Run2016GH]

DoubleEG_Run2016                  = getSample('DoubleEG',         'Run2016',           (35.9)*1000)
DoubleMuon_Run2016                = getSample('DoubleMuon',       'Run2016',           (35.9)*1000)
SingleElectron_Run2016            = getSample('SingleElectron',   'Run2016',           (35.9)*1000)
SingleMuon_Run2016                = getSample('SingleMuon',       'Run2016',           (35.9)*1000)
MuonEG_Run2016                    = getSample('MuonEG',           'Run2016',           (35.9)*1000)
allSamples_Data25ns += [MuonEG_Run2016, DoubleEG_Run2016, DoubleMuon_Run2016, SingleElectron_Run2016, SingleMuon_Run2016]

Run2016 = Sample.combine("Run2016", [MuonEG_Run2016, DoubleEG_Run2016, DoubleMuon_Run2016, SingleElectron_Run2016, SingleMuon_Run2016], texName = "Run2016")
Run2016.lumi = (35.9)*1000
allSamples_Data25ns.append(Run2016)

Run2016BCD = Sample.combine("Run2016BCD", [MuonEG_Run2016BCD, DoubleEG_Run2016BCD, DoubleMuon_Run2016BCD, SingleElectron_Run2016BCD, SingleMuon_Run2016BCD], texName = "Run2016BCD")
Run2016BCD.lumi = (35.9)*1000
allSamples_Data25ns.append(Run2016BCD)

Run2016EF = Sample.combine("Run2016EF", [MuonEG_Run2016EF, DoubleEG_Run2016EF, DoubleMuon_Run2016EF, SingleElectron_Run2016EF, SingleMuon_Run2016EF], texName = "Run2016EF")
Run2016EF.lumi = (35.9)*1000
allSamples_Data25ns.append(Run2016EF)

Run2016GH = Sample.combine("Run2016GH", [MuonEG_Run2016GH, DoubleEG_Run2016GH, DoubleMuon_Run2016GH, SingleElectron_Run2016GH, SingleMuon_Run2016GH], texName = "Run2016GH")
Run2016GH.lumi = (35.9)*1000
allSamples_Data25ns.append(Run2016GH)

for s in allSamples_Data25ns:
  s.color   = ROOT.kBlack
  s.isData  = True
