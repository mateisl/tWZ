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
    import tWZ.samples.nanoAODv6_private_locations as locations
    directory_ = locations.data_2017

logger.info("Loading data samples from directory %s", directory_)

dirs = {}
for (run, version) in [('B',''),('C',''),('D',''),('E',''),('F','')]:
    runTag = 'Run2017' + run + '_25Oct2019' + version
    dirs["DoubleEG_Run2017"         + run + version ] = ["DoubleEG_"          + runTag ]
    dirs["DoubleMuon_Run2017"       + run + version ] = ["DoubleMuon_"        + runTag ]
    dirs["MuonEG_Run2017"           + run + version ] = ["MuonEG_"            + runTag ]
    dirs["SingleMuon_Run2017"       + run + version ] = ["SingleMuon_"        + runTag ]
    dirs["SingleElectron_Run2017"   + run + version ] = ["SingleElectron_"    + runTag ]

def merge(pd, totalRunName, listOfRuns):
    dirs[pd + '_' + totalRunName] = []
    for run in listOfRuns: dirs[pd + '_' + totalRunName].extend(dirs[pd + '_' + run])

for pd in ['MuonEG', 'DoubleMuon',  'SingleElectron', 'SingleMuon']:  #'DoubleEG',
    merge(pd, 'Run2017',    ['Run2017B', 'Run2017D', 'Run2017C', 'Run2017E', 'Run2017F']) 
    merge(pd, 'Run2017CDE', ['Run2017C', 'Run2017D', 'Run2017E'])  
    merge(pd, 'Run2017BCDE',['Run2017B', 'Run2017D', 'Run2017C', 'Run2017E'])

for pd in ['DoubleEG']:
    merge(pd, 'Run2017',    ['Run2017B', 'Run2017C', 'Run2017E', 'Run2017F']) #'Run2017D',
    merge(pd, 'Run2017CDE', ['Run2017C', 'Run2017E']) #'Run2017D', 
    merge(pd, 'Run2017BCDE',['Run2017B', 'Run2017C', 'Run2017E']) #'Run2017D',

for key in dirs:
    dirs[key] = [ os.path.join( directory_, dir) for dir in dirs[key]]

def getSample(pd, runName, lumi):
    sample      = Sample.fromDirectory(name=(pd + '_' + runName), treeName="Events", texName=(pd + ' (' + runName + ')'), directory=dirs[pd + '_' + runName])
    sample.lumi = lumi
    return sample

allSamples_Data25ns = []

DoubleEG_Run2017                = getSample('DoubleEG',         'Run2017',       (41.5)*1000)
DoubleEG_Run2017B               = getSample('DoubleEG',         'Run2017B',      (4.823)*1000)
DoubleEG_Run2017CDE             = getSample('DoubleEG',         'Run2017CDE',    (9.664+4.252+9.278)*1000)
DoubleEG_Run2017BCDE            = getSample('DoubleEG',         'Run2017BCDE',   (4.823+9.664+4.252+9.278)*1000)
DoubleEG_Run2017F               = getSample('DoubleEG',         'Run2017F',      (13.540)*1000)
allSamples_Data25ns += [DoubleEG_Run2017, DoubleEG_Run2017B, DoubleEG_Run2017CDE, DoubleEG_Run2017BCDE, DoubleEG_Run2017F]

DoubleMuon_Run2017              = getSample('DoubleMuon',       'Run2017',       (41.5)*1000)
DoubleMuon_Run2017B             = getSample('DoubleMuon',       'Run2017B',      (4.823)*1000)
DoubleMuon_Run2017CDE           = getSample('DoubleMuon',       'Run2017CDE',    (9.664+4.252+9.278)*1000)
DoubleMuon_Run2017BCDE          = getSample('DoubleMuon',       'Run2017BCDE',   (4.823+9.664+4.252+9.278)*1000)
DoubleMuon_Run2017F             = getSample('DoubleMuon',       'Run2017F',      (13.540)*1000)
allSamples_Data25ns += [DoubleMuon_Run2017, DoubleMuon_Run2017B, DoubleMuon_Run2017CDE, DoubleMuon_Run2017BCDE, DoubleMuon_Run2017F]

MuonEG_Run2017                  = getSample('MuonEG',           'Run2017',       (41.5)*1000)
MuonEG_Run2017B                 = getSample('MuonEG',           'Run2017B',      (4.823)*1000)
MuonEG_Run2017CDE               = getSample('MuonEG',           'Run2017CDE',    (9.664+4.252+9.278)*1000)
MuonEG_Run2017BCDE              = getSample('MuonEG',           'Run2017BCDE',   (4.823+9.664+4.252+9.278)*1000)
MuonEG_Run2017F                 = getSample('MuonEG',           'Run2017F',      (13.540)*1000)
allSamples_Data25ns += [MuonEG_Run2017, MuonEG_Run2017B, MuonEG_Run2017CDE, MuonEG_Run2017BCDE, MuonEG_Run2017F]

SingleMuon_Run2017              = getSample('SingleMuon',       'Run2017',       (41.5)*1000)
SingleMuon_Run2017B             = getSample('SingleMuon',       'Run2017B',      (4.823)*1000)
SingleMuon_Run2017CDE           = getSample('SingleMuon',       'Run2017CDE',    (9.664+4.252+9.278)*1000)
SingleMuon_Run2017BCDE          = getSample('SingleMuon',       'Run2017BCDE',   (4.823+9.664+4.252+9.278)*1000)
SingleMuon_Run2017F             = getSample('SingleMuon',       'Run2017F',      (13.540)*1000)
allSamples_Data25ns += [SingleMuon_Run2017, SingleMuon_Run2017B, SingleMuon_Run2017CDE, SingleMuon_Run2017BCDE, SingleMuon_Run2017F]

SingleElectron_Run2017          = getSample('SingleElectron',   'Run2017',       (41.5)*1000)
SingleElectron_Run2017B         = getSample('SingleElectron',   'Run2017B',      (4.823)*1000)
SingleElectron_Run2017CDE       = getSample('SingleElectron',   'Run2017CDE',    (9.664+4.252+9.278)*1000)
SingleElectron_Run2017BCDE      = getSample('SingleElectron',   'Run2017BCDE',   (4.823+9.664+4.252+9.278)*1000)
SingleElectron_Run2017F         = getSample('SingleElectron',   'Run2017F',      (13.540)*1000)
allSamples_Data25ns += [SingleElectron_Run2017, SingleElectron_Run2017B, SingleElectron_Run2017CDE, SingleElectron_Run2017BCDE, SingleElectron_Run2017F]

Run2017 = Sample.combine("Run2017", [MuonEG_Run2017, DoubleEG_Run2017, DoubleMuon_Run2017, SingleMuon_Run2017, SingleMuon_Run2017], texName = "Run2017")
Run2017.lumi  = (41.5)*1000
Run2017B = Sample.combine("Run2017B", [MuonEG_Run2017B, DoubleEG_Run2017B, DoubleMuon_Run2017B, SingleElectron_Run2017B, SingleMuon_Run2017B], texName = "Run2017B")
Run2017B.lumi = (4.823)*1000
Run2017CDE = Sample.combine("Run2017CDE", [MuonEG_Run2017CDE, DoubleEG_Run2017CDE, DoubleMuon_Run2017CDE, SingleElectron_Run2017CDE, SingleMuon_Run2017CDE], texName = "Run2017CDE")
Run2017CDE.lumi = (9.664+4.252+9.278)*1000
Run2017BCDE = Sample.combine("Run2017BCDE", [MuonEG_Run2017BCDE, DoubleEG_Run2017BCDE, DoubleMuon_Run2017BCDE, SingleElectron_Run2017BCDE, SingleMuon_Run2017BCDE], texName = "Run2017BCDE")
Run2017BCDE.lumi = (4.823+9.664+4.252+9.278)*1000
Run2017F = Sample.combine("Run2017F", [MuonEG_Run2017F, DoubleEG_Run2017F, DoubleMuon_Run2017F, SingleElectron_Run2017F, SingleMuon_Run2017F], texName = "Run2017F")
Run2017F.lumi = (13.540)*1000

allSamples_Data25ns += [Run2017, Run2017B, Run2017CDE, Run2017BCDE, Run2017F]

for s in allSamples_Data25ns:
  s.color   = ROOT.kBlack
  s.isData  = True


