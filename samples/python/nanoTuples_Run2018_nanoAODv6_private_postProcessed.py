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
    directory_ = locations.data_2018

logger.info("Loading data samples from directory %s", directory_)

dirs = {}
for (run, version) in [('A',''), ('B',''), ('C',''), ('D','')]:
    runTag = 'Run2018' + run + '_25Oct2019' + version
    dirs["EGamma_Run2018"           + run + version ] = ["EGamma_"            + runTag ]
    dirs["DoubleMuon_Run2018"       + run + version ] = ["DoubleMuon_"        + runTag ]
    dirs["MuonEG_Run2018"           + run + version ] = ["MuonEG_"            + runTag ]
    dirs["SingleMuon_Run2018"       + run + version ] = ["SingleMuon_"        + runTag ]

def merge(pd, totalRunName, listOfRuns):
    dirs[pd + '_' + totalRunName] = []
    for run in listOfRuns: dirs[pd + '_' + totalRunName].extend(dirs[pd + '_' + run])

for pd in ['MuonEG', 'DoubleMuon', 'EGamma', 'SingleMuon']:
    merge(pd, 'Run2018',    ['Run2018A', 'Run2018B', 'Run2018C', 'Run2018D'])

for key in dirs:
    dirs[key] = [ os.path.join( directory_, dir) for dir in dirs[key]]


def getSample(pd, runName, lumi):
    sample      = Sample.fromDirectory(name=(pd + '_' + runName), treeName="Events", texName=(pd + ' (' + runName + ')'), directory=dirs[pd + '_' + runName])
    sample.lumi = lumi
    return sample

allSamples_Data25ns = []
EGamma_Run2018A                  = getSample('EGamma',           'Run2018A',       ((14.00)*1000))
DoubleMuon_Run2018A              = getSample('DoubleMuon',       'Run2018A',       ((14.00)*1000))
MuonEG_Run2018A                  = getSample('MuonEG',           'Run2018A',       ((14.00)*1000))
SingleMuon_Run2018A              = getSample('SingleMuon',       'Run2018A',       ((14.00)*1000))
allSamples_Data25ns += [MuonEG_Run2018A, EGamma_Run2018A, DoubleMuon_Run2018A, SingleMuon_Run2018A]

EGamma_Run2018B                  = getSample('EGamma',           'Run2018B',       ((7.10)*1000))
DoubleMuon_Run2018B              = getSample('DoubleMuon',       'Run2018B',       ((7.10)*1000))
MuonEG_Run2018B                  = getSample('MuonEG',           'Run2018B',       ((7.10)*1000))
SingleMuon_Run2018B              = getSample('SingleMuon',       'Run2018B',       ((7.10)*1000))
allSamples_Data25ns += [MuonEG_Run2018B, EGamma_Run2018B, DoubleMuon_Run2018B, SingleMuon_Run2018B]

EGamma_Run2018C                  = getSample('EGamma',           'Run2018C',       ((6.94)*1000))
DoubleMuon_Run2018C              = getSample('DoubleMuon',       'Run2018C',       ((6.94)*1000))
MuonEG_Run2018C                  = getSample('MuonEG',           'Run2018C',       ((6.94)*1000))
SingleMuon_Run2018C              = getSample('SingleMuon',       'Run2018C',       ((6.94)*1000))
allSamples_Data25ns += [MuonEG_Run2018C, EGamma_Run2018C, DoubleMuon_Run2018C, SingleMuon_Run2018C]

EGamma_Run2018D                  = getSample('EGamma',           'Run2018D',       ((31.93)*1000))
DoubleMuon_Run2018D              = getSample('DoubleMuon',       'Run2018D',       ((31.93)*1000))
MuonEG_Run2018D                  = getSample('MuonEG',           'Run2018D',       ((31.93)*1000))
SingleMuon_Run2018D              = getSample('SingleMuon',       'Run2018D',       ((31.93)*1000))
allSamples_Data25ns += [MuonEG_Run2018D, EGamma_Run2018D, DoubleMuon_Run2018D, SingleMuon_Run2018D]

EGamma_Run2018                   = getSample('EGamma',           'Run2018',        ((59.97)*1000))
DoubleMuon_Run2018               = getSample('DoubleMuon',       'Run2018',        ((59.97)*1000))
MuonEG_Run2018                   = getSample('MuonEG',           'Run2018',        ((59.97)*1000))
SingleMuon_Run2018               = getSample('SingleMuon',       'Run2018',        ((59.97)*1000))
allSamples_Data25ns += [MuonEG_Run2018, EGamma_Run2018, DoubleMuon_Run2018, SingleMuon_Run2018]

Run2018A = Sample.combine("Run2018A", [MuonEG_Run2018A, EGamma_Run2018A, DoubleMuon_Run2018A, SingleMuon_Run2018A], texName = "Run2018")
Run2018A.lumi = (14.00)*1000
allSamples_Data25ns.append( Run2018A )

Run2018B = Sample.combine("Run2018B", [MuonEG_Run2018B, EGamma_Run2018B, DoubleMuon_Run2018B, SingleMuon_Run2018B], texName = "Run2018B")
Run2018B.lumi = (7.10)*1000
allSamples_Data25ns.append( Run2018B )

Run2018C = Sample.combine("Run2018C", [MuonEG_Run2018C, EGamma_Run2018C, DoubleMuon_Run2018C, SingleMuon_Run2018C], texName = "Run2018C")
Run2018C.lumi = (6.94)*1000
allSamples_Data25ns.append( Run2018C )

Run2018ABC = Sample.combine("Run2018ABC", [Run2018A, Run2018B, Run2018C], texName = "Run2018ABC")
Run2018ABC.lumi = (14.00+7.10+6.94)*1000
allSamples_Data25ns.append( Run2018ABC )

Run2018D = Sample.combine("Run2018D", [MuonEG_Run2018D, EGamma_Run2018D, DoubleMuon_Run2018D, SingleMuon_Run2018D], texName = "Run2018D")
Run2018D.lumi = (31.93)*1000
allSamples_Data25ns.append( Run2018D )

Run2018 = Sample.combine("Run2018", [MuonEG_Run2018, EGamma_Run2018, DoubleMuon_Run2018, SingleMuon_Run2018], texName = "Run2018")
Run2018.lumi = (59.97)*1000
allSamples_Data25ns.append( Run2018 )

for s in allSamples_Data25ns:
  s.color   = ROOT.kBlack
  s.isData  = True


