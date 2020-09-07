import copy, os, sys
from RootTools.core.Sample import Sample
import ROOT

# Logging
import logging
logger = logging.getLogger(__name__)

from tWZ.samples.color import color

# Data directory
try:
    directory_ = sys.modules['__main__'].directory_
except:
    import tWZ.samples.nanoAODv6_private_locations as locations
    directory_ = locations.fakes_2016

logger.info("Loading MC samples from directory %s", directory_)

def make_dirs( dirs ):
    return [ os.path.join( directory_, dir_ ) for dir_ in dirs ]

dirs = {}

# diboson exclusive
dirs['QCD_Mu']               = ["QCD_Mu_pt1000toInf_comb", "QCD_Mu_pt120to170", "QCD_Mu_pt15to20", "QCD_Mu_pt170to300_comb", "QCD_Mu_pt20to30", "QCD_Mu_pt300to470_comb", "QCD_Mu_pt30to50", "QCD_Mu_pt470to600_comb", "QCD_Mu_pt50to80", "QCD_Mu_pt600to800_comb", "QCD_Mu_pt800to1000_comb", "QCD_Mu_pt80to120_comb"]
QCD_Mu = Sample.fromDirectory(name="QCD_Mu", treeName="Events", isData=False, color=color.QCD, texName="QCD(#mu)", directory=make_dirs( dirs['QCD_Mu']))

# TT
dirs['TTbar']           = ['TTbar']
TTbar  = Sample.fromDirectory(name="TTbar", treeName="Events", isData=False, color=color.TTJets, texName="t/t#bar{t}", directory=make_dirs(dirs['TTbar']))

# WJets
dirs['WJetsToLNu']           = ['WJetsToLNu_comb']
WJetsToLNu  = Sample.fromDirectory(name="WJetsToLNu", treeName="Events", isData=False, color=color.WJets, texName="W+jets", directory=make_dirs(dirs['WJetsToLNu']))

# DY
dirs['DY_LO']              = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO_ext1_comb'] #,'DYJetsToLL_M10to50_LO']
DY_LO  = Sample.fromDirectory(name="DY", treeName="Events", isData=False, color=color.DY, texName="DY(LO)+Jets", directory=make_dirs(dirs['DY_LO']))

for (run, version) in [('B','_ver2'),('C',''),('D',''),('E',''),('F',''),('G',''),('H','')]: # no event that passes json in B_ver1
    runTag = 'Run2016' + run + '_25Oct2019' + version
    #dirs["DoubleEG_Run2016"         + run + version ] = ["DoubleEG_"          + runTag ]
    dirs["DoubleMuon_Run2016"       + run + version ] = ["DoubleMuon_"        + runTag ]
    #dirs["SingleElectron_Run2016"   + run + version ] = ["SingleElectron_"    + runTag ]
    #dirs["SingleMuon_Run2016"       + run + version ] = ["SingleMuon_"        + runTag ]
    #dirs["MuonEG_Run2016"           + run + version ] = ["MuonEG_"            + runTag ]

def merge(pd, totalRunName, listOfRuns):
    dirs[pd + '_' + totalRunName] = []
    for run in listOfRuns: dirs[pd + '_' + totalRunName].extend(dirs[pd + '_' + run])

for pd in [ 'DoubleMuon' ]:#, 'DoubleEG', 'SingleElectron', 'SingleMuon']:
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

DoubleMuon_Run2016   = getSample('DoubleMuon',       'Run2016',           (35.9)*1000)
allSamples_Data25ns += [ DoubleMuon_Run2016 ]

for s in allSamples_Data25ns:
  s.color   = ROOT.kBlack
  s.isData  = True
