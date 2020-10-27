import copy, os, sys
from RootTools.core.Sample import Sample
import ROOT

# Logging
import logging
logger = logging.getLogger(__name__)

from tWZ.samples.color import color

# Data directory
import tWZ.samples.nanoAODv6_private_locations as locations

logger.info("Loading MC samples from directory %s and %s", locations.fakes_mu_2016, locations.fakes_ele_2016)

def make_dirs( base, dirs ):
    return [ os.path.join( base, dir_ ) for dir_ in dirs ]

dirs = {}

# mu enriched QCD (few bins missing!) 
dirs['QCD_mu'] = ["QCD_Mu_pt1000toInf_comb", "QCD_Mu_pt120to170", "QCD_Mu_pt15to20", "QCD_Mu_pt170to300_comb", "QCD_Mu_pt20to30", "QCD_Mu_pt300to470_comb", "QCD_Mu_pt30to50", "QCD_Mu_pt470to600_comb", "QCD_Mu_pt50to80", "QCD_Mu_pt600to800_comb", "QCD_Mu_pt800to1000_comb", "QCD_Mu_pt80to120_comb"]
QCD_mu = Sample.fromDirectory(name="QCD_mu", treeName="Events", isData=False, color=color.QCD, texName="QCD(#mu)", directory=make_dirs( locations.fakes_mu_2016, dirs['QCD_mu']))

# pt binned inclusive QCD (mu)
#dirs['QCD_pt'] = ["QCD_pt15to30", "QCD_pt30to50", "QCD_pt50to80", "QCD_pt80to120_comb", "QCD_pt120to170_comb", "QCD_pt170to300_comb", "QCD_pt300to470_comb", "QCD_pt470to600", "QCD_pt600to800_comb", "QCD_pt1400to1800_comb", "QCD_pt1800to2400_comb", "QCD_pt2400to3200_comb",]
#QCD_pt_mu = Sample.fromDirectory(name="QCD_pt_mu", treeName="Events", isData=False, color=color.QCD, texName="QCD p_{T} (#mu)", directory=make_dirs( locations.fakes_mu_2016, dirs['QCD_pt']))

# EM enriched QCD (bcToE missing!)
dirs['QCD_ele'] = [ "QCD_Ele_pt20to30", "QCD_Ele_pt30to50_comb", "QCD_Ele_pt50to80_comb", "QCD_Ele_pt80to120_comb", "QCD_Ele_pt120to170_comb", "QCD_Ele_pt170to300", "QCD_Ele_pt300toInf" ] # "QCD_Ele_pt20to30"

QCD_ele = Sample.fromDirectory(name="QCD_ele", treeName="Events", isData=False, color=color.QCD, texName="QCD(e)", directory=make_dirs( locations.fakes_ele_2016, dirs['QCD_ele']))

# pt binned inclusive QCD (ele)
#QCD_pt_ele = Sample.fromDirectory(name="QCD_pt_ele", treeName="Events", isData=False, color=color.QCD, texName="QCD p_{T} (e)", directory=make_dirs( locations.fakes_ele_2016, dirs['QCD_pt']))

# TT
dirs['TTbar'] = ['TTbar']
TTbar_mu   = Sample.fromDirectory(name="TTbar_mu", treeName="Events", isData=False, color=color.TTJets, texName="t/t#bar{t}", directory=make_dirs(locations.fakes_mu_2016, dirs['TTbar']))
TTbar_ele  = Sample.fromDirectory(name="TTbar_ele", treeName="Events", isData=False, color=color.TTJets, texName="t/t#bar{t}", directory=make_dirs(locations.fakes_ele_2016, dirs['TTbar']))

# WJets
dirs['WJetsToLNu'] = ['WJetsToLNu_comb']
WJetsToLNu_mu  = Sample.fromDirectory(name="WJetsToLNu_mu", treeName="Events", isData=False, color=color.WJets, texName="W+jets", directory=make_dirs(locations.fakes_mu_2016, dirs['WJetsToLNu']))
WJetsToLNu_ele = Sample.fromDirectory(name="WJetsToLNu_ele", treeName="Events", isData=False, color=color.WJets, texName="W+jets", directory=make_dirs(locations.fakes_ele_2016, dirs['WJetsToLNu']))

# DY
dirs['DY_LO'] = ['DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO_ext1_comb'] #,'DYJetsToLL_M10to50_LO']
DY_LO_mu  = Sample.fromDirectory(name="DY_mu", treeName="Events", isData=False, color=color.DY, texName="DY(LO)+Jets", directory=make_dirs(locations.fakes_mu_2016, dirs['DY_LO']))
DY_LO_ele = Sample.fromDirectory(name="DY_ele", treeName="Events", isData=False, color=color.DY, texName="DY(LO)+Jets", directory=make_dirs(locations.fakes_ele_2016, dirs['DY_LO']))

for (run, version) in [('B','_ver2'),('C',''),('D',''),('E',''),('F',''),('G',''),('H','')]: # no event that passes json in B_ver1
    runTag = 'Run2016' + run + '_25Oct2019' + version
    dirs["DoubleEG_Run2016"         + run + version ] = ["DoubleEG_"          + runTag ]
    dirs["DoubleMuon_Run2016"       + run + version ] = ["DoubleMuon_"        + runTag ]
    #dirs["SingleElectron_Run2016"   + run + version ] = ["SingleElectron_"    + runTag ]
    #dirs["SingleMuon_Run2016"       + run + version ] = ["SingleMuon_"        + runTag ]
    #dirs["MuonEG_Run2016"           + run + version ] = ["MuonEG_"            + runTag ]

def merge(pd, totalRunName, listOfRuns):
    dirs[pd + '_' + totalRunName] = []
    for run in listOfRuns: dirs[pd + '_' + totalRunName].extend(dirs[pd + '_' + run])

for pd in [ 'DoubleMuon', 'DoubleEG']:#, 'SingleElectron', 'SingleMuon']:
    merge(pd, 'Run2016BCD',    ['Run2016B_ver2', 'Run2016C', 'Run2016D'])
    merge(pd, 'Run2016BCDEFG', ['Run2016BCD', 'Run2016E', 'Run2016F', 'Run2016G'])
    merge(pd, 'Run2016EF', ['Run2016E', 'Run2016F'])
    merge(pd, 'Run2016GH', ['Run2016G', 'Run2016H'])
    merge(pd, 'Run2016',       ['Run2016BCDEFG', 'Run2016H'])

def getSample(base, pd, runName, lumi):
    sample      = Sample.fromDirectory(name=(pd + '_' + runName), treeName="Events", texName=(pd + ' (' + runName + ')'), directory=map( lambda d: os.path.join( base, d), dirs[pd + '_' + runName]) )
    sample.lumi = lumi
    return sample

l_tot = (5.883+2.646+4.353+4.050+3.124+7.554+8.494+0.217) # old normtag

allSamples_Data25ns = []

DoubleMuon_Run2016   = getSample(locations.fakes_mu_2016,  'DoubleMuon',       'Run2016',           (35.9)*1000)
DoubleEG_Run2016     = getSample(locations.fakes_ele_2016, 'DoubleEG',       'Run2016',           (35.9)*1000)
allSamples_Data25ns += [ DoubleMuon_Run2016, DoubleEG_Run2016]

for s in allSamples_Data25ns:
  s.color   = ROOT.kBlack
  s.isData  = True
