#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT, os
ROOT.gROOT.SetBatch(True)
import itertools
import copy
import array
import operator
from math                                import sqrt, cos, sin, pi, atan2, cosh

# RootTools
from RootTools.core.standard             import *

# tWZ
from tWZ.Tools.user                      import plot_directory
from tWZ.Tools.cutInterpreter            import cutInterpreter
from tWZ.Tools.objectSelection           import cbEleIdFlagGetter, vidNestedWPBitMapNamingList
from tWZ.Tools.objectSelection           import mu_string, ele_string, getJets, getGoodElectrons, getGoodMuons, eleSelector, muonSelector 
# Analysis
from Analysis.Tools.helpers              import deltaPhi, deltaR
from Analysis.Tools.puProfileCache       import *
from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons
from Analysis.Tools.puReweighting        import getReweightingFunction
from Analysis.Tools.metFilters               import getFilterCut

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--noData',         action='store_true', default=False, help='also plot data?')
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--dataMCScaling',  action='store_true', help='Data MC scaling?', )
argParser.add_argument('--plot_directory', action='store', default='tWZ_v3')
argParser.add_argument('--era',            action='store', type=str, default="2016")
argParser.add_argument('--selection',      action='store', default='trilepMini0p12-minDLmass12-onZ1-njet4p-btag2p')
argParser.add_argument('--year',        action='store',                     type=int,                                               help="Which year?" )
argParser.add_argument('--nanoAODv4',   action='store_true',                                                                        help="Run on nanoAODv4?" )
argParser.add_argument('--samples',     action='store',         nargs='*',  type=str, default=['TTZToLLNuNu_ext'],                  help="List of samples to be post-processed, given as CMG component name" )
argParser.add_argument('--event',                       action='store',     type=int, default=-1,                                   help="Just process event no")
argParser.add_argument('--skim',        action='store',         nargs='?',  type=str, default='trilep',                             help="Skim conditions to be applied for post-processing" )
argParser.add_argument('--triggerSelection',            action='store_true',                                                        help="Trigger selection?" ) 
argParser.add_argument('--LHEHTCut',    action='store',         nargs='?',  type=int, default=-1,                                   help="LHE cut." )



args = argParser.parse_args()

options = argParser.parse_args()

# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"
if args.noData:                       args.plot_directory += "_noData"

logger.info( "Working in era %s", args.era)

#tWZ_sample = TWZ if args.nominalSignal else yt_TWZ_filter

from tWZ.samples.nanoTuples_RunII_nanoAODv4_postProcessed import *

if args.era == "Run2016":
    mc = [Summer16.TWZ, Summer16.TTZ, Summer16.TTX_rare, Summer16.TZQ, Summer16.WZ, Summer16.triBoson, Summer16.ZZ, Summer16.nonprompt_3l]
elif args.era == "Run2017":
    mc = [Fall17.TWZ, Fall17.TTZ, Fall17.TTX_rare, Fall17.TZQ, Fall17.WZ, Fall17.triBoson, Fall17.ZZ, Fall17.nonprompt_3l]
elif args.era == "Run2018":
    mc = [Autumn18.TWZ, Autumn18.TTZ, Autumn18.TTX_rare, Autumn18.TZQ, Autumn18.WZ, Autumn18.triBoson, Autumn18.ZZ, Autumn18.nonprompt_3l]
elif args.era == "RunII":
    mc = [TWZ, TTZ, TTX_rare, TZQ, WZ, triBoson, ZZ, nonprompt_3l]

# data sample
try:
  data_sample = eval(args.era)
except Exception as e:
  logger.error( "Didn't find %s", args.era )
  raise e

lumi_scale                 = data_sample.lumi/1000
data_sample.scale          = 1.
for sample in mc:
    sample.scale           = 1 # Scale MCs individually with lumi

if args.small:
    for sample in mc + [data_sample]:
        sample.normalization = 1.
        #sample.reduceFiles( factor = 40 )
        sample.reduceFiles( to=1)
        sample.scale /= sample.normalization

# Text on the plots
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)
tex.SetTextAlign(11) # align right

def charge(pdgId):
    return -pdgId/abs(pdgId)

def drawObjects( plotData, dataMCScale, lumi_scale ):
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if plotData else 'CMS Simulation'), 
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    if "mt2ll100" in args.selection and args.noData: lines += [(0.55, 0.5, 'M_{T2}(ll) > 100 GeV')] # Manually put the mt2ll > 100 GeV label
    return [tex.DrawLatex(*l) for l in lines] 

def drawPlots(plots, mode, dataMCScale):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, mode + ("_log" if log else ""), args.selection)
    for plot in plots:
      if not max(l.GetMaximum() for l in sum(plot.histos,[])): continue # Empty plot
      if not args.noData: 
        if mode == "all": plot.histos[1][0].legendText = "Data"
        if mode == "SF":  plot.histos[1][0].legendText = "Data (SF)"

      _drawObjects = []

      if isinstance( plot, Plot):
          plotting.draw(plot,
            plot_directory = plot_directory_,
            ratio = {'yRange':(0.1,1.9)} if not args.noData else None,
            logX = False, logY = log, sorting = True,
            yRange = (0.03, "auto") if log else (0.001, "auto"),
            scaling = {0:1} if args.dataMCScaling else {},
            legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
            drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ) + _drawObjects,
            copyIndexPHP = True, extensions = ["png"],
          )

#Samples: Load samples
maxN = 1 if options.small else None
if options.small:
    options.job = 0
    options.nJobs = 10000 # set high to just run over 1 input file

if options.nanoAODv4:
    if options.year == 2016:
        from Samples.nanoAOD.Summer16_private_legacy_v1 import allSamples as mcSamples
        from Samples.nanoAOD.Run2016_17Jul2018_private  import allSamples as dataSamples
        allSamples = mcSamples + dataSamples
    elif options.year == 2017:
        from Samples.nanoAOD.Fall17_private_legacy_v1   import allSamples as mcSamples
        from Samples.nanoAOD.Run2017_31Mar2018_private  import allSamples as dataSamples
        allSamples = mcSamples + dataSamples
    elif options.year == 2018:
        from Samples.nanoAOD.Autumn18_private_legacy_v1 import allSamples as mcSamples
        from Samples.nanoAOD.Run2018_17Sep2018_private  import allSamples as dataSamples
        allSamples = mcSamples + dataSamples
    else:
        raise NotImplementedError
else:
    if options.year == 2016:
        from Samples.nanoAOD.Summer16_nanoAODv6         import allSamples as mcSamples
        from Samples.nanoAOD.Run2016_nanoAODv6          import allSamples as dataSamples
        allSamples = mcSamples + dataSamples
    elif options.year == 2017:
        from Samples.nanoAOD.Fall17_nanoAODv6           import allSamples as mcSamples
        from Samples.nanoAOD.Run2017_nanoAODv6          import allSamples as dataSamples
        allSamples = mcSamples + dataSamples
    elif options.year == 2018:
        from Samples.nanoAOD.Autumn18_nanoAODv6         import allSamples as mcSamples
        from Samples.nanoAOD.Run2018_nanoAODv6          import allSamples as dataSamples
        allSamples = mcSamples + dataSamples
    else:
        raise NotImplementedError

samples = []
for selectedSamples in options.samples:
    for sample in allSamples:
        if selectedSamples == sample.name:
            samples.append(sample)

#if len(samples)==0:
 #   logger.info( "No samples found. Was looking for %s. Exiting" % options.samples )
  #  sys.exit(-1)

isData = False not in [s.isData for s in samples]
isMC   =  True not in [s.isData for s in samples]

# Flags 
isDiLep         = options.skim.lower().startswith('dilep')
isTriLep        = options.skim.lower().startswith('trilep')
isSingleLep     = options.skim.lower().startswith('singlelep')
isSmall         = options.skim.lower().count('small')
isInclusive     = options.skim.lower().count('inclusive') 

# Skim condition
skimConds = []

if options.event > 0:
    skimConds.append( "event == %s"%options.event )

if isDiLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.4) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.4)>=2" )
if isTriLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)&&Electron_pfRelIso03_all<0.4) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5&&Muon_pfRelIso03_all<0.4)>=2 && Sum$(Electron_pt>10&&abs(Electron_eta)<2.5)+Sum$(Muon_pt>10&&abs(Muon_eta)<2.5)>=3" )
elif isSingleLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.5) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5)>=1" )

if isInclusive:
    skimConds.append('(1)')
    isSingleLep = True #otherwise no lepton variables?!
    isDiLep     = True
    isTriLep    = True

# Trigger selection
from tWZ.Tools.triggerSelector import triggerSelector
ts           = triggerSelector(options.year)
triggerCond  = ts.getSelection(options.samples[0] if isData else "MC")
treeFormulas = {"triggerDecision": {'string':triggerCond} }

            
if isData and options.triggerSelection:
    logger.info("Sample will have the following trigger skim: %s"%triggerCond)
    skimConds.append( triggerCond )

# apply MET filter
skimConds.append( getFilterCut(options.year, isData=isData, ignoreJSON=True, skipWeight=True) )

# LHE cut (DY samples)
if options.LHEHTCut>0:
    sample.name+="_lheHT"+str(options.LHEHTCut)
    logger.info( "Adding upper LHE cut at %f", options.LHEHTCut )
    skimConds.append( "LHE_HTIncoming<%f"%options.LHEHTCut )

#
# Read variables and sequences
#
jetMCInfo   = ['genJetIdx/I','hadronFlavour/I']

jetVars         = ['pt/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'rawFactor/F', 'eta/F', 'phi/F', 'jetId/I', 'btagDeepB/F', 'btagDeepFlavB/F', 'btagCSVV2/F', 'area/F', 'pt_nom/F', 'corr_JER/F'] 
if isMC:
    jetVars     += jetMCInfo
    jetVars     += ['pt_jesTotalUp/F', 'pt_jesTotalDown/F', 'pt_jerUp/F', 'pt_jerDown/F', 'corr_JER/F', 'corr_JEC/F']
jetVarNames     = [x.split('/')[0] for x in jetVars]
print jetVarNames

#ppvariables

read_variables =  [ 'MET_pt/F', 'MET_phi/F', 'run/I', 'luminosityBlock/I', 'event/l', 'PV_npvs/I', 'PV_npvsGood/I'] 
if options.year == 2017:
    read_variables += [ 'METFixEE2017_pt/F', 'METFixEE2017_phi/F', 'METFixEE2017_pt_nom/F', 'METFixEE2017_phi_nom/F']
    if isMC:
        read_variables += [ 'METFixEE2017_pt_jesTotalUp/F', 'METFixEE2017_pt_jesTotalDown/F', 'METFixEE2017_pt_jerUp/F', 'METFixEE2017_pt_jerDown/F', 'METFixEE2017_pt_unclustEnDown/F', 'METFixEE2017_pt_unclustEnUp/F', 'METFixEE2017_phi_jesTotalUp/F', 'METFixEE2017_phi_jesTotalDown/F', 'METFixEE2017_phi_jerUp/F', 'METFixEE2017_phi_jerDown/F', 'METFixEE2017_phi_unclustEnDown/F', 'METFixEE2017_phi_unclustEnUp/F', 'METFixEE2017_pt_jer/F', 'METFixEE2017_phi_jer/F']
else:
    read_variables += [ 'MET_pt_nom/F', 'MET_phi_nom/F' ]
    if isMC:
        read_variables += [ 'MET_pt_jesTotalUp/F', 'MET_pt_jesTotalDown/F', 'MET_pt_jerUp/F', 'MET_pt_jerDown/F', 'MET_pt_unclustEnDown/F', 'MET_pt_unclustEnUp/F', 'MET_phi_jesTotalUp/F', 'MET_phi_jesTotalDown/F', 'MET_phi_jerUp/F', 'MET_phi_jerDown/F', 'MET_phi_unclustEnDown/F', 'MET_phi_unclustEnUp/F', 'MET_pt_jer/F', 'MET_phi_jer/F']
if isMC:
    read_variables += [ 'GenMET_pt/F', 'GenMET_phi/F' ]

read_variables += [ TreeVariable.fromString('nPhoton/I'),
                    VectorTreeVariable.fromString('Photon[pt/F,eta/F,phi/F,mass/F,cutBased/I,pdgId/I]') if (options.year == 2016) else VectorTreeVariable.fromString('Photon[pt/F,eta/F,phi/F,mass/F,cutBasedBitmap/I,pdgId/I]') ]

#new_variables = [ 'weight/F', 'triggerDecision/I', 'year/I']
if isMC:
    read_variables += [ TreeVariable.fromString('Pileup_nTrueInt/F') ]
    # reading gen particles for top pt reweighting
    read_variables.append( 'nGenPart/I' )
    read_variables.append( 'nISR/I' )
    read_variables.append( 'GenPart[pt/F,mass/F,phi/F,eta/F,pdgId/I,genPartIdxMother/I,status/I,statusFlags/I]') # default nMax is 100, which would lead to corrupt values in this case
    read_variables.append( 'genWeight/F' )
    read_variables.append( 'nGenJet/I' )
    read_variables.append( 'GenJet[pt/F,eta/F,phi/F]'  )
#    new_variables.extend([ 'reweightTopPt/F', 'reweight_nISR/F', 'reweight_nISRUp/F', 'reweight_nISRDown/F', 'reweightPU/F','reweightPUUp/F','reweightPUDown/F', 'reweightPUVUp/F','reweightPUVVUp/F', 'reweightPUVDown/F', 'reweightL1Prefire/F', 'reweightL1PrefireUp/F', 'reweightL1PrefireDown/F'])
#    if not options.skipGenLepMatching:
#        new_variables.append( TreeVariable.fromString( 'nGenLep/I' ) )
#        new_variables.append( 'GenLep[%s]'% ( ','.join(genLepVars) ) )
#    if options.doCRReweighting:
#        new_variables.append('reweightCR/F')

read_variables += [\
    'nElectron/I',
    'Electron[pt/F,eta/F,phi/F,pdgId/I,cutBased/I,miniPFRelIso_all/F,pfRelIso03_all/F,sip3d/F,lostHits/b,mvaFall17V2Iso_WP80/O,mvaFall17V2Iso_WP90/O,convVeto/O,dxy/F,dz/F,charge/I,deltaEtaSC/F,vidNestedWPBitmap/I,mvaTTH/F]',
    'nMuon/I',
    'Muon[pt/F,eta/F,phi/F,pdgId/I,mediumId/O,miniPFRelIso_all/F,pfRelIso03_all/F,sip3d/F,dxy/F,dz/F,charge/I,mvaTTH/F]',
    'nJet/I',
    'Jet[%s]'% ( ','.join(jetVars) )  ]

#new_variables += [\
#    'overlapRemoval/I','nlep/I',
#    'JetGood[%s]'% ( ','.join(jetVars+['index/I']) + ',genPt/F' ),
#    'met_pt/F', 'met_phi/F', 'met_pt_min/F'
#]
#
#if sample.isData: new_variables.extend( ['jsonPassed/I','isData/I'] )
#new_variables.extend( ['nBTag/I', 'm3/F', 'minDLmass/F'] )
#
#new_variables.append( 'lep[%s]'% ( ','.join(lepVars) ) )
#
#if isTriLep or isDiLep or isSingleLep:
#    new_variables.extend( ['nGoodMuons/I', 'nGoodElectrons/I', 'nGoodLeptons/I' ] )
#    new_variables.extend( ['l1_pt/F', 'l1_eta/F', 'l1_phi/F', 'l1_pdgId/I', 'l1_index/I', 'l1_jetPtRelv2/F', 'l1_jetPtRatiov2/F', 'l1_miniRelIso/F', 'l1_relIso03/F', 'l1_dxy/F', 'l1_dz/F', 'l1_mIsoWP/I', 'l1_eleIndex/I', 'l1_muIndex/I' ] )
#    new_variables.extend( ['mlmZ_mass/F'])
#    if isMC: 
#        new_variables.extend(['reweightLeptonSF/F', 'reweightLeptonSFUp/F', 'reweightLeptonSFDown/F'])
#if isTriLep or isDiLep:
#    new_variables.extend( ['l2_pt/F', 'l2_eta/F', 'l2_phi/F', 'l2_pdgId/I', 'l2_index/I', 'l2_jetPtRelv2/F', 'l2_jetPtRatiov2/F', 'l2_miniRelIso/F', 'l2_relIso03/F', 'l2_dxy/F', 'l2_dz/F', 'l2_mIsoWP/I', 'l2_eleIndex/I', 'l2_muIndex/I' ] )
#    if isMC: new_variables.extend( \
#        [   'genZ1_pt/F', 'genZ1_eta/F', 'genZ1_phi/F',
#            'genZ2_pt/F', 'genZ1_eta/F', 'genZ1_phi/F',  
#            'reweightTrigger/F', 'reweightTriggerUp/F', 'reweightTriggerDown/F',
#            'reweightLeptonTrackingSF/F',
#         ] )
#if isTriLep:
#    new_variables.extend( ['l3_pt/F', 'l3_eta/F', 'l3_phi/F', 'l3_pdgId/I', 'l3_index/I', 'l3_jetPtRelv2/F', 'l3_jetPtRatiov2/F', 'l3_miniRelIso/F', 'l3_relIso03/F', 'l3_dxy/F', 'l3_dz/F', 'l3_mIsoWP/I', 'l3_eleIndex/I', 'l3_muIndex/I' ] )
#new_variables.extend( ['nPhotonGood/I','photon_pt/F','photon_eta/F','photon_phi/F','photon_idCutBased/I'] )
#if isMC: new_variables.extend( ['photon_genPt/F', 'photon_genEta/F', 'genZ_mass/F'] )
#new_variables.extend( ['photonJetdR/F','photonLepdR/F'] )
#
### ttZ related variables
#new_variables.extend( ['Z1_l1_index/I', 'Z1_l2_index/I', 'Z2_l1_index/I', 'Z2_l2_index/I', 'nonZ1_l1_index/I', 'nonZ1_l2_index/I'] )
#for i in [1,2]:
#    new_variables.extend( ['Z%i_pt/F'%i, 'Z%i_eta/F'%i, 'Z%i_phi/F'%i, 'Z%i_lldPhi/F'%i, 'Z%i_lldR/F'%i,  'Z%i_mass/F'%i, 'Z%i_cosThetaStar/F'%i] )
#
#if options.checkTTGJetsOverlap:
#    new_variables.extend( ['TTGJetsEventType/I'] )
#
#if addSystematicVariations:
#    for var in ['jesTotalUp', 'jesTotalDown', 'jerUp', 'jer', 'jerDown', 'unclustEnUp', 'unclustEnDown']:
#        if not var.startswith('unclust'):
#            new_variables.extend( ['nJetGood_'+var+'/I', 'nBTag_'+var+'/I'] )
#        new_variables.extend( ['met_pt_'+var+'/F', 'met_phi_'+var+'/F'] )
#
## Btag weights Method 1a
#for var in btagEff.btagWeightNames:
#    if var!='MC':
#        new_variables.append('reweightBTag_'+var+'/F')


reader = sample.treeReader( \
    variables = read_variables,
    selectionString = "&&".join(skimConds) 
    )   

eleSelector_ = eleSelector( "tightMiniIso02", year = options.year )
muSelector_  = muonSelector("tightMiniIso02", year = options.year )

#def getjets(event, sample):
r = reader.event

# get leptons before jets in order to clean jets
electrons_pt10  = getGoodElectrons(r, ele_selector = eleSelector_)
muons_pt10      = getGoodMuons    (r, mu_selector  = muSelector_ )

for e in electrons_pt10:
    e['pdgId']      = int( -11*e['charge'] )
    e['eleIndex']   = e['index']
    e['muIndex']    = -1
for m in muons_pt10:
    m['pdgId']      = int( -13*m['charge'] )
    m['muIndex']    = m['index']
    m['eleIndex']   = -1

# make list of leptons
leptons = electrons_pt10+muons_pt10
leptons.sort(key = lambda p:-p['pt'])

# now get jets, cleaned against good leptons
all_jets     = getJets(r, jetVars=jetVarNames)
clean_jets,_ = cleanJetsAndLeptons( all_jets, leptons )

jets         = filter(lambda j:j['pt']>30, clean_jets)

# Filling jets
maxNJet = 100 
store_jets = jets #if not options.keepAllJets else soft_jets + jets
store_jets = store_jets[:maxNJet]
store_jets.sort( key = lambda j:-j['pt'])
event.nJetGood   = len(store_jets)

for iJet, jet in enumerate(store_jets):
    event.JetGoodwoetacut_index[iJet] = jet['index']
    for b in jetVarNames:
        getattr(event, "JetGoodwoetacut_"+b)[iJet] = jet[b]
    event.abseta = abs( event.Jet['eta'])

sequence       = []

def getWpt( event, sample):

    # get the lepton and met
    lepton  = ROOT.TLorentzVector()
    met     = ROOT.TLorentzVector()
    lepton.SetPtEtaPhiM(event.lep_pt[event.nonZ1_l1_index], event.lep_eta[event.nonZ1_l1_index], event.lep_phi[event.nonZ1_l1_index], 0)
    met.SetPtEtaPhiM(event.met_pt, 0, event.met_phi, 0)

    # get the W boson candidate
    W   = lepton + met
    event.W_pt = W.Pt()

sequence.append( getWpt )

def getM3l( event, sample ):
    # get the invariant mass of the 3l system
    l = []
    for i in range(3):
        l.append(ROOT.TLorentzVector())
        l[i].SetPtEtaPhiM(event.lep_pt[i], event.lep_eta[i], event.lep_phi[i],0)
    event.M3l = (l[0] + l[1] + l[2]).M()

sequence.append( getM3l )

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l2_eta/F", "l2_phi/F", 
    "JetGood[pt/F,eta/F,phi/F]",
    "lep[pt/F,eta/F,phi/F,pdgId/I,muIndex/I,eleIndex/I]",
    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I", 
    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
]

read_variables_MC = ['reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F']

# define 3l selections

def getLeptonSelection( mode ):
    if   mode=="mumumu": return "Sum$({mu_string})==3&&Sum$({ele_string})==0".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="mumue":  return "Sum$({mu_string})==2&&Sum$({ele_string})==1".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="muee":   return "Sum$({mu_string})==1&&Sum$({ele_string})==2".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="eee":    return "Sum$({mu_string})==0&&Sum$({ele_string})==3".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=='all':    return "Sum$({mu_string})+Sum$({ele_string})==3".format(mu_string=mu_string,ele_string=ele_string)

# Getter functor for lepton quantities
def lep_getter( branch, index, abs_pdg = None, functor = None, debug=False):
    if functor is not None:
        if abs_pdg == 13:
            def func_( event, sample ):
                if debug:
                    print "Returning", "Muon_%s"%branch, index, abs_pdg, "functor", functor, "result",
                    print functor(getattr( event, "Muon_%s"%branch )[event.lep_muIndex[index]]) if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
                return functor(getattr( event, "Muon_%s"%branch )[event.lep_muIndex[index]]) if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
        else:
            def func_( event, sample ):
                return functor(getattr( event, "Electron_%s"%branch )[event.lep_eleIndex[index]]) if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
    else:
        if abs_pdg == 13:
            def func_( event, sample ):
                if debug:
                    print "Returning", "Muon_%s"%branch, index, abs_pdg, "functor", functor, "result",
                    print getattr( event, "Muon_%s"%branch )[event.lep_muIndex[index]] if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
                return getattr( event, "Muon_%s"%branch )[event.lep_muIndex[index]] if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
        else:
            def func_( event, sample ):
                return getattr( event, "Electron_%s"%branch )[event.lep_eleIndex[index]] if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
    return func_

#mu0_charge   = lep_getter("pdgId", 0, 13, functor = charge)
#ele0_charge = lep_getter("pdgId", 0, 11, functor = charge)
#mu1_charge   = lep_getter("pdgId", 1, 13, functor = charge)
#ele1_charge = lep_getter("pdgId", 1, 11, functor = charge)
#def test(event, sample):
#    mu0_ch  = mu0_charge(event, sample)
#    ele0_ch = ele0_charge(event, sample)
#    mu1_ch  = mu1_charge(event, sample)
#    ele1_ch = ele1_charge(event, sample)
#    print "mu0_ch",mu0_ch, "ele0_ch",ele0_ch, "mu1_ch",mu1_ch, "ele1_ch",ele1_ch
#
#sequence.append( test )

# 3l trainign variables

def make_training_observables_3l(event, sample):

    event.nonZ1l1_Z1_deltaPhi = deltaPhi(event.lep_phi[event.nonZ1_l1_index], event.Z1_phi)
    event.nonZ1l1_Z1_deltaEta = abs(event.lep_eta[event.nonZ1_l1_index] - event.Z1_eta)
    event.nonZ1l1_Z1_deltaR   = deltaR({'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet0_Z1_deltaR      = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet0_nonZ1l1_deltaR = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    event.jet1_Z1_deltaR      = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet1_nonZ1l1_deltaR = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})

sequence.append( make_training_observables_3l )

yields     = {}
allPlots   = {}
allModes   = ['mumumu','mumue','muee', 'eee']
for i_mode, mode in enumerate(allModes):
    yields[mode] = {}
    if not args.noData:
        data_sample.texName = "data"
        data_sample.setSelectionString([getLeptonSelection(mode)])
        data_sample.name           = "data"
        data_sample.style          = styles.errorStyle(ROOT.kBlack)
        lumi_scale                 = data_sample.lumi/1000

    weight_ = lambda event, sample: event.weight if sample.isData else event.weight*lumi_year[event.year]/1000.

    for sample in mc: sample.style = styles.fillStyle(sample.color)
    
    for sample in mc:
      sample.read_variables = read_variables_MC 
      sample.setSelectionString([getLeptonSelection(mode)])
      sample.weight = lambda event, sample: event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger#*event.reweightLeptonSF

    #yt_TWZ_filter.scale = lumi_scale * 1.07314

    if not args.noData:
      stack = Stack(mc, data_sample)
    else:
      stack = Stack(mc)

    # Use some defaults
    Plot.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = cutInterpreter.cutString(args.selection))

    plots = []


    plots.append(Plot(
      name = 'abseta', 
      texX = 'abs(#eta)', 
      texY = 'Number of Events',
      attribute = lambda event, sample: event.abseta_eta,
      binning=[20, 0, 5],
    ))


#    plots.append(Plot(
#      name = 'yield', texX = '', texY = 'Number of Events',
#      attribute = lambda event, sample: 0.5 + i_mode,
#      binning=[4, 0, 4],
#    ))
#
#    plots.append(Plot(
#      name = 'nVtxs', texX = 'vertex multiplicity', texY = 'Number of Events',
#      attribute = TreeVariable.fromString( "PV_npvsGood/I" ),
#      binning=[50,0,50],
#      addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = TreeVariable.fromString( "met_pt/F" ),
#        binning=[400/20,0,400],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        texX = '#phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
#        attribute = TreeVariable.fromString( "met_phi/F" ),
#        binning=[10,-pi,pi],
#    ))
#
#    plots.append(Plot(
#        name = "Z1_pt",
#        texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = TreeVariable.fromString( "Z1_pt/F" ),
#        binning=[20,0,400],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'Z1_pt_coarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events / 50 GeV',
#        attribute = TreeVariable.fromString( "Z1_pt/F" ),
#        binning=[16,0,800],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'Z1_pt_superCoarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
#        attribute = TreeVariable.fromString( "Z1_pt/F" ),
#        binning=[3,0,600],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'Z1_pt_coarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events / 50 GeV',
#        attribute = TreeVariable.fromString( "Z1_pt/F" ),
#        binning=[16,0,800],
#    ))
#
#    plots.append(Plot(
#        name = 'Z1_pt_superCoarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
#        attribute = TreeVariable.fromString( "Z1_pt/F" ),
#        binning=[3,0,600],
#    ))
#
#    plots.append(Plot(
#        name = "M3l",
#        texX = 'M(3l) (GeV)', texY = 'Number of Events',
#        attribute = lambda event, sample:event.M3l,
#        binning=[25,0,500],
#    ))
#
#    plots.append(Plot(
#        name = "dPhiZJet",
#        texX = '#Delta#phi(Z,j1)', texY = 'Number of Events',
#        attribute = lambda event, sample: deltaPhi(event.Z1_phi, event.JetGood_phi[0]),
#        binning=[20,0,pi],
#    ))
#
#    plots.append(Plot(
#        name = "l1_Z1_pt",
#        texX = 'p_{T}(l_{1,Z}) (GeV)', texY = 'Number of Events / 10 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.Z1_l1_index],
#        binning=[30,0,300],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = "l1_Z1_pt_coarse",
#        texX = 'p_{T}(l_{1,Z}) (GeV)', texY = 'Number of Events / 40 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.Z1_l1_index],
#        binning=[10,0,400],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'l1_Z1_pt_ext', texX = 'p_{T}(l_{1,Z}) (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.Z1_l1_index],
#        binning=[20,40,440],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = "l2_Z1_pt",
#        texX = 'p_{T}(l_{2,Z}) (GeV)', texY = 'Number of Events / 10 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.Z1_l2_index],
#        binning=[20,0,200],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#      texX = 'p_{T}(leading l) (GeV)', texY = 'Number of Events / 20 GeV',
#      name = 'lep1_pt', attribute = lambda event, sample: event.lep_pt[0],
#      binning=[400/20,0,400],
#    ))
#
#    plots.append(Plot(
#      texX = 'p_{T}(subleading l) (GeV)', texY = 'Number of Events / 10 GeV',
#      name = 'lep2_pt', attribute = lambda event, sample: event.lep_pt[1],
#      binning=[200/10,0,200],
#    ))
#
#    plots.append(Plot(
#      texX = 'p_{T}(trailing l) (GeV)', texY = 'Number of Events / 10 GeV',
#      name = 'lep3_pt', attribute = lambda event, sample: event.lep_pt[2],
#      binning=[150/10,0,150],
#    ))
#
#    plots.append(Plot(
#        name = "l2_Z1_pt_coarse",
#        texX = 'p_{T}(l_{2,Z}) (GeV)', texY = 'Number of Events / 10 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.Z1_l2_index],
#        binning=[10,0,200],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'l2_Z1_pt_ext', texX = 'p_{T}(l_{2,Z}) (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.Z1_l2_index],
#        binning=[20,0,400],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'lnonZ1_pt',
#        texX = 'p_{T}(l_{1,extra}) (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.nonZ1_l1_index],
#        binning=[15,0,300],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'lnonZ1_pt_coarse',
#        texX = 'p_{T}(l_{1,extra}) (GeV)', texY = 'Number of Events / 60 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.nonZ1_l1_index],
#        binning=[3,0,180],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'lnonZ1_charge',
#        texX = 'Charge(l_{1,extra})', texY = 'Number of Events',
#        attribute = lambda event, sample:-event.lep_pdgId[event.nonZ1_l1_index]/abs(event.lep_pdgId[event.nonZ1_l1_index]),
#        binning=[2,-1,1],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'lnonZ1_eta',
#        texX = '#eta(l_{1,extra})', texY = 'Number of Events',
#        attribute = lambda event, sample: event.lep_eta[event.nonZ1_l1_index],
#        binning=[20,-3,3],
#    ))
#
#    plots.append(Plot(
#        texX = 'M(ll) (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = TreeVariable.fromString( "Z1_mass/F" ),
#        binning=[10,81,101],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = "Z1_mass_wide",
#        texX = 'M(ll) (GeV)', texY = 'Number of Events / 2 GeV',
#        attribute = TreeVariable.fromString( "Z1_mass/F" ),
#        binning=[50,20,120],
#        addOverFlowBin='upper',
#    )) 
#
#    plots.append(Plot(
#        name = "Z1_cosThetaStar", texX = 'cos#theta(l-)', texY = 'Number of Events / 0.2',
#        attribute = lambda event, sample:event.Z1_cosThetaStar,
#        binning=[10,-1,1],
#    ))
#
#    plots.append(Plot(
#        name = "Z2_mass_wide",
#        texX = 'M(ll) of 2nd OSDL pair', texY = 'Number of Events / 2 GeV',
#        attribute = TreeVariable.fromString( "Z2_mass/F" ),
#        binning=[60,0,120],
#        addOverFlowBin='upper',
#    )) 
#
#    plots.append(Plot(
#        name = "minDLmass",
#        texX = 'min mass of all DL pairs', texY = 'Number of Events / 2 GeV',
#        attribute = TreeVariable.fromString( "minDLmass/F" ),
#        binning=[60,0,120],
#        addOverFlowBin='upper',
#    )) 
#
#    plots.append(Plot(
#        texX = '#Delta#phi(Z_{1}(ll))', texY = 'Number of Events',
#        attribute = TreeVariable.fromString( "Z1_lldPhi/F" ),
#        binning=[10,0,pi],
#    ))
#
#    plots.append(Plot(
#        texX = '#Delta R(Z_{1}(ll))', texY = 'Number of Events',
#        attribute = TreeVariable.fromString( "Z1_lldR/F" ),
#        binning=[10,0,6],
#    ))
#
#    plots.append(Plot(
#      texX = 'N_{jets}', texY = 'Number of Events',
#      attribute = TreeVariable.fromString( "nJetGood/I" ), #nJetSelected
#      binning=[8,-0.5,7.5],
#    ))
#
#    plots.append(Plot(
#      texX = 'N_{b-tag}', texY = 'Number of Events',
#      attribute = TreeVariable.fromString( "nBTag/I" ), #nJetSelected
#      binning=[4,-0.5,3.5],
#    ))
#
#    plots.append(Plot(
#      texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events / 30 GeV',
#      name = 'jet0_pt', attribute = lambda event, sample: event.JetGood_pt[0],
#      binning=[600/30,0,600],
#    ))
#
#    plots.append(Plot(
#      texX = 'p_{T}(subleading jet) (GeV)', texY = 'Number of Events / 30 GeV',
#      name = 'jet1_pt', attribute = lambda event, sample: event.JetGood_pt[1],
#      binning=[600/30,0,600],
#    ))
#
#    plots.append(Plot(
#        name = "W_pt",
#        texX = 'p_{T}(W) (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = lambda event, sample:event.W_pt,
#        binning=[20,0,400],
#    ))
#
    # 3l training variables
#
#    plots.append(Plot(
#      texX = '#Delta\#phi(nonZ-l_{1}, Z_{1})', texY = 'Number of Events',
#      name = 'nonZ1l1_Z1_deltaPhi', attribute = lambda event, sample: event.nonZ1l1_Z1_deltaPhi,
#      binning=[20,0,pi],
#    ))
#    plots.append(Plot(
#      texX = '#Delta#eta(nonZ-l_{1}, Z_{1})', texY = 'Number of Events',
#      name = 'nonZ1l1_Z1_deltaEta', attribute = lambda event, sample: event.nonZ1l1_Z1_deltaEta,
#      binning=[20,0,6],
#    ))
#    plots.append(Plot(
#      texX = '#Delta R(nonZ-l_{1}, Z_{1})', texY = 'Number of Event',
#      name = 'nonZ1l1_Z1_deltaR', attribute = lambda event, sample: event.nonZ1l1_Z1_deltaR,
#      binning=[20,0,6],
#    ))
#
#    plots.append(Plot(
#      texX = '#Delta R(jet_{0}, Z_{1})', texY = 'Number of Events',
#      name = 'jet0_Z1_deltaR', attribute = lambda event, sample: event.jet0_Z1_deltaR,
#      binning=[20,0,6],
#    ))
#    plots.append(Plot(
#      texX = '#Delta R(jet_{0}, nonZ-l_{1})', texY = 'Number of Events',
#      name = 'jet0_nonZ1l1_deltaR', attribute = lambda event, sample: event.jet0_nonZ1l1_deltaR,
#      binning=[20,0,6],
#    ))
#    plots.append(Plot(
#      texX = '#Delta R(jet_{1}, Z_{1})', texY = 'Number of Events',
#      name = 'jet1_Z1_deltaR', attribute = lambda event, sample: event.jet1_Z1_deltaR,
#      binning=[20,0,6],
#    ))
#    plots.append(Plot(
#      texX = '#Delta R(jet_{1}, nonZ-l_{1})', texY = 'Number of Events',
#      name = 'jet1_nonZ1l1', attribute = lambda event, sample: event.jet1_nonZ1l1_deltaR,
#      binning=[20,0,6],
#    ))
#    
#    for index in range(3):
#        for abs_pdg in [11, 13]:
#            lep_name = "mu" if abs_pdg==13 else "ele"
#            plots.append(Plot(
#              texX = 'p_{T}(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_pt'%(lep_name, index), attribute = lep_getter("pt", index, abs_pdg),
#              binning=[400/20,0,400],
#            ))
#            plots.append(Plot(
#              texX = '#eta(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_eta'%(lep_name, index), attribute = lep_getter("eta", index, abs_pdg),
#              binning=[30,-3,3],
#            ))
#            plots.append(Plot(
#              texX = '#phi(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_phi'%(lep_name, index), attribute = lep_getter("phi", index, abs_pdg),
#              binning=[30,-pi,pi],
#            ))
#            plots.append(Plot(
#              texX = 'dxy(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_dxy'%(lep_name, index), attribute = lep_getter("dxy", index, abs_pdg, functor = lambda x: abs(x)),
#              binning=[50,0,0.05],
#            ))
#            plots.append(Plot(
#              texX = 'dz(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_dz'%(lep_name, index), attribute = lep_getter("dz", index, abs_pdg, functor = lambda x: abs(x)),
#              binning=[50,0,0.05],
#            ))
#            plots.append(Plot(
#              texX = 'IP_{3D}(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_ip3d'%(lep_name, index), attribute = lep_getter("ip3d", index, abs_pdg, functor = lambda x: abs(x)),
#              binning=[50,0,0.05],
#            ))
#            plots.append(Plot(
#              texX = '#sigma(IP)_{3D}(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_sip3d'%(lep_name, index), attribute = lep_getter("sip3d", index, abs_pdg, functor = lambda x: abs(x)),
#              binning=[40,0,8],
#            ))
#            plots.append(Plot(
#              texX = 'jetRelIso(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_jetRelIso'%(lep_name, index), attribute = lep_getter("jetRelIso", index, abs_pdg),
#              binning=[50,-.15,0.5],
#            ))
#            plots.append(Plot(
#              texX = 'miniPFRelIso_all(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_miniPFRelIso_all'%(lep_name, index), attribute = lep_getter("miniPFRelIso_all", index, abs_pdg),
#              binning=[50,0,.5],
#            ))
#            plots.append(Plot(
#              texX = 'pfRelIso03_all(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_pfRelIso03_all'%(lep_name, index), attribute = lep_getter("pfRelIso03_all", index, abs_pdg),
#              binning=[50,0,.5],
#            ))
#            plots.append(Plot(
#              texX = 'mvaTTH(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_mvaTTH'%(lep_name, index), attribute = lep_getter("mvaTTH", index, abs_pdg),
#              binning=[24,-1.2,1.2],
#            ))
#            plots.append(Plot(
#              texX = 'charge(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_charge'%(lep_name, index), attribute = lep_getter("pdgId", index, abs_pdg, functor = charge),
#              binning=[3,-1,2],
#            ))
#            if lep_name == "mu":
#                plots.append(Plot(
#                  texX = 'segmentComp(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#                  name = '%s%i_segmentComp'%(lep_name, index), attribute = lep_getter("segmentComp", index, abs_pdg),
#                  binning=[50,0,1],
#                ))
#                plots.append(Plot(
#                  texX = 'nStations(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#                  name = '%s%i_nStations'%(lep_name, index), attribute = lep_getter("nStations", index, abs_pdg),
#                  binning=[10,0,10],
#                ))
#                plots.append(Plot(
#                  texX = 'nTrackerLayers(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#                  name = '%s%i_nTrackerLayers'%(lep_name, index), attribute = lep_getter("nTrackerLayers", index, abs_pdg),
#                  binning=[20,0,20],
#                ))
#            if lep_name == "ele":
#                for cbIdFlag in vidNestedWPBitMapNamingList:
#                    plots.append(Plot(
#                      texX = '%s(%s_{%i}) (GeV)'%(cbIdFlag, lep_name, index), texY = 'Number of Events',
#                      name = '%s%i_%s_Flag'%(lep_name, index, cbIdFlag), attribute = lep_getter("vidNestedWPBitmap", index, abs_pdg, functor = cbEleIdFlagGetter(cbIdFlag)),
#                      binning=[5,0,5],
#                    ))

    plotting.fill(plots, read_variables = read_variables, sequence = sequence)

    # Get normalization yields from yield histogram
    for plot in plots:
      if plot.name == "yield":
        for i, l in enumerate(plot.histos):
          for j, h in enumerate(l):
            yields[mode][plot.stack[i][j].name] = h.GetBinContent(h.FindBin(0.5+i_mode))
            h.GetXaxis().SetBinLabel(1, "#mu#mu#mu")
            h.GetXaxis().SetBinLabel(2, "#mu#mue")
            h.GetXaxis().SetBinLabel(3, "#muee")
            h.GetXaxis().SetBinLabel(4, "eee")
      if plot.name.endswith("_Flag"):
        for i, l in enumerate(plot.histos):
          for j, h in enumerate(l):
            h.GetXaxis().SetBinLabel(1, "fail")
            h.GetXaxis().SetBinLabel(2, "veto")
            h.GetXaxis().SetBinLabel(3, "loose")
            h.GetXaxis().SetBinLabel(4, "medium")
            h.GetXaxis().SetBinLabel(5, "tight")
        
    if args.noData: yields[mode]["data"] = 0

    yields[mode]["MC"] = sum(yields[mode][s.name] for s in mc)
    dataMCScale        = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')

    drawPlots(plots, mode, dataMCScale)
    allPlots[mode] = plots


# Add the different channels into SF and all
for mode in ["comb1","comb2","all"]:
    yields[mode] = {}
    for y in yields[allModes[0]]:
        try:    yields[mode][y] = sum(yields[c][y] for c in ['eee','muee','mumue', 'mumumu'])
        except: yields[mode][y] = 0
    dataMCScale = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')
    
    for plot in allPlots['mumumu']:
        if mode=="comb1":
            tmp = allPlots['mumue']
        elif mode=="comb2":
            tmp = allPlots['muee']
        else:
            tmp = allPlots['eee']
        for plot2 in (p for p in tmp if p.name == plot.name):
            for i, j in enumerate(list(itertools.chain.from_iterable(plot.histos))):
                for k, l in enumerate(list(itertools.chain.from_iterable(plot2.histos))):
                    if i==k:
                        j.Add(l)
    
    if mode == "all": drawPlots(allPlots['mumumu'], mode, dataMCScale)

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
