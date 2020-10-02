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
from tWZ.Tools.objectSelection           import lepString

from Analysis.Tools.helpers              import deltaPhi, deltaR, getCollection, getObjDict
from Analysis.Tools.puProfileCache       import *
from Analysis.Tools.puReweighting        import getReweightingFunction
from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons
import Analysis.Tools.syncer 


# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--noData',         action='store_true', default=False, help='also plot data?')
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--dataMCScaling',  action='store_true', help='Data MC scaling?', )
argParser.add_argument('--plot_directory', action='store', default='tWZ_v3')
argParser.add_argument('--era',            action='store', type=str, default="Run2016")
argParser.add_argument('--selection',      action='store', default='trilepM-minDLmass12-onZ1-njet4p-btag2p')
#argParser.add_argument('--nanoAODv4',   default=True, action='store_true',                                                                        help="Run on nanoAODv4?" )
argParser.add_argument('--samples',        action='store',         nargs='*',  type=str, default=['TTZToLLNuNu_ext'],                  help="List of samples to be post-processed, given as CMG component name" )
#flagg for parton selection
argParser.add_argument('--partonweight',   action='store_true', help='weight partons?', )
argParser.add_argument('--partonweightprivate', action='store_true', help='weight partons?', )
argParser.add_argument('--normalize',          action='store_true', default=False,                                                                   help="Normalize to 1" )
argParser.add_argument('--privateTWZttZ',  action='store_true', help='twzDR and ttZ', )
argParser.add_argument('--scaled',         action='store_true', help='scaling', )
#argParser.add_argument('--mva',            action='store', type=str )

args = argParser.parse_args()
options = argParser.parse_args()

# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"
if args.noData:                       args.plot_directory += "_noData"
if args.partonweightprivate:          args.plot_directory += "_weightedpartons(private)"
if args.privateTWZttZ:                args.plot_directory += "_privateTWZDR_ttZ"
if args.normalize:                    args.plot_directory += "_normalize"
if args.scaled:                       args.plot_directory += "_scaled"
logger.info( "Working in era %s", args.era)

from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed import *

#priavate tWZ match samples
#TWZ match samples 
if args.partonweightprivate: 
    tWZ_ud_match = copy.deepcopy( Summer16.TWZ_NLO_DR  )
    tWZ_ud_match.name = "tWZ_ud_match"
    tWZ_ud_match.texName = "private tWZ (fwd u/d)"
    
    tWZ_gluon_match = copy.deepcopy( Summer16.TWZ_NLO_DR  )
    tWZ_gluon_match.name = "tWZ_gluon_match"
    tWZ_gluon_match.texName = "private tWZ (fwd gluon)"
    
    tWZ_other_match = copy.deepcopy( Summer16.TWZ_NLO_DR  )
    tWZ_other_match.name = "tWZ_other_match"
    tWZ_other_match.texName = "private tWZ (fwd others)"


if args.era == "Run2016":
    if args.partonweightprivate: 
        mc = [tWZ_gluon_match, tWZ_ud_match, tWZ_other_match, Summer16.TWZ_NLO_DR, Summer16.TTZ]
    elif args.privateTWZttZ: 
        mc = [Summer16.TWZ_NLO_DR, Summer16.TTZ]
    else: 
        mc = [Summer16.TTZ, Summer16.TTX_rare, Summer16.TZQ, Summer16.WZ, Summer16.triBoson, Summer16.ZZ, Summer16.nonprompt_3l, Summer16.TWZ_NLO_DR]
elif args.era == "Run2017":
    mc = [Fall17.TWZ_NLO_DR, Fall17.TTZ, Fall17.TTX_rare, Fall17.TZQ, Fall17.WZ, Fall17.triBoson, Fall17.ZZ, Fall17.nonprompt_3l]
elif args.era == "Run2018":
    mc = [Autumn18.TWZ_NLO_DR, Autumn18.TTZ, Autumn18.TTX_rare, Autumn18.TZQ, Autumn18.WZ, Autumn18.triBoson, Autumn18.ZZ, Autumn18.nonprompt_3l]
elif args.era == "RunII":
    mc = [TWZ_NLO_DR, TTZ, TTX_rare, TZQ, WZ, triBoson, ZZ, nonprompt_3l]

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
    for sample in mc:# + [data_sample]:
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
        
      if args.normalize:
        for comp in plot.histos:
            for h in comp:
                if h.Integral()!=0.: h.Scale(1./h.Integral())  
      
      if isinstance( plot, Plot):
            plotting.draw(plot,
            plot_directory = plot_directory_,
            ratio = {'yRange':(0.1,1.9)} if not args.noData else None,
            #ratio = {'yRange': (0.1, 1.9), 'histos':[(1,0)], 'texY':'Ratio'} if args.scaled else None,
            #ratio = {'yRange': (0.1, 1.9), 'histos':[(1,0),(2,0),(3,0),(4,0)], 'texY':'Ratio'} if args.scaled else None,
            logX = False, logY = log, sorting = True,
            yRange = (0.03, "auto") if log else (0.001, "auto"),
            scaling =  { 1:0 } if args.scaled else {}, 
            #scaling =  { i+1:0 for i in range(4) } if args.scaled else {}, 
            legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
            drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ) + _drawObjects,
            copyIndexPHP = True, extensions = ["png"],
          )

# Read variables and sequences

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

#jets without eta cut
jetVars          = ['pt/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'rawFactor/F', 'eta/F', 'phi/F', 'jetId/I', 'btagDeepB/F', 'btagDeepFlavB/F', 'btagCSVV2/F', 'area/F'] 
jetVarNames      = [x.split('/')[0] for x in jetVars]

lepVars         = ['pt/F','eta/F','phi/F','pdgId/I','muIndex/I','eleIndex/I','cutBased/I','miniPFRelIso_all/F','pfRelIso03_all/F','mvaFall17V2Iso_WP90/O', 'sip3d/F','lostHits/I','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','eleIndex/I','muIndex/I']
lepVarNames     = [x.split('/')[0] for x in lepVars]

bVars      = ['pt/F','eta/F','phi/F']
bVarNames  = [x.split('/')[0] for x in bVars]

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I", 
    "l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
    "l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
    "JetGood[pt/F,eta/F,phi/F]", 
    "nJet/I", 
    "nlep/I", "lep[%s]"%(",".join(lepVars)),
    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I", 
    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,pdgId/I,vidNestedWPBitmap/I]",
]

#b tagger
b_tagger = "DeepJet"

# read only for data:
read_variables_data = [ "Jet[%s]"%(",".join(jetVars))] 

# read only for MC:
read_variables_MC = ['reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F']
read_variables_MC.append( "Jet[%s,genJetIdx/I]"%(",".join(jetVars)) )
read_variables_MC.append( "GenJet[pt/F,eta/F,phi/F,hadronFlavour/b,partonFlavour/I]")

# define 3l selections

#def getbjets( event, sample ):
#    #jets einlesen (in case of MC also reat the index of the genjet)
#    alljets     = getCollection( event, 'Jet', jetVarNames + (['genJetIdx'] if not sample.isData else []), 'nJet')
#    event.jets_no_eta         = filter(lambda j:j['pt']>30, clean_jets)
#    #bjets filtern, nur 2016
#    #print jetgood
#    bJets = filter(lambda j:isBJet(j, tagger=b_tagger, year=2016) and abs(j['eta'])<=2.4    , alljets)
#    #print bJets
#    for i in range(len(bJets)) :
#        event.bjet_pt  = bJets[i]['pt']
#        event.bjet_eta = bJets[i]['eta']
#        event.bjet_phi = bJets[i]['phi']
# 
#    event.bjet_Z1_deltaR      = deltaR({'eta':event.bjet_eta, 'phi':event.bjet_phi}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
#    event.bjet_nonZ1l1_deltaR = deltaR({'eta':event.bjet_eta, 'phi':event.bjet_phi}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
#
#sequence.append( getbjets )

def getjets( event, sample ):
    #jets einlesen (in case of MC also reat the index of the genjet)
    alljets   = getCollection( event, 'Jet', jetVarNames + (['genJetIdx'] if not sample.isData else []), 'nJet')  
    alljets.sort( key = lambda j: -j['pt'] )
    leptons    = getCollection(event, "lep", lepVarNames, 'nlep') 
    # clean against good leptons
    clean_jets,_ = cleanJetsAndLeptons( alljets, leptons )
    jets  = filter(lambda j:j['pt']>30, clean_jets)
    # filter pt, but not eta (I store the list of jets in the event because I want to use it in the next function)
    event.jets_no_eta         = filter(lambda j:j['pt']>30, clean_jets)
    # very nice, Rosmarie. Here are all the python built on functions, so you get an idea what you can do: https://docs.python.org/2.7/library/functions.html
    #print (event.jets_no_eta)
    if event.jets_no_eta: 
        event.maxEta_of_pt30jets  = max( [ abs(j['eta']) for j in event.jets_no_eta ] )
        #print (event.maxEta_of_pt30jets)
#    #bjets filtern( nur 2016 )
#    bJets = filter(lambda j:isBJet(j, tagger=b_tagger, year=2016) and abs(j['eta'])<=2.4    , jets)
#    # print bJets
#    #print bJets
#    for i in range(len(bJets)) :
#        event.bjet_pt  = bJets[i]['pt']
#        event.bjet_eta = bJets[i]['eta']
#        event.bjet_phi = bJets[i]['phi']
#    event.bjet_maxbtagdeepb  =  max( [ j['btagDeepB'] for j in bJets ] ) 
#    
#   # event.bjet_2maxbtagdeepb =  max( [ j['btagDeepB'] for j in bJets ] )
#
#    event.bjet_Z1_deltaR      = deltaR({'eta':event.bjet_eta, 'phi':event.bjet_phi}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
#    event.bjet_nonZ1l1_deltaR = deltaR({'eta':event.bjet_eta, 'phi':event.bjet_phi}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
#
sequence.append( getjets )

def genJetStuff( event, sample ):
    # only do something in simulation
    if sample.isData:
        return
    # let's now match jets with genjets. It's been done for us in the nanoAOD. 
    # here you have all nanoAOD branches: https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html#Muon    
    # Jet_genJetIdx is the index of the genjet that matches the jet. It only makes sense in simulation. 
    # See above and below how I changed the logic of when we read this collection (look at the read_variables_MC/data).

    # If you look up the max function you learn that you can give it a method that tells python how it should do the comparison. Let's use it to get the highest eta *jet* 
    if event.jets_no_eta:
        max_eta_jet = max( event.jets_no_eta, key = lambda j:abs(j['eta']) )
        #print (max_eta_jet)
        #print (max_eta_jet['genJetIdx'])
    # Let's see if it has a gen match:
    # Set default values
    event.partonsinfwdjets = -8
    if max_eta_jet['genJetIdx']>=0:
        # we have a gen jet index (negative number means nothing found). 
        # Let's not read all the genjets, just the one we need. Looking up getCollection you'll find it only loops getObjDict, i.e. a method that makes a dictionary from event.Collection_branch
        # the "prefix" argument of getObjDict needs an underscore ... we could have done that better. 
        max_eta_genjet = getObjDict( event, "GenJet_", ["pt", "eta", "phi", "hadronFlavour", "partonFlavour"], max_eta_jet['genJetIdx'] )
        # You'll find that GenJet_hadronFlavour is a UChar_t, i.e. to save disk space we're not using a 32 bit integer but an 8 bit character (it explains the /b above = Byte)
        # We change the byte to a normal number (the function ord is in the link above)
        #max_eta_genjet['hadronFlavour'] = ord(max_eta_genjet['hadronFlavour'])
        # Genjets are clustered from all stable particles after the shower (pythia) and hadronisation. Neutrinos are not taken into account. Look at the gen jet flavour: 
        # partonFlavor: This flavor is obtained by including the generator partons (from the hard scatter) in the jet clustering but with infinitesimal momentum (otherwise the genjet would be changed)
        #               The partonflavor is the pdgId of the genparton that gets clustered with the jet. E.g. 5 means b-jet. 
        #               Negative numbers mean antiparticles.
        #               The hadronFlavour is the same, but instead of partons, the generated hadrons are used (e.g. after shower+hadronisation). The most important difference is gluon splitting:
        #               For strongly interacting particles, we can radiate a gluon and gluon can produce a b/bbar pair. Such b quarks from the shower are not what we want to use when we decide
        #               whether a jet was originally a b-jet because they don't originate from the top but are rather produced inside jets. 
        #               Thus, we mostly use partonFlavour (i.e. hard scatter) while e.g. the performance measurements of b-tagging are done with the hadronFlavor (after all, it's a true b)
        # Look at the numbers. Here is the dictionary: http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
        #print max_eta_jet['genJetIdx'], max_eta_genjet
        # partonFlavour number 
        #print max_eta_genjet['partonFlavour']
        #in event schreiben 
        
        event.partonsinfwdjets =  max_eta_genjet['partonFlavour']
        event.ud_match    =  max_eta_genjet['partonFlavour'] in [ 1, 2, -1, -2]
        event.gluon_match =  max_eta_genjet['partonFlavour'] in [ 21 ]
        event.other_match =  max_eta_genjet['partonFlavour'] not in [ 1, 2, -1, -2, 21 ]
    else:
        event.ud_match    = 0
        event.gluon_match = 0
        event.other_match = 1

sequence.append( genJetStuff )

#do the same for min. eta
#def genjetmineta( event, sample ):    
#    if sample.isData:
#        return
#    min_eta_jet = min( event.jets_no_eta, key = lambda j:abs(j['eta']) )
#    event.partonsinfwdjetsmineta = -8
#    if min_eta_jet['genJetIdx']>=0:
#        min_eta_genjet = getObjDict( event, "GenJet_", ["pt", "eta", "phi", "hadronFlavour", "partonFlavour"], min_eta_jet['genJetIdx'] )
#        event.partonsinfwdjetsmineta =  min_eta_genjet['partonFlavour']
#
#sequence.append( genjetmineta )

#DeltaR 
#def deltaRfwdjeteta(event, sample):
#    #get highest eta jet(most fwd) 
#    event.max_eta_jet = max( event.jets_no_eta, key = lambda j:abs(j['eta']) )
#    #use deltaR function to get deltaR 
#    event.maxeta_Z_deltaR      = deltaR({'eta':event.max_eta_jet['eta'], 'phi':event.max_eta_jet['phi']}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
#    
#sequence.append( deltaRfwdjeteta )

#DeltaEta
def deltaEtaZll( event,sample ):
    event.Z1_lldEta =  event.lep_eta[event.Z1_l2_index] - event.lep_eta[event.Z1_l1_index]
#    print event.Z1_lldEta

sequence.append( deltaEtaZll )

mu_string  = lepString('mu','VL')
ele_string = lepString('ele','VL')
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

#MVA
#WZ 
from Analysis.TMVA.Reader    import Reader

from tWZ.MVA.MVA_TWZ_3l      import mva_variables #, all_mlp_np40, all_mlp_np5s0c3e0c5, all_mlp_ncnc1s0c3e0c5 
from tWZ.MVA.MVA_TWZ_3l      import sequence as mva_sequence
from tWZ.MVA.MVA_TWZ_3l      import read_variables as mva_read_variables

import tWZ.MVA.MVA_TWZ_3l as tWZ_vs_ttZ
import tWZ.MVA.MVA_TWZ_WZ_3l as tWZ_vs_WZ

mvas = [tWZ_vs_ttZ.all_mlp_ncnc1s0c3e0c5, tWZ_vs_WZ.wz_all_mlp_ncnc1s0c3e0c5, tWZ_vs_WZ.wz_all_mlp_np5s0c3e0c5, tWZ_vs_WZ.wz_all_mlp_np40 ]
wzmvas = [tWZ_vs_WZ.wz_all_mlp_ncnc1s0c3e0c5, tWZ_vs_WZ.wz_all_mlp_np5s0c3e0c5, tWZ_vs_WZ.wz_all_mlp_np40 ]

from tWZ.Tools.user          import mva_directory

sequence.extend( mva_sequence )
read_variables.extend( mva_read_variables )

mva_reader = Reader(
    mva_variables     = mva_variables,
    weight_directory  = os.path.join( mva_directory, "TWZ_3l" ),
    label             = "TWZ_3l")

def makeDiscriminator( mva ):
    def _getDiscriminator( event, sample ):
        kwargs = {name:func(event,None) for name, func in mva_variables.iteritems()}
        setattr( event, mva['name'], mva_reader.evaluate(mva['name'], **kwargs))
    return _getDiscriminator

def discriminator_getter(name):
    def _disc_getter( event, sample ):
        return getattr( event, name )
    return _disc_getter

for mva in mvas:
    mva_reader.addMethod(method=mva)
    sequence.append( makeDiscriminator(mva) )

def get_wzmva( threshold, weightmva ):
    def wzmva(event, sample):
        mva = weightmva   
        passweight = discriminator_getter(mva['name'])(event, sample)
        event.wzpassweight = int(passweight > threshold)
        #print event.wzpassweight  
    return wzmva

thresholds = {'cut0p85':0.85, 'cut0p9':0.9, 'cut0p95':0.95}

for wzmva in wzmvas: 
    for cut in thresholds: 
        sequence.append( get_wzmva(threshold=thresholds[cut], weightmva=wzmva) ) 
        print ('passweight'+wzmva['name']+cut) 
        #event. 
        passweight =  lambda event, sample: event.wzpassweight #if sample.isData else event.wzpassweight*lumi_year[event.year]/1000.

weight_ = lambda event, sample: event.weight if sample.isData else event.weight*lumi_year[event.year]/1000.

#weight_wzmva_pass_ = lambda event, sample: event.wzpassweight #if sample.isData else event.wzpassweight*lumi_year[event.year]/1000.


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
        data_sample.read_variables = read_variables_data
        lumi_scale                 = data_sample.lumi/1000

    for sample in mc: sample.style = styles.fillStyle(sample.color)
    #for sample in mc: sample.style = styles.lineStyle(sample.color)
    
    if args.partonweightprivate: 
        tWZ_ud_match.style    = styles.lineStyle(ROOT.kBlue)
        tWZ_gluon_match.style = styles.lineStyle(ROOT.kYellow)
        tWZ_other_match.style = styles.lineStyle(ROOT.kGreen)
        TTZ.style             = styles.lineStyle(1)       
        TWZ_NLO_DR.style      = styles.lineStyle(ROOT.kCyan)

    for sample in mc:
      sample.read_variables = read_variables_MC 
      sample.setSelectionString([getLeptonSelection(mode)])
      # if you are going to use sample e.g. in sample.name, 'sample' must be the loop iterator or otherwise defined
      #if   sample.name == "tWZ_ud_match"    : sample.weight = lambda event, sample: event.ud_match*event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger
      #elif sample.name == "tWZ_gluon_match" : sample.weight = lambda event, sample: event.gluon_match*event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger
      #elif sample.name == "tWZ_other_match" : sample.weight = lambda event, sample: event.other_match*event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger
      #else: sample.weight = lambda event, sample: event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger

    #yt_TWZ_filter.scale = lumi_scale * 1.07314
    
    if args.partonweightprivate:
        stack = Stack([Summer16.TWZ_NLO_DR],[tWZ_ud_match],[tWZ_gluon_match],[tWZ_other_match],[Summer16.TTZ])
    if args.scaled: 
        stack = Stack()
    elif args.privateTWZttZ:
        mc = ([Summer16.TWZ_NLO_DR], [Summer16.TTZ])
    else:
        if not args.noData:
          stack = Stack(mc, data_sample)
        else:
          stack = Stack(mc)
    
    for sample in mc+[data_sample]:
        sample.weight = weight_

    # Use some defaults
    Plot.setDefaults(stack = stack, weight = None, selectionString = cutInterpreter.cutString(args.selection))

    plots = []

    for mva in mvas:
        for cut in thresholds: 
            plots.append(Plot(
                texX = 'MVA_{3l}'+'_'+cut, texY = 'Number of Events',
                name = mva['name']+'_'+cut, attribute = discriminator_getter(mva['name']),
                binning=[25, 0, 1],
                weight = passweight , 
            ))
            
            plots.append(Plot(
                texX = 'MVA_{3l}'+'_'+cut, texY = 'Number of Events',
                name = mva['name']+'_'+cut+'_coarse', attribute = discriminator_getter(mva['name']),
                binning=[10, 0, 1],
                weight = passweight, 
            ))

    for mva in mvas:
        plots.append(Plot(
            texX = 'MVA{3l}', texY = 'Number of Events',
            name = mva['name'], attribute = discriminator_getter(mva['name']),
            binning=[25, 0, 1],
        ))

    for mva in mvas:
        plots.append(Plot(
            texX = 'MVA_{3l}', texY = 'Number of Events',
            name = mva['name']+'_coarse', attribute = discriminator_getter(mva['name']),
            binning=[10, 0, 1],
        ))

#    for mva in mvas:
#        plots.append(Plot(
#            texX = 'MVA_pass{3l}', texY = 'Number of Events',
#            name = mva['name']+'_WZ0p95', attribute = discriminator_getter(mva['name']),
#            binning=[25, 0, 1],
#            weight =  weight_wzmva_pass_ , 
#        ))
#
#    for mva in mvas:
#        plots.append(Plot(
#            texX = 'MVA_pass{3l}', texY = 'Number of Events',
#            name = mva['name']+'_WZ0p95_coarse', attribute = discriminator_getter(mva['name']),
#            binning=[10, 0, 1],
#            weight = weight_wzmva_pass_ , 
#        ))

#    plots.append(Plot(
#      name = 'maxpt',
#      texX = 'maxpt',
#      texY = 'Number of Events',
#      attribute = lambda event, sample: event.max_pt_jet, 
#      binning=[20, 0, 400],
#    ))
#
#    plots.append(Plot(
#      name = 'maxabseta',
#      texX = 'abs(#eta)_max',
#      texY = 'Number of Events',
#      attribute = lambda event, sample: event.maxEta_of_pt30jets,
#      binning=[20, 0, 5],
#    ))
#
#    plots.append(Plot(
#      name = 'partons in jets mineta',
#      texX = 'partons in jets (mineta)',
#      texY = 'Number of Events',
#      attribute = lambda event, sample: event.partonsinfwdjetsmineta,
#      binning=[28, -6, 22],
#    ))
#
#    plots.append(Plot(
#      name = 'partons in fwd jets maxeta',
#      texX = 'partons in fwd jets (maxeta)',
#      texY = 'Number of Events',
#      attribute = lambda event, sample: event.partonsinfwdjets,
#      binning=[28, -6, 22],
#    ))

    plots.append(Plot(
      name = 'yield', texX = '', texY = 'Number of Events',
      attribute = lambda event, sample: 0.5 + i_mode,
      binning=[4, 0, 4],
    ))

#    plots.append(Plot(
#      texX = '#Delta R(maxetajet, Z_{1})', texY = 'Number of Events',
#      name = 'maxptjet_Z1_deltaR', attribute = lambda event, sample: event.maxeta_Z_deltaR,
#      binning=[20,0,6],
#    ))
#
#    plots.append(Plot(
#        name = "M2l",
#        texX = 'M(2l) (GeV)', texY = 'Number of Events',
#        attribute = lambda event, sample:event.M2l,
#        binning=[25,0,500],
#    ))
#
#    plots.append(Plot(
#        name = "m2l",
#        texX = 'M(2l) (GeV)', texY = 'Number of Events',
#        attribute = lambda event, sample:event.M2l,
#        binning=[25,50,150],
#    ))
#
#    plots.append(Plot(
#        name = "M(2l)",
#        texX = 'M(2l) (GeV)', texY = 'Number of Events',
#        attribute = lambda event, sample:event.M2l,
#        binning=[25,80,100],
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
#        name = 'Deltaeta_Z1_ll',
#        texX = '#Delta#eta(Z_{1}(ll))', texY = 'Number of Events',
#        attribute = lambda event, sample: event.Z1_lldEta,
#        binning=[10,-3,3],
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

   #bjets 
#
#    plots.append(Plot(
#      texX = 'p_{T}(b jet) (GeV)', texY = 'Number of Events / 30 GeV',
#      name = 'bjet_pt', attribute = lambda event, sample: event.bjet_pt,
#      binning=[600/30,0,600],
#    ))
#
#    
#    plots.append(Plot(
#      texX = '#Delta R(bjet, Z_{1})', texY = 'Number of Events',
#      name = 'bjet_Z1_deltaR', attribute = lambda event, sample: event.bjet_Z1_deltaR,
#      binning=[20,0,6],
#    ))
#
#    plots.append(Plot(
#      texX = '#Delta R(bjet, nonZ-l_{1})', texY = 'Number of Events',
#      name = 'bjet_nonZ1l1_deltaR', attribute = lambda event, sample: event.bjet_nonZ1l1_deltaR,
#      binning=[20,0,6],
#    ))
#
#    plots.append(Plot(
#      name = 'btagDeepB',
#      texX = 'btagDeepB',
#      texY = 'Number of Events',
#      attribute = lambda event, sample: event.bjet_maxbtagdeepb,
#      binning=[20, 0, 5],
#    ))
#
   # 3l training variables

#    plots.append(Plot(
#      texX = '#Delta\#phi(nonZ-l_{1}, Z_{1})', texY = 'Number of Events',
#      name = 'nonZ1l1_Z1_deltaPhi', attribute = lambda event, sample: event.nonZ1l1_Z1_deltaPhi,
#      binning=[20,0,pi],
#    ))
#
#    plots.append(Plot(
#      texX = '#Delta#eta(nonZ-l_{1}, Z_{1})', texY = 'Number of Events',
#      name = 'nonZ1l1_Z1_deltaEta', attribute = lambda event, sample: event.nonZ1l1_Z1_deltaEta,
#      binning=[20,0,6],
#    ))
#
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
#
#    plots.append(Plot(
#      texX = '#Delta R(jet_{0}, nonZ-l_{1})', texY = 'Number of Events',
#      name = 'jet0_nonZ1l1_deltaR', attribute = lambda event, sample: event.jet0_nonZ1l1_deltaR,
#      binning=[20,0,6],
#    ))
#
#    plots.append(Plot(
#      texX = '#Delta R(jet_{1}, Z_{1})', texY = 'Number of Events',
#      name = 'jet1_Z1_deltaR', attribute = lambda event, sample: event.jet1_Z1_deltaR,
#      binning=[20,0,6],
#    ))
#
#    plots.append(Plot(
#      texX = '#Delta R(jet_{1}, nonZ-l_{1})', texY = 'Number of Events',
#      name = 'jet1_nonZ1l1', attribute = lambda event, sample: event.jet1_nonZ1l1_deltaR,
#      binning=[20,0,6],
#    ))
#    
  #  for index in range(3):
  #      for abs_pdg in [11, 13]:
  #          lep_name = "mu" if abs_pdg==13 else "ele"
  #          plots.append(Plot(
  #            texX = 'p_{T}(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #            name = '%s%i_pt'%(lep_name, index), attribute = lep_getter("pt", index, abs_pdg),
  #            binning=[400/20,0,400],
  #          ))
  #          plots.append(Plot(
  #            texX = '#eta(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #            name = '%s%i_eta'%(lep_name, index), attribute = lep_getter("eta", index, abs_pdg),
  #            binning=[30,-3,3],
  #          ))
  #          plots.append(Plot(
  #            texX = '#phi(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #            name = '%s%i_phi'%(lep_name, index), attribute = lep_getter("phi", index, abs_pdg),
  #            binning=[30,-pi,pi],
  #          ))
  #          plots.append(Plot(
  #            texX = 'dxy(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #            name = '%s%i_dxy'%(lep_name, index), attribute = lep_getter("dxy", index, abs_pdg, functor = lambda x: abs(x)),
  #            binning=[50,0,0.05],
  #          ))
  #          plots.append(Plot(
  #            texX = 'dz(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #            name = '%s%i_dz'%(lep_name, index), attribute = lep_getter("dz", index, abs_pdg, functor = lambda x: abs(x)),
  #            binning=[50,0,0.05],
  #          ))
  #          plots.append(Plot(
  #            texX = 'IP_{3D}(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #            name = '%s%i_ip3d'%(lep_name, index), attribute = lep_getter("ip3d", index, abs_pdg, functor = lambda x: abs(x)),
  #            binning=[50,0,0.05],
  #          ))
  #          plots.append(Plot(
  #            texX = '#sigma(IP)_{3D}(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #            name = '%s%i_sip3d'%(lep_name, index), attribute = lep_getter("sip3d", index, abs_pdg, functor = lambda x: abs(x)),
  #            binning=[40,0,8],
  #          ))
  #          plots.append(Plot(
  #            texX = 'jetRelIso(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #            name = '%s%i_jetRelIso'%(lep_name, index), attribute = lep_getter("jetRelIso", index, abs_pdg),
  #            binning=[50,-.15,0.5],
  #          ))
  #          plots.append(Plot(
  #            texX = 'miniPFRelIso_all(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #            name = '%s%i_miniPFRelIso_all'%(lep_name, index), attribute = lep_getter("miniPFRelIso_all", index, abs_pdg),
  #            binning=[50,0,.5],
  #          ))
  #          plots.append(Plot(
  #            texX = 'pfRelIso03_all(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #            name = '%s%i_pfRelIso03_all'%(lep_name, index), attribute = lep_getter("pfRelIso03_all", index, abs_pdg),
  #            binning=[50,0,.5],
  #          ))
  #          plots.append(Plot(
  #            texX = 'mvaTTH(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #            name = '%s%i_mvaTTH'%(lep_name, index), attribute = lep_getter("mvaTTH", index, abs_pdg),
  #            binning=[24,-1.2,1.2],
  #          ))
  #          plots.append(Plot(
  #            texX = 'charge(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #            name = '%s%i_charge'%(lep_name, index), attribute = lep_getter("pdgId", index, abs_pdg, functor = charge),
  #            binning=[3,-1,2],
  #          ))
  #          if lep_name == "mu":
  #              plots.append(Plot(
  #                texX = 'segmentComp(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #                name = '%s%i_segmentComp'%(lep_name, index), attribute = lep_getter("segmentComp", index, abs_pdg),
  #                binning=[50,0,1],
  #              ))
  #              plots.append(Plot(
  #                texX = 'nStations(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #                name = '%s%i_nStations'%(lep_name, index), attribute = lep_getter("nStations", index, abs_pdg),
  #                binning=[10,0,10],
  #              ))
  #              plots.append(Plot(
  #                texX = 'nTrackerLayers(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
  #                name = '%s%i_nTrackerLayers'%(lep_name, index), attribute = lep_getter("nTrackerLayers", index, abs_pdg),
  #                binning=[20,0,20],
  #              ))
  #          if lep_name == "ele":
  #              for cbIdFlag in vidNestedWPBitMapNamingList:
  #                  plots.append(Plot(
  #                    texX = '%s(%s_{%i}) (GeV)'%(cbIdFlag, lep_name, index), texY = 'Number of Events',
  #                    name = '%s%i_%s_Flag'%(lep_name, index, cbIdFlag), attribute = lep_getter("vidNestedWPBitmap", index, abs_pdg, functor = cbEleIdFlagGetter(cbIdFlag)),
  #                    binning=[5,0,5],
  #                  ))

#newplots 

    #plots.append(Plot(
    #   texX = 'mvaTOP(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
    #   name = '%s%i_mvaTOP'%(lep_name, index), attribute = lep_getter("mvaTOP", index, abs_pdg),
    #   binning=[24,-1.2,1.2],
    # ))

#    plots.append(Plot(
#       name = 'l1_pt',
#       texX = 'p_{T}(l_{1}) (GeV)', texY = 'Number of Events / 20 GeV',
#       attribute = lambda event, sample:event.l1_pt,
#       binning=[15,0,300],
#       addOverFlowBin='upper',
#   ))
#
#    plots.append(Plot(
#       name = 'l1_eta',
#       texX = '#eta(l_{1})', texY = 'Number of Events',
#       attribute = lambda event, sample: event.l1_eta,
#       binning=[20,-3,3],
#   ))
#
#    plots.append(Plot(
#       name = 'l1_mvaTOP',
#       texX = 'MVA_{TOP}(l_{1})', texY = 'Number of Events',
#       attribute = lambda event, sample: event.l1_mvaTOP,
#       binning=[20,-1,1],
#   ))
#
#    plots.append(Plot(
#       name = 'l1_mvaTOPWP',
#       texX = 'MVA_{TOP}(l_{1}) WP', texY = 'Number of Events',
#       attribute = lambda event, sample: event.l1_mvaTOPWP,
#       binning=[5,0,5],
#   ))
#
#    plots.append(Plot(
#       name = 'l2_pt',
#       texX = 'p_{T}(l_{2}) (GeV)', texY = 'Number of Events / 20 GeV',
#       attribute = lambda event, sample:event.l2_pt,
#       binning=[15,0,300],
#       addOverFlowBin='upper',
#   ))
#
#    plots.append(Plot(
#       name = 'l2_eta',
#       texX = '#eta(l_{2})', texY = 'Number of Events',
#       attribute = lambda event, sample: event.l2_eta,
#       binning=[20,-3,3],
#   ))
#
#    plots.append(Plot(
#       name = 'l2_mvaTOP',
#       texX = 'MVA_{TOP}(l_{2})', texY = 'Number of Events',
#       attribute = lambda event, sample: event.l2_mvaTOP,
#       binning=[20,-1,1],
#   ))
#
#    plots.append(Plot(
#       name = 'l2_mvaTOPWP',
#       texX = 'MVA_{TOP}(l_{1}) WP', texY = 'Number of Events',
#       attribute = lambda event, sample: event.l2_mvaTOPWP,
#       binning=[5,0,5],
#   ))
#
#    plots.append(Plot(
#       name = 'l3_pt',
#       texX = 'p_{T}(l_{3}) (GeV)', texY = 'Number of Events / 20 GeV',
#       attribute = lambda event, sample:event.l3_pt,
#       binning=[15,0,300],
#       addOverFlowBin='upper',
#   ))
#
#    plots.append(Plot(
#       name = 'l3_eta',
#       texX = '#eta(l_{3})', texY = 'Number of Events',
#       attribute = lambda event, sample: event.l3_eta,
#       binning=[20,-3,3],
#   ))
#
#    plots.append(Plot(
#       name = 'l3_mvaTOP',
#       texX = 'MVA_{TOP}(l_{3})', texY = 'Number of Events',
#       attribute = lambda event, sample: event.l3_mvaTOP,
#       binning=[20,-1,1],
#   ))
#
#    plots.append(Plot(
#       name = 'l3_mvaTOPWP',
#       texX = 'MVA_{TOP}(l_{1}) WP', texY = 'Number of Events',
#       attribute = lambda event, sample: event.l3_mvaTOPWP,
#       binning=[5,0,5],
#   ))

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
