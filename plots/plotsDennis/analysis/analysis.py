#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT, os
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas() # do this to avoid version conflict in png.h with keras import ...
c1.Draw()
c1.Print('delete.png')
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
from tWZ.Tools.helpers          import getCollection

# Analysis
from Analysis.Tools.helpers              import deltaPhi, deltaR
from Analysis.Tools.puProfileCache       import *
from Analysis.Tools.puReweighting        import getReweightingFunction
from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons

import Analysis.Tools.syncer
import numpy as np

################################################################################
# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--noData',         action='store_true', default=False, help='also plot data?')
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--dataMCScaling',  action='store_true', help='Data MC scaling?', )
argParser.add_argument('--plot_directory', action='store', default='tWZ_v1')
argParser.add_argument('--era',            action='store', type=str, default="Run2016")
argParser.add_argument('--selection',      action='store', default='trilepT-minDLmass12-onZ1-njet4p-btag1')
args = argParser.parse_args()

################################################################################
# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"
if args.noData:                       args.plot_directory += "_noData"

logger.info( "Working in era %s", args.era)
if args.dataMCScaling:
    logger.info( "Data/MC scaling active")
else:
    logger.info( "Data/MC scaling not active")

################################################################################
# Define the MC samples
from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed import *

if args.era == "Run2016":
    mc = [Summer16.TWZ_NLO_DR, Summer16.TTZ, Summer16.TTX_rare, Summer16.TZQ, Summer16.WZ, Summer16.triBoson, Summer16.ZZ, Summer16.nonprompt_3l]
elif args.era == "Run2017":
    mc = [Fall17.TWZ_NLO_DR, Fall17.TTZ, Fall17.TTX_rare, Fall17.TZQ, Fall17.WZ, Fall17.triBoson, Fall17.ZZ, Fall17.nonprompt_3l]
elif args.era == "Run2018":
    mc = [Autumn18.TWZ_NLO_DR, Autumn18.TTZ, Autumn18.TTX_rare, Autumn18.TZQ, Autumn18.WZ, Autumn18.triBoson, Autumn18.ZZ, Autumn18.nonprompt_3l]
elif args.era == "RunII":
    mc = [TWZ_NLO_DR, TTZ, TTX_rare, TZQ, WZ, triBoson, ZZ, nonprompt_3l]

################################################################################
# Define the data sample
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
        sample.reduceFiles( to = 1 )
        #sample.reduceFiles( to=1)
        sample.scale /= sample.normalization

################################################################################
# Text on the plots
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)
tex.SetTextAlign(11) # align right

################################################################################
# Functions needed specifically for this analysis routine
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

def getBJetindex( event ):
    maxscore = 0.0
    index = -1
    for i in range(event.nJetGood):
        btagscore = event.JetGood_btagDeepB[i]
        if btagscore > maxscore:
            maxscore = btagscore
            index = i
    return index


def getWhad( event ):
    # get Whad from min dijet mass
    # exclude highest b tag from jet selection
    index_btag = getBJetindex(event)
    min_mass = 1000.
    Whad = ROOT.TLorentzVector()
    for i in range(event.nJetGood):
        if i == index_btag:
            continue
        jetindex1 = event.JetGood_index[i]
        jet1  = ROOT.TLorentzVector()
        jet1.SetPtEtaPhiM(event.Jet_pt[jetindex1], event.Jet_eta[jetindex1], event.Jet_phi[jetindex1], event.Jet_mass[jetindex1])

        for j in range(event.nJetGood):
            if j == index_btag or j == i:
                continue
            jetindex2 = event.JetGood_index[j]
            jet2  = ROOT.TLorentzVector()
            jet2.SetPtEtaPhiM(event.Jet_pt[jetindex2], event.Jet_eta[jetindex2], event.Jet_phi[jetindex2], event.Jet_mass[jetindex2])

            dijet = jet1+jet2
            if dijet.M() < min_mass:
                min_mass = dijet.M()
                Whad = dijet

    return Whad

def getWlep( event ):
    Wlep = ROOT.TLorentzVector()
    lepton  = ROOT.TLorentzVector()
    met     = ROOT.TLorentzVector()
    lepton.SetPtEtaPhiM(event.lep_pt[event.nonZ1_l1_index], event.lep_eta[event.nonZ1_l1_index], event.lep_phi[event.nonZ1_l1_index], 0)
    met.SetPtEtaPhiM(event.met_pt, 0, event.met_phi, 0)

    lepton_pT = ROOT.TVector3(lepton.Px(), lepton.Py(), 0)
    neutrino_pT = ROOT.TVector3(met.Px(), met.Py(), 0)

    mass_w = 80.399
    mu = mass_w * mass_w / 2 + lepton_pT * neutrino_pT
    A = - (lepton_pT * lepton_pT)
    B = mu * lepton.Pz()
    C = mu * mu - lepton.E() * lepton.E() * (neutrino_pT * neutrino_pT)
    discriminant = B * B - A * C
    if discriminant <= 0:
        # Take only real part of the solution for pz:
        neutrino = ROOT.TLorentzVector()
        neutrino.SetPxPyPzE(met.Px(),met.Py(),-B / A,0)
        neutrino.SetE(neutrino.P())
    else:
        discriminant = sqrt(discriminant)
        neutrino1 = ROOT.TLorentzVector()
        neutrino1.SetPxPyPzE(met.Px(),met.Py(),(-B - discriminant) / A,0)
        neutrino1.SetE(neutrino1.P())
        neutrino2 = ROOT.TLorentzVector()
        neutrino2.SetPxPyPzE(met.Px(),met.Py(),(-B + discriminant) / A,0)
        neutrino2.SetE(neutrino2.P())
        if neutrino1.E() > neutrino2.E():
            neutrino = neutrino1
        else:
            neutrino = neutrino2

    Wlep = lepton + neutrino
    return Wlep

def getTopHypos( event ):
    Wh = getWhad( event )
    Wl = getWlep( event )
    bjetindex = event.JetGood_index[getBJetindex(event)]
    bJet = ROOT.TLorentzVector()
    bJet.SetPtEtaPhiM(event.Jet_pt[bjetindex], event.Jet_eta[bjetindex], event.Jet_phi[bjetindex], event.Jet_mass[bjetindex])

    toph = Wh + bJet
    topl = Wl + bJet

    return toph, topl

def get2TopHypos( event ):
    topA1 = ROOT.TLorentzVector()
    topA2 = ROOT.TLorentzVector()
    topB1 = ROOT.TLorentzVector()
    topB2 = ROOT.TLorentzVector()
    # Get Ws
    Wh = getWhad( event )
    Wl = getWlep( event )
    # Get highest btag
    bjetindex = event.JetGood_index[getBJetindex(event)]
    bJet1 = ROOT.TLorentzVector()
    bJet1.SetPtEtaPhiM(event.Jet_pt[bjetindex], event.Jet_eta[bjetindex], event.Jet_phi[bjetindex], event.Jet_mass[bjetindex])

    # now search for second highest bjet
    secondmax = 0
    index = -1
    bJet2 = ROOT.TLorentzVector()
    for i in range(event.nJetGood):
        if i == bjetindex:
            continue
        btagscore = event.JetGood_btagDeepB[i]
        if btagscore > secondmax:
            secondmax = btagscore
            index = i
    if index != -1:
        bjet2index = event.JetGood_index[index]
        bJet2.SetPtEtaPhiM(event.Jet_pt[bjet2index], event.Jet_eta[bjet2index], event.Jet_phi[bjet2index], event.Jet_mass[bjet2index])
    else:
        return topA1, topA2

    # construct tops
    topA1 = Wh + bJet1
    topA2 = Wl + bJet2
    topB1 = Wh + bJet2
    topB2 = Wl + bJet1

    # get hypothesis closest to mtop
    mtop = 172.5
    diffA1 = abs(topA1.M()-mtop)
    diffA2 = abs(topA2.M()-mtop)
    diffA = diffA1 if diffA1<diffA2 else diffA2
    diffB1 = abs(topB1.M()-mtop)
    diffB2 = abs(topB2.M()-mtop)
    diffB = diffB1 if diffB1<diffB2 else diffB2

    if diffA < diffB:
        return topA1, topA2
    else:
        return topB1, topB2

################################################################################
# Define sequences
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
    event.m3l = (l[0] + l[1] + l[2]).M()

sequence.append( getM3l )

def getBJetPt( event, sample ):
    index = getBJetindex(event)
    if index != -1:
        event.BJetPt = event.JetGood_pt[index]
        event.BJetIndex = index
    else:
        event.BJetPt = float('nan')
        event.BJetIndex = index

sequence.append( getBJetPt )

def getBadTop( event, sample):
    top1, top2  = get2TopHypos(event)
    mtop = 172.5
    if abs(top1.M()-172.5) < abs(top2.M()-172.5):
        event.BadTopMass = top2.M()
    else:
        event.BadTopMass = top1.M()

sequence.append( getBadTop )

def getWhadMass( event, sample):
    event.WhadMass = getWhad( event ).M()

sequence.append( getWhadMass )

def getWlepMass( event, sample):
    event.WlepMass = getWlep( event ).M()

sequence.append( getWlepMass )

def getTopProperties( event, sample):
    tophad, toplep = getTopHypos( event )

    if abs(tophad.M()-172.5) < abs(toplep.M()-172.5):
        top = tophad
        event.TopHadMass = top.M()
        event.TopLepMass = float('nan')

    else:
        top = toplep
        event.TopLepMass = top.M()
        event.TopHadMass = float('nan')


    event.TopMass = top.M()
    event.TopPt = top.Pt()

sequence.append( getTopProperties )

def getMaxEtaJet( event, sample ):
    maxeta = 0
    index = -1
    for i in range(event.nJet):
        if abs(event.Jet_eta[i]) > maxeta and event.Jet_pt[i] > 30:
            maxeta = abs(event.Jet_eta[i])
            index = i
    if index != -1:
        event.maxeta = event.Jet_eta[index]
    else:
        event.maxeta = 0.0

sequence.append( getMaxEtaJet )

def getDRTopZ( event, sample ):
    tophad, toplep = getTopHypos( event )

    if abs(tophad.M()-172.5) < abs(toplep.M()-172.5):
        top = tophad
    else:
        top = toplep

    Z = ROOT.TLorentzVector()
    Z.SetPtEtaPhiM(event.Z1_pt, event.Z1_eta, event.Z1_phi, event.Z1_mass)
    event.DRTopZ = top.DeltaR(Z)


sequence.append( getDRTopZ )

def get2TopChi2(event,sample):
    tophad, toplep = getTopHypos( event )
    # Values from UHH2
    Mtlep_mean  = 174.
    Mtlep_sigma =  18.
    Mthad_mean  = 181.
    Mthad_sigma =  15.
    chi2 = pow((tophad.M()-Mthad_mean)/Mthad_sigma,2)+pow((toplep.M()-Mtlep_mean)/Mtlep_sigma,2)
    event.chi2 = chi2
sequence.append( get2TopChi2 )

def getDeltaMaxEta(event,sample):
    jet1 = ROOT.TLorentzVector()
    # First get jet with max eta
    maxeta = 0
    i_maxeta = -1
    for i in range(event.nJet):
        if abs(event.Jet_eta[i]) > maxeta and event.Jet_pt[i] > 30:
            maxeta = abs(event.Jet_eta[i])
            found_jet1 = True
            i_maxeta = i
    jet1.SetPtEtaPhiM(event.Jet_pt[i_maxeta], event.Jet_eta[i_maxeta], event.Jet_phi[i_maxeta], event.Jet_mass[i_maxeta])
    # get second jet which is not the btag and gives max m_{ij}
    i_bjet = getBJetindex(event)
    maxdijetmass = 0
    i_maxdijetmass = -1
    for i in range(event.nJetGood):
        if i==i_bjet: continue
        jetindex = event.JetGood_index[i]
        jetcandidate = ROOT.TLorentzVector()
        jetcandidate.SetPtEtaPhiM(event.Jet_pt[jetindex], event.Jet_eta[jetindex], event.Jet_phi[jetindex], event.Jet_mass[jetindex])
        dijet = jet1+jetcandidate
        if dijet.M() > maxdijetmass:
            maxdijetmass=dijet.M()
            i_maxdijetmass=jetindex

    if i_maxeta!=-1 and i_maxdijetmass!=-1:
        event.deltamaxeta = abs(event.Jet_eta[i_maxdijetmass]-maxeta)
    else:
        event.deltamaxeta = float('nan')

sequence.append( getDeltaMaxEta )

def getAngles(event, sample=None):
    event.nonZ1_l1_Z1_deltaPhi = deltaPhi(event.lep_phi[event.nonZ1_l1_index], event.Z1_phi)
    event.Z1_j1_deltaPhi       = deltaPhi(event.Z1_phi, event.JetGood_phi[0])
    event.nonZ1_l1_Z1_deltaEta = abs(event.lep_eta[event.nonZ1_l1_index] - event.Z1_eta)
    event.nonZ1_l1_Z1_deltaR   = deltaR({'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet0_Z1_deltaR       = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet0_nonZ1_l1_deltaR = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    event.jet1_Z1_deltaR       = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet1_nonZ1_l1_deltaR = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    event.jet2_Z1_deltaR       = deltaR({'eta':event.JetGood_eta[2], 'phi':event.JetGood_phi[2]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet2_nonZ1_l1_deltaR = deltaR({'eta':event.JetGood_eta[2], 'phi':event.JetGood_phi[2]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    i_bjet = getBJetindex(event)
    event.bJet_Z1_deltaR      = deltaR({'eta':event.JetGood_eta[i_bjet], 'phi':event.JetGood_phi[i_bjet]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.bJet_nonZ1l1_deltaR = deltaR({'eta':event.JetGood_eta[i_bjet], 'phi':event.JetGood_phi[i_bjet]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
sequence.append( getAngles )


def forwardJets( event, sample=None ):
    jetVars          = ['pt/F', 'eta/F', 'phi/F', 'btagDeepB/F', 'jetId/I', 'btagDeepFlavB/F', 'mass/F']
    lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','cutBased/I','miniPFRelIso_all/F','pfRelIso03_all/F','mvaFall17V2Iso_WP90/O', 'mvaTOP/F', 'sip3d/F','lostHits/I','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','eleIndex/I','muIndex/I']
    jetVarNames      = [x.split('/')[0] for x in jetVars]
    lepVarNames      = [x.split('/')[0] for x in lepVars]
    #jets einlesen (in case of MC also reat the index of the genjet)
    alljets   = getCollection( event, 'Jet', jetVarNames, 'nJet')
    alljets.sort( key = lambda j: -j['pt'] )
    leptons   = getCollection(event, "lep", lepVarNames, 'nlep')
    # clean against good leptons
    clean_jets,_ = cleanJetsAndLeptons( alljets, leptons )
    # filter pt, but not eta
    jets_no_eta         = filter(lambda j:j['pt']>30, clean_jets)
    if jets_no_eta:
        event.maxAbsEta_of_pt30jets = max( [ abs(j['eta']) for j in jets_no_eta ])
    else:
        event.maxAbsEta_of_pt30jets = -1
sequence.append( forwardJets )

def getTTbar( event, sample):
    top1, top2  = get2TopHypos(event)
    mtop = 172.5
    if abs(top1.M()-172.5) < abs(top2.M()-172.5):
        event.MTop1 = top1.M()
        event.MTop2 = top2.M()
    else:
        event.MTop1 = top2.M()
        event.MTop2 = top1.M()

    event.MTop_average = (top1.M()+top2.M())/2

sequence.append( getTTbar )

################################################################################
# Read variables

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",  "nJet/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I",
    "l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
    "l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
    "JetGood[pt/F,eta/F,phi/F,area/F,btagDeepB/F,index/I]",
    "Jet[pt/F,eta/F,phi/F,mass/F]",
    "lep[pt/F,eta/F,phi/F,pdgId/I,muIndex/I,eleIndex/I]",
    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I",
    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
]

read_variables_MC = [
    'reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F',
    "genZ1_pt/F", "genZ1_eta/F", "genZ1_phi/F",
]

################################################################################
# MVA
import tWZ.MVA.configs as configs
config = configs.tWZ_3l
config_topReco = configs.tWZ_3l_topReco
config_lstm = configs.tWZ_3l
read_variables += config.read_variables
read_variables += config_topReco.read_variables
read_variables += config_lstm.read_variables


# Add sequence that computes the MVA inputs
def make_mva_inputs( event, sample ):
    for mva_variable, func in config.mva_variables:
        setattr( event, mva_variable, func(event, sample) )
    for mva_variable, func in config_topReco.mva_variables:
        setattr( event, mva_variable, func(event, sample) )
    for mva_variable, func in config_lstm.mva_variables:
        setattr( event, mva_variable, func(event, sample) )
sequence.append( make_mva_inputs )

# load models
from keras.models import load_model

models = [
    ("tWZ_3l", False, load_model("/mnt/hephy/cms/dennis.schwarz/tWZ/models/tWZ_3l_ttz/tWZ_3l/multiclass_model.h5")),
    ("tWZ_3l_topReco", False, load_model("/mnt/hephy/cms/dennis.schwarz/tWZ/models/tWZ_3l_ttz_topReco/tWZ_3l_topReco/multiclass_model.h5")),
    ("tWZ_3l_lstm", True, load_model("/mnt/hephy/cms/dennis.schwarz/tWZ/models/tWZ_3l_ttz_LSTM_LSTM/tWZ_3l/multiclass_model.h5")),
]

def keras_predict( event, sample ):
    # get model inputs assuming lstm
    for name, has_lstm, model in models:
        if name == "tWZ_3l": cfg = config
        elif name == "tWZ_3l_topReco": cfg = config_topReco
        elif name == "tWZ_3l_lstm": cfg = config_lstm
        flat_variables, lstm_jets = cfg.predict_inputs( event, sample, jet_lstm = True)
        # print name, has_lstm, flat_variables, lstm_jets
        prediction = model.predict( flat_variables if not has_lstm else [flat_variables, lstm_jets] )
        for i_val, val in enumerate( prediction[0] ):
            setattr( event, name+'_'+cfg.training_samples[i_val].name, val)
sequence.append( keras_predict )

################################################################################
# define 3l selections
mu_string  = lepString('mu','VL')
ele_string = lepString('ele','VL')
def getLeptonSelection( mode ):
    if   mode=="mumumu": return "Sum$({mu_string})==3&&Sum$({ele_string})==0".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="mumue":  return "Sum$({mu_string})==2&&Sum$({ele_string})==1".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="muee":   return "Sum$({mu_string})==1&&Sum$({ele_string})==2".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="eee":    return "Sum$({mu_string})==0&&Sum$({ele_string})==3".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=='all':    return "Sum$({mu_string})+Sum$({ele_string})==3".format(mu_string=mu_string,ele_string=ele_string)

################################################################################
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


################################################################################
# Set up channels and values for plotting
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

    ################################################################################
    # Now define the plots

    plots = []

    plots.append(Plot(
      name = 'yield', texX = '', texY = 'Number of Events',
      attribute = lambda event, sample: 0.5 + i_mode,
      binning=[4, 0, 4],
    ))

    plots.append(Plot(
      name = 'nVtxs', texX = 'vertex multiplicity', texY = 'Number of Events',
      attribute = TreeVariable.fromString( "PV_npvsGood/I" ),
      binning=[50,0,50],
      addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'l1_pt',
        texX = 'p_{T}(l_{1}) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = lambda event, sample:event.l1_pt,
        binning=[15,0,300],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'l1_eta',
        texX = '#eta(l_{1})', texY = 'Number of Events',
        attribute = lambda event, sample: event.l1_eta,
        binning=[20,-3,3],
    ))

    plots.append(Plot(
        name = 'l1_mvaTOP',
        texX = 'MVA_{TOP}(l_{1})', texY = 'Number of Events',
        attribute = lambda event, sample: event.l1_mvaTOP,
        binning=[20,-1,1],
    ))

    plots.append(Plot(
        name = 'l1_mvaTOPWP',
        texX = 'MVA_{TOP}(l_{1}) WP', texY = 'Number of Events',
        attribute = lambda event, sample: event.l1_mvaTOPWP,
        binning=[5,0,5],
    ))

    plots.append(Plot(
        name = 'l2_pt',
        texX = 'p_{T}(l_{2}) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = lambda event, sample:event.l2_pt,
        binning=[15,0,300],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'l2_eta',
        texX = '#eta(l_{2})', texY = 'Number of Events',
        attribute = lambda event, sample: event.l2_eta,
        binning=[20,-3,3],
    ))

    plots.append(Plot(
        name = 'l2_mvaTOP',
        texX = 'MVA_{TOP}(l_{2})', texY = 'Number of Events',
        attribute = lambda event, sample: event.l2_mvaTOP,
        binning=[20,-1,1],
    ))

    plots.append(Plot(
        name = 'l2_mvaTOPWP',
        texX = 'MVA_{TOP}(l_{1}) WP', texY = 'Number of Events',
        attribute = lambda event, sample: event.l2_mvaTOPWP,
        binning=[5,0,5],
    ))

    plots.append(Plot(
        name = 'l3_pt',
        texX = 'p_{T}(l_{3}) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = lambda event, sample:event.l3_pt,
        binning=[15,0,300],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'l3_eta',
        texX = '#eta(l_{3})', texY = 'Number of Events',
        attribute = lambda event, sample: event.l3_eta,
        binning=[20,-3,3],
    ))

    plots.append(Plot(
        name = 'l3_mvaTOP',
        texX = 'MVA_{TOP}(l_{3})', texY = 'Number of Events',
        attribute = lambda event, sample: event.l3_mvaTOP,
        binning=[20,-1,1],
    ))

    plots.append(Plot(
        name = 'l3_mvaTOPWP',
        texX = 'MVA_{TOP}(l_{1}) WP', texY = 'Number of Events',
        attribute = lambda event, sample: event.l3_mvaTOPWP,
        binning=[5,0,5],
    ))

    plots.append(Plot(
        texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "met_pt/F" ),
        binning=[400/20,0,400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        texX = '#phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "met_phi/F" ),
        binning=[10,-pi,pi],
    ))

    plots.append(Plot(
        name = "Z1_pt",
        texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "Z1_pt/F" ),
        binning=[20,0,400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'Z1_pt_coarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events / 50 GeV',
        attribute = TreeVariable.fromString( "Z1_pt/F" ),
        binning=[16,0,800],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'Z1_pt_superCoarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
        attribute = TreeVariable.fromString( "Z1_pt/F" ),
        binning=[3,0,600],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'Z1_pt_coarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events / 50 GeV',
        attribute = TreeVariable.fromString( "Z1_pt/F" ),
        binning=[16,0,800],
    ))

    plots.append(Plot(
        name = 'Z1_pt_superCoarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
        attribute = TreeVariable.fromString( "Z1_pt/F" ),
        binning=[3,0,600],
    ))

    plots.append(Plot(
        name = "M3l",
        texX = 'M(3l) (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample:event.m3l,
        binning=[25,0,500],
    ))

    plots.append(Plot(
        name = "dPhiZJet",
        texX = '#Delta#phi(Z,j1)', texY = 'Number of Events',
        attribute = lambda event, sample: deltaPhi(event.Z1_phi, event.JetGood_phi[0]),
        binning=[20,0,pi],
    ))

    plots.append(Plot(
        name = "l1_Z1_pt",
        texX = 'p_{T}(l_{1,Z}) (GeV)', texY = 'Number of Events / 10 GeV',
        attribute = lambda event, sample:event.lep_pt[event.Z1_l1_index],
        binning=[30,0,300],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "l1_Z1_pt_coarse",
        texX = 'p_{T}(l_{1,Z}) (GeV)', texY = 'Number of Events / 40 GeV',
        attribute = lambda event, sample:event.lep_pt[event.Z1_l1_index],
        binning=[10,0,400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'l1_Z1_pt_ext', texX = 'p_{T}(l_{1,Z}) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = lambda event, sample:event.lep_pt[event.Z1_l1_index],
        binning=[20,40,440],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "l2_Z1_pt",
        texX = 'p_{T}(l_{2,Z}) (GeV)', texY = 'Number of Events / 10 GeV',
        attribute = lambda event, sample:event.lep_pt[event.Z1_l2_index],
        binning=[20,0,200],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
      texX = 'p_{T}(leading l) (GeV)', texY = 'Number of Events / 20 GeV',
      name = 'lep1_pt', attribute = lambda event, sample: event.lep_pt[0],
      binning=[400/20,0,400],
    ))

    plots.append(Plot(
      texX = 'p_{T}(subleading l) (GeV)', texY = 'Number of Events / 10 GeV',
      name = 'lep2_pt', attribute = lambda event, sample: event.lep_pt[1],
      binning=[200/10,0,200],
    ))

    plots.append(Plot(
      texX = 'p_{T}(trailing l) (GeV)', texY = 'Number of Events / 10 GeV',
      name = 'lep3_pt', attribute = lambda event, sample: event.lep_pt[2],
      binning=[150/10,0,150],
    ))

    plots.append(Plot(
        name = "l2_Z1_pt_coarse",
        texX = 'p_{T}(l_{2,Z}) (GeV)', texY = 'Number of Events / 10 GeV',
        attribute = lambda event, sample:event.lep_pt[event.Z1_l2_index],
        binning=[10,0,200],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'l2_Z1_pt_ext', texX = 'p_{T}(l_{2,Z}) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = lambda event, sample:event.lep_pt[event.Z1_l2_index],
        binning=[20,0,400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'lnonZ1_pt',
        texX = 'p_{T}(l_{1,extra}) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = lambda event, sample:event.lep_pt[event.nonZ1_l1_index],
        binning=[15,0,300],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'lnonZ1_pt_coarse',
        texX = 'p_{T}(l_{1,extra}) (GeV)', texY = 'Number of Events / 60 GeV',
        attribute = lambda event, sample:event.lep_pt[event.nonZ1_l1_index],
        binning=[3,0,180],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'lnonZ1_charge',
        texX = 'Charge(l_{1,extra})', texY = 'Number of Events',
        attribute = lambda event, sample:-event.lep_pdgId[event.nonZ1_l1_index]/abs(event.lep_pdgId[event.nonZ1_l1_index]),
        binning=[2,-1,1],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'lnonZ1_eta',
        texX = '#eta(l_{1,extra})', texY = 'Number of Events',
        attribute = lambda event, sample: event.lep_eta[event.nonZ1_l1_index],
        binning=[20,-3,3],
    ))

    plots.append(Plot(
        texX = 'M(ll) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "Z1_mass/F" ),
        binning=[10,81,101],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "Z1_mass_wide",
        texX = 'M(ll) (GeV)', texY = 'Number of Events / 2 GeV',
        attribute = TreeVariable.fromString( "Z1_mass/F" ),
        binning=[50,20,120],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "Z1_cosThetaStar", texX = 'cos#theta(l-)', texY = 'Number of Events / 0.2',
        attribute = lambda event, sample:event.Z1_cosThetaStar,
        binning=[10,-1,1],
    ))

    plots.append(Plot(
        name = "Z2_mass_wide",
        texX = 'M(ll) of 2nd OSDL pair', texY = 'Number of Events / 2 GeV',
        attribute = TreeVariable.fromString( "Z2_mass/F" ),
        binning=[60,0,120],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "minDLmass",
        texX = 'min mass of all DL pairs', texY = 'Number of Events / 2 GeV',
        attribute = TreeVariable.fromString( "minDLmass/F" ),
        binning=[60,0,120],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        texX = '#Delta#phi(Z_{1}(ll))', texY = 'Number of Events',
        attribute = TreeVariable.fromString( "Z1_lldPhi/F" ),
        binning=[10,0,pi],
    ))

    plots.append(Plot(
        texX = '#Delta R(Z_{1}(ll))', texY = 'Number of Events',
        attribute = TreeVariable.fromString( "Z1_lldR/F" ),
        binning=[10,0,6],
    ))

    plots.append(Plot(
      texX = 'p_{T}^{gen}(Z1) (GeV)', texY = 'Number of Events / 30 GeV',
      name = 'genZ1_pt', attribute = lambda event, sample: event.genZ1_pt if not sample.isData else float('nan'),
      binning=[600/30,0,600],
    ))

    plots.append(Plot(
      texX = 'p_{T}(highest b-tagged jet) (GeV)', texY = 'Number of Events / 30 GeV',
      name = 'btag_pt',
      attribute = lambda event, sample: event.BJetPt,
      binning=[600/30,0,600],
    ))

    plots.append(Plot(
      texX = 'Index of b-tagged jet (GeV)', texY = 'Number of Events',
      name = 'btag_index',
      attribute = lambda event, sample: event.BJetIndex,
      binning=[8,-1.5, 6.5],
    ))

    plots.append(Plot(
      texX = 'm_{had W} (GeV)', texY = 'Number of Events / 10 GeV',
      name = 'Whad_mass',
      attribute = lambda event, sample: event.WhadMass,
      binning=[300/10,0,300],
    ))

    plots.append(Plot(
      texX = 'm_{lep W} (GeV)', texY = 'Number of Events / 10 GeV',
      name = 'Wlep_mass',
      attribute = lambda event, sample: event.WlepMass,
      binning=[300/10,0,300],
    ))

    plots.append(Plot(
      texX = 'm_{second top} (GeV)', texY = 'Number of Events / 20 GeV',
      name = 'BadTop_mass',
      attribute = lambda event, sample: event.BadTopMass,
      binning=[500/20,0,500],
    ))

    plots.append(Plot(
      texX = 'm_{had top} (GeV)', texY = 'Number of Events / 20 GeV',
      name = 'TopHad_mass',
      attribute = lambda event, sample: event.TopHadMass,
      binning=[500/20,0,500],
    ))

    plots.append(Plot(
      texX = 'm_{lep top} (GeV)', texY = 'Number of Events / 20 GeV',
      name = 'TopLep_mass',
      attribute = lambda event, sample: event.TopLepMass,
      binning=[500/20,0,500],
    ))

    plots.append(Plot(
      texX = 'm_{top} (GeV)', texY = 'Number of Events / 20 GeV',
      name = 'Top_mass',
      attribute = lambda event, sample: event.TopMass,
      binning=[500/20,0,500],
    ))

    plots.append(Plot(
      texX = 'm_{top1} (GeV)', texY = 'Number of Events / 20 GeV',
      name = 'MTop1',
      attribute = lambda event, sample: event.MTop1,
      binning=[500/20,0,500],
    ))

    plots.append(Plot(
      texX = 'm_{top2} (GeV)', texY = 'Number of Events / 20 GeV',
      name = 'MTop2',
      attribute = lambda event, sample: event.MTop2,
      binning=[500/20,0,500],
    ))

    plots.append(Plot(
      texX = '(m_{top1}+m_{top2})/2 (GeV)', texY = 'Number of Events / 20 GeV',
      name = 'MTop_average',
      attribute = lambda event, sample: event.MTop_average,
      binning=[500/20,0,500],
    ))

    plots.append(Plot(
      texX = 'p_{T}(top) (GeV)', texY = 'Number of Events / 30 GeV',
      name = 'Top_pt',
      attribute = lambda event, sample: event.TopPt,
      binning=[600/30,0,600],
    ))

    plots.append(Plot(
      texX = 'N_{jets}', texY = 'Number of Events',
      attribute = TreeVariable.fromString( "nJetGood/I" ), #nJetSelected
      binning=[8,-0.5,7.5],
    ))

    plots.append(Plot(
      texX = 'N_{b-tag}', texY = 'Number of Events',
      attribute = TreeVariable.fromString( "nBTag/I" ), #nJetSelected
      binning=[4,-0.5,3.5],
    ))

    plots.append(Plot(
      texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events / 30 GeV',
      name = 'jet0_pt', attribute = lambda event, sample: event.JetGood_pt[0],
      binning=[600/30,0,600],
    ))

    plots.append(Plot(
      texX = 'p_{T}(subleading jet) (GeV)', texY = 'Number of Events / 30 GeV',
      name = 'jet1_pt', attribute = lambda event, sample: event.JetGood_pt[1],
      binning=[600/30,0,600],
    ))

    plots.append(Plot(
        name = 'maxeta',
        texX = 'max #eta(jet)', texY = 'Number of Events',
        attribute = lambda event, sample: event.maxeta,
        binning=[20,-5,5],
    ))

    plots.append(Plot(
        name = 'jet0_eta',
        texX = '#eta(leading jet)', texY = 'Number of Events',
        attribute = lambda event, sample: event.JetGood_eta[0],
        binning=[20,-3,3],
    ))

    plots.append(Plot(
        name = 'jet1_eta',
        texX = '#eta(subleading jet)', texY = 'Number of Events',
        attribute = lambda event, sample: event.JetGood_eta[1],
        binning=[20,-3,3],
    ))

    plots.append(Plot(
        name = "W_pt",
        texX = 'p_{T}(W) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = lambda event, sample:event.W_pt,
        binning=[20,0,400],
    ))

    plots.append(Plot(
        name = "DR_TopZ",
        texX = '#Delta R (top, Z)', texY = 'Number of Events',
        attribute = lambda event, sample:event.DRTopZ,
        binning=[20,0,8],
    ))

    plots.append(Plot(
        name = "chi2",
        texX = '#chi^{2}', texY = 'Number of Events',
        attribute = lambda event, sample:event.chi2,
        binning=[50,0,200],
    ))

    plots.append(Plot(
        name = "DeltaMaxEta",
        texX = '#Delta #eta_{j}', texY = 'Number of Events',
        attribute = lambda event, sample:event.deltamaxeta,
        binning=[20,0,7],
    ))


    # MVA plot
    for name, has_lstm, model in models:
        #print has_lstm, flat_variables, lstm_jets
        if name == "tWZ_3l": cfg = config
        elif name == "tWZ_3l_topReco": cfg = config_topReco
        elif name == "tWZ_3l_lstm": cfg = config_lstm
        for i_tr_s, tr_s in enumerate( cfg.training_samples ):
            disc_name = name+'_'+cfg.training_samples[i_tr_s].name
            plots.append(Plot(
                texX = disc_name, texY = 'Number of Events',
                name = "MVA_"+disc_name,
                attribute = lambda event, sample, disc_name=disc_name: getattr( event, disc_name ),
                binning=[50, 0, 1],
            ))
            plots.append(Plot(
                texX = disc_name, texY = 'Number of Events',
                name = "MVA_"+disc_name+"_ZOOM",
                attribute = lambda event, sample, disc_name=disc_name: getattr( event, disc_name ),
                binning=[10, 0, 1],
            ))




    plotting.fill(plots, read_variables = read_variables, sequence = sequence)

    ################################################################################
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


################################################################################
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

    if mode == "all":
        drawPlots(allPlots['mumumu'], mode, dataMCScale)
        # Write MVA score in root file
        outfile = ROOT.TFile('MVA_score_'+args.era+'.root', 'recreate')
        outfile.cd()
        for plot in plots:
            if "MVA_tWZ" in plot.name:
                model = "MVA_tWZ_3l"
                if "topReco" in plot.name: model = "MVA_tWZ_3l_topReco"
                elif "lstm" in plot.name: model = "MVA_tWZ_3l_lstm"

                for i, l in enumerate(plot.histos):
                    for j, h in enumerate(l):
                        histname = h.GetName()
                        cropped_name = ""
                        if model+"_TWZ_NLO_DR" in histname:
                            node_name = "TWZnode"
                            cropped_name = histname.replace(model+"_TWZ_NLO_DR", "")
                        elif model+"_TTZ" in histname:
                            cropped_name = histname.replace(model+"_TTZ", "")
                            node_name = "TTZnode"
                        else:
                            print "[ERROR] cannot identify node of MVA plot"
                            continue
                        if "TWZ_NLO_DR" in cropped_name: process = "tWZ"
                        elif "TTZ" in cropped_name: process = "ttZ"
                        elif "TTX_rare" in cropped_name: process = "ttX"
                        elif "TZQ" in cropped_name: process = "tZq"
                        elif "WZ" in cropped_name: process = "WZ"
                        elif "ZZ" in cropped_name: process = "ZZ"
                        elif "triBoson" in cropped_name: process = "triBoson"
                        elif "nonprompt" in cropped_name: process = "nonprompt"
                        elif "data" in cropped_name: process = "data"
                        h.Write(model+"__"+process+"__"+node_name)
        outfile.Close()


logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
