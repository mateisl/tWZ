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
from math                                import sqrt, cos, sin, pi, atan2, cosh, exp

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
from Analysis.Tools.WeightInfo           import WeightInfo

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
argParser.add_argument('--plot_directory', action='store', default='tWZ_EFT_v1')
argParser.add_argument('--era',            action='store', type=str, default="Run2018")
argParser.add_argument('--selection',      action='store', default='trilepT-minDLmass12-onZ1-njet4p-deepjet1')
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
import tWZ.samples.nanoTuples_Autumn18_nanoAODv6_private_SMEFTsim_fast_postProcessed as SMEFTsim_fast

if args.era == "Run2016":
    mc = [Summer16.TWZ_NLO_DR, Summer16.TTZ, Summer16.TTX_rare, Summer16.TZQ, Summer16.WZ, Summer16.triBoson, Summer16.ZZ, Summer16.nonprompt_3l]
elif args.era == "Run2017":
    mc = [Fall17.TWZ_NLO_DR, Fall17.TTZ, Fall17.TTX_rare, Fall17.TZQ, Fall17.WZ, Fall17.triBoson, Fall17.ZZ, Fall17.nonprompt_3l]
elif args.era == "Run2018":
    mc = [Autumn18.TWZ_NLO_DR, Autumn18.TTZ, Autumn18.TTX_rare, Autumn18.TZQ, Autumn18.WZ, Autumn18.triBoson, Autumn18.ZZ, Autumn18.nonprompt_3l]
    sample_eft = SMEFTsim_fast.ttZ01j
elif args.era == "RunII":
    mc = [TWZ_NLO_DR, TTZ, TTX_rare, TZQ, WZ, triBoson, ZZ, nonprompt_3l]

################################################################################
# EFT reweight

# WeightInfo
w = WeightInfo(sample_eft.reweight_pkl)
w.set_order(2)

# define which Wilson coefficients to plot
#cHq1Re11 cHq1Re22 cHq1Re33 cHq3Re11 cHq3Re22 cHq3Re33 cHuRe11 cHuRe22 cHuRe33 cHdRe11 cHdRe22 cHdRe33 cHudRe11 cHudRe22 cHudRe33

WCs = [
   # ('cHq3Re11', 1.0, ROOT.kCyan),
   # ('cHq3Re22', 1.0, ROOT.kMagenta),
   # ('cHq3Re33', 1.0, ROOT.kBlue),
    ('cHq1Re11', 2.0, ROOT.kRed),
    ('cHq1Re22', 2.0, ROOT.kGreen),
    ('cHq1Re33', 2.0, ROOT.kOrange),
    # ('cHuRe11',  2.0, ROOT.kCyan),
    # ('cHuRe22',  2.0, ROOT.kMagenta),
    # ('cHuRe33',  2.0, ROOT.kBlue),
    # ('cHdRe11',  2.0, ROOT.kViolet-9),
    # ('cHdRe22',  2.0, ROOT.kGray),
    # ('cHdRe33',  2.0, ROOT.kAzure+10),
]

params =  [ ]

for i_wc, (WC, WCval, color) in enumerate(WCs):
    params.append ({'legendText':'%s=%3.2f'%(WC, WCval), 'color':color,  'WC':{WC:WCval} })

for i_param, param in enumerate(params):
    param['sample']   = sample_eft
    sample_eft.weight = lambda event, sample: event.weight*lumi_year[event.year]/1000
    param['style']    = styles.lineStyle( param['color'] )

# Creating a list of weights
weight = []
# Add MC weights
weight_mc = []
for sample in mc:
    weight_ = lambda event, sample: 1. # Add event.weight and lumi weight to sample.weight later
    weight_mc.append(weight_)
weight.append(weight_mc)
# Add data weight
if not args.noData: weight.append([lambda event, sample: event.weight])
# Add EFT weight
weight_eft = []
for param in params:
    eft_weight = w.get_weight_func(**param['WC'])
    weight.append([eft_weight])

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

for param in params:
    param['sample'].scale = 1

if args.small:
    for sample in mc + [data_sample]:
        sample.normalization = 1.
        sample.reduceFiles( to = 1 )
        sample.scale /= sample.normalization
    for param in params:
        param['sample'].normalization = 1.
        param['sample'].reduceFiles( to = 1 )
        param['sample'].scale /= sample.normalization
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
        btagscore = event.JetGood_btagDeepFlavB[i]
        if btagscore > maxscore:
            maxscore = btagscore
            index = i
    return index

def getDeepJetsWP(disc,year):
    WP_L = {2016:0.0614, 2017:0.0521, 2018:0.0494}
    WP_M = {2016:0.3093, 2017:0.3033, 2018:0.2770}
    WP_T = {2016:0.7221, 2017:0.7489, 2018:0.7264}
    wp = 0
    if disc > WP_L[year]: wp = 1
    if disc > WP_M[year]: wp = 2
    if disc > WP_T[year]: wp = 3
    return wp




## top reco functions ##############################################################
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
    neutrinos = []
    if discriminant <= 0:
        # Take only real part of the solution for pz:
        neutrino = ROOT.TLorentzVector()
        neutrino.SetPxPyPzE(met.Px(),met.Py(),-B / A,0)
        neutrino.SetE(neutrino.P())
        neutrinos.append(neutrino)
    else:
        discriminant = sqrt(discriminant)
        neutrino1 = ROOT.TLorentzVector()
        neutrino1.SetPxPyPzE(met.Px(),met.Py(),(-B - discriminant) / A,0)
        neutrino1.SetE(neutrino1.P())
        neutrino2 = ROOT.TLorentzVector()
        neutrino2.SetPxPyPzE(met.Px(),met.Py(),(-B + discriminant) / A,0)
        neutrino2.SetE(neutrino2.P())
        if neutrino1.E() > neutrino2.E():
            neutrinos.append(neutrino1)
            neutrinos.append(neutrino2)
        else:
            neutrinos.append(neutrino2)
            neutrinos.append(neutrino1)

    Wleps = []
    for neu in neutrinos:
        Wlep = lepton + neu
        Wleps.append([Wlep, lepton, neu])
    return Wleps

def calculate_chi2(toplep, tophad, Whad, blep_disc, bhad_disc, mode):
    Mtlep_mean   = 171.
    Mtlep_sigma  =  16.
    Mthad_mean   = 171.
    Mthad_sigma  =  17.
    MWhad_mean   =  83.
    MWhad_sigma  =  11.
    Mtdiff_sigma = sqrt(pow(Mtlep_sigma,2)+pow(Mthad_sigma,2))
    b_disc_mean  = 1.0
    b_disc_sigma = 0.4
    toplepterm  = pow((toplep.M()-Mtlep_mean)/Mtlep_sigma,2)
    tophadterm  = pow((tophad.M()-Mthad_mean)/Mthad_sigma,2)
    Whadterm    = pow((  Whad.M()-MWhad_mean)/MWhad_sigma,2)
    topdiffterm = pow((toplep.M()-tophad.M())/Mtdiff_sigma,2)
    blepterm    = pow((blep_disc-b_disc_mean)/b_disc_sigma,2)
    bhadterm    = pow((bhad_disc-b_disc_mean)/b_disc_sigma,2)
    if mode=="topdiff":
        chi2 = topdiffterm+Whadterm
    elif mode=="btag":
        chi2 = topdiffterm+Whadterm+blepterm+bhadterm
    else:
        chi2 = toplepterm+tophadterm+Whadterm
    return chi2

def getTopHypos(event, Njetsmax):
    # Set maximum number of jets
    if event.nJetGood < Njetsmax:
        Njetsmax = event.nJetGood
    if Njetsmax < 4: return []
    # Find all possible solutions for the 4 missing jets
    # (blep,bhad,Wdecay1, Wdecay2)
    jet_permutations = []
    for i in range(Njetsmax):
        for j in range(Njetsmax):
            if j == i:
                continue
            for k in range(Njetsmax):
                if k==i or k==j:
                    continue
                for l in range(Njetsmax):
                    if l==i or l==j or l==k:
                        continue
                    jet_i = ROOT.TLorentzVector()
                    jet_j = ROOT.TLorentzVector()
                    jet_k = ROOT.TLorentzVector()
                    jet_l = ROOT.TLorentzVector()
                    jetidx_i = event.JetGood_index[i]
                    jetidx_j = event.JetGood_index[j]
                    jetidx_k = event.JetGood_index[k]
                    jetidx_l = event.JetGood_index[l]
                    jet_i.SetPtEtaPhiM(event.Jet_pt[jetidx_i], event.Jet_eta[jetidx_i], event.Jet_phi[jetidx_i], event.Jet_mass[jetidx_i])
                    jet_j.SetPtEtaPhiM(event.Jet_pt[jetidx_j], event.Jet_eta[jetidx_j], event.Jet_phi[jetidx_j], event.Jet_mass[jetidx_j])
                    jet_k.SetPtEtaPhiM(event.Jet_pt[jetidx_k], event.Jet_eta[jetidx_k], event.Jet_phi[jetidx_k], event.Jet_mass[jetidx_k])
                    jet_l.SetPtEtaPhiM(event.Jet_pt[jetidx_l], event.Jet_eta[jetidx_l], event.Jet_phi[jetidx_l], event.Jet_mass[jetidx_l])
                    jet_permutations.append([[jet_i, jet_j, jet_k, jet_l], [event.JetGood_btagDeepFlavB[i],event.JetGood_btagDeepFlavB[j]]])
    # Get Wlep from lepton + MET
    Wleps   = getWlep(event)
    # build hypotheses
    hypotheses = []
    for Wlep, lepton, neutrino in Wleps:
        for permutation, btag_disc in jet_permutations:
            hypo = {}
            hypo['toplep']       = Wlep + permutation[0]
            hypo['Wlep']         = Wlep
            hypo['lepton']       = lepton
            hypo['neutrino']     = neutrino
            hypo['blep']         = permutation[0]
            hypo['blep_disc']    = btag_disc[0]
            hypo['tophad']       = permutation[1] + permutation[2] + permutation[3]
            hypo['Whad']         = permutation[2] + permutation[3]
            hypo['bhad']         = permutation[1]
            hypo['bhad_disc']    = btag_disc[1]
            hypo['WhadDecay1']   = permutation[2]
            hypo['WhadDecay2']   = permutation[3]
            hypo['chi2']         = calculate_chi2(hypo['toplep'],hypo['tophad'],hypo['Whad'],hypo['blep_disc'],hypo['bhad_disc'],"normal")
            hypotheses.append(hypo)
    return hypotheses

################################################################################
# Define sequences
sequence       = []

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

def getWpt( event, sample=None):
    # get the lepton and met
    lepton  = ROOT.TLorentzVector()
    met     = ROOT.TLorentzVector()
    lepton.SetPtEtaPhiM(event.lep_pt[event.nonZ1_l1_index], event.lep_eta[event.nonZ1_l1_index], event.lep_phi[event.nonZ1_l1_index], 0)
    met.SetPtEtaPhiM(event.met_pt, 0, event.met_phi, 0)
    # get the W boson candidate
    W   = lepton + met
    event.W_pt = W.Pt()
sequence.append( getWpt )

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

## Top Reco sequence ###########################################################
def TopReco(event, sample):
    hypotheses=getTopHypos(event, 6)
    chi2min = 10000
    foundhypo = False
    for hypo in hypotheses:
        if hypo['chi2']<chi2min:
            chi2min = hypo['chi2']
            hypo_selected = hypo
            foundhypo = True
    if foundhypo:
        if hypo_selected['toplep'].M() > hypo_selected['tophad'].M():
            mtop_hi = hypo_selected['toplep'].M()
            mtop_lo = hypo_selected['tophad'].M()
        else:
            mtop_hi = hypo_selected['tophad'].M()
            mtop_lo = hypo_selected['toplep'].M()

    # Also get Z
    Z = ROOT.TLorentzVector()
    Z.SetPtEtaPhiM(event.Z1_pt, event.Z1_eta, event.Z1_phi, event.Z1_mass)

    LeptonNoZ = ROOT.TLorentzVector()
    LeptonNoZ.SetPtEtaPhiM(event.lep_pt[event.nonZ1_l1_index], event.lep_eta[event.nonZ1_l1_index], event.lep_phi[event.nonZ1_l1_index], 0)

    # observables for mtop1 and mtop2 closest to mtop
    event.mtop_average = (hypo_selected['toplep'].M()+hypo_selected['tophad'].M())/2 if foundhypo else -1
    event.mtoplep = hypo_selected['toplep'].M() if foundhypo else -1
    event.mtophad = hypo_selected['tophad'].M() if foundhypo else -1
    event.mtop_hi = mtop_hi if foundhypo else -1
    event.mtop_lo = mtop_lo if foundhypo else -1
    event.mtop_diff = (mtop_hi-mtop_lo)/(mtop_hi+mtop_lo) if foundhypo else -1
    event.pt_diff = abs(hypo_selected['toplep'].Pt()-hypo_selected['tophad'].Pt()) if foundhypo else -1
    event.mW_lep = hypo_selected['Wlep'].M() if foundhypo else -1
    event.mW_had = hypo_selected['Whad'].M() if foundhypo else -1
    event.ptW_lep = hypo_selected['Wlep'].Pt() if foundhypo else -1
    event.ptW_had = hypo_selected['Whad'].Pt() if foundhypo else -1
    event.chi2 = hypo_selected['chi2'] if foundhypo else -1
    event.pgof = exp(-0.5*hypo_selected['chi2']) if foundhypo else -1
    event.hypofound = 1 if foundhypo else 0
    event.dR_tops = hypo_selected['toplep'].DeltaR(hypo_selected['tophad']) if foundhypo else -1
    event.dR_Ws = hypo_selected['Wlep'].DeltaR(hypo_selected['Whad']) if foundhypo else -1
    event.dR_bottoms = hypo_selected['bhad'].DeltaR(hypo_selected['blep']) if foundhypo else -1
    event.dR_toplep_Z = hypo_selected['toplep'].DeltaR(Z) if foundhypo else -1
    event.dR_tophad_Z = hypo_selected['tophad'].DeltaR(Z) if foundhypo else -1
    event.dR_toplep_LepNoZ = hypo_selected['toplep'].DeltaR(LeptonNoZ) if foundhypo else -1
    event.dR_tophad_LepNoZ = hypo_selected['tophad'].DeltaR(LeptonNoZ) if foundhypo else -1
    event.dR_Wlep_LepNoZ = hypo_selected['Wlep'].DeltaR(LeptonNoZ) if foundhypo else -1
    event.dR_Whad_LepNoZ = hypo_selected['Whad'].DeltaR(LeptonNoZ) if foundhypo else -1
    event.dR_blep_LepNoZ = hypo_selected['blep'].DeltaR(LeptonNoZ) if foundhypo else -1
    event.dR_bhad_LepNoZ = hypo_selected['bhad'].DeltaR(LeptonNoZ) if foundhypo else -1
    event.blep_disc = hypo_selected['blep_disc'] if foundhypo else -1
    event.bhad_disc = hypo_selected['bhad_disc'] if foundhypo else -1
sequence.append(TopReco)

def getDeepJetTags(event, sample):
    nDeepJetLoose = 0
    nDeepJetMedium = 0
    nDeepJetTight = 0
    for i in range(event.nJetGood):
        idx_jet = event.JetGood_index[i]
        disc = event.Jet_btagDeepFlavB[idx_jet]
        wp = getDeepJetsWP(disc,event.year)
        if wp>=1: nDeepJetLoose+=1
        if wp>=2: nDeepJetMedium+=1
        if wp>=3: nDeepJetTight+=1
    event.nLoose  = nDeepJetLoose
    event.nMedium = nDeepJetMedium
    event.nTight  = nDeepJetTight

sequence.append(getDeepJetTags)

################################################################################
# Read variables

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",  "nJet/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I",
    "l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
    "l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
    "JetGood[pt/F,eta/F,phi/F,area/F,btagDeepB/F,btagDeepFlavB/F,index/I]",
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

read_variables_eft = [
    "np/I", VectorTreeVariable.fromString("p[C/F]",nMax=200)
]

################################################################################
# MVA
import tWZ.MVA.configs as configs
config = configs.tWZ_3l
config_topReco = configs.tWZ_3l_topReco
read_variables += config.read_variables
read_variables += config_topReco.read_variables


# Add sequence that computes the MVA inputs
def make_mva_inputs( event, sample ):
    for mva_variable, func in config.mva_variables:
        setattr( event, mva_variable, func(event, sample) )
    for mva_variable, func in config_topReco.mva_variables:
        setattr( event, mva_variable, func(event, sample) )
sequence.append( make_mva_inputs )

# load models
from keras.models import load_model

models = [
    ("tWZ_3l", False, load_model("/mnt/hephy/cms/dennis.schwarz/tWZ/models/tWZ_3l_ttz/tWZ_3l/multiclass_model.h5")),
    ("tWZ_3l_topReco", False, load_model("/mnt/hephy/cms/dennis.schwarz/tWZ/models/tWZ_3l_ttz_topReco/tWZ_3l_topReco/multiclass_model.h5")),
]

def keras_predict( event, sample ):
    # get model inputs assuming lstm
    for name, has_lstm, model in models:
        if name == "tWZ_3l": cfg = config
        elif name == "tWZ_3l_topReco": cfg = config_topReco
        flat_variables, lstm_jets = cfg.predict_inputs( event, sample, jet_lstm = True)
        # print name, has_lstm, flat_variables, lstm_jets
        prediction = model.predict( flat_variables if not has_lstm else [flat_variables, lstm_jets] )
        for i_val, val in enumerate( prediction[0] ):
            setattr( event, name+'_'+cfg.training_samples[i_val].name, val)
sequence.append( keras_predict )

################################################################################
# define 3l selections
mu_string  = lepString('mu','VL') + "&&lep_mediumId"
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

    for sample in mc: sample.style = styles.fillStyle(sample.color)

    for sample in mc:
      sample.read_variables = read_variables_MC
      sample.setSelectionString([getLeptonSelection(mode)])
      sample.weight = lambda event, sample: event.weight*lumi_year[event.year]/1000*event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger#*event.reweightLeptonSF

    for param in params:
      param['sample'].read_variables = read_variables_MC + read_variables_eft
      param['sample'].setSelectionString([getLeptonSelection(mode)])
      # param['sample'].weight = lambda event, sample: event.weight*lumi_year[event.year]/1000*event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger   #*event.reweightLeptonSF

    if not args.noData:
      stack = Stack(mc, data_sample, *[ [ param['sample'] ] for param in params ])
      noneftidxs = [0,1]
    else:
      stack = Stack(mc, *[ [ param['sample'] ] for param in params ])
      noneftidxs = [0]
    # Use some defaults
    Plot.setDefaults(stack = stack, weight = weight, selectionString = cutInterpreter.cutString(args.selection))

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
        name = 'Nbtag_DeepJet_L',
        texX = 'Number of b-tagged jets (loose)', texY = 'Number of Events',
        attribute = lambda event, sample: event.nLoose,
        binning=[11,-0.5, 10.5],
    ))

    plots.append(Plot(
        name = 'Nbtag_DeepJet_M',
        texX = 'Number of b-tagged jets (medium)', texY = 'Number of Events',
        attribute = lambda event, sample: event.nMedium,
        binning=[11,-0.5, 10.5],
    ))

    plots.append(Plot(
        name = 'Nbtag_DeepJet_T',
        texX = 'Number of b-tagged jets (tight)', texY = 'Number of Events',
        attribute = lambda event, sample: event.nTight,
        binning=[11,-0.5, 10.5],
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
        name = "DeltaMaxEta",
        texX = '#Delta #eta_{j}', texY = 'Number of Events',
        attribute = lambda event, sample:event.deltamaxeta,
        binning=[20,0,7],
    ))

    #### Top reco plots ########################################################
    plots.append(Plot(
        name = 'm_top_average',
        texX = '(m_{t1}+m_{t2})/2', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_average,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_lep',
        texX = 'm_{tlep}', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtoplep,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_had',
        texX = 'm_{thad}', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtophad,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_hi',
        texX = 'max(m_{t1},m_{t2})', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_hi,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_lo',
        texX = 'min(m_{t1},m_{t2})', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_lo,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_diff',
        texX = '|m_{t1}-m_{t2}|/|m_{t1}+m_{t2}|', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_diff,
        binning=[50, 0.0, 1.0],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'pt_top_diff',
        texX = '|p_{T}^{t1}-p_{T}^{t2}|', texY = 'Number of Events',
        attribute = lambda event, sample: event.pt_diff,
        binning=[40, 0.0, 400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'mW_lep',
        texX = 'm_{Wlep}', texY = 'Number of Events',
        attribute = lambda event, sample: event.mW_lep,
        binning=[40, 0.0, 200],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'mW_had',
        texX = 'm_{Whad}', texY = 'Number of Events',
        attribute = lambda event, sample: event.mW_had,
        binning=[40, 0.0, 200],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'ptW_lep',
        texX = 'p_{T}^{Wlep}', texY = 'Number of Events',
        attribute = lambda event, sample: event.ptW_lep,
        binning=[40, 0.0, 400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'ptW_had',
        texX = 'p_{T}^{Whad}', texY = 'Number of Events',
        attribute = lambda event, sample: event.ptW_had,
        binning=[40, 0.0, 400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'chi2',
        texX = '#chi^{2}', texY = 'Number of Events',
        attribute = lambda event, sample: event.chi2,
        binning=[40, 0.0, 80.],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'p_gof',
        texX = 'P_{gof}', texY = 'Number of Events',
        attribute = lambda event, sample: event.pgof,
        binning=[50, 0., 1.],
    ))

    plots.append(Plot(
        name = 'dR_tops',
        texX = '#Delta R(top_{lep},top_{had})', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_tops,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
    name = 'dR_Ws',
    texX = '#Delta R(W_{lep},W_{had})', texY = 'Number of Events',
    attribute = lambda event, sample: event.dR_Ws,
    binning=[35, 0.0, 7],
    addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'dR_bottoms',
        texX = '#Delta R(b_{lep},b_{had})', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_bottoms,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'dR_toplep_Z',
        texX = '#Delta R(top_{lep},Z)', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_toplep_Z,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'dR_tophad_Z',
        texX = '#Delta R(top_{had},Z)', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_tophad_Z,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'dR_toplep_LepNoZ',
        texX = '#Delta R(top_{lep},lep (not Z))', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_toplep_LepNoZ,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'dR_tophad_LepNoZ',
        texX = '#Delta R(top_{had},lep (no Z))', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_tophad_LepNoZ,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'dR_Wlep_LepNoZ',
        texX = '#Delta R(W_{lep},lep (not Z))', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_Wlep_LepNoZ,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'dR_Whad_LepNoZ',
        texX = '#Delta R(W_{had},lep (no Z))', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_Whad_LepNoZ,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'dR_blep_LepNoZ',
        texX = '#Delta R(b_{lep},lep (not Z))', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_blep_LepNoZ,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'dR_bhad_LepNoZ',
        texX = '#Delta R(b_{had},lep (no Z))', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_bhad_LepNoZ,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    # MVA plot
    for name, has_lstm, model in models:
        #print has_lstm, flat_variables, lstm_jets
        if name == "tWZ_3l": cfg = config
        elif name == "tWZ_3l_topReco": cfg = config_topReco
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


    for plot in plots:
        for idx, histo_list in enumerate(plot.histos):
            if idx in noneftidxs:
                continue
            # Get number of EFT sample (idx 0=MCstack, 1=data, 2=EFT0, 3=EFT1,...)
            # So, identify how many stacks are not EFT
            n_noneft = len(noneftidxs)
            histo_list[0].legendText = params[idx-n_noneft]['legendText']
            histo_list[0].style = params[idx-n_noneft]['style']

    # Now plot
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
        plot_dir = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, mode, args.selection)
        outfile = ROOT.TFile(plot_dir+'/MVA_score.root', 'recreate')
        outfile.cd()
        for plot in allPlots['mumumu']:
            if "MVA_tWZ" in plot.name and not "_ZOOM" in plot.name:
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
