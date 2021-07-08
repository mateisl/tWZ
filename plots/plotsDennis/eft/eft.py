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
import numpy as np
from math                                import sqrt, cos, sin, pi, atan2, cosh, exp

# RootTools
from RootTools.core.standard             import *

# tWZ
from tWZ.Tools.user                      import plot_directory
from tWZ.Tools.cutInterpreter            import cutInterpreter
from tWZ.Tools.objectSelection           import cbEleIdFlagGetter, vidNestedWPBitMapNamingList
from tWZ.Tools.objectSelection           import lepString
# Analysis
from Analysis.Tools.helpers              import deltaPhi, deltaR
from Analysis.Tools.puProfileCache       import *
from Analysis.Tools.puReweighting        import getReweightingFunction
from Analysis.Tools.WeightInfo           import WeightInfo
import Analysis.Tools.syncer
import numpy as np

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--sample',         action='store', default='ttZ01j')
argParser.add_argument('--selection',      action='store', default='trilepL-minDLmass12-onZ1-njet3p-deepjet1p')
argParser.add_argument('--plot_directory', action='store', default='eft')

args = argParser.parse_args()

# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"

#tWZ_sample = TWZ if args.nominalSignal else yt_TWZ_filter
import tWZ.samples.nanoTuples_Autumn18_nanoAODv6_private_SMEFTsim_fast_postProcessed as SMEFTsim_fast

sample = getattr( SMEFTsim_fast, args.sample)

# WeightInfo
w = WeightInfo(sample.reweight_pkl)
w.set_order(2)

# define which Wilson coefficients to plot
#cHq1Re11 cHq1Re22 cHq1Re33 cHq3Re11 cHq3Re22 cHq3Re33 cHuRe11 cHuRe22 cHuRe33 cHdRe11 cHdRe22 cHdRe33 cHudRe11 cHudRe22 cHudRe33

WCs = [
#    ('cHq3Re11', 1.0, ROOT.kCyan),
#    ('cHq3Re22', 1.0, ROOT.kMagenta),
#    ('cHq3Re33', 1.0, ROOT.kBlue),
    ('cHq1Re11', 2.0, ROOT.kRed),
    ('cHq1Re22', 2.0, ROOT.kGreen),
    ('cHq1Re33', 2.0, ROOT.kOrange),
    ('cHuRe11',  2.0, ROOT.kCyan),
    ('cHuRe22',  2.0, ROOT.kMagenta),
    ('cHuRe33',  2.0, ROOT.kBlue),
    ('cHdRe11',  2.0, ROOT.kViolet-9),
    ('cHdRe22',  2.0, ROOT.kGray),
    ('cHdRe33',  2.0, ROOT.kAzure+10),
]

params =  [ ]

for i_wc, (WC, WCval, color) in enumerate(WCs):
    params.append ({'legendText':'%s=%3.2f'%(WC, WCval), 'color':color,  'WC':{WC:WCval} })
params +=  [ {'legendText':'SM',  'color':ROOT.kBlack, 'WC':{}} ]

lumi_scale = 137
for i_param, param in enumerate(params):
    param['sample'] = sample
    sample.weight   = lambda event, sample: event.weight*lumi_scale
    param['style']  = styles.lineStyle( param['color'] )

stack = Stack(*[ [ param['sample'] ] for param in params ] )
weight= [ [ w.get_weight_func(**param['WC']) ] for param in params ]

if args.small:
    sample.reduceFiles( to = 1 )

# Text on the plots
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)
tex.SetTextAlign(11) # align right

# load BIT
import sys
sys.path.append('/users/dennis.schwarz/CMSSW_10_6_0/src/BIT/')
from  BoostedInformationTree import BoostedInformationTree
bit = BoostedInformationTree.load('/users/dennis.schwarz/CMSSW_10_6_0/src/tWZ/BIT/python/BIT_ttZ.pkl')


def charge(pdgId):
    return -pdgId/abs(pdgId)

def drawObjects( plotData, lumi_scale ):
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if plotData else 'CMS Simulation'),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)'% ( lumi_scale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines]

def drawPlots(plots):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, 'eft',  args.plot_directory, sample.name, ("log" if log else "lin"), args.selection)
    for plot in plots:
      if not max(l.GetMaximum() for l in sum(plot.histos,[])): continue # Empty plot

      _drawObjects = []

      if isinstance( plot, Plot):
          plotting.draw(plot,
            plot_directory = plot_directory_,
            #ratio = {'yRange':(0.1,1.9)} if not args.noData else None,
            logX = False, logY = log, sorting = False,
            yRange = (0.03, "auto") if log else (0.001, "auto"),
            #scaling = {0:1} if args.dataMCScaling else {},
            legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
            drawObjects = drawObjects( False, lumi_scale ) + _drawObjects,
            copyIndexPHP = True, extensions = ["png", "pdf"],
          )

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

## BIT predict #################################################################
def BITpredict(event,sample):
    feature_list = []
    feature_list.append(event.Z1_pt)
    feature_list.append(event.Z1_phi)
    feature_list.append(event.Z1_eta)
    feature_list.append(event.Z1_mass)
    feature_list.append(event.lep_pt[event.Z1_l1_index])
    feature_list.append(event.lep_phi[event.Z1_l1_index])
    feature_list.append(event.lep_eta[event.Z1_l1_index])
    feature_list.append(event.lep_pt[event.Z1_l2_index])
    feature_list.append(event.lep_phi[event.Z1_l2_index])
    feature_list.append(event.lep_eta[event.Z1_l2_index])
    event.BIT = bit.predict(np.array(feature_list))
    # print event.BIT
sequence.append(BITpredict)


read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I",
    "l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
    "JetGood[pt/F,eta/F,phi/F,area/F,btagDeepB/F,btagDeepFlavB/F,index/I]",
    "Jet[pt/F,eta/F,phi/F,mass/F]",
    "lep[pt/F,eta/F,phi/F,pdgId/I,muIndex/I,eleIndex/I]",
    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I",
    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
    "np/I", VectorTreeVariable.fromString("p[C/F]",nMax=200)
]

read_variables_MC = ['reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F']

# define 3l selections
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


# Use some defaults
Plot.setDefaults(stack = stack, weight = weight, selectionString = "("+cutInterpreter.cutString(args.selection)+")&&("+getLeptonSelection('all')+")")

plots = []

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
    texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
    attribute = TreeVariable.fromString( "met_pt/F" ),
    binning=[400/20,0,400],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = "Z1_eta",
    texX = '#eta(Z_{1})', texY = 'Number of Events / 20 GeV',
    attribute = TreeVariable.fromString( "Z1_eta/F" ),
    binning=[20,-3,3],
    addOverFlowBin='upper',
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
    name = 'BIT_1', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT,
    binning=[50,-1,1],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'BIT_2', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT,
    binning=[50,-0.5,0.5],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'BIT_3', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT,
    binning=[50,-0.1,0.1],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'BIT_4', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT,
    binning=[50,-0.01,0.01],
    addOverFlowBin='upper',
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

plotting.fill(plots, read_variables = read_variables, sequence = sequence)

for plot in plots:
    for i_h, hl in enumerate(plot.histos):
        # dress up
        hl[0].legendText = params[i_h]['legendText']
        hl[0].style = params[i_h]['style']

drawPlots(plots)


logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
