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
from tWZ.Tools.helpers                   import getCollection, cosThetaStarNew, getTheta, gettheta, getphi

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
argParser.add_argument('--small',          action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',       action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--dataMCScaling',  action='store_true', help='Data MC scaling?', )
argParser.add_argument('--plot_directory', action='store', default='EFT_SYS_v1')
argParser.add_argument('--era',            action='store', type=str, default="Run2018")
argParser.add_argument('--selection',      action='store', default='trilepVL-minDLmass12-onZ1-onZ2-nLeptons4')
argParser.add_argument('--sys',            action='store', default='central')
argParser.add_argument('--nicePlots',      action='store_true', default=False)
argParser.add_argument('--twoD',           action='store_true', default=False)
args = argParser.parse_args()

################################################################################
# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"
if args.noData:                       args.plot_directory += "_noData"
if args.sys is not 'central':         args.plot_directory += "_%s" %(args.sys)


logger.info( "Working in era %s", args.era)
if args.dataMCScaling:
    logger.info( "Data/MC scaling active")
else:
    logger.info( "Data/MC scaling not active")

if args.nicePlots:
    logger.info( "Only draw the plots")
else:
    logger.info( "Only saving into root file")

if args.twoD:
    logger.info( "Create EFT points in 2D")
else:
    logger.info( "Create EFT points in 1D")
################################################################################
# Possible SYS variations
variations = [
    "Trigger_UP", "Trigger_DOWN",
    "LepID_UP", "LepID_DOWN",
    "BTag_b_UP", "BTag_b_DOWN",
    "BTag_l_UP", "BTag_l_DOWN",
    "PU_UP", "PU_DOWN",
    "JES_UP", "JES_DOWN"
]

jet_variations = {
    "JES_UP": "jesTotalUp",
    "JES_DOWN": "jesTotalDown",
}
################################################################################
# Check if we know the variation else don't use data!
if args.sys is not 'central' and args.sys not in variations:
    raise RuntimeError( "Variation %s not among the known: %s", args.sys, ",".join( variations ) )
else:
    args.noData = True


################################################################################
# Selection modifier
def jetSelectionModifier( sys, returntype = "func"):
    #Need to make sure all jet variations of the following observables are in the ntuple
    variiedJetObservables = ['nJetGood', 'nBTag', 'met_pt']
    if returntype == "func":
        def changeCut_( string ):
            for s in variiedJetObservables:
                string = string.replace(s, s+'_'+sys)
            return string
        return changeCut_
    elif returntype == "list":
        return [ v+'_'+sys for v in variiedJetObservables ]

def metSelectionModifier( sys, returntype = 'func'):
    #Need to make sure all MET variations of the following observables are in the ntuple
    variiedMetObservables = ['met_pt']
    if returntype == "func":
        def changeCut_( string ):
            for s in variiedMetObservables:
                string = string.replace(s, s+'_'+sys)
            return string
        return changeCut_
    elif returntype == "list":
        return [ v+'_'+sys for v in variiedMetObservables ]

################################################################################
# Add a selection selectionModifier

if args.sys in jet_variations.keys():
    selectionModifier = jetSelectionModifier(jet_variations[args.sys])
else:
    selectionModifier = None

################################################################################
# Define the MC samples
from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed import *
import tWZ.samples.nanoTuples_Autumn18_nanoAODv6_private_SMEFTsim_fast_postProcessed as SMEFTsim_fast

if args.era == "Run2016":
    mc = [Summer16.TWZ_NLO_DR, Summer16.TTZ, Summer16.TTX_rare, Summer16.TZQ, Summer16.WZ, Summer16.triBoson, Summer16.ZZ, Summer16.nonprompt_3l]
elif args.era == "Run2017":
    mc = [Fall17.TWZ_NLO_DR, Fall17.TTZ, Fall17.TTX_rare, Fall17.TZQ, Fall17.WZ, Fall17.triBoson, Fall17.ZZ, Fall17.nonprompt_3l]
elif args.era == "Run2018":
    mc = [Autumn18.TWZ_NLO_DR, Autumn18.TTX_rare, Autumn18.TZQ, Autumn18.triBoson, Autumn18.nonprompt_3l]
    samples_eft = [SMEFTsim_fast.ZZ, SMEFTsim_fast.WZ, SMEFTsim_fast.ttZ01j]
    if args.nicePlots:
        mc = [Autumn18.TWZ_NLO_DR, Autumn18.TTX_rare, Autumn18.TZQ, Autumn18.triBoson, Autumn18.nonprompt_3l, SMEFTsim_fast.ZZ, SMEFTsim_fast.WZ, SMEFTsim_fast.ttZ01j]
        signal_idx = [5, 6, 7] # INDEX OF signal sample
elif args.era == "RunII":
    mc = [TWZ_NLO_DR, TTZ, TTX_rare, TZQ, WZ, triBoson, ZZ, nonprompt_3l]

################################################################################
# EFT reweight

# WeightInfo
eftweights = []
for sample in samples_eft:
    w = WeightInfo(sample.reweight_pkl)
    w.set_order(2)
    eftweights.append(w)

# define which Wilson coefficients to plot
#cHq1Re11 cHq1Re22 cHq1Re33 cHq3Re11 cHq3Re22 cHq3Re33 cHuRe11 cHuRe22 cHuRe33 cHdRe11 cHdRe22 cHdRe33 cHudRe11 cHudRe22 cHudRe33

# WCs = [
#    # ('cHq3Re11', 1.0, ROOT.kCyan),
#    # ('cHq3Re22', 1.0, ROOT.kMagenta),
#    # ('cHq3Re33', 1.0, ROOT.kBlue),
#    ('cHq1Re11', -2.0, ROOT.kRed),
#     ('cHq1Re11', 2.0, ROOT.kRed),
#     ('cHq1Re22', -2.0, ROOT.kGreen+2),
#     ('cHq1Re22', 2.0, ROOT.kGreen+2),
#     ('cHq1Re33', -2.0, ROOT.kOrange-3),
#     ('cHq1Re33', 2.0, ROOT.kOrange-3),
#     # ('cHuRe11',  2.0, ROOT.kCyan),
#     # ('cHuRe22',  2.0, ROOT.kMagenta),
#     # ('cHuRe33',  2.0, ROOT.kBlue),
#     # ('cHdRe11',  2.0, ROOT.kViolet-9),
#     # ('cHdRe22',  2.0, ROOT.kGray),
#     # ('cHdRe33',  2.0, ROOT.kAzure+10),
# ]

minval = -10.0
maxval = 10.0
Npoints = 51

if args.nicePlots:
    Npoints = 0

WCs = []
WC_setup = [
    ('cHq1Re11', ROOT.kRed),
    ('cHq1Re22', ROOT.kGreen+2),
    ('cHq1Re33', ROOT.kOrange-3),
    ('cHq3Re11', ROOT.kCyan),
    ('cHq3Re22', ROOT.kMagenta),
    ('cHq3Re33', ROOT.kBlue),
]
for i_wc, (WCname, color) in enumerate(WC_setup):
    for i in range(Npoints):
        value = minval + ((maxval-minval)/(Npoints-1))*i
        WCs.append( (WCname, value, color) )

params =  []
for i_sample, sample in enumerate(samples_eft):
    for i_wc, (WC, WCval, color) in enumerate(WCs):
        params.append({'legendText':'%s=%3.2f'%(WC, WCval), 'color':color,  'WC':{WC:WCval} , 'sample': sample, 'i_sample': i_sample})

#### 2D scan
if args.twoD:
    minval = -4.0
    maxval = 4.0
    Npoints = 11
    params = []
    for i in range(Npoints):
        value1 = minval + ((maxval-minval)/(Npoints-1))*i
        for j in range(Npoints):
            value2 = minval + ((maxval-minval)/(Npoints-1))*j
            for i_sample, sample in enumerate(samples_eft):
                params.append({'legendText':'%i,%i'%(i, j), 'color':ROOT.kRed,  'WC':{'cHq1Re11':value1, 'cHq3Re11':value2} , 'sample': sample, 'i_sample': i_sample})

####

for i_param, param in enumerate(params):
    param['style']    = styles.lineStyle( param['color'] )

# Creating a list of weights
plotweights = []
weights_SM = []
# Add MC weights
weight_mc = []
for sample in mc:
    weight_ = lambda event, sample: 1. # Add event.weight and lumi weight to sample.weight later
    weight_mc.append(weight_)
plotweights.append(weight_mc)

# Add data weight
if not args.noData:
    plotweights.append([lambda event, sample: event.weight])
# Add EFT weight
for param in params:
    i_sample = param['i_sample']
    eft_weight = eftweights[i_sample].get_weight_func(**param['WC'])
    plotweights.append([eft_weight])

################################################################################
# Define the data sample
try:
  data_sample = eval(args.era)
except Exception as e:
  logger.error( "Didn't find %s", args.era )
  raise e

lumi_scale                 = data_sample.lumi/1000
data_sample.scale          = 1.

# Set up MC sample
for sample in mc:
    sample.scale           = 1

for param in params:
    param['sample'].scale = 1

if args.small:
    for sample in mc + [data_sample]:
        sample.normalization = 1.
        sample.reduceFiles( to = 1 )
        sample.scale /= sample.normalization
    for param in params[i_sample]:
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
            n_stacks=len(plot.histos)
            if isinstance( plot, Plot):
                plotting.draw(plot,
                  plot_directory = plot_directory_,
                  ratio = {'histos': [[i+1,0] for i in range(n_stacks-1)], 'yRange':(0.1,1.9)} if not args.noData else None,
                  logX = False, logY = log, sorting = True,
                  yRange = (0.03, "auto") if log else (0.001, "auto"),
                  scaling = {0:1} if args.dataMCScaling else {},
                  legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
                  drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ) + _drawObjects,
                  copyIndexPHP = True, extensions = ["png"],
                )

def getDeepJetsWP(disc,year):
    WP_L = {2016:0.0614, 2017:0.0521, 2018:0.0494}
    WP_M = {2016:0.3093, 2017:0.3033, 2018:0.2770}
    WP_T = {2016:0.7221, 2017:0.7489, 2018:0.7264}
    wp = 0
    if disc > WP_L[year]: wp = 1
    if disc > WP_M[year]: wp = 2
    if disc > WP_T[year]: wp = 3
    return wp

def getClosestBJetindex( event, object, minBtagValue ):
    minDR = 100
    closestjet = ROOT.TLorentzVector()
    for i in range(event.nJetGood):
        btagscore = event.JetGood_btagDeepFlavB[i]
        if btagscore > minBtagValue:
            jet = ROOT.TLorentzVector()
            jetidx = event.JetGood_index[i]
            jet.SetPtEtaPhiM(event.Jet_pt[jetidx], event.Jet_eta[jetidx], event.Jet_phi[jetidx], event.Jet_mass[jetidx])
            if object.DeltaR(jet) < minDR:
                minDR = object.DeltaR(jet)
                closestjet = jet
    return closestjet

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
################################################################################
# Define sequences
sequence       = []

def getMlb(sample, event):
    lepton = ROOT.TLorentzVector()
    if not 'nLeptons4' in args.selection:
        lepton.SetPtEtaPhiM(event.lep_pt[event.nonZ1_l1_index], event.lep_eta[event.nonZ1_l1_index], event.lep_phi[event.nonZ1_l1_index], 0)
    else:
        lepton.SetPtEtaPhiM(0,0,0,0)

    # Get closest tight b:
    minBtagValue = 0.7221
    if   event.year == 2017: minBtagValue = 0.7489
    elif event.year == 2018: minBtagValue = 0.7264
    closestjet = ROOT.TLorentzVector()
    closestjet = getClosestBJetindex( event, lepton, minBtagValue )
    combination = closestjet + lepton
    event.mlb = combination.M()
sequence.append( getMlb )

def getCosThetaStar(sample, event):
    lepton = ROOT.TLorentzVector()
    Z = ROOT.TLorentzVector()
    Z.SetPtEtaPhiM(event.Z1_pt, event.Z1_eta, event.Z1_phi, event.Z1_mass)
    idx_l1 = event.Z1_l1_index
    idx_l2 = event.Z1_l2_index
    if abs(event.lep_pdgId[idx_l1]) == 11:
        lepmass = 0.0005
    elif abs(event.lep_pdgId[idx_l1]) == 13:
        lepmass = 0.1
    elif abs(event.lep_pdgId[idx_l1]) == 15:
        lepmass = 1.777
    else:
        print "Z does not decay into leptons"
        lepmass = 0

    charge_l1 = event.lep_pdgId[idx_l1]/abs(event.lep_pdgId[idx_l1])
    charge_l2 = event.lep_pdgId[idx_l2]/abs(event.lep_pdgId[idx_l2])
    if charge_l1 < 0 and charge_l2 > 0:
        lepton.SetPtEtaPhiM(event.lep_pt[idx_l1], event.lep_eta[idx_l1], event.lep_phi[idx_l1], lepmass)


    event.cosThetaStar = cosThetaStarNew(lepton, Z)
    # boostvector = Z.BoostVector()
    # lepton_newSys = ROOT.TLorentzVector()
    # lepton_newSys = lepton
    # lepton_newSys.Boost(-boostvector)
    # event.cosThetaStar = lepton_newSys.CosTheta()
sequence.append( getCosThetaStar )

def getDiBosonAngles(sample, event):
    lep1,lep2,boson = ROOT.TLorentzVector(),ROOT.TLorentzVector(),ROOT.TLorentzVector()
    # Get lep1 and lep2 from Z boson
    idx_l1 = event.Z1_l1_index
    idx_l2 = event.Z1_l2_index
    if abs(event.lep_pdgId[idx_l1]) == 11:
        lepmass = 0.0005
    elif abs(event.lep_pdgId[idx_l1]) == 13:
        lepmass = 0.1
    elif abs(event.lep_pdgId[idx_l1]) == 15:
        lepmass = 1.777
    else:
        print "Z does not decay into leptons"
        lepmass = 0
    lep1.SetPtEtaPhiM(event.lep_pt[idx_l1], event.lep_eta[idx_l1], event.lep_phi[idx_l1], lepmass)
    lep2.SetPtEtaPhiM(event.lep_pt[idx_l2], event.lep_eta[idx_l2], event.lep_phi[idx_l2], lepmass)

    # For 2nd boson distinguish between cases
    if "nLeptons4" in args.selection:
        boson.SetPtEtaPhiM(event.Z2_pt, event.Z2_eta, event.Z2_phi, event.Z2_mass)
    elif "deepjet0" in args.selection:
        Wleps = getWlep(event)
        if len(Wleps) == 1:
            boson = Wleps[0][0]
        elif len(Wleps) == 2:
            boson = Wleps[0][0] if Wleps[0][0].Pt()>Wleps[1][0].Pt() else Wleps[1][0]
        else:
            print "[ERROR] getWlep should not return more than 2 options"
    else:
        boson.SetPtEtaPhiM(0,0,0,0)

    event.Theta = getTheta(lep1, lep2, boson)
    event.theta = gettheta(lep1, lep2, boson)
    event.phi = getphi(lep1, lep2, boson)
sequence.append( getDiBosonAngles )
################################################################################
# Read variables

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",  "nJet/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I",
    "l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
    "l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
    "JetGood[pt/F,eta/F,phi/F,area/F,btagDeepB/F,btagDeepFlavB/F,index/I]",
    "Jet[pt/F,eta/F,phi/F,mass/F,btagDeepFlavB/F]",
    "lep[pt/F,eta/F,phi/F,pdgId/I,muIndex/I,eleIndex/I]",
    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I",
    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I,mediumId/O,tightId/O,isPFcand/B,isTracker/B,isGlobal/B]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
]

if "nLeptons4" in args.selection:
    read_variables.append("Z2_phi/F", "Z2_pt/F", "Z2_mass/F", "Z2_eta/F")

read_variables_MC = [
    "weight/F", 'reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F',
    "genZ1_pt/F", "genZ1_eta/F", "genZ1_phi/F",
    "Muon[genPartFlav/I]"
]

read_variables_eft = [
    "np/I", VectorTreeVariable.fromString("p[C/F]",nMax=200)
]


################################################################################
# MVA

################################################################################
# define 3l selections
if "lepVeto" in args.selection:
    mu_string  = lepString('mu','VL')
else:
    mu_string  = lepString('mu','VL') + "&&lep_mediumId"
ele_string = lepString('ele','VL')
def getLeptonSelection( mode ):
    if   mode=="mumumu":   return "Sum$({mu_string})==3&&Sum$({ele_string})==0".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="mumue":    return "Sum$({mu_string})==2&&Sum$({ele_string})==1".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="muee":     return "Sum$({mu_string})==1&&Sum$({ele_string})==2".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="eee":      return "Sum$({mu_string})==0&&Sum$({ele_string})==3".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="mumumumu": return "Sum$({mu_string})==4&&Sum$({ele_string})==0".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="mumuee":   return "Sum$({mu_string})==2&&Sum$({ele_string})==2".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="eeee":     return "Sum$({mu_string})==0&&Sum$({ele_string})==4".format(mu_string=mu_string,ele_string=ele_string)

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
allPlots_SM= {}
allModes   = ['mumumu', 'mumue', 'muee', 'eee']
if 'nLeptons4' in args.selection:
    allModes = ['mumumumu','mumuee','eeee']

print "Working on channels:", allModes

for i_mode, mode in enumerate(allModes):
    yields[mode] = {}
    if not args.noData:
        data_sample.texName = "data"
        data_sample.setSelectionString([getLeptonSelection(mode)])
        data_sample.name           = "data"
        data_sample.style          = styles.errorStyle(ROOT.kBlack)
        lumi_scale                 = data_sample.lumi/1000

    for sample in mc: sample.style = styles.fillStyle(sample.color)

    ###### SYS #################################################################
    if args.sys in jet_variations:
        new_variables = ['%s/F'%v for v in jetSelectionModifier(jet_variations[args.sys],'list')]
        read_variables_MC += new_variables
        read_variables    += new_variables

    weightnames = ['weight', 'reweightBTag_SF', 'reweightPU', 'reweightL1Prefire' , 'reweightTrigger'] # 'reweightLeptonSF'
    sys_weights = {
        "BTag_b_UP"     : ('reweightBTag_SF','reweightBTag_SF_b_Up'),
        "BTag_b_DOWN"   : ('reweightBTag_SF','reweightBTag_SF_b_Down'),
        "BTag_l_UP"     : ('reweightBTag_SF','reweightBTag_SF_l_Up'),
        "BTag_l_DOWN"   : ('reweightBTag_SF','reweightBTag_SF_l_Down'),
        'Trigger_UP'    : ('reweightTrigger','reweightTriggerUp'),
        'Trigger_DOWN'  : ('reweightTrigger','reweightTriggerDown'),
        'LepID_UP'      : ('reweightLeptonSF','reweightLeptonSFUp'),
        'LepID_DOWN'    : ('reweightLeptonSF','reweightLeptonSFDown'),
        'PU_UP'         : ('reweightPU','reweightPUUp'),
        'PU_DOWN'       : ('reweightPU','reweightPUDown'),
    }

    if args.sys in sys_weights:
        oldname, newname = sys_weights[args.sys]
        for i, weight in enumerate(weightnames):
            if weight == oldname:
                weightnames[i] = newname
                read_variables_MC += ['%s/F'%(newname)]

    getters = map( operator.attrgetter, weightnames)
    def weight_function( event, sample):
        # Calculate weight, this becomes: w = event.weightnames[0]*event.weightnames[1]*...
        w = reduce(operator.mul, [g(event) for g in getters], 1)
        # Get Lumi weight and multiply to weight
        lumi_weight = lumi_year[event.year]/1000.
        w *= lumi_weight
        return w


    for sample in mc:
        sample.read_variables = read_variables_MC
        sample.setSelectionString([getLeptonSelection(mode)])
        sample.weight = weight_function

    for param in params:
        param['sample'].read_variables = read_variables_MC + read_variables_eft
        param['sample'].setSelectionString([getLeptonSelection(mode)])
        param['sample'].weight = weight_function

    if not args.noData:
        stack = Stack(mc, data_sample, *[ [ param['sample'] ] for param in params ])
        noneftidxs = [0,1]
        if args.nicePlots:
            stack = Stack(mc, data_sample)
    else:
        stack = Stack(mc, *[ [ param['sample'] ] for param in params ])
        noneftidxs = [0]
        if args.nicePlots:
            stack = Stack(mc)

    # Use some defaults
    selection_string = selectionModifier(cutInterpreter.cutString(args.selection)) if selectionModifier is not None else cutInterpreter.cutString(args.selection)
    Plot.setDefaults(stack = stack, weight = plotweights, selectionString = selection_string)

    ################################################################################
    # Now define the plots

    plots = []

    plots.append(Plot(
      name = 'yield', texX = '', texY = 'Number of Events',
      attribute = lambda event, sample: 0.5 + i_mode,
      binning=[4, 0, 4],
    ))

    plots.append(Plot(
        name = "Z1_pt",
        texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events / 160 GeV',
        attribute = TreeVariable.fromString( "Z1_pt/F" ),
        binning=[5, 0, 800],
    ))

    if args.nicePlots:
        plots.append(Plot(
            name = "M_lb",
            texX = 'm_{lb} (GeV)', texY = 'Number of Events / 20 GeV',
            attribute = lambda event, sample: event.mlb,
            binning=[20, 0, 400],
        ))

        plots.append(Plot(
            name = "CosThetaStar",
            texX = 'cos #theta^{*}', texY = 'Number of Events',
            attribute = lambda event, sample: event.cosThetaStar,
            binning=[20, -1, 1],
        ))

        plots.append(Plot(
            name = "CosThetaStar_old",
            texX = 'cos #theta^{*}', texY = 'Number of Events',
            attribute = lambda event, sample: event.Z1_cosThetaStar,
            binning=[20, -1, 1],
        ))

        plots.append(Plot(
            name = "theta",
            texX = '#theta', texY = 'Number of Events',
            attribute = lambda event, sample: event.theta,
            binning=[20, 0, pi],
        ))

        plots.append(Plot(
            name = "Theta",
            texX = '#Theta', texY = 'Number of Events',
            attribute = lambda event, sample: event.Theta,
            binning=[20, 0, pi],
        ))

        plots.append(Plot(
            name = "phi",
            texX = '#phi', texY = 'Number of Events',
            attribute = lambda event, sample: event.phi,
            binning=[20, -pi, pi],
        ))

    plotting.fill(plots, read_variables = read_variables, sequence = sequence)


    ################################################################################
    # Get normalization yields from yield histogram
    for plot in plots:
      if plot.name == "yield":
        for i, l in enumerate(plot.histos):
          for j, h in enumerate(l):
            yields[mode][plot.stack[i][j].name] = h.GetBinContent(h.FindBin(0.5+i_mode))
            if 'nLeptons4' in args.selection:
                h.GetXaxis().SetBinLabel(1, "#mu#mu#mu#mu")
                h.GetXaxis().SetBinLabel(2, "#mu#muee")
                h.GetXaxis().SetBinLabel(3, "eeee")
            else:
                h.GetXaxis().SetBinLabel(1, "#mu#mu#mu")
                h.GetXaxis().SetBinLabel(2, "#mu#mue")
                h.GetXaxis().SetBinLabel(3, "#muee")
                h.GetXaxis().SetBinLabel(4, "eee")

    if args.noData: yields[mode]["data"] = 0

    yields[mode]["MC"] = sum(yields[mode][s.name] for s in mc  + samples_eft)
    dataMCScale        = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')


    # if args.nicePlots:
    #     for plot in plots:
    #         for idx, histo_list in enumerate(plot.histos):
    #             if idx in noneftidxs:
    #                 continue
    #             # Get number of EFT sample (idx 0=MCstack, 1=data, 2=EFT0, 3=EFT1,...)
    #             # So, identify how many stacks are not EFT
    #             n_noneft = len(noneftidxs)
    #             histo_list[0].legendText = params[idx-n_noneft]['legendText']
    #             histo_list[0].style = params[idx-n_noneft]['style']
    #             # Now also add other SM samples (backgrounds)
    #             for idx_sm, histo_list_sm in enumerate(plot.histos):
    #                 # Only add SM MC
    #                 if idx_sm == 0:
    #                     for i in range(len(histo_list_sm)):
    #                         # Also don't add the signal again
    #                         if i in signal_idx:
    #                             continue
    #                         histo_list[0].Add(histo_list_sm[i])

    if args.nicePlots: drawPlots(plots, mode, dataMCScale)

    allPlots[mode] = plots


################################################################################
# Add all different channels
yields["all"] = {}
for y in yields[allModes[0]]:
    try:    yields["all"][y] = sum(yields[c][y] for c in allModes)
    except: yields["all"][y] = 0
dataMCScale = yields["all"]["data"]/yields["all"]["MC"] if yields["all"]["MC"] != 0 else float('nan')

allPlots["all"] = allPlots[allModes[0]]
for plot in allPlots['all']:
    for i_mode,mode in enumerate(allModes):
        if i_mode == 0:
            continue
        tmp = allPlots[mode]
        for plot2 in (p for p in tmp if p.name == plot.name):
            for i, j in enumerate(list(itertools.chain.from_iterable(plot.histos))):
                for k, l in enumerate(list(itertools.chain.from_iterable(plot2.histos))):
                    if i==k:
                        j.Add(l)


if args.nicePlots:
    drawPlots(allPlots['all'], "all", dataMCScale)

if not args.nicePlots:
    # Write Result Hist in root file
    plot_dir = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, "all", args.selection)
    if not os.path.exists(plot_dir):
        try:
            os.makedirs(plot_dir)
        except:
            print 'Could not crete', plot_dir
    outfilename = plot_dir+'/Results.root'
    if args.twoD: outfilename = plot_dir+'/Results_twoD.root'
    outfile = ROOT.TFile(outfilename, 'recreate')
    outfile.cd()
    for plot in allPlots['all']:
        if plot.name == "Z1_pt":
            for idx, histo_list in enumerate(plot.histos):
                for j, h in enumerate(histo_list):
                    histname = h.GetName()
                    if "TWZ_NLO_DR" in histname: process = "tWZ"
                    elif "TTZ" in histname: process = "ttZ"
                    elif "ttZ01j" in histname: process = "ttZ"
                    elif "TTX_rare" in histname: process = "ttX"
                    elif "TZQ" in histname: process = "tZq"
                    elif "WZ" in histname: process = "WZ"
                    elif "ZZ" in histname: process = "ZZ"
                    elif "triBoson" in histname: process = "triBoson"
                    elif "nonprompt" in histname: process = "nonprompt"
                    elif "data" in histname: process = "data"
                    # Also add a string for the eft signal samples
                    n_noneft = len(noneftidxs)
                    if idx not in noneftidxs: h.Write("Z1_pt__"+process+"__"+params[idx-n_noneft]['legendText'])
                    else: h.Write("Z1_pt__"+process)
    outfile.Close()



logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
