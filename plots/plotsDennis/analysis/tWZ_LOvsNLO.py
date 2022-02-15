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
argParser.add_argument('--small',          action='store_true', help='Run only on a small subset of the data?', )
argParser.add_argument('--dataMCScaling',  action='store_true', help='Data MC scaling?', )
argParser.add_argument('--plot_directory', action='store', default='tWZ_LOvsNLO_v1')
argParser.add_argument('--era',            action='store', type=str, default="Run2018")
argParser.add_argument('--selection',      action='store', default='trilepT-minDLmass12-onZ1')
argParser.add_argument('--nicePlots',      action='store_true', default=False)
args = argParser.parse_args()

################################################################################
# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)


################################################################################
# Some info messages
if args.small:                        args.plot_directory += "_small"

logger.info( "Working in era %s", args.era)
if args.dataMCScaling:
    logger.info( "Data/MC scaling active")
else:
    logger.info( "Data/MC scaling not active")

if args.nicePlots:
    logger.info( "Only draw the plots")
else:
    logger.info( "Only saving into root file")


################################################################################
# Define the MC samples
from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed import *
import tWZ.samples.nanoTuples_Autumn18_nanoAODv6_private_SMEFTsim_fast_postProcessed as SMEFTsim_fast

mc1 = [SMEFTsim_fast.tWZToLL01j_lepWFilter]
mc2 = [Autumn18.TWZ_NLO_DR]
mc3 = [Autumn18.TWZ_NLO_DS]

################################################################################
# Creating a list of weights
plotweights = []
# Add MC weights
weight_mc1 = []
for sample in mc1:
    weight_ = lambda event, sample: 1. # Add event.weight and lumi weight to sample.weight later
    weight_mc1.append(weight_)
plotweights.append(weight_mc1)

weight_mc2 = []
for sample in mc2:
    weight_ = lambda event, sample: 1. # Add event.weight and lumi weight to sample.weight later
    weight_mc2.append(weight_)
plotweights.append(weight_mc2)

weight_mc3 = []
for sample in mc3:
    weight_ = lambda event, sample: 1. # Add event.weight and lumi weight to sample.weight later
    weight_mc3.append(weight_)
plotweights.append(weight_mc3)

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
for sample in mc1+mc2+mc3:
    sample.scale           = 1

if args.small:
    for sample in mc1 + mc2 + mc3 + [data_sample]:
        sample.normalization = 1.
        sample.reduceFiles( to = 1 )
        sample.scale /= sample.normalization
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
    if "mt2ll100" in args.selection: lines += [(0.55, 0.5, 'M_{T2}(ll) > 100 GeV')] # Manually put the mt2ll > 100 GeV label
    return [tex.DrawLatex(*l) for l in lines]

def drawPlots(plots, mode, dataMCScale):
    for log in [False, True]:
        plot_directory_ = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, mode + ("_log" if log else ""), args.selection)
        for plot in plots:
            if not max(l.GetMaximum() for l in sum(plot.histos,[])): continue # Empty plot
            _drawObjects = []
            n_stacks=len(plot.histos)
            if isinstance( plot, Plot):
                plotting.draw(plot,
                  plot_directory = plot_directory_,
                  ratio = None,
                  logX = False, logY = log, sorting = True,
                  yRange = (0.03, "auto") if log else (0.001, "auto"),
                  scaling = {0:1} if args.dataMCScaling else {},
                  legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
                  drawObjects = drawObjects( False, dataMCScale , lumi_scale ) + _drawObjects,
                  copyIndexPHP = True, extensions = ["png", "pdf", "root"],
                )

def getFinalstates(event):
    W,Z,top,b = -1, -1, -1, -1
    for i in range(event.nGenPart):
        if i > 10: break
        if event.GenPart_genPartIdxMother[i] in [0,1]:
            if abs(event.GenPart_pdgId[i]) == 6:
                top = i
            elif abs(event.GenPart_pdgId[i]) == 23:
                Z = i
            elif abs(event.GenPart_pdgId[i]) == 24:
                W = i
    return W,Z,top,b

def getbfromtop(event):
    for i in range(event.nGenPart):
        i_mother = event.GenPart_genPartIdxMother[i]
        if i_mother > 0:
            if abs(event.GenPart_pdgId[i]) == 5 and abs(event.GenPart_pdgId[i_mother]) == 6:
                return i
    return -1
################################################################################
# Define sequences
sequence       = []

def getBquarks(event,sample):
    W,Z,b,bTop = ROOT.TLorentzVector(),ROOT.TLorentzVector(),ROOT.TLorentzVector(),ROOT.TLorentzVector()
    i_W, i_Z, i_top, i_b = getFinalstates(event)
    foundW, foundZ, foundb, foundbtop = False, False, False, False
    if i_b > 0:
        b.SetPtEtaPhiM(event.GenPart_pt[i_b], event.GenPart_eta[i_b], event.GenPart_phi[i_b], event.GenPart_mass[i_b])
        foundb = True
    if i_top > 0:
        i_bTop = getbfromtop(event)
        if i_bTop > 0:
            bTop.SetPtEtaPhiM(event.GenPart_pt[i_bTop], event.GenPart_eta[i_bTop], event.GenPart_phi[i_bTop], event.GenPart_mass[i_bTop])
            foundbtop = True
    if i_W > 0:
        W.SetPtEtaPhiM(event.GenPart_pt[i_W], event.GenPart_eta[i_W], event.GenPart_phi[i_W], event.GenPart_mass[i_W])
        foundW = True
    if i_Z > 0:
        Z.SetPtEtaPhiM(event.GenPart_pt[i_Z], event.GenPart_eta[i_Z], event.GenPart_phi[i_Z], event.GenPart_mass[i_Z])
        foundZ = True
        
    event.b_pt = b.Pt() if foundb else float('nan')
    event.b_eta = b.Eta()  if foundb else float('nan')
    event.bTop_pt = bTop.Pt()  if foundbtop else float('nan')
    event.bTop_eta = bTop.Eta() if foundbtop else float('nan')
    event.b_diff_pt = bTop.Pt() - b.Pt() if (foundb and foundbtop) else float('nan')
    event.W_pt = W.Pt() if foundW else float('nan')
    event.W_eta = W.Eta() if foundW else float('nan')
    event.Z_eta = Z.Eta() if foundZ else float('nan')            
    event.Z_pt = Z.Pt() if foundZ else float('nan')

    pass_b_veto = True
    if foundb:
        if event.b_pt > 30 or abs(event.b_eta)<2.5:
            pass_b_veto = False
            
    pass_b_tag = False
    if foundbtop:
        if event.bTop_pt > 30 and abs(event.bTop_eta)<2.5:
            pass_b_tag = True

    event.pass_b_veto = pass_b_veto
    event.pass_b_tag = pass_b_tag
sequence.append(getBquarks)


def getMWZ(event,sample):
    mWZ = -1
    i_W, i_Z, i_top, i_b = getFinalstates(event)
    W,Z = ROOT.TLorentzVector(),ROOT.TLorentzVector()
    if i_W > 0 and i_Z > 0:
        W.SetPtEtaPhiM(event.GenPart_pt[i_W], event.GenPart_eta[i_W], event.GenPart_phi[i_W], event.GenPart_mass[i_W])
        Z.SetPtEtaPhiM(event.GenPart_pt[i_Z], event.GenPart_eta[i_Z], event.GenPart_phi[i_Z], event.GenPart_mass[i_Z])
        WZ = W+Z
        mWZ = WZ.M()
    event.mWZ = mWZ
sequence.append(getMWZ)

def countStatus(event,sample):
    nIn = 0
    nInter = 0
    nOut = 0
    for i in range(event.nGenPart):
        if   event.GenPart_status[i] == 21: nIn += 1
        elif event.GenPart_status[i] == 22: nInter += 1
        elif event.GenPart_status[i] == 23: nOut += 1
    event.nIn = nIn if nIn<10 else 10
    event.nInter = nInter if nInter<10 else 10
    event.nOut = nOut if nOut<10 else 10
sequence.append(countStatus)
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
    read_variables.append("Z2_phi/F")
    read_variables.append("Z2_pt/F")
    read_variables.append("Z2_eta/F")
    read_variables.append("Z2_mass/F")

read_variables_MC = [
    "weight/F", 'reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F',
    "genZ1_pt/F", "genZ1_eta/F", "genZ1_phi/F",
    "Muon[genPartFlav/I]",
    VectorTreeVariable.fromString( "GenPart[pt/F,mass/F,phi/F,eta/F,pdgId/I,genPartIdxMother/I,status/I,statusFlags/I]", nMax=1000),
    'nGenPart/I',
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


    for sample in mc1: sample.style = styles.fillStyle(sample.color)
    for sample in mc2: sample.style = styles.lineStyle(ROOT.kAzure+7)
    for sample in mc3: sample.style = styles.lineStyle(ROOT.kBlue)

    ###### Weights #################################################################
    weightnames = ['weight', 'reweightBTag_SF', 'reweightPU', 'reweightL1Prefire' , 'reweightTrigger'] # 'reweightLeptonSF'

    getters = map( operator.attrgetter, weightnames)
    def weight_function( event, sample):
        # Calculate weight, this becomes: w = event.weightnames[0]*event.weightnames[1]*...
        w = reduce(operator.mul, [g(event) for g in getters], 1)
        # Get Lumi weight and multiply to weight
        lumi_weight = lumi_year[event.year]/1000.
        w *= lumi_weight
        return w


    for sample in mc1:
        sample.read_variables = read_variables_MC
        sample.setSelectionString([getLeptonSelection(mode)])
        sample.weight = weight_function
    for sample in mc2:
        sample.read_variables = read_variables_MC
        sample.setSelectionString([getLeptonSelection(mode)])
        sample.weight = weight_function
    for sample in mc3:
        sample.read_variables = read_variables_MC
        sample.setSelectionString([getLeptonSelection(mode)])
        sample.weight = weight_function

    stack = Stack(mc1, mc2, mc3)

    # Use some defaults
    selection_string = cutInterpreter.cutString(args.selection)
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
      name = 'nIn', texX = '', texY = 'Number of incomming particles',
      attribute = lambda event, sample: event.nIn,
      binning=[11, -0.5, 10.5],
    ))

    plots.append(Plot(
      name = 'nInter', texX = '', texY = 'Number of intermediate particles',
      attribute = lambda event, sample: event.nInter,
      binning=[11, -0.5, 10.5],
    ))
    
    plots.append(Plot(
      name = 'nOut', texX = '', texY = 'Number of outgoing particles',
      attribute = lambda event, sample: event.nOut,
      binning=[11, -0.5, 10.5],
    ))
    
    plots.append(Plot(
        name = "gen_Z_eta",
        texX = 'Z #eta', texY = 'Number of Events',
        attribute = lambda event, sample: event.Z_eta,
        binning=[50, -5, 5],
    ))
    
    plots.append(Plot(
        name = "gen_W_eta",
        texX = 'W #eta', texY = 'Number of Events',
        attribute = lambda event, sample: event.W_eta,
        binning=[50, -5, 5],
    ))    
        
    plots.append(Plot(
        name = "gen_W_pt",
        texX = 'p_{T}(W) (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.W_pt,
        binning=[30, 0, 700],
    ))
            
    plots.append(Plot(
        name = "Z1_pt",
        texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
        attribute = TreeVariable.fromString( "Z1_pt/F" ),
        binning=[30, 0, 700],
    ))

    plots.append(Plot(
        name = "l1_pt",
        texX = 'p_{T}(l_{1}) (GeV)', texY = 'Number of Events',
        attribute = TreeVariable.fromString( "l1_pt/F" ),
        binning=[30, 0, 500],
    ))

    plots.append(Plot(
        name = "gen_Z1_pt",
        texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
        attribute = TreeVariable.fromString( "genZ1_pt/F" ),
        binning=[30, 0, 700],
    ))

    plots.append(Plot(
        name = "gen_Z_pt_alt",
        texX = 'p_{T}(Z) (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.Z_pt,
        binning=[30, 0, 700],
    ))

    plots.append(Plot(
        name = "gen_m_WZ",
        texX = 'm_{WZ} (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mWZ,
        binning=[30, 0, 1400],
    ))

    plots.append(Plot(
        name = "gen_Z1_pt_bveto",
        texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.genZ1_pt if event.pass_b_veto else float('nan'),
        binning=[30, 0, 700],
    ))

    plots.append(Plot(
        name = "gen_m_WZ_bveto",
        texX = 'm_{WZ} (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mWZ  if event.pass_b_veto else float('nan'),
        binning=[30, 0, 1400],
    ))

    plots.append(Plot(
        name = "gen_Z1_pt_bveto_btag",
        texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.genZ1_pt if (event.pass_b_veto and event.pass_b_tag) else float('nan'),
        binning=[30, 0, 700],
    ))

    plots.append(Plot(
        name = "gen_m_WZ_bveto_btag",
        texX = 'm_{WZ} (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mWZ  if (event.pass_b_veto and event.pass_b_tag) else float('nan'),
        binning=[30, 0, 1400],
    ))

    plots.append(Plot(
        name = "gen_bveto",
        texX = 'pass b veto', texY = 'Number of Events',
        attribute = lambda event, sample: event.pass_b_veto,
        binning=[2, -0.5, 1.5],
    ))

    plots.append(Plot(
        name = "gen_btag",
        texX = 'pass b tag', texY = 'Number of Events',
        attribute = lambda event, sample: event.pass_b_tag,
        binning=[2, -0.5, 1.5],
    ))
    
    plots.append(Plot(
        name = "gen_b_pt",
        texX = 'b p_{T} (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.b_pt,
        binning=[50, 0, 500],
    ))

    plots.append(Plot(
        name = "gen_b_eta",
        texX = 'b #eta', texY = 'Number of Events',
        attribute = lambda event, sample: event.b_eta,
        binning=[50, -5, 5],
    ))

    plots.append(Plot(
        name = "gen_bTop_pt",
        texX = 'b_{top} p_{T} (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.bTop_pt,
        binning=[50, 0, 500],
    ))

    plots.append(Plot(
        name = "gen_bTop_eta",
        texX = 'b_{top} #eta', texY = 'Number of Events',
        attribute = lambda event, sample: event.bTop_eta,
        binning=[50, -5, 5],
    ))

    plots.append(Plot(
        name = "gen_b_diff_pt",
        texX = 'p_{T}(b_{top}) - p_{T}(b)', texY = 'Number of Events',
        attribute = lambda event, sample: event.b_diff_pt,
        binning=[50, -250, 250],
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

    yields[mode]["data"] = 0

    yields[mode]["MC"] = sum(yields[mode][s.name] for s in mc1+mc2+mc3)
    dataMCScale        = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')

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
    plotstowrite = [
        "gen_Z_pt_alt",
        "gen_Z1_pt",
        "gen_W_pt",
        "gen_m_WZ",
        "gen_Z1_pt_bveto",
        "gen_m_WZ_bveto",
        "gen_Z1_pt_bveto_btag",
        "gen_m_WZ_bveto_btag",
        "l1_pt",
        "Z1_pt",
    ]
    # Write Result Hist in root file
    plot_dir = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, "all", args.selection)
    if not os.path.exists(plot_dir):
        try:
            os.makedirs(plot_dir)
        except:
            print 'Could not crete', plot_dir
    outfilename = plot_dir+'/Results.root'
    outfile = ROOT.TFile(outfilename, 'recreate')
    outfile.cd()
    for plot in allPlots['all']:
        if plot.name in plotstowrite:
            for idx, histo_list in enumerate(plot.histos):
                for j, h in enumerate(histo_list):
                    histname = h.GetName()
                    if "TWZ_NLO_DR" in histname: process = "tWZ_NLO_DR"
                    elif "TWZ_NLO_DS" in histname: process = "tWZ_NLO_DS"
                    elif "tWZToLL01j_lepWFilter" in histname: process = "tWZ"
                    elif "TTZ" in histname: process = "ttZ"
                    elif "ttZ01j_lepWFilter" in histname: process = "ttZ"
                    elif "ttZ01j" in histname: process = "ttZ_old"
                    elif "TTX_rare" in histname: process = "ttX"
                    elif "TZQ" in histname: process = "tZq"
                    elif "WZ" in histname: process = "WZ"
                    elif "ZZ" in histname: process = "ZZ"
                    elif "triBoson" in histname: process = "triBoson"
                    elif "nonprompt" in histname: process = "nonprompt"
                    elif "data" in histname: process = "data"
                    h.Write(plot.name+"__"+process)
    outfile.Close()



logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
