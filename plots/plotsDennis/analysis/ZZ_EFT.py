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
argParser.add_argument('--SMEFTsim',     action='store_true', default=False, help='use FastSim samples?')
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--dataMCScaling',  action='store_true', help='Data MC scaling?', )
argParser.add_argument('--plot_directory', action='store', default='ZZ_EFT_v1')
argParser.add_argument('--era',            action='store', type=str, default="Run2018")
argParser.add_argument('--selection',      action='store', default='trilepVL-minDLmass12-onZ1-onZ2-nLeptons4')
args = argParser.parse_args()

################################################################################
# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"
if args.noData:                       args.plot_directory += "_noData"
if args.SMEFTsim:                     args.plot_directory += "_SMEFTsim"

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
    sample_eft = SMEFTsim_fast.ZZ
    if args.SMEFTsim:
        mc = [Autumn18.TWZ_NLO_DR, Autumn18.TTZ, Autumn18.TTX_rare, Autumn18.TZQ, Autumn18.WZ, Autumn18.triBoson, SMEFTsim_fast.ZZ, Autumn18.nonprompt_3l]
elif args.era == "RunII":
    mc = [TWZ_NLO_DR, TTZ, TTX_rare, TZQ, WZ, triBoson, ZZ, nonprompt_3l]

signal_idx = 6 # INDEX OF ZZ sample
################################################################################
# EFT reweight

# WeightInfo
w = WeightInfo(sample_eft.reweight_pkl)
w.set_order(2)

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
WCs = []
WC_setup = [
    ('cHq1Re11', ROOT.kRed),
    ('cHq1Re22', ROOT.kGreen+2),
    ('cHq1Re33', ROOT.kOrange-3),
]
for i_wc, (WCname, color) in enumerate(WC_setup):
    for i in range(Npoints):
        value = minval + ((maxval-minval)/(Npoints-1))*i
        WCs.append( (WCname, value, color) )

params =  [ ]

for i_wc, (WC, WCval, color) in enumerate(WCs):
    params.append({'legendText':'%s=%3.2f'%(WC, WCval), 'color':color,  'WC':{WC:WCval} })
# params.append({'legendText':'SM', 'color':ROOT.kBlack,  'WC':{} })

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

# Set up MC sample
for sample in mc:
    sample.scale           = 1 # Scale MCs individually with lumi

for param in params:
    param['sample'].scale = 1

if args.small:
    for sample in mc + [data_sample]:
        sample.normalization = 1.
        sample.reduceFiles( to = 1 )
        #sample.reduceFiles( to=1)
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
            n_stacks=len(plot.histos)
            if isinstance( plot, Plot):
                plotting.draw(plot,
                  plot_directory = plot_directory_,
                  ratio = {'histos': [[i+1,0] for i in range(n_stacks-1)], 'yRange':(0.1,1.9)},
                  logX = False, logY = log, sorting = True,
                  yRange = (0.03, "auto") if log else (0.001, "auto"),
                  scaling = {0:1} if args.dataMCScaling else {},
                  legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
                  drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ) + _drawObjects,
                  copyIndexPHP = True, extensions = ["png"],
                )


################################################################################
# Define sequences
sequence       = []


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

read_variables_MC = [
    'reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F',
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
mu_string  = lepString('mu','VL') + "&&lep_mediumId"
ele_string = lepString('ele','VL')
def getLeptonSelection( mode ):
    if   mode=="mumumumu": return "Sum$({mu_string})==4&&Sum$({ele_string})==0".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="mumuee":  return "Sum$({mu_string})==2&&Sum$({ele_string})==2".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="eeee":   return "Sum$({mu_string})==0&&Sum$({ele_string})==4".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=='all':    return "Sum$({mu_string})+Sum$({ele_string})==4".format(mu_string=mu_string,ele_string=ele_string)

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
allModes   = ['mumumumu','mumuee','eeee']
for i_mode, mode in enumerate(allModes):
    yields[mode] = {}
    if not args.noData:
        data_sample.texName = "data"
        data_sample.setSelectionString([getLeptonSelection(mode)])
        data_sample.name           = "data"
        data_sample.style          = styles.errorStyle(ROOT.kBlack)
        lumi_scale                 = data_sample.lumi/1000

    # weight_ = lambda event, sample: event.weight if sample.isData else event.weight*lumi_year[event.year]/1000.

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
      texX = 'N_{jets}', texY = 'Number of Events',
      attribute = TreeVariable.fromString( "nJetGood/I" ), #nJetSelected
      binning=[8,-0.5,7.5],
    ))

    plots.append(Plot(
      texX = 'N_{leptons}', texY = 'Number of Events',
      attribute = TreeVariable.fromString( "nGoodLeptons/I" ),
      binning=[8,-0.5,7.5],
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
        name = "Z1_pt",
        texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "Z1_pt/F" ),
        binning=[20,0,400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "Z2_pt",
        texX = 'p_{T}(Z_{2}) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "Z2_pt/F" ),
        binning=[20,0,400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "Z1_mass",
        texX = 'm(Z_{1}) (GeV)', texY = 'Number of Events / 10 GeV',
        attribute = TreeVariable.fromString( "Z1_mass/F" ),
        binning=[20,0,200],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "Z2_mass",
        texX = 'm(Z_{2}) (GeV)', texY = 'Number of Events / 10 GeV',
        attribute = TreeVariable.fromString( "Z2_mass/F" ),
        binning=[20,0,200],
        addOverFlowBin='upper',
    ))




    plotting.fill(plots, read_variables = read_variables, sequence = sequence)

    ################################################################################
    # Get normalization yields from yield histogram
    for plot in plots:
      if plot.name == "yield":
        for i, l in enumerate(plot.histos):
          for j, h in enumerate(l):
            yields[mode][plot.stack[i][j].name] = h.GetBinContent(h.FindBin(0.5+i_mode))
            h.GetXaxis().SetBinLabel(1, "#mu#mu#mu#mu")
            h.GetXaxis().SetBinLabel(2, "#mu#muee")
            h.GetXaxis().SetBinLabel(3, "eeee")

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
            # Now also add other SM samples (backgrounds)
            for idx_sm, histo_list_sm in enumerate(plot.histos):
                # Only add SM MC
                if idx_sm == 0:
                    for i in range(len(histo_list_sm)):
                        # Also don't add the signal again
                        if i == signal_idx:
                            continue
                        histo_list[0].Add(histo_list_sm[i])

    drawPlots(plots, mode, dataMCScale)

    allPlots[mode] = plots


################################################################################
# Add the different channels into SF and all
for mode in ["comb1","all"]:
    yields[mode] = {}
    for y in yields[allModes[0]]:
        try:    yields[mode][y] = sum(yields[c][y] for c in ['eeee','mumuee','mumumumu'])
        except: yields[mode][y] = 0
    dataMCScale = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')

    for plot in allPlots['mumumumu']:
        if mode=="comb1":
            tmp = allPlots['mumuee']
        else:
            tmp = allPlots['eeee']
        for plot2 in (p for p in tmp if p.name == plot.name):
            for i, j in enumerate(list(itertools.chain.from_iterable(plot.histos))):
                for k, l in enumerate(list(itertools.chain.from_iterable(plot2.histos))):
                    if i==k:
                        j.Add(l)

    if mode == "all":
        drawPlots(allPlots['mumumumu'], mode, dataMCScale)
        # Write Result Hist in root file
        plot_dir = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, mode, args.selection)
        outfile = ROOT.TFile(plot_dir+'/Results.root', 'recreate')
        outfile.cd()
        for plot in allPlots['mumumumu']:
            if plot.name == "Z1_pt":
                for idx, histo_list in enumerate(plot.histos):
                    for j, h in enumerate(histo_list):
                        histname = h.GetName()
                        if "TWZ_NLO_DR" in histname: process = "tWZ"
                        elif "TTZ" in histname: process = "ttZ"
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
