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
argParser.add_argument('--plot_directory', action='store', default='tWZ_CUTFLOW_v1')
argParser.add_argument('--selection',      action='store', default='trilepT-minDLmass12-onZ1-njet4p-btag1')
argParser.add_argument('--era',            action='store', type=str, default="Run2016")
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
def drawObjects( plotData, dataMCScale, lumi_scale ):
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if plotData else 'CMS Simulation'),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    if "mt2ll100" in args.selection and args.noData: lines += [(0.55, 0.5, 'M_{T2}(ll) > 100 GeV')] # Manually put the mt2ll > 100 GeV label
    return [tex.DrawLatex(*l) for l in lines]

def drawPlots(plots, mode, dataMCScale):
  for log in [False]:
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
            copyIndexPHP = True, extensions = ["png","root"],
          )
################################################################################
# Define sequences
sequence       = []


################################################################################
# Read variables

read_variables = [
    "weight/F", "year/I", "PV_npvsGood/I"
]

read_variables_MC = [
    'reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F',
]

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
hists = {}

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
      # sample.weight = lambda event, sample: event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger#*event.reweightLeptonSF

    #yt_TWZ_filter.scale = lumi_scale * 1.07314

    nominalMCWeights = ["weight", "reweightPU", "reweightTrigger", "reweightBTag_SF", "reweightL1Prefire"]
    weightstring = "*".join(nominalMCWeights)
    lumistring = "("+"+".join( [ "%s*(year==%s)"%(l/1000, y) for (y,l) in lumi_year.iteritems() ])+")"
    weightstring += "*"+lumistring

    weightstring_data = "weight"

    Nevents = {}
    for sample in mc:
        yieldlist = [
            sample.getYieldFromDraw(selectionString = None, weightString = weightstring)['val'],
            sample.getYieldFromDraw(selectionString = cutInterpreter.cutString("trilepT"), weightString = weightstring)['val'],
            sample.getYieldFromDraw(selectionString = cutInterpreter.cutString("trilepT-minDLmass12"), weightString = weightstring)['val'],
            sample.getYieldFromDraw(selectionString = cutInterpreter.cutString("trilepT-minDLmass12-onZ1"), weightString = weightstring)['val'],
            sample.getYieldFromDraw(selectionString = cutInterpreter.cutString("trilepT-minDLmass12-onZ1-njet4p"), weightString = weightstring)['val'],
            sample.getYieldFromDraw(selectionString = cutInterpreter.cutString("trilepT-minDLmass12-onZ1-njet4p-btag1"), weightString = weightstring)['val']
        ]
        Nevents[sample.name] = yieldlist

    yieldlist = [
        data_sample.getYieldFromDraw(selectionString = None, weightString = weightstring_data)['val'],
        data_sample.getYieldFromDraw(selectionString = cutInterpreter.cutString("trilepT"), weightString = weightstring_data)['val'],
        data_sample.getYieldFromDraw(selectionString = cutInterpreter.cutString("trilepT-minDLmass12"), weightString = weightstring_data)['val'],
        data_sample.getYieldFromDraw(selectionString = cutInterpreter.cutString("trilepT-minDLmass12-onZ1"), weightString = weightstring_data)['val'],
        data_sample.getYieldFromDraw(selectionString = cutInterpreter.cutString("trilepT-minDLmass12-onZ1-njet4p"), weightString = weightstring_data)['val'],
        data_sample.getYieldFromDraw(selectionString = cutInterpreter.cutString("trilepT-minDLmass12-onZ1-njet4p-btag1"), weightString = weightstring_data)['val']
    ]
    Nevents[data_sample.name] = yieldlist

    for sample in mc:
        Ncuts = len(Nevents[sample.name])
        hist = ROOT.TH1F("yield_"+sample.name+"_"+mode, "cut", Ncuts, 0, Ncuts)
        for i in range(Ncuts):
            hist.SetBinContent(i+1, Nevents[sample.name][i])
        hists[sample.name+"_"+mode] = hist

    Ncuts = len(Nevents[data_sample.name])
    hist = ROOT.TH1F("yield_"+data_sample.name+"_"+mode, "cut", Ncuts, 0, Ncuts)
    for i in range(Ncuts):
        hist.SetBinContent(i+1, Nevents[data_sample.name][i])
    hists[data_sample.name+"_"+mode] = hist

################################################################################
# Write file
outfile = ROOT.TFile('Cutflow.root', 'recreate')
outfile.cd()
hist_background = hists["TTZ_mumumu"].Clone("yield_background_all")
hist_background.Reset()
for name in ["TWZ_NLO_DR","TTZ","TTX_rare","TZQ","WZ","triBoson","ZZ","nonprompt_3l","data"]:
    hist_tot = hists[name+"_mumumu"].Clone("yield_"+name+"_all")
    hist_tot.Reset()
    for mode in ['mumumu','mumue','muee', 'eee']:
        hist_tot.Add(hists[name+"_"+mode])
        hists[name+"_"+mode].Write()
        if name != "data" and name != "TWZ_NLO_DR":
            hist_background.Add(hists[name+"_"+mode])
    hist_tot.Write()
hist_background.Write()
outfile.Close()


logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
