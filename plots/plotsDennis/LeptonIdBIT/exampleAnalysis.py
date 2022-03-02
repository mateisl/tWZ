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

import Analysis.Tools.syncer
import numpy as np

################################################################################
# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory', action='store', default='LeptonID')
argParser.add_argument('--selection',      action='store', default='trilep')
argParser.add_argument('--era',            action='store', type=str, default="Run2018")
args = argParser.parse_args()

################################################################################
# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

################################################################################
# Define the MC samples
sample_directory = "/scratch-cbe/users/dennis.schwarz/LeptonID/nanoTuples/2018/"
TTZ = Sample.fromDirectory(name="TTZ", treeName="Events", isData=False, color=ROOT.kAzure, texName="ttZ", directory=sample_directory+"TTZ")

mc = [TTZ]
lumi_scale = 60

################################################################################
# Text on the plots
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)
tex.SetTextAlign(11) # align right

################################################################################
# Functions needed specifically for this analysis routine

def drawObjects( plotData, lumi_scale ):
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if plotData else 'CMS Simulation'),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) '% ( lumi_scale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines]

def drawPlots(plots, mode):
    for log in [False, True]:
        plot_directory_ = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, mode+ ("_log" if log else ""), args.selection)
        for plot in plots:
            if not max(l.GetMaximum() for l in sum(plot.histos,[])): continue # Empty plot

            _drawObjects = []
            n_stacks=len(plot.histos)
            plotData=False
            if isinstance( plot, Plot):
                plotting.draw(plot,
                  plot_directory = plot_directory_,
                  ratio =  None,
                  logX = False, logY = log, sorting = True,
                  yRange = (0.03, "auto") if log else (0.001, "auto"),
                  scaling = {},
                  legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
                  drawObjects = drawObjects( plotData , lumi_scale ) + _drawObjects,
                  copyIndexPHP = True, extensions = ["png", "pdf", "root"],
                )

################################################################################
# Define sequences
sequence       = []

def mysequence( event, sample ):
    print event.l1_pt
sequence.append( mysequence )

################################################################################
# Read variables

read_variables = [
    "weight/F", "year/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F",
    "l2_pt/F", "l2_eta/F" , "l2_phi/F",
    "l3_pt/F", "l3_eta/F" , "l3_phi/F",
    "Muon[pt/F,eta/F,phi/F,pdgId/I,mediumId/O]",
    "Electron[pt/F,eta/F,phi/F,pdgId/I,mvaFall17V2Iso_WP80/O]",
    
]

read_variables_MC = [
    # 'reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F',
    # "genZ1_pt/F", "genZ1_eta/F", "genZ1_phi/F",
]


################################################################################
# define 3l selections
mu_string  = "Muon_pt>10&&abs(Muon_eta)<2.4&&abs(Muon_pdgId)==13&&Muon_mediumId"
ele_string  = "Electron_pt>10&&abs(Electron_eta)<2.4&&abs(Electron_pdgId)==11&&mvaFall17V2Iso_WP80"

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
allPlots   = {}
allModes   = ['mumumu']
for i_mode, mode in enumerate(allModes):

    weight_ = lambda event, sample: event.weight if sample.isData else event.weight#*lumi_year[event.year]/1000.

    for sample in mc: sample.style = styles.fillStyle(sample.color)

    for sample in mc:
      sample.read_variables = read_variables_MC
      sample.setSelectionString([getLeptonSelection(mode)])
      sample.weight = lambda event, sample: 1#*event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger*event.reweightLeptonSF

    stack = Stack(mc)

    # Use some defaults
    Plot.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = cutInterpreter.cutString(args.selection))
    
    ################################################################################
    # Now define the plots
    
    plots = []
    
    plots.append(Plot(
        name = "Lepton1_pt",
        texX = 'Leading lepton p_{T}', texY = 'Number of Events',
        attribute = lambda event, sample: event.l1_pt,
        binning=[40, 0., 200.],
    ))
    
    plotting.fill(plots, read_variables = read_variables, sequence = sequence)

    drawPlots(plots, mode)
    
logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
