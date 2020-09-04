#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT, os
ROOT.gROOT.SetBatch(True)
#from math                                import sqrt, cos, sin, pi, atan2, cosh

# RootTools
from RootTools.core.standard             import *

# tWZ
from tWZ.Tools.user                      import plot_directory
#from tWZ.Tools.cutInterpreter            import cutInterpreter
#from tWZ.Tools.objectSelection           import lepString 
# Analysis
from Analysis.Tools.helpers              import deltaPhi, deltaR
#from Analysis.Tools.puProfileCache       import *
#from Analysis.Tools.puReweighting        import getReweightingFunction
import Analysis.Tools.syncer

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
#argParser.add_argument('--noData',         action='store_true', default=False, help='also plot data?')
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
#argParser.add_argument('--dataMCScaling',  action='store_true', help='Data MC scaling?', )
argParser.add_argument('--plot_directory', action='store', default='tWZ_v3')
argParser.add_argument('--era',            action='store', type=str, default="Run2016")
#argParser.add_argument('--selection',      action='store', default='trilep-minDLmass12-onZ1-njet4p-btag2p')
args = argParser.parse_args()

# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"
#if args.noData:                       args.plot_directory += "_noData"

#logger.info( "Working in era %s", args.era)

#tWZ_sample = TWZ if args.nominalSignal else yt_TWZ_filter
#from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed import *

from Samples.nanoAOD.Run2016_private_nanoAODv6 import *

if args.era == "Run2016":
    data_sample = Sample.combine( "DoubleMuon_Run2016", [
        DoubleMuon_Run2016B_25Oct2019_ver2,
        DoubleMuon_Run2016C_25Oct2019,
        DoubleMuon_Run2016D_25Oct2019,
        DoubleMuon_Run2016E_25Oct2019,
        DoubleMuon_Run2016F_25Oct2019,
        DoubleMuon_Run2016G_25Oct2019,
        DoubleMuon_Run2016H_25Oct2019,
        ])
    triggers    = ["HLT_Mu3_PFJet40", "HLT_Mu8", "HLT_Mu17"]#, "HLT_Mu27"] HLT_Mu27 is actually in SingleMuon!
    #mc = [Summer16.TWZ_NLO_DR, Summer16.TTZ, Summer16.TTX_rare, Summer16.TZQ, Summer16.WZ, Summer16.triBoson, Summer16.ZZ, Summer16.nonprompt_3l]
#elif args.era == "Run2017":
#    mc = [Fall17.TWZ_NLO_DR, Fall17.TTZ, Fall17.TTX_rare, Fall17.TZQ, Fall17.WZ, Fall17.triBoson, Fall17.ZZ, Fall17.nonprompt_3l]
#elif args.era == "Run2018":
#    mc = [Autumn18.TWZ_NLO_DR, Autumn18.TTZ, Autumn18.TTX_rare, Autumn18.TZQ, Autumn18.WZ, Autumn18.triBoson, Autumn18.ZZ, Autumn18.nonprompt_3l]
#elif args.era == "RunII":
#    mc = [TWZ_NLO_DR, TTZ, TTX_rare, TZQ, WZ, triBoson, ZZ, nonprompt_3l]

triggerSelection = '('+"||".join(triggers)+')'

#lumi_scale                 = data_sample.lumi/1000
data_sample.scale          = 1.
#for sample in mc:
#    sample.scale           = 1 # Scale MCs individually with lumi

if args.small:
    for sample in [data_sample]:
        sample.normalization = 1.
        sample.reduceFiles( factor = 10 )
        #sample.reduceFiles( to=1)
        sample.scale /= sample.normalization

# Text on the plots
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)
tex.SetTextAlign(11) # align right

def drawObjects():
    lines = [
      (0.15, 0.95, 'CMS Preliminary'), 
      #(0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines] 

def drawPlots(plots):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, ("log" if log else "lin"))#, args.selection)
    for plot in plots:
      if not max(l.GetMaximum() for l in sum(plot.histos,[])): continue # Empty plot

      _drawObjects = []

      if isinstance( plot, Plot):
          plotting.draw(plot,
            plot_directory = plot_directory_,
            #ratio = {'yRange':(0.1,1.9)} if not args.noData else None,
            logX = False, logY = log, sorting = True,
            yRange = (0.03, "auto") if log else (0.001, "auto"),
            #scaling = {0:1} if args.dataMCScaling else {},
            legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
            drawObjects = drawObjects() + _drawObjects,
            copyIndexPHP = True, extensions = ["png"],
          )
            
# Read variables and sequences
sequence       = []

read_variables = [
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
    ]

allPlots   = {}

data_sample.texName = "data"
data_sample.setSelectionString(triggerSelection)
data_sample.name    = "data"
data_sample.style   = styles.errorStyle(ROOT.kBlack)
#lumi_scale          = data_sample.lumi/1000

#weight_ = lambda event, sample: event.weight if sample.isData else event.weight*lumi_year[event.year]/1000.

#for sample in mc: sample.style = styles.fillStyle(sample.color)
#for sample in mc:
#  sample.read_variables = read_variables_MC 
#  sample.setSelectionString([getLeptonSelection(mode)])
#  sample.weight = lambda event, sample: event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger#*event.reweightLeptonSF

#stack = Stack([data_sample])
#
## Use some defaults
#Plot.setDefaults(stack = stack)
#
#plots = []
#
#plots.append(Plot(
#  name = 'nVtxs', texX = 'vertex multiplicity', texY = 'Number of Events',
#  attribute = TreeVariable.fromString( "PV_npvsGood/I" ),
#  binning=[50,0,50],
#  addOverFlowBin='upper',
#))
#
#plotting.fill(plots, read_variables = read_variables, sequence = sequence)
#
#drawPlots(plots)
