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
argParser.add_argument('--plot_directory', action='store', default='BIT_training')

args = argParser.parse_args()

features = []
weights = []
diff_weights = []
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


# Read variables and sequences
sequence       = []

def WriteTrainingData(event,sample):
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
    feature_list.append(event.nBTag)
    feature_list.append(event.nJetGood)
    features.append(feature_list)

    weight = w.get_weight_func(cHq1Re33=0)
    weight_diff_cHq1Re11 = w.get_diff_weight_func('cHq1Re33', cHq1Re33=0)

    weights.append(weight(event,sample))
    diff_weights.append(weight_diff_cHq1Re11(event,sample))
sequence.append(WriteTrainingData)

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
    name = "n_Jet",
    texX = 'Number of jets', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "nJetGood/I" ),
    binning=[11,-0.5,10.5],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = "n_bTags",
    texX = 'Number of b tagged jets', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "nBTag/I" ),
    binning=[11,-0.5,10.5],
    addOverFlowBin='upper',
))


plotting.fill(plots, read_variables = read_variables, sequence = sequence)

for plot in plots:
    for i_h, hl in enumerate(plot.histos):
        # dress up
        hl[0].legendText = params[i_h]['legendText']
        hl[0].style = params[i_h]['style']

drawPlots(plots)

####### BIT training
import sys
sys.path.append('/users/dennis.schwarz/CMSSW_10_6_0/src/BIT/')
from  BoostedInformationTree import BoostedInformationTree

learning_rate = 0.02
n_trees       = 100
learning_rate = 0.2
max_depth     = 11
min_size      = 50

# split into training and validation
# validation_fraction = 0.3
# Nevents = len(weights)
# middle_index = Nevents//(1/validation_fraction) # divide and round down
#
# training_features = features[middle_index:]
# training_weights = weights[middle_index:]
# training_diff_weights = diff_weights[middle_index:]
#
# validation_weights = weights[:middle_index]
# validation_diff_weights = diff_weights[:middle_index]
# validation_features = features[middle_index:]

# for now use all data for training
training_features = features
training_weights = weights
training_diff_weights = diff_weights

bit = BoostedInformationTree(
        training_features = np.array(training_features),
        training_weights      = np.array(training_weights),
        training_diff_weights = np.array(training_diff_weights),
        learning_rate = learning_rate,
        n_trees = n_trees,
        max_depth=max_depth,
        min_size=min_size,
        split_method='vectorized_split_and_weight_sums',
        weights_update_method='vectorized')

bit.boost()
bit.save('BIT_ttZ_cHq1Re33_0.pkl')

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
