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
   # ('cHq3Re11', 1.0, ROOT.kCyan),
   # ('cHq3Re22', 1.0, ROOT.kMagenta),
   # ('cHq3Re33', 1.0, ROOT.kBlue),
    ('cHq1Re11', 0.0, ROOT.kRed),
    # ('cHq1Re22', 2.0, ROOT.kGreen),
    # ('cHq1Re33', 2.0, ROOT.kOrange),
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
bit_cHq1Re11_0 = BoostedInformationTree.load('/users/dennis.schwarz/CMSSW_10_6_0/src/tWZ/BIT/BIT_ttZ_cHq1Re11_0.pkl')
bit_cHq1Re11_2 = BoostedInformationTree.load('/users/dennis.schwarz/CMSSW_10_6_0/src/tWZ/BIT/BIT_ttZ_cHq1Re11_2.pkl')
bit_cHq1Re33_0 = BoostedInformationTree.load('/users/dennis.schwarz/CMSSW_10_6_0/src/tWZ/BIT/BIT_ttZ_cHq1Re33_0.pkl')
bit_cHq1Re11_0_doublediff = BoostedInformationTree.load('/users/dennis.schwarz/CMSSW_10_6_0/src/tWZ/BIT/BIT_ttZ_cHq1Re11_0_doublediff.pkl')
bins_BIT = [-1., 0.0, 0.2, 0.4, 0.8, 1.] # for plots

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
    feature_list.append(event.nBTag)
    feature_list.append(event.nJetGood)
    event.BIT_cHq1Re11_0 = bit_cHq1Re11_0.predict(np.array(feature_list))
    event.BIT_cHq1Re11_2 = bit_cHq1Re11_2.predict(np.array(feature_list))
    event.BIT_cHq1Re33_0 = bit_cHq1Re33_0.predict(np.array(feature_list))
    event.BIT_cHq1Re11_0_doublediff = bit_cHq1Re11_0_doublediff.predict(np.array(feature_list))
    # print event.BIT
sequence.append(BITpredict)

################################################################################

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I", "minDLmass/F",
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
    "np/I", VectorTreeVariable.fromString("p[C/F]",nMax=200)
]

read_variables_MC = ['reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F']

# define 3l selections
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
    elif mode=='all':
        if args.sample == "ZZ":
            return "Sum$({mu_string})+Sum$({ele_string})==4".format(mu_string=mu_string,ele_string=ele_string)
        else:
            return "Sum$({mu_string})+Sum$({ele_string})==3".format(mu_string=mu_string,ele_string=ele_string)

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
    name = "n_Jet",
    texX = 'Number of jets', texY = 'Number of Events',
    attribute = lambda event, sample: event.nJetGood,
    binning=[11,-0.5,10.5],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = "n_btag",
    texX = 'Number of jets', texY = 'Number of Events',
    attribute = lambda event, sample: event.nBTag,
    binning=[11,-0.5,10.5],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'BIT_cHq1Re11_0_A', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT_cHq1Re11_0,
    binning=[50,-5,5],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'BIT_cHq1Re11_0_B', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT_cHq1Re11_0,
    binning=[50,-1.0,1.0],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'BIT_cHq1Re11_0_C', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT_cHq1Re11_0,
    binning=[50,-0.1,0.1],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'BIT_cHq1Re11_0_doublediff', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT_cHq1Re11_0_doublediff,
    binning=[50,-10,100],
    addOverFlowBin='upper',
))




plots.append(Plot(
    name = 'BIT_cHq1Re11_2_A', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT_cHq1Re11_2,
    binning=[50,-1,1],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'BIT_cHq1Re11_2_A2', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT_cHq1Re11_2,
    binning=[10,-1,1],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'BIT_cHq1Re11_2_B', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT_cHq1Re11_2,
    binning=[50,-0.5,0.5],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'BIT_cHq1Re11_2_C', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT_cHq1Re11_2,
    binning=[50,-0.1,0.1],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'BIT_cHq1Re33_0_A', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT_cHq1Re33_0,
    binning=[50,-1,1],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'BIT_cHq1Re33_0_B', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT_cHq1Re33_0,
    binning=[50,-0.5,0.5],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'BIT_cHq1Re33_0_C', texX = 'BIT score', texY = 'Number of Events',
    attribute = lambda event, sample: event.BIT_cHq1Re33_0,
    binning=[50,-0.1,0.1],
    addOverFlowBin='upper',
))

### binned in BIT score
# Z pt
plots.append(Plot(
    name = "Z1_pt_bin_1",
    texX = 'p_{T}(Z_{1}) (-1.0 < score < 0.2) [GeV]', texY = 'Number of Events / 20 GeV',
    attribute = lambda event, sample: event.Z1_pt if event.BIT_cHq1Re11_2 < 0.2 else float('nan'),
    binning=[30,0,600],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = "Z1_pt_bin_2",
    texX = 'p_{T}(Z_{1}) (0.2 < score < 0.6) [GeV]', texY = 'Number of Events / 20 GeV',
    attribute = lambda event, sample: event.Z1_pt if event.BIT_cHq1Re11_2 > 0.2 and event.BIT_cHq1Re11_2 < 0.6 else float('nan'),
    binning=[30,0,600],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = "Z1_pt_bin_3",
    texX = 'p_{T}(Z_{1}) (0.6 < score < 1.0) [GeV]', texY = 'Number of Events / 20 GeV',
    attribute = lambda event, sample: event.Z1_pt if event.BIT_cHq1Re11_2 > 0.6 else float('nan'),
    binning=[30,0,600],
    addOverFlowBin='upper',
))

# Z eta
plots.append(Plot(
    name = "Z1_eta_bin_1",
    texX = '#eta(Z_{1}) (-1.0 < score < 0.2) [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.Z1_eta if event.BIT_cHq1Re11_2 < 0.2 else float('nan'),
    binning=[20, -3, 3],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = "Z1_eta_bin_2",
    texX = '#eta(Z_{1}) (0.2 < score < 0.6) [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.Z1_eta if event.BIT_cHq1Re11_2 > 0.2 and event.BIT_cHq1Re11_2 < 0.6 else float('nan'),
    binning=[20, -3, 3],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = "Z1_eta_bin_3",
    texX = '#eta(Z_{1}) (0.6 < score < 1.0) [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.Z1_eta if event.BIT_cHq1Re11_2 > 0.6 else float('nan'),
    binning=[20, -3, 3],
    addOverFlowBin='upper',
))


    # plots.append(Plot(
    #     name = "l1_Z_pt_bin"+str(i),
    #     texX = 'p_{T}(l_{1}^Z) '+score_string+' [GeV]', texY = 'Number of Events / 20 GeV',
    #     attribute = lambda event, sample: event.lep_pt[event.Z1_l1_index] if event.BIT_cHq1Re11_2>lo and event.BIT_cHq1Re11_2<hi else -1,
    #     binning=[20,0,400],
    #     addOverFlowBin='upper',
    # ))
    # plots.append(Plot(
    #     name = "n_Jet_bin"+str(i),
    #     texX = 'Number of jets '+score_string, texY = 'Number of Events',
    #     attribute = lambda event, sample: event.nJetGood if event.BIT_cHq1Re11_2>lo and event.BIT_cHq1Re11_2<hi else -1,
    #     binning=[11,-0.5,10.5],
    #     addOverFlowBin='upper',
    # ))
    # plots.append(Plot(
    #     name = "n_btag_bin"+str(i),
    #     texX = 'Number of jets '+score_string, texY = 'Number of Events',
    #     attribute = lambda event, sample: event.nBTag if event.BIT_cHq1Re11_2>lo and event.BIT_cHq1Re11_2<hi else -1,
    #     binning=[11,-0.5,10.5],
    #     addOverFlowBin='upper',
    # ))

plotting.fill(plots, read_variables = read_variables, sequence = sequence)

for plot in plots:
    for i_h, hl in enumerate(plot.histos):
        # dress up
        hl[0].legendText = params[i_h]['legendText']
        hl[0].style = params[i_h]['style']

drawPlots(plots)

# getCoeffList
selString_ptZ = cutInterpreter.cutString(args.selection)+"&&Z1_pt>=500"

coeffList_ptZ = w.getCoeffListFromDraw(sample, selectionString=selString_ptZ, weightString = "weight*137" )
values = [-4, -3, -2, -1, 1, 2, 3, 4]
h_cHq1Re11_ptZ = ROOT.TH1F("yields_cHq1Re11_ptZ", "cHq1Re11", len(values), values[0]-0.5, values[-1]+0.5)
# h_cHq1Re33_ptZ = ROOT.TH1F("yields_cHq1Re33_ptZ", "cHq1Re33", len(values), values[0]-0.5, values[-1]+0.5)

for i in range(len(values)):
    h_cHq1Re11_ptZ.SetBinContent(i+1, w.get_weight_yield(coeffList_ptZ, cHq1Re11=values[i]))
    # h_cHq1Re33_ptZ.SetBinContent(i+1, w.get_weight_yield(coeffList_ptZ, cHq1Re33=values[i]))

dir = os.path.join(plot_directory, 'eft',  args.plot_directory, sample.name, "lin", args.selection)

c3 = ROOT.TCanvas()
h_cHq1Re11_ptZ.Draw("HIST")
h_cHq1Re11_ptZ.Fit('pol2')
c3.Print(dir+"/yields_cHq1Re11_ptZ500.pdf")

# c4 = ROOT.TCanvas()
# h_cHq1Re33_ptZ.Draw("HIST")
# h_cHq1Re33_ptZ.Fit('pol2')
# c4.Print(dir+"/yields_cHq1Re33_ptZ500.pdf")

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
