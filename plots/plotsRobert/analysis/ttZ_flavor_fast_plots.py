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
argParser.add_argument('--selection',      action='store', default='trilepL-minDLmass12-onZ1-njet3p-btag1p')
argParser.add_argument('--plot_directory', action='store', default='SMEFTsim_fast')

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
    #sample.normalization = 1.
    sample.reduceFiles( to = 1 )
    #sample.reduceFiles( to=1)
    #sample.scale /= sample.normalization

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
    plot_directory_ = os.path.join(plot_directory,  args.plot_directory, sample.name, ("log" if log else "lin"), args.selection)
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

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I", 
    "l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
#    "l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
    "JetGood[pt/F,eta/F,phi/F]",
    "lep[pt/F,eta/F,phi/F,pdgId/I,muIndex/I,eleIndex/I]",
    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I", 
    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
    "np/I", VectorTreeVariable.fromString("p[C/F]",nMax=200)
]

read_variables_MC = ['reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F']
# define 3l selections

##MVA
#import TMB.MVA.configs as configs
#config = configs.ttZ_3l_flavor
#read_variables += config.read_variables
#
## Add sequence that computes the MVA inputs
#def make_mva_inputs( event, sample ):
#    for mva_variable, func in config.mva_variables:
#        setattr( event, mva_variable, func(event, sample) )
#sequence.extend( config.sequence )
#sequence.append( make_mva_inputs )

# load models
from keras.models import load_model

#if args.onlyMVA is not None:
#    has_lstm = ('LSTM' in args.onlyMVA)
#    name = args.onlyMVA.split('/')[-4]
#    models = [ (name, has_lstm, load_model(args.onlyMVA) ), ]
#else:
models = [
#    ("ttZ_3l_flavor", False,  load_model("/mnt/hephy/cms/robert.schoefbeck/TMB/models/ttZ_3l_flavor/ttZ_3l_flavor/multiclass_model.h5")),
#    ("ttZ_3l_flavor_LSTM", True,  load_model("/mnt/hephy/cms/robert.schoefbeck/TMB/models/ttZ_3l_flavor_LSTM/ttZ_3l_flavor/multiclass_model.h5")),
]

def keras_predict( event, sample ):

    # get model inputs assuming lstm
    flat_variables, lstm_jets = config.predict_inputs( event, sample, jet_lstm = True)
    for name, has_lstm, model in models:
        #print has_lstm, flat_variables, lstm_jets
        prediction = model.predict( flat_variables if not has_lstm else [flat_variables, lstm_jets] )
        for i_val, val in enumerate( prediction[0] ):
            setattr( event, name+'_'+config.training_samples[i_val].name, val)

#sequence.append( keras_predict )

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

## 3l trainign variables
#def make_training_observables_3l(event, sample):
#
#    event.nonZ1l1_Z1_deltaPhi = deltaPhi(event.lep_phi[event.nonZ1_l1_index], event.Z1_phi)
#    event.nonZ1l1_Z1_deltaEta = abs(event.lep_eta[event.nonZ1_l1_index] - event.Z1_eta)
#    event.nonZ1l1_Z1_deltaR   = deltaR({'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
#    event.jet0_Z1_deltaR      = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
#    event.jet0_nonZ1l1_deltaR = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
#    event.jet1_Z1_deltaR      = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
#    event.jet1_nonZ1l1_deltaR = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
#
#sequence.append( make_training_observables_3l )

# Use some defaults
Plot.setDefaults(stack = stack, weight = weight, selectionString = "("+cutInterpreter.cutString(args.selection)+")&&("+getLeptonSelection('all')+")")

plots = []

for name, has_lstm, model in models:
    #print has_lstm, flat_variables, lstm_jets
    for i_tr_s, tr_s in enumerate( config.training_samples ):
        disc_name = name+'_'+config.training_samples[i_tr_s].name
        plots.append(Plot(
            texX = disc_name, texY = 'Number of Events',
            name = disc_name, attribute = lambda event, sample, disc_name=disc_name: getattr( event, disc_name ),
            binning=[50, 0, 1],
        ))
    for i_tr_s, tr_s in enumerate( config.training_samples ):
        disc_name = name+'_'+config.training_samples[i_tr_s].name
        plots.append(Plot(
            texX = disc_name, texY = 'Number of Events',
            name = 'coarse_'+disc_name, attribute = lambda event, sample, disc_name=disc_name: getattr( event, disc_name ),
            binning=[10, 0, 1],
        ))

#for mva in mvas:
#    plots.append(Plot(
#        texX = 'MVA_{3l}', texY = 'Number of Events',
#        name = mva['name'], attribute = discriminator_getter(mva['name']),
#        binning=[25, 0, 1],
#    ))

#for mva in mvas:
#    plots.append(Plot(
#        texX = 'MVA_{3l}', texY = 'Number of Events',
#        name = mva['name']+'_coarse', attribute = discriminator_getter(mva['name']),
#        binning=[10, 0, 1],
#    ))

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

#    plots.append(Plot(
#        name = 'l3_pt',
#        texX = 'p_{T}(l_{3}) (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = lambda event, sample:event.l3_pt,
#        binning=[15,0,300],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'l3_eta',
#        texX = '#eta(l_{3})', texY = 'Number of Events',
#        attribute = lambda event, sample: event.l3_eta,
#        binning=[20,-3,3],
#    ))
#
#    plots.append(Plot(
#        name = 'l3_mvaTOP',
#        texX = 'MVA_{TOP}(l_{3})', texY = 'Number of Events',
#        attribute = lambda event, sample: event.l3_mvaTOP,
#        binning=[20,-1,1],
#    ))
#
#    plots.append(Plot(
#        name = 'l3_mvaTOPWP',
#        texX = 'MVA_{TOP}(l_{1}) WP', texY = 'Number of Events',
#        attribute = lambda event, sample: event.l3_mvaTOPWP,
#        binning=[5,0,5],
#    ))

plots.append(Plot(
    texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
    attribute = TreeVariable.fromString( "met_pt/F" ),
    binning=[400/20,0,400],
    addOverFlowBin='upper',
))

#    plots.append(Plot(
#        texX = '#phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
#        attribute = TreeVariable.fromString( "met_phi/F" ),
#        binning=[10,-pi,pi],
#    ))

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
    name = 'Z1_pt_coarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events / 50 GeV',
    attribute = TreeVariable.fromString( "Z1_pt/F" ),
    binning=[16,0,800],
))

plots.append(Plot(
    name = 'Z1_pt_superCoarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "Z1_pt/F" ),
    binning=[3,0,600],
))

#    plots.append(Plot(
#        name = "M3l",
#        texX = 'M(3l) (GeV)', texY = 'Number of Events',
#        attribute = lambda event, sample:event.M3l,
#        binning=[25,0,500],
#    ))
#
#    plots.append(Plot(
#        name = "dPhiZJet",
#        texX = '#Delta#phi(Z,j1)', texY = 'Number of Events',
#        attribute = lambda event, sample: deltaPhi(event.Z1_phi, event.JetGood_phi[0]),
#        binning=[20,0,pi],
#    ))
#
#    plots.append(Plot(
#        name = "l1_Z1_pt",
#        texX = 'p_{T}(l_{1,Z}) (GeV)', texY = 'Number of Events / 10 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.Z1_l1_index],
#        binning=[30,0,300],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = "l1_Z1_pt_coarse",
#        texX = 'p_{T}(l_{1,Z}) (GeV)', texY = 'Number of Events / 40 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.Z1_l1_index],
#        binning=[10,0,400],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'l1_Z1_pt_ext', texX = 'p_{T}(l_{1,Z}) (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.Z1_l1_index],
#        binning=[20,40,440],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = "l2_Z1_pt",
#        texX = 'p_{T}(l_{2,Z}) (GeV)', texY = 'Number of Events / 10 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.Z1_l2_index],
#        binning=[20,0,200],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#      texX = 'p_{T}(leading l) (GeV)', texY = 'Number of Events / 20 GeV',
#      name = 'lep1_pt', attribute = lambda event, sample: event.lep_pt[0],
#      binning=[400/20,0,400],
#    ))
#
#    plots.append(Plot(
#      texX = 'p_{T}(subleading l) (GeV)', texY = 'Number of Events / 10 GeV',
#      name = 'lep2_pt', attribute = lambda event, sample: event.lep_pt[1],
#      binning=[200/10,0,200],
#    ))
#
#    plots.append(Plot(
#      texX = 'p_{T}(trailing l) (GeV)', texY = 'Number of Events / 10 GeV',
#      name = 'lep3_pt', attribute = lambda event, sample: event.lep_pt[2],
#      binning=[150/10,0,150],
#    ))
#
#    plots.append(Plot(
#        name = "l2_Z1_pt_coarse",
#        texX = 'p_{T}(l_{2,Z}) (GeV)', texY = 'Number of Events / 10 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.Z1_l2_index],
#        binning=[10,0,200],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'l2_Z1_pt_ext', texX = 'p_{T}(l_{2,Z}) (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.Z1_l2_index],
#        binning=[20,0,400],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'lnonZ1_pt',
#        texX = 'p_{T}(l_{1,extra}) (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.nonZ1_l1_index],
#        binning=[15,0,300],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'lnonZ1_pt_coarse',
#        texX = 'p_{T}(l_{1,extra}) (GeV)', texY = 'Number of Events / 60 GeV',
#        attribute = lambda event, sample:event.lep_pt[event.nonZ1_l1_index],
#        binning=[3,0,180],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'lnonZ1_charge',
#        texX = 'Charge(l_{1,extra})', texY = 'Number of Events',
#        attribute = lambda event, sample:-event.lep_pdgId[event.nonZ1_l1_index]/abs(event.lep_pdgId[event.nonZ1_l1_index]),
#        binning=[2,-1,1],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'lnonZ1_eta',
#        texX = '#eta(l_{1,extra})', texY = 'Number of Events',
#        attribute = lambda event, sample: event.lep_eta[event.nonZ1_l1_index],
#        binning=[20,-3,3],
#    ))
#
#    plots.append(Plot(
#        texX = 'M(ll) (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = TreeVariable.fromString( "Z1_mass/F" ),
#        binning=[10,81,101],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = "Z1_mass_wide",
#        texX = 'M(ll) (GeV)', texY = 'Number of Events / 2 GeV',
#        attribute = TreeVariable.fromString( "Z1_mass/F" ),
#        binning=[50,20,120],
#        addOverFlowBin='upper',
#    )) 
#
#    plots.append(Plot(
#        name = "Z1_cosThetaStar", texX = 'cos#theta(l-)', texY = 'Number of Events / 0.2',
#        attribute = lambda event, sample:event.Z1_cosThetaStar,
#        binning=[10,-1,1],
#    ))
#
#    plots.append(Plot(
#        name = "Z2_mass_wide",
#        texX = 'M(ll) of 2nd OSDL pair', texY = 'Number of Events / 2 GeV',
#        attribute = TreeVariable.fromString( "Z2_mass/F" ),
#        binning=[60,0,120],
#        addOverFlowBin='upper',
#    )) 
#
#    plots.append(Plot(
#        name = "minDLmass",
#        texX = 'min mass of all DL pairs', texY = 'Number of Events / 2 GeV',
#        attribute = TreeVariable.fromString( "minDLmass/F" ),
#        binning=[60,0,120],
#        addOverFlowBin='upper',
#    )) 
#
#    plots.append(Plot(
#        texX = '#Delta#phi(Z_{1}(ll))', texY = 'Number of Events',
#        attribute = TreeVariable.fromString( "Z1_lldPhi/F" ),
#        binning=[10,0,pi],
#    ))
#
#    plots.append(Plot(
#        texX = '#Delta R(Z_{1}(ll))', texY = 'Number of Events',
#        attribute = TreeVariable.fromString( "Z1_lldR/F" ),
#        binning=[10,0,6],
#    ))
#
#    plots.append(Plot(
#      texX = 'N_{jets}', texY = 'Number of Events',
#      attribute = TreeVariable.fromString( "nJetGood/I" ), #nJetSelected
#      binning=[8,-0.5,7.5],
#    ))
#
#    plots.append(Plot(
#      texX = 'N_{b-tag}', texY = 'Number of Events',
#      attribute = TreeVariable.fromString( "nBTag/I" ), #nJetSelected
#      binning=[4,-0.5,3.5],
#    ))
#
#    plots.append(Plot(
#      texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events / 30 GeV',
#      name = 'jet0_pt', attribute = lambda event, sample: event.JetGood_pt[0],
#      binning=[600/30,0,600],
#    ))
#
#    plots.append(Plot(
#      texX = 'p_{T}(subleading jet) (GeV)', texY = 'Number of Events / 30 GeV',
#      name = 'jet1_pt', attribute = lambda event, sample: event.JetGood_pt[1],
#      binning=[600/30,0,600],
#    ))
#
#    plots.append(Plot(
#        name = "W_pt",
#        texX = 'p_{T}(W) (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = lambda event, sample:event.W_pt,
#        binning=[20,0,400],
#    ))
#
#    # 3l training variables
#
#    plots.append(Plot(
#      texX = '#Delta\#phi(nonZ-l_{1}, Z_{1})', texY = 'Number of Events',
#      name = 'nonZ1l1_Z1_deltaPhi', attribute = lambda event, sample: event.nonZ1l1_Z1_deltaPhi,
#      binning=[20,0,pi],
#    ))
#    plots.append(Plot(
#      texX = '#Delta#eta(nonZ-l_{1}, Z_{1})', texY = 'Number of Events',
#      name = 'nonZ1l1_Z1_deltaEta', attribute = lambda event, sample: event.nonZ1l1_Z1_deltaEta,
#      binning=[20,0,6],
#    ))
#    plots.append(Plot(
#      texX = '#Delta R(nonZ-l_{1}, Z_{1})', texY = 'Number of Event',
#      name = 'nonZ1l1_Z1_deltaR', attribute = lambda event, sample: event.nonZ1l1_Z1_deltaR,
#      binning=[20,0,6],
#    ))
#
#    plots.append(Plot(
#      texX = '#Delta R(jet_{0}, Z_{1})', texY = 'Number of Events',
#      name = 'jet0_Z1_deltaR', attribute = lambda event, sample: event.jet0_Z1_deltaR,
#      binning=[20,0,6],
#    ))
#    plots.append(Plot(
#      texX = '#Delta R(jet_{0}, nonZ-l_{1})', texY = 'Number of Events',
#      name = 'jet0_nonZ1l1_deltaR', attribute = lambda event, sample: event.jet0_nonZ1l1_deltaR,
#      binning=[20,0,6],
#    ))
#    plots.append(Plot(
#      texX = '#Delta R(jet_{1}, Z_{1})', texY = 'Number of Events',
#      name = 'jet1_Z1_deltaR', attribute = lambda event, sample: event.jet1_Z1_deltaR,
#      binning=[20,0,6],
#    ))
#    plots.append(Plot(
#      texX = '#Delta R(jet_{1}, nonZ-l_{1})', texY = 'Number of Events',
#      name = 'jet1_nonZ1l1', attribute = lambda event, sample: event.jet1_nonZ1l1_deltaR,
#      binning=[20,0,6],
#    ))
#    
#    for index in range(3):
#        for abs_pdg in [11, 13]:
#            lep_name = "mu" if abs_pdg==13 else "ele"
#            plots.append(Plot(
#              texX = 'p_{T}(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_pt'%(lep_name, index), attribute = lep_getter("pt", index, abs_pdg),
#              binning=[400/20,0,400],
#            ))
#            plots.append(Plot(
#              texX = '#eta(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_eta'%(lep_name, index), attribute = lep_getter("eta", index, abs_pdg),
#              binning=[30,-3,3],
#            ))
#            plots.append(Plot(
#              texX = '#phi(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_phi'%(lep_name, index), attribute = lep_getter("phi", index, abs_pdg),
#              binning=[30,-pi,pi],
#            ))
#            plots.append(Plot(
#              texX = 'dxy(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_dxy'%(lep_name, index), attribute = lep_getter("dxy", index, abs_pdg, functor = lambda x: abs(x)),
#              binning=[50,0,0.05],
#            ))
#            plots.append(Plot(
#              texX = 'dz(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_dz'%(lep_name, index), attribute = lep_getter("dz", index, abs_pdg, functor = lambda x: abs(x)),
#              binning=[50,0,0.05],
#            ))
#            plots.append(Plot(
#              texX = 'IP_{3D}(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_ip3d'%(lep_name, index), attribute = lep_getter("ip3d", index, abs_pdg, functor = lambda x: abs(x)),
#              binning=[50,0,0.05],
#            ))
#            plots.append(Plot(
#              texX = '#sigma(IP)_{3D}(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_sip3d'%(lep_name, index), attribute = lep_getter("sip3d", index, abs_pdg, functor = lambda x: abs(x)),
#              binning=[40,0,8],
#            ))
#            plots.append(Plot(
#              texX = 'jetRelIso(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_jetRelIso'%(lep_name, index), attribute = lep_getter("jetRelIso", index, abs_pdg),
#              binning=[50,-.15,0.5],
#            ))
#            plots.append(Plot(
#              texX = 'miniPFRelIso_all(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_miniPFRelIso_all'%(lep_name, index), attribute = lep_getter("miniPFRelIso_all", index, abs_pdg),
#              binning=[50,0,.5],
#            ))
#            plots.append(Plot(
#              texX = 'pfRelIso03_all(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_pfRelIso03_all'%(lep_name, index), attribute = lep_getter("pfRelIso03_all", index, abs_pdg),
#              binning=[50,0,.5],
#            ))
#            plots.append(Plot(
#              texX = 'mvaTTH(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_mvaTTH'%(lep_name, index), attribute = lep_getter("mvaTTH", index, abs_pdg),
#              binning=[24,-1.2,1.2],
#            ))
#            plots.append(Plot(
#              texX = 'mvaTOP(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_mvaTOP'%(lep_name, index), attribute = lep_getter("mvaTOP", index, abs_pdg),
#              binning=[24,-1.2,1.2],
#            ))
#            plots.append(Plot(
#              texX = 'charge(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#              name = '%s%i_charge'%(lep_name, index), attribute = lep_getter("pdgId", index, abs_pdg, functor = charge),
#              binning=[3,-1,2],
#            ))
#            if lep_name == "mu":
#                plots.append(Plot(
#                  texX = 'segmentComp(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#                  name = '%s%i_segmentComp'%(lep_name, index), attribute = lep_getter("segmentComp", index, abs_pdg),
#                  binning=[50,0,1],
#                ))
#                plots.append(Plot(
#                  texX = 'nStations(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#                  name = '%s%i_nStations'%(lep_name, index), attribute = lep_getter("nStations", index, abs_pdg),
#                  binning=[10,0,10],
#                ))
#                plots.append(Plot(
#                  texX = 'nTrackerLayers(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
#                  name = '%s%i_nTrackerLayers'%(lep_name, index), attribute = lep_getter("nTrackerLayers", index, abs_pdg),
#                  binning=[20,0,20],
#                ))
#            if lep_name == "ele":
#                for cbIdFlag in vidNestedWPBitMapNamingList:
#                    plots.append(Plot(
#                      texX = '%s(%s_{%i}) (GeV)'%(cbIdFlag, lep_name, index), texY = 'Number of Events',
#                      name = '%s%i_%s_Flag'%(lep_name, index, cbIdFlag), attribute = lep_getter("vidNestedWPBitmap", index, abs_pdg, functor = cbEleIdFlagGetter(cbIdFlag)),
#                      binning=[5,0,5],
#                    ))

plotting.fill(plots, read_variables = read_variables, sequence = sequence)

for plot in plots:
    for i_h, hl in enumerate(plot.histos):
        # dress up
        hl[0].legendText = params[i_h]['legendText']
        hl[0].style = params[i_h]['style']

drawPlots(plots)


logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
