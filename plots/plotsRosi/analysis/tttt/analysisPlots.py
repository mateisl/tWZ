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
c1.Print('/tmp/delete.png')

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
from tWZ.Tools.objectSelection           import lepString, isBJet
from tWZ.Tools.mt2Calculator             import mt2Calculator

from Analysis.Tools.helpers              import deltaPhi, deltaR, getCollection, getObjDict
from Analysis.Tools.puProfileCache       import *
from Analysis.Tools.puReweighting        import getReweightingFunction
from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons
import Analysis.Tools.syncer 

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--noData',         action='store_true', default=False, help='also plot data?')
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--dataMCScaling',  action='store_true', help='Data MC scaling?', )
argParser.add_argument('--plot_directory', action='store', default='tttt_test')
argParser.add_argument('--era',            action='store', type=str, default="Run2016")
argParser.add_argument('--selection',      action='store', default='trilepM-offZ1-minDLmass12-njet4p-btag2p')
#argParser.add_argument('--nanoAODv4',   default=True, action='store_true',                                                                        help="Run on nanoAODv4?" )
argParser.add_argument('--samples',        action='store',         nargs='*',  type=str, default=['TTZToLLNuNu_ext'],                  help="List of samples to be post-processed, given as CMG component name" )
#flagg for parton selection
argParser.add_argument('--normalize',      action='store_true', default=False,                                                                   help="Normalize to 1" )
argParser.add_argument('--scaled',         action='store_true', help='scaling', )
#argParser.add_argument('--mva',            action='store', type=str )

args = argParser.parse_args()
options = argParser.parse_args()

# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"
if args.noData:                       args.plot_directory += "_noData"
if args.normalize:                    args.plot_directory += "_normalize"
if args.scaled:                       args.plot_directory += "_scaled"
logger.info( "Working in era %s", args.era)

from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed import *

if args.era == "Run2016":
    mc = [Summer16.TTW, Summer16.TTZ, Summer16.TTTT , Summer16.nonprompt_3l ] #Summer16.TTX_rare, Summer16.TZQ, Summer16.WZ, Summer16.triBoson, Summer16.ZZ, Summer16.nonprompt_3l]
elif args.era == "Run2017":
    mc = [Fall17.TTW, Fall17.TTZ, Fall17.TTTT , Fall17.nonprompt_3l ] #, Fall17.TTX_rare, Fall17.TZQ, Fall17.WZ, Fall17.triBoson, Fall17.ZZ, Fall17.nonprompt_3l]
elif args.era == "Run2018":
    mc = [Autumn18.TTW, Autumn18.TTZ, Autumn18.TTTT , Autumn18.nonprompt_3l ] #, Autumn18.TTX_rare, Autumn18.TZQ, Autumn18.WZ, Autumn18.triBoson, Autumn18.ZZ, Autumn18.nonprompt_3l]
elif args.era == "RunII":
    #mc = [TTTT, TTZ, TTW, WZ ]
    mc = [TTW, TTZ, TTTT , nonprompt_3l, WZ ]
# data sample
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
    for sample in mc:# + [data_sample]:
        sample.normalization = 1.
        #sample.reduceFiles( factor = 40 )
        sample.reduceFiles( to=1)
        sample.scale /= sample.normalization

# Text on the plots
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)
tex.SetTextAlign(11) # align right

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
        
      if args.normalize:
        for comp in plot.histos:
            for h in comp:
                if h.Integral()!=0.: h.Scale(1./h.Integral())  
      
      if isinstance( plot, Plot):
            plotting.draw(plot,
            plot_directory = plot_directory_,
            ratio = {'yRange':(0.1,1.9)} if not args.noData else None,
            #ratio = {'yRange': (0.1, 1.9), 'histos':[(1,0)], 'texY':'Ratio'} if args.scaled else None,
            #ratio = {'yRange': (0.1, 1.9), 'histos':[(1,0),(2,0),(3,0),(4,0)], 'texY':'Ratio'} if args.scaled else None,
            logX = False, logY = log, sorting = False,
            yRange = (0.03, "auto") if log else (0.001, "auto"),
            scaling =  { 1:0 } if args.scaled else {}, 
            #scaling =  { i+1:0 for i in range(4) } if args.scaled else {}, 
            legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
            drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ) + _drawObjects,
            copyIndexPHP = True, extensions = ["png", "pdf", "root"],
          )

# Read variables and sequences

sequence       = []

def getmt2(event, sample): 
    mt2Calculator.reset()
    #set met 
    mt2Calculator.setMet(event.met_pt, event.met_phi )
    #set leptons
    mt2Calculator.setLeptons(event.l1_pt, event.l1_eta, event.l1_phi, event.l2_pt, event.l2_eta, event.l2_phi )

    setattr(event, "mt2ll", mt2Calculator.mt2ll() )

sequence.append(getmt2)

def getDeltaR(event, sample):
    event.jets     = [getObjDict(event, 'JetGood_', jetVarNames, i) for i in range(int(event.nJetGood))]
    bjets          = filter(lambda j:isBJet(j, year=event.year) and abs(j['eta'])<=2.4    , event.jets)
    event.minDRbjets = min( [ deltaR(b1, b2) for i, b1 in enumerate(bjets[:-1]) for b2 in bjets[i+1:]  ]) 
sequence.append(getDeltaR)

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


#jets without eta cut
jetVars          = ['pt/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'rawFactor/F', 'eta/F', 'phi/F', 'jetId/I', 'btagDeepB/F', 'btagDeepFlavB/F', 'btagCSVV2/F', 'area/F'] 
jetVarNames      = [x.split('/')[0] for x in jetVars]

lepVars         = ['pt/F','eta/F','phi/F','pdgId/I','muIndex/I','eleIndex/I','cutBased/I','miniPFRelIso_all/F','pfRelIso03_all/F','mvaFall17V2Iso_WP90/O', 'sip3d/F','lostHits/I','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','eleIndex/I','muIndex/I']
lepVarNames     = [x.split('/')[0] for x in lepVars]

bVars      = ['pt/F','eta/F','phi/F']
bVarNames  = [x.split('/')[0] for x in bVars]

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I", 
    "l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
    "l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
    "JetGood[pt/F,eta/F,phi/F]", 
    "nJet/I", 
    "nlep/I", "lep[%s]"%(",".join(lepVars)),
    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I", 
    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,pdgId/I,vidNestedWPBitmap/I]",
]

#b tagger
b_tagger = "DeepJet"

# read only for data:
read_variables_data = [ "Jet[%s]"%(",".join(jetVars))] 

# read only for MC:
read_variables_MC = ['reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F']
read_variables_MC.append( "Jet[%s,genJetIdx/I]"%(",".join(jetVars)) )
read_variables_MC.append( "GenJet[pt/F,eta/F,phi/F,hadronFlavour/b,partonFlavour/I]")

# define 3l selections

#DeltaEta
#def deltaEtaZll( event,sample ):
#    event.Z1_lldEta =  event.lep_eta[event.Z1_l2_index] - event.lep_eta[event.Z1_l1_index]
##    print event.Z1_lldEta
#
#sequence.append( deltaEtaZll )

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

#MVA

import TMB.MVA.configs as configs
config = configs.tttt_3l
read_variables += config.read_variables
sequence += config.sequence
# Add sequence that computes the MVA inputs
def make_mva_inputs( event, sample ):
    for mva_variable, func in config.mva_variables:
        setattr( event, mva_variable, func(event, sample) )
sequence.append( make_mva_inputs )


#load models

from keras.models import load_model

models = [
     #("FI_ctZ_BSM_TTG",      False, load_model("/mnt/hephy/cms/robert.schoefbeck/TMB/models/ctZ_BSM_TTG/ttG_WG/FI_ctZ_BSM/regression_model.h5")),
     #("tttt_ttw_ttz_nonprompt_LSTM",   True , load_model("/mnt/hephy/cms/rosmarie.schoefbeck/TMB/models/tttt_3l_ttw_ttz_nonprompt_LSTM/tttt_3l/regression_model.h5")),
     #("tttt_ttw_ttz_nonprompt", False, load_model("/mnt/hephy/cms/rosmarie.schoefbeck/TMB/models/tttt_3l_ttw_ttz_nonprompt/tttt_3l/regression_model.h5")),
     ("tttt_ttw_ttz_nonprompt_LSTM", True, load_model("/mnt/hephy/cms/rosmarie.schoefbeck/TMB/models/tttt_3l_ttw_ttz_nonprompt_v3_LSTM/tttt_3l/multiclass_model.h5")),
     ("tttt_ttw_ttz_nonprompt", False, load_model("/mnt/hephy/cms/rosmarie.schoefbeck/TMB/models/tttt_3l_ttw_ttz_nonprompt_v3/tttt_3l/multiclass_model.h5")),
     #("tttt_ttw_ttz_nonprompt_LSTM",   True , load_model("/mnt/hephy/cms/robert.schoefbeck/TMB/models/tttt_3l_ttw_ttz_nonprompt_v2_LSTM/tttt_3l/regression_model.h5")),
     #("tttt_ttw_ttz_nonprompt", False, load_model("/mnt/hephy/cms/robert.schoefbeck/TMB/models/tttt_3l_ttw_ttz_nonprompt_v2/tttt_3l/regression_model.h5")),
#     ("TTTT_Multiclass_LSTM",    True, load_model("/mnt/hephy/cms/rosmarie.schoefbeck/TMB/models/tttt_3l_test_LSTM/tttt_3l/regression_model.h5")),
#     ("Multiclass_TTTT",  False, load_model("/mnt/hephy/cms/rosmarie.schoefbeck/TMB/models/tttt_3l_test/tttt_3l/regression_model.h5")),
]

def keras_predict( event, sample ):

    # get model inputs assuming lstm
    flat_variables, lstm_jets = config.predict_inputs( event, sample, jet_lstm = True)
    for name, has_lstm, model in models:
        #print has_lstm, flat_variables, lstm_jets
        prediction = model.predict( flat_variables if not has_lstm else [flat_variables, lstm_jets] )
        for i_pred, value in enumerate(prediction[0]):  
           setattr( event, name+'_'+config.training_samples[i_pred].name, value )

sequence.append( keras_predict )


weight_ = lambda event, sample: event.weight if sample.isData else event.weight*lumi_year[event.year]/1000.


yields     = {}
allPlots   = {}
allModes   = ['mumumu','mumue','muee', 'eee']
for i_mode, mode in enumerate(allModes):
    yields[mode] = {}
    if not args.noData:
        data_sample.texName = "data"
        data_sample.setSelectionString([getLeptonSelection(mode)])
        data_sample.name           = "data"
        data_sample.style          = styles.errorStyle(ROOT.kBlack)
        data_sample.read_variables = read_variables_data 
        lumi_scale                 = data_sample.lumi/1000

    for sample in mc: 
        sample.style = styles.fillStyle(sample.color)
        #sample.style = styles.lineStyle(sample.color)
        sample.read_variables = read_variables_MC 
        sample.setSelectionString([getLeptonSelection(mode)])
        sample.weight = lambda event, sample: event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger#*event.reweightLeptonSF

    
    if args.scaled: 
        stack = Stack()
    else:
        if not args.noData:
          stack = Stack(mc, data_sample)
        else:
          stack = Stack(mc)
    
#    for sample in mc+[data_sample]:
#        sample.weight = weight_
#
    # Use some defaults
    Plot.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = cutInterpreter.cutString(args.selection))

    plots = []

    for name, has_lstm, model in models:
        for i_tr_s, tr_s in enumerate( config.training_samples ):
            disc_name = name+'_'+config.training_samples[i_tr_s].name
            plots.append(Plot(
                texX = disc_name, texY = 'Number of Events',
                name = disc_name, attribute = lambda event, sample, disc_name=disc_name: getattr( event, disc_name ),
                binning=[50, 0, 1],
            ))

            plots.append(Plot(
                texX = disc_name+"coarse", texY = 'Number of Events',
                name = disc_name+"coarse", attribute = lambda event, sample, disc_name=disc_name: getattr( event, disc_name ),
                binning=[18, 0, 1],
            ))

            plots.append(Plot(
                texX = disc_name+"medcoarse", texY = 'Number of Events',
                name = disc_name+"medcoarse", attribute = lambda event, sample, disc_name=disc_name: getattr( event, disc_name ),
                binning=[32, 0, 1],
            ))

    plots.append(Plot(
        texX = "mt2ll" , texY = 'Number of Events',
        name = "mt2ll", attribute = lambda event, sample: getattr( event, "mt2ll" ),
        binning=[20, 0, 300],
    ))

    plots.append(Plot(
        texX = "minDR_bjets" , 
        texY = 'Number of Events',
        name = "minDR_bjets_coarse", 
        attribute = lambda event, sample: event.minDRbjets ,
        binning=[10, 0, 5],
    ))

    plots.append(Plot(
        texX = "minDR_bjets" , 
        texY = 'Number of Events',
        name = "minDR_bjets", 
        attribute = lambda event, sample: event.minDRbjets ,
        binning=[20, 0, 5],
    ))

#    plots.append(Plot(
#        texX = 'Multiclass_TTTT_p0', 
#        texY = 'Number of Events',
#        name = 'multiclass_tttt_p0',
#        attribute = lambda event, sample: event.Multiclass_TTTT_p0,
#        binning=[50,0,1],
#    ))
#
#    plots.append(Plot(
#        texX = 'Multiclass_TTTT_p1', 
#        texY = 'Number of Events',
#        name = 'multiclass_tttt_p1', attribute = lambda event, sample: event.Multiclass_TTTT_p1,
#        binning=[50,0,1],
#    ))
#
#    plots.append(Plot(
#        texX = 'Multiclass_TTTT_p2', 
#        texY = 'Number of Events',
#        name = 'multiclass_tttt_p2', attribute = lambda event, sample: event.Multiclass_TTTT_p2,
#        binning=[50,0,1],
#    ))
#    plots.append(Plot(
#      name = 'maxpt',
#      texX = 'maxpt',
#      texY = 'Number of Events',
#      attribute = lambda event, sample: event.max_pt_jet, 
#      binning=[20, 0, 400],
#    ))
#
#    plots.append(Plot(
#      name = 'maxabseta',
#      texX = 'abs(#eta)_max',
#      texY = 'Number of Events',
#      attribute = lambda event, sample: event.maxEta_of_pt30jets,
#      binning=[20, 0, 5],
#    ))
#
#    plots.append(Plot(
#      name = 'partons in jets mineta',
#      texX = 'partons in jets (mineta)',
#      texY = 'Number of Events',
#      attribute = lambda event, sample: event.partonsinfwdjetsmineta,
#      binning=[28, -6, 22],
#    ))
#
#    plots.append(Plot(
#      name = 'partons in fwd jets maxeta',
#      texX = 'partons in fwd jets (maxeta)',
#      texY = 'Number of Events',
#      attribute = lambda event, sample: event.partonsinfwdjets,
#      binning=[28, -6, 22],
#    ))

    plots.append(Plot(
      name = 'yield', texX = '', texY = 'Number of Events',
      attribute = lambda event, sample: 0.5 + i_mode,
      binning=[4, 0, 4],
    ))

#    plots.append(Plot(
#      texX = '#Delta R(maxetajet, Z_{1})', texY = 'Number of Events',
#      name = 'maxptjet_Z1_deltaR', attribute = lambda event, sample: event.maxeta_Z_deltaR,
#      binning=[20,0,6],
#    ))
#
#    plots.append(Plot(
#        name = "M2l",
#        texX = 'M(2l) (GeV)', texY = 'Number of Events',
#        attribute = lambda event, sample:event.M2l,
#        binning=[25,0,500],
#    ))
#
#    plots.append(Plot(
#        name = "m2l",
#        texX = 'M(2l) (GeV)', texY = 'Number of Events',
#        attribute = lambda event, sample:event.M2l,
#        binning=[25,50,150],
#    ))
#

#    plots.append(Plot(
#        name = "Deltaphi_l12",
#        texX = '#Delta#phi(l_12)', texY = 'Number of Events',
#        attribute = lambda event, sample: acos(cos(event.l1_phi-event.l2_phi)),
#        binning=[10,0,pi],
#    ))
#
#    plots.append(Plot(
#        name = "Deltaphi_j12",
#        texX = '#Delta#phi(j_12)', texY = 'Number of Events',
#        attribute = lambda event, sample: acos(cos(event.JetGood_phi[0]-event.JetGood_phi[1])) if event.nJetGood >=2 else 0,
#        binning=[10,0,pi],
#    ))

    plots.append(Plot(
        name = 'Deltaeta_l12',
        texX = '#Delta#eta(l12)', texY = 'Number of Events',
        attribute = lambda event, sample: event.l1_eta-event.l2_eta,
        binning=[10,-3,3],
    ))

    plots.append(Plot(
        name = 'Deltaeta_j12',
        texX = '#Delta#eta(j12)', texY = 'Number of Events',
        attribute = lambda event, sample: event.JetGood_eta[0] - event.JetGood_eta[1] if event.nJetGood >=2 else -10,
        binning=[10,-3,3],
    ))

     
    plots.append(Plot(
        name = "ht",
        texX = 'ht (GeV)', texY = 'Number of Events ',
        attribute = lambda event, sample: sum( [j['pt'] for j in event.jets] ),
        binning=[20,0,200],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "htb",
        texX = 'htb (GeV)', texY = 'Number of Events ',
        attribute = lambda event, sample: sum( [j['pt'] for j in event.bJets] ),
        binning=[20,0,400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "ht_ratio",
        texX = 'ht_ratio (GeV)', texY = 'Number of Events ',
        attribute = lambda event, sample: sum( [j['pt'] for j in event.jets[:4]])/ sum( [j['pt'] for j in event.jets ]) if event.nJetGood>=4 else 1 ,
        binning=[20,0,2],
        addOverFlowBin='upper',
    ))
     
    plots.append(Plot(
        name = "mT_l1",
        texX = 'mT_l1 (GeV)', texY = 'Number of Events ',
        attribute = lambda event, sample:sqrt(2*event.l1_pt*event.met_pt*(1-cos(event.l1_phi-event.met_phi))),
        binning=[20,0,300],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "mT_l2",
        texX = 'mT_l2 (GeV)', texY = 'Number of Events ',
        attribute = lambda event, sample:sqrt(2*event.l2_pt*event.met_pt*(1-cos(event.l2_phi-event.met_phi))),
        binning=[20,0,300],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "mT_l3",
        texX = 'mT_l3 (GeV)', texY = 'Number of Events ',
        attribute = lambda event, sample:sqrt(2*event.l3_pt*event.met_pt*(1-cos(event.l3_phi-event.met_phi))),
        binning=[20,0,300],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "mlj_11",
        texX = 'mlj_11 (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: sqrt(event.l1_pt*event.jets[0]['pt']*cosh(event.l1_eta-event.jets[0]['eta'])-cos(event.l1_phi-event.jets[0]['phi'])) if event.nJetGood >=1 else 0 ,
        binning=[25,0,500],
    ))

    plots.append(Plot(
        name = "mlj_12",
        texX = 'mlj_12 (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: sqrt(event.l1_pt*event.jets[1]['pt']*cosh(event.l1_eta-event.jets[1]['eta'])-cos(event.l1_phi-event.jets[1]['phi'])) if event.nJetGood >=2 else 0 ,
        binning=[25,0,500],
    ))

    plots.append(Plot(
        name = "ml_12",
        texX = 'ml_12 (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: sqrt(2*event.l1_pt*event.l2_pt*(cosh(event.l1_eta-event.l2_eta)-cos(event.l1_phi-event.l2_phi))) ,
        binning=[25,0,500],
    ))

    plots.append(Plot(
        name = "mj_12",
        texX = 'mj_12 (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: sqrt(event.jets[0]['pt']*event.jets[1]['pt']*cosh(event.jets[0]['eta']-event.jets[1]['eta'])-cos(event.jets[0]['phi']-event.jets[1]['phi']))  if event.nJetGood >=2 else 0 ,
        binning=[25,0,500],
    ))

    plots.append(Plot(
      name = 'nVtxs', texX = 'vertex multiplicity', texY = 'Number of Events',
      attribute = TreeVariable.fromString( "PV_npvsGood/I" ),
      binning=[50,0,50],
      addOverFlowBin='upper',
    ))

    plots.append(Plot(
        texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "met_pt/F" ),
        binning=[400/20,0,400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        texX = '#phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "met_phi/F" ),
        binning=[10,-pi,pi],
    ))

#    plots.append(Plot(
#        name = 'Z1_pt_coarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events / 50 GeV',
#        attribute = TreeVariable.fromString( "Z1_pt/F" ),
#        binning=[16,0,800],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'Z1_pt_superCoarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
#        attribute = TreeVariable.fromString( "Z1_pt/F" ),
#        binning=[3,0,600],
#        addOverFlowBin='upper',
#    ))
#
#    plots.append(Plot(
#        name = 'Z1_pt_coarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events / 50 GeV',
#        attribute = TreeVariable.fromString( "Z1_pt/F" ),
#        binning=[16,0,800],
#    ))
#
#    plots.append(Plot(
#        name = 'Z1_pt_superCoarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
#        attribute = TreeVariable.fromString( "Z1_pt/F" ),
#        binning=[3,0,600],
#    ))

    plots.append(Plot(
        name = "M3l",
        texX = 'M(3l) (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample:event.M3l,
        binning=[25,0,500],
    ))

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

    plots.append(Plot(
      texX = 'p_{T}(leading l) (GeV)', texY = 'Number of Events / 20 GeV',
      name = 'lep1_pt', attribute = lambda event, sample: event.lep_pt[0],
      binning=[400/20,0,400],
    ))

    plots.append(Plot(
      texX = 'p_{T}(subleading l) (GeV)', texY = 'Number of Events / 10 GeV',
      name = 'lep2_pt', attribute = lambda event, sample: event.lep_pt[1],
      binning=[200/10,0,200],
    ))

    plots.append(Plot(
      texX = 'p_{T}(trailing l) (GeV)', texY = 'Number of Events / 10 GeV',
      name = 'lep3_pt', attribute = lambda event, sample: event.lep_pt[2],
      binning=[150/10,0,150],
    ))

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
    plots.append(Plot(
        texX = 'M(ll) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "Z1_mass/F" ),
        binning=[10,81,101],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "Z1_mass_wide",
        texX = 'M(ll) (GeV)', texY = 'Number of Events / 2 GeV',
        attribute = TreeVariable.fromString( "Z1_mass/F" ),
        binning=[50,20,120],
        addOverFlowBin='upper',
    )) 

#    plots.append(Plot(
#        name = "Z1_cosThetaStar", texX = 'cos#theta(l-)', texY = 'Number of Events / 0.2',
#        attribute = lambda event, sample:event.Z1_cosThetaStar,
#        binning=[10,-1,1],
#    ))

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
#        name = 'Deltaeta_Z1_ll',
#        texX = '#Delta#eta(Z_{1}(ll))', texY = 'Number of Events',
#        attribute = lambda event, sample: event.Z1_lldEta,
#        binning=[10,-3,3],
#    ))

    plots.append(Plot(
      texX = 'N_{jets}', texY = 'Number of Events',
      attribute = TreeVariable.fromString( "nJetGood/I" ), #nJetSelected
      binning=[8,-0.5,13.5],
    ))

    plots.append(Plot(
      texX = 'N_{b-tag}', texY = 'Number of Events',
      attribute = TreeVariable.fromString( "nBTag/I" ), #nJetSelected
      binning=[4,-0.5,10.5],
    ))

    plots.append(Plot(
      texX = 'N_{lep}', texY = 'Number of Events',
      attribute = TreeVariable.fromString( "nlep/I" ), #nJetSelected
      binning=[4,-0.5,10.5],
    ))


    plots.append(Plot(
      texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events ',
      name = 'jet0_pt', attribute = lambda event, sample: event.JetGood_pt[0],
      binning=[10,0,300],
    ))

    plots.append(Plot(
      texX = 'p_{T}(3) (GeV)', texY = 'Number of Events ',
      name = 'jet3_pt', attribute = lambda event, sample: event.JetGood_pt[2] if event.nJetGood >=4 else 0,
      binning=[10,0,300],
    ))

    plots.append(Plot(
      texX = 'p_{T}(4) (GeV)', texY = 'Number of Events',
      name = 'jet4_pt', attribute = lambda event, sample: event.JetGood_pt[2] if event.nJetGood >=5 else 0,
      binning=[10,0,300],
    ))

    plots.append(Plot(
      texX = 'p_{T}(5) (GeV)', texY = 'Number of Events',
      name = 'jet5_pt', attribute = lambda event, sample: event.JetGood_pt[2] if event.nJetGood >=6 else 0,
      binning=[10,0,300],
    ))

    plots.append(Plot(
      texX = 'p_{T}(6) (GeV)', texY = 'Number of Events',
      name = 'jet6_pt', attribute = lambda event, sample: event.JetGood_pt[2] if event.nJetGood >=7 else 0,
      binning=[10,0,300],
    ))

    plots.append(Plot(
      texX = 'p_{T}(7) (GeV)', texY = 'Number of Events',
      name = 'jet7_pt', attribute = lambda event, sample: event.JetGood_pt[2] if event.nJetGood >=8 else 0,
      binning=[10,0,300],
    ))

    plots.append(Plot(
      texX = 'p_{T}(subleading jet) (GeV)', texY = 'Number of Events ',
      name = 'jet1_pt', attribute = lambda event, sample: event.JetGood_pt[1],
      binning=[10,0,400],
    ))

    plots.append(Plot(
      texX = 'p_{T}(2) (GeV)', texY = 'Number of Events ',
      name = 'jet2_pt', attribute = lambda event, sample: event.JetGood_pt[2],
      binning=[10,0,300],
    ))

    plots.append(Plot(
        name = "W_pt",
        texX = 'p_{T}(W) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = lambda event, sample:event.W_pt,
        binning=[20,0,400],
    ))

 #bjets 

#    plots.append(Plot(
#      texX = 'p_{T}(b jet) (GeV)', texY = 'Number of Events / 30 GeV',
#      name = 'bjet_pt', attribute = lambda event, sample: event.bjet_pt,
#      binning=[600/30,0,600],
#    ))
#
#    
#    plots.append(Plot(
#      texX = '#Delta R(bjet, Z_{1})', texY = 'Number of Events',
#      name = 'bjet_Z1_deltaR', attribute = lambda event, sample: event.bjet_Z1_deltaR,
#      binning=[20,0,6],
#    ))
#
#    plots.append(Plot(
#      texX = '#Delta R(bjet, nonZ-l_{1})', texY = 'Number of Events',
#      name = 'bjet_nonZ1l1_deltaR', attribute = lambda event, sample: event.bjet_nonZ1l1_deltaR,
#      binning=[20,0,6],
#    ))
#
#    plots.append(Plot(
#      name = 'btagDeepB',
#      texX = 'btagDeepB',
#      texY = 'Number of Events',
#      attribute = lambda event, sample: event.bjet_maxbtagdeepb,
#      binning=[20, 0, 5],
#    ))

   # 3l training variables


#    plots.append(Plot(
#      texX = '#Delta#eta(nonZ-l_{1}, Z_{1})', texY = 'Number of Events',
#      name = 'nonZ1l1_Z1_deltaEta', attribute = lambda event, sample: event.nonZ1l1_Z1_deltaEta,
#      binning=[20,0,6],
#    ))
#
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
#
#    plots.append(Plot(
#      texX = '#Delta R(jet_{0}, nonZ-l_{1})', texY = 'Number of Events',
#      name = 'jet0_nonZ1l1_deltaR', attribute = lambda event, sample: event.jet0_nonZ1l1_deltaR,
#      binning=[20,0,6],
#    ))
#
#    plots.append(Plot(
#      texX = '#Delta R(jet_{1}, Z_{1})', texY = 'Number of Events',
#      name = 'jet1_Z1_deltaR', attribute = lambda event, sample: event.jet1_Z1_deltaR,
#      binning=[20,0,6],
#    ))
#
#    plots.append(Plot(
#      texX = '#Delta R(jet_{1}, nonZ-l_{1})', texY = 'Number of Events',
#      name = 'jet1_nonZ1l1', attribute = lambda event, sample: event.jet1_nonZ1l1_deltaR,
#      binning=[20,0,6],
#    ))
#    
#newplots 

    #plots.append(Plot(
    #   texX = 'mvaTOP(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
    #   name = '%s%i_mvaTOP'%(lep_name, index), attribute = lep_getter("mvaTOP", index, abs_pdg),
    #   binning=[24,-1.2,1.2],
    # ))

    plots.append(Plot(
       name = 'l1_pt',
       texX = 'p_{T}(l_{1}) (GeV)', texY = 'Number of Events GeV',
       attribute = lambda event, sample:event.l1_pt,
       binning=[15,0,300],
   ))

    plots.append(Plot(
       name = 'l1_eta',
       texX = '#eta(l_{1})', texY = 'Number of Events',
       attribute = lambda event, sample: event.l1_eta,
       binning=[20,-3,3],
   ))

    plots.append(Plot(
       name = 'l2_eta',
       texX = '#eta(l_{2})', texY = 'Number of Events',
       attribute = lambda event, sample: event.l2_eta,
       binning=[20,-3,3],
   ))

    plots.append(Plot(
       name = 'l3_eta',
       texX = '#eta(l_{3})', texY = 'Number of Events',
       attribute = lambda event, sample: event.l3_eta,
       binning=[20,-3,3],
   ))

    plots.append(Plot(
       name = 'jet0_eta',
       texX = '#eta(jet_{0})', texY = 'Number of Events',
       attribute = lambda event, sample: event.JetGood_eta[0] if event.nJetGood >=1 else -10,
       binning=[20,-3,3],
   ))

    plots.append(Plot(
       name = 'jet1_eta',
       texX = '#eta(jet_{1})', texY = 'Number of Events',
       attribute = lambda event, sample: event.JetGood_eta[1] if event.nJetGood >=2 else -10,
       binning=[20,-3,3],
   ))

    plots.append(Plot(
       name = 'jet2_eta',
       texX = '#eta(jet_{2})', texY = 'Number of Events',
       attribute = lambda event, sample: event.JetGood_eta[2] if event.nJetGood >=3 else -10,
       binning=[20,-3,3],
   ))

    plots.append(Plot(
      name = 'jet0_btagDeepB',
      texX = 'jet0_btagDeepB',
      texY = 'Number of Events',
      attribute = lambda event, sample: event.JetGood_btagDeepB[0]   if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10,
      binning=[20, 0, 1],
    ))

    plots.append(Plot(
      name = 'jet1_btagDeepB',
      texX = 'jet1_btagDeepB',
      texY = 'Number of Events',
      attribute = lambda event, sample: event.JetGood_btagDeepB[1]   if (event.nJetGood >=2 and event.JetGood_btagDeepB[1]>-10) else -10,
      binning=[20, 0, 1],
    ))

    plots.append(Plot(
      name = 'jet2_btagDeepB',
      texX = 'jet2_btagDeepB',
      texY = 'Number of Events',
      attribute = lambda event, sample: event.JetGood_btagDeepB[2]   if (event.nJetGood >=3 and event.JetGood_btagDeepB[2]>-10) else -10,
      binning=[20, 0, 1],
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
       texX = 'p_{T}(l_{2}) (GeV)', texY = 'Number of Events',
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
       texX = 'MVA_{TOP}(l_{2}) WP', texY = 'Number of Events',
       attribute = lambda event, sample: event.l2_mvaTOPWP,
       binning=[5,0,5],
   ))

    plots.append(Plot(
       name = 'l3_pt',
       texX = 'p_{T}(l_{3}) (GeV)', texY = 'Number of Events',
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
       texX = 'MVA_{TOP}(l_{3}) WP', texY = 'Number of Events',
       attribute = lambda event, sample: event.l3_mvaTOPWP,
       binning=[5,0,5],
   ))

    plotting.fill(plots, read_variables = read_variables, sequence = sequence)

    # Get normalization yields from yield histogram
    for plot in plots:
      if plot.name == "yield":
        for i, l in enumerate(plot.histos):
          for j, h in enumerate(l):
            yields[mode][plot.stack[i][j].name] = h.GetBinContent(h.FindBin(0.5+i_mode))
            h.GetXaxis().SetBinLabel(1, "#mu#mu#mu")
            h.GetXaxis().SetBinLabel(2, "#mu#mue")
            h.GetXaxis().SetBinLabel(3, "#muee")
            h.GetXaxis().SetBinLabel(4, "eee")
      if plot.name.endswith("_Flag"):
        for i, l in enumerate(plot.histos):
          for j, h in enumerate(l):
            h.GetXaxis().SetBinLabel(1, "fail")
            h.GetXaxis().SetBinLabel(2, "veto")
            h.GetXaxis().SetBinLabel(3, "loose")
            h.GetXaxis().SetBinLabel(4, "medium")
            h.GetXaxis().SetBinLabel(5, "tight")
        
    if args.noData: yields[mode]["data"] = 0

    yields[mode]["MC"] = sum(yields[mode][s.name] for s in mc)
    dataMCScale        = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')

    drawPlots(plots, mode, dataMCScale)
    allPlots[mode] = plots


# Add the different channels into SF and all
for mode in ["comb1","comb2","all"]:
    yields[mode] = {}
    for y in yields[allModes[0]]:
        try:    yields[mode][y] = sum(yields[c][y] for c in ['eee','muee','mumue', 'mumumu'])
        except: yields[mode][y] = 0
    dataMCScale = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')
    
    for plot in allPlots['mumumu']:
        if mode=="comb1":
            tmp = allPlots['mumue']
        elif mode=="comb2":
            tmp = allPlots['muee']
        else:
            tmp = allPlots['eee']
        for plot2 in (p for p in tmp if p.name == plot.name):
            for i, j in enumerate(list(itertools.chain.from_iterable(plot.histos))):
                for k, l in enumerate(list(itertools.chain.from_iterable(plot2.histos))):
                    if i==k:
                        j.Add(l)
    
    if mode == "all": drawPlots(allPlots['mumumu'], mode, dataMCScale)

import pickle
pickle.dump( {p.name: p.histos for p in allPlots['mumumu'] if "TTTT" in p.name}, file( os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era,  args.selection+'.pkl'), 'w' ))

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
