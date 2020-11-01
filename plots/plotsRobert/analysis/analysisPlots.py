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
import Analysis.Tools.syncer
import numpy as np

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--noData',         action='store_true', default=False, help='also plot data?')
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--dataMCScaling',  action='store_true', help='Data MC scaling?', )
argParser.add_argument('--plot_directory', action='store', default='tWZ_v3')
argParser.add_argument('--era',            action='store', type=str, default="Run2016")
argParser.add_argument('--selection',      action='store', default='trilepT-minDLmass12-onZ1-njet4p-btag1')
args = argParser.parse_args()

# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"
if args.noData:                       args.plot_directory += "_noData"

logger.info( "Working in era %s", args.era)

#tWZ_sample = TWZ if args.nominalSignal else yt_TWZ_filter
from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed import *

if args.era == "Run2016":
    mc = [Summer16.TWZ_NLO_DR, Summer16.TTZ, Summer16.TTX_rare, Summer16.TZQ, Summer16.WZ, Summer16.triBoson, Summer16.ZZ, Summer16.nonprompt_3l]
elif args.era == "Run2017":
    mc = [Fall17.TWZ_NLO_DR, Fall17.TTZ, Fall17.TTX_rare, Fall17.TZQ, Fall17.WZ, Fall17.triBoson, Fall17.ZZ, Fall17.nonprompt_3l]
elif args.era == "Run2018":
    mc = [Autumn18.TWZ_NLO_DR, Autumn18.TTZ, Autumn18.TTX_rare, Autumn18.TZQ, Autumn18.WZ, Autumn18.triBoson, Autumn18.ZZ, Autumn18.nonprompt_3l]
elif args.era == "RunII":
    mc = [TWZ_NLO_DR, TTZ, TTX_rare, TZQ, WZ, triBoson, ZZ, nonprompt_3l]

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
    for sample in mc + [data_sample]:
        sample.normalization = 1.
        sample.reduceFiles( to = 1 )
        #sample.reduceFiles( to=1)
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

      if isinstance( plot, Plot):
          plotting.draw(plot,
            plot_directory = plot_directory_,
            ratio = {'yRange':(0.1,1.9)} if not args.noData else None,
            logX = False, logY = log, sorting = True,
            yRange = (0.03, "auto") if log else (0.001, "auto"),
            scaling = {0:1} if args.dataMCScaling else {},
            legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
            drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ) + _drawObjects,
            copyIndexPHP = True, extensions = ["png"],
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
    "l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
    "JetGood[pt/F,eta/F,phi/F]",
    "lep[pt/F,eta/F,phi/F,pdgId/I,muIndex/I,eleIndex/I]",
    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I", 
    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
]

read_variables_MC = ['reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F']
# define 3l selections

#MVA
#from Analysis.TMVA.Reader    import Reader
from tWZ.MVA.MVA_TWZ_3l_LTS  import mva_variables#, mlp1, mlp_tanh
from tWZ.MVA.MVA_TWZ_3l_LTS  import sequence as mva_sequence
from tWZ.MVA.MVA_TWZ_3l_LTS  import read_variables as mva_read_variables
from tWZ.Tools.user          import mva_directory

sequence.extend( mva_sequence )
read_variables.extend( mva_read_variables )

#mva_reader = Reader(
#    mva_variables     = mva_variables,
#    weight_directory  = os.path.join( mva_directory, "TWZ_3l" ),
#    label             = "TWZ_3l")
#
#def makeDiscriminator( mva ):
#    def _getDiscriminator( event, sample ):
#        kwargs = {name:func(event,None) for name, func in mva_variables.iteritems()}
#        setattr( event, mva['name'], mva_reader.evaluate(mva['name'], **kwargs))
#    return _getDiscriminator
#
def discriminator_getter(name):
    def _disc_getter( event, sample ):
        return getattr( event, name )
    return _disc_getter

#mvas = [mlp_tanh]
#for mva in mvas:
#    mva_reader.addMethod(method=mva)
#    sequence.append( makeDiscriminator(mva) )

from ML.models.tWZ.tWZ_multiclass import variables as keras_varnames
from ML.models.tWZ.tWZ_multiclass import output_specification as keras_output_specification 
from ML.models.tWZ.tWZ_multiclass import model     as keras_multiclass

def keras_predict( event, sample ):
    #print np.array([[mva_variables[varname](event, sample) for varname in keras_varnames]])
    event.keras_multiclass_prediction = keras_multiclass.predict( np.array([[mva_variables[varname](event, sample) for varname in keras_varnames]])) 
    event.keras_multiclass_prediction = event.keras_multiclass_prediction[0]
    #for i in xrange( len(keras_output_specification) ):
    #    print "keras_multiclass_"+keras_output_specification[i]
    #    print  event.keras_multiclass_prediction, i
    #    print  event.keras_multiclass_prediction[i]
    for i in xrange( len(keras_output_specification) ):
        setattr( event, "keras_multiclass_"+keras_output_specification[i], event.keras_multiclass_prediction[i] ) 
    #print event.keras_multiclass_prediction
sequence.append( keras_predict )

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

#mu0_charge   = lep_getter("pdgId", 0, 13, functor = charge)
#ele0_charge = lep_getter("pdgId", 0, 11, functor = charge)
#mu1_charge   = lep_getter("pdgId", 1, 13, functor = charge)
#ele1_charge = lep_getter("pdgId", 1, 11, functor = charge)
#def test(event, sample):
#    mu0_ch  = mu0_charge(event, sample)
#    ele0_ch = ele0_charge(event, sample)
#    mu1_ch  = mu1_charge(event, sample)
#    ele1_ch = ele1_charge(event, sample)
#    print "mu0_ch",mu0_ch, "ele0_ch",ele0_ch, "mu1_ch",mu1_ch, "ele1_ch",ele1_ch
#
#sequence.append( test )

# 3l trainign variables

def make_training_observables_3l(event, sample):

    event.nonZ1l1_Z1_deltaPhi = deltaPhi(event.lep_phi[event.nonZ1_l1_index], event.Z1_phi)
    event.nonZ1l1_Z1_deltaEta = abs(event.lep_eta[event.nonZ1_l1_index] - event.Z1_eta)
    event.nonZ1l1_Z1_deltaR   = deltaR({'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet0_Z1_deltaR      = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet0_nonZ1l1_deltaR = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    event.jet1_Z1_deltaR      = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet1_nonZ1l1_deltaR = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})

sequence.append( make_training_observables_3l )

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
        lumi_scale                 = data_sample.lumi/1000

    weight_ = lambda event, sample: event.weight if sample.isData else event.weight*lumi_year[event.year]/1000.

    for sample in mc: sample.style = styles.fillStyle(sample.color)
    
    for sample in mc:
      sample.read_variables = read_variables_MC 
      sample.setSelectionString([getLeptonSelection(mode)])
      sample.weight = lambda event, sample: event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger#*event.reweightLeptonSF

    #yt_TWZ_filter.scale = lumi_scale * 1.07314

    if not args.noData:
      stack = Stack(mc, data_sample)
    else:
      stack = Stack(mc)

    # Use some defaults
    Plot.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = cutInterpreter.cutString(args.selection))

    plots = []

    plots.append(Plot(
      name = 'yield', texX = '', texY = 'Number of Events',
      attribute = lambda event, sample: 0.5 + i_mode,
      binning=[4, 0, 4],
    ))

    for output in keras_output_specification:
        plots.append(Plot(
            texX = 'keras multiclass '+output, texY = 'Number of Events',
            name = 'keras_multiclass_'+output, attribute = discriminator_getter('keras_multiclass_'+output),
            binning=[50, 0, 1],
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
        texX = '#phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "met_phi/F" ),
        binning=[10,-pi,pi],
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

    plots.append(Plot(
        name = "M3l",
        texX = 'M(3l) (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample:event.M3l,
        binning=[25,0,500],
    ))

    plots.append(Plot(
        name = "dPhiZJet",
        texX = '#Delta#phi(Z,j1)', texY = 'Number of Events',
        attribute = lambda event, sample: deltaPhi(event.Z1_phi, event.JetGood_phi[0]),
        binning=[20,0,pi],
    ))

    plots.append(Plot(
        name = "l1_Z1_pt",
        texX = 'p_{T}(l_{1,Z}) (GeV)', texY = 'Number of Events / 10 GeV',
        attribute = lambda event, sample:event.lep_pt[event.Z1_l1_index],
        binning=[30,0,300],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "l1_Z1_pt_coarse",
        texX = 'p_{T}(l_{1,Z}) (GeV)', texY = 'Number of Events / 40 GeV',
        attribute = lambda event, sample:event.lep_pt[event.Z1_l1_index],
        binning=[10,0,400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'l1_Z1_pt_ext', texX = 'p_{T}(l_{1,Z}) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = lambda event, sample:event.lep_pt[event.Z1_l1_index],
        binning=[20,40,440],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = "l2_Z1_pt",
        texX = 'p_{T}(l_{2,Z}) (GeV)', texY = 'Number of Events / 10 GeV',
        attribute = lambda event, sample:event.lep_pt[event.Z1_l2_index],
        binning=[20,0,200],
        addOverFlowBin='upper',
    ))

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

    plots.append(Plot(
        name = "l2_Z1_pt_coarse",
        texX = 'p_{T}(l_{2,Z}) (GeV)', texY = 'Number of Events / 10 GeV',
        attribute = lambda event, sample:event.lep_pt[event.Z1_l2_index],
        binning=[10,0,200],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'l2_Z1_pt_ext', texX = 'p_{T}(l_{2,Z}) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = lambda event, sample:event.lep_pt[event.Z1_l2_index],
        binning=[20,0,400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'lnonZ1_pt',
        texX = 'p_{T}(l_{1,extra}) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = lambda event, sample:event.lep_pt[event.nonZ1_l1_index],
        binning=[15,0,300],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'lnonZ1_pt_coarse',
        texX = 'p_{T}(l_{1,extra}) (GeV)', texY = 'Number of Events / 60 GeV',
        attribute = lambda event, sample:event.lep_pt[event.nonZ1_l1_index],
        binning=[3,0,180],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'lnonZ1_charge',
        texX = 'Charge(l_{1,extra})', texY = 'Number of Events',
        attribute = lambda event, sample:-event.lep_pdgId[event.nonZ1_l1_index]/abs(event.lep_pdgId[event.nonZ1_l1_index]),
        binning=[2,-1,1],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'lnonZ1_eta',
        texX = '#eta(l_{1,extra})', texY = 'Number of Events',
        attribute = lambda event, sample: event.lep_eta[event.nonZ1_l1_index],
        binning=[20,-3,3],
    ))

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

    plots.append(Plot(
        name = "Z1_cosThetaStar", texX = 'cos#theta(l-)', texY = 'Number of Events / 0.2',
        attribute = lambda event, sample:event.Z1_cosThetaStar,
        binning=[10,-1,1],
    ))

    plots.append(Plot(
        name = "Z2_mass_wide",
        texX = 'M(ll) of 2nd OSDL pair', texY = 'Number of Events / 2 GeV',
        attribute = TreeVariable.fromString( "Z2_mass/F" ),
        binning=[60,0,120],
        addOverFlowBin='upper',
    )) 

    plots.append(Plot(
        name = "minDLmass",
        texX = 'min mass of all DL pairs', texY = 'Number of Events / 2 GeV',
        attribute = TreeVariable.fromString( "minDLmass/F" ),
        binning=[60,0,120],
        addOverFlowBin='upper',
    )) 

    plots.append(Plot(
        texX = '#Delta#phi(Z_{1}(ll))', texY = 'Number of Events',
        attribute = TreeVariable.fromString( "Z1_lldPhi/F" ),
        binning=[10,0,pi],
    ))

    plots.append(Plot(
        texX = '#Delta R(Z_{1}(ll))', texY = 'Number of Events',
        attribute = TreeVariable.fromString( "Z1_lldR/F" ),
        binning=[10,0,6],
    ))

    plots.append(Plot(
      texX = 'N_{jets}', texY = 'Number of Events',
      attribute = TreeVariable.fromString( "nJetGood/I" ), #nJetSelected
      binning=[8,-0.5,7.5],
    ))

    plots.append(Plot(
      texX = 'N_{b-tag}', texY = 'Number of Events',
      attribute = TreeVariable.fromString( "nBTag/I" ), #nJetSelected
      binning=[4,-0.5,3.5],
    ))

    plots.append(Plot(
      texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events / 30 GeV',
      name = 'jet0_pt', attribute = lambda event, sample: event.JetGood_pt[0],
      binning=[600/30,0,600],
    ))

    plots.append(Plot(
      texX = 'p_{T}(subleading jet) (GeV)', texY = 'Number of Events / 30 GeV',
      name = 'jet1_pt', attribute = lambda event, sample: event.JetGood_pt[1],
      binning=[600/30,0,600],
    ))

    plots.append(Plot(
        name = "W_pt",
        texX = 'p_{T}(W) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = lambda event, sample:event.W_pt,
        binning=[20,0,400],
    ))

    # 3l training variables

    plots.append(Plot(
      texX = '#Delta\#phi(nonZ-l_{1}, Z_{1})', texY = 'Number of Events',
      name = 'nonZ1l1_Z1_deltaPhi', attribute = lambda event, sample: event.nonZ1l1_Z1_deltaPhi,
      binning=[20,0,pi],
    ))
    plots.append(Plot(
      texX = '#Delta#eta(nonZ-l_{1}, Z_{1})', texY = 'Number of Events',
      name = 'nonZ1l1_Z1_deltaEta', attribute = lambda event, sample: event.nonZ1l1_Z1_deltaEta,
      binning=[20,0,6],
    ))
    plots.append(Plot(
      texX = '#Delta R(nonZ-l_{1}, Z_{1})', texY = 'Number of Event',
      name = 'nonZ1l1_Z1_deltaR', attribute = lambda event, sample: event.nonZ1l1_Z1_deltaR,
      binning=[20,0,6],
    ))

    plots.append(Plot(
      texX = '#Delta R(jet_{0}, Z_{1})', texY = 'Number of Events',
      name = 'jet0_Z1_deltaR', attribute = lambda event, sample: event.jet0_Z1_deltaR,
      binning=[20,0,6],
    ))
    plots.append(Plot(
      texX = '#Delta R(jet_{0}, nonZ-l_{1})', texY = 'Number of Events',
      name = 'jet0_nonZ1l1_deltaR', attribute = lambda event, sample: event.jet0_nonZ1l1_deltaR,
      binning=[20,0,6],
    ))
    plots.append(Plot(
      texX = '#Delta R(jet_{1}, Z_{1})', texY = 'Number of Events',
      name = 'jet1_Z1_deltaR', attribute = lambda event, sample: event.jet1_Z1_deltaR,
      binning=[20,0,6],
    ))
    plots.append(Plot(
      texX = '#Delta R(jet_{1}, nonZ-l_{1})', texY = 'Number of Events',
      name = 'jet1_nonZ1l1', attribute = lambda event, sample: event.jet1_nonZ1l1_deltaR,
      binning=[20,0,6],
    ))
    
    for index in range(3):
        for abs_pdg in [11, 13]:
            lep_name = "mu" if abs_pdg==13 else "ele"
            plots.append(Plot(
              texX = 'p_{T}(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
              name = '%s%i_pt'%(lep_name, index), attribute = lep_getter("pt", index, abs_pdg),
              binning=[400/20,0,400],
            ))
            plots.append(Plot(
              texX = '#eta(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
              name = '%s%i_eta'%(lep_name, index), attribute = lep_getter("eta", index, abs_pdg),
              binning=[30,-3,3],
            ))
            plots.append(Plot(
              texX = '#phi(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
              name = '%s%i_phi'%(lep_name, index), attribute = lep_getter("phi", index, abs_pdg),
              binning=[30,-pi,pi],
            ))
            plots.append(Plot(
              texX = 'dxy(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
              name = '%s%i_dxy'%(lep_name, index), attribute = lep_getter("dxy", index, abs_pdg, functor = lambda x: abs(x)),
              binning=[50,0,0.05],
            ))
            plots.append(Plot(
              texX = 'dz(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
              name = '%s%i_dz'%(lep_name, index), attribute = lep_getter("dz", index, abs_pdg, functor = lambda x: abs(x)),
              binning=[50,0,0.05],
            ))
            plots.append(Plot(
              texX = 'IP_{3D}(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
              name = '%s%i_ip3d'%(lep_name, index), attribute = lep_getter("ip3d", index, abs_pdg, functor = lambda x: abs(x)),
              binning=[50,0,0.05],
            ))
            plots.append(Plot(
              texX = '#sigma(IP)_{3D}(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
              name = '%s%i_sip3d'%(lep_name, index), attribute = lep_getter("sip3d", index, abs_pdg, functor = lambda x: abs(x)),
              binning=[40,0,8],
            ))
            plots.append(Plot(
              texX = 'jetRelIso(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
              name = '%s%i_jetRelIso'%(lep_name, index), attribute = lep_getter("jetRelIso", index, abs_pdg),
              binning=[50,-.15,0.5],
            ))
            plots.append(Plot(
              texX = 'miniPFRelIso_all(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
              name = '%s%i_miniPFRelIso_all'%(lep_name, index), attribute = lep_getter("miniPFRelIso_all", index, abs_pdg),
              binning=[50,0,.5],
            ))
            plots.append(Plot(
              texX = 'pfRelIso03_all(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
              name = '%s%i_pfRelIso03_all'%(lep_name, index), attribute = lep_getter("pfRelIso03_all", index, abs_pdg),
              binning=[50,0,.5],
            ))
            plots.append(Plot(
              texX = 'mvaTTH(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
              name = '%s%i_mvaTTH'%(lep_name, index), attribute = lep_getter("mvaTTH", index, abs_pdg),
              binning=[24,-1.2,1.2],
            ))
            plots.append(Plot(
              texX = 'mvaTOP(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
              name = '%s%i_mvaTOP'%(lep_name, index), attribute = lep_getter("mvaTOP", index, abs_pdg),
              binning=[24,-1.2,1.2],
            ))
            plots.append(Plot(
              texX = 'charge(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
              name = '%s%i_charge'%(lep_name, index), attribute = lep_getter("pdgId", index, abs_pdg, functor = charge),
              binning=[3,-1,2],
            ))
            if lep_name == "mu":
                plots.append(Plot(
                  texX = 'segmentComp(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
                  name = '%s%i_segmentComp'%(lep_name, index), attribute = lep_getter("segmentComp", index, abs_pdg),
                  binning=[50,0,1],
                ))
                plots.append(Plot(
                  texX = 'nStations(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
                  name = '%s%i_nStations'%(lep_name, index), attribute = lep_getter("nStations", index, abs_pdg),
                  binning=[10,0,10],
                ))
                plots.append(Plot(
                  texX = 'nTrackerLayers(%s_{%i}) (GeV)'%(lep_name, index), texY = 'Number of Events',
                  name = '%s%i_nTrackerLayers'%(lep_name, index), attribute = lep_getter("nTrackerLayers", index, abs_pdg),
                  binning=[20,0,20],
                ))
            if lep_name == "ele":
                for cbIdFlag in vidNestedWPBitMapNamingList:
                    plots.append(Plot(
                      texX = '%s(%s_{%i}) (GeV)'%(cbIdFlag, lep_name, index), texY = 'Number of Events',
                      name = '%s%i_%s_Flag'%(lep_name, index, cbIdFlag), attribute = lep_getter("vidNestedWPBitmap", index, abs_pdg, functor = cbEleIdFlagGetter(cbIdFlag)),
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

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
