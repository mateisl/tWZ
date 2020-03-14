#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT, os
ROOT.gROOT.SetBatch(True)
import itertools
import copy
import array
import operator

from math                                import sqrt, cos, sin, pi, atan2, cosh
from RootTools.core.standard             import *
from tWZ.Tools.user                      import plot_directory
#from tWZ.tools.helpers         import deltaPhi, deltaR, add_histos
from Analysis.Tools.metFilters           import getFilterCut
#from tWZ.tools.cutInterpreter            import cutInterpreter
from Analysis.Tools.puProfileCache import *

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--noData',         action='store_true', default=False, help='also plot data?')
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--dataMCScaling',  action='store_true', help='Data MC scaling?', )
argParser.add_argument('--plot_directory', action='store', default='tWZ_v0')
argParser.add_argument('--era',            action='store', type=str, default="2016")
argParser.add_argument('--selection',      action='store', default='lepSel-njet2p-btag0')
args = argParser.parse_args()

#
# Logger
#
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"
if args.noData:                       args.plot_directory += "_noData"
#
# Make samples, will be searched for in the postProcessing directory
#
from Analysis.Tools.puReweighting import getReweightingFunction

if "2016" in args.era:
    year = 2016
elif "2017" in args.era:
    year = 2017
elif "2018" in args.era:
    year = 2018

logger.info( "Working in year %i", year )

#tWZ_sample = TWZ if args.nominalSignal else yt_TWZ_filter

if year == 2016:
    #from StopsDilepton.samples.nanoTuples_Run2016_17Jul2018_postProcessed import *
    from StopsDilepton.samples.nanoTuples_Summer16_postProcessed import *
    mc             = [ [ tWZ_sample, TTZtoLLNuNu, TTX_rare_TWZ, TZQ, WZ_amcatnlo, rare, ZZ, nonprompt_TWZ_3l]
elif year == 2017:
    from StopsDilepton.samples.nanoTuples_Fall17_postProcessed import *
    from StopsDilepton.samples.nanoTuples_Run2017_31Mar2018_postProcessed import *
    mc             = [ Top_pow_17, TTXNoZ_17, TTZ_17, multiBoson_17, DY_HT_LO_17]
elif year == 2018:
    from StopsDilepton.samples.nanoTuples_Run2018_PromptReco_postProcessed import *
    from StopsDilepton.samples.nanoTuples_Autumn18_postProcessed import *
    mc             = [ Top_pow_18, TTXNoZ_18, TTZ_18, multiBoson_18, DY_HT_LO_18]


try:
  data_sample = eval(args.era)
except Exception as e:
  logger.error( "Didn't find %s", args.era )
  raise e

lumi_scale                 = data_sample.lumi/1000
data_sample.scale          = 1.
for sample in mc:
    sample.scale          = lumi_scale

#if args.small:
#    for sample in mc + [data_sample]:
#        sample.normalization = 1.
#        #sample.reduceFiles( factor = 40 )
#        sample.reduceFiles( to=1)
#        sample.scale /= sample.normalization
#
##
## Text on the plots
##
#tex = ROOT.TLatex()
#tex.SetNDC()
#tex.SetTextSize(0.04)
#tex.SetTextAlign(11) # align right
#
#def drawObjects( plotData, dataMCScale, lumi_scale ):
#    lines = [
#      (0.15, 0.95, 'CMS Preliminary' if plotData else 'CMS Simulation'), 
#      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
#    ]
#    if "mt2ll100" in args.selection and args.noData: lines += [(0.55, 0.5, 'M_{T2}(ll) > 100 GeV')] # Manually put the mt2ll > 100 GeV label
#    return [tex.DrawLatex(*l) for l in lines] 
#
#def drawPlots(plots, mode, dataMCScale):
#  for log in [False, True]:
#    plot_directory_ = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, mode + ("_log" if log else ""), args.selection)
#    for plot in plots:
#      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
#      if not args.noData: 
#        if mode == "all": plot.histos[1][0].legendText = "Data"
#        if mode == "SF":  plot.histos[1][0].legendText = "Data (SF)"
#
#      _drawObjects = []
#
#      if ("u_para" in plot.name or "u_perp" in plot.name) and not args.noData:
#          h_mc   = plot.histos_added[0][0].Clone()
#          h_data = plot.histos_added[1][0].Clone()
#          if h_mc.Integral()>0:
#              h_mc.Scale(h_data.Integral()/h_mc.Integral())
#          q_mc   = tuple(get_quantiles( h_mc ))
#          q_data = tuple(get_quantiles( h_data ))
#          median_shift = q_data[2]-q_mc[2]
#          sigma1_ratio = (q_data[3]-q_data[1])/(q_mc[3]-q_mc[1]) if q_mc[3]-q_mc[1]!=0 else 0
#          sigma2_ratio = (q_data[4]-q_data[0])/(q_mc[4]-q_mc[0]) if q_mc[4]-q_mc[0]!=0 else 0
#
#          _drawObjects.append( tex.DrawLatex(0.22, 0.62, '#Delta(med): %+3.1f   1#sigma: %4.3f  2#sigma  %4.3f' % ( median_shift, sigma1_ratio, sigma2_ratio) ) )
#
#      if isinstance( plot, Plot):
#          plotting.draw(plot,
#            plot_directory = plot_directory_,
#            ratio = {'yRange':(0.1,1.9)} if not args.noData else None,
#            logX = False, logY = log, sorting = not (args.splitMET or args.splitMETSig or args.splitNvtx), # and args.sorting is not None,
#            yRange = (0.03, "auto") if log else (0.001, "auto"),
#            scaling = {0:1} if args.dataMCScaling else {},
#            legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
#            drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ) + _drawObjects,
#            copyIndexPHP = True, extensions = ["png", "pdf"],
#          )
#            
##
## Read variables and sequences
##
#read_variables = ["weight/F", "l1_pt/F", "dl_phi/F", "dl_pt/F", "l2_pt/F", "l1_eta/F" , "l1_phi/F", "l2_eta/F", "l2_phi/F", "JetGood[pt/F,eta/F,phi/F]", "dl_mass/F", "dl_eta/F", "dl_mt2ll/F", "dl_mt2bb/F", "dl_mt2blbl/F", "met_pt/F", "met_phi/F", "MET_significance/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I", "RawMET_pt/F", "RawMET_phi/F"]
#read_variables += ["l2_eleIndex/I"]
#read_variables+= ["event/l", "luminosityBlock/I", "run/I"]
#if "2017" in args.era:
#    read_variables.append( "MET_pt_min/F" ) 
#
#sequence = []
#
##
## Loop over channels
##
#yields     = {}
#allPlots   = {}
#allModes   = ['mumu','mue','ee']
#for index, mode in enumerate(allModes):
#  yields[mode] = {}
#
#  data_sample.setSelectionString([getFilterCut(isData=True, year=year), getLeptonSelection(mode)])
#  data_sample.name           = "data"
#  data_sample.read_variables = ["event/I","run/I", "reweightHEM/F"]
#  data_sample.style          = styles.errorStyle(ROOT.kBlack)
#  weight_ = lambda event, sample: event.weight*event.reweightHEM
#
#  #data_sample_filtered = copy.deepcopy( data_sample )
#  #data_sample_filtered.style = styles.errorStyle(ROOT.kRed)
#  #data_sample_filtered.weight = lambda event, sample: event.weight*event.passes_veto
#  #data_sample_filtered.name   += "_filtered"
#  #data_sample_filtered.texName+= " (filtered)"
#
#  for sample in mc + signals:
#    sample.read_variables = ['reweightPU/F', 'reweightL1Prefire/F', 'Pileup_nTrueInt/F', 'reweightDilepTrigger/F','reweightLeptonSF/F','reweightBTag_SF/F', 'reweightLeptonTrackingSF/F', 'GenMET_pt/F', 'GenMET_phi/F', 'reweightHEM/F', 'reweightLeptonHit0SF/F', 'reweightLeptonSip3dSF/F']
#    # Need individual pu reweighting functions for each sample in 2017, so nTrueInt_puRW is only defined here
#
#    if args.rwHit0:
#        weight_Hit0  = lambda event: 1 
#    else:
#        weight_Hit0  = operator.attrgetter( 'reweightLeptonHit0SF' ) 
#    if args.rwSip3d:
#        weight_sip3d = lambda event: 1 
#    else:
#        weight_sip3d = operator.attrgetter( 'reweightLeptonSip3dSF' ) 
#
#    if args.reweightPU and args.reweightPU not in ["noPUReweighting", "nvtx"]:
#        sample.read_variables.append( 'reweightPU/F' if args.reweightPU=='Central' else 'reweightPU%s/F'%args.reweightPU )
#
#    if args.reweightPU == "noPUReweighting":
#        sample.weight         = lambda event, sample: event.reweightDilepTrigger*event.reweightLeptonSF*event.reweightBTag_SF*event.reweightLeptonTrackingSF*event.reweightL1Prefire
#    elif args.reweightPU == "nvtx":
#        sample.weight         = lambda event, sample: nvtx_puRW(event.PV_npvsGood) * event.reweightDilepTrigger*event.reweightLeptonSF*event.reweightBTag_SF*event.reweightLeptonTrackingSF*event.reweightL1Prefire
#    elif args.reweightPU:
#        pu_getter = operator.attrgetter( 'reweightPU' if args.reweightPU=='Central' else 'reweightPU%s'%args.reweightPU )
#        sample.weight         = lambda event, sample: pu_getter(event) * event.reweightDilepTrigger*event.reweightLeptonSF*event.reweightBTag_SF*event.reweightLeptonTrackingSF*event.reweightL1Prefire*weight_sip3d(event)*weight_Hit0(event)
#    else: #default
#        sample.weight         = lambda event, sample: event.reweightPU*event.reweightDilepTrigger*event.reweightLeptonSF*event.reweightBTag_SF*event.reweightLeptonTrackingSF*event.reweightL1Prefire*weight_sip3d(event)*weight_Hit0(event)
#
#    sample.setSelectionString([getFilterCut(isData=False, year=year), getLeptonSelection(mode)])
#
#  if args.splitMET:
#    mc_ = splitMetMC(mc)
#  elif args.splitMETSig:
#    mc_ = splitMetSigMC(mc)
#  elif args.splitNvtx:
#    mc_ = splitNvtxMC(mc)
#  else:
#    mc_ = mc 
#
#  for sample in mc_: sample.style = styles.fillStyle(sample.color)
#
#  if not args.noData:
#    stack = Stack(mc_, data_sample)#, data_sample_filtered)
#  else:
#    stack = Stack(mc_)
#
#  stack.extend( [ [s] for s in signals ] )
#
#  # Use some defaults
#  Plot  .setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = cutInterpreter.cutString(args.selection), addOverFlowBin='upper', histo_class=ROOT.TH1D)
#  Plot2D.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = cutInterpreter.cutString(args.selection), histo_class=ROOT.TH2D)
#  
#  plots   = []
#  plots2D = []
#  plots.append(Plot(
#    name = 'yield', texX = 'yield', texY = 'Number of Events',
#    attribute = lambda event, sample: 0.5 + index,
#    binning=[3, 0, 3],
#  ))
#
#  plots.append(Plot(
#    name = 'PV_npvsGood', texX = 'N_{PV} (good)', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "PV_npvsGood/I" ),
#    binning=[100,0,100],
#  ))
#    
#  plots.append(Plot(
#    name = 'PV_npvs', texX = 'N_{PV} (total)', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "PV_npvs/I" ),
#    binning=[100,0,100],
#  ))
#
#  plots.append(Plot(
#      texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
#      attribute = TreeVariable.fromString( "met_pt/F" ),
#      binning=[400/20,0,400],
#  ))
#
#  plots2D.append(Plot2D(
#    name = 'allJets_occupancy_pt30',
#    stack = stack,
#    attribute = (
#      lambda event, sample: [ j['eta'] for j in event.jets if j['pt']>30 ],
#      lambda event, sample: [ j['phi'] for j in event.jets if j['pt']>30 ],
#    ),
#    texX = '#eta (all jets)', texY = '#phi',
#    binning=[52, -5.2, 5.2, 32, -pi, pi],
#  ))
#
#  plots.append(Plot(
#    name = 'allJets_eta_pt30',
#    stack = stack,
#    attribute = lambda event, sample: [ j['eta'] for j in event.jets if j['pt']>30 ],
#    texX = '#eta (all jets)', 
#    binning=[52, -5.2, 5.2],
#  ))
#
#  plots.append(Plot(
#    name = 'allJets_eta_pt50',
#    stack = stack,
#    attribute = lambda event, sample: [ j['eta'] for j in event.jets if j['pt']>50 ],
#    texX = '#eta (all jets)', 
#    binning=[52, -5.2, 5.2],
#  ))
#
#  plots.append(Plot(
#    name = 'allJets_eta_pt100',
#    stack = stack,
#    attribute =  lambda event, sample: [ j['eta'] for j in event.jets if j['pt']>100 ],
#    texX = '#eta (all jets)', 
#    binning=[52, -5.2, 5.2],
#  ))
#
#  if "2017" in args.era:
#    plots.append(Plot(
#        texX = 'min E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
#        attribute = TreeVariable.fromString( "MET_pt_min/F" ),
#        binning=[400/20,0,400],
#    ))
#    plots.append(Plot(name = "MET_pt_min_delta",
#        texX = '#Delta min E_{T}^{miss} (GeV)', texY = 'Number of Events / 10 GeV',
#        attribute = lambda event, sample: event.met_pt - event.MET_pt_min,
#        binning=[200/10,0,200],
#    ))
#
##Sum$(abs(Jet_eta)>2.6&&abs(Jet_eta)<3.1&&Jet_pt<50)==0
#
#  #plots.append(Plot( name = "met_pt_raw",
#  #    texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
#  #    attribute = TreeVariable.fromString( "RawMET_pt/F" ),
#  #    binning=[400/20,0,400],
#  #))
#
#  plots.append(Plot(
#      texX = 'E_{T}^{miss} significance', texY = 'Number of Events',
#      attribute = TreeVariable.fromString( "MET_significance/F" ),
#      binning=[34,0,102],
#  ))
#
#  plots.append(Plot(
#      texX = '#phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
#      attribute = TreeVariable.fromString( "met_phi/F" ),
#      binning=[10,-pi,pi],
#  ))
#
#  #plots.append(Plot( name = "met_phi_corr",
#  #    texX = '#phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
#  #    attribute = lambda event, sample: event.met_phi_corr,
#  #    binning=[10,-pi,pi],
#  #))
#
#  #plots.append(Plot( name = "met_phi_raw",
#  #    texX = 'raw #phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
#  #    attribute = TreeVariable.fromString( "RawMET_phi/F" ),
#  #    binning=[10,-pi,pi],
#  #))
#
#  #plots.append(Plot( name = "met_phi_raw_corr",
#  #    texX = 'raw #phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
#  #    attribute = lambda event, sample: event.RawMET_phi_corr,
#  #    binning=[10,-pi,pi],
#  #))
#
#  if not args.blinded:
#    plots.append(Plot(
#      texX = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV',
#      attribute = TreeVariable.fromString( "dl_mt2ll/F" ),
#      binning=[300/20,0,300],
#    ))
#
#    plots.append(Plot( name = "dl_mt2ll_coarse",
#      texX = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV',
#      attribute = TreeVariable.fromString( "dl_mt2ll/F" ),
#      binning=Binning.fromThresholds([0,20,40,60,80,100,140,240,340]),
#    ))
#    #plots.append(Plot( name = "dl_mt2ll_raw",
#    #  texX = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV',
#    #  attribute = lambda event, sample: event.dl_mt2ll_raw,
#    #  binning=[300/20,0,300],
#    #))
#
#  #plots.append(Plot( name = "met_pt_corr",
#  #    texX = 'corr E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
#  #    attribute = lambda event, sample: event.met_pt_corr,
#  #    binning=[400/20,0,400],
#  #))
#
#  #plots.append(Plot( name = "met_pt_raw_corr",
#  #    texX = 'corr E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
#  #    attribute = lambda event, sample: event.RawMET_pt_corr,
#  #    binning=[400/20,0,400],
#  #))
#    
#  #plots.append(Plot( name = "dl_mt2ll_corr",
#  #    texX = 'corr M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV',
#  #    attribute = lambda event, sample: event.dl_mt2ll_corr,
#  #    binning=[300/20,0,300],
#  #))
#
#  #plots.append(Plot( name = "dl_mt2ll_raw_corr",
#  #    texX = 'raw corr M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV',
#  #    attribute = lambda event, sample: event.dl_mt2ll_raw_corr,
#  #    binning=[300/20,0,300],
#  #))
#
#  plots.append(Plot( name = "qT",
#    texX = 'q_{T} (GeV)', texY = 'Number of Events / 50 GeV',
#    attribute = lambda event, sample: sqrt((event.l1_pt*cos(event.l1_phi) + event.l2_pt*cos(event.l2_phi) + event.met_pt*cos(event.met_phi))**2 + (event.l1_pt*sin(event.l1_phi) + event.l2_pt*sin(event.l2_phi) + event.met_pt*sin(event.met_phi))**2),
#    binning= [1000/50,0,1000]),
#  )
#
#  plots.append(Plot(
#    texX = 'number of jets', texY = 'Number of Events',
#    attribute = TreeVariable.fromString('nJetGood/I'),
#    binning=[14,0,14],
#  ))
#  if "2017" in args.era:  
#    plots.append(Plot( name = "nJet_EE",
#      texX = 'number of jets', texY = 'Number of Events',
#      attribute = lambda event, sample: event.nJet_EE,
#      binning=[14,0,14],
#    ))
#    plots.append(Plot( name = "nJet_EE_pt30To50",
#      texX = 'number of jets', texY = 'Number of Events',
#      attribute = lambda event, sample: event.nJet_EE_pt30To50,
#      binning=[14,0,14],
#    ))
#    plots.append(Plot( name = "nJet_EE_ptTo50",
#      texX = 'number of jets', texY = 'Number of Events',
#      attribute = lambda event, sample: event.nJet_EE_ptTo50,
#      binning=[14,0,14],
#    ))
#    plots.append(Plot( name = "nJet_EE_ptTo40",
#      texX = 'number of jets', texY = 'Number of Events',
#      attribute = lambda event, sample: event.nJet_EE_ptTo40,
#      binning=[14,0,14],
#    ))
#    plots.append(Plot( name = "nJet_EE_ptTo30",
#      texX = 'number of jets', texY = 'Number of Events',
#      attribute = lambda event, sample: event.nJet_EE_ptTo30,
#      binning=[14,0,14],
#    ))
#    plots.append(Plot( name = "nJet_EE_pt50",
#      texX = 'number of jets', texY = 'Number of Events',
#      attribute = lambda event, sample: event.nJet_EE_pt50,
#      binning=[14,0,14],
#    ))
#    plots.append(Plot( name = "nJet_EE_pt40",
#      texX = 'number of jets', texY = 'Number of Events',
#      attribute = lambda event, sample: event.nJet_EE_pt40,
#      binning=[14,0,14],
#    ))
#    plots.append(Plot( name = "nJet_EE_pt30",
#      texX = 'number of jets', texY = 'Number of Events',
#      attribute = lambda event, sample: event.nJet_EE_pt30,
#      binning=[14,0,14],
#    ))
#    plots.append(Plot( name = "badJetE",
#      texX = 'badEEJetEnergy', texY = 'Number of Events',
#      attribute = lambda event, sample: event.badJetE,
#      binning=[40,0,400],
#    ))
#    plots.append(Plot( name = "badJetPt",
#      texX = 'badEEJetPt', texY = 'Number of Events',
#      attribute = lambda event, sample: event.badJetPt,
#      binning=[40,0,200],
#    ))
#
#  plots.append(Plot(
#    texX = 'number of medium b-tags (CSVM)', texY = 'Number of Events',
#    attribute = TreeVariable.fromString('nBTag/I'),
#    binning=[8,0,8],
#  ))
#
#  #plots.append(Plot(
#  #  texX = 'H_{T} (GeV)', texY = 'Number of Events / 25 GeV',
#  #  attribute = TreeVariable.fromString( "ht/F" ),
#  #  binning=[500/25,0,600],
#  #))
#
#  plots.append(Plot(
#    texX = 'm(ll) of leading dilepton (GeV)', texY = 'Number of Events / 4 GeV',
#    attribute = TreeVariable.fromString( "dl_mass/F" ),
#    binning=[200/4,0,200],
#  ))
#
#  plots.append(Plot(
#    texX = 'p_{T}(ll) (GeV)', texY = 'Number of Events / 10 GeV',
#    attribute = TreeVariable.fromString( "dl_pt/F" ),
#    binning=[20,0,400],
#  ))
#
#  plots.append(Plot(
#      texX = '#eta(ll) ', texY = 'Number of Events',
#      name = 'dl_eta', attribute = lambda event, sample: abs(event.dl_eta), read_variables = ['dl_eta/F'],
#      binning=[10,0,3],
#  ))
#
#  plots.append(Plot(
#    texX = '#phi(ll)', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "dl_phi/F" ),
#    binning=[10,-pi,pi],
#  ))
#
#  plots.append(Plot(
#    texX = 'Cos(#Delta#phi(ll, E_{T}^{miss}))', texY = 'Number of Events',
#    name = 'cosZMetphi',
#    attribute = lambda event, sample: cos( event.dl_phi - event.met_phi ), 
#    read_variables = ["met_phi/F", "dl_phi/F"],
#    binning = [10,-1,1],
#  ))
#
#  plots.append(Plot(
#    texX = 'p_{T}(l_{1}) (GeV)', texY = 'Number of Events / 15 GeV',
#    attribute = TreeVariable.fromString( "l1_pt/F" ),
#    binning=[20,0,300],
#  ))
#
#  plots.append(Plot(
#    texX = '#eta(l_{1})', texY = 'Number of Events',
#    name = 'l1_eta', attribute = lambda event, sample: abs(event.l1_eta), read_variables = ['l1_eta/F'],
#    binning=[15,0,3],
#  ))
#
#  plots.append(Plot(
#    texX = '#phi(l_{1})', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "l1_phi/F" ),
#    binning=[10,-pi,pi],
#  ))
#
#  plots.append(Plot(
#    texX = 'I_{mini}(l_{1})', texY = 'Number of Events',
#    name = 'l1_miniRelIso', attribute = lambda event, sample: event.l1_miniRelIso, read_variables = ['l1_miniRelIso/F'],
#    binning=[20,0,.5],
#  ))
#
#  plots.append(Plot(
#    texX = 'pdgId(l1)', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "l1_pdgId/I" ),
#    binning=[30,-15,15],
#  ))
#
#  plots.append(Plot(
#    texX = 'p_{T}(l_{2}) (GeV)', texY = 'Number of Events / 15 GeV',
#    attribute = TreeVariable.fromString( "l2_pt/F" ),
#    binning=[20,0,300],
#  ))
#
#  plots.append(Plot(
#    texX = '#eta(l_{2})', texY = 'Number of Events',
#    name = 'l2_eta', attribute = lambda event, sample: abs(event.l2_eta), read_variables = ['l2_eta/F'],
#    binning=[15,0,3],
#  ))
#
#  plots.append(Plot(
#    texX = '#phi(l_{2})', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "l2_phi/F" ),
#    binning=[10,-pi,pi],
#  ))
#  
#  plots.append(Plot(
#    texX = 'I_{mini}(l_{2})', texY = 'Number of Events',
#    name = 'l2_miniRelIso', attribute = lambda event, sample: event.l2_miniRelIso, read_variables = ['l2_miniRelIso/F'],
#    binning=[20,0,.5],
#  ))
#
#  plots.append(Plot(
#    texX = 'pdgId(l2)', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "l2_pdgId/I" ),
#    binning=[30,-15,15],
#  ))
#
#  # plot trailing lepton quantities
#  plots2D.append(Plot2D(
#    name = 'eff_num_trailing_ele_sip3Dlt4',
#    stack = stack,
#    attribute = (
#      lambda event, sample: event.l2_eta if event.l2_eleIndex>0 and event.Electron_sip3d[event.l2_eleIndex]<4 else float('nan'),
#      lambda event, sample: event.l2_pt  if event.l2_eleIndex>0 and event.Electron_sip3d[event.l2_eleIndex]<4 else float('nan'),
#    ),
#    texX = 'trailing lepton #eta', texY = 'trailing lepton p_{T} (GeV)',
#    binning=[Binning.fromThresholds([-2.5, -2, -1.566, -1.444, -0.8, 0, 0.8, 1.444, 1.566, 2.0, 2.5]), Binning.fromThresholds([20, 35, 50, 100, 200])],
#  ))
#  plots2D.append(Plot2D(
#    name = 'eff_den_trailing_ele_sip3Dlt4', #keep the cut in the name, but this is the denominator, so remove it from the lambda function
#    stack = stack,
#    attribute = (
#      lambda event, sample: event.l2_eta if event.l2_eleIndex>0 else float('nan'),
#      lambda event, sample: event.l2_pt  if event.l2_eleIndex>0 else float('nan'),
#    ),
#    texX = 'trailing lepton #eta', texY = 'trailing lepton p_{T} (GeV)',
#    binning=[Binning.fromThresholds([-2.5, -2, -1.566, -1.444, -0.8, 0, 0.8, 1.444, 1.566, 2.0, 2.5]), Binning.fromThresholds([20, 35, 50, 100, 200])],
#  ))
#
#  # Plots only when at least one jet:
#  if args.selection.count('njet2') or args.selection.count('njet1') or args.selection.count('njet01'):
#
#    plots2D.append(Plot2D(
#      name = 'leading_jet_occ',
#      stack = stack,
#      attribute = (
#        lambda event, sample: event.JetGood_eta[0],
#        lambda event, sample: event.JetGood_phi[0],
#      ),
#      texX = '#eta(leading jet) (GeV)', texY = '#phi(leading jet) (GeV)',
#      binning=[16, -3.0, 3.0, 10, -pi, pi],
#    ))
#
#    plots.append(Plot(
#      texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events / 30 GeV',
#      name = 'jet1_pt', attribute = lambda event, sample: event.JetGood_pt[0],
#      binning=[600/30,0,600],
#    ))
#
#    plots.append(Plot(
#      texX = '#eta(leading jet) (GeV)', texY = 'Number of Events',
#      name = 'jet1_eta', attribute = lambda event, sample: event.JetGood_eta[0],
#      binning=[20,-3,3],
#    ))
#
#    plots.append(Plot(
#      texX = '#phi(leading jet) (GeV)', texY = 'Number of Events',
#      name = 'jet1_phi', attribute = lambda event, sample: event.JetGood_phi[0],
#      binning=[20,-pi,pi],
#    ))
#
#    plots.append(Plot(
#      name = 'cosMetJet1phi',
#      texX = 'Cos(#Delta#phi(E_{T}^{miss}, leading jet))', texY = 'Number of Events',
#      attribute = lambda event, sample: cos( event.met_phi - event.JetGood_phi[0]), 
#      read_variables = ["met_phi/F", "JetGood[phi/F]"],
#      binning = [10,-1,1],
#    ))
#    
#    plots.append(Plot(
#      name = 'cosMetJet1phi_smallBinning',
#      texX = 'Cos(#Delta#phi(E_{T}^{miss}, leading jet))', texY = 'Number of Events',
#      attribute = lambda event, sample: cos( event.met_phi - event.JetGood_phi[0] ) , 
#      read_variables = ["met_phi/F", "JetGood[phi/F]"],
#      binning = [20,-1,1],
#    ))
#
#    plots.append(Plot(
#      name = 'cosZJet1phi',
#      texX = 'Cos(#Delta#phi(Z, leading jet))', texY = 'Number of Events',
#      attribute = lambda event, sample: cos( event.dl_phi - event.JetGood_phi[0] ),
#      read_variables =  ["dl_phi/F", "JetGood[phi/F]"],
#      binning = [10,-1,1],
#    ))
#    plots.append(Plot( name = "u_para", 
#      texX = "u_{#parallel} (GeV)", texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: - event.met_pt*cos(event.met_phi-event.dl_phi),
#      binning=[80, -200,200],
#    ))
#    plots.append(Plot( name = "u_perp", 
#      texX = "u_{#perp} (GeV)", texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: - event.met_pt*cos(event.met_phi-(event.dl_phi-pi/2)),
#      binning=[80, -200,200],
#    ))
#    #plots.append(Plot( name = "u_para_raw", 
#    #  texX = "u_{#parallel} (GeV)", texY = 'Number of Events / 30 GeV',
#    #  attribute = lambda event, sample: - event.RawMET_pt*cos(event.RawMET_phi-event.dl_phi),
#    #  binning=[80, -200,200],
#    #))
#    #plots.append(Plot( name = "u_perp_raw", 
#    #  texX = "u_{#perp} (GeV)", texY = 'Number of Events / 30 GeV',
#    #  attribute = lambda event, sample: - event.RawMET_pt*cos(event.RawMET_phi-(event.dl_phi-pi/2)),
#    #  binning=[80, -200,200],
#    #))
#    #plots.append(Plot( name = "u_para_corr", 
#    #  texX = "u_{#parallel} corr. (GeV)", texY = 'Number of Events / 30 GeV',
#    #  attribute = lambda event, sample: - event.met_pt_corr*cos(event.met_phi_corr-event.dl_phi),
#    #  binning=[80, -200,200],
#    #))
#    #plots.append(Plot( name = "u_perp_corr", 
#    #  texX = "u_{#perp} corr. (GeV)", texY = 'Number of Events / 30 GeV',
#    #  attribute = lambda event, sample: - event.met_pt_corr*cos(event.met_phi_corr-(event.dl_phi-pi/2)),
#    #  binning=[80, -200,200],
#    #))
#    #plots.append(Plot( name = "u_para_raw_corr", 
#    #  texX = "u_{#parallel} corr. (GeV)", texY = 'Number of Events / 30 GeV',
#    #  attribute = lambda event, sample: - event.RawMET_pt_corr*cos(event.RawMET_phi_corr-event.dl_phi),
#    #  binning=[80, -200,200],
#    #))
#    #plots.append(Plot( name = "u_perp_raw_corr", 
#    #  texX = "u_{#perp} corr. (GeV)", texY = 'Number of Events / 30 GeV',
#    #  attribute = lambda event, sample: - event.RawMET_pt_corr*cos(event.RawMET_phi_corr-(event.dl_phi-pi/2)),
#    #  binning=[80, -200,200],
#    #))
#
#    if args.plotUPara:
#        # u_para u_perp closure plots
#        u_para_binning   =  [ i*20 for i in range(-10, 11) ]
#        qt_binning    = [0, 50, 100, 150, 200, 300 ]
#        qt_bins = [ (qt_binning[i],qt_binning[i+1]) for i in range(len(qt_binning)-1) ]
#        var_binning   = [ pi*(i-5)/5. for i in range(0,11) ]
#        var_bins      = [ (var_binning[i],var_binning[i+1]) for i in range(len(var_binning)-1) ]
#        #for var_bin in var_bins:
#        #    for qt_bin in qt_bins:
#        #        postfix = "phill_%3.2f_%3.2f_qt_%i_%i"%( var_bin[0], var_bin[1], qt_bin[0], qt_bin[1] )
#        #        plots.append(Plot( name = "u_para_" + postfix, 
#        #          texX = "u_{#parallel} (GeV)", texY = 'Number of Events / 30 GeV',
#        #          attribute = lambda event, sample: - event.met_pt*cos(event.met_phi-event.dl_phi),
#        #          weight = recoil_weight(var_bin, qt_bin),
#        #          binning=[80, -200,200],
#        #        ))
#        #        plots.append(Plot( name = "u_perp_" + postfix, 
#        #          texX = "u_{#perp} (GeV)", texY = 'Number of Events / 30 GeV',
#        #          attribute = lambda event, sample: - event.met_pt*cos(event.met_phi-(event.dl_phi-pi/2)),
#        #          weight = recoil_weight(var_bin, qt_bin),
#        #          binning=[80, -200,200],
#        #        ))
#        #        plots.append(Plot( name = "u_para_corr_" + postfix, 
#        #          texX = "u_{#parallel} corr. (GeV)", texY = 'Number of Events / 30 GeV',
#        #          attribute = lambda event, sample: - event.met_pt_corr*cos(event.met_phi_corr-event.dl_phi),
#        #          weight = recoil_weight(var_bin, qt_bin),
#        #          binning=[80, -200,200],
#        #        ))
#        #        plots.append(Plot( name = "u_perp_corr_" + postfix, 
#        #          texX = "u_{#perp} corr. (GeV)", texY = 'Number of Events / 30 GeV',
#        #          attribute = lambda event, sample: - event.met_pt_corr*cos(event.met_phi_corr-(event.dl_phi-pi/2)),
#        #          weight = recoil_weight(var_bin, qt_bin),
#        #          binning=[80, -200,200],
#        #        ))
#
##  # Plots only when at least two jets:
#  if args.selection.count('njet2'):
#    plots.append(Plot(
#      texX = 'p_{T}(2nd leading jet) (GeV)', texY = 'Number of Events / 30 GeV',
#      name = 'jet2_pt', attribute = lambda event, sample: event.JetGood_pt[1],
#      binning=[600/30,0,600],
#    ))
#
#    plots.append(Plot(
#      texX = '#eta(2nd leading jet) (GeV)', texY = 'Number of Events',
#      name = 'jet2_eta', attribute = lambda event, sample: event.JetGood_eta[1],
#      binning=[20,-3,3],
#    ))
#
#    plots.append(Plot(
#      texX = '#phi(2nd leading jet) (GeV)', texY = 'Number of Events',
#      name = 'jet2_phi', attribute = lambda event, sample: event.JetGood_phi[1],
#      binning=[20,-pi,pi],
#    ))
#
#    plots.append(Plot(
#      name = 'cosMetJet2phi',
#      texX = 'Cos(#Delta#phi(E_{T}^{miss}, second jet))', texY = 'Number of Events',
#      attribute = lambda event, sample: cos( event.met_phi - event.JetGood_phi[1] ) , 
#      read_variables = ["met_phi/F", "JetGood[phi/F]"],
#      binning = [10,-1,1],
#    ))
#    
#    plots.append(Plot(
#      name = 'cosMetJet2phi_smallBinning',
#      texX = 'Cos(#Delta#phi(E_{T}^{miss}, second jet))', texY = 'Number of Events',
#      attribute = lambda event, sample: cos( event.met_phi - event.JetGood_phi[1] ) , 
#      read_variables = ["met_phi/F", "JetGood[phi/F]"],
#      binning = [20,-1,1],
#    ))
#
#    plots.append(Plot(
#      name = 'cosZJet2phi',
#      texX = 'Cos(#Delta#phi(Z, 2nd leading jet))', texY = 'Number of Events',
#      attribute = lambda event, sample: cos( event.dl_phi - event.JetGood_phi[0] ),
#      read_variables = ["dl_phi/F", "JetGood[phi/F]"],
#      binning = [10,-1,1],
#    ))
#
#    plots.append(Plot(
#      name = 'cosJet1Jet2phi',
#      texX = 'Cos(#Delta#phi(leading jet, 2nd leading jet))', texY = 'Number of Events',
#      attribute = lambda event, sample: cos( event.JetGood_phi[1] - event.JetGood_phi[0] ) ,
#      read_variables =  ["JetGood[phi/F]"],
#      binning = [10,-1,1],
#    ))
#
#    plots.append(Plot(
#      texX = 'M_{T2}(bb) (GeV)', texY = 'Number of Events / 30 GeV',
#      attribute = TreeVariable.fromString( "dl_mt2bb/F" ),
#      binning=[420/30,70,470],
#    ))
#
#    plots.append(Plot(
#      texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 30 GeV',
#      attribute = TreeVariable.fromString( "dl_mt2blbl/F" ),
#      binning=[420/30,0,400],
#
#    ))
#
#    plots.append(Plot( name = "dl_mt2blbl_coarse",       # SR binning of MT2ll
#      texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 30 GeV',
#      attribute = TreeVariable.fromString( "dl_mt2blbl/F" ),
#      binning=Binning.fromThresholds([0,20,40,60,80,100,120,140,160,200,250,300,350]),
#    ))
#
#    plots.append(Plot( name = "dl_mt2blbl_bj_dR",       # SR binning of MT2ll
#      texX = 'dR(bj0,bj1))', texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: event.mt2blbl_bj_dR,
#      binning=[12,0,6],
#    ))
#
#    plots.append(Plot( name = "dl_mt2blbl_bj_dPhi",       # SR binning of MT2ll
#      texX = 'dPhi(bj0,bj1))', texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: event.mt2blbl_bj_dPhi,
#      binning=[12,0,pi],
#    ))
#
#    plots.append(Plot( name = "dl_mt2blbl_bj_dEta",       # SR binning of MT2ll
#      texX = 'dEta(bj0,bj1))', texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: event.mt2blbl_bj_dEta,
#      binning=[12,-6,6],
#    ))
#
#    plots.append(Plot( name = "dl_mt2blbl_bj_mass",       # SR binning of MT2ll
#      texX = 'mass(bj0,bj1))', texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: event.mt2blbl_bj_mass,
#      binning=[12,0,600],
#    ))
#
#    plots.append(Plot( name = "dl_mt2blbl_bj0_pt",       # SR binning of MT2ll
#      texX = 'bj0(pt)', texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: event.bj0['pt'] if event.bj0 is not None else float('nan'),
#      binning=[12,0,300],
#    ))
#    plots.append(Plot( name = "dl_mt2blbl_bj1_pt",       # SR binning of MT2ll
#      texX = 'bj1(pt)', texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: event.bj1['pt'] if event.bj1 is not None else float('nan'),
#      binning=[12,0,300],
#    ))
#    plots.append(Plot( name = "dl_mt2blbl_bj0_phi",       # SR binning of MT2ll
#      texX = 'bj0(phi)', texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: event.bj0['phi'] if event.bj0 is not None else float('nan'),
#      binning=[12,-pi,pi],
#    ))
#    plots.append(Plot( name = "dl_mt2blbl_bj1_phi",       # SR binning of MT2ll
#      texX = 'bj1(phi)', texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: event.bj1['phi'] if event.bj1 is not None else float('nan'),
#      binning=[12,-pi,pi],
#    ))
#    plots.append(Plot( name = "dl_mt2blbl_bj0_eta",       # SR binning of MT2ll
#      texX = 'bj0(eta)', texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: event.bj0['eta'] if event.bj0 is not None else float('nan'),
#      binning=[12,-3,3],
#    ))
#    plots.append(Plot( name = "dl_mt2blbl_bj1_eta",       # SR binning of MT2ll
#      texX = 'bj1(eta)', texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: event.bj1['eta'] if event.bj1 is not None else float('nan'),
#      binning=[12,-3,3],
#    ))
#    plots.append(Plot( name = "dl_mt2blbl_bj0_DeepCSV",       # SR binning of MT2ll
#      texX = 'bj0(DeepCSV)', texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: event.bj0['btagDeepB'] if event.bj0 is not None else float('nan'),
#      binning=[10,0,1],
#    ))
#    plots.append(Plot( name = "dl_mt2blbl_bj1_DeepCSV",       # SR binning of MT2ll
#      texX = 'bj1(DeepCSV)', texY = 'Number of Events / 30 GeV',
#      attribute = lambda event, sample: event.bj1['btagDeepB'] if event.bj1 is not None else float('nan'),
#      binning=[10,0,1],
#    ))
#
#   
#  plotting.fill(plots + plots2D, read_variables = read_variables, sequence = sequence)
#
#  # Get normalization yields from yield histogram
#  for plot in plots:
#    if plot.name == "yield":
#      for i, l in enumerate(plot.histos):
#        for j, h in enumerate(l):
#          yields[mode][plot.stack[i][j].name] = h.GetBinContent(h.FindBin(0.5+index))
#          h.GetXaxis().SetBinLabel(1, "#mu#mu")
#          h.GetXaxis().SetBinLabel(2, "e#mu")
#          h.GetXaxis().SetBinLabel(3, "ee")
#  if args.noData: yields[mode]["data"] = 0
#
#  yields[mode]["MC"] = sum(yields[mode][s.name] for s in mc_)
#  dataMCScale        = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')
#
#  drawPlots(plots + plots2D, mode, dataMCScale)
#  allPlots[mode] = plots + plots2D
#
## Add the different channels into SF and all
#for mode in ["SF","all"]:
#  yields[mode] = {}
#  for y in yields[allModes[0]]:
#    try:    yields[mode][y] = sum(yields[c][y] for c in (['ee','mumu'] if mode=="SF" else ['ee','mumu','mue']))
#    except: yields[mode][y] = 0
#  dataMCScale = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')
#
#  for plot in allPlots['mumu']:
#    for plot2 in (p for p in (allPlots['ee'] if mode=="SF" else allPlots["mue"]) if p.name == plot.name):  #For SF add EE, second round add EMu for all
#      for i, j in enumerate(list(itertools.chain.from_iterable(plot.histos))):
#	for k, l in enumerate(list(itertools.chain.from_iterable(plot2.histos))):
#	  if i==k:
#	    j.Add(l)
#
#  drawPlots(allPlots['mumu'], mode, dataMCScale)
#
#
#logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
#
