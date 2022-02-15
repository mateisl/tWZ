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
from tWZ.Tools.helpers                   import getCollection, cosThetaStarNew, getTheta, gettheta, getphi

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
argParser.add_argument('--small',          action='store_true', help='Run only on a small subset of the data?', )
argParser.add_argument('--dataMCScaling',  action='store_true', help='Data MC scaling?', )
argParser.add_argument('--plot_directory', action='store', default='EFT_GEN_v1')
argParser.add_argument('--era',            action='store', type=str, default="Run2018")
argParser.add_argument('--selection',      action='store', default='onZ1')
argParser.add_argument('--nicePlots',      action='store_true', default=False)
argParser.add_argument('--twoD',           action='store_true', default=False)
args = argParser.parse_args()

################################################################################
# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)


################################################################################
# Some info messages
if args.small:                        args.plot_directory += "_small"

logger.info( "Working in era %s", args.era)
if args.dataMCScaling:
    logger.info( "Data/MC scaling active")
else:
    logger.info( "Data/MC scaling not active")

if args.nicePlots:
    logger.info( "Only draw the plots")
else:
    logger.info( "Only saving into root file")

if args.twoD:
    logger.info( "Create EFT points in 2D")
else:
    logger.info( "Create EFT points in 1D")

################################################################################
# Define the MC samples
from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed import *
import tWZ.samples.nanoTuples_Autumn18_nanoAODv6_private_SMEFTsim_fast_postProcessed as SMEFTsim_fast

if args.era == "Run2016":
    mc = [Summer16.TWZ_NLO_DR, Summer16.TTZ, Summer16.TTX_rare, Summer16.TZQ, Summer16.WZ, Summer16.triBoson, Summer16.ZZ, Summer16.nonprompt_3l]
elif args.era == "Run2017":
    mc = [Fall17.TWZ_NLO_DR, Fall17.TTZ, Fall17.TTX_rare, Fall17.TZQ, Fall17.WZ, Fall17.triBoson, Fall17.ZZ, Fall17.nonprompt_3l]
elif args.era == "Run2018":
    mc = []
    samples_eft = [SMEFTsim_fast.tWZToLL01j_lepWFilter, SMEFTsim_fast.ttZ01j_lepWFilter]
    if args.nicePlots:
        mc = [SMEFTsim_fast.tWZToLL01j_lepWFilter, SMEFTsim_fast.ttZ01j_lepWFilter]
        samples_eft = []
        # mc = [SMEFTsim_fast.tWZToLL01j_lepWFilter]
elif args.era == "RunII":
    mc = [TWZ_NLO_DR, TTZ, TTX_rare, TZQ, WZ, triBoson, ZZ, nonprompt_3l]

################################################################################
# EFT reweight

# WeightInfo
eftweights = []
for sample in samples_eft:
    w = WeightInfo(sample.reweight_pkl)
    w.set_order(2)
    eftweights.append(w)

Npoints = 5

if args.nicePlots:
    Npoints = 0

WCs = []
WC_setup = [
    ('cHq1Re11', ROOT.kRed),
    ('cHq1Re22', ROOT.kGreen+2),
    ('cHq1Re33', ROOT.kOrange-3),
    ('cHq3Re11', ROOT.kCyan),
    ('cHq3Re22', ROOT.kMagenta),
    ('cHq3Re33', ROOT.kBlue),
    ('cHuRe11',  ROOT.kCyan),
    ('cHuRe22',  ROOT.kMagenta),
    ('cHuRe33',  ROOT.kBlue),
    ('cHdRe11',  ROOT.kViolet-9),
    ('cHdRe22',  ROOT.kGray),
    ('cHdRe33',  ROOT.kAzure+10),
    ('cW'     ,  ROOT.kBlue+4),
    ('cWtil'  ,  ROOT.kGreen-3),
]
for i_wc, (WCname, color) in enumerate(WC_setup):
    for i in range(Npoints):
        minval = -2.
        maxval = 2.
        # if 'cHq3Re11' in WCname:
        #     minval = -0.2
        #     maxval = 0.2
        value = minval + ((maxval-minval)/(Npoints-1))*i
        WCs.append( (WCname, value, color) )

params =  []
for i_sample, sample in enumerate(samples_eft):
    for i_wc, (WC, WCval, color) in enumerate(WCs):
        params.append({'legendText':'%s=%3.4f'%(WC, WCval), 'color':color,  'WC':{WC:WCval} , 'sample': sample, 'i_sample': i_sample})

#### 2D scan
if args.twoD:
    minval1  = -4.0
    maxval1  = 4.0
    minval2  = -4.0
    maxval2  = 4.0
    Npoints1 = 21
    Npoints2 = 21
    WC1  = 'cHq1Re1122'
    WC1a = 'cHq1Re11'
    WC1b = 'cHq1Re22'
    WC2  = 'cHq1Re33'
    params = []
    for i in range(Npoints1):
        value1 = minval1 + ((maxval1-minval1)/(Npoints1-1))*i
        for j in range(Npoints2):
            value2 = minval2 + ((maxval2-minval2)/(Npoints2-1))*j
            for i_sample, sample in enumerate(samples_eft):
                params.append({'legendText':'%s=%3.4f, %s=%3.4f'%(WC1,value1,WC2,value2), 'color':ROOT.kRed,  'WC':{WC1a:value1, WC1b:value1, WC2:value2} , 'sample': sample, 'i_sample': i_sample})
####

for i_param, param in enumerate(params):
    param['style']    = styles.lineStyle( param['color'] )

# Creating a list of weights
plotweights = []
weights_SM = []
# Add MC weights
weight_mc = []
for sample in mc:
    weight_ = lambda event, sample: 1. # Add event.weight and lumi weight to sample.weight later
    weight_mc.append(weight_)
plotweights.append(weight_mc)

# Add EFT weight
for param in params:
    i_sample = param['i_sample']
    eft_weight = eftweights[i_sample].get_weight_func(**param['WC'])
    plotweights.append([eft_weight])

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
    sample.scale           = 1

for param in params:
    param['sample'].scale = 1

if args.small:
    for sample in mc + [data_sample]:
        sample.normalization = 1.
        sample.reduceFiles( to = 1 )
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
    if "mt2ll100" in args.selection: lines += [(0.55, 0.5, 'M_{T2}(ll) > 100 GeV')] # Manually put the mt2ll > 100 GeV label
    return [tex.DrawLatex(*l) for l in lines]

def drawPlots(plots, mode, dataMCScale):
    for log in [False, True]:
        plot_directory_ = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, mode + ("_log" if log else ""), args.selection)
        for plot in plots:
            if not max(l.GetMaximum() for l in sum(plot.histos,[])): continue # Empty plot
            _drawObjects = []
            n_stacks=len(plot.histos)
            if isinstance( plot, Plot):
                plotting.draw(plot,
                  plot_directory = plot_directory_,
                  ratio = None,
                  logX = False, logY = log, sorting = True,
                  yRange = (0.03, "auto") if log else (0.001, "auto"),
                  scaling = {0:1} if args.dataMCScaling else {},
                  legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
                  drawObjects = drawObjects( False, dataMCScale , lumi_scale ) + _drawObjects,
                  copyIndexPHP = True, extensions = ["png", "pdf", "root"],
                )

def getFinalstates(event):
    W,Z,top,b = -1, -1, -1, -1
    # print "--------------------------"
    for i in range(event.nGenPart):
        # if abs(event.GenPart_pdgId[i]) == 5:
        #     print i, abs(event.GenPart_pdgId[i]), event.GenPart_genPartIdxMother[i], event.GenPart_status[i]
        if i > 10: break
        if event.GenPart_genPartIdxMother[i] in [0,1]:
            if abs(event.GenPart_pdgId[i]) == 6:
                top = i
            elif abs(event.GenPart_pdgId[i]) == 23:
                Z = i
            elif abs(event.GenPart_pdgId[i]) == 24:
                W = i
    return W,Z,top,b

def getbfromtop(event):
    for i in range(event.nGenPart):
        i_mother = event.GenPart_genPartIdxMother[i]
        if i_mother > 0:
            if abs(event.GenPart_pdgId[i]) == 5 and abs(event.GenPart_pdgId[i_mother]) == 6:
                return i
    return -1
################################################################################
# Define sequences
sequence       = []

def getBquarks(event,sample):
    W,Z,b,bTop = ROOT.TLorentzVector(),ROOT.TLorentzVector(),ROOT.TLorentzVector(),ROOT.TLorentzVector()
    i_W, i_Z, i_top, i_b = getFinalstates(event)
    foundW, foundZ, foundb, foundbtop = False, False, False, False
    if i_b > 0:
        b.SetPtEtaPhiM(event.GenPart_pt[i_b], event.GenPart_eta[i_b], event.GenPart_phi[i_b], event.GenPart_mass[i_b])
        foundb = True
    if i_top > 0:
        i_bTop = getbfromtop(event)
        if i_bTop > 0:
            bTop.SetPtEtaPhiM(event.GenPart_pt[i_bTop], event.GenPart_eta[i_bTop], event.GenPart_phi[i_bTop], event.GenPart_mass[i_bTop])
            foundbtop = True
    if i_W > 0:
        W.SetPtEtaPhiM(event.GenPart_pt[i_W], event.GenPart_eta[i_W], event.GenPart_phi[i_W], event.GenPart_mass[i_W])
        foundW = True
    if i_Z > 0:
        Z.SetPtEtaPhiM(event.GenPart_pt[i_Z], event.GenPart_eta[i_Z], event.GenPart_phi[i_Z], event.GenPart_mass[i_Z])
        foundZ = True
                
    event.b_pt = b.Pt() if foundb else float('nan')
    event.b_eta = b.Eta()  if foundb else float('nan')
    event.bTop_pt = bTop.Pt()  if foundbtop else float('nan')
    event.bTop_eta = bTop.Eta() if foundbtop else float('nan')
    event.b_diff_pt = bTop.Pt() - b.Pt() if (foundb and foundbtop) else float('nan')
    event.W_pt = W.Pt() if foundW else float('nan')
    event.W_eta = W.Eta() if foundW else float('nan')
    event.Z_eta = Z.Eta() if foundZ else float('nan')
    
    
    pass_b_veto = True
    if foundb:
        if event.b_pt > 30 or abs(event.b_eta)<2.5:
            pass_b_veto = False
            
    pass_b_tag = False
    if foundbtop:
        if event.bTop_pt > 30 and abs(event.bTop_eta)<2.5:
            pass_b_tag = True

    event.pass_b_veto = pass_b_veto
    event.pass_b_tag = pass_b_tag
sequence.append(getBquarks)


def getMWZ(event,sample):
    mWZ = -1
    i_W, i_Z, i_top, i_b = getFinalstates(event)
    W,Z = ROOT.TLorentzVector(),ROOT.TLorentzVector()
    if i_W > 0 and i_Z > 0:
        W.SetPtEtaPhiM(event.GenPart_pt[i_W], event.GenPart_eta[i_W], event.GenPart_phi[i_W], event.GenPart_mass[i_W])
        Z.SetPtEtaPhiM(event.GenPart_pt[i_Z], event.GenPart_eta[i_Z], event.GenPart_phi[i_Z], event.GenPart_mass[i_Z])
        WZ = W+Z
        mWZ = WZ.M()
    event.mWZ = mWZ
sequence.append(getMWZ)
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
    "weight/F", 'reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F',
    "genZ1_pt/F", "genZ1_eta/F", "genZ1_phi/F",
    "Muon[genPartFlav/I]",
    VectorTreeVariable.fromString( "GenPart[pt/F,mass/F,phi/F,eta/F,pdgId/I,genPartIdxMother/I,status/I,statusFlags/I]", nMax=1000),
    'nGenPart/I',
]

read_variables_eft = [
    "np/I", VectorTreeVariable.fromString("p[C/F]",nMax=200)
]


################################################################################
# Set up channels and values for plotting
yields     = {}
allPlots   = {}
allModes   = ['allchannels']

print "Working on channels:", allModes

for i_mode, mode in enumerate(allModes):
    yields[mode] = {}

    for sample in mc: sample.style = styles.fillStyle(sample.color)

    ###### Weights #################################################################
    weightnames = ['weight', 'reweightBTag_SF', 'reweightPU', 'reweightL1Prefire' , 'reweightTrigger'] # 'reweightLeptonSF'

    getters = map( operator.attrgetter, weightnames)
    def weight_function( event, sample):
        # Calculate weight, this becomes: w = event.weightnames[0]*event.weightnames[1]*...
        w = reduce(operator.mul, [g(event) for g in getters], 1)
        # Get Lumi weight and multiply to weight
        lumi_weight = lumi_year[event.year]/1000.
        w *= lumi_weight
        return w


    for sample in mc:
        sample.read_variables = read_variables_MC
        sample.weight = weight_function

    for param in params:
        param['sample'].read_variables = read_variables_MC + read_variables_eft
        param['sample'].weight = weight_function

    stack = Stack(mc, *[ [ param['sample'] ] for param in params ])
    noneftidxs = [0]
    if args.nicePlots:
        stack = Stack(mc)

    # Use some defaults
    selection_string = cutInterpreter.cutString(args.selection)
    Plot.setDefaults(stack = stack, weight = plotweights, selectionString = selection_string)

    ################################################################################
    # Now define the plots

    plots = []

    plots.append(Plot(
      name = 'yield', texX = '', texY = 'Number of Events',
      attribute = lambda event, sample: 0.5 + i_mode,
      binning=[4, 0, 4],
    ))

    plots.append(Plot(
        name = "gen_Z_eta",
        texX = 'Z #eta', texY = 'Number of Events',
        attribute = lambda event, sample: event.Z_eta,
        binning=[50, -5, 5],
    ))
    
    plots.append(Plot(
        name = "gen_W_eta",
        texX = 'W #eta', texY = 'Number of Events',
        attribute = lambda event, sample: event.W_eta,
        binning=[50, -5, 5],
    ))    
    
    plots.append(Plot(
        name = "gen_W_pt",
        texX = 'p_{T}(W) (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.W_pt,
        binning=[30, 0, 700],
    ))

    plots.append(Plot(
        name = "gen_Z1_pt",
        texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
        attribute = TreeVariable.fromString( "genZ1_pt/F" ),
        binning=[30, 0, 700],
    ))

    plots.append(Plot(
        name = "gen_m_WZ",
        texX = 'm_{WZ} (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mWZ,
        binning=[30, 0, 1400],
    ))

    plots.append(Plot(
        name = "gen_Z1_pt_bveto",
        texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.genZ1_pt if event.pass_b_veto else float('nan'),
        binning=[30, 0, 700],
    ))

    plots.append(Plot(
        name = "gen_m_WZ_bveto",
        texX = 'm_{WZ} (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mWZ  if event.pass_b_veto else float('nan'),
        binning=[30, 0, 1400],
    ))

    plots.append(Plot(
        name = "gen_Z1_pt_bveto_btag",
        texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.genZ1_pt if (event.pass_b_veto and event.pass_b_tag) else float('nan'),
        binning=[30, 0, 700],
    ))

    plots.append(Plot(
        name = "gen_m_WZ_bveto_btag",
        texX = 'm_{WZ} (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mWZ  if (event.pass_b_veto and event.pass_b_tag) else float('nan'),
        binning=[30, 0, 1400],
    ))

    plots.append(Plot(
        name = "gen_bveto",
        texX = 'pass b veto', texY = 'Number of Events',
        attribute = lambda event, sample: event.pass_b_veto,
        binning=[2, -0.5, 1.5],
    ))

    plots.append(Plot(
        name = "gen_btag",
        texX = 'pass b tag', texY = 'Number of Events',
        attribute = lambda event, sample: event.pass_b_tag,
        binning=[2, -0.5, 1.5],
    ))

    plots.append(Plot(
        name = "gen_b_pt",
        texX = 'b p_{T} (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.b_pt,
        binning=[50, 0, 500],
    ))

    plots.append(Plot(
        name = "gen_b_eta",
        texX = 'b #eta', texY = 'Number of Events',
        attribute = lambda event, sample: event.b_eta,
        binning=[50, -5, 5],
    ))

    plots.append(Plot(
        name = "gen_bTop_pt",
        texX = 'b_{top} p_{T} (GeV)', texY = 'Number of Events',
        attribute = lambda event, sample: event.bTop_pt,
        binning=[50, 0, 500],
    ))

    plots.append(Plot(
        name = "gen_bTop_eta",
        texX = 'b_{top} #eta', texY = 'Number of Events',
        attribute = lambda event, sample: event.bTop_eta,
        binning=[50, -5, 5],
    ))

    plots.append(Plot(
        name = "gen_b_diff_pt",
        texX = 'p_{T}(b_{top}) - p_{T}(b)', texY = 'Number of Events',
        attribute = lambda event, sample: event.b_diff_pt,
        binning=[50, -250, 250],
    ))

    plotting.fill(plots, read_variables = read_variables, sequence = sequence)


    ################################################################################
    # Get normalization yields from yield histogram
    for plot in plots:
        if plot.name == "yield":
            for i, l in enumerate(plot.histos):
                for j, h in enumerate(l):
                    yields[mode][plot.stack[i][j].name] = h.GetBinContent(h.FindBin(0.5+i_mode))
                    h.GetXaxis().SetBinLabel(1, "all")


    yields[mode]["data"] = 0

    yields[mode]["MC"] = sum(yields[mode][s.name] for s in mc  + samples_eft)
    dataMCScale        = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')

    if args.nicePlots: drawPlots(plots, mode, dataMCScale)

    allPlots[mode] = plots


################################################################################
# Add all different channels
yields["all"] = {}
for y in yields[allModes[0]]:
    try:    yields["all"][y] = sum(yields[c][y] for c in allModes)
    except: yields["all"][y] = 0
dataMCScale = yields["all"]["data"]/yields["all"]["MC"] if yields["all"]["MC"] != 0 else float('nan')

allPlots["all"] = allPlots[allModes[0]]
for plot in allPlots['all']:
    for i_mode,mode in enumerate(allModes):
        if i_mode == 0:
            continue
        tmp = allPlots[mode]
        for plot2 in (p for p in tmp if p.name == plot.name):
            for i, j in enumerate(list(itertools.chain.from_iterable(plot.histos))):
                for k, l in enumerate(list(itertools.chain.from_iterable(plot2.histos))):
                    if i==k:
                        j.Add(l)


if args.nicePlots:
    drawPlots(allPlots['all'], "all", dataMCScale)

if not args.nicePlots:
    plotstowrite = [
        "gen_Z1_pt",
        "gen_m_WZ",
        "gen_Z1_pt_bveto",
        "gen_m_WZ_bveto",
        "gen_Z1_pt_bveto_btag",
        "gen_m_WZ_bveto_btag",
    ]
    # Write Result Hist in root file
    plot_dir = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, "all", args.selection)
    if not os.path.exists(plot_dir):
        try:
            os.makedirs(plot_dir)
        except:
            print 'Could not crete', plot_dir
    outfilename = plot_dir+'/Results.root'
    if args.twoD: outfilename = plot_dir+'/Results_twoD.root'
    outfile = ROOT.TFile(outfilename, 'recreate')
    outfile.cd()
    for plot in allPlots['all']:
        if plot.name in plotstowrite:
            for idx, histo_list in enumerate(plot.histos):
                for j, h in enumerate(histo_list):
                    histname = h.GetName()
                    if "TWZ_NLO_DR" in histname: process = "tWZ"
                    elif "tWZToLL01j_lepWFilter" in histname: process = "tWZ"
                    elif "TTZ" in histname: process = "ttZ"
                    elif "ttZ01j" in histname: process = "ttZ"
                    elif "ttZ01j_lepWFilter" in histname: process = "ttZ"
                    elif "TTX_rare" in histname: process = "ttX"
                    elif "TZQ" in histname: process = "tZq"
                    elif "WZ" in histname: process = "WZ"
                    elif "ZZ" in histname: process = "ZZ"
                    elif "triBoson" in histname: process = "triBoson"
                    elif "nonprompt" in histname: process = "nonprompt"
                    elif "data" in histname: process = "data"
                    # Also add a string for the eft signal samples
                    n_noneft = len(noneftidxs)
                    if idx not in noneftidxs:
                        h.Write(plot.name+"__"+process+"__"+params[idx-n_noneft]['legendText'])
                    else:
                        h.Write(plot.name+"__"+process)
    outfile.Close()



logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
