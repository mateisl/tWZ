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

################################################################################
# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',          action='store_true', help='Run only on a small subset of the data?', )
argParser.add_argument('--plot_directory', action='store', default='tWZ_GEN_v1')
argParser.add_argument('--era',            action='store', type=str, default="Run2016")
argParser.add_argument('--selection',      action='store', default='trilepT-minDLmass12-onZ1-njet4p-btag1')
args = argParser.parse_args()

################################################################################
# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"

logger.info( "Working in era %s", args.era)

################################################################################
# Define the MC samples
from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed import *

lumi_scale = 137.4

if args.era == "Run2016":
    mc = [Summer16.TWZ_NLO_DR]
    lumi_scale = 35.9
elif args.era == "Run2017":
    mc = [Fall17.TWZ_NLO_DR]
    lumi_scale = 41.5
elif args.era == "Run2018":
    mc = [Autumn18.TWZ_NLO_DR]
    lumi_scale = 60.0
elif args.era == "RunII":
    mc = [TWZ_NLO_DR]
for sample in mc:
    sample.scale           = 1 # Scale MCs individually with lumi


if args.small:
    for sample in mc:
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
def charge(pdgId):
    return -pdgId/abs(pdgId)

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

      if isinstance( plot, Plot):
          plotting.draw(plot,
            plot_directory = plot_directory_,
            ratio = None,
            logX = False, logY = log, sorting = True,
            yRange = (0.03, "auto") if log else (0.001, "auto"),
            scaling = {},
            legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
            drawObjects = drawObjects( False, dataMCScale , lumi_scale ) + _drawObjects,
            copyIndexPHP = True, extensions = ["png"],
          )


def getZDecay(event):
    Z1 = ROOT.TLorentzVector()
    Z2 = ROOT.TLorentzVector()
    NZdecays = 0
    for i in range(event.nGenPart):
        i_mother = event.GenPart_genPartIdxMother[i]
        if i_mother!=-1 and abs(event.GenPart_pdgId[i_mother]) == 23:
            # Go on if daughter is still a Z
            if abs(event.GenPart_pdgId[i]) == 23:
                continue
            else:
                NZdecays+=1
                if NZdecays == 1:
                    Z1.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_mass[i])
                elif NZdecays == 2:
                    Z2.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_mass[i])

    if NZdecays != 2:
        print "[ERROR] Found %s Z boson decay products" % ( NZdecays )
        return []
    else:
        if Z1.Pt() > Z2.Pt():
            return [Z1, Z2]
        else:
            return [Z2, Z1]

def getWorigin(index, event):
    while True:
        i_mother = event.GenPart_genPartIdxMother[index]
        if abs(event.GenPart_pdgId[i_mother]) == 24:
            index = i_mother
        else:
            return i_mother


def getWDecay(event):
    W1 = ROOT.TLorentzVector()
    W2 = ROOT.TLorentzVector()
    NWdecays = 0
    stati = []
    ids = []

    # Check which particle have a W as mother
    for i in range(event.nGenPart):
        i_mother = event.GenPart_genPartIdxMother[i]
        # Mother has to be a W
        if i_mother!=-1 and abs(event.GenPart_pdgId[i_mother]) == 24:
            # Continue if mother comes from top
            i_origin = getWorigin(i_mother, event)
            if abs(event.GenPart_pdgId[i_origin]) == 6:
                continue
            # Particle itself should not be a W or radiated Photon
            if abs(event.GenPart_pdgId[i]) == 24 or abs(event.GenPart_pdgId[i]) == 22:
                continue
            NWdecays+=1
            if NWdecays==1:
                W1.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_mass[i])
            elif NWdecays==2:
                W2.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_mass[i])


    if NWdecays != 2:
        print "[ERROR] Found %s W boson decay products" % ( NWdecays )
        return []
    else:
        if W1.Pt() > W2.Pt():
            return [W1, W2]
        else:
            return [W2, W1]

def getTopDecay(event):
    bquark = ROOT.TLorentzVector()
    Wdecay1 = ROOT.TLorentzVector()
    Wdecay2 = ROOT.TLorentzVector()
    Windex = -1
    Nb = 0
    NW = 0
    NWdaughter = 0
    # Get b and W
    for i in range(event.nGenPart):
        i_mother = event.GenPart_genPartIdxMother[i]
        if i_mother != -1 and abs(event.GenPart_pdgId[i_mother]) == 6:
            if abs(event.GenPart_pdgId[i]) == 5:
                Nb+=1
                bquark.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_mass[i])
            elif abs(event.GenPart_pdgId[i]) == 24:
                NW+=1
                Windex = i
    # Get W decay
    Wdecays = []
    for i in range(event.nGenPart):
        i_mother = event.GenPart_genPartIdxMother[i]
        if i_mother == Windex:
            if abs(event.GenPart_pdgId[i]) == 24: # If daughter is also a W, continue follow the decay
                Windex = i
                continue
            NWdaughter+=1
            if NWdaughter == 1:
                Wdecay1.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_mass[i])
            elif NWdaughter == 2:
                Wdecay2.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_mass[i])

    valid = True
    if Nb != 1:
        print "[ERROR] Found %s b quarks" % ( Nb )
        valid = False
    if NW != 1:
        print "[ERROR] Found %s W bosons" % ( NW )
        valid = False
    if NWdaughter != 2:
        print "[ERROR] Found %s W daughters" % ( NWdaughter )
        valid = False

    if Wdecay2.Pt() > Wdecay1.Pt():
        store = Wdecay1
        Wdecay1 = Wdecay2
        Wdecay2 = store

    if valid:
        return [bquark, Wdecay1, Wdecay2]
    else:
        return []

def getMaxDeltaR(particles):
    maxDR = 0
    for i in range(len(particles)):
        for j in range(len(particles)):
            if i!=j and particles[i].DeltaR(particles[j]) > maxDR:
                maxDR = particles[i].DeltaR(particles[j])
    return maxDR


################################################################################
# Define sequences
sequence       = []

def printGenParticles(event, sample):
    quarkIDs = [1,2,3,4,5] #exclude top
    leptonIDs = [11,12,13,14,15,16]
    statusvetos = [21,22,71]
    Nb = 0
    # print '----------------------'
    for i in range(event.nGenPart):
        status = event.GenPart_status[i]
        id = event.GenPart_pdgId[i]
        if status in statusvetos: continue
        if abs(id) in quarkIDs or abs(id) in leptonIDs:
            # print "Index: %s, ID: %s, Status %s" %(i, id, status)
            if abs(id) == 5: Nb+=1

    event.NBquarks = Nb
    event.test = [1,2,3]

sequence.append(printGenParticles)




def getZDecayProperties(event, sample):
    products = getZDecay(event)
    if len(products) < 2:
        event.ZDecay_Z1_pt = -1
        event.ZDecay_Z2_pt = -1
        event.ZDecay_DR = -1
    else:
        event.ZDecay_Z1_pt = products[0].Pt()
        event.ZDecay_Z2_pt = products[1].Pt()
        event.ZDecay_DR = products[0].DeltaR(products[1])

sequence.append(getZDecayProperties)

def getWDecayProperties(event, sample):
    products = getWDecay(event)
    if len(products) < 2:
        event.WDecay_W1_pt = -1
        event.WDecay_W2_pt = -1
        event.WDecay_DR = -1
    else:
        event.WDecay_W1_pt = products[0].Pt()
        event.WDecay_W2_pt = products[1].Pt()
        event.WDecay_DR = products[0].DeltaR(products[1])

sequence.append(getWDecayProperties)

def getTopDecayProperties(event, sample):
    products = getTopDecay(event)
    if len(products) < 3:
        event.TopDecay_b_pt = -1
        event.TopDecay_W1_pt = -1
        event.TopDecay_W2_pt = -1
        event.TopDecay_maxDR = -1
    else:
        event.TopDecay_b_pt = products[0].Pt()
        event.TopDecay_W1_pt = products[1].Pt()
        event.TopDecay_W2_pt = products[2].Pt()
        event.TopDecay_maxDR = getMaxDeltaR(products)

sequence.append(getTopDecayProperties)

def getPDGIDs(event, sample):
    ids = []
    for i in range(event.nGenPart):
        ids.append(event.GenPart_pdgId[i])
    event.ids = ids
sequence.append(getPDGIDs)

def getOrientations(event, sample):
    TopProducts = getTopDecay(event)
    WProducts = getWDecay(event)
    ZProducts = getZDecay(event)
    if len(TopProducts)<3 or len(WProducts)<2 or len(ZProducts)<2:
        event.DR_TopW = -1
        event.DR_TopZ = -1
        event.DR_WZ = -1
        event.mTop = -1
        event.mW = -1
        event.mZ = -1
    else:
        Top = TopProducts[0]+TopProducts[1]+TopProducts[2]
        W = WProducts[0]+WProducts[1]
        Z = ZProducts[0]+ZProducts[1]
        event.DR_TopW = Top.DeltaR(W)
        event.DR_TopZ = Top.DeltaR(Z)
        event.DR_WZ = W.DeltaR(Z)
        event.mTop = Top.M()
        event.mW = W.M()
        event.mZ = Z.M()
sequence.append(getOrientations)

################################################################################
# Read variables

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",
    "lep[pt/F,eta/F,phi/F,pdgId/I,muIndex/I,eleIndex/I]",
]

read_variables_MC = [
    'reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F',
    "genZ1_pt/F", "genZ1_eta/F", "genZ1_phi/F",
    "GenJet[pt/F,eta/F,phi/F,partonFlavour/I,hadronFlavour/I]",
    VectorTreeVariable.fromString( "GenPart[pt/F,mass/F,phi/F,eta/F,pdgId/I,genPartIdxMother/I,status/I,statusFlags/I]", nMax=1000),
    'nGenPart/I',
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
for i_mode, mode in enumerate(allModes):
    yields[mode] = {}

    weight_ = lambda event, sample: event.weight*lumi_year[event.year]/1000.

    for sample in mc: sample.style = styles.fillStyle(sample.color)

    for sample in mc:
      sample.read_variables = read_variables_MC
      sample.setSelectionString([getLeptonSelection(mode)])
      sample.weight = lambda event, sample: event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger#*event.reweightLeptonSF

    #yt_TWZ_filter.scale = lumi_scale * 1.07314

    stack = Stack(mc)

    # Use some defaults
    Plot.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = cutInterpreter.cutString(args.selection))

    ################################################################################
    # Now define the plots

    plots = []


    plots.append(Plot(
        name = "N_GenParts",
        texX = 'Number of generator particles', texY = 'Number of Events',
        attribute = lambda event, sample: event.nGenPart,
        binning=[50,0,200],
    ))

    plots.append(Plot(
        name = "N_Bquarks",
        texX = 'Number of generated b', texY = 'Number of Events',
        attribute = lambda event, sample: event.NBquarks,
        binning=[11,-0.5,10.5],
    ))

    plots.append(Plot(
        name = "TopDecay_b_pt",
        texX = 'p_{T}(b from top) GeV', texY = 'Number of Events',
        attribute = lambda event, sample: event.TopDecay_b_pt,
        binning=[50,0,400],
    ))

    plots.append(Plot(
        name = "TopDecay_W1_pt",
        texX = 'p_{T}(Leading W decay from top) GeV', texY = 'Number of Events',
        attribute = lambda event, sample: event.TopDecay_W1_pt,
        binning=[50,0,400],
    ))

    plots.append(Plot(
        name = "TopDecay_W2_pt",
        texX = 'p_{T}(Sub-leading W decay from top) GeV', texY = 'Number of Events',
        attribute = lambda event, sample: event.TopDecay_W2_pt,
        binning=[50,0,400],
    ))

    plots.append(Plot(
        name = "TopDecay_maxDR",
        texX = 'max #Delta R between top decays', texY = 'Number of Events',
        attribute = lambda event, sample: event.TopDecay_maxDR,
        binning=[50,0,10],
    ))

    plots.append(Plot(
        name = "m_top",
        texX = 'm_{top}', texY = 'Number of Events',
        attribute = lambda event, sample: event.mTop,
        binning=[50,0,200],
    ))

    plots.append(Plot(
        name = "ZDecay_Z1_pt",
        texX = 'p_{T}(Leading Z decay) GeV', texY = 'Number of Events',
        attribute = lambda event, sample: event.ZDecay_Z1_pt,
        binning=[50,0,400],
    ))

    plots.append(Plot(
        name = "ZDecay_Z2_pt",
        texX = 'p_{T}(Sub-leading Z decay) GeV', texY = 'Number of Events',
        attribute = lambda event, sample: event.ZDecay_Z2_pt,
        binning=[50,0,400],
    ))

    plots.append(Plot(
        name = "ZDecay_DR",
        texX = '#Delta R between Z decays', texY = 'Number of Events',
        attribute = lambda event, sample: event.ZDecay_DR,
        binning=[50,0,10],
    ))

    plots.append(Plot(
        name = "m_Z",
        texX = 'm_{Z}', texY = 'Number of Events',
        attribute = lambda event, sample: event.mZ,
        binning=[50,0,200],
    ))

    plots.append(Plot(
        name = "WDecay_W1_pt",
        texX = 'p_{T}(Leading W decay) GeV', texY = 'Number of Events',
        attribute = lambda event, sample: event.WDecay_W1_pt,
        binning=[50,0,400],
    ))

    plots.append(Plot(
        name = "WDecay_W2_pt",
        texX = 'p_{T}(Sub-leading W decay) GeV', texY = 'Number of Events',
        attribute = lambda event, sample: event.WDecay_W2_pt,
        binning=[50,0,400],
    ))

    plots.append(Plot(
        name = "WDecay_DR",
        texX = '#Delta R between W decays', texY = 'Number of Events',
        attribute = lambda event, sample: event.WDecay_DR,
        binning=[50,0,10],
    ))

    plots.append(Plot(
        name = "m_W",
        texX = 'm_{W}', texY = 'Number of Events',
        attribute = lambda event, sample: event.mW,
        binning=[50,0,200],
    ))

    plots.append(Plot(
        name = "DR_TopW",
        texX = '#Delta R(top,W)', texY = 'Number of Events',
        attribute = lambda event, sample: event.DR_TopW,
        binning=[50,0,10],
    ))

    plots.append(Plot(
        name = "DR_TopZ",
        texX = '#Delta R(top,Z)', texY = 'Number of Events',
        attribute = lambda event, sample: event.DR_TopZ,
        binning=[50,0,10],
    ))

    plots.append(Plot(
        name = "DR_WZ",
        texX = '#Delta R(W,Z)', texY = 'Number of Events',
        attribute = lambda event, sample: event.DR_WZ,
        binning=[50,0,10],
    ))

    plots.append(Plot(
        name = "PDGIDs",
        texX = 'pdg ID', texY = 'Number of Events',
        attribute = lambda event, sample: event.ids,
        binning=[51,-25,25],
    ))

    plotting.fill(plots, read_variables = read_variables, sequence = sequence)

    ################################################################################
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


    # yields[mode]["MC"] = sum(yields[mode][s.name] for s in mc)

    drawPlots(plots, mode, 1.0)

    allPlots[mode] = plots

    # This is how you would wirte a root file
    # outfile = ROOT.TFile('test.root', 'recreate')
    # outfile.cd()
    # for plot in plots:
    #     print "--------------------"
    #     print plot.name
    #     for i, l in enumerate(plot.histos):
    #         for j, h in enumerate(l):
    #             h.Write()
    # outfile.Close()
################################################################################
# Add the different channels into SF and all
for mode in ["comb1","comb2","all"]:
    yields[mode] = {}
    for y in yields[allModes[0]]:
        try:    yields[mode][y] = sum(yields[c][y] for c in ['eee','muee','mumue', 'mumumu'])
        except: yields[mode][y] = 0

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

    if mode == "all": drawPlots(allPlots['mumumu'], mode, 1.0)

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
