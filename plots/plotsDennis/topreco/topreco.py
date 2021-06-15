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
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
argParser.add_argument('--plot_directory', action='store', default='tWZ_topreco_v1')
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

lumi_scale = 0

if args.era == "Run2016":
    mc = [Summer16.TWZ_NLO_DR, Summer16.TTZ]
    lumi_scale = 35.9
elif args.era == "Run2017":
    mc = [Fall17.TWZ_NLO_DR, Fall17.TTZ]
    lumi_scale = 41.5
elif args.era == "Run2018":
    mc = [Autumn18.TWZ_NLO_DR, Autumn18.TTZ]
    lumi_scale = 60.0
elif args.era == "RunII":
    mc = [TWZ_NLO_DR, TTZ]
    lumi_scale = 137.4

for sample in mc:
    sample.scale           = 1 # Scale MCs individually with lumi

if args.small:
    for sample in mc:
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
def charge(pdgId):
    return -pdgId/abs(pdgId)

def drawObjects( plotData, lumi_scale ):
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if plotData else 'CMS Simulation'),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) '% ( lumi_scale) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines]

def drawPlots(plots, mode):
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
            drawObjects = drawObjects( False, lumi_scale ) + _drawObjects,
            copyIndexPHP = True, extensions = ["png"],
          )


def getGenWs(event):
    indices = []
    for i in range(event.nGenPart):
        if abs(event.GenPart_pdgId[i]) == 24:
            indices.append(i)

    Whads = []
    Wleps = []
    quarkIDs = [1,2,3,4,5] #exclude top
    leptonIDs = [11,12,13,14,15,16]
    for i_mother in indices:
        W = ROOT.TLorentzVector()
        W.SetPtEtaPhiM(event.GenPart_pt[i_mother], event.GenPart_eta[i_mother], event.GenPart_phi[i_mother], event.GenPart_mass[i_mother])
        for i in range(event.nGenPart):
            if i_mother == event.GenPart_genPartIdxMother[i]:
                if abs(event.GenPart_pdgId[i]) in quarkIDs:
                    Whads.append(W)
                    break # prevent double counting
                elif abs(event.GenPart_pdgId[i]) in leptonIDs:
                    Wleps.append(W)
                    break # prevent double counting
    return [Wleps, Whads]

def getBJetindex( event ):
    maxscore = 0.0
    index = -1
    for i in range(event.nJetGood):
        btagscore = event.JetGood_btagDeepB[i]
        if btagscore > maxscore:
            maxscore = btagscore
            index = i
    return index


def getWhad( event ):
    min_dM = 1000.
    mW = 80.4
    Whad = ROOT.TLorentzVector()
    for i in range(event.nJetGood):
        jetindex1 = event.JetGood_index[i]
        jet1  = ROOT.TLorentzVector()
        jet1.SetPtEtaPhiM(event.Jet_pt[jetindex1], event.Jet_eta[jetindex1], event.Jet_phi[jetindex1], event.Jet_mass[jetindex1])
        for j in range(event.nJetGood):
            if j == i:
                continue
            jetindex2 = event.JetGood_index[j]
            jet2  = ROOT.TLorentzVector()
            jet2.SetPtEtaPhiM(event.Jet_pt[jetindex2], event.Jet_eta[jetindex2], event.Jet_phi[jetindex2], event.Jet_mass[jetindex2])

            dijet = jet1+jet2
            if abs(dijet.M()-mW) < min_dM:
                min_dM = abs(dijet.M()-mW)
                Whad = dijet

    return Whad

def getWlep( event ):
    Wlep = ROOT.TLorentzVector()
    lepton  = ROOT.TLorentzVector()
    met     = ROOT.TLorentzVector()
    lepton.SetPtEtaPhiM(event.lep_pt[event.nonZ1_l1_index], event.lep_eta[event.nonZ1_l1_index], event.lep_phi[event.nonZ1_l1_index], 0)
    met.SetPtEtaPhiM(event.met_pt, 0, event.met_phi, 0)

    lepton_pT = ROOT.TVector3(lepton.Px(), lepton.Py(), 0)
    neutrino_pT = ROOT.TVector3(met.Px(), met.Py(), 0)

    mass_w = 80.399
    mu = mass_w * mass_w / 2 + lepton_pT * neutrino_pT
    A = - (lepton_pT * lepton_pT)
    B = mu * lepton.Pz()
    C = mu * mu - lepton.E() * lepton.E() * (neutrino_pT * neutrino_pT)
    discriminant = B * B - A * C
    neutrinos = []
    if discriminant <= 0:
        # Take only real part of the solution for pz:
        neutrino = ROOT.TLorentzVector()
        neutrino.SetPxPyPzE(met.Px(),met.Py(),-B / A,0)
        neutrino.SetE(neutrino.P())
        neutrinos.append(neutrino)
    else:
        discriminant = sqrt(discriminant)
        neutrino1 = ROOT.TLorentzVector()
        neutrino1.SetPxPyPzE(met.Px(),met.Py(),(-B - discriminant) / A,0)
        neutrino1.SetE(neutrino1.P())
        neutrino2 = ROOT.TLorentzVector()
        neutrino2.SetPxPyPzE(met.Px(),met.Py(),(-B + discriminant) / A,0)
        neutrino2.SetE(neutrino2.P())
        if neutrino1.E() > neutrino2.E():
            neutrinos.append(neutrino1)
            neutrinos.append(neutrino2)
        else:
            neutrinos.append(neutrino2)
            neutrinos.append(neutrino1)

    Wleps = []
    for neu in neutrinos:
        Wlep = lepton + neu
        Wleps.append([Wlep, lepton, neu])
    return Wleps

def getGenTop_indices(event, mode):
    top_indices = []
    for i in range(event.nGenPart):
        # Get tops after production
        if mode == "production":
            if abs(event.GenPart_pdgId[i]) == 6:
                i_mother = event.GenPart_genPartIdxMother[i]
                if i_mother<0 or abs(event.GenPart_pdgId[i_mother]) != 6:
                    top_indices.append(i)
        # Get tops before decay (as mother of particle)
        elif mode == "decay":
            i_mother = event.GenPart_genPartIdxMother[i]
            if i_mother > 0 and abs(event.GenPart_pdgId[i_mother]) == 6:
                if abs(event.GenPart_pdgId[i]) != 6: # particle itself is not a top
                    if i_mother not in top_indices:
                        top_indices.append(i_mother)
        else:
            print '[ERROR] No top mode given! Unable to find top indices'
    return top_indices

def getGenBottomFromTop_indices(event):
    bot_indices = []
    for i in range(event.nGenPart):
        if abs(event.GenPart_pdgId[i]) == 5:
            i_mother = event.GenPart_genPartIdxMother[i]
            if i_mother<0 or abs(event.GenPart_pdgId[i_mother]) != 5:
                if i_mother>0 and abs(event.GenPart_pdgId[i_mother])==6:
                    bot_indices.append(i)
    return bot_indices

def getGenW_indices(event):
    W_indices = []
    for i in range(event.nGenPart):
        if abs(event.GenPart_pdgId[i]) == 24:
            i_mother = event.GenPart_genPartIdxMother[i]
            if i_mother<0 or abs(event.GenPart_pdgId[i_mother]) != 24:
                    W_indices.append(i)
    return W_indices

def getGenWDecay_indices(event):
    Wdecay_indices = []
    Wdecay_mother_indices = []
    NWdecays = 0
    # Check which particle have a W as mother
    for i in range(event.nGenPart):
        i_mother = event.GenPart_genPartIdxMother[i]
        # Mother has to be a W
        if i_mother!=-1 and abs(event.GenPart_pdgId[i_mother]) == 24:
            # Particle itself should not be a W or radiated Photon
            if abs(event.GenPart_pdgId[i]) == 24 or abs(event.GenPart_pdgId[i]) == 22:
                continue
            Wdecay_indices.append(i)

    for i in Wdecay_indices:
        Wdecay_mother_indices.append(event.GenPart_genPartIdxMother[i])

    Wdecay_indices_sorted = []
    used_idx = []
    for i in range(len(Wdecay_indices)):
        if i in used_idx: continue
        for j in range(len(Wdecay_indices)):
            if i==j or j in used_idx:
                continue
            if Wdecay_mother_indices[i] == Wdecay_mother_indices[j]:
                Wdecay_indices_sorted.append( [ Wdecay_indices[i],Wdecay_indices[j] ] )
                used_idx.append(i)
                used_idx.append(j)
    return Wdecay_indices_sorted
    # returns [ [W1decay1, W1decay2], [W2decay1, W2decay2], ... ]

def getGenWhadDecays(event):
    quarkIDs = [1,2,3,4,5] #exclude top
    GenWDecay_idx = getGenWDecay_indices(event)

    WDecay1_had = ROOT.TLorentzVector()
    WDecay2_had = ROOT.TLorentzVector()

    decay_pairs = []

    for idx in GenWDecay_idx:
        id = event.GenPart_pdgId[idx[0]]
        if abs(id) in quarkIDs:
            WDecay1_had.SetPtEtaPhiM(event.GenPart_pt[idx[0]], event.GenPart_eta[idx[0]], event.GenPart_phi[idx[0]], event.GenPart_mass[idx[0]])
            WDecay2_had.SetPtEtaPhiM(event.GenPart_pt[idx[1]], event.GenPart_eta[idx[1]], event.GenPart_phi[idx[1]], event.GenPart_mass[idx[1]])
            decay_pairs.append([WDecay1_had, WDecay2_had])

    return decay_pairs

def findMotherTop(event,idx):
    topidx = -1
    motherindex = event.GenPart_genPartIdxMother[idx]
    if motherindex<0: return -1
    # If particle is a b and comes from top, return top idx
    if abs(event.GenPart_pdgId[idx]) == 5 and abs(event.GenPart_pdgId[motherindex]) == 6:
        return motherindex
    # If particle comes from a W but is not a W, follow W
    if abs(event.GenPart_pdgId[motherindex]) == 24 and abs(event.GenPart_pdgId[idx]) != 24:
        foundtop = False
        Wmotherindex = event.GenPart_genPartIdxMother[motherindex]
        while(not foundtop):
            # If no top is found, leave loop
            if Wmotherindex < 0:
                foundtop = True
            # If top is found, return top idx
            elif abs(event.GenPart_pdgId[Wmotherindex]) == 6:
                return Wmotherindex
            # Else, go to mother and repeat
            else:
                Wmotherindex = event.GenPart_genPartIdxMother[Wmotherindex]
    return -1

def getTopDecays(event):
    alltops = getGenTop_indices(event, "decay")
    alldecays = [] # [ [topidx, decayindx], [], ...]
    for i in range(event.nGenPart):
        topidx = findMotherTop(event, i)
        if topidx > 0:
            alldecays.append([topidx, i])

    decays_sorted = []
    for topindex in alltops:
        decaysonetop = []
        for topidx, decayidx in alldecays:
            if topindex == topidx:
                decaysonetop.append(decayidx)
        decays_sorted.append(decaysonetop)
    return decays_sorted

def getMatchingJets(event, parton_idxs):
    jets = []
    selectedjet_idxs = []
    for i in parton_idxs:
        parton = ROOT.TLorentzVector()
        parton.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_mass[i])
        min_dR = 0.4
        found_match = False
        selectedjet = ROOT.TLorentzVector()
        for j in range(event.nJetGood):
            if j in selectedjet_idxs:
                continue
            jetidx = event.JetGood_index[j]
            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(event.Jet_pt[jetidx], event.Jet_eta[jetidx], event.Jet_phi[jetidx], event.Jet_mass[jetidx])
            if parton.DeltaR(jet) < min_dR:
                min_dR = parton.DeltaR(jet)
                selectedjet = jet
                selectedjet_idx = j
                found_match = True
        if found_match:
            jets.append(selectedjet)
            selectedjet_idxs.append(selectedjet_idx)
    return jets

def isHadronicDecay(event, idxs):
    quarkIDs = [1,2,3,4,5]
    leptonIDs = [11,12,13,14,15,16]
    isHadronic = True
    for i in idxs:
        if abs(event.GenPart_pdgId[i]) not in quarkIDs:
            isHadronic = False
    return isHadronic


def calculate_chi2(toplep, tophad, Whad, mode):
    Mtlep_mean   = 174.
    Mtlep_sigma  =  18.
    Mthad_mean   = 181.
    Mthad_sigma  =  15.
    MWhad_mean   =  80.399
    MWhad_sigma  =  12.
    Mtdiff_sigma = sqrt(pow(Mtlep_sigma,2)+pow(Mthad_sigma,2))
    toplepterm  = pow((toplep.M()-Mtlep_mean)/Mtlep_sigma,2)
    tophadterm  = pow((tophad.M()-Mthad_mean)/Mthad_sigma,2)
    Whadterm    = pow((  Whad.M()-MWhad_mean)/MWhad_sigma,2)
    topdiffterm = pow((toplep.M()-tophad.M())/Mtdiff_sigma,2)
    if mode=="topdiff":
        chi2 = topdiffterm+Whadterm
    else:
        chi2 = toplepterm+tophadterm+Whadterm
    return chi2

def getTopHypos(event, Njetsmax):
    # Set maximum number of jets
    if event.nJetGood < Njetsmax:
        Njetsmax = event.nJetGood

    # Find all possible solutions for the 4 missing jets
    # (blep,bhad,Wdecay1, Wdecay2)
    jet_permutations = []
    for i in range(Njetsmax):
        for j in range(Njetsmax):
            if j == i:
                continue
            for k in range(Njetsmax):
                if k==i or k==j:
                    continue
                for l in range(Njetsmax):
                    if l==i or l==j or l==k:
                        continue
                    jet_i = ROOT.TLorentzVector()
                    jet_j = ROOT.TLorentzVector()
                    jet_k = ROOT.TLorentzVector()
                    jet_l = ROOT.TLorentzVector()
                    jetidx_i = event.JetGood_index[i]
                    jetidx_j = event.JetGood_index[j]
                    jetidx_k = event.JetGood_index[k]
                    jetidx_l = event.JetGood_index[l]
                    jet_i.SetPtEtaPhiM(event.Jet_pt[jetidx_i], event.Jet_eta[jetidx_i], event.Jet_phi[jetidx_i], event.Jet_mass[jetidx_i])
                    jet_j.SetPtEtaPhiM(event.Jet_pt[jetidx_j], event.Jet_eta[jetidx_j], event.Jet_phi[jetidx_j], event.Jet_mass[jetidx_j])
                    jet_k.SetPtEtaPhiM(event.Jet_pt[jetidx_k], event.Jet_eta[jetidx_k], event.Jet_phi[jetidx_k], event.Jet_mass[jetidx_k])
                    jet_l.SetPtEtaPhiM(event.Jet_pt[jetidx_l], event.Jet_eta[jetidx_l], event.Jet_phi[jetidx_l], event.Jet_mass[jetidx_l])
                    jet_permutations.append([jet_i, jet_j, jet_k, jet_l])

    # Get Wlep from lepton + MET
    Wleps   = getWlep(event)
    # build hypotheses
    hypotheses = []
    for Wlep, lepton, neutrino in Wleps:
        for permutation in jet_permutations:
            hypo = {}
            hypo['toplep'] = Wlep + permutation[0]
            hypo['Wlep']   = Wlep
            hypo['tophad'] = permutation[1] + permutation[2] + permutation[3]
            hypo['Whad']   = permutation[2] + permutation[3]
            hypo['chi2']   = calculate_chi2(hypo['toplep'], hypo['tophad'], hypo['Whad'], "normal")
            hypo['chi2_topdiff']   = calculate_chi2(hypo['toplep'], hypo['tophad'], hypo['Whad'], "topdiff")
            hypotheses.append(hypo)

    return hypotheses

################################################################################
# Define sequences
sequence       = []

def MatchWs(event, sample):
    Whads_gen = []
    Wleps_gen = []
    Wleps_gen, Whads_gen = getGenWs(event)
    event.NWlep = len(Wleps_gen)
    event.NWhad = len(Whads_gen)

    Wlep_reco = getWlep(event)[0][0]
    Whad_reco = getWhad(event)

    if len(Wleps_gen)!=1 or len(Whads_gen)!=1:
        return

    event.dR_Wlep_gen_reco = Wlep_reco.DeltaR(Wleps_gen[0])
    event.dR_Whad_gen_reco = Whad_reco.DeltaR(Whads_gen[0])

    event.mWlep = Wlep_reco.M()
    event.mWhad = Whad_reco.M()
sequence.append(MatchWs)

def MatchWlepDecay(event, sample):
    dR_lepton = -1
    dR_neutrino = -1
    lepton_reco = getWlep(event)[0][1]
    neutrino_reco = getWlep(event)[0][2]
    W_idxs = getGenWDecay_indices(event)
    Wlep_idxs = []
    for Wdecay1,Wdecay2 in W_idxs:
        if not isHadronicDecay(event, [Wdecay1,Wdecay2]):
            Wlep_idxs.append([Wdecay1,Wdecay2])
    event.NWlep_matching = len(Wlep_idxs)
    if len(Wlep_idxs) != 1:
        return
    else:
        lepton_gen = ROOT.TLorentzVector()
        neutrino_gen = ROOT.TLorentzVector()
        if abs(event.GenPart_pdgId[Wdecay1]) in [11,13,15]:
            lepton_gen.SetPtEtaPhiM(event.GenPart_pt[Wdecay1], event.GenPart_eta[Wdecay1], event.GenPart_phi[Wdecay1], event.GenPart_mass[Wdecay1])
            neutrino_gen.SetPtEtaPhiM(event.GenPart_pt[Wdecay2], event.GenPart_eta[Wdecay2], event.GenPart_phi[Wdecay2], event.GenPart_mass[Wdecay2])
        else:
            lepton_gen.SetPtEtaPhiM(event.GenPart_pt[Wdecay2], event.GenPart_eta[Wdecay2], event.GenPart_phi[Wdecay2], event.GenPart_mass[Wdecay2])
            neutrino_gen.SetPtEtaPhiM(event.GenPart_pt[Wdecay1], event.GenPart_eta[Wdecay1], event.GenPart_phi[Wdecay1], event.GenPart_mass[Wdecay1])
        dR_lepton = lepton_reco.DeltaR(lepton_gen)
        dR_neutrino = neutrino_reco.DeltaR(neutrino_gen)
    event.dR_lepton_gen_reco = dR_lepton
    event.dR_neutrino_gen_reco = dR_neutrino
sequence.append(MatchWlepDecay)

def BestWhadReco(event, sample):
    decays = getGenWhadDecays(event)
    if len(decays) != 1:
        event.mWhad_matching = float('nan')
        event.maxDR_Whaddecays = float('nan')
        return

    WDecay1_had = decays[0][0]
    WDecay2_had = decays[0][1]

    min_dR = 1000
    jet1 = ROOT.TLorentzVector()
    i_jet1 = -1
    for i in range(event.nJetGood):
        jetindex = event.JetGood_index[i]
        jet_tmp  = ROOT.TLorentzVector()
        jet_tmp.SetPtEtaPhiM(event.Jet_pt[jetindex], event.Jet_eta[jetindex], event.Jet_phi[jetindex], event.Jet_mass[jetindex])
        if WDecay1_had.DeltaR(jet_tmp) < min_dR:
            min_dR = WDecay1_had.DeltaR(jet_tmp)
            jet1 = jet_tmp
            i_jet1 = i

    min_dR = 1000 # reset min_dR
    jet2 = ROOT.TLorentzVector()
    for i in range(event.nJetGood):
        if i == i_jet1:
            continue
        jetindex = event.JetGood_index[i]
        jet_tmp  = ROOT.TLorentzVector()
        jet_tmp.SetPtEtaPhiM(event.Jet_pt[jetindex], event.Jet_eta[jetindex], event.Jet_phi[jetindex], event.Jet_mass[jetindex])
        if WDecay2_had.DeltaR(jet_tmp) < min_dR:
            min_dR = WDecay2_had.DeltaR(jet_tmp)
            jet2 = jet_tmp

    Wreco = jet1+jet2
    event.mWhad_matching = Wreco.M()
    event.maxDR_Whaddecays = WDecay1_had.DeltaR(jet1) if WDecay1_had.DeltaR(jet1) > WDecay2_had.DeltaR(jet2) else WDecay2_had.DeltaR(jet2)
sequence.append(BestWhadReco)

def Matching(event, sample):
    top_idxs = getGenTop_indices(event, "production")
    bot_idxs = getGenBottomFromTop_indices(event)
    event.Ntops = len(top_idxs)
    event.Nbots = len(bot_idxs)
    event.mtophad_matched = -1

    topdecay_idxs = getTopDecays(event)
    Nt=0
    Nthad=0
    Ntlep=0

    numberoftopdecays = []
    for topdecays in topdecay_idxs:
        numberoftopdecays.append(len(topdecays))
        if len(topdecays) == 3:
            Nt+=1
            if isHadronicDecay(event, topdecays):
                Nthad+=1
                matched_jets = getMatchingJets(event, topdecays)
                event.Nmatched_tophad = len(matched_jets)
                if len(matched_jets) == 3:
                    top = matched_jets[0]+matched_jets[1]+matched_jets[2]
                    event.mtophad_matched = top.M()
            else:
                Ntlep+=1
    event.Ntops_matching = Nt
    event.Ntopshad_matching = Nthad
    event.Ntopslep_matching = Ntlep
    event.Ntopdecays = numberoftopdecays
sequence.append(Matching)

def TopReco(event, sample):
    hypotheses=getTopHypos(event, 5)
    chi2min = 10000
    chi2min_topdiff = 10000
    hypo_selected = hypotheses[0]
    hypo_selected_topdiff = hypotheses[0]
    foundhypo = False
    foundhypo_topdiff = False
    for hypo in hypotheses:
        if hypo['chi2']<chi2min:
            chi2min = hypo['chi2']
            hypo_selected = hypo
            foundhypo = True
        if hypo['chi2_topdiff']<chi2min:
            chi2min_topdiff = hypo['chi2_topdiff']
            hypo_selected_topdiff = hypo
            foundhypo_topdiff = True

    event.mtop_average = (hypo_selected['toplep'].M()+hypo_selected['tophad'].M())/2 if foundhypo else -1
    event.chi2 = hypo_selected['chi2'] if foundhypo else -1
    event.pgof = exp(-0.5*hypo_selected['chi2']) if foundhypo else -1
    event.hypofound = 1 if foundhypo else 0
    event.mtop_average_topdiff = (hypo_selected_topdiff['toplep'].M()+hypo_selected_topdiff['tophad'].M())/2 if foundhypo_topdiff else -1
    event.chi2_topdiff = hypo_selected_topdiff['chi2_topdiff'] if foundhypo_topdiff else -1
    event.pgof_topdiff = exp(-0.5*hypo_selected_topdiff['chi2']) if foundhypo_topdiff else -1
    event.hypofound_topdiff = 1 if foundhypo_topdiff else 0


    # Match with Gen Level
    top_idxs = getGenTop_indices(event, "decay")
    tops = []
    for idx in top_idxs:
        top = ROOT.TLorentzVector()
        top.SetPtEtaPhiM(event.GenPart_pt[idx], event.GenPart_eta[idx], event.GenPart_phi[idx], event.GenPart_mass[idx])
        tops.append(top)

    if len(tops) == 2:
        dR_lep_1 = hypo_selected['toplep'].DeltaR(tops[0])
        dR_lep_2 = hypo_selected['toplep'].DeltaR(tops[1])
        dR_had_1 = hypo_selected['tophad'].DeltaR(tops[0])
        dR_had_2 = hypo_selected['tophad'].DeltaR(tops[1])
        if dR_lep_1+dR_had_2 < dR_lep_2+dR_had_1:
            toplep_gen = tops[0]
            tophad_gen = tops[1]
        else:
            toplep_gen = tops[1]
            tophad_gen = tops[0]
        event.dR_toplep_gen_reco = hypo_selected['toplep'].DeltaR(toplep_gen) if foundhypo else -1
        event.dR_tophad_gen_reco = hypo_selected['tophad'].DeltaR(tophad_gen) if foundhypo else -1
        event.reso_mass_toplep = (hypo_selected['toplep'].M()-toplep_gen.M())/toplep_gen.M() if foundhypo else float('nan')
        event.reso_mass_tophad = (hypo_selected['tophad'].M()-tophad_gen.M())/tophad_gen.M() if foundhypo else float('nan')
        event.reso_pt_toplep = (hypo_selected['toplep'].Pt()-toplep_gen.Pt())/toplep_gen.Pt() if foundhypo else float('nan')
        event.reso_pt_tophad = (hypo_selected['tophad'].Pt()-tophad_gen.Pt())/tophad_gen.Pt() if foundhypo else float('nan')
    else:
        event.dR_toplep_gen_reco = float('nan')
        event.dR_tophad_gen_reco = float('nan')
        event.reso_mass_toplep = float('nan')
        event.reso_mass_tophad = float('nan')
        event.reso_pt_toplep = float('nan')
        event.reso_pt_tophad = float('nan')


sequence.append(TopReco)




def Zmother(event, sample):
    for i in range(event.nGenPart):
        if abs(event.GenPart_pdgId[i]) == 23:
            i_mother = event.GenPart_genPartIdxMother[i]
            pdgid_mother = event.GenPart_pdgId[i_mother]
            if abs(pdgid_mother) == 23:
                continue
            else:
                event.Zmother = abs(pdgid_mother)
sequence.append(Zmother)

################################################################################
# Read variables

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",  "nJet/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I",
    "l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
    "l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
    "JetGood[pt/F,eta/F,phi/F,area/F,btagDeepB/F,index/I]",
    "Jet[pt/F,eta/F,phi/F,mass/F]",
    "lep[pt/F,eta/F,phi/F,pdgId/I,muIndex/I,eleIndex/I]",
    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I",
    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
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
allPlots   = {}
allModes   = ['mumumu','mumue','muee', 'eee']
for i_mode, mode in enumerate(allModes):

    weight_ = lambda event, sample: event.weight if sample.isData else event.weight*lumi_year[event.year]/1000.

    for sample in mc: sample.style = styles.fillStyle(sample.color)

    for sample in mc:
      sample.read_variables = read_variables_MC
      sample.setSelectionString([getLeptonSelection(mode)])
      sample.weight = lambda event, sample: event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger#*event.reweightLeptonSF

    stack = Stack(mc)

    # Use some defaults
    Plot.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = cutInterpreter.cutString(args.selection))

    ################################################################################
    # Now define the plots

    plots = []


    plots.append(Plot(
        name = 'N_tops',
        texX = 'number of generated tops', texY = 'Number of Events',
        attribute = lambda event, sample: event.Ntops,
        binning=[11,-0.5,10.5],
    ))

    plots.append(Plot(
        name = 'N_tops_matching',
        texX = 'number of generated tops from decay products', texY = 'Number of Events',
        attribute = lambda event, sample: event.Ntops_matching,
        binning=[11,-0.5,10.5],
    ))

    plots.append(Plot(
        name = 'N_tops_had_matching',
        texX = 'number of generated had. tops from decay products', texY = 'Number of Events',
        attribute = lambda event, sample: event.Ntopshad_matching,
        binning=[11,-0.5,10.5],
    ))

    plots.append(Plot(
        name = 'N_tops_lep_matching',
        texX = 'number of generated lep. tops from decay products', texY = 'Number of Events',
        attribute = lambda event, sample: event.Ntopslep_matching,
        binning=[11,-0.5,10.5],
    ))

    plots.append(Plot(
        name = 'N_W_lep_matching',
        texX = 'number of generated lep. W from decay products', texY = 'Number of Events',
        attribute = lambda event, sample: event.NWlep_matching,
        binning=[11,-0.5,10.5],
    ))

    plots.append(Plot(
        name = 'N_jets_mached_tophad',
        texX = 'number of jets matched to had. top', texY = 'Number of Events',
        attribute = lambda event, sample: event.Nmatched_tophad,
        binning=[11,-0.5,10.5],
    ))

    plots.append(Plot(
        name = 'm_top_matched',
        texX = 'm(top_{had}) after matching', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtophad_matched,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'N_top_decays',
        texX = 'number of generated top decay products', texY = 'Number of Events',
        attribute = lambda event, sample: event.Ntopdecays,
        binning=[11,-0.5,10.5],
    ))


    plots.append(Plot(
        name = 'N_bottoms',
        texX = 'number of generated bottoms', texY = 'Number of Events',
        attribute = lambda event, sample: event.Nbots,
        binning=[11,-0.5,10.5],
    ))

    plots.append(Plot(
        name = 'N_Wlep',
        texX = 'number of generated leptonic Ws', texY = 'Number of Events',
        attribute = lambda event, sample: event.NWlep,
        binning=[11,-0.5,10.5],
    ))

    plots.append(Plot(
        name = 'N_Wlep',
        texX = 'number of generated leptonic Ws', texY = 'Number of Events',
        attribute = lambda event, sample: event.NWlep,
        binning=[11,-0.5,10.5],
    ))

    plots.append(Plot(
        name = 'N_Whad',
        texX = 'number of generated hadronic Ws', texY = 'Number of Events',
        attribute = lambda event, sample: event.NWhad,
        binning=[11,-0.5,10.5],
    ))

    plots.append(Plot(
        name = 'm_Wlep',
        texX = 'm(W_{lep})', texY = 'Number of Events',
        attribute = lambda event, sample: event.mWlep,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_Whad',
        texX = 'm(W_{had})', texY = 'Number of Events',
        attribute = lambda event, sample: event.mWhad,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'dR_lepton_gen_reco',
        texX = '#Delta R(#ell^{reco}, #ell^{gen})', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_lepton_gen_reco,
        binning=[35, 0.0, 7],
    ))

    plots.append(Plot(
        name = 'dR_neutrino_gen_reco',
        texX = '#Delta R(#nu^{reco}, #nu^{gen})', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_neutrino_gen_reco,
        binning=[35, 0.0, 7],
    ))

    plots.append(Plot(
        name = 'dR_Wlep_gen_reco',
        texX = '#Delta R(W_{lep}^{reco}, W_{lep}^{gen})', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_Wlep_gen_reco,
        binning=[35, 0.0, 7],
    ))

    plots.append(Plot(
        name = 'dR_Whad_gen_reco',
        texX = '#Delta R(W_{had}^{reco}, W_{had}^{gen})', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_Whad_gen_reco,
        binning=[35, 0.0, 7],
    ))

    plots.append(Plot(
        name = 'maxDR_Whaddecays',
        texX = 'max #Delta R(W_{had decay}^{reco}, W_{had decay}^{gen}) after matching', texY = 'Number of Events',
        attribute = lambda event, sample: event.maxDR_Whaddecays,
        binning=[35, 0.0, 7],
    ))

    plots.append(Plot(
        name = 'mWhad_matching',
        texX = 'm(W_{had}) after matching', texY = 'Number of Events',
        attribute = lambda event, sample: event.mWhad_matching,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = "Zmother",
        texX = 'pdg ID of Z mother', texY = 'Number of Events',
        attribute = lambda event, sample: event.Zmother,
        binning=[26,-0.5, 25.5],
    ))


    plots.append(Plot(
        name = 'm_top_average',
        texX = '(m_{t1}+m_{t2})/2', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_average,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_average_topdiff',
        texX = '(m_{t1}+m_{t2})/2 (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_average_topdiff,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'chi2',
        texX = '#chi^{2}', texY = 'Number of Events',
        attribute = lambda event, sample: event.chi2,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'p_gof',
        texX = 'P_{gof}', texY = 'Number of Events',
        attribute = lambda event, sample: event.pgof,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'chi2_topdiff',
        texX = '#chi^{2} (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.chi2_topdiff,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'p_gof_topdiff',
        texX = 'P_{gof} (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.pgof_topdiff,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'hypofound',
        texX = 'found ttbar hypothesis', texY = 'Number of Events',
        attribute = lambda event, sample: event.hypofound,
        binning=[2, -0.5, 1.5],
    ))

    plots.append(Plot(
        name = 'hypofound_topdiff',
        texX = 'found ttbar hypothesis (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.hypofound_topdiff,
        binning=[2, -0.5, 1.5],
    ))

    plots.append(Plot(
        name = 'dR_toplep_gen_reco',
        texX = '#Delta R(top_{lep}^{reco},top_{lep}^{gen})', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_toplep_gen_reco,
        binning=[35, 0.0, 7],
    ))

    plots.append(Plot(
        name = 'reso_mass_toplep',
        texX = '[M(top_{lep}^{reco})-M(top_{lep}^{gen})]/M(top_{lep}^{gen})', texY = 'Number of Events',
        attribute = lambda event, sample: event.reso_mass_toplep,
        binning=[40, -2.0, 2.0],
    ))

    plots.append(Plot(
        name = 'reso_pt_toplep',
        texX = '[p_{T}(top_{lep}^{reco})-p_{T}(top_{lep}^{gen})]/p_{T}(top_{lep}^{gen})', texY = 'Number of Events',
        attribute = lambda event, sample: event.reso_pt_toplep,
        binning=[40, -2.0, 2.0],
    ))

    plots.append(Plot(
        name = 'dR_tophad_gen_reco',
        texX = '#Delta R(top_{had}^{reco},top_{had}^{gen})', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_tophad_gen_reco,
        binning=[35, 0.0, 7],
    ))

    plots.append(Plot(
        name = 'reso_mass_tophad',
        texX = '[M(top_{had}^{reco})-M(top_{had}^{gen})]/M(top_{had}^{gen})', texY = 'Number of Events',
        attribute = lambda event, sample: event.reso_mass_tophad,
        binning=[40, -2.0, 2.0],
    ))

    plots.append(Plot(
        name = 'reso_pt_tophad',
        texX = '[p_{T}(top_{had}^{reco})-p_{T}(top_{had}^{gen})]/p_{T}(top_{had}^{gen})', texY = 'Number of Events',
        attribute = lambda event, sample: event.reso_pt_tophad,
        binning=[40, -2.0, 2.0],
    ))


    plotting.fill(plots, read_variables = read_variables, sequence = sequence)

    ################################################################################
    drawPlots(plots, mode)
    allPlots[mode] = plots


################################################################################
# Add the different channels into SF and all
for mode in ["comb1","comb2","all"]:
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

    if mode == "all":
        drawPlots(allPlots['mumumu'], mode)
        # Write MVA score in root file
        outfile = ROOT.TFile('MVA_score_'+args.era+'.root', 'recreate')
        outfile.cd()
        for plot in plots:
            if "MVA_tWZ" in plot.name:
                model = "MVA_tWZ_3l"
                if "topReco" in plot.name: model = "MVA_tWZ_3l_topReco"
                elif "lstm" in plot.name: model = "MVA_tWZ_3l_lstm"

                for i, l in enumerate(plot.histos):
                    for j, h in enumerate(l):
                        histname = h.GetName()
                        cropped_name = ""
                        if model+"_TWZ_NLO_DR" in histname:
                            node_name = "TWZnode"
                            cropped_name = histname.replace(model+"_TWZ_NLO_DR", "")
                        elif model+"_TTZ" in histname:
                            cropped_name = histname.replace(model+"_TTZ", "")
                            node_name = "TTZnode"
                        else:
                            print "[ERROR] cannot identify node of MVA plot"
                            continue
                        if "TWZ_NLO_DR" in cropped_name: process = "tWZ"
                        elif "TTZ" in cropped_name: process = "ttZ"
                        elif "TTX_rare" in cropped_name: process = "ttX"
                        elif "TZQ" in cropped_name: process = "tZq"
                        elif "WZ" in cropped_name: process = "WZ"
                        elif "ZZ" in cropped_name: process = "ZZ"
                        elif "triBoson" in cropped_name: process = "triBoson"
                        elif "nonprompt" in cropped_name: process = "nonprompt"
                        elif "data" in cropped_name: process = "data"
                        h.Write(model+"__"+process+"__"+node_name)
        outfile.Close()


logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
