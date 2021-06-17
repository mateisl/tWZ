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
argParser.add_argument('--small',          action='store_true', help='Run only on a small subset of the data?', )
argParser.add_argument('--shape',          action='store_true', help='Compare MC shapes?', )
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

if args.shape:                        args.plot_directory += "_shape"
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

def drawPlots(plots, mode, compareShape):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, mode + ("_log" if log else ""), args.selection)
    for plot in plots:
      if not max(l.GetMaximum() for l in sum(plot.histos,[])): continue # Empty plot

      _drawObjects = []

      if isinstance( plot, Plot):
          plotting.draw(plot,
            plot_directory = plot_directory_,
            ratio = {'yRange':(0.1,1.9)} if compareShape else None,
            logX = False, logY = log, sorting = True,
            yRange = (0.03, "auto") if log else (0.001, "auto"),
            scaling = {0:1} if compareShape else {},
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

def matchBjets(event):
    all_idxs = getTopDecays(event)
    bjet_idxs = []
    for idxs in all_idxs:
        for idx in idxs:
            if abs(event.GenPart_pdgId[idx])==5:
                dRmin = 0.4
                bjet_idx = -1
                parton = ROOT.TLorentzVector()
                parton.SetPtEtaPhiM(event.GenPart_pt[idx], event.GenPart_eta[idx], event.GenPart_phi[idx], event.GenPart_mass[idx])
                for i in range(event.nJetGood):
                    jetindex = event.JetGood_index[i]
                    jet  = ROOT.TLorentzVector()
                    jet.SetPtEtaPhiM(event.Jet_pt[jetindex], event.Jet_eta[jetindex], event.Jet_phi[jetindex], event.Jet_mass[jetindex])
                    if jet.DeltaR(parton) < dRmin:
                        dRmin = jet.DeltaR(parton)
                        bjet_idx = i
                if bjet_idx >= 0:
                    bjet_idxs.append(bjet_idx)
    return bjet_idxs

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
                if abs(event.GenPart_pdgId[i]) == 5 or abs(event.GenPart_pdgId[i]) == 24: # particle itself is a W or b
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
        if topidx > 0 and abs(event.GenPart_pdgId[i]) != 22:
            alldecays.append([topidx, i])

    decays_sorted = []
    for topindex in alltops:
        decaysonetop = []
        for topidx, decayidx in alldecays:
            if topindex == topidx:
                decaysonetop.append(decayidx)
        if len(decaysonetop)==3:
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


def calculate_chi2(toplep, tophad, Whad, blep_disc, bhad_disc, mode):
    Mtlep_mean   = 174.
    Mtlep_sigma  =  18.
    Mthad_mean   = 181.
    Mthad_sigma  =  15.
    MWhad_mean   =  80.399
    MWhad_sigma  =  12.
    Mtdiff_sigma = sqrt(pow(Mtlep_sigma,2)+pow(Mthad_sigma,2))
    b_disc_mean  = 1.0
    b_disc_sigma = 0.4
    toplepterm  = pow((toplep.M()-Mtlep_mean)/Mtlep_sigma,2)
    tophadterm  = pow((tophad.M()-Mthad_mean)/Mthad_sigma,2)
    Whadterm    = pow((  Whad.M()-MWhad_mean)/MWhad_sigma,2)
    topdiffterm = pow((toplep.M()-tophad.M())/Mtdiff_sigma,2)
    blepterm    = pow((blep_disc-b_disc_mean)/b_disc_sigma,2)
    bhadterm    = pow((bhad_disc-b_disc_mean)/b_disc_sigma,2)
    if mode=="topdiff":
        chi2 = topdiffterm+Whadterm
    elif mode=="btag":
        chi2 = topdiffterm+Whadterm+blepterm+bhadterm
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
                    jet_permutations.append([[jet_i, jet_j, jet_k, jet_l], [event.JetGood_btagDeepB[i],event.JetGood_btagDeepB[j]]])

    # Get Wlep from lepton + MET
    Wleps   = getWlep(event)
    # build hypotheses
    hypotheses = []
    for Wlep, lepton, neutrino in Wleps:
        for permutation, btag_disc in jet_permutations:
            hypo = {}
            hypo['toplep']       = Wlep + permutation[0]
            hypo['Wlep']         = Wlep
            hypo['lepton']       = lepton
            hypo['neutrino']     = neutrino
            hypo['blep']         = permutation[0]
            hypo['blep_disc']    = btag_disc[0]
            hypo['tophad']       = permutation[1] + permutation[2] + permutation[3]
            hypo['Whad']         = permutation[2] + permutation[3]
            hypo['bhad']         = permutation[1]
            hypo['bhad_disc']    = btag_disc[1]
            hypo['WhadDecay1']   = permutation[2]
            hypo['WhadDecay2']   = permutation[3]
            hypo['chi2']         = calculate_chi2(hypo['toplep'],hypo['tophad'],hypo['Whad'],hypo['bhad_disc'],hypo['blep_disc'],"normal")
            hypo['chi2_topdiff'] = calculate_chi2(hypo['toplep'],hypo['tophad'],hypo['Whad'],hypo['bhad_disc'],hypo['blep_disc'],"topdiff")
            hypo['chi2_btag']    = calculate_chi2(hypo['toplep'],hypo['tophad'],hypo['Whad'],hypo['bhad_disc'],hypo['blep_disc'],"btag")
            hypotheses.append(hypo)
    return hypotheses

def getMatchCategory(event, hypo):
    all_idxs = getTopDecays(event)
    Whadrecos = [ hypo['WhadDecay1'], hypo['WhadDecay2'] ]
    NWhaddecay_match = 0
    MatchedLepton = False
    MatchedBlep = False
    MatchedBhad = False
    Ntophad = 0
    Ntoplep = 0
    for topdecay_idxs in all_idxs:
        if isHadronicDecay(event, topdecay_idxs):
            Ntophad += 1
            for idx in topdecay_idxs:
                if abs(event.GenPart_pdgId[idx])==5:
                    bhad = ROOT.TLorentzVector()
                    bhad.SetPtEtaPhiM(event.GenPart_pt[idx], event.GenPart_eta[idx], event.GenPart_phi[idx], event.GenPart_mass[idx])
                    if bhad.DeltaR(hypo['bhad']) < 0.4:
                        MatchedBhad = True
                else:
                    parton = ROOT.TLorentzVector()
                    parton.SetPtEtaPhiM(event.GenPart_pt[idx], event.GenPart_eta[idx], event.GenPart_phi[idx], event.GenPart_mass[idx])
                    for hadreco in Whadrecos:
                        if parton.DeltaR(hadreco) < 0.4:
                            NWhaddecay_match+=1
        else:
            Ntoplep += 1
            for idx in topdecay_idxs:
                if abs(event.GenPart_pdgId[idx])==11 or abs(event.GenPart_pdgId[idx])==13:
                    lepton = ROOT.TLorentzVector()
                    lepton.SetPtEtaPhiM(event.GenPart_pt[idx], event.GenPart_eta[idx], event.GenPart_phi[idx], event.GenPart_mass[idx])
                    if lepton.DeltaR(hypo['lepton']) < 0.4:
                        MatchedLepton = True
                elif abs(event.GenPart_pdgId[idx])==5:
                    blep = ROOT.TLorentzVector()
                    blep.SetPtEtaPhiM(event.GenPart_pt[idx], event.GenPart_eta[idx], event.GenPart_phi[idx], event.GenPart_mass[idx])
                    if blep.DeltaR(hypo['blep']) < 0.4:
                        MatchedBlep = True

    category = -1
    # First, all cases where semilep ttbar is matchable
    if Ntoplep==1 and Ntophad==1:
        if MatchedBlep and MatchedLepton and MatchedBhad and NWhaddecay_match==2:
            category = 1
        else:
            category = 0
    elif Ntoplep==0 and Ntophad==1:
        if MatchedBhad and NWhaddecay_match==2:
            category = 3
        else:
            category = 2
    elif Ntoplep==1 and Ntophad==0:
        if MatchedBlep and MatchedLepton:
            category = 5
        else:
            category = 4
    else:
        category = 6
    return category

################################################################################
# Define sequences
sequence       = []

# def MatchWs(event, sample):
#     Whads_gen = []
#     Wleps_gen = []
#     Wleps_gen, Whads_gen = getGenWs(event)
#     event.NWlep = len(Wleps_gen)
#     event.NWhad = len(Whads_gen)
#
#     Wlep_reco = getWlep(event)[0][0]
#     Whad_reco = getWhad(event)
#
#     if len(Wleps_gen)!=1 or len(Whads_gen)!=1:
#         return
#
#     event.dR_Wlep_gen_reco = Wlep_reco.DeltaR(Wleps_gen[0])
#     event.dR_Whad_gen_reco = Whad_reco.DeltaR(Whads_gen[0])
#
#     event.mWlep = Wlep_reco.M()
#     event.mWhad = Whad_reco.M()
# sequence.append(MatchWs)

# def MatchWlepDecay(event, sample):
#     dR_lepton = -1
#     dR_neutrino = -1
#     lepton_reco = getWlep(event)[0][1]
#     neutrino_reco = getWlep(event)[0][2]
#     W_idxs = getGenWDecay_indices(event)
#     Wlep_idxs = []
#     for Wdecay1,Wdecay2 in W_idxs:
#         if not isHadronicDecay(event, [Wdecay1,Wdecay2]):
#             Wlep_idxs.append([Wdecay1,Wdecay2])
#     event.NWlep_matching = len(Wlep_idxs)
#     if len(Wlep_idxs) != 1:
#         return
#     else:
#         lepton_gen = ROOT.TLorentzVector()
#         neutrino_gen = ROOT.TLorentzVector()
#         if abs(event.GenPart_pdgId[Wdecay1]) in [11,13,15]:
#             lepton_gen.SetPtEtaPhiM(event.GenPart_pt[Wdecay1], event.GenPart_eta[Wdecay1], event.GenPart_phi[Wdecay1], event.GenPart_mass[Wdecay1])
#             neutrino_gen.SetPtEtaPhiM(event.GenPart_pt[Wdecay2], event.GenPart_eta[Wdecay2], event.GenPart_phi[Wdecay2], event.GenPart_mass[Wdecay2])
#         else:
#             lepton_gen.SetPtEtaPhiM(event.GenPart_pt[Wdecay2], event.GenPart_eta[Wdecay2], event.GenPart_phi[Wdecay2], event.GenPart_mass[Wdecay2])
#             neutrino_gen.SetPtEtaPhiM(event.GenPart_pt[Wdecay1], event.GenPart_eta[Wdecay1], event.GenPart_phi[Wdecay1], event.GenPart_mass[Wdecay1])
#         dR_lepton = lepton_reco.DeltaR(lepton_gen)
#         dR_neutrino = neutrino_reco.DeltaR(neutrino_gen)
#     event.dR_lepton_gen_reco = dR_lepton
#     event.dR_neutrino_gen_reco = dR_neutrino
# sequence.append(MatchWlepDecay)

# def BestWhadReco(event, sample):
#     decays = getGenWhadDecays(event)
#     if len(decays) != 1:
#         event.mWhad_matching = float('nan')
#         event.maxDR_Whaddecays = float('nan')
#         return
#
#     WDecay1_had = decays[0][0]
#     WDecay2_had = decays[0][1]
#
#     min_dR = 1000
#     jet1 = ROOT.TLorentzVector()
#     i_jet1 = -1
#     for i in range(event.nJetGood):
#         jetindex = event.JetGood_index[i]
#         jet_tmp  = ROOT.TLorentzVector()
#         jet_tmp.SetPtEtaPhiM(event.Jet_pt[jetindex], event.Jet_eta[jetindex], event.Jet_phi[jetindex], event.Jet_mass[jetindex])
#         if WDecay1_had.DeltaR(jet_tmp) < min_dR:
#             min_dR = WDecay1_had.DeltaR(jet_tmp)
#             jet1 = jet_tmp
#             i_jet1 = i
#
#     min_dR = 1000 # reset min_dR
#     jet2 = ROOT.TLorentzVector()
#     for i in range(event.nJetGood):
#         if i == i_jet1:
#             continue
#         jetindex = event.JetGood_index[i]
#         jet_tmp  = ROOT.TLorentzVector()
#         jet_tmp.SetPtEtaPhiM(event.Jet_pt[jetindex], event.Jet_eta[jetindex], event.Jet_phi[jetindex], event.Jet_mass[jetindex])
#         if WDecay2_had.DeltaR(jet_tmp) < min_dR:
#             min_dR = WDecay2_had.DeltaR(jet_tmp)
#             jet2 = jet_tmp
#
#     Wreco = jet1+jet2
#     event.mWhad_matching = Wreco.M()
#     event.maxDR_Whaddecays = WDecay1_had.DeltaR(jet1) if WDecay1_had.DeltaR(jet1) > WDecay2_had.DeltaR(jet2) else WDecay2_had.DeltaR(jet2)
# sequence.append(BestWhadReco)

def Matching(event, sample):
    top_idxs = getGenTop_indices(event, "production")
    bot_idxs = getGenBottomFromTop_indices(event)
    event.Ntops = len(top_idxs)
    event.Nbots = len(bot_idxs)
    event.mtophad_matched = -1
    event.mtophadreduced_matched = -1
    event.mWhad_matched = -1

    topdecay_idxs = getTopDecays(event)
    Nt=0
    Nthad=0
    Ntlep=0
    numberoftopdecays = []
    for topdecays in topdecay_idxs:
        numberoftopdecays.append(len(topdecays))
        Wdecays = []
        for decay in topdecays:
            if abs(event.GenPart_pdgId[decay]) != 5:
                Wdecays.append(decay)
        if len(topdecays) == 3:
            Nt+=1
            if isHadronicDecay(event, topdecays):
                Nthad+=1
                matched_jets = getMatchingJets(event, topdecays)
                event.Nmatched_tophad = len(matched_jets)
                if len(matched_jets) == 3:
                    top = matched_jets[0]+matched_jets[1]+matched_jets[2]
                    event.mtophad_matched = top.M()
                elif len(matched_jets) == 2:
                    top = matched_jets[0]+matched_jets[1]
                    event.mtophadreduced_matched = top.M()
            else:
                Ntlep+=1
        if len(Wdecays) == 2:
            if isHadronicDecay(event, Wdecays):
                matched_jets_W = getMatchingJets(event, Wdecays)
                if len(matched_jets_W) == 2:
                    W = matched_jets[0]+matched_jets[1]
                    event.mWhad_matched = W.M()

    event.Ntops_matching = Nt
    event.Ntopshad_matching = Nthad
    event.Ntopslep_matching = Ntlep
    event.Ntopdecays = numberoftopdecays
sequence.append(Matching)

def TopReco(event, sample):
    hypotheses=getTopHypos(event, 6)
    chi2min = 10000
    chi2min_topdiff = 10000
    chi2min_btag = 10000
    hypo_selected = hypotheses[0]
    hypo_selected_topdiff = hypotheses[0]
    hypo_selected_btag = hypotheses[0]
    foundhypo = False
    foundhypo_topdiff = False
    foundhypo_btag = False
    for hypo in hypotheses:
        if hypo['chi2']<chi2min:
            chi2min = hypo['chi2']
            hypo_selected = hypo
            foundhypo = True
        if hypo['chi2_topdiff']<chi2min_topdiff:
            chi2min_topdiff = hypo['chi2_topdiff']
            hypo_selected_topdiff = hypo
            foundhypo_topdiff = True
        if hypo['chi2_btag']<chi2min_btag:
            chi2min_btag = hypo['chi2_btag']
            hypo_selected_btag = hypo
            foundhypo_btag = True

    if foundhypo:
        if hypo_selected['toplep'].M() > hypo_selected['tophad'].M():
            mtop_hi = hypo_selected['toplep'].M()
            mtop_lo = hypo_selected['tophad'].M()
        else:
            mtop_hi = hypo_selected['tophad'].M()
            mtop_lo = hypo_selected['toplep'].M()
    if foundhypo_topdiff:
        if hypo_selected_topdiff['toplep'].M() > hypo_selected_topdiff['tophad'].M():
            mtop_hi_topdiff = hypo_selected_topdiff['toplep'].M()
            mtop_lo_topdiff = hypo_selected_topdiff['tophad'].M()
        else:
            mtop_hi_topdiff = hypo_selected_topdiff['tophad'].M()
            mtop_lo_topdiff = hypo_selected_topdiff['toplep'].M()
    if foundhypo_btag:
        if hypo_selected_btag['toplep'].M() > hypo_selected_btag['tophad'].M():
            mtop_hi_btag = hypo_selected_btag['toplep'].M()
            mtop_lo_btag = hypo_selected_btag['tophad'].M()
        else:
            mtop_hi_btag = hypo_selected_btag['tophad'].M()
            mtop_lo_btag = hypo_selected_btag['toplep'].M()

    # Also get Z
    Z = ROOT.TLorentzVector()
    Z.SetPtEtaPhiM(event.Z1_pt, event.Z1_eta, event.Z1_phi, event.Z1_mass)

    # observables for mtop1 and mtop2 closest to mtop
    event.mtop_average = (hypo_selected['toplep'].M()+hypo_selected['tophad'].M())/2 if foundhypo else -1
    event.mtoplep = hypo_selected['toplep'].M() if foundhypo else -1
    event.mtophad = hypo_selected['tophad'].M() if foundhypo else -1
    event.mtop_hi = mtop_hi if foundhypo else -1
    event.mtop_lo = mtop_lo if foundhypo else -1
    event.mtop_diff = (mtop_hi-mtop_lo)/(mtop_hi+mtop_lo) if foundhypo else -1
    event.pt_diff = abs(hypo_selected['toplep'].Pt()-hypo_selected['tophad'].Pt()) if foundhypo else -1
    event.dR_tops = hypo_selected['toplep'].DeltaR(hypo_selected['tophad']) if foundhypo else -1
    event.mW_lep = hypo_selected['Wlep'].M() if foundhypo else -1
    event.mW_had = hypo_selected['Whad'].M() if foundhypo else -1
    event.ptW_lep = hypo_selected['Wlep'].Pt() if foundhypo else -1
    event.ptW_had = hypo_selected['Whad'].Pt() if foundhypo else -1
    event.chi2 = hypo_selected['chi2'] if foundhypo else -1
    event.pgof = exp(-0.5*hypo_selected['chi2']) if foundhypo else -1
    event.hypofound = 1 if foundhypo else 0
    event.dR_toplep_Z = hypo_selected['toplep'].DeltaR(Z) if foundhypo else 0
    event.dR_tophad_Z = hypo_selected['tophad'].DeltaR(Z) if foundhypo else 0

    # observables for mtop1 closest to mtop2
    event.mtop_average_topdiff = (hypo_selected_topdiff['toplep'].M()+hypo_selected_topdiff['tophad'].M())/2 if foundhypo_topdiff else -1
    event.mtoplep_topdiff = hypo_selected_topdiff['toplep'].M() if foundhypo_topdiff else -1
    event.mtophad_topdiff = hypo_selected_topdiff['tophad'].M() if foundhypo_topdiff else -1
    event.mtop_hi_topdiff = mtop_hi_topdiff if foundhypo_topdiff else -1
    event.mtop_lo_topdiff = mtop_lo_topdiff if foundhypo_topdiff else -1
    event.mtop_diff_topdiff = (mtop_hi_topdiff-mtop_lo_topdiff)/(mtop_hi_topdiff+mtop_lo_topdiff) if foundhypo_topdiff else -1
    event.pt_diff_topdiff = abs(hypo_selected_topdiff['toplep'].Pt()-hypo_selected_topdiff['tophad'].Pt()) if foundhypo_topdiff else -1
    event.dR_tops_topdiff = hypo_selected_topdiff['toplep'].DeltaR(hypo_selected_topdiff['tophad']) if foundhypo_topdiff else -1
    event.mW_lep_topdiff = hypo_selected_topdiff['Wlep'].M() if foundhypo_topdiff else -1
    event.mW_had_topdiff = hypo_selected_topdiff['Whad'].M() if foundhypo_topdiff else -1
    event.ptW_lep_topdiff = hypo_selected_topdiff['Wlep'].Pt() if foundhypo_topdiff else -1
    event.ptW_had_topdiff = hypo_selected_topdiff['Whad'].Pt() if foundhypo_topdiff else -1
    event.chi2_topdiff = hypo_selected_topdiff['chi2_topdiff'] if foundhypo_topdiff else -1
    event.pgof_topdiff = exp(-0.5*hypo_selected_topdiff['chi2']) if foundhypo_topdiff else -1
    event.hypofound_topdiff = 1 if foundhypo_topdiff else 0
    event.dR_toplep_Z_topdiff = hypo_selected_topdiff['toplep'].DeltaR(Z) if foundhypo_topdiff else 0
    event.dR_tophad_Z_topdiff = hypo_selected_topdiff['tophad'].DeltaR(Z) if foundhypo_topdiff else 0

    # observables for mtop1 closest to mtop2 + btag
    event.mtop_average_btag = (hypo_selected_btag['toplep'].M()+hypo_selected_btag['tophad'].M())/2 if foundhypo_btag else -1
    event.mtoplep_btag = hypo_selected_btag['toplep'].M() if foundhypo_btag else -1
    event.mtophad_btag = hypo_selected_btag['tophad'].M() if foundhypo_btag else -1
    event.mtop_hi_btag = mtop_hi_btag if foundhypo_btag else -1
    event.mtop_lo_btag = mtop_lo_btag if foundhypo_btag else -1
    event.mtop_diff_btag = (mtop_hi_btag-mtop_lo_btag)/(mtop_hi_btag+mtop_lo_btag) if foundhypo_btag else -1
    event.pt_diff_btag = abs(hypo_selected_btag['toplep'].Pt()-hypo_selected_btag['tophad'].Pt()) if foundhypo_btag else -1
    event.dR_tops_btag = hypo_selected_btag['toplep'].DeltaR(hypo_selected_btag['tophad']) if foundhypo_btag else -1
    event.mW_lep_btag = hypo_selected_btag['Wlep'].M() if foundhypo_btag else -1
    event.mW_had_btag = hypo_selected_btag['Whad'].M() if foundhypo_btag else -1
    event.ptW_lep_btag = hypo_selected_btag['Wlep'].Pt() if foundhypo_btag else -1
    event.ptW_had_btag = hypo_selected_btag['Whad'].Pt() if foundhypo_btag else -1
    event.chi2_btag = hypo_selected_btag['chi2_btag'] if foundhypo_btag else -1
    event.pgof_btag = exp(-0.5*hypo_selected_btag['chi2']) if foundhypo_btag else -1
    event.hypofound_btag = 1 if foundhypo_btag else 0
    event.dR_toplep_Z_btag = hypo_selected_btag['toplep'].DeltaR(Z) if foundhypo_btag else 0
    event.dR_tophad_Z_btag = hypo_selected_btag['tophad'].DeltaR(Z) if foundhypo_btag else 0

    # Match hypo with b tag in event
    bjet_goodjet_idx = getBJetindex(event)
    bjet_idx = event.JetGood_index[bjet_goodjet_idx]
    bjet = ROOT.TLorentzVector()
    bjet.SetPtEtaPhiM(event.Jet_pt[bjet_idx], event.Jet_eta[bjet_idx], event.Jet_phi[bjet_idx], event.Jet_mass[bjet_idx])
    bmatch = 0
    if bjet.DeltaR(hypo_selected['blep']) < 0.01:
        bmatch = 1
    elif bjet.DeltaR(hypo_selected['bhad']) < 0.01:
        bmatch = 2
    bmatch_topdiff = 0
    if bjet.DeltaR(hypo_selected_topdiff['blep']) < 0.01:
        bmatch_topdiff = 1
    elif bjet.DeltaR(hypo_selected_topdiff['bhad']) < 0.01:
        bmatch_topdiff = 2
    bmatch_btag = 0
    if bjet.DeltaR(hypo_selected_btag['blep']) < 0.01:
        bmatch_btag = 1
    elif bjet.DeltaR(hypo_selected_btag['bhad']) < 0.01:
        bmatch_btag = 2
    event.bmatch = bmatch
    event.bmatch_topdiff = bmatch_topdiff
    event.bmatch_btag = bmatch_btag

    # Match with Gen Level
    event.category = getMatchCategory(event, hypo_selected) if foundhypo else -1
    event.category_topdiff = getMatchCategory(event, hypo_selected_topdiff) if foundhypo_topdiff else -1
    event.category_btag = getMatchCategory(event, hypo_selected_btag) if foundhypo_btag else -1
sequence.append(TopReco)

def getMatchedBjetProperties(event, sample):
    btag_values = []
    jet_idxs = matchBjets(event)
    for idx in jet_idxs:
        btag_values.append(event.JetGood_btagDeepB[idx])
    event.btag_disc_matched = btag_values
sequence.append(getMatchedBjetProperties)

# def printGenParticles(event,sample):
#     print '-------'
#     for i in range(event.nGenPart):
#         if abs(event.GenPart_pdgId[i]) == 5:
#             print event.GenPart_status[i]
# sequence.append(printGenParticles)

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

    if args.shape:
        mc[0].style = styles.lineStyle(mc[0].color, width=2)
        mc[1].style = styles.fillStyle(mc[1].color)
    else:
        for sample in mc: sample.style = styles.fillStyle(sample.color)

    for sample in mc:
      sample.read_variables = read_variables_MC
      sample.setSelectionString([getLeptonSelection(mode)])
      sample.weight = lambda event, sample: event.reweightBTag_SF*event.reweightPU*event.reweightL1Prefire*event.reweightTrigger#*event.reweightLeptonSF

    if args.shape:
        stack = Stack([mc[1]], [mc[0]])
    else:
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
        name = 'm_W_matched',
        texX = 'm(W_{had}) after matching', texY = 'Number of Events',
        attribute = lambda event, sample: event.mWhad_matched,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_reduced_matched',
        texX = 'm(reduced top_{had}) after matching', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtophadreduced_matched,
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
        name = 'btag_disc_matched',
        texX = 'b tag score of matched b jets', texY = 'Number of Events',
        attribute = lambda event, sample: event.btag_disc_matched,
        binning=[25,0.0,1.0],
    ))

    plots.append(Plot(
        name = "Zmother",
        texX = 'pdg ID of Z mother', texY = 'Number of Events',
        attribute = lambda event, sample: event.Zmother,
        binning=[26,-0.5, 25.5],
    ))


    ############################################################################
    ## Top Reco Plots

    plots.append(Plot(
        name = 'm_top_average',
        texX = '(m_{t1}+m_{t2})/2', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_average,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_lep',
        texX = 'm_{tlep}', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtoplep,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_had',
        texX = 'm_{thad}', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtophad,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_hi',
        texX = 'max(m_{t1},m_{t2})', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_hi,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_lo',
        texX = 'min(m_{t1},m_{t2})', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_lo,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_diff',
        texX = '|m_{t1}-m_{t2}|/|m_{t1}+m_{t2}|', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_diff,
        binning=[50, 0.0, 1.0],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'pt_top_diff',
        texX = '|p_{T}^{t1}-p_{T}^{t2}|', texY = 'Number of Events',
        attribute = lambda event, sample: event.pt_diff,
        binning=[40, 0.0, 400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'dR_tops',
        texX = '#Delta R(top_{lep},top_{had})', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_tops,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'mW_lep',
        texX = 'm_{Wlep}', texY = 'Number of Events',
        attribute = lambda event, sample: event.mW_lep,
        binning=[40, 0.0, 200],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'mW_had',
        texX = 'm_{Whad}', texY = 'Number of Events',
        attribute = lambda event, sample: event.mW_had,
        binning=[40, 0.0, 200],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'ptW_lep',
        texX = 'p_{T}^{Wlep}', texY = 'Number of Events',
        attribute = lambda event, sample: event.ptW_lep,
        binning=[40, 0.0, 400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'ptW_had',
        texX = 'p_{T}^{Whad}', texY = 'Number of Events',
        attribute = lambda event, sample: event.ptW_had,
        binning=[40, 0.0, 400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'bmatch',
        texX = ' ', texY = 'Number of Events',
        attribute = lambda event, sample: event.bmatch,
        binning=[3, -0.5, 2.5],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'chi2',
        texX = '#chi^{2}', texY = 'Number of Events',
        attribute = lambda event, sample: event.chi2,
        binning=[40, 0.0, 80.],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'p_gof',
        texX = 'P_{gof}', texY = 'Number of Events',
        attribute = lambda event, sample: event.pgof,
        binning=[50, 0., 1.],
    ))

    plots.append(Plot(
        name = 'hypofound',
        texX = 'found ttbar hypothesis', texY = 'Number of Events',
        attribute = lambda event, sample: event.hypofound,
        binning=[2, -0.5, 1.5],
    ))

    plots.append(Plot(
        name = 'category',
        texX = 'matching category', texY = 'Number of Events',
        attribute = lambda event, sample: event.category,
        binning=[7, -0.5, 6.5],
    ))

    plots.append(Plot(
        name = 'dR_toplep_Z',
        texX = '#Delta R(top_{lep},Z)', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_toplep_Z,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'dR_tophad_Z',
        texX = '#Delta R(top_{had},Z)', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_tophad_Z,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))
    ############################################################################
    ## Top Reco Plots
    plots.append(Plot(
        name = 'm_top_average_topdiff',
        texX = '(m_{t1}+m_{t2})/2 (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_average_topdiff,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_lep_topdiff',
        texX = 'm_{tlep} (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtoplep_topdiff,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_had_topdiff',
        texX = 'm_{thad} (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtophad_topdiff,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_hi_topdiff',
        texX = 'max(m_{t1},m_{t2}) (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_hi_topdiff,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_lo_topdiff',
        texX = 'min(m_{t1},m_{t2}) (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_lo_topdiff,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_diff_topdiff',
        texX = '|m_{t1}-m_{t2}|/|m_{t1}+m_{t2}| (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_diff_topdiff,
        binning=[50, 0.0, 1.0],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'pt_top_diff_topdiff',
        texX = '|p_{T}^{t1}-p_{T}^{t2}| (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.pt_diff_topdiff,
        binning=[40, 0.0, 400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'dR_tops_topdiff',
        texX = '#Delta R(top_{lep},top_{had}) (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_tops_topdiff,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'mW_lep_topdiff',
        texX = 'm_{Wlep} (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mW_lep_topdiff,
        binning=[40, 0.0, 200],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'mW_had_topdiff',
        texX = 'm_{Whad} (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mW_had_topdiff,
        binning=[40, 0.0, 200],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'ptW_lep_topdiff',
        texX = 'p_{T}^{Wlep} (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.ptW_lep_topdiff,
        binning=[40, 0.0, 400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'ptW_had_topdiff',
        texX = 'p_{T}^{Whad} (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.ptW_had_topdiff,
        binning=[40, 0.0, 400],
        addOverFlowBin='upper',
    ))


    plots.append(Plot(
        name = 'bmatch_topdiff',
        texX = ' ', texY = 'Number of Events',
        attribute = lambda event, sample: event.bmatch_topdiff,
        binning=[3, -0.5, 2.5],
        addOverFlowBin='upper',
    ))


    plots.append(Plot(
        name = 'chi2_topdiff',
        texX = '#chi^{2} (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.chi2_topdiff,
        binning=[40, 0.0, 80.],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'p_gof_topdiff',
        texX = 'P_{gof} (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.pgof_topdiff,
        binning=[50, 0., 1.],
    ))



    plots.append(Plot(
        name = 'hypofound_topdiff',
        texX = 'found ttbar hypothesis (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.hypofound_topdiff,
        binning=[2, -0.5, 1.5],
    ))

    plots.append(Plot(
        name = 'category_topdiff',
        texX = 'matching category (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.category_topdiff,
        binning=[7, -0.5, 6.5],
    ))

    plots.append(Plot(
        name = 'dR_toplep_Z_topdiff',
        texX = '#Delta R(top_{lep},Z) (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_toplep_Z_topdiff,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'dR_tophad_Z_topdiff',
        texX = '#Delta R(top_{had},Z) (top symmetry)', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_tophad_Z_topdiff,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    ############################################################################
    ## Top Reco Plots
    plots.append(Plot(
        name = 'm_top_average_btag',
        texX = '(m_{t1}+m_{t2})/2 (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_average_btag,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_lep_btag',
        texX = 'm_{tlep} (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtoplep_btag,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_had_btag',
        texX = 'm_{thad} (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtophad_btag,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_hi_btag',
        texX = 'max(m_{t1},m_{t2}) (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_hi_btag,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_lo_btag',
        texX = 'min(m_{t1},m_{t2}) (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_lo_btag,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_top_diff_btag',
        texX = '|m_{t1}-m_{t2}|/|m_{t1}+m_{t2}| (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtop_diff_btag,
        binning=[50, 0.0, 1.0],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'pt_top_diff_btag',
        texX = '|p_{T}^{t1}-p_{T}^{t2}| (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.pt_diff_btag,
        binning=[40, 0.0, 400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'dR_tops_btag',
        texX = '#Delta R(top_{lep},top_{had}) (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_tops_btag,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'mW_lep_btag',
        texX = 'm_{Wlep} (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mW_lep_btag,
        binning=[40, 0.0, 200],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'mW_had_btag',
        texX = 'm_{Whad} (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.mW_had_btag,
        binning=[40, 0.0, 200],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'ptW_lep_btag',
        texX = 'p_{T}^{Wlep} (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.ptW_lep_btag,
        binning=[40, 0.0, 400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'ptW_had_btag',
        texX = 'p_{T}^{Whad} (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.ptW_had_btag,
        binning=[40, 0.0, 400],
        addOverFlowBin='upper',
    ))


    plots.append(Plot(
        name = 'bmatch_btag',
        texX = ' ', texY = 'Number of Events',
        attribute = lambda event, sample: event.bmatch_btag,
        binning=[3, -0.5, 2.5],
        addOverFlowBin='upper',
    ))


    plots.append(Plot(
        name = 'chi2_btag',
        texX = '#chi^{2} (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.chi2_btag,
        binning=[40, 0.0, 80.],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'p_gof_btag',
        texX = 'P_{gof} (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.pgof_btag,
        binning=[50, 0., 1.],
    ))



    plots.append(Plot(
        name = 'hypofound_btag',
        texX = 'found ttbar hypothesis (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.hypofound_btag,
        binning=[2, -0.5, 1.5],
    ))

    plots.append(Plot(
        name = 'category_btag',
        texX = 'matching category (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.category_btag,
        binning=[7, -0.5, 6.5],
    ))

    plots.append(Plot(
        name = 'dR_toplep_Z_btag',
        texX = '#Delta R(top_{lep},Z) (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_toplep_Z_btag,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plots.append(Plot(
        name = 'dR_tophad_Z_btag',
        texX = '#Delta R(top_{had},Z) (top symmetry + b discriminator)', texY = 'Number of Events',
        attribute = lambda event, sample: event.dR_tophad_Z_btag,
        binning=[35, 0.0, 7],
        addOverFlowBin = 'upper',
    ))

    plotting.fill(plots, read_variables = read_variables, sequence = sequence)

    ################################################################################
    drawPlots(plots, mode, args.shape)
    allPlots[mode] = plots

    for plot in plots:
      if "bmatch" in plot.name:
        for i, l in enumerate(plot.histos):
          for j, h in enumerate(l):
            h.GetXaxis().SetBinLabel(1, "not matched")
            h.GetXaxis().SetBinLabel(2, "b_{lep}")
            h.GetXaxis().SetBinLabel(3, "b_{had}")

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
        drawPlots(allPlots['mumumu'], mode, args.shape)
        # Write matched hists in root file
        plot_dir = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era)
        outfile = ROOT.TFile(plot_dir+'/MatchHists.root', 'recreate')
        outfile.cd()
        for plot in plots:
            if "_matched" in plot.name:
                for i, l in enumerate(plot.histos):
                    for j, h in enumerate(l):
                        histname = h.GetName()
                        if "TWZ_NLO_DR" in histname: process = "tWZ"
                        elif "TTZ" in histname: process = "ttZ"
                        h.Write(plot.name+"__"+process)
        outfile.Close()




logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
