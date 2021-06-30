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
argParser.add_argument('--selection',      action='store', default='trilepT-minDLmass12-onZ1-njet4p-deepjet1')
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


def getBJetindex( event ):
    maxscore = 0.0
    index = -1
    for i in range(event.nJetGood):
        btagscore = event.JetGood_btagDeepB[i]
        if btagscore > maxscore:
            maxscore = btagscore
            index = i
    return index

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
        min_dR = 0.2
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
    isHadronic = True
    for i in idxs:
        if abs(event.GenPart_pdgId[i]) not in quarkIDs:
            isHadronic = False
    return isHadronic


def calculate_chi2(toplep, tophad, Whad, blep_disc, bhad_disc, mode):
    Mtlep_mean   = 171.
    Mtlep_sigma  =  16.
    Mthad_mean   = 171.
    Mthad_sigma  =  17.
    MWhad_mean   =  83.
    MWhad_sigma  =  11.
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
            hypo['chi2']         = calculate_chi2(hypo['toplep'],hypo['tophad'],hypo['Whad'],hypo['blep_disc'],hypo['bhad_disc'],"normal")
            hypotheses.append(hypo)
    return hypotheses

def getMatchCategory(event, hypo):
    leptonIDs = [11,13,15]
    neutrinoIDs = [12,14,16]
    matchdR = 0.4
    all_idxs = getTopDecays(event)
    categories = []
    for topdecay_idxs in all_idxs:
        # If decay is hadronic, check jets, each parton
        if isHadronicDecay(event, topdecay_idxs):
            hadmatch = 1
            if len(topdecay_idxs) != 3: hadmatch = 4
            elif not isMatchable(event,[topdecay_idxs]): hadmatch = 3
            else:
                for idx in topdecay_idxs:
                    parton = ROOT.TLorentzVector()
                    parton.SetPtEtaPhiM(event.GenPart_pt[idx], event.GenPart_eta[idx], event.GenPart_phi[idx], event.GenPart_mass[idx])
                    if parton.DeltaR(hypo['bhad']) > matchdR and parton.DeltaR(hypo['WhadDecay1']) > matchdR and parton.DeltaR(hypo['WhadDecay2']) > matchdR:
                        hadmatch=2
            categories.append(hadmatch)
        else:
            lepmatch = 1
            if len(topdecay_idxs) != 3: lepmatch = 4
            elif not isMatchable(event,[topdecay_idxs]): lepmatch = 3
            else:
                for idx in topdecay_idxs:
                    parton = ROOT.TLorentzVector()
                    parton.SetPtEtaPhiM(event.GenPart_pt[idx], event.GenPart_eta[idx], event.GenPart_phi[idx], event.GenPart_mass[idx])
                    if abs(event.GenPart_pdgId[idx]) == 5:
                        if parton.DeltaR(hypo['blep']) > matchdR:
                            lepmatch=2
                    elif abs(event.GenPart_pdgId[idx]) in leptonIDs:
                        if parton.DeltaR(hypo['lepton']) > matchdR:
                            lepmatch=2
                    # elif abs(event.GenPart_pdgId[idx]) in neutrinoIDs:
                    #     if parton.DeltaR(hypo['neutrino']) > matchdR:
                    #         lepmatch=2
            categories.append(lepmatch+4)

    return categories

def getBGenJets(event):
    genjets = []
    for i in range(event.nGenJet):
        if abs(event.GenJet_partonFlavour[i]) == 5:
            genjet = ROOT.TLorentzVector()
            genjet.SetPtEtaPhiM(event.GenJet_pt[i],event.GenJet_eta[i],event.GenJet_phi[i],event.GenJet_mass[i])
            genjets.append(genjet)
    return genjets

def isMatchable(event,topdecay_idxs):
    leptonIDs = [11,13,15]
    neutrinoIDs = [12,14,16]
    quarkIDs = [1,2,3,4,5]
    matchdR = 0.4
    foundlepton = True
    foundjets = True
    for topdecays in topdecay_idxs:
        for idx in topdecays:
            parton = ROOT.TLorentzVector()
            parton.SetPtEtaPhiM(event.GenPart_pt[idx], event.GenPart_eta[idx], event.GenPart_phi[idx], event.GenPart_mass[idx])
            if abs(event.GenPart_pdgId[idx]) in leptonIDs:
                recolepton = ROOT.TLorentzVector()
                recolepton.SetPtEtaPhiM(event.lep_pt[event.nonZ1_l1_index], event.lep_eta[event.nonZ1_l1_index], event.lep_phi[event.nonZ1_l1_index], 0)
                if parton.DeltaR(recolepton) > matchdR: foundlepton = False
            elif abs(event.GenPart_pdgId[idx]) in quarkIDs:
                foundjet = False
                for i in range(event.nJetGood):
                    jetidx = event.JetGood_index[i]
                    recojet = ROOT.TLorentzVector()
                    recojet.SetPtEtaPhiM(event.Jet_pt[jetidx], event.Jet_eta[jetidx], event.Jet_phi[jetidx], event.Jet_mass[jetidx])
                    if parton.DeltaR(recojet) < 0.4: foundjet = True
                if not foundjet: foundjets = False
    ismatchable = (foundlepton and foundjets)
    return ismatchable





################################################################################
# Define sequences
sequence       = []

def Matching(event, sample):
    top_idxs = getGenTop_indices(event, "production")
    bot_idxs = getGenBottomFromTop_indices(event)
    event.Ntops = len(top_idxs)
    event.Nbots = len(bot_idxs)
    event.mtophad_matched = -1
    event.mtoplep_matched = -1
    event.mtophadreduced_matched = -1
    event.mWhad_matched = -1

    topdecay_idxs = getTopDecays(event)
    event.ismatchable = isMatchable(event,topdecay_idxs)
    Nt=0
    Nthad=0
    Ntlep=0
    numberoftopdecays = []
    for topdecays in topdecay_idxs:
        numberoftopdecays.append(len(topdecays))
        # save W decay partons in additional list
        Wdecays = []
        for decay in topdecays:
            if abs(event.GenPart_pdgId[decay]) != 5:
                Wdecays.append(decay)
        # Firts look at top decay
        if len(topdecays) == 3:
            Nt+=1
            # match hadronic top and reconstruct mass
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
            # match leptonic top and reconstruct mass
            else:
                Ntlep+=1
                bparton = []
                lepton_gen = ROOT.TLorentzVector()
                neutrino_gen = ROOT.TLorentzVector()
                leptonIDs = [11,13,15]
                neutrinoIDs = [12,14,16]
                foundb = False
                foundlep = False
                foundneu = False
                for idx in topdecays:
                    if abs(event.GenPart_pdgId[idx]) == 5:
                        bparton.append(idx)
                        foundb = True
                    elif abs(event.GenPart_pdgId[idx]) in leptonIDs:
                        lepton_gen.SetPtEtaPhiM(event.GenPart_pt[idx],event.GenPart_eta[idx],event.GenPart_phi[idx],event.GenPart_mass[idx])
                        foundlep = True
                    elif abs(event.GenPart_pdgId[idx]) in neutrinoIDs:
                        neutrino_gen.SetPtEtaPhiM(event.GenPart_pt[idx],event.GenPart_eta[idx],event.GenPart_phi[idx],event.GenPart_mass[idx])
                        foundneu = True
                if len(bparton) == 1:
                    matched_jets_blep = getMatchingJets(event, bparton)
                    if len(matched_jets_blep) == 1 and foundb and foundlep and foundneu:
                        Wlephypos = getWlep(event)
                        for Wlep, lepton, neutrino in Wlephypos:
                            if lepton.DeltaR(lepton_gen) < 0.3 and neutrino.DeltaR(neutrino_gen) < 0.3:
                                toplep = lepton+neutrino+matched_jets_blep[0]
                                event.mtoplep_matched = toplep.M()
        # match hadronic W and reconstruct mass
        if len(Wdecays) == 2:
            if isHadronicDecay(event, Wdecays):
                matched_jets_W = getMatchingJets(event, Wdecays)
                if len(matched_jets_W) == 2:
                    W = matched_jets_W[0]+matched_jets_W[1]
                    event.mWhad_matched = W.M()

    event.Ntops_matching = Nt
    event.Ntopshad_matching = Nthad
    event.Ntopslep_matching = Ntlep
    event.Ntopdecays = numberoftopdecays
sequence.append(Matching)

def TopReco(event, sample):
    hypotheses=getTopHypos(event, 6)
    chi2min = 10000
    foundhypo = False
    for hypo in hypotheses:
        if hypo['chi2']<chi2min:
            chi2min = hypo['chi2']
            hypo_selected = hypo
            foundhypo = True
    if foundhypo:
        if hypo_selected['toplep'].M() > hypo_selected['tophad'].M():
            mtop_hi = hypo_selected['toplep'].M()
            mtop_lo = hypo_selected['tophad'].M()
        else:
            mtop_hi = hypo_selected['tophad'].M()
            mtop_lo = hypo_selected['toplep'].M()

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
    event.dR_toplep_Z = hypo_selected['toplep'].DeltaR(Z) if foundhypo else -1
    event.dR_tophad_Z = hypo_selected['tophad'].DeltaR(Z) if foundhypo else -1
    event.blep_disc = hypo_selected['blep_disc'] if foundhypo else -1
    event.bhad_disc = hypo_selected['bhad_disc'] if foundhypo else -1

    # Match hypo with b genjet in event
    genBjets = getBGenJets(event)
    bmatch         = 0
    for genjet in genBjets:
        if   genjet.DeltaR(hypo_selected['blep']) < 0.4: bmatch += 1
        elif genjet.DeltaR(hypo_selected['bhad']) < 0.4: bmatch += 2

    event.bmatch         = bmatch         if bmatch < 4 else 4

    # Match with Gen Level
    event.category = getMatchCategory(event, hypo_selected) if foundhypo else -1

sequence.append(TopReco)

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

def bjetMatch(event, sample):
    btagscoreList = []
    genBjets = getBGenJets(event)
    for genjet in genBjets:
        dRmin = 0.4
        i_match = -1
        for i in range(event.nJetGood):
            jetidx = event.JetGood_index[i]
            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(event.Jet_pt[jetidx], event.Jet_eta[jetidx], event.Jet_phi[jetidx], event.Jet_mass[jetidx])
            if jet.DeltaR(genjet) < dRmin:
                dRmin = jet.DeltaR(genjet)
                i_match = i
        if i_match != -1:
            btagscoreList.append(event.JetGood_btagDeepB[i_match])
    event.btagscore_matched = btagscoreList
sequence.append(bjetMatch)

def genJetFlavor(event, sample):
    nbjet = 0
    nlightjet = 0
    for i in range(event.nGenJet):
        if abs(event.GenJet_partonFlavour[i]) == 5:
            nbjet += 1
        elif abs(event.GenJet_partonFlavour[i]) in [1,2,3,4]:
            nlightjet += 1
    event.nbjet = nbjet
    event.nlightjet = nlightjet

sequence.append(genJetFlavor)

################################################################################
# Read variables

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",  "nJet/I",
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
]

read_variables_MC = [
    'reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F',
    "genZ1_pt/F", "genZ1_eta/F", "genZ1_phi/F",
    "nGenJet/I",
    "GenJet[pt/F,eta/F,phi/F,partonFlavour/I,hadronFlavour/I,mass/F]",
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
        name = "N_GenJets_b",
        texX = 'Number of b GenJets', texY = 'Number of Events',
        attribute = lambda event, sample: event.nbjet,
        binning=[11,-0.5,10.5],
    ))

    plots.append(Plot(
        name = "N_GenJets_light",
        texX = 'Number of light GenJets', texY = 'Number of Events',
        attribute = lambda event, sample: event.nlightjet,
        binning=[11,-0.5,10.5],
    ))

    plots.append(Plot(
        name = 'isMatchable',
        texX = 'top decay is matchable', texY = 'Number of Events',
        attribute = lambda event, sample: event.ismatchable,
        binning=[2,-0.5,1.5],
    ))

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
        name = 'm_tophad_matched',
        texX = 'm(top_{had}) after matching', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtophad_matched,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_toplep_matched',
        texX = 'm(top_{lep}) after matching', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtoplep_matched,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_Whad_matched',
        texX = 'm(W_{had}) after matching', texY = 'Number of Events',
        attribute = lambda event, sample: event.mWhad_matched,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'm_tophad_reduced_matched',
        texX = 'm(reduced top_{had}) after matching', texY = 'Number of Events',
        attribute = lambda event, sample: event.mtophadreduced_matched,
        binning=[40, 0.0, 400],
    ))

    plots.append(Plot(
        name = 'btagscore_matched',
        texX = 'b tag score after matching', texY = 'Number of Events',
        attribute = lambda event, sample: event.btagscore_matched,
        binning=[25, 0.0, 1.0],
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
        binning=[5, -0.5, 4.5],
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
        binning=[8, 0.5, 8.5],
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

    plots.append(Plot(
        name = 'blep_disc',
        texX = 'b tag score of b_{lep} jet', texY = 'Number of Events',
        attribute = lambda event, sample: event.blep_disc,
        binning=[40, 0.0, 1.0],
    ))

    plots.append(Plot(
        name = 'bhad_disc',
        texX = 'b tag score of b_{had} jet', texY = 'Number of Events',
        attribute = lambda event, sample: event.bhad_disc,
        binning=[40, 0.0, 1.0],
    ))

    plotting.fill(plots, read_variables = read_variables, sequence = sequence)

    ################################################################################
    drawPlots(plots, mode, args.shape)
    allPlots[mode] = plots

    for plot in plots:
      if "bmatch" in plot.name:
        for i, l in enumerate(plot.histos):
          for j, h in enumerate(l):
            h.GetXaxis().SetBinLabel(1, "no match")
            h.GetXaxis().SetBinLabel(2, "b_{lep}")
            h.GetXaxis().SetBinLabel(3, "b_{had}")
            h.GetXaxis().SetBinLabel(4, "b_{lep}+b_{had}")
            h.GetXaxis().SetBinLabel(5, " ")
      if "category" in plot.name:
        for i, l in enumerate(plot.histos):
          for j, h in enumerate(l):
            h.GetXaxis().SetBinLabel(1, "m_{had}")
            h.GetXaxis().SetBinLabel(2, "u_{had}")
            h.GetXaxis().SetBinLabel(3, "n_{had}")
            h.GetXaxis().SetBinLabel(4, "e_{had}")
            h.GetXaxis().SetBinLabel(5, "m_{lep}")
            h.GetXaxis().SetBinLabel(6, "u_{lep}")
            h.GetXaxis().SetBinLabel(7, "n_{lep}")
            h.GetXaxis().SetBinLabel(8, "e_{lep}")

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
