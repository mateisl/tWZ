#!/usr/bin/env python

# Standard imports
from operator                   import attrgetter
from math                       import pi, sqrt, cosh, cos, acos, exp
import ROOT, os

# RootTools
from RootTools.core.standard     import *

# helpers
from Analysis.Tools.helpers              import deltaPhi, deltaR, getObjDict
from tWZ.Tools.helpers          import getCollection


from tWZ.Tools.objectSelection  import isBJet, isAnalysisJet
from Analysis.Tools.WeightInfo       import WeightInfo

import logging
logger = logging.getLogger(__name__)

from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons

jetVars          = ['pt/F', 'eta/F', 'phi/F', 'btagDeepB/F', 'jetId/I', 'btagDeepFlavB/F', 'mass/F']

jetVarNames      = [x.split('/')[0] for x in jetVars]

lstm_jets_maxN   = 10
lstm_jetVars     = ['pt/F', 'eta/F', 'phi/F', 'btagDeepFlavB/F', 'btagDeepFlavC/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'muEF/F', 'puId/F', 'qgl/F']
lstm_jetVarNames = [x.split('/')[0] for x in lstm_jetVars]

lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','cutBased/I','miniPFRelIso_all/F','pfRelIso03_all/F','mvaFall17V2Iso_WP90/O', 'mvaTOP/F', 'sip3d/F','lostHits/I','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','eleIndex/I','muIndex/I']
lepVarNames      = [x.split('/')[0] for x in lepVars]

# Training variables
read_variables = [\
    "year/I",
    "nBTag/I",
    "nJetGood/I",
    "met_pt/F", "met_phi/F",
    "lep[pt/F,eta/F,phi/F]",
    "Z1_pt/F", "Z1_eta/F", "Z1_phi/F", "Z1_mass/F", "Z1_cosThetaStar/F",
    "nonZ1_l1_index/I",
    "Jet[%s]"%(",".join(jetVars)),
    "nJet/I",
    "JetGood[pt/F,eta/F,phi/F,btagDeepB/F,index/I,btagDeepFlavB/F]",
    "nlep/I",
    "Z1_lldPhi/F",
    "Z1_lldR/F",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I",
    "l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
    "l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
]

b_tagger = "DeepJet"


## Some Functions ##############################################################
def getBJetindex( event ):
    maxscore = 0.0
    index = -1
    for i in range(event.nJetGood):
        btagscore = event.JetGood_btagDeepFlavB[i]
        if btagscore > maxscore:
            maxscore = btagscore
            index = i
    return index

## top reco functions ##############################################################
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
    if Njetsmax < 4: return[]
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
                    jet_permutations.append([[jet_i, jet_j, jet_k, jet_l], [event.JetGood_btagDeepFlavB[i],event.JetGood_btagDeepFlavB[j]]])

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

## Sequence ####################################################################
sequence = []

def getM3l( event, sample=None):
    # get the invariant mass of the 3l system
    l = []
    for i in range(3):
        l.append(ROOT.TLorentzVector())
        l[i].SetPtEtaPhiM(event.lep_pt[i], event.lep_eta[i], event.lep_phi[i],0)
    event.m3l = (l[0] + l[1] + l[2]).M()
sequence.append( getM3l )

from tWZ.Tools.objectSelection import isBJet
def make_jets( event, sample ):
    event.jets     = [getObjDict(event, 'JetGood_', jetVarNames, i) for i in range(int(event.nJetGood))]
sequence.append( make_jets )

def getAngles(event, sample=None):
    event.nonZ1_l1_Z1_deltaPhi = deltaPhi(event.lep_phi[event.nonZ1_l1_index], event.Z1_phi)
    event.Z1_j1_deltaPhi       = deltaPhi(event.Z1_phi, event.JetGood_phi[0])
    event.nonZ1_l1_Z1_deltaEta = abs(event.lep_eta[event.nonZ1_l1_index] - event.Z1_eta)
    event.nonZ1_l1_Z1_deltaR   = deltaR({'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet0_Z1_deltaR       = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet0_nonZ1_l1_deltaR = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    event.jet1_Z1_deltaR       = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet1_nonZ1_l1_deltaR = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    event.jet2_Z1_deltaR       = deltaR({'eta':event.JetGood_eta[2], 'phi':event.JetGood_phi[2]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet2_nonZ1_l1_deltaR = deltaR({'eta':event.JetGood_eta[2], 'phi':event.JetGood_phi[2]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    i_bjet = getBJetindex(event)
    event.bJet_Z1_deltaR      = deltaR({'eta':event.JetGood_eta[i_bjet], 'phi':event.JetGood_phi[i_bjet]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.bJet_nonZ1l1_deltaR = deltaR({'eta':event.JetGood_eta[i_bjet], 'phi':event.JetGood_phi[i_bjet]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
sequence.append( getAngles )


def forwardJets( event, sample=None ):
    #jets einlesen (in case of MC also reat the index of the genjet)
    alljets   = getCollection( event, 'Jet', jetVarNames, 'nJet')
    alljets.sort( key = lambda j: -j['pt'] )
    leptons   = getCollection(event, "lep", lepVarNames, 'nlep')
    # clean against good leptons
    clean_jets,_ = cleanJetsAndLeptons( alljets, leptons )
    # filter pt, but not eta
    jets_no_eta         = filter(lambda j:j['pt']>30, clean_jets)
    if jets_no_eta:
        event.maxAbsEta_of_pt30jets = max( [ abs(j['eta']) for j in jets_no_eta ])
    else:
        event.maxAbsEta_of_pt30jets = -1

sequence.append( forwardJets )


## Top Reco sequence ###########################################################
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

    LeptonNoZ = ROOT.TLorentzVector()
    LeptonNoZ.SetPtEtaPhiM(event.lep_pt[event.nonZ1_l1_index], event.lep_eta[event.nonZ1_l1_index], event.lep_phi[event.nonZ1_l1_index], 0)

    # observables for mtop1 and mtop2 closest to mtop
    event.mtop_average = (hypo_selected['toplep'].M()+hypo_selected['tophad'].M())/2 if foundhypo else -1
    event.mtoplep = hypo_selected['toplep'].M() if foundhypo else -1
    event.mtophad = hypo_selected['tophad'].M() if foundhypo else -1
    event.mtop_hi = mtop_hi if foundhypo else -1
    event.mtop_lo = mtop_lo if foundhypo else -1
    event.mtop_diff = (mtop_hi-mtop_lo)/(mtop_hi+mtop_lo) if foundhypo else -1
    event.pt_diff = abs(hypo_selected['toplep'].Pt()-hypo_selected['tophad'].Pt()) if foundhypo else -1
    event.mW_lep = hypo_selected['Wlep'].M() if foundhypo else -1
    event.mW_had = hypo_selected['Whad'].M() if foundhypo else -1
    event.ptW_lep = hypo_selected['Wlep'].Pt() if foundhypo else -1
    event.ptW_had = hypo_selected['Whad'].Pt() if foundhypo else -1
    event.chi2 = hypo_selected['chi2'] if foundhypo else -1
    event.pgof = exp(-0.5*hypo_selected['chi2']) if foundhypo else -1
    event.hypofound = 1 if foundhypo else 0
    event.dR_tops = hypo_selected['toplep'].DeltaR(hypo_selected['tophad']) if foundhypo else -1
    event.dR_Ws = hypo_selected['Wlep'].DeltaR(hypo_selected['Whad']) if foundhypo else -1
    event.dR_bottoms = hypo_selected['bhad'].DeltaR(hypo_selected['blep']) if foundhypo else -1
    event.dR_toplep_Z = hypo_selected['toplep'].DeltaR(Z) if foundhypo else -1
    event.dR_tophad_Z = hypo_selected['tophad'].DeltaR(Z) if foundhypo else -1
    event.dR_toplep_LepNoZ = hypo_selected['toplep'].DeltaR(LeptonNoZ) if foundhypo else -1
    event.dR_tophad_LepNoZ = hypo_selected['tophad'].DeltaR(LeptonNoZ) if foundhypo else -1
    event.dR_Wlep_LepNoZ = hypo_selected['Wlep'].DeltaR(LeptonNoZ) if foundhypo else -1
    event.dR_Whad_LepNoZ = hypo_selected['Whad'].DeltaR(LeptonNoZ) if foundhypo else -1
    event.dR_blep_LepNoZ = hypo_selected['blep'].DeltaR(LeptonNoZ) if foundhypo else -1
    event.dR_bhad_LepNoZ = hypo_selected['bhad'].DeltaR(LeptonNoZ) if foundhypo else -1
    event.blep_disc = hypo_selected['blep_disc'] if foundhypo else -1
    event.bhad_disc = hypo_selected['bhad_disc'] if foundhypo else -1
sequence.append(TopReco)
################################################################################

all_mva_variables = {

# global event properties
     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
     "mva_met_pt"                :(lambda event, sample: event.met_pt),
     "mva_m3l"                   :(lambda event, sample: event.m3l),
     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
     "mva_nBTag"                 :(lambda event, sample: event.nBTag),

# top reconstruction
     "mva_mtop_average"          :(lambda event, sample: event.mtop_average),
     "mva_mtoplep"               :(lambda event, sample: event.mtoplep),
     "mva_mtophad"               :(lambda event, sample: event.mtophad),
     "mva_mtop_hi"               :(lambda event, sample: event.mtop_hi),
     "mva_mtop_lo"               :(lambda event, sample: event.mtop_lo),
     "mva_mtop_diff"             :(lambda event, sample: event.mtop_diff),
     "mva_mW_lep"                :(lambda event, sample: event.mW_lep),
     "mva_mW_had"                :(lambda event, sample: event.mW_had),
     "mva_ptW_lep"               :(lambda event, sample: event.ptW_lep),
     "mva_ptW_had"               :(lambda event, sample: event.ptW_had),
     "mva_chi2"                  :(lambda event, sample: event.chi2),
     "mva_pgof"                  :(lambda event, sample: event.pgof),
     "mva_dR_tops"               :(lambda event, sample: event.dR_tops),
     "mva_dR_toplep_Z"           :(lambda event, sample: event.dR_toplep_Z),
     "mva_dR_tophad_Z"           :(lambda event, sample: event.dR_tophad_Z),
     "mva_dR_Ws"                 :(lambda event, sample: event.dR_Ws),
     "mva_dR_bottoms"            :(lambda event, sample: event.dR_bottoms),
     "mva_dR_toplep_LepNoZ"      :(lambda event, sample: event.dR_toplep_LepNoZ),
     "mva_dR_tophad_LepNoZ"      :(lambda event, sample: event.dR_tophad_LepNoZ),
     "mva_dR_Wlep_LepNoZ"        :(lambda event, sample: event.dR_Wlep_LepNoZ),
     "mva_dR_Whad_LepNoZ"        :(lambda event, sample: event.dR_Whad_LepNoZ),
     "mva_dR_blep_LepNoZ"        :(lambda event, sample: event.dR_blep_LepNoZ),
     "mva_dR_bhad_LepNoZ"        :(lambda event, sample: event.dR_bhad_LepNoZ),



# jet kinmatics
     "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
     "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
     "mva_jet0_btagDeepFlavB"    :(lambda event, sample: event.JetGood_btagDeepFlavB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepFlavB[0]>-10) else -10),
     "mva_jet1_pt"               :(lambda event, sample: event.JetGood_pt[1]          if event.nJetGood >=2 else 0),
     "mva_jet1_eta"              :(lambda event, sample: event.JetGood_eta[1]         if event.nJetGood >=2 else -10),
     "mva_jet1_btagDeepFlavB"    :(lambda event, sample: event.JetGood_btagDeepFlavB[1] if (event.nJetGood >=2 and event.JetGood_btagDeepFlavB[1]>-10) else -10),
     "mva_jet2_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=3 else 0),
     "mva_jet2_eta"              :(lambda event, sample: event.JetGood_eta[2]         if event.nJetGood >=3 else -10),
     "mva_jet2_btagDeepFlavB"    :(lambda event, sample: event.JetGood_btagDeepFlavB[2]   if (event.nJetGood >=3 and event.JetGood_btagDeepFlavB[1]>-10) else -10),

# Z1 kinematics
     "mva_Z1_pt"                 :(lambda event, sample: event.Z1_pt),
     "mva_Z1_eta"                :(lambda event, sample: event.Z1_eta),
     "mva_Z1_cosThetaStar"       :(lambda event, sample: event.Z1_cosThetaStar),

# extra lepton kinematics
     "mva_lnonZ1_pt"             :(lambda event, sample: event.lep_pt[event.nonZ1_l1_index]),
     "mva_lnonZ1_eta"            :(lambda event, sample: event.lep_eta[event.nonZ1_l1_index]),

# Z1 vs. other objects
     "mva_nonZ1_l1_Z1_deltaPhi"  :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nlep >= 1 else -1 ),
     "mva_nonZ1_l1_Z1_deltaR"    :(lambda event, sample: event.nonZ1_l1_Z1_deltaR),

     "mva_jet0_Z1_deltaR"        :(lambda event, sample: event.jet0_Z1_deltaR         if event.nJetGood >=1 else -1),
     "mva_jet1_Z1_deltaR"        :(lambda event, sample: event.jet1_Z1_deltaR         if event.nJetGood >=2 else -1),
     "mva_jet2_Z1_deltaR"        :(lambda event, sample: event.jet2_Z1_deltaR         if event.nJetGood >=3 else -1),

# nonZ1_l1 vs. other objects
     "mva_jet0_nonZl1_deltaR"    :(lambda event, sample: event.jet0_nonZ1_l1_deltaR    if event.nJetGood >=1 else -1),
     "mva_jet1_nonZl1_deltaR"    :(lambda event, sample: event.jet1_nonZ1_l1_deltaR    if event.nJetGood >=2 else -1),
     "mva_bJet_Z1_deltaR"        :(lambda event, sample: event.bJet_Z1_deltaR),
     "mva_bJet_non_Z1l1_deltaR"  :(lambda event, sample: event.bJet_nonZ1l1_deltaR),
     "mva_maxAbsEta_of_pt30jets" :(lambda event, sample: event.maxAbsEta_of_pt30jets),

#leptonMVA
     "mva_l1_mvaTOP"             :(lambda event, sample: event.l1_mvaTOP),
     "mva_l2_mvaTOP"             :(lambda event, sample: event.l2_mvaTOP),
     "mva_l3_mvaTOP"             :(lambda event, sample: event.l3_mvaTOP),

     "mva_l1_mvaTOPWP"           :(lambda event, sample: event.l1_mvaTOPWP),
     "mva_l2_mvaTOPWP"           :(lambda event, sample: event.l2_mvaTOPWP),
     "mva_l3_mvaTOPWP"           :(lambda event, sample: event.l3_mvaTOPWP),


}

def lstm_jets(event, sample):
    jets = [ getObjDict( event, 'Jet_', lstm_jetVarNames, event.JetGood_index[i] ) for i in range(int(event.nJetGood)) ]
    #jets = filter( jet_vector_var['selector'], jets )
    return jets

# for the filler
mva_vector_variables    =   {
    #"mva_Jet":  {"name":"Jet", "vars":lstm_jetVars, "varnames":lstm_jetVarNames, "selector": (lambda jet: True), 'maxN':10}
    "mva_Jet":  {"func":lstm_jets, "name":"Jet", "vars":lstm_jetVars, "varnames":lstm_jetVarNames}
}

## Using all variables
mva_variables_ = all_mva_variables.keys()
mva_variables_.sort()
mva_variables  = [ (key, value) for key, value in all_mva_variables.iteritems() if key in mva_variables_ ]

import numpy as np
import operator

# make predictions to be used with keras.predict
def predict_inputs( event, sample, jet_lstm = False):
    flat_variables = np.array([[getattr( event, mva_variable) for mva_variable, _ in mva_variables]])
    if jet_lstm:
        lstm_jets_maxN = 10 #remove after retraining
        jet_vector_var = mva_vector_variables["mva_Jet"]
        jets = mva_vector_variables["mva_Jet"]["func"](event,sample=None)
        jets =  [ [ operator.itemgetter(varname)(jet) for varname in lstm_jetVarNames] for jet in jets[:lstm_jets_maxN] ]
        # zero padding
        jets += [ [0.]*len(lstm_jetVarNames)]*(max(0, lstm_jets_maxN-len(jets)))
        jets = np.array([jets])

        return [ flat_variables, jets ]
    else:
        return   flat_variables

#define training samples for multiclassification
from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed import *
#use only Summer16
training_samples = [ Summer16.TWZ_NLO_DR, Summer16.TTZ ]

assert len(training_samples)==len(set([s.name for s in training_samples])), "training_samples names are not unique!"

# training selection
from tWZ.Tools.cutInterpreter import cutInterpreter
# selectionString = cutInterpreter.cutString( 'trilepT-onZ1-btag1-njet3p' )
selectionString = cutInterpreter.cutString( 'trilepT-minDLmass12-onZ1-njet4p-deepjet1' )
