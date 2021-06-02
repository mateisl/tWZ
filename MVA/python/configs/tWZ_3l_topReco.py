#!/usr/bin/env python

# Standard imports
from operator                   import attrgetter
from math                       import pi, sqrt, cosh, cos, acos
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
    "JetGood[pt/F,eta/F,phi/F,btagDeepB/F,index/I]",
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
        btagscore = event.JetGood_btagDeepB[i]
        if btagscore > maxscore:
            maxscore = btagscore
            index = i
    return index

def getWhad( event ):
    # get Whad from min dijet mass
    # exclude highest b tag from jet selection
    index_btag = getBJetindex(event)
    min_mass = 1000.
    Whad = ROOT.TLorentzVector()
    for i in range(event.nJetGood):
        if i == index_btag:
            continue
        jetindex1 = event.JetGood_index[i]
        jet1  = ROOT.TLorentzVector()
        jet1.SetPtEtaPhiM(event.Jet_pt[jetindex1], event.Jet_eta[jetindex1], event.Jet_phi[jetindex1], event.Jet_mass[jetindex1])

        for j in range(event.nJetGood):
            if j == index_btag or j == i:
                continue
            jetindex2 = event.JetGood_index[j]
            jet2  = ROOT.TLorentzVector()
            jet2.SetPtEtaPhiM(event.Jet_pt[jetindex2], event.Jet_eta[jetindex2], event.Jet_phi[jetindex2], event.Jet_mass[jetindex2])

            dijet = jet1+jet2
            if dijet.M() < min_mass:
                min_mass = dijet.M()
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
    if discriminant <= 0:
        # Take only real part of the solution for pz:
        neutrino = ROOT.TLorentzVector()
        neutrino.SetPxPyPzE(met.Px(),met.Py(),-B / A,0)
        neutrino.SetE(neutrino.P())
    else:
        discriminant = sqrt(discriminant)
        neutrino1 = ROOT.TLorentzVector()
        neutrino1.SetPxPyPzE(met.Px(),met.Py(),(-B - discriminant) / A,0)
        neutrino1.SetE(neutrino1.P())
        neutrino2 = ROOT.TLorentzVector()
        neutrino2.SetPxPyPzE(met.Px(),met.Py(),(-B + discriminant) / A,0)
        neutrino2.SetE(neutrino2.P())
        if neutrino1.E() > neutrino2.E():
            neutrino = neutrino1
        else:
            neutrino = neutrino2

    Wlep = lepton + neutrino
    return Wlep

def getTopHypos( event ):
    Wh = getWhad( event )
    Wl = getWlep( event )
    bjetindex = event.JetGood_index[getBJetindex(event)]
    bJet = ROOT.TLorentzVector()
    bJet.SetPtEtaPhiM(event.Jet_pt[bjetindex], event.Jet_eta[bjetindex], event.Jet_phi[bjetindex], event.Jet_mass[bjetindex])

    toph = Wh + bJet
    topl = Wl + bJet

    return toph, topl

def get2TopHypos( event ):
    topA1 = ROOT.TLorentzVector()
    topA2 = ROOT.TLorentzVector()
    topB1 = ROOT.TLorentzVector()
    topB2 = ROOT.TLorentzVector()
    # Get Ws
    Wh = getWhad( event )
    Wl = getWlep( event )
    # Get highest btag
    bjetindex = event.JetGood_index[getBJetindex(event)]
    bJet1 = ROOT.TLorentzVector()
    bJet1.SetPtEtaPhiM(event.Jet_pt[bjetindex], event.Jet_eta[bjetindex], event.Jet_phi[bjetindex], event.Jet_mass[bjetindex])

    # now search for second highest bjet
    secondmax = 0
    index = -1
    bJet2 = ROOT.TLorentzVector()
    for i in range(event.nJetGood):
        if i == bjetindex:
            continue
        btagscore = event.JetGood_btagDeepB[i]
        if btagscore > secondmax:
            secondmax = btagscore
            index = i
    if index != -1:
        bjet2index = event.JetGood_index[index]
        bJet2.SetPtEtaPhiM(event.Jet_pt[bjet2index], event.Jet_eta[bjet2index], event.Jet_phi[bjet2index], event.Jet_mass[bjet2index])
    else:
        return topA1, topA2

    # construct tops
    topA1 = Wh + bJet1
    topA2 = Wl + bJet2
    topB1 = Wh + bJet2
    topB2 = Wl + bJet1

    # get hypothesis closest to mtop
    mtop = 172.5
    diffA = pow(topA1.M()-mtop,2) + pow(topA2.M()-mtop,2)
    diffB = pow(topB1.M()-mtop,2) + pow(topB2.M()-mtop,2)

    if diffA < diffB:
        return topA1, topA2
    else:
        return topB1, topB2


## Sequence ####################################################################
sequence = []

def getTTbar( event, sample):
    top1, top2  = get2TopHypos(event)
    mtop = 172.5
    if abs(top1.M()-172.5) < abs(top2.M()-172.5):
        event.MTop1 = top1.M()
        event.MTop2 = top2.M()
    else:
        event.MTop1 = top2.M()
        event.MTop2 = top1.M()

sequence.append( getTTbar )

def getWhadMass( event, sample):
    event.WhadMass = getWhad( event ).M()

sequence.append( getWhadMass )

def getWpt( event, sample=None):
    # get the lepton and met
    lepton  = ROOT.TLorentzVector()
    met     = ROOT.TLorentzVector()
    lepton.SetPtEtaPhiM(event.lep_pt[event.nonZ1_l1_index], event.lep_eta[event.nonZ1_l1_index], event.lep_phi[event.nonZ1_l1_index], 0)
    met.SetPtEtaPhiM(event.met_pt, 0, event.met_phi, 0)
    # get the W boson candidate
    W   = lepton + met
    event.W_pt = W.Pt()
sequence.append( getWpt )

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
    event.bJets    = filter(lambda j:isBJet(j, year=event.year) and abs(j['eta'])<=2.4    , event.jets)
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


def getDeltaMaxEta(event,sample):
    jet1 = ROOT.TLorentzVector()
    # First get jet with max eta
    maxeta = 0
    i_maxeta = -1
    for i in range(event.nJet):
        if abs(event.Jet_eta[i]) > maxeta and event.Jet_pt[i] > 30:
            maxeta = abs(event.Jet_eta[i])
            found_jet1 = True
            i_maxeta = i
    jet1.SetPtEtaPhiM(event.Jet_pt[i_maxeta], event.Jet_eta[i_maxeta], event.Jet_phi[i_maxeta], event.Jet_mass[i_maxeta])
    # get second jet which is not the btag and gives max m_{ij}
    i_bjet = getBJetindex(event)
    maxdijetmass = 0
    i_maxdijetmass = -1
    for i in range(event.nJetGood):
        if i==i_bjet: continue
        jetindex = event.JetGood_index[i]
        jetcandidate = ROOT.TLorentzVector()
        jetcandidate.SetPtEtaPhiM(event.Jet_pt[jetindex], event.Jet_eta[jetindex], event.Jet_phi[jetindex], event.Jet_mass[jetindex])
        dijet = jet1+jetcandidate
        if dijet.M() > maxdijetmass:
            maxdijetmass=dijet.M()
            i_maxdijetmass=jetindex

    if i_maxeta!=-1 and i_maxdijetmass!=-1:
        event.deltamaxeta = abs(event.Jet_eta[i_maxdijetmass]-maxeta)
    else:
        event.deltamaxeta = float('nan')

sequence.append( getDeltaMaxEta )

def getDRTopZ( event, sample ):
    tophad, toplep = getTopHypos( event )

    if abs(tophad.M()-172.5) < abs(toplep.M()-172.5):
        top = tophad
    else:
        top = toplep

    Z = ROOT.TLorentzVector()
    Z.SetPtEtaPhiM(event.Z1_pt, event.Z1_eta, event.Z1_phi, event.Z1_mass)
    event.DRTopZ = top.DeltaR(Z)


sequence.append( getDRTopZ )

def get2TopChi2(event,sample):
    tophad, toplep = getTopHypos( event )
    # Values from UHH2
    Mtlep_mean  = 174.
    Mtlep_sigma =  18.
    Mthad_mean  = 181.
    Mthad_sigma =  15.
    chi2 = pow((tophad.M()-Mthad_mean)/Mthad_sigma,2)+pow((toplep.M()-Mtlep_mean)/Mtlep_sigma,2)
    event.chi2 = chi2
sequence.append( get2TopChi2 )

all_mva_variables = {

# global event properties
     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
     "mva_met_pt"                :(lambda event, sample: event.met_pt),
     "mva_m3l"                   :(lambda event, sample: event.m3l),
     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
     "mva_nBTag"                 :(lambda event, sample: event.nBTag),

# top reconstruction
     "mva_mtop1"                 :(lambda event, sample: event.MTop1),
     "mva_mtop2"                 :(lambda event, sample: event.MTop1),
     "mva_mWhad"                 :(lambda event, sample: event.WhadMass),
     "mva_dR_top_Z"              :(lambda event, sample: event.DRTopZ),
     "mva_chi2"                  :(lambda event, sample: event.chi2),

# jet kinmatics
     "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
     "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),
     "mva_jet1_pt"               :(lambda event, sample: event.JetGood_pt[1]          if event.nJetGood >=2 else 0),
     "mva_jet1_eta"              :(lambda event, sample: event.JetGood_eta[1]         if event.nJetGood >=2 else -10),
     "mva_jet1_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[1] if (event.nJetGood >=2 and event.JetGood_btagDeepB[1]>-10) else -10),
     "mva_jet2_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=3 else 0),
     "mva_jet2_eta"              :(lambda event, sample: event.JetGood_eta[2]         if event.nJetGood >=3 else -10),
     "mva_jet2_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[2]   if (event.nJetGood >=3 and event.JetGood_btagDeepB[1]>-10) else -10),

# Z1 kinematics
     "mva_Z1_pt"                 :(lambda event, sample: event.Z1_pt),
     "mva_Z1_eta"                :(lambda event, sample: event.Z1_eta),
     "mva_Z1_cosThetaStar"       :(lambda event, sample: event.Z1_cosThetaStar),

# extra lepton kinematics
     "mva_lnonZ1_pt"             :(lambda event, sample: event.lep_pt[event.nonZ1_l1_index]),
     "mva_lnonZ1_eta"            :(lambda event, sample: event.lep_eta[event.nonZ1_l1_index]),

# leptonic W
     "mva_W_pt"                  :(lambda event, sample: event.W_pt),

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

# Delta eta
     "mva_deltamaxeta"           :(lambda event, sample: event.deltamaxeta),

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
selectionString = cutInterpreter.cutString( 'trilepT-onZ1-btag1-njet3p' )
# selectionString = cutInterpreter.cutString( 'trilepT-minDLmass12-onZ1-njet4p-btag1' )
