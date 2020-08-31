#!/usr/bin/env python

# Standard imports
from operator                   import attrgetter
from math                       import pi, sqrt, cosh
import ROOT, os

# RootTools
from RootTools.core.standard     import *

# helpers
from tWZ.Tools.helpers          import deltaPhi, deltaR2, deltaR, getCollection, getObjDict
from tWZ.Tools.objectSelection  import isBJet, isAnalysisJet

# Logger
import tWZ.Tools.logger as logger
logger = logger.get_logger("INFO", logFile = None )

from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--variables',          action='store', type=str,   default='all')
argParser.add_argument('--overwrite',          action='store_true')
argParser.add_argument('--mva',                action='store', type=str,   default='mlp')
argParser.add_argument('--NTrees',             action='store', type=int, default=250)
argParser.add_argument('--maxdepth',           action='store', type=int, default=1)
argParser.add_argument('--ncuts',              action='store', type=int, default=50)
argParser.add_argument('--layer',              action='store', type=str, default='N')
argParser.add_argument('--sampling',           action='store', type=float, default=1)
argParser.add_argument('--epoch',              action='store', type=float, default=1)

args = argParser.parse_args()

jetVars          = ['pt/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'rawFactor/F', 'eta/F', 'phi/F', 'jetId/I', 'btagDeepB/F', 'btagDeepFlavB/F', 'btagCSVV2/F', 'area/F']
jetVarNames      = [x.split('/')[0] for x in jetVars]

lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','cutBased/I','miniPFRelIso_all/F','pfRelIso03_all/F','mvaFall17V2Iso_WP90/O', 'mvaTTH/F', 'sip3d/F','lostHits/I','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','eleIndex/I','muIndex/I']
lepVarNames      = [x.split('/')[0] for x in lepVars]

# Training variables
read_variables = [\
                    "nBTag/I",
                    "nJetGood/I",
                    "met_pt/F", "met_phi/F",
                    "lep[pt/F,eta/F,phi/F]", 
                    "Z1_pt/F", "Z1_eta/F", "Z1_phi/F", "Z1_mass/F", "Z1_cosThetaStar/F",
                    "nonZ1_l1_index/I",
                    "Jet[%s]"%(",".join(jetVars)),
                    "nJet/I",
                    "JetGood[pt/F,eta/F,phi/F,btagDeepB/F]",
                    "nlep/I",
                    "Z1_lldPhi/F",
                    "Z1_lldR/F",
                    ]

#b tagger
b_tagger = "DeepCSV"
#def flavorBin( event, sample=None):
#    event.flavorBin = 0
#
#    if      event.nMuons_tight_3l==3 and event.nElectrons_tight_3l==0: event.flavorBin = 1
#    elif    event.nMuons_tight_3l==2 and event.nElectrons_tight_3l==1: event.flavorBin = 2
#    elif    event.nMuons_tight_3l==1 and event.nElectrons_tight_3l==2: event.flavorBin = 3
#    elif    event.nMuons_tight_3l==0 and event.nElectrons_tight_3l==3: event.flavorBin = 4
#sequence.append( flavorBin )

# sequence 
sequence = []

def getbjets( event, sample=None ):
#    #bjets filtern( nur 2016 )
    alljets   = getCollection( event, 'Jet', jetVarNames , 'nJet')  
    goodjets = filter( isAnalysisJet, alljets ) 
    bJets = filter(lambda j:isBJet(j, tagger=b_tagger, year=2016) and abs(j['eta'])<=2.4, goodjets)
    for i in range(len(bJets)) :
        event.bjet_pt  = bJets[i]['pt']
        event.bjet_eta = bJets[i]['eta']
        event.bjet_phi = bJets[i]['phi']

    event.bjet_Z1_deltaR      = deltaR({'eta':event.bjet_eta, 'phi':event.bjet_phi}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.bjet_nonZ1l1_deltaR = deltaR({'eta':event.bjet_eta, 'phi':event.bjet_phi}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
sequence.append( getbjets )

def getM3l( event, sample=None):
    # get the invariant mass of the 3l system
    l = []
    for i in range(3):
        l.append(ROOT.TLorentzVector())
        l[i].SetPtEtaPhiM(event.lep_pt[i], event.lep_eta[i], event.lep_phi[i],0)
    event.M3l = (l[0] + l[1] + l[2]).M()
sequence.append( getM3l )

def getDeltaPhi(event, sample=None):
    event.nonZ1_l1_Z1_deltaPhi = deltaPhi(event.lep_phi[event.nonZ1_l1_index], event.Z1_phi)
    event.Z1_j1_deltaPhi       = deltaPhi(event.Z1_phi, event.JetGood_phi[0]) 
sequence.append( getDeltaPhi )

def getDeltaEta(event, sample=None):
    event.nonZ1_l1_Z1_deltaEta = abs(event.lep_eta[event.nonZ1_l1_index] - event.Z1_eta)
sequence.append( getDeltaEta )

def getDeltaR(event, sample=None):
    event.nonZ1_l1_Z1_deltaR   = deltaR({'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet0_Z1_deltaR       = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet0_nonZ1_l1_deltaR = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    event.jet1_Z1_deltaR       = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet1_nonZ1_l1_deltaR = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    event.jet2_Z1_deltaR       = deltaR({'eta':event.JetGood_eta[2], 'phi':event.JetGood_phi[2]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet2_nonZ1_l1_deltaR = deltaR({'eta':event.JetGood_eta[2], 'phi':event.JetGood_phi[2]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})

    event.jet0_l1_deltaR        = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.lep_eta[1], 'phi':event.lep_phi[1]})
    event.jet0_l2_deltaR        = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.lep_eta[2], 'phi':event.lep_phi[2]})
    event.jet1_l1_deltaR        = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.lep_eta[1], 'phi':event.lep_phi[1]})
    event.jet1_l2_deltaR        = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.lep_eta[2], 'phi':event.lep_phi[2]})

    event.Z1_l1_deltaR          = deltaR({'eta':event.Z1_eta, 'phi':event.Z1_phi},{'eta':event.lep_eta[1], 'phi':event.lep_phi[1]})
sequence.append( getDeltaR )

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

def getjets( event, sample=None ):
    #jets einlesen (in case of MC also reat the index of the genjet)
    alljets   = getCollection( event, 'Jet', jetVarNames, 'nJet')  
    alljets.sort( key = lambda j: -j['pt'] )
    leptons   = getCollection(event, "lep", lepVarNames, 'nlep') 
    # clean against good leptons
    clean_jets,_ = cleanJetsAndLeptons( alljets, leptons )
    # filter pt, but not eta 
    jets_no_eta         = filter(lambda j:j['pt']>30, clean_jets)
    if jets_no_eta: 
        event.maxEta_of_pt30jets  = max( [ abs(j['eta']) for j in jets_no_eta ])   
sequence.append( getjets )

## met, ht, nonZ1_pt/eta, Z1_pt, nJet, nBTag, lep1_eta
#mva_variables =  {
#                    "met_pt"    :attrgetter("met_pt"), # copy
#                    "ht"        :attrgetter("ht"), # copy
#                    "lnonZ1_pt" :(lambda event: event.lep_pt[event.nonZ1_l1_index_4l]),
#                    "lnonZ1_eta":(lambda event: event.lep_eta[event.nonZ1_l1_index_4l]),
#                    "Z1_pt_4l"  :attrgetter("Z1_pt_4l"),
##                    "lep1_pt"   :(lambda event: event.lep_pt[0]),
##                    "lep2_pt"   :(lambda event: event.lep_pt[1]),
#                    "lep1_eta"  :(lambda event: event.lep_eta[0]),
##                    "lep2_eta"  :(lambda event: event.lep_eta[1]),
#                    "nJetGood":attrgetter("nJetGood"),
#                    "nBTag"     :attrgetter("nBTag"),      
##                    "yield"     :(lambda event: event.flavorBin),
##                    "jet1_pt"   :(lambda event: event.jet_pt[0]),
##                    "nLepLoose":(lambda event: event.nlep),
#                    #"myvar1" :(lambda event: event.nBTag), # calculate on the fly
#                    #"myvar2" :(lambda event: event.myFancyVar), # from sequence

if args.variables == 'tommy':
    mva_variables = {
                     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                     "mva_met_pt"                :(lambda event, sample: event.met_pt),
                     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                     "mva_nBTag"                 :(lambda event, sample: event.nBTag),

                     "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
                     "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
                     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),
                     "mva_jet1_pt"               :(lambda event, sample: event.JetGood_pt[1]          if event.nJetGood >=2 else 0),
                     "mva_jet1_eta"              :(lambda event, sample: event.JetGood_eta[1]         if event.nJetGood >=2 else -10),
                     "mva_jet1_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[1] if (event.nJetGood >=2 and event.JetGood_btagDeepB[1]>-10) else -10),
                     "mva_jet2_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=3 else 0),
                     "mva_jet2_eta"              :(lambda event, sample: event.JetGood_eta[2]         if event.nJetGood >=3 else -10),
                     "mva_jet2_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[2]   if (event.nJetGood >=3 and event.JetGood_btagDeepB[1]>-10) else -10),

                     "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
                     "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
                     "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),
                     "mva_Z1_mass"                :(lambda event, sample: event.Z1_mass),

                     "mva_jet0_Z1_deltaR"         :(lambda event, sample: event.jet0_Z1_deltaR         if event.nJetGood >=1 else -1),
                     "mva_jet1_Z1_deltaR"         :(lambda event, sample: event.jet1_Z1_deltaR         if event.nJetGood >=2 else -1),

                     "mva_jet0_l1_deltaR"         :(lambda event, sample: event.jet0_l1_deltaR         if event.nJetGood >=1 else -1),
                     "mva_jet1_l1_deltaR"         :(lambda event, sample: event.jet1_l1_deltaR         if event.nJetGood >=2 else -1),
                     "mva_Z1_l1_deltaR"           :(lambda event, sample: event.Z1_l1_deltaR),

                     "mva_l1_pt"                  :(lambda event, sample: event.lep_pt[1]),
                     "mva_l2_pt"                  :(lambda event, sample: event.lep_pt[2]),

                    }
if args.variables == 'alljets': 
    mva_variables = {
                     "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
                     "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
                     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),     
                     "mva_jet1_pt"               :(lambda event, sample: event.JetGood_pt[1]          if event.nJetGood >=2 else 0),
                     "mva_jet1_eta"              :(lambda event, sample: event.JetGood_eta[1]         if event.nJetGood >=2 else -10),
                     "mva_jet1_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[1] if (event.nJetGood >=2 and event.JetGood_btagDeepB[1]>-10) else -10),
                     "mva_jet2_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=3 else 0),
                     "mva_jet2_eta"              :(lambda event, sample: event.JetGood_eta[2]         if event.nJetGood >=3 else -10),
                     "mva_jet2_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[2]   if (event.nJetGood >=3 and event.JetGood_btagDeepB[1]>-10) else -10),
                     #bjets
                     "mva_bjet_pt"                :(lambda event, sample: event.bjet_pt),
                     "mva_bjet_Z1_deltaR"         :(lambda event, sample: event.bjet_Z1_deltaR),
                     "mva_bjet_non_Z1l1_deltaR"   :(lambda event, sample: event.bjet_nonZ1l1_deltaR),
                    }

elif args.variables == 'rand': 
      mva_variables = {
                      "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                      "mva_met_pt"                :(lambda event, sample: event.met_pt),
                      "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                      "mva_nBTag"                 :(lambda event, sample: event.nBTag),
#                      "mva_flavorBin"             :(lambda event, sample: event.flavorBin),
#                      "mva_nlep"                  :(lambda event, sample: event.nlep),                        
                      }

elif args.variables == 'Z1': 
      mva_variables = {
                      "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
                      "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
                      "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),
#                      "mva_Z1_mass"                :(lambda event, sample: event.Z1_mass),
                      }

elif args.variables == 'dRdE':
      mva_variables = {
#                       "mva_nonZl1_Z_deltaEta"     :(lambda event, sample: event.nonZl1_Z_deltaEta),
                       "mva_nonZ1_l1_Z1_deltaR"     :(lambda event, sample: event.nonZ1_l1_Z1_deltaR),

                       "mva_jet0_Z1_deltaR"         :(lambda event, sample: event.jet0_Z1_deltaR         if event.nJetGood >=1 else -1),
                       "mva_jet0_nonZl1_deltaR"     :(lambda event, sample: event.jet0_nonZ1_l1_deltaR    if event.nJetGood >=1 else -1),
                       "mva_jet1_Z1_deltaR"         :(lambda event, sample: event.jet1_Z1_deltaR         if event.nJetGood >=2 else -1),
                       "mva_jet1_nonZl1_deltaR"     :(lambda event, sample: event.jet1_nonZ1_l1_deltaR    if event.nJetGood >=2 else -1),
#                       "mva_jet2_Z1_deltaR"         :(lambda event, sample: event.jet2_Z1_deltaR         if event.nJetGood >=3 else -1),
#                       "mva_jet2_nonZl1_deltaR"    :(lambda event, sample: event.jet2_nonZl1_deltaR    if event.nJetGood >=3 else -1),
                       }

elif args.variables == 'randZWdphilnonZfwdjet' :
    mva_variables = {
                     "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
                     "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
#                     "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),

                     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                     "mva_met_pt"                :(lambda event, sample: event.met_pt),
                     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                     "mva_nBTag"                 :(lambda event, sample: event.nBTag),

                     "mva_W_pt"                   :(lambda event, sample: event.W_pt),

                     "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),
                     "mva_nonZ1_l1_Z1_deltaPhi"   :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nlep >= 2 else -1 ),
#                     "mva_Z_ll_deltaPhi"          :(lambda event, (TreeVariable.fromString( "Z1_lldPhi/F" ))),

                     "mva_lnonZ1_pt"              :(lambda event, sample: event.lep_pt[event.nonZ1_l1_index]),
                     "mva_lnonZ1_eta"             :(lambda event, sample: event.lep_eta[event.nonZ1_l1_index]),

                     "mva_jetsnoetacutptg30"      :(lambda event, sample: event.maxEta_of_pt30jets),
                    }

elif args.variables == 'randZWdphibjet' :
    mva_variables = {
                     "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
                     "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
                     "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),
                     
                     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                     "mva_met_pt"                :(lambda event, sample: event.met_pt),
                     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                     "mva_nBTag"                 :(lambda event, sample: event.nBTag),
                     
                     "mva_W_pt"                   :(lambda event, sample: event.W_pt),
                     
                     "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),
                     "mva_nonZ1_l1_Z1_deltaPhi"   :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nlep >= 2 else -1 ),

                    #bjets
                     "mva_bjet_pt"                :(lambda event, sample: event.bjet_pt),
                     "mva_bjet_Z1_deltaR"         :(lambda event, sample: event.bjet_Z1_deltaR),
                     "mva_bjet_non_Z1l1_deltaR"   :(lambda event, sample: event.bjet_nonZ1l1_deltaR),
                    }


elif args.variables == 'randZWdphijet' :
    mva_variables = {
                     "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
                     "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
                     "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),
                     
                     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                     "mva_met_pt"                :(lambda event, sample: event.met_pt),
                     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                     "mva_nBTag"                 :(lambda event, sample: event.nBTag),
                     
                     "mva_W_pt"                   :(lambda event, sample: event.W_pt),
                     
                     "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),
                     "mva_nonZ1_l1_Z1_deltaPhi"   :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nlep >= 2 else -1 ),

                     "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
                     "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
                     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),
                     "mva_jet1_pt"               :(lambda event, sample: event.JetGood_pt[1]          if event.nJetGood >=2 else 0),
                     "mva_jet1_eta"              :(lambda event, sample: event.JetGood_eta[1]         if event.nJetGood >=2 else -10),
                     "mva_jet1_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[1] if (event.nJetGood >=2 and event.JetGood_btagDeepB[1]>-10) else -10),
                     "mva_jet2_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=3 else 0),
                     "mva_jet2_eta"              :(lambda event, sample: event.JetGood_eta[2]         if event.nJetGood >=3 else -10),
                     "mva_jet2_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[2]   if (event.nJetGood >=3 and event.JetGood_btagDeepB[1]>-10) else -10),
                    }

elif args.variables == 'randZWdphijet0' :
    mva_variables = {
                     "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
                     "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
                     "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),

                     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                     "mva_met_pt"                :(lambda event, sample: event.met_pt),
                     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                     "mva_nBTag"                 :(lambda event, sample: event.nBTag),

                     "mva_W_pt"                   :(lambda event, sample: event.W_pt),

                     "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),
                     "mva_nonZ1_l1_Z1_deltaPhi"   :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nlep >= 2 else -1 ),

                     "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
                     "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
                     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),
                    }

elif args.variables == 'randZWdphidR' :
    mva_variables = {
                     "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
                     "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
                     "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),
                     
                     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                     "mva_met_pt"                :(lambda event, sample: event.met_pt),
                     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                     "mva_nBTag"                 :(lambda event, sample: event.nBTag),
                     
                     "mva_W_pt"                   :(lambda event, sample: event.W_pt),
                     
                     "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),
                     "mva_nonZ1_l1_Z1_deltaPhi"   :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nlep >= 2 else -1 ),
                     
                     "mva_jet0_Z1_deltaR"         :(lambda event, sample: event.jet0_Z1_deltaR         if event.nJetGood >=1 else -1),
                     "mva_jet0_nonZl1_deltaR"     :(lambda event, sample: event.jet0_nonZ1_l1_deltaR    if event.nJetGood >=1 else -1),
                     "mva_jet1_Z1_deltaR"         :(lambda event, sample: event.jet1_Z1_deltaR         if event.nJetGood >=2 else -1),
                     "mva_jet1_nonZl1_deltaR"     :(lambda event, sample: event.jet1_nonZ1_l1_deltaR    if event.nJetGood >=2 else -1),
                     "mva_jet2_Z1_deltaR"         :(lambda event, sample: event.jet2_Z1_deltaR         if event.nJetGood >=3 else -1),
                    }

elif args.variables == 'randZWdphi' :
    mva_variables = {
                     "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
                     "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
                     "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),
                     
                     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                     "mva_met_pt"                :(lambda event, sample: event.met_pt),
                     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                     "mva_nBTag"                 :(lambda event, sample: event.nBTag),
                     
                     "mva_W_pt"                   :(lambda event, sample: event.W_pt),

                     "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),
                     "mva_nonZ1_l1_Z1_deltaPhi"   :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nlep >= 2 else -1 ),
                     
#                     "mva_jet0_Z1_deltaR"         :(lambda event, sample: event.jet0_Z1_deltaR         if event.nJetGood >=1 else -1),
#                     "mva_jet0_nonZl1_deltaR"     :(lambda event, sample: event.jet0_nonZ1_l1_deltaR    if event.nJetGood >=1 else -1),
                    }

elif args.variables == 'randZWdphidRnonZ' :
    mva_variables = {
                     "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
                     "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
                     "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),
                     
                     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                     "mva_met_pt"                :(lambda event, sample: event.met_pt),
                     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                     "mva_nBTag"                 :(lambda event, sample: event.nBTag),
                     
                     "mva_W_pt"                   :(lambda event, sample: event.W_pt),

                     "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),
                     "mva_nonZ1_l1_Z1_deltaPhi"   :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nlep >= 2 else -1 ),
                     
                     "mva_jet0_nonZl1_deltaR"     :(lambda event, sample: event.jet0_nonZ1_l1_deltaR    if event.nJetGood >=1 else -1),
                     "mva_jet1_nonZl1_deltaR"     :(lambda event, sample: event.jet1_nonZ1_l1_deltaR    if event.nJetGood >=2 else -1),
                    }


elif args.variables == 'randZWdphidRnonZjet0bjet' :
    mva_variables = {
                     "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
                     "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
                     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),

                     "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
                     "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
                     "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),

                     "mva_jet0_nonZl1_deltaR"     :(lambda event, sample: event.jet0_nonZ1_l1_deltaR    if event.nJetGood >=1 else -1),
                     "mva_jet1_nonZl1_deltaR"     :(lambda event, sample: event.jet1_nonZ1_l1_deltaR    if event.nJetGood >=2 else -1),
#
                     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                     "mva_met_pt"                :(lambda event, sample: event.met_pt),
                     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                     "mva_nBTag"                 :(lambda event, sample: event.nBTag),

                     "mva_W_pt"                   :(lambda event, sample: event.W_pt),

                     "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),
                     "mva_nonZ1_l1_Z1_deltaPhi"   :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nlep >= 2 else -1 ),

                    #bjets
                     "mva_bjet_pt"                :(lambda event, sample: event.bjet_pt),
                     "mva_bjet_Z1_deltaR"         :(lambda event, sample: event.bjet_Z1_deltaR),
                     "mva_bjet_nonZ1l1_deltaR"   :(lambda event, sample: event.bjet_nonZ1l1_deltaR),
                    }

elif  args.variables == 'bjetstest':
    mva_variables = {
                    #bjets
                    "mva_bjet_pt"                :(lambda event, sample: event.bjet_pt),
                    "mva_bjet_Z1_deltaR"         :(lambda event, sample: event.bjet_Z1_deltaR), 
                    "mva_bjet_non_Z1l1_deltaR"   :(lambda event, sample: event.bjet_nonZ1l1_deltaR), 
                    }

elif args.variables == 'highrankjet0bjet' :
    mva_variables = {
                     "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
                     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),

                     "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),

                     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                     "mva_met_pt"                :(lambda event, sample: event.met_pt),
                     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                     "mva_nBTag"                 :(lambda event, sample: event.nBTag),

                     "mva_W_pt"                   :(lambda event, sample: event.W_pt),

                     "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),

                     "mva_bjet_pt"                :(lambda event, sample: event.bjet_pt),
                    }

elif args.variables == 'highrankjet01bjet' :
    mva_variables = {
                     "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
                     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),
                     "mva_jet1_pt"               :(lambda event, sample: event.JetGood_pt[1]          if event.nJetGood >=2 else 0),
                     "mva_jet1_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[1] if (event.nJetGood >=2 and event.JetGood_btagDeepB[0]>-10) else -10),

                     "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),

                     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                     "mva_met_pt"                :(lambda event, sample: event.met_pt),
                     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                     "mva_nBTag"                 :(lambda event, sample: event.nBTag),

                     "mva_W_pt"                   :(lambda event, sample: event.W_pt),

                     "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),

                     "mva_bjet_pt"                :(lambda event, sample: event.bjet_pt),
                    }

elif args.variables == 'highrankjet0btagbjet' :
    mva_variables = {
                     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),
                     "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),

                     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                     "mva_met_pt"                :(lambda event, sample: event.met_pt),
                     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                     "mva_nBTag"                 :(lambda event, sample: event.nBTag),

                     "mva_W_pt"                   :(lambda event, sample: event.W_pt),

                     "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),

                     "mva_bjet_pt"                :(lambda event, sample: event.bjet_pt),
                    }

elif args.variables == 'highrankbjetdRnonZ' :
    mva_variables = {
                     "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),

                     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                     "mva_met_pt"                :(lambda event, sample: event.met_pt),
                     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                     "mva_nBTag"                 :(lambda event, sample: event.nBTag),

                     "mva_W_pt"                   :(lambda event, sample: event.W_pt),

                     "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),

                     "mva_bjet_pt"                :(lambda event, sample: event.bjet_pt),

                     "mva_jet0_nonZl1_deltaR"     :(lambda event, sample: event.jet0_nonZ1_l1_deltaR    if event.nJetGood >=1 else -1),
                     "mva_jet1_nonZl1_deltaR"     :(lambda event, sample: event.jet1_nonZ1_l1_deltaR    if event.nJetGood >=2 else -1),
                    }

elif args.variables == 'plotvar' :
    mva_variables = {
                     "mva_bjet_Z1_deltaR"         :(lambda event, sample: event.bjet_Z1_deltaR),
                     "mva_W_pt"                   :(lambda event, sample: event.W_pt),
                     "lnonZ1_pt"                  :(lambda event, sample: event.lep_pt[event.nonZ1_l1_index]),
                     "lnonZ1_eta"                 :(lambda event, sample: event.lep_eta[event.nonZ1_l1_index]),
                     "mva_jet1_nonZl1_deltaR"     :(lambda event, sample: event.jet1_nonZ1_l1_deltaR    if event.nJetGood >=2 else -1),
                    }

elif  args.variables == 'original': 
    mva_variables = {
                    "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                    "mva_met_pt"                :(lambda event, sample: event.met_pt),
                    "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                    "mva_nBTag"                 :(lambda event, sample: event.nBTag),
    
                    "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
                    "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
                    "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),              
                    "mva_jet1_pt"               :(lambda event, sample: event.JetGood_pt[1]          if event.nJetGood >=2 else 0),
                    "mva_jet1_eta"              :(lambda event, sample: event.JetGood_eta[1]         if event.nJetGood >=2 else -10),
                    "mva_jet1_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[1] if (event.nJetGood >=2 and event.JetGood_btagDeepB[1]>-10) else -10),
    
                    "mva_nonZ1_l1_pt"           :(lambda event, sample: event.lep_pt[event.nonZ1_l1_index]),
                   
                    "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
                    "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
                    "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),
#                    "mva_Z1_mass"                :(lambda event, sample: event.Z1_mass),
    
                    "mva_nonZ1_l1_Z1_deltaR"     :(lambda event, sample: event.nonZ1_l1_Z1_deltaR),
     
                    "mva_jet0_Z1_deltaR"         :(lambda event, sample: event.jet0_Z1_deltaR         if event.nJetGood >=1 else -1),
                    "mva_jet0_nonZl1_deltaR"     :(lambda event, sample: event.jet0_nonZ1_l1_deltaR    if event.nJetGood >=1 else -1),
                    "mva_jet1_Z1_deltaR"         :(lambda event, sample: event.jet1_Z1_deltaR         if event.nJetGood >=2 else -1),
                    "mva_jet1_nonZl1_deltaR"     :(lambda event, sample: event.jet1_nonZ1_l1_deltaR    if event.nJetGood >=2 else -1),            
                    }

if args.mva == "bdt": 
    Bezeichnung = args.mva+"_NTrees"+str(args.NTrees)+"_maxDepth"+str(args.maxdepth)+"_nCuts"+str(args.ncuts)
    NTrees      = args.NTrees
    MaxDepth    = args.maxdepth
    nCuts       = args.ncuts 
    
    Bezeichnung = {
    "type"                : ROOT.TMVA.Types.kBDT,
    "name"                : str(Bezeichnung),
    "color"               : ROOT.kGreen,
    "options"             : ["!H","!V","NTrees="+str(NTrees),"BoostType=Grad","Shrinkage=0.20","UseBaggedBoost","GradBaggingFraction=0.5","SeparationType=GiniIndex","nCuts="+str(nCuts),"PruneMethod=NoPruning","MaxDepth="+str(MaxDepth)],
    }

elif args.mva == "mlp":
    sampling        = str(args.sampling)
    epoch           = str(args.epoch)
    layer           = args.layer
    sampling        = sampling.replace(".","c")
    epoch           = epoch.replace(".","c")
    layer           = layer.replace("p","+")
    layer           = layer.replace("m","-")
    layer           = layer.replace("c",",")

    Bezeichnung     = args.mva+args.layer+"_sampling"+sampling+"_epoch"+epoch
    
    Bezeichnung = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : str(Bezeichnung),
    "layers"              : layer,
    "color"               : ROOT.kMagenta+3,
    "options"             : ["!H","!V","VarTransform=Norm,Deco","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.02", "DecayRate=0.01","Sampling="+str(args.sampling),"SamplingEpoch="+str(args.epoch),"ConvergenceTests=1","CreateMVAPdfs=True","TestRate=10" ],
    }
