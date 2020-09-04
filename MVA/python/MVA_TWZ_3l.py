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
b_tagger = "DeepJet"
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
    if len(bJets)>=1:
        
        #bjet_pt  = bJets[0]['pt']
        bjet_eta = bJets[0]['eta']
        bjet_phi = bJets[0]['phi']
    
        event.bjet_Z1_deltaR      = deltaR({'eta':bjet_eta, 'phi':bjet_phi}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
        event.bjet_nonZ1l1_deltaR = deltaR({'eta':bjet_eta, 'phi':bjet_phi}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    else:
        event.bjet_Z1_deltaR      = -1 
        event.bjet_nonZ1l1_deltaR = -1 

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

mva_variables = {
                 "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
                 "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
                 "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),
                 
                 "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
                 "mva_met_pt"                :(lambda event, sample: event.met_pt),
                 "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                 "mva_nBTag"                 :(lambda event, sample: event.nBTag),
                 
                 "mva_W_pt"                   :(lambda event, sample: event.W_pt),

                 "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=1 else -1),
                 "mva_nonZ1_l1_Z1_deltaPhi"   :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nlep >= 1 else -1 ),
                 
#                     "mva_jet0_Z1_deltaR"         :(lambda event, sample: event.jet0_Z1_deltaR         if event.nJetGood >=1 else -1),
#                     "mva_jet0_nonZl1_deltaR"     :(lambda event, sample: event.jet0_nonZ1_l1_deltaR    if event.nJetGood >=1 else -1),
                }

bdt1 = {
    "type"                : ROOT.TMVA.Types.kBDT,
    "name"                : 'bdt1',
    "color"               : ROOT.kGreen,
    "options"             : ["!H","!V","NTrees=250","BoostType=Grad","Shrinkage=0.20","UseBaggedBoost","GradBaggingFraction=0.5","SeparationType=GiniIndex","nCuts=50","PruneMethod=NoPruning","MaxDepth=1"],
}
mlp1 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'mlp1',
    "layers"              : "N,N+7",
    "color"               : ROOT.kMagenta+3,
    "options"             : ["!H","!V","VarTransform=Norm,Deco","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.02", "DecayRate=0.01","Sampling=0.3","SamplingEpoch=0.8","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=10" ],
    }
