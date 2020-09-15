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

def getbJets( event, sample=None ):
#    #bJets filtern( nur 2016 )
    alljets   = getCollection( event, 'Jet', jetVarNames , 'nJet')  
    goodjets = filter( isAnalysisJet, alljets ) 
    bJets = filter(lambda j:isBJet(j, tagger=b_tagger, year=2016) and abs(j['eta'])<=2.4, goodjets)

    if len(bJets)>=1:
        event.bJet = bJets[0]
        event.bJet_Z1_deltaR      = deltaR({'eta':bJets[0]['eta'], 'phi':bJets[0]['phi']}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
        event.bJet_nonZ1l1_deltaR = deltaR({'eta':bJets[0]['eta'], 'phi':bJets[0]['phi']}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    else:
        event.bJet_Z1_deltaR      = -1 
        event.bJet_nonZ1l1_deltaR = -1 

sequence.append( getbJets )

def getM3l( event, sample=None):
    # get the invariant mass of the 3l system
    l = []
    for i in range(3):
        l.append(ROOT.TLorentzVector())
        l[i].SetPtEtaPhiM(event.lep_pt[i], event.lep_eta[i], event.lep_phi[i],0)
    event.m3l = (l[0] + l[1] + l[2]).M()
sequence.append( getM3l )

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

sequence.append( getAngles )

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

all_mva_variables = {

# global event properties     
     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
     "mva_met_pt"                :(lambda event, sample: event.met_pt),
     "mva_m3l"                   :(lambda event, sample: event.m3l),
     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
     "mva_nBTag"                 :(lambda event, sample: event.nBTag),

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
                }

#varcon1
#mva_variables_ = {
#    "mva_Z1_pt",
#    "mva_Z1_eta",
#    "mva_Z1_cosThetaStar",
#    "mva_ht",
#    "mva_met_pt",
#    "mva_nJetGood",
#    "mva_nBTag",
#    "mva_W_pt",
#    "mva_Z1_j1_deltaPhi",
#    "mva_nonZ1_l1_Z1_deltaPhi",
#}
#
#varcon2
#mva_variables_ = {
#    "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
#    "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
#    "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),
#    "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
#    "mva_met_pt"                :(lambda event, sample: event.met_pt),
#    "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
#    "mva_nBTag"                 :(lambda event, sample: event.nBTag),
#    "mva_W_pt"                   :(lambda event, sample: event.W_pt),
#    "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),
#    "mva_nonZ1_l1_Z1_deltaPhi"   :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nlep >= 2 else -1 ),
#    "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
#    "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
#    "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),
#}
 
#varcon3
#mva_variables_ = {
#    "mva_Z1_pt",
#    "mva_Z1_eta",
#    "mva_Z1_cosThetaStar",
#    "mva_ht",
#    "mva_met_pt",
#    "mva_nJetGood",
#    "mva_nBTag",
#    "mva_W_pt",
#    "mva_Z1_j1_deltaPhi",
#    "mva_nonZ1_l1_Z1_deltaPhi",
#    "mva_jet0_nonZl1_deltaR",
#    "mva_jet1_nonZl1_deltaR",
#}

#varcon4
#mva_variables_ = {
#    "mva_Z1_pt"                  :(lambda event, sample: event.Z1_pt),
#    "mva_Z1_eta"                 :(lambda event, sample: event.Z1_eta),
#    "mva_Z1_cosThetaStar"        :(lambda event, sample: event.Z1_cosThetaStar),
#    "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
#    "mva_met_pt"                :(lambda event, sample: event.met_pt),
#    "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
#    "mva_nBTag"                 :(lambda event, sample: event.nBTag),
#    "mva_W_pt"                   :(lambda event, sample: event.W_pt),
#    "mva_Z_j1_deltaPhi"          :(lambda event, sample: event.Z1_j1_deltaPhi           if event.nJetGood >=2 else -1),
#    "mva_nonZ1_l1_Z1_deltaPhi"   :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nlep >= 2 else -1 ),
#    "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
#    "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
#    "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),
#    "mva_jet0_nonZl1_deltaR",
#    "mva_jet1_nonZl1_deltaR",
#}

# Using all variables
mva_variables_ = all_mva_variables.keys()

mva_variables = {key:value for key, value in all_mva_variables.iteritems() if key in mva_variables_}

#bdt
all_bdt_nTr1000_maxD1_mNS5 = {
    "type"                : ROOT.TMVA.Types.kBDT,
    "name"                : 'all_bdt_Tr1000_maxD1_mNS5',
    "color"               : ROOT.kGreen,
    "options"             : ["!H","!V","NTrees=220", "BoostType=AdaBoost","SeparationType=GiniIndex","nCuts=20","PruneMethod=NoPruning","MaxDepth=1", "MinNodeSize=5"],
}

all_bdt_nTr1000_maxD1_mNS10 = {
    "type"                : ROOT.TMVA.Types.kBDT,
    "name"                : 'all_bdt_nTr1000_maxD1_mNS10',
    "color"               : ROOT.kGreen+4,
    "options"             : ["!H","!V","NTrees=1000","BoostType=AdaBoost","SeparationType=GiniIndex","nCuts=20","PruneMethod=NoPruning","MaxDepth=1", "MinNodeSize=10"],
}

all_bdt_nTr1000_maxD3_mNS10 = {
    "type"                : ROOT.TMVA.Types.kBDT,
    "name"                : 'all_bdt_nTr1000_maxD3_mNS10',
    "color"               : ROOT.kGreen-5,
    "options"             : ["!H","!V","NTrees=1000","BoostType=AdaBoost","SeparationType=GiniIndex","nCuts=20","PruneMethod=NoPruning","MaxDepth=3", "MinNodeSize=10"],
}

all_bdt_nTr1000_maxD4_mNS20 = {
    "type"                : ROOT.TMVA.Types.kBDT,
    "name"                : 'all_bdt_nTr1000_maxD4_mNS20',
    "color"               : ROOT.kGreen-10,
    "options"             : ["!H","!V","NTrees=1000","BoostType=AdaBoost","SeparationType=GiniIndex","nCuts=20","PruneMethod=NoPruning","MaxDepth=4", "MinNodeSize=20"],
}

#mlp
all_mlp_np5s0c5e0c8 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_np5s0c5e0c8',
    "layers"              : "N+5",
    "color"               : ROOT.kMagenta+3,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.5","SamplingEpoch=0.8","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }
all_mlp_np5s0c3e0c5 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_np5s0c3e0c5',
    "layers"              : "N+5",
    "color"               : ROOT.kMagenta-5,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.3","SamplingEpoch=0.5","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }
all_mlp_np5s0c5e1 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_np5s0c5e1',
    "layers"              : "N+5",
    "color"               : ROOT.kMagenta-10,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.5","SamplingEpoch=1","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }
all_mlp_np7s0c3e0c5 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_np7s0c3e0c5',
    "layers"              : "N+7",
    "color"               : ROOT.kBlue+3,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.3","SamplingEpoch=0.5","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }
all_mlp_np7s0c5e1 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_np7s0c5e1',
    "layers"              : "N+7",
    "color"               : ROOT.kBlue-4,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.5","SamplingEpoch=1","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }
all_mlp_np7s0c5e0c8 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_np7s0c5e0c8',
    "layers"              : "N+7",
    "color"               : ROOT.kBlue-5,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.5","SamplingEpoch=0.8","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }
all_mlp_ncnc1s0c3e0c5 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_ncnc1s0c3e0c5',
    "layers"              : "N,N,1",
    "color"               : ROOT.kBlue-10,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.3","SamplingEpoch=0.5","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }
all_mlp_ncnc1s0c5e1 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_ncnc1s0c5e1',
    "layers"              : "N,N,1",
    "color"               : ROOT.kRed,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.5","SamplingEpoch=1","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }
all_mlp_ncnc1s0c5e0c8 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_ncnc1s0c5e0c8',
    "layers"              : "N,N,1",
    "color"               : ROOT.kRed+3,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.5","SamplingEpoch=0.8","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }
all_mlp_ncnp5c1s0c3e0c5 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_ncnp5c1s0c3e0c5',
    "layers"              : "N,N+5,1",
    "color"               : ROOT.kRed-5,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.3","SamplingEpoch=0.5","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }
all_mlp_ncnp5c1s0c5e1 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_ncnp5c1s0c5e1',
    "layers"              : "N,N+5,1",
    "color"               : ROOT.kRed-9,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.5","SamplingEpoch=1","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }
all_mlp_ncnp5c1s0c5e0c8 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_ncnp5c1s0c5e0c8',
    "layers"              : "N,N+5,1",
    "color"               : ROOT.kOrange-6,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.5","SamplingEpoch=0.8","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }

all_mlp_np20 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_np20',
    "layers"              : "N+20",
    "color"               : ROOT.kYellow+3,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.3","SamplingEpoch=0.5","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }
all_mlp_np30 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_np30',
    "layers"              : "N+30",
    "color"               : ROOT.kGray,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.3","SamplingEpoch=0.5","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }
all_mlp_np40 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_np40',
    "layers"              : "N+40",
    "color"               : ROOT.kCyan+3,
    "options"             : ["!H","!V","VarTransform=Norm","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.04", "DecayRate=0.01","Sampling=0.3","SamplingEpoch=0.5","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=5" ],
    }

#oldconfig
all_mlp_oldconfig_ncnp5c1s0c3e0c5 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_oldconfig_ncnp5c1s0c3e0c5',
    "layers"              : "N,N+5,1",
    "color"               : ROOT.kCyan-6,
    "options"             : ["!H","!V","VarTransform=Norm,Deco","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.02", "DecayRate=0.01","Sampling=0.3","SamplingEpoch=0.5","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=10" ],
    }
all_mlp_oldconfig_ncnp5c1s0c3e0c3 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_oldconfig_ncnp5c1s0c3e0c3',
    "layers"              : "N,N+5,1",
    "color"               : ROOT.kBlack,
    "options"             : ["!H","!V","VarTransform=Norm,Deco","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.02", "DecayRate=0.01","Sampling=0.3","SamplingEpoch=0.3","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=10" ],
    }
all_mlp_oldconfig_np7s0c3e0c8 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_oldconfig_np7s0c3e0c8',
    "layers"              : "N+7",
    "color"               : ROOT.kYellow,
    "options"             : ["!H","!V","VarTransform=Norm,Deco","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.02", "DecayRate=0.01","Sampling=0.3","SamplingEpoch=0.8","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=10" ],
    }
all_mlp_oldconfig_np7c1s0c5e0c5 = {
    "type"                : ROOT.TMVA.Types.kMLP,
    "name"                : 'all_mlp_oldconfig_np7c1s0c5e0c5',
    "layers"              : "N+7",
    "color"               : ROOT.kGreen-8,
    "options"             : ["!H","!V","VarTransform=Norm,Deco","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.02", "DecayRate=0.01","Sampling=0.5","SamplingEpoch=0.5","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=10" ],
    }
