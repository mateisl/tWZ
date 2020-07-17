#!/usr/bin/env python

# Standard imports
from operator                    import attrgetter
from math import pi, sqrt, cosh, cos
import ROOT

#from tWZ.Tools.objectSelection import jetSelector, getParticles 
from Analysis.Tools.helpers import deltaPhi, deltaR2, deltaR

# Logger
import Analysis.Tools.logger as logger
logger = logger.get_logger("INFO", logFile = None )

# Training variables
read_variables = [\
    "year/I",
    "nBTag/I",
    "nJetGood/I", "met_pt/F", 
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l2_eta/F", "l2_phi/F", 
    "JetGood[pt/F,eta/F,phi/F]", "bJets[pt/F,eta/F,phi/F]",
    "nJet/I", 
    "nlep/I",# "lep[%s]"%(",".join(lepVars)),
    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I", 
    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "jets[pt/F,eta/F,phi/F,btagDeepB/F,jetId/I]",
    #"ht/F", "mT/F", "LeptonTight0_pt/F","LeptonTight0_phi/F",
    #"LeptonTight0_eta/F", "m3/F", "Jet[pt/F,eta/F,phi/F,btagDeepB/F,jetId/I]",
    #"PhotonGood0_mvaID/F",
    #"PhotonGood0_pt/F", "PhotonGood0_eta/F","PhotonGood0_phi/F",
    #"nJet/I" 
    ]

# sequence 
sequence = []

#recoJetSel = {key:jetSelector(key) for key in [2016, 2017, 2018]}

def makeGoodJets(event, sample=None):
    # read all jets
    allJets     = getParticles( event, collVars=["pt","eta","phi","jetId","btagDeepB"], coll="Jet" )
    # selection
    JetGood = list( filter( lambda j: recoJetSel[event.year](j), allJets ) )
    JetGood = filter( lambda j :deltaR(j, {'eta':event.PhotonGood0_eta, 'phi':event.PhotonGood0_phi})>0.4, JetGood )
    JetGood = filter( lambda j :deltaR(j, {'eta':event.LeptonTight0_eta, 'phi':event.LeptonTight0_phi})>0.4, JetGood )
    event.JetGood  = JetGood
    event.nJetGood = len(JetGood)

sequence.append( makeGoodJets )


mva_variables = {
#                "mva_year"                  :(lambda event, sample: event.year),
                #"mva_PhotonGood0_mvaID"     :(lambda event, sample: event.PhotonGood0_mvaID),
                #"mva_ht"                    :(lambda event, sample: event.ht),
                "mva_met_pt"                :(lambda event, sample: event.met_pt),
                "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
                "mva_nBTag"                 :(lambda event, sample: event.nBTag),
                #"mva_mT"                    :(lambda event, sample: event.mT),
                #"mva_m3"                    :(lambda event, sample: event.m3),
                #"mva_dRlg"                  :(lambda event, sample: deltaR({'eta':event.PhotonGood0_eta, 'phi':event.PhotonGood0_phi},{'eta':event.LeptonTight0_eta,'phi':event.LeptonTight0_phi})),
                #"mva_mlg"                   :(lambda event, sample: sqrt(2*event.LeptonTight0_pt*event.PhotonGood0_pt*(cosh(event.LeptonTight0_eta-event.PhotonGood0_eta)-cos(event.LeptonTight0_phi-event.PhotonGood0_phi)))),

                "mva_jet0_pt"               :(lambda event, sample: event.JetGood[0]['pt']   if event.nJetGood >=1 else 0),
                "mva_jet0_eta"              :(lambda event, sample: event.JetGood[0]['eta']   if event.nJetGood >=1 else -10),
                "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood[0]['btagDeepB'] if (event.nJetGood >=1 and event.JetGood[0]['btagDeepB']>-10) else -10),              
                "mva_jet1_pt"               :(lambda event, sample: event.JetGood[1]['pt']   if event.nJetGood >=2 else 0),
                "mva_jet1_eta"              :(lambda event, sample: event.JetGood[1]['eta']   if event.nJetGood >=2 else -10),
                "mva_jet1_btagDeepB"        :(lambda event, sample: event.JetGood[1]['btagDeepB'] if (event.nJetGood >=2 and event.JetGood[1]['btagDeepB']>-10) else -10),
                "mva_jet2_pt"               :(lambda event, sample: event.JetGood[2]['pt']   if event.nJetGood >=3 else 0),
                "mva_jet2_eta"              :(lambda event, sample: event.JetGood[2]['eta']   if event.nJetGood >=3 else -10),
                "mva_jet2_btagDeepB"        :(lambda event, sample: event.JetGood[2]['btagDeepB']   if (event.nJetGood >=3 and event.JetGood[2]['btagDeepB']>-10) else -10),
                "mva_jet3_pt"               :(lambda event, sample: event.JetGood[3]['pt']   if event.nJetGood >=4 else 0),
                "mva_jet3_eta"              :(lambda event, sample: event.JetGood[3]['eta']   if event.nJetGood >=4 else -10),
                "mva_jet3_btagDeepB"        :(lambda event, sample: event.JetGood[3]['btagDeepB']   if (event.nJetGood >=4 and event.JetGood[3]['btagDeepB']>-10) else -10),
    
                }

bdt1 = {
"type"                : ROOT.TMVA.Types.kBDT,
"name"                : "bdt1",
"color"               : ROOT.kGreen,
"options"             : ["!H","!V","NTrees=250","BoostType=Grad","Shrinkage=0.20","UseBaggedBoost","GradBaggingFraction=0.5","SeparationType=GiniIndex","nCuts=50","PruneMethod=NoPruning","MaxDepth=1"],
}

mlp1 = {
"type"                : ROOT.TMVA.Types.kMLP,
"name"                : "mlp1",
"layers"              : "N+7",
"color"               : ROOT.kRed+5,
"options"             : ["!H","!V","VarTransform=Norm,Deco","NeuronType=sigmoid","NCycles=10000","TrainingMethod=BP","LearningRate=0.03", "DecayRate=0.01","Sampling=0.3","SamplingEpoch=0.8","ConvergenceTests=1","CreateMVAPdfs=True","TestRate=10" ],
}
