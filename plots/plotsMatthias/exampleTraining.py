#!/usr/bin/env python
''' BIT test script for learning nuisances 
'''

# Standard imports and batch mode
import ROOT, os
import sys
ROOT.gROOT.SetBatch(True)
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

import Analysis.Tools.syncer
import numpy as np

################################################################################
# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory', action='store', default='LeptonID')
argParser.add_argument('--selection',      action='store', default='trilep')
argParser.add_argument('--era',            action='store', type=str, default="Run2018")
args = argParser.parse_args()

################################################################################
# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

################################################################################
# Define the MC samples
from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed import *

sample_directory = "/scratch-cbe/users/dennis.schwarz/LeptonID/nanoTuples/2018/"
sample = Sample.fromDirectory(name="TTZ", treeName="Events", isData=False, color=ROOT.kAzure+4, texName="ttZ", directory=sample_directory+"TTZ")

maxEvents = 100

################################################################################
# Text on the plots
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)
tex.SetTextAlign(11) # align right


################################################################################
# Define sequences
sequence       = []

def buildZ( event, sample ):
    # This is how you can construct your own observables
    # In this example I build a Z boson from two leptons
    # I choose those two leptons that result in a mass close to mZ
    lepton1  = ROOT.TLorentzVector()
    lepton2  = ROOT.TLorentzVector()
    lepton3  = ROOT.TLorentzVector()
    lepton1.SetPtEtaPhiM(event.l1_pt, event.l1_eta, event.l1_phi, 0)
    lepton2.SetPtEtaPhiM(event.l2_pt, event.l2_eta, event.l2_phi, 0)
    lepton3.SetPtEtaPhiM(event.l3_pt, event.l3_eta, event.l3_phi, 0)

    mZ = 91.2
    Z1 = lepton1+lepton2
    Z2 = lepton1+lepton3
    Z3 = lepton2+lepton3
    Zmass = float('nan')
    mindiff = 1000
    for Z_candidate in [Z1,Z2,Z3]:
        if abs(Z_candidate.M()-mZ)<mindiff:
            mindiff = abs(Z_candidate.M()-mZ)
            Zmass = Z_candidate.M()
            
    event.Zmass = Zmass

sequence.append( buildZ )

################################################################################
# Read variables

read_variables = [
    "weight/F", "year/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F",
    "l2_pt/F", "l2_eta/F" , "l2_phi/F",
    "l3_pt/F", "l3_eta/F" , "l3_phi/F",
    "Muon[pt/F,eta/F,phi/F,pdgId/I,mediumId/O]",
    "Electron[pt/F,eta/F,phi/F,pdgId/I,mvaFall17V2Iso_WP80/O]",
    "Lep1SF/F", "Lep1SF_up/F", "Lep1SF_down/F",
    "Lep2SF/F", "Lep2SF_up/F", "Lep2SF_down/F",
    "Lep3SF/F", "Lep3SF_up/F", "Lep3SF_down/F",
]

r = sample.treeReader( variables = read_variables, sequence = sequence)
r.start()

counter = 0
while r.run():
    counter+=1
    print r.event.Zmass
    if counter>maxEvents: break 

# Boosting
sys.path.insert(0,os.path.expandvars("$CMSSW_BASE/src/BIT"))
from BoostedInformationTree import BoostedInformationTree

#training_features = np.array([[l1_pt, l2_pt, l3_pt, l1_eta, l2...], 
#                                [...], []])

bit = BoostedInformationTree(
        training_features     = training_features,
        training_weights      = training_weights,
        training_diff_weights = training_diff_weights,
        split_method          = 'vectorized_split_and_weight_sums',
        weights_update_method = 'vectorized',
        #bagging_fraction      = args.bagging_fraction,
        #**config.bit_cfg[derivative],
            ) 

bit.boost()

bit.save("v0.pkl")

