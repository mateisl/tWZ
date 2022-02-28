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
from tWZ.Tools.helpers                   import getCollection
from tWZ.Tools.leptonSF_BIT              import leptonSF_BIT

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
argParser.add_argument('--noData',         action='store_true', default=False, help='also plot data?')
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--dataMCScaling',  action='store_true', help='Data MC scaling?', )
argParser.add_argument('--plot_directory', action='store', default='tWZ_v1')
argParser.add_argument('--era',            action='store', type=str, default="Run2018")
argParser.add_argument('--selection',      action='store', default='triMuon-vetoElec')
args = argParser.parse_args()

################################################################################
# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"
if args.noData:                       args.plot_directory += "_noData"

logger.info( "Working in era %s", args.era)
if args.dataMCScaling:
    logger.info( "Data/MC scaling active")
else:
    logger.info( "Data/MC scaling not active")

################################################################################
# Define the MC samples
from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed import *

if args.era == "Run2016":
    mc = [Summer16.TWZ_NLO_DR, Summer16.TTZ, Summer16.TTX_rare, Summer16.TZQ, Summer16.WZ, Summer16.triBoson, Summer16.ZZ, Summer16.nonprompt_3l]
elif args.era == "Run2017":
    mc = [Fall17.TWZ_NLO_DR, Fall17.TTZ, Fall17.TTX_rare, Fall17.TZQ, Fall17.WZ, Fall17.triBoson, Fall17.ZZ, Fall17.nonprompt_3l]
elif args.era == "Run2018":
    # mc = [Autumn18.TWZ_NLO_DR, Autumn18.TTZ, Autumn18.TTX_rare, Autumn18.TZQ, Autumn18.WZ, Autumn18.triBoson, Autumn18.ZZ, Autumn18.nonprompt_3l]
    mc = [Autumn18.TTZ]
elif args.era == "RunII":
    mc = [TWZ_NLO_DR, TTZ, TTX_rare, TZQ, WZ, triBoson, ZZ, nonprompt_3l]


mu2018SF = leptonSF_BIT( year=2018, ID="medium" )

################################################################################
# Define sequences
sequence       = []

def getLeptonSF( event, sample ):
    
    #
    idx_l1 = event.l1_index
    event.Lep1SF = mu2018SF.getSF( pdgId=event.lep_pdgId[idx_l1], pt=event.lep_pt[idx_l1], eta=((event.lep_eta[idx_l1]+event.lep_deltaEtaSC[idx_l1]) if abs(event.lep_pdgId[idx_l1])==11 else event.lep_eta[idx_l1]) )
    event.Lep1SF_up = mu2018SF.getSF( pdgId=event.lep_pdgId[idx_l1], pt=event.lep_pt[idx_l1], eta=((event.lep_eta[idx_l1]+event.lep_deltaEtaSC[idx_l1]) if abs(event.lep_pdgId[idx_l1])==11 else event.lep_eta[idx_l1]), sigma = +1 )
    event.Lep1SF_down = mu2018SF.getSF( pdgId=event.lep_pdgId[idx_l1], pt=event.lep_pt[idx_l1], eta=((event.lep_eta[idx_l1]+event.lep_deltaEtaSC[idx_l1]) if abs(event.lep_pdgId[idx_l1])==11 else event.lep_eta[idx_l1]), sigma = -1 )

    idx_l2 = event.l2_index
    event.Lep2SF = mu2018SF.getSF( pdgId=event.lep_pdgId[idx_l2], pt=event.lep_pt[idx_l2], eta=((event.lep_eta[idx_l2]+event.lep_deltaEtaSC[idx_l2]) if abs(event.lep_pdgId[idx_l2])==11 else event.lep_eta[idx_l2]) )
    event.Lep2SF_up = mu2018SF.getSF( pdgId=event.lep_pdgId[idx_l2], pt=event.lep_pt[idx_l2], eta=((event.lep_eta[idx_l2]+event.lep_deltaEtaSC[idx_l2]) if abs(event.lep_pdgId[idx_l2])==11 else event.lep_eta[idx_l2]), sigma = +1 )
    event.Lep2SF_down = mu2018SF.getSF( pdgId=event.lep_pdgId[idx_l2], pt=event.lep_pt[idx_l2], eta=((event.lep_eta[idx_l2]+event.lep_deltaEtaSC[idx_l2]) if abs(event.lep_pdgId[idx_l2])==11 else event.lep_eta[idx_l2]), sigma = -1 )

    idx_l3 = event.l3_index
    event.Lep3SF = mu2018SF.getSF( pdgId=event.lep_pdgId[idx_l3], pt=event.lep_pt[idx_l3], eta=((event.lep_eta[idx_l3]+event.lep_deltaEtaSC[idx_l3]) if abs(event.lep_pdgId[idx_l3])==11 else event.lep_eta[idx_l3]) )
    event.Lep3SF_up = mu2018SF.getSF( pdgId=event.lep_pdgId[idx_l3], pt=event.lep_pt[idx_l3], eta=((event.lep_eta[idx_l3]+event.lep_deltaEtaSC[idx_l3]) if abs(event.lep_pdgId[idx_l3])==11 else event.lep_eta[idx_l3]), sigma = +1 )
    event.Lep3SF_down = mu2018SF.getSF( pdgId=event.lep_pdgId[idx_l3], pt=event.lep_pt[idx_l3], eta=((event.lep_eta[idx_l3]+event.lep_deltaEtaSC[idx_l3]) if abs(event.lep_pdgId[idx_l3])==11 else event.lep_eta[idx_l3]), sigma = -1 )


sequence.append( getLeptonSF )

################################################################################
# Read variables

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",  "nJet/I",
    "nElectron/I", "nMuon/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I",
    "l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
    "l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
    "JetGood[pt/F,eta/F,phi/F,area/F,btagDeepB/F,btagDeepFlavB/F,index/I]",
    "Jet[pt/F,eta/F,phi/F,mass/F]",
    "lep[pt/F,eta/F,phi/F,pdgId/I,muIndex/I,eleIndex/I,deltaEtaSC/F,mediumId/I]",
    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I",
    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I,mediumId/O]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I,mvaFall17V2Iso_WP80/O]",
]

read_variables_MC = [
    'reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F', 'reweightTrigger/F',
    "genZ1_pt/F", "genZ1_eta/F", "genZ1_phi/F",
]


################################################################################
# define 3l selections
mu_string  = lepString('mu','VL') + "&&lep_mediumId"
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

store_variables = [
    "weight/F", "year/I",
    "nElectron/I", "nMuon/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F",
    "l2_pt/F", "l2_eta/F" , "l2_phi/F",
    "l3_pt/F", "l3_eta/F" , "l3_phi/F",
    "Lep1SF/F", "Lep1SF_up/F", "Lep1SF_down/F",
    "Lep2SF/F", "Lep2SF_up/F", "Lep2SF_down/F",
    "Lep3SF/F", "Lep3SF_up/F", "Lep3SF_down/F",
    ]
    
store_MuonVars = ["pt/F","eta/F","phi/F","pdgId/I","mediumId/O"]
store_MuonVars_names = [x.split('/')[0] for x in store_MuonVars]
store_variables_muon = ['Muon[%s]'% ( ','.join(store_MuonVars))]    

store_ElectronVars = ["pt/F","eta/F","phi/F","pdgId/I","mvaFall17V2Iso_WP80/O"]
store_ElectronVars_names = [x.split('/')[0] for x in store_ElectronVars]
store_variables_elec = ['Electron[%s]'% ( ','.join(store_ElectronVars))]  

for sample in mc:
    print "Running on", sample.name
    r = sample.treeReader( variables = read_variables, sequence = sequence, selectionString = cutInterpreter.cutString(args.selection) )

    def filler( event ):
        for x in store_variables:
            y = TreeVariable.fromString(x) if type(x)==type("") else x
            setattr(event,"%s"%(y.name),getattr( r.event,"%s"%(y.name))) 
        for imu in range(r.event.nMuon):
            for var in store_MuonVars_names:
                getattr(event, "Muon_"+var)[imu] = getattr(r.event, "Muon_"+var)[imu]
        for iel in range(r.event.nElectron):
            for var in store_ElectronVars_names:
                getattr(event, "Electron_"+var)[iel] = getattr(r.event, "Electron_"+var)[iel]

    
    maker = TreeMaker(
        sequence  = [ filler ],
        variables = [  TreeVariable.fromString(x) if type(x)==type("") else x for x in store_variables+store_variables_muon+store_variables_elec ],
        treeName = "Events"
        )
    maker.start()
    sample.chain.SetBranchStatus("*",1)

    count=0
    count_threshold = 1000
    r.start()
    while r.run():
        count+=1
        maker.run()
        if count >= count_threshold:
            count_threshold += 1000
            print "Processed", count, "events"
    

    sample_directory = "/scratch-cbe/users/dennis.schwarz/LeptonID/nanoTuples/2018/"
    outfile = ROOT.TFile(sample_directory+sample.name+"/Ntuple.root", "recreate")
    outfile.cd()
    maker.tree.Write()
    outfile.Close()

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
