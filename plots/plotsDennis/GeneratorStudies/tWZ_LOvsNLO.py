import ROOT
import Analysis.Tools.syncer
from tWZ.Tools.user                              import plot_directory
from MyRootTools.plotter.Plotter                 import Plotter


def BuildZCandidate(leptons):
    minDiffMass = 100000
    Z = ROOT.TLorentzVector()
    idx1 = -1
    idx2 = -1
    for i, lep_i in enumerate(leptons):
        for j, lep_j in enumerate(leptons):
            if i==j:
                continue
            Zcand = lep_i+lep_j
            if abs(Zcand.M()-91.2)<minDiffMass:
                Z = Zcand
                minDiffMass = abs(Zcand.M()-91.2)
                idx1 = i
                idx2 = j
    return Z, idx1, idx2


nanoAODnames = [
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_0.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_1.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_2.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_3.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_4.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_5.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_6.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_7.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_8.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_9.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_10.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_11.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_12.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_13.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_14.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_15.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_16.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_17.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_18.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_19.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_20.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_21.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_22.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_23.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_24.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_25.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_26.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_27.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_28.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_29.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_30.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_31.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_32.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_33.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_34.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_35.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_36.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_37.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_38.root",
	"/eos/vbc/group/cms/robert.schoefbeck/nanoAODSim_fast_private/dennis/tWZToLL01j_lepWFilter/tWZtoLL01j_lepWFilter_schoef-Autumn18-mAODv1-10222-58e8664f142fdf477807c95cd1ce2e2a_USER/NANOAODSIMoutput_39.root",
]
################################################################################

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--small',    action='store_true', help='Run only on a small subset of the data?', default=False)
args = argParser.parse_args()

NeventsMax = 100000
if args.small:
    NeventsMax = 1000
    print "SMALL option, only processing %s events per sample"%(NeventsMax)
files = [
    ("LO + jet lepWFilter", ["tWZToLL01j_lepWFilter.root"], ROOT.kBlack),
    # ("LO + jet lepWFilter nanoAOD", nanoAODnames, ROOT.kGreen),
    # ("LO + jet", ["/users/dennis.schwarz/CMSSW_10_6_20/src/tWZ01j_newMatching/tWZtoLL01j-vec_200k.root"], ROOT.kAzure+7),
    ("LO + jet (m_{b}=0)", ["/users/dennis.schwarz/CMSSW_10_6_20/src/tWZnoB/tWZnoB_200k.root"], 13),
    # ("LO + jet (m_{b}=0) new gen", ["/users/dennis.schwarz/CMSSW_10_6_20/src/tWZnoB/tWZnob_newGen_200k.root"], ROOT.kAzure+7),
    ("NLO DR", ["tWZ_NLO_DR.root"], ROOT.kRed),
    # ("NLO DR tlep Whad", ["tWZ_tlep_Whad_NLO_DR.root"], ROOT.kRed),
    # ("NLO DR tlep Wlep", ["tWZ_tlep_Wlep_NLO_DR.root"], ROOT.kBlue),
    # ("NLO DR thad Wlep", ["tWZ_thad_Wlep_NLO_DR.root"], ROOT.kGreen),
    # ("NLO DS", ["tWZ_NLO_DS.root"], ROOT.kGreen-2),
]

p_njet = Plotter("Njet")
p_njet.plot_dir = plot_directory+"/LOvsNLO/"
p_njet.xtitle = "Number of jets"
p_njet.ytitle = "a.u."
p_njet.showLumi = False

p_njet30 = Plotter("Njet30")
p_njet30.plot_dir = plot_directory+"/LOvsNLO/"
p_njet30.xtitle = "Number of jets"
p_njet30.ytitle = "a.u."
p_njet30.showLumi = False

p_njet40 = Plotter("Njet40")
p_njet40.plot_dir = plot_directory+"/LOvsNLO/"
p_njet40.xtitle = "Number of jets"
p_njet40.ytitle = "a.u."
p_njet40.showLumi = False

p_njet50 = Plotter("Njet50")
p_njet50.plot_dir = plot_directory+"/LOvsNLO/"
p_njet50.xtitle = "Number of jets"
p_njet50.ytitle = "a.u."
p_njet50.showLumi = False

p_ptjet = Plotter("Ptjet")
p_ptjet.plot_dir = plot_directory+"/LOvsNLO/"
p_ptjet.xtitle = "Jet p_{T}"
p_ptjet.ytitle = "a.u."
p_ptjet.showLumi = False

p_etajet = Plotter("etajet")
p_etajet.plot_dir = plot_directory+"/LOvsNLO/"
p_etajet.xtitle = "Jet #eta"
p_etajet.ytitle = "a.u."
p_etajet.showLumi = False

p_dRlepjet = Plotter("dRlepjet")
p_dRlepjet.plot_dir = plot_directory+"/LOvsNLO/"
p_dRlepjet.xtitle = "min #Delta R(lep,jet)"
p_dRlepjet.ytitle = "a.u."
p_dRlepjet.showLumi = False

p_dRneujet = Plotter("dRneujet")
p_dRneujet.plot_dir = plot_directory+"/LOvsNLO/"
p_dRneujet.xtitle = "min #Delta R(neutrino,jet)"
p_dRneujet.ytitle = "a.u."
p_dRneujet.showLumi = False

p_nlep_matched = Plotter("NleptonsMatched")
p_nlep_matched.plot_dir = plot_directory+"/LOvsNLO/"
p_nlep_matched.xtitle = "Number of leptons (Matched to W/Z)"
p_nlep_matched.ytitle = "a.u."
p_nlep_matched.showLumi = False

p_nWlep = Plotter("NleptonsFromW")
p_nWlep.plot_dir = plot_directory+"/LOvsNLO/"
p_nWlep.xtitle = "Number of leptons from W"
p_nWlep.ytitle = "a.u."
p_nWlep.showLumi = False

p_nZlep = Plotter("NleptonsFromZ")
p_nZlep.plot_dir = plot_directory+"/LOvsNLO/"
p_nZlep.xtitle = "Number of leptons from Z"
p_nZlep.ytitle = "a.u."
p_nZlep.showLumi = False

p_charge = Plotter("ChargeLepNonZ")
p_charge.plot_dir = plot_directory+"/LOvsNLO/"
p_charge.xtitle = "Charge of non Z lepton"
p_charge.ytitle = "a.u."
p_charge.showLumi = False
p_charge.yfactor = 2.5

p_chargeSum = Plotter("ChargeZ")
p_chargeSum.plot_dir = plot_directory+"/LOvsNLO/"
p_chargeSum.xtitle = "Charge of reconstructed Z"
p_chargeSum.ytitle = "a.u."
p_chargeSum.showLumi = False

p_ptlep = Plotter("Ptleptons")
p_ptlep.plot_dir = plot_directory+"/LOvsNLO/"
p_ptlep.xtitle = "Lepton p_{T}"
p_ptlep.ytitle = "a.u."
p_ptlep.showLumi = False

p_Zmass = Plotter("Zmass")
p_Zmass.plot_dir = plot_directory+"/LOvsNLO/"
p_Zmass.xtitle = "m_{Z}"
p_Zmass.ytitle = "a.u."
p_Zmass.showLumi = False

p_ptZ = Plotter("PtZ")
p_ptZ.plot_dir = plot_directory+"/LOvsNLO/"
p_ptZ.xtitle = "Z p_{T}"
p_ptZ.ytitle = "a.u."
p_ptZ.showLumi = False

p_ptZ_filter = Plotter("PtZ_filter")
p_ptZ_filter.plot_dir = plot_directory+"/LOvsNLO/"
p_ptZ_filter.xtitle = "Z p_{T}"
p_ptZ_filter.ytitle = "a.u."
p_ptZ_filter.showLumi = False

p_nBjets = Plotter("NBjets")
p_nBjets.plot_dir = plot_directory+"/LOvsNLO/"
p_nBjets.xtitle = "Number of b jets"
p_nBjets.ytitle = "a.u."
p_nBjets.showLumi = False

p_genweight = Plotter("GenWeight")
p_genweight.plot_dir = plot_directory+"/LOvsNLO/"
p_genweight.xtitle = "Gen weight"
p_genweight.ytitle = "a.u."
p_genweight.showLumi = False

h_njets = {}
h_njets30 = {}
h_njets40 = {}
h_njets50 = {}
h_ptjets = {}
h_etajets = {}
h_dRlepjet = {}
h_dRneujet = {}
h_nleps_matched = {}
h_nWlep = {}
h_nZlep = {}
h_charge = {}
h_chargeSum = {}
h_ptleps = {}
h_ptZ = {}
h_ptZ_filter = {}
h_Zmass = {}
h_nBjets = {}
h_genweight = {}

leptonIds = [11,13]
neutrinoIds = [12,14,16]

for (name,filenames,color) in files:
    print "Now reading %s" %(name)
    h_njets[name] = ROOT.TH1F("njets_"+name, "number of jets", 21, -0.5, 20.5)
    h_njets30[name] = ROOT.TH1F("njets30_"+name, "number of jets", 21, -0.5, 20.5)
    h_njets40[name] = ROOT.TH1F("njets40_"+name, "number of jets", 21, -0.5, 20.5)
    h_njets50[name] = ROOT.TH1F("njets50_"+name, "number of jets", 21, -0.5, 20.5)
    h_ptjets[name] = ROOT.TH1F("ptjets_"+name, "jet p_{T}", 30, 0, 300)
    h_etajets[name] = ROOT.TH1F("etajets_"+name, "jet #eta", 30, -5., 5.)
    h_dRlepjet[name] = ROOT.TH1F("dRlepjet_"+name, "#Delta R(lepton,jet)", 30, 0., 2.)
    h_dRneujet[name] = ROOT.TH1F("dRneujet_"+name, "#Delta R(lepton,jet)", 30, 0., 2.)
    h_nleps_matched[name] = ROOT.TH1F("nleptons_matched_"+name, "number of leptons", 11, -0.5, 10.5)
    h_nWlep[name] = ROOT.TH1F("nleptons_fromW_"+name, "number of leptons from W", 11, -0.5, 10.5)
    h_nZlep[name] = ROOT.TH1F("nleptons_fromZ_"+name, "number of leptons from Z", 11, -0.5, 10.5)
    h_charge[name] = ROOT.TH1F("charge_"+name, "charge of non Z leptons", 2, -1.5, 1.5)
    h_chargeSum[name] = ROOT.TH1F("chargeSum_"+name, "charge of Z", 3, -1.5, 1.5)
    h_ptleps[name] = ROOT.TH1F("ptleps_"+name, "lepton p_{T}", 30, 0, 300)
    h_ptZ[name] = ROOT.TH1F("ptZ_"+name, "Z p_{t}", 30, 0, 500)
    h_ptZ_filter[name] = ROOT.TH1F("ptZ_filter_"+name, "Z p_{t}", 30, 0, 500)
    h_Zmass[name] = ROOT.TH1F("ZMass_"+name, "m_{Z}", 8, 70.0, 110)
    h_nBjets[name] = ROOT.TH1F("nBjets_"+name, "number of b jets", 11, -0.5, 10.5)
    h_genweight[name] = ROOT.TH1F("genweight_"+name, "gen weight", 50, -10, 10)
    eventcounter = 0
    for filename in filenames:
        file = ROOT.TFile(filename)
        tree = file.Get("Events")
        if tree.GetEntry(0)<=0:
            print 'EMPTY TREE'
        for event in tree:
            if eventcounter%1000==0: print "  - %i events processed" %(eventcounter)
            if eventcounter > NeventsMax: break
            eventcounter+=1
            # Lepton selection
            leptons = []
            neutrinos = []
            lep_ids = []
            neutrinos = []
            N_Wleps = 0
            N_Zleps = 0
            allMothers = []
            for i in range(event.nGenPart):
                allMothers.append(event.GenPart_genPartIdxMother[i])
            for i in range(event.nGenPart):
                # Find a lepton that is not the mother of any particle
                if abs(event.GenPart_pdgId[i]) in leptonIds:
                    if i not in allMothers:
                        lep  = ROOT.TLorentzVector()
                        lep.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], 0)
                        foundMother = False
                        idx = i
                        # Go up the decay chain until a mother is found that has a different ID
                        while not foundMother:
                            motherIdx = event.GenPart_genPartIdxMother[idx]
                            pdgID = event.GenPart_pdgId[idx]
                            # print idx, pdgID, motherIdx, event.GenPart_pdgId[motherIdx]
                            if motherIdx>0 and pdgID != event.GenPart_pdgId[motherIdx]:
                                foundMother = True 
                                if abs(event.GenPart_pdgId[motherIdx]) == 24:
                                    N_Wleps+=1
                                elif abs(event.GenPart_pdgId[motherIdx]) == 23:
                                    N_Zleps+=1
                            else:
                                idx = motherIdx
                                if idx < 0:
                                    break
                        leptons.append(lep)
                        lep_ids.append(i)
                #####
                if abs(event.GenPart_pdgId[i]) in neutrinoIds:
                    lep  = ROOT.TLorentzVector()
                    lep.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], 0)
                    neutrinos.append(lep)
            # if len(leptons) != 3: continue
            if N_Wleps+N_Zleps != 3: continue
            Z, idx1, idx2 = BuildZCandidate(leptons)
            if abs(Z.M()-91.2)>10: continue
            for i, idx in enumerate(lep_ids):
                if i != idx1 and i != idx2:
                    nonZlep_idx = idx
                    h_charge[name].Fill(event.GenPart_pdgId[nonZlep_idx]/abs(event.GenPart_pdgId[nonZlep_idx]), event.genWeight)
            charge1 = event.GenPart_pdgId[lep_ids[idx1]]/abs(event.GenPart_pdgId[lep_ids[idx1]])
            charge2 = event.GenPart_pdgId[lep_ids[idx2]]/abs(event.GenPart_pdgId[lep_ids[idx2]])
            h_chargeSum[name].Fill(charge1+charge2, event.genWeight)
            h_ptZ[name].Fill(Z.Pt(), event.genWeight)
            if N_Wleps >= 1:
                h_ptZ_filter[name].Fill(Z.Pt(), event.genWeight)
            h_Zmass[name].Fill(Z.M(), event.genWeight)
            h_nleps_matched[name].Fill(len(leptons), event.genWeight)
            h_nWlep[name].Fill(N_Wleps, event.genWeight)
            h_nZlep[name].Fill(N_Zleps, event.genWeight)
            h_genweight[name].Fill(event.genWeight)
            for lep in leptons:
                h_ptleps[name].Fill(lep.Pt(), event.genWeight)
            # Count jets
            h_njets[name].Fill(event.nGenJet, event.genWeight)
            Njets30 = 0
            Njets40 = 0
            Njets50 = 0
            Nbjets = 0
            mindRlep = 1000
            mindRneu = 1000
            for i in range(event.nGenJet):
                jet  = ROOT.TLorentzVector()
                jet.SetPtEtaPhiM(event.GenJet_pt[i], event.GenJet_eta[i], event.GenJet_phi[i], event.GenJet_mass[i])
                for neu in neutrinos:
                    if jet.DeltaR(neu) < 0.4:
                        jet = jet-neu
                h_etajets[name].Fill(jet.Eta(), event.genWeight)
                h_ptjets[name].Fill(jet.Pt(), event.genWeight)
                if jet.Pt() > 30 and abs(jet.Eta()) < 2.4:
                    Njets30 += 1
                    if event.GenJet_partonFlavour[i] == 5:
                        Nbjets +=1
                    for lep in leptons:
                        if jet.DeltaR(lep) < mindRlep:
                            mindRlep = jet.DeltaR(lep)
                    for neu in neutrinos:
                        if jet.DeltaR(neu) < mindRneu:
                            mindRneu = jet.DeltaR(neu)
                if jet.Pt() > 40 and abs(jet.Eta()) < 2.4:
                    Njets40 += 1
                if jet.Pt() > 50 and abs(jet.Eta()) < 2.4:
                    Njets50 += 1

            h_njets30[name].Fill(Njets30, event.genWeight)
            h_njets40[name].Fill(Njets40, event.genWeight)
            h_njets50[name].Fill(Njets50, event.genWeight)
            h_nBjets[name].Fill(Nbjets, event.genWeight)
            h_dRlepjet[name].Fill(mindRlep, event.genWeight)
            h_dRneujet[name].Fill(mindRneu, event.genWeight)

            
    h_njets[name].Scale(1/h_njets[name].Integral())
    h_njets30[name].Scale(1/h_njets30[name].Integral())
    h_njets40[name].Scale(1/h_njets40[name].Integral())
    h_njets50[name].Scale(1/h_njets50[name].Integral())
    h_etajets[name].Scale(1/h_etajets[name].Integral())
    h_ptjets[name].Scale(1/h_ptjets[name].Integral())
    h_nleps_matched[name].Scale(1/h_nleps_matched[name].Integral())
    h_nWlep[name].Scale(1/h_nWlep[name].Integral())
    h_nZlep[name].Scale(1/h_nZlep[name].Integral())
    h_charge[name].Scale(1/h_charge[name].Integral())
    h_chargeSum[name].Scale(1/h_chargeSum[name].Integral())
    h_ptleps[name].Scale(1/h_ptleps[name].Integral())
    h_ptZ[name].Scale(1/h_ptZ[name].Integral())
    h_ptZ_filter[name].Scale(1/h_ptZ_filter[name].Integral())
    h_Zmass[name].Scale(1/h_Zmass[name].Integral())
    h_nBjets[name].Scale(1/h_nBjets[name].Integral())
    h_dRlepjet[name].Scale(1/h_dRlepjet[name].Integral())
    h_dRneujet[name].Scale(1/h_dRneujet[name].Integral())
    h_genweight[name].Scale(1/h_genweight[name].Integral())

    p_njet.addSignal(h_njets[name], name, color)
    p_njet30.addSignal(h_njets30[name], name, color)
    p_njet40.addSignal(h_njets40[name], name, color)
    p_njet50.addSignal(h_njets50[name], name, color)
    p_etajet.addSignal(h_etajets[name], name, color)
    p_ptjet.addSignal(h_ptjets[name], name, color)
    p_nlep_matched.addSignal(h_nleps_matched[name], name, color)
    p_nWlep.addSignal(h_nWlep[name], name, color)
    p_nZlep.addSignal(h_nZlep[name], name, color)
    p_charge.addSignal(h_charge[name], name, color)
    p_chargeSum.addSignal(h_chargeSum[name], name, color)
    p_ptlep.addSignal(h_ptleps[name], name, color)
    p_ptZ.addSignal(h_ptZ[name], name, color)
    p_ptZ_filter.addSignal(h_ptZ_filter[name], name, color)
    p_Zmass.addSignal(h_Zmass[name], name, color)
    p_nBjets.addSignal(h_nBjets[name], name, color)
    p_dRlepjet.addSignal(h_dRlepjet[name], name, color)
    p_dRneujet.addSignal(h_dRneujet[name], name, color)
    p_genweight.addSignal(h_genweight[name], name, color)

    file.Close()

p_njet.draw()
p_njet30.draw()
p_njet40.draw()
p_njet50.draw()
p_ptjet.draw()
p_etajet.draw()
p_dRlepjet.draw()
p_dRneujet.draw()
p_nlep_matched.draw()
p_nWlep.draw()
p_nZlep.draw()
p_charge.draw()
p_chargeSum.draw()
p_ptlep.draw()
p_ptZ.draw()
p_ptZ_filter.draw()
p_Zmass.draw()
p_nBjets.draw()
p_genweight.draw()
