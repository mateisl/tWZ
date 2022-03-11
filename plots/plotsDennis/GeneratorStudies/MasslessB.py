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
################################################################################


files = [
    # ("LO + jet", "/users/dennis.schwarz/CMSSW_10_6_20/src/tWZ01j/tWZ01j-vec_10000events.root", ROOT.kBlue),
    ("LO + jet", "/users/dennis.schwarz/CMSSW_10_6_20/src/tWZ01j_newMatching/tWZ01j-vec_newMatching_50000events.root", ROOT.kAzure+7),
    # ("LO + jet (new matching, noEtaMax)", "/users/dennis.schwarz/CMSSW_10_6_20/src/tWZ01j_newMatching/tWZ01j-vec_newMatching_noExternalDecays_xqcut10_qcut30_noEtaMax.root", 798),
    # ("LO + jet (new matching nJetMax=0)", "/users/dennis.schwarz/CMSSW_10_6_20/src/tWZ01j_newMatching/tWZ01j-vec_newMatching_nJetMax0.root", ROOT.kAzure+2),
    ("NLO DR", "tWZ_NLO_DR.root", ROOT.kRed),
    ("NLO DS", "tWZ_NLO_DS.root", ROOT.kGreen-2),
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

p_nlep = Plotter("Nleptons")
p_nlep.plot_dir = plot_directory+"/LOvsNLO/"
p_nlep.xtitle = "Number of leptons"
p_nlep.ytitle = "a.u."
p_nlep.showLumi = False

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

p_nBjets = Plotter("NBjets")
p_nBjets.plot_dir = plot_directory+"/LOvsNLO/"
p_nBjets.xtitle = "Number of b jets"
p_nBjets.ytitle = "a.u."
p_nBjets.showLumi = False

h_njets = {}
h_njets30 = {}
h_njets40 = {}
h_njets50 = {}
h_ptjets = {}
h_etajets = {}
h_dRlepjet = {}
h_dRneujet = {}
h_nleps_matched = {}
h_nleps = {}
h_charge = {}
h_chargeSum = {}
h_ptleps = {}
h_Zmass = {}
h_nBjets = {}

leptonIds = [11,13]
neutrinoIds = [12,14,16]
allleptonIds = [11,13,15]

for (name,filename,color) in files:
    print "Now reading %s" %(name)
    h_njets[name] = ROOT.TH1F("njets_"+name, "number of jets", 21, -0.5, 20.5)
    h_njets30[name] = ROOT.TH1F("njets30_"+name, "number of jets", 21, -0.5, 20.5)
    h_njets40[name] = ROOT.TH1F("njets40_"+name, "number of jets", 21, -0.5, 20.5)
    h_njets50[name] = ROOT.TH1F("njets50_"+name, "number of jets", 21, -0.5, 20.5)
    h_ptjets[name] = ROOT.TH1F("ptjets_"+name, "jet p_{T}", 30, 0, 300)
    h_etajets[name] = ROOT.TH1F("etajets_"+name, "jet #eta", 30, -5., 5.)
    h_dRlepjet[name] = ROOT.TH1F("dRlepjet_"+name, "#Delta R(lepton,jet)", 30, 0., 2.)
    h_dRneujet[name] = ROOT.TH1F("dRneujet_"+name, "#Delta R(lepton,jet)", 30, 0., 2.)
    h_nleps[name] = ROOT.TH1F("nleptons_"+name, "number of leptons", 11, -0.5, 10.5)
    h_nleps_matched[name] = ROOT.TH1F("nleptons_matched_"+name, "number of leptons", 11, -0.5, 10.5)
    h_charge[name] = ROOT.TH1F("charge_"+name, "charge of non Z leptons", 2, -1.5, 1.5)
    h_chargeSum[name] = ROOT.TH1F("chargeSum_"+name, "charge of Z", 3, -1.5, 1.5)
    h_ptleps[name] = ROOT.TH1F("ptleps_"+name, "lepton p_{T}", 30, 0, 300)
    h_Zmass[name] = ROOT.TH1F("ZMass_"+name, "m_{Z}", 8, 70.0, 110)
    h_nBjets[name] = ROOT.TH1F("nBjets_"+name, "number of b jets", 11, -0.5, 10.5)
    eventcounter = 0
    file = ROOT.TFile(filename)
    tree = file.Get("Events")
    if tree.GetEntry(0)<=0:
        print 'EMPTY TREE'
    for event in tree:
        if eventcounter%1000==0: print "  - %i events processed" %(eventcounter)
        eventcounter+=1
        # Lepton selection
        leptons_all = []
        leptons = []
        neutrinos = []
        lep_ids = []
        allleptons = []
        neutrinos = []
        for i in range(event.nGenPart):
            if abs(event.GenPart_pdgId[i]) in allleptonIds:
                lep  = ROOT.TLorentzVector()
                lep.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], 0)
                allleptons.append(lep)
            if abs(event.GenPart_pdgId[i]) in neutrinoIds:
                lep  = ROOT.TLorentzVector()
                lep.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], 0)
                neutrinos.append(lep)
            if abs(event.GenPart_pdgId[i]) in leptonIds and event.GenPart_status[i]==1:
                lep  = ROOT.TLorentzVector()
                lep.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], 0)
                leptons_all.append(lep)
                motherIdx = event.GenPart_genPartIdxMother[i]
                if motherIdx >= 0:
                    if abs(event.GenPart_pdgId[motherIdx]) == 24 or abs(event.GenPart_pdgId[motherIdx]) == 23:
                        leptons.append(lep)
                        lep_ids.append(i)
        if len(leptons) != 3: continue
        Z, idx1, idx2 = BuildZCandidate(leptons)
        if abs(Z.M()-91.2)>10: continue
        for i, idx in enumerate(lep_ids):
            if i != idx1 and i != idx2:
                nonZlep_idx = idx
                h_charge[name].Fill(event.GenPart_pdgId[nonZlep_idx]/abs(event.GenPart_pdgId[nonZlep_idx]))
        charge1 = event.GenPart_pdgId[lep_ids[idx1]]/abs(event.GenPart_pdgId[lep_ids[idx1]])
        charge2 = event.GenPart_pdgId[lep_ids[idx2]]/abs(event.GenPart_pdgId[lep_ids[idx2]])
        h_chargeSum[name].Fill(charge1+charge2)
        h_Zmass[name].Fill(Z.M())
        h_nleps[name].Fill(len(leptons_all))
        h_nleps_matched[name].Fill(len(leptons))
        for lep in leptons:
            h_ptleps[name].Fill(lep.Pt())
        # Count jets
        h_njets[name].Fill(event.nGenJet)
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
            h_etajets[name].Fill(jet.Eta())
            h_ptjets[name].Fill(jet.Pt())
            if jet.Pt() > 30 and abs(jet.Eta()) < 2.4:
                Njets30 += 1
                if event.GenJet_partonFlavour[i] == 5:
                    Nbjets +=1
                for lep in allleptons:
                    if jet.DeltaR(lep) < mindRlep:
                        mindRlep = jet.DeltaR(lep)
                for neu in neutrinos:
                    if jet.DeltaR(neu) < mindRneu:
                        mindRneu = jet.DeltaR(neu)
            if jet.Pt() > 40 and abs(jet.Eta()) < 2.4:
                Njets40 += 1
            if jet.Pt() > 50 and abs(jet.Eta()) < 2.4:
                Njets50 += 1

        h_njets30[name].Fill(Njets30)
        h_njets40[name].Fill(Njets40)
        h_njets50[name].Fill(Njets50)
        h_nBjets[name].Fill(Nbjets)
        h_dRlepjet[name].Fill(mindRlep)
        h_dRneujet[name].Fill(mindRneu)

        ####
        if eventcounter > 50000: break
    h_njets[name].Scale(1/h_njets[name].Integral())
    h_njets30[name].Scale(1/h_njets30[name].Integral())
    h_njets40[name].Scale(1/h_njets40[name].Integral())
    h_njets50[name].Scale(1/h_njets50[name].Integral())
    h_etajets[name].Scale(1/h_etajets[name].Integral())
    h_ptjets[name].Scale(1/h_ptjets[name].Integral())
    h_nleps[name].Scale(1/h_nleps[name].Integral())
    h_nleps_matched[name].Scale(1/h_nleps_matched[name].Integral())
    h_charge[name].Scale(1/h_charge[name].Integral())
    h_chargeSum[name].Scale(1/h_chargeSum[name].Integral())
    h_ptleps[name].Scale(1/h_ptleps[name].Integral())
    h_Zmass[name].Scale(1/h_Zmass[name].Integral())
    h_nBjets[name].Scale(1/h_nBjets[name].Integral())
    h_dRlepjet[name].Scale(1/h_dRlepjet[name].Integral())
    h_dRneujet[name].Scale(1/h_dRneujet[name].Integral())

    p_njet.addSignal(h_njets[name], name, color)
    p_njet30.addSignal(h_njets30[name], name, color)
    p_njet40.addSignal(h_njets40[name], name, color)
    p_njet50.addSignal(h_njets50[name], name, color)
    p_etajet.addSignal(h_etajets[name], name, color)
    p_ptjet.addSignal(h_ptjets[name], name, color)
    p_nlep.addSignal(h_nleps[name], name, color)
    p_nlep_matched.addSignal(h_nleps_matched[name], name, color)
    p_charge.addSignal(h_charge[name], name, color)
    p_chargeSum.addSignal(h_chargeSum[name], name, color)
    p_ptlep.addSignal(h_ptleps[name], name, color)
    p_Zmass.addSignal(h_Zmass[name], name, color)
    p_nBjets.addSignal(h_nBjets[name], name, color)
    p_dRlepjet.addSignal(h_dRlepjet[name], name, color)
    p_dRneujet.addSignal(h_dRneujet[name], name, color)

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
p_nlep.draw()
p_charge.draw()
p_chargeSum.draw()
p_ptlep.draw()
p_Zmass.draw()
p_nBjets.draw()
