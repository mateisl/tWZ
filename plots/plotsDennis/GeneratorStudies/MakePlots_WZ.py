#!/usr/bin/env python

import ROOT, os
import Analysis.Tools.syncer
from tWZ.Tools.user                      import plot_directory
from tWZ.Tools.helpers                   import cosThetaStarNew, cosThetaStarTop

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")

args = argParser.parse_args()

path = "/users/dennis.schwarz/CMSSW_10_6_20/src/"

histograms = {
    "PtWdecays_LHE": ROOT.TH1F("PtWdecays_LHE", "PtWdecays_LHE", 20, 0, 400),
    "AWdecays_LHE": ROOT.TH1F("AWdecays_LHE", "AWdecays_LHE", 10, 0, 1),
    "CosThetaStar_LHE": ROOT.TH1F("CosThetaStar_LHE", "CosThetaStar_LHE", 10, -1, 1),
    "WMass_LHE": ROOT.TH1F("WMass_LHE", "WMass_LHE", 20, 0, 100),
}

xtitles = {
    "PtWdecays_LHE": "LHE W decay products p_{T}",
    "AWdecays_LHE": "LHE W decay products |#frac{p_{T,1}-p_{T,2}}{p_{T,1}+p_{T,2}}|",
    "CosThetaStar_LHE": "cos #theta* W",
    "WMass_LHE": "m_{W}",
}

setups = [
    # {
    # "filename": path+"/WZ-decay_00/WZ-decay-vec_cHq3Re11_0p0_10000events.root",
    # "legendtext": "decay + single insertion, cHq3Re11 = 0.0 ",
    # "color": ROOT.kBlack,
    # },
    {
    "filename": path+"/WZ-decay_01/WZ-decay-vec_cHq3Re11_0p1_10000events.root",
    "legendtext": "decay + single insertion, cHq3Re11 = 0.1 ",
    "color": ROOT.kRed,
    },
    {
    "filename": path+"/WZ-decay_05/WZ-decay-vec_cHq3Re11_0p5_10000events.root",
    "legendtext": "decay + single insertion, cHq3Re11 = 0.5 ",
    "color": ROOT.kGreen,
    },
    {
    "filename": path+"/WZ-decay_10/WZ-decay-vec_cHq3Re11_1p0_10000events.root",
    "legendtext": "decay + single insertion, cHq3Re11 = 1.0 ",
    "color": ROOT.kAzure+7,
    },
]

# Setup the hists for each file
counter = 0
for setup in setups:
    counter += 1
    for key in histograms:
        setup[key] = histograms[key].Clone(key+str(counter))



for setup in setups:
    print "reading %s ..." %(setup["filename"])
    eventcounter = 0
    file = ROOT.TFile(setup["filename"])
    tree = file.Get("Events")
    if tree.GetEntry(0)<=0:
        print 'EMPTY TREE'
    for event in tree:
        if eventcounter%1000==0: print "  - %i events processed" %(eventcounter)
        # in WZ the W decay products have indices 2 (lepton) and 3 (neutrino)
        setup["PtWdecays_LHE"].Fill(event.LHEPart_pt[2])
        setup["PtWdecays_LHE"].Fill(event.LHEPart_pt[3])
        setup["AWdecays_LHE"].Fill(abs( (event.LHEPart_pt[2]-event.LHEPart_pt[3])/(event.LHEPart_pt[2]+event.LHEPart_pt[3])) )
        lepton = ROOT.TLorentzVector()
        if abs(event.LHEPart_pdgId[2]) == 11:
            lepmass = 0.0005
        elif abs(event.LHEPart_pdgId[2]) == 13:
            lepmass = 0.1
        elif abs(event.LHEPart_pdgId[2]) == 15:
            lepmass = 1.777
        lepton.SetPtEtaPhiM(event.LHEPart_pt[2], event.LHEPart_eta[2], event.LHEPart_phi[2], lepmass)
        neutrino = ROOT.TLorentzVector()
        neutrino.SetPtEtaPhiM(event.LHEPart_pt[3], event.LHEPart_eta[3], event.LHEPart_phi[3], 0.0)
        W = lepton+neutrino
        setup["CosThetaStar_LHE"].Fill(cosThetaStarNew(lepton,W))
        setup["WMass_LHE"].Fill(W.M())

        eventcounter +=1

# Draw the plots
for key in histograms:
    c = ROOT.TCanvas("c"+key, "c"+key, 600, 600)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)
    leg = ROOT.TLegend(.5,.6,.85,.85)
    if "AWdecays" in key:     leg = ROOT.TLegend(.2,.2,.55,.55)
    if "CosThetaStar" in key: leg = ROOT.TLegend(.2,.2,.55,.55)

    for i,setup in enumerate(setups):
        setup[key].SetLineColor(setup["color"])
        setup[key].SetMarkerColor(setup["color"])
        setup[key].SetMarkerStyle(20)
        setup[key].SetLineWidth(2)
        setup[key].GetXaxis().SetTitle(xtitles[key])
        setup[key].GetYaxis().SetTitle("Events")
        setup[key].GetYaxis().SetRangeUser(0, 1.4*setup[key].GetMaximum())
        if i==0:
            setup[key].Draw("E1")
        else:
            setup[key].Draw("E1 SAME")
        leg.AddEntry(setup[key], setup["legendtext"], "l")
    leg.Draw()
    ROOT.gPad.RedrawAxis()
    c.Print(os.path.join(plot_directory, "GeneratorTest_WZ_cHq3Re11_"+key+".pdf"))
    c.Delete()

Analysis.Tools.syncer.sync()
