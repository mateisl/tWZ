#!/usr/bin/env python

import ROOT, os
import Analysis.Tools.syncer
from tWZ.Tools.user                      import plot_directory
from tWZ.Tools.helpers                   import cosThetaStarNew, cosThetaStarTop

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--mode',      action='store', default='default')
argParser.add_argument('--WC',      action='store', default='cHq3Re33')

args = argParser.parse_args()

leptons = [11,13,15]
neutrinos = [12,14,16]

path = "/users/dennis.schwarz/CMSSW_10_6_20/src/"

histograms = {
    "Njets": ROOT.TH1F("Njets", "Njets", 16, -0.5, 15.5),
    "Ptjets": ROOT.TH1F("Ptjets", "Ptjets", 20, 0, 400),
    "PtWdecays": ROOT.TH1F("PtWdecays", "PtWdecays", 20, 0, 400),
    "AWdecays": ROOT.TH1F("AWdecays", "AWdecays", 10, 0, 1),
    "NWlep": ROOT.TH1F("NWlep", "NWlep", 6, -0.5, 5.5),
    "PtWdecays_LHE": ROOT.TH1F("PtWdecays_LHE", "PtWdecays_LHE", 20, 0, 400),
    "AWdecays_LHE": ROOT.TH1F("AWdecays_LHE", "AWdecays_LHE", 10, 0, 1),
    "CosThetaStar_LHE": ROOT.TH1F("CosThetaStar_LHE", "CosThetaStar_LHE", 10, -1, 1),
    "CosThetaStarTop_LHE": ROOT.TH1F("CosThetaStarTop_LHE", "CosThetaStarTop_LHE", 10, -1, 1),
    "TopMass_LHE": ROOT.TH1F("TopMass_LHE", "TopMass_LHE", 20, 0, 250),
}

xtitles = {
    "Njets": "Number of jets",
    "Ptjets": "Jet p_{T}",
    "PtWdecays": "W decay products p_{T}",
    "AWdecays": "W decay products |#frac{p_{T,1}-p_{T,2}}{p_{T,1}+p_{T,2}}|",
    "NWlep": "Number of leptonic W decays",
    "PtWdecays_LHE": "LHE W decay products p_{T}",
    "AWdecays_LHE": "LHE W decay products |#frac{p_{T,1}-p_{T,2}}{p_{T,1}+p_{T,2}}|",
    "CosThetaStar_LHE": "cos #theta* W",
    "CosThetaStarTop_LHE": "cos #theta* top",
    "TopMass_LHE": "m_{top}",
}

if args.mode == "default":
    fname = ""
    setups = [
        {
        "filename": path+"/tt-decay_00/tt-decay-vec_cHq3Re33_0p0_10000events.root",
        "legendtext": "decay + single insertion, cHq3Re33 = 0.0 ",
        "color": ROOT.kBlack,
        },
        {
        "filename": path+"/tt-decay_01/tt-decay-vec_cHq3Re33_0p1_10000events.root",
        "legendtext": "decay + single insertion, cHq3Re33 = 0.1 ",
        "color": ROOT.kRed,
        },
        {
        "filename": path+"/tt-decay_10/tt-decay-vec_cHq3Re33_1p0_10000events.root",
        "legendtext": "decay + single insertion, cHq3Re33 = 1.0 ",
        "color": ROOT.kAzure+7,
        },
    ]
    if args.WC == "ctW":
        setups = [
            {
            "filename": path+"/tt-decay_00/tt-decay-vec_cHq3Re33_0p0_10000events.root",
            "legendtext": "decay + single insertion, ctW = 0.0 ",
            "color": ROOT.kBlack,
            },
            {
            "filename": path+"/tt-decay-ctW01/tt-decay-vec_ctW_0p1_10000events.root",
            "legendtext": "decay + single insertion, ctW = 0.1",
            "color": ROOT.kRed,
            },
            {
            "filename": path+"/tt-decay-ctW05/tt-decay-vec_ctW_0p5_10000events.root",
            "legendtext": "decay + single insertion, ctW = 0.5 ",
            "color": ROOT.kGreen,
            },
            {
            "filename": path+"/tt-decay-ctW10/tt-decay-vec_ctW_1p0_10000events.root",
            "legendtext": "decay + single insertion, ctW = 1.0 ",
            "color": ROOT.kAzure+7,
            },
        ]
    if args.WC == "ctW_1":
        setups = [
            {
            "filename": path+"/tt-decay_00/tt-decay-vec_cHq3Re33_0p0_10000events.root",
            "legendtext": "decay + single insertion, ctW = 0.0 ",
            "color": ROOT.kBlack,
            },
            {
            "filename": path+"/tt-decay-ctW01/tt-decay-vec_ctW_0p1_10000events_1.root",
            "legendtext": "decay + single insertion, ctW = 0.1",
            "color": ROOT.kRed,
            },
            {
            "filename": path+"/tt-decay-ctW05/tt-decay-vec_ctW_0p5_10000events_1.root",
            "legendtext": "decay + single insertion, ctW = 0.5 ",
            "color": ROOT.kGreen,
            },
            {
            "filename": path+"/tt-decay-ctW10/tt-decay-vec_ctW_1p0_10000events_1.root",
            "legendtext": "decay + single insertion, ctW = 1.0 ",
            "color": ROOT.kAzure+7,
            },
        ]
    elif args.WC == "cHtb":
        setups = [
            {
            "filename": path+"/tt-decay_00/tt-decay-vec_cHq3Re33_0p0_10000events.root",
            "legendtext": "decay + single insertion, cHtb = 0.0 ",
            "color": ROOT.kBlack,
            },
            {
            "filename": path+"/tt-decay-cHtb01/tt-decay-vec_cHtb_0p1_10000events.root",
            "legendtext": "decay + single insertion, cHtb = 0.1",
            "color": ROOT.kRed,
            },
            {
            "filename": path+"/tt-decay-cHtb05/tt-decay-vec_cHtb_0p5_10000events.root",
            "legendtext": "decay + single insertion, cHtb = 0.5 ",
            "color": ROOT.kGreen,
            },
            {
            "filename": path+"/tt-decay-cHtb10/tt-decay-vec_cHtb_1p0_10000events.root",
            "legendtext": "decay + single insertion, cHtb = 1.0 ",
            "color": ROOT.kAzure+7,
            },
        ]
    elif args.WC == "cHq3Re11":
        setups = [
            {
            "filename": path+"/tt-decay-cHq3Re1100/tt-decay-vec_cHq3Re11_0p0_10000events.root",
            "legendtext": "decay + single insertion, cHq3Re11 = 0.0 ",
            "color": ROOT.kBlack,
            },
            {
            "filename": path+"/tt-decay-cHq3Re1101/tt-decay-vec_cHq3Re11_0p1_10000events.root",
            "legendtext": "decay + single insertion, cHq3Re11 = 0.1",
            "color": ROOT.kRed,
            },
            {
            "filename": path+"/tt-decay-cHq3Re1105/tt-decay-vec_cHq3Re11_0p5_10000events.root",
            "legendtext": "decay + single insertion, cHq3Re11 = 0.5 ",
            "color": ROOT.kGreen,
            },
            {
            "filename": path+"/tt-decay-cHq3Re1110/tt-decay-vec_cHq3Re11_1p0_10000events.root",
            "legendtext": "decay + single insertion, cHq3Re11 = 1.0 ",
            "color": ROOT.kAzure+7,
            },
        ]
else:
    if args.mode == "0p0":
        dname = "00"
        fname = "0p0"
        wcstring = "cHq3Re33 = 0.0"
    elif args.mode == "0p1":
        dname = "01"
        fname = "0p1"
        wcstring = "cHq3Re33 = 0.1"
    elif args.mode == "1p0":
        dname = "10"
        fname = "1p0"
        wcstring = "cHq3Re33 = 1.0"
    setups = [
        {
        "filename": path+"/tt-decay_"+dname+"/tt-decay-vec_cHq3Re33_"+fname+"_10000events.root",
        "legendtext": "decay + single insertion, "+wcstring,
        "color": ROOT.kRed,
        },
        {
        "filename": path+"/tt-decay-double_"+dname+"/tt-decay-double-vec_cHq3Re33_"+fname+"_10000events.root",
        "legendtext": "decay + double insertion, "+wcstring,
        "color": ROOT.kAzure+7,
        },
        {
        "filename": path+"/tt_"+dname+"/tt-vec_cHq3Re33_"+fname+"_10000events.root",
        "legendtext": "only top width, "+wcstring,
        "color": ROOT.kBlack,
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
        NJets = len(event.GenJet_pt)
        NJets_pass = 0
        GenJetPt = []
        for i in range(NJets):
            GenJetPt.append(event.GenJet_pt[i])
            if event.GenJet_pt[i] > 30 and abs(event.GenJet_eta[i])<2.5:
                NJets_pass += 1
                setup["Ptjets"].Fill(event.GenJet_pt[i])
            setup["Njets"].Fill(NJets_pass)
        # in this ttbar sample the LHC particles at i=2 and 5 are the leptons
        # and 3 and 6 are the neutrinos from the W decay
        setup["PtWdecays_LHE"].Fill(event.LHEPart_pt[2])
        setup["PtWdecays_LHE"].Fill(event.LHEPart_pt[5])
        setup["PtWdecays_LHE"].Fill(event.LHEPart_pt[3])
        setup["PtWdecays_LHE"].Fill(event.LHEPart_pt[6])
        setup["AWdecays_LHE"].Fill(abs( (event.LHEPart_pt[2]-event.LHEPart_pt[3])/(event.LHEPart_pt[2]+event.LHEPart_pt[3])) )
        setup["AWdecays_LHE"].Fill(abs( (event.LHEPart_pt[5]-event.LHEPart_pt[6])/(event.LHEPart_pt[5]+event.LHEPart_pt[6])) )
        leptons = [ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector()]
        lepton_idx = [2,3,5,6]
        for i, idx in enumerate(lepton_idx):
            if abs(event.LHEPart_pdgId[idx]) == 11:
                lepmass = 0.0005
            elif abs(event.LHEPart_pdgId[idx]) == 13:
                lepmass = 0.1
            elif abs(event.LHEPart_pdgId[idx]) == 15:
                lepmass = 1.777
            else:
                lepmass = 0
            leptons[i].SetPtEtaPhiM(event.LHEPart_pt[idx], event.LHEPart_eta[idx], event.LHEPart_phi[idx], lepmass)
        W1 = leptons[0]+leptons[1]
        W2 = leptons[2]+leptons[3]
        setup["CosThetaStar_LHE"].Fill(cosThetaStarNew(leptons[0],W1))
        setup["CosThetaStar_LHE"].Fill(cosThetaStarNew(leptons[2],W2))
        b1 = ROOT.TLorentzVector()
        b2 = ROOT.TLorentzVector()
        b1.SetPtEtaPhiM(event.LHEPart_pt[4], event.LHEPart_eta[4], event.LHEPart_phi[4], 4.18)
        b2.SetPtEtaPhiM(event.LHEPart_pt[7], event.LHEPart_eta[7], event.LHEPart_phi[7], 4.18)
        top1 = b1+W1
        top2 = b2+W2
        setup["CosThetaStarTop_LHE"].Fill(cosThetaStarTop(leptons[0],W1,top1))
        setup["CosThetaStarTop_LHE"].Fill(cosThetaStarTop(leptons[2],W2,top2))
        setup["TopMass_LHE"].Fill(top1.M())
        setup["TopMass_LHE"].Fill(top2.M())
        eventcounter +=1

# Draw the plots
for key in histograms:
    c = ROOT.TCanvas("c"+key, "c"+key, 600, 600)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetLeftMargin(0.17)
    ROOT.gPad.SetBottomMargin(0.2)

    leg = ROOT.TLegend(.5,.6,.85,.85)
    if "AWdecays" in key:     leg = ROOT.TLegend(.2,.2,.55,.55)
    if "TopMass" in key:      leg = ROOT.TLegend(.2,.2,.55,.55)
    if "CosThetaStar" in key: leg = ROOT.TLegend(.2,.2,.55,.55)

    for i,setup in enumerate(setups):
        setup[key].SetLineColor(setup["color"])
        setup[key].SetMarkerColor(setup["color"])
        setup[key].SetMarkerStyle(20)
        setup[key].SetLineWidth(2)
        setup[key].GetXaxis().SetTitle(xtitles[key])
        setup[key].GetXaxis().SetTitleOffset(1.3)
        setup[key].GetYaxis().SetTitleOffset(1.3)
        setup[key].GetXaxis().SetNdivisions(505)
        setup[key].GetYaxis().SetNdivisions(505)
        setup[key].GetYaxis().SetTitle("Events")
        setup[key].GetYaxis().SetRangeUser(0, 1.4*setup[key].GetMaximum())
        if i==0:
            setup[key].Draw("E1")
        else:
            setup[key].Draw("E1 SAME")
        leg.AddEntry(setup[key], setup["legendtext"], "l")
    leg.Draw()
    ROOT.gPad.RedrawAxis()
    c.Print(os.path.join(plot_directory, "GeneratorTest_"+args.WC+"_"+fname+key+".pdf"))
    c.Delete()

Analysis.Tools.syncer.sync()
