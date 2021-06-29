from RootTools.core.standard import *
import Analysis.Tools.syncer
import os
import ROOT
import array
import numpy as np
from math                                import sqrt
from tWZ.Tools.user                      import plot_directory



ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/tWZ/Tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()

dirname = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/tWZ_DeepJet_v1_noData/RunII/"
filename = "NjetNbjet.root"

signal = "tWZ"
plotnames = ["binNrL", "binNrM", "binNrT"]
backgrounds = ["ttZ", "ttX", "tZq", "WZ", "ZZ", "triBoson", "nonprompt"]

f = ROOT.TFile.Open(dirname+filename)

Njetcuts = [2,3,4]
Nbjetcuts = [1,2]
Nbins = len(Njetcuts)*len(Nbjetcuts)
hists = {}
for plot in plotnames:
    sigcontent = {}
    bkgcontent = {}
    # get signal and backgrounds
    sig = f.Get(plot+"__"+signal)
    bkg = f.Get(plot+"__"+backgrounds[0])
    if len(backgrounds)>1:
        for i in range(len(backgrounds)):
            if i==0: continue
            else:
                tmp = f.Get(plot+"__"+backgrounds[i])
                bkg.Add(tmp)
    for Njet in range(10):
        sigcontent[Njet] = {}
        bkgcontent[Njet] = {}
        for Nbjet in range(10):
            bin = 10*Njet+Nbjet+1
            sigcontent[Njet][Nbjet] = sig.GetBinContent(bin)
            bkgcontent[Njet][Nbjet] = bkg.GetBinContent(bin)

    binnumber = 1
    hist=ROOT.TH1F('','',Nbins,0,Nbins)
    for Njetcut in Njetcuts:
        for Nbjetcut in Nbjetcuts:
            sigsum = 0
            bkgsum = 0
            for Njet in range(10):
                for Nbjet in range(10):
                    if Njet<Njetcut:continue
                    if Nbjet != Nbjetcut: continue
                    sigsum += sigcontent[Njet][Nbjet]
                    bkgsum += bkgcontent[Njet][Nbjet]
            signifi = sigsum/sqrt(sigsum+bkgsum)
            binlabel = 'Njet#geq'+str(Njetcut)+',Nbjet='+str(Nbjetcut)
            hist.SetBinContent(binnumber,signifi)
            hist.GetXaxis().SetBinLabel(binnumber,binlabel)
            binnumber+=1

    hists[plot] = hist



ROOT.gStyle.SetOptStat(0)
colors=[ROOT.kGreen-7, ROOT.kAzure+7, ROOT.kRed-2]
c1 = ROOT.TCanvas()
ROOT.gPad.SetRightMargin(0.15)
hists["binNrL"].SetLineColor(ROOT.kGreen-7)
hists["binNrM"].SetLineColor(ROOT.kAzure+7)
hists["binNrT"].SetLineColor(ROOT.kRed-7)
hists["binNrL"].SetLineWidth(2)
hists["binNrM"].SetLineWidth(2)
hists["binNrT"].SetLineWidth(2)
hists["binNrL"].GetYaxis().SetTitle("S/#sqrt{S+B}")
hists["binNrL"].GetYaxis().SetRangeUser(0.0, 2.0)
hists["binNrL"].Draw("HIST")
hists["binNrM"].Draw("HIST SAME")
hists["binNrT"].Draw("HIST SAME")
leg=ROOT.TLegend(.5,.6,.85,.85)
leg.AddEntry(hists["binNrL"], "loose", "l")
leg.AddEntry(hists["binNrM"], "medium", "l")
leg.AddEntry(hists["binNrT"], "tight", "l")
leg.Draw("SAME")
c1.SetTitle("")
c1.RedrawAxis()
c1.Print(os.path.join(plot_directory, "NjetNbjet.pdf"))
Analysis.Tools.syncer.sync()
