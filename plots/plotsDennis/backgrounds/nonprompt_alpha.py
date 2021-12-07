#!/usr/bin/env python

import ROOT, os
from math                                        import sqrt
import Analysis.Tools.syncer
from tWZ.Tools.user                              import plot_directory
from MyRootTools.backgroundAlpha.backgroundAlpha import backgroundAlpha
from MyRootTools.plotter.Plotter                 import Plotter


dir_CR = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1/Run2018/all/trilepVetoT-minDLmass12-onZ1-deepjet0-met60/"
dir_SR = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1_noData/Run2018/all/trilepT-minDLmass12-onZ1-deepjet0-met60/"

f_CR = ROOT.TFile(dir_CR+"Results.root")
f_SR = ROOT.TFile(dir_SR+"Results.root")

histname = "Z1_pt"
backgrounds = ["tWZ", "ttX", "tZq", "triBoson"]
signals = ["WZ", "ZZ", "ttZ"]

rebin = 2

nonprompt_CR = f_CR.Get(histname+"__nonprompt")
nonprompt_SR = f_SR.Get(histname+"__nonprompt")
data_CR = f_CR.Get(histname+"__data")
nonprompt_CR.Rebin(rebin)
nonprompt_SR.Rebin(rebin)
data_CR.Rebin(rebin)

estimator = backgroundAlpha()
estimator.setFitRange(0, 300)
estimator.setFitFormula("[0]*exp([1]*x)")
estimator.setFitFormula("[0]+[1]*x")
estimator.setMCCR(nonprompt_CR)
estimator.setMCSR(nonprompt_SR)
estimator.setDataCR(data_CR)
for bkg in backgrounds:
    hist = f_CR.Get(histname+"__"+bkg)
    hist.Rebin(rebin)
    estimator.subtractBackgroundCR(hist)
for sig in signals:
    hist = f_CR.Get(histname+"__"+sig+"__cHq1Re11=0.0000")
    hist.Rebin(rebin)
    estimator.subtractBackgroundCR(hist)
prediction = estimator.getPrediction()
fitfunctions = estimator.getAlphaFunctions()
alpha_hist = estimator.getAlphaHist()


c=ROOT.TCanvas()
alpha_hist.Draw("HIST")
colors = [ROOT.kRed, ROOT.kAzure, ROOT.kGreen]
for i, function in enumerate(fitfunctions):
    function.SetLineColor(colors[i])
    function.Draw("SAME")
c.Print(plot_directory+"/backgrounds/Alpha.pdf")

# Plot prediction
p2 = Plotter("Nonprompt_alpha")
p2.plot_dir = plot_directory+"/backgrounds/"
p2.addSignal(prediction, "Data driven", ROOT.kBlack)
p2.addSignal(nonprompt_SR, "MC prediction", ROOT.kRed)
p2.draw()
