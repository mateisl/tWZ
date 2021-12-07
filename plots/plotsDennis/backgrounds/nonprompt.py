#!/usr/bin/env python

import ROOT, os
from math                                        import sqrt
import Analysis.Tools.syncer
from tWZ.Tools.user                              import plot_directory
from MyRootTools.backgroundABCD.backgroundABCD   import BackgroundABCD
from MyRootTools.plotter.Plotter                 import Plotter


# Directories
dir_A = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1/Run2018/all/trilepVetoT-minDLmass12-offZ1-deepjet0-met60/"
dir_B = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1/Run2018/all/trilepT-minDLmass12-offZ1-deepjet0-met60/"
dir_C = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1/Run2018/all/trilepVetoT-minDLmass12-onZ1-deepjet0-met60/"
dir_D = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1/Run2018/all/trilepT-minDLmass12-onZ1-deepjet0-met60/"

f_A = ROOT.TFile(dir_A+"Results.root")
f_B = ROOT.TFile(dir_B+"Results.root")
f_C = ROOT.TFile(dir_C+"Results.root")
f_D = ROOT.TFile(dir_C+"Results.root")

histname = "Z1_pt"
backgrounds = ["tWZ", "ttX", "tZq", "triBoson"]
signals = ["WZ", "ZZ", "ttZ"]

nonprompt_A = f_A.Get(histname+"__nonprompt")
nonprompt_B = f_B.Get(histname+"__nonprompt")
data_C = f_C.Get(histname+"__data")
nonprompt_D = f_D.Get(histname+"__nonprompt")

estimater = BackgroundABCD()
estimater.setSampleRegionA(nonprompt_A)
estimater.setSampleRegionB(nonprompt_B)
estimater.setSampleRegionC(data_C)
for bkg in backgrounds:
    hist = f_C.Get(histname+"__"+bkg)
    estimater.addBackgroundRegionC(hist)
for sig in signals:
    hist = f_C.Get(histname+"__"+sig+"__cHq1Re11=0.0000")
    estimater.addBackgroundRegionC(hist)
tf = estimater.getTransferFactor()
prediction, prediction_up, prediction_down = estimater.getBackgroundPrediction()

# Plot transfer factor
p = Plotter("Nonprompt_TF")
p.plot_dir = plot_directory+"/backgrounds/"
p.addSignal(tf, "Transfer factor", ROOT.kRed)
p.setCustomYRange(0, 2.)
p.draw()

# Plot prediction
p2 = Plotter("Nonprompt_prediction")
p2.plot_dir = plot_directory+"/backgrounds/"
p2.addBackground(prediction, "Data driven", 15)
p2.addSystematic(prediction_up, prediction_down, "TFuncert", "Data driven")
p2.addSignal(nonprompt_D, "MC prediction", ROOT.kRed)
p2.draw()
