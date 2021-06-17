from RootTools.core.standard import *
import Analysis.Tools.syncer
import os
import ROOT
import array
import numpy as np
from tWZ.Tools.user                      import plot_directory

# ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/tWZ/Tools/scripts/tdrstyle.C")
# ROOT.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
dirname = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/tWZ_topreco_v1_shape/Run2016/"
filename = "MatchHists.root"
processes = ["ttZ", "tWZ"]
plots = ["m_top_matched", "m_W_matched"]
axislabels = {
    "m_top_matched" : "m_{t} [GeV]",
    "m_W_matched"   : "m_{W} [GeV]",
}

file = ROOT.TFile.Open(dirname+filename)

for process in processes:
    for plot in plots:
        hist = file.Get(plot+"__"+process)
        range = [100,250] if 'm_top' in plot else [40,150]
        fit = ROOT.TF1("fit","gaus",range[0], range[1])
        hist.Fit('fit','QR')
        mean = fit.GetParameter(1)
        sigma = fit.GetParameter(2)

        c1 = ROOT.TCanvas("c", "c", 600, 600)
        ROOT.gPad.SetTickx(1)
        ROOT.gPad.SetTicky(1)
        ROOT.gPad.SetBottomMargin(.12)
        ROOT.gPad.SetLeftMargin(.14)

        hist.GetYaxis().SetRangeUser(0, 1.2*hist.GetMaximum())
        hist.SetTitle("")
        hist.GetYaxis().SetTitle("a.u")
        hist.GetXaxis().SetTitle(axislabels[plot])
        hist.GetXaxis().SetNdivisions(505)
        hist.SetLineColor(ROOT.kAzure-7)
        hist.SetFillColor(ROOT.kAzure-7)
        fit.SetLineWidth(2)
        fit.SetLineColor(ROOT.kRed)

        hist.Draw("HIST")
        fit.Draw("SAME")

        leg = ROOT.TLegend(0.55, 0.65, 0.9, 0.85)
        leg.SetFillStyle(0)
        leg.SetShadowColor(ROOT.kWhite)
        leg.SetBorderSize(0)
        leg.AddEntry( hist, process, "f")
        leg.AddEntry( fit, "Fit", "l")
        leg.AddEntry(0, "#mu = %.2f" %(mean), "");
        leg.AddEntry(0, "#sigma = %.2f" %(sigma), "");
        leg.Draw()
        c1.SetTitle("")
        c1.RedrawAxis()
        c1.Print(os.path.join(plot_directory, "MatchingHists_"+plot+"_"+process+".pdf"))
Analysis.Tools.syncer.sync()
