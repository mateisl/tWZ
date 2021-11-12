from RootTools.core.standard import *
import Analysis.Tools.syncer
import os
import ROOT
import array
import numpy as np
from tWZ.Tools.user                      import plot_directory

def getAUC2( xvalues, yvalues ):
    # The TGraph.Integral function only returns the size of the area enclosed by the TGraph.
    # So, we add points at (x,0) for every existing value of (x,y).
    # This has to happen in reverse order, so a graph with values
    # x = [1, 2, 3],          y = [1, 2, 3] is converted to
    # x = [1, 2, 3, 3, 2, 1], y = [1, 2, 3, 0, 0, 0]
    # Then we can use the Integral function
    newx = []
    newy = []
    for x in xvalues:
        newx.append(x) # set same x values
    for y in yvalues:
        newy.append(y) # set same y values
    for x in reversed(xvalues):
        newx.append(x) # x values in revers
        newy.append(0) # y is all zeros
    # create TGraph and return integral
    newgraph = ROOT.TGraph(len(newx), array.array('d',newx), array.array('d',newy))
    return newgraph.Integral()

def getAUC( graph ):
    stepsize = 0.001
    integral = 0.0
    x_lo = 0
    x_hi = x_lo + stepsize
    while x_hi < 1.0:
        int_tmp = (x_hi-x_lo)*0.5*(graph.Eval(x_lo)+graph.Eval(x_hi))
        integral += int_tmp
        x_lo += stepsize
        x_hi += stepsize
    return integral


ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/tWZ/Tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()

dirname = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/tWZ_v3_noData/Run2016/all/trilepT-minDLmass12-onZ1-njet4p-deepjet1/"
roc = {}
auc = {}
filename = "MVA_score.root"
f = ROOT.TFile.Open(dirname+filename)
models = ["MVA_tWZ_3l", "MVA_tWZ_3l_topReco"]
legends = {
    "MVA_tWZ_3l":           "MVA 3l",
    "MVA_tWZ_3l_topReco":   "MVA 3l + top reco",
}

for model in models:
    node = "TWZnode"
    signal = "tWZ"
    # backgrounds = ["ttZ", "ttX", "tZq", "WZ", "ZZ", "triBoson", "nonprompt"]
    backgrounds = ["ttZ"]
    print "Getting histograms from model %s" %(model)
    # get signal and backgrounds
    print "  reading", model+"__"+signal+"__"+node, "..."
    sig = f.Get(model+"__"+signal+"__"+node)
    print "  reading", model+"__"+backgrounds[0]+"__"+node, "..."
    bkg = f.Get(model+"__"+backgrounds[0]+"__"+node)
    if len(backgrounds)>1:
        for i in range(len(backgrounds)):
            if i==0:
                continue
            else:
                print "  reading", model+"__"+backgrounds[i]+"__"+node, "..."
                tmp = f.Get(model+"__"+backgrounds[i]+"__"+node)
                bkg.Add(tmp)
    # normalize to unit area
    sig.Scale(1./sig.Integral())
    bkg.Scale(1./bkg.Integral())

    sig_eff = []
    bkg_eff = []
    for i_bin in reversed(range(1,sig.GetNbinsX()+1)):
        sig_eff.append( sig.Integral(i_bin, sig.GetNbinsX()))
        bkg_eff.append( bkg.Integral(i_bin, sig.GetNbinsX()))

    roc[model] = ROOT.TGraph(len(sig_eff), array.array('d',bkg_eff), array.array('d',sig_eff))
    auc[model] = getAUC(roc[model])

ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/tWZ/Tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()

colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen, ROOT.kAzure-7, 798]
styles = [1, 7, 2, 1, 1, 1, 1]

ROOT.gStyle.SetOptStat(0)
c1 = ROOT.TCanvas()
leg = ROOT.TLegend(0.33, 0.2, 0.9, 0.35)
leg.SetFillStyle(0)
leg.SetShadowColor(ROOT.kWhite)
leg.SetBorderSize(0)
firstkey = roc.keys()[0]
roc[firstkey].Draw("AL")
roc[firstkey].GetXaxis().SetTitle("ttZ background efficiency")
roc[firstkey].GetYaxis().SetTitle("tWZ signal efficiency")
roc[firstkey].GetXaxis().SetRangeUser(0,1)
roc[firstkey].GetYaxis().SetRangeUser(0,1)
icol = 0
for model in models:
    roc[model].SetLineColor(colors[icol])
    roc[model].SetLineWidth(2)
    roc[model].SetLineStyle(styles[icol])
    roc[model].SetMarkerStyle(0)
    roc[model].Draw("L SAME")
    leg.AddEntry( roc[model], legends[model]+", auc = %.2f" %(auc[model]), "l")
    icol += 1
leg.Draw()
diagonal = ROOT.TLine(0,0,1,1)
diagonal.SetLineWidth(2)
diagonal.SetLineColor(13)
diagonal.SetLineStyle(7)
diagonal.Draw("SAME")
c1.SetTitle("")
c1.RedrawAxis()
c1.Print(os.path.join(plot_directory, "roc.pdf"))

##### MVA scores
for model in models:
    h_tWZ = f.Get(model+"__tWZ__"+node)
    h_ttZ = f.Get(model+"__ttZ__"+node)
    h_tWZ.Scale(1./h_tWZ.Integral())
    h_ttZ.Scale(1./h_ttZ.Integral())
    c2 = ROOT.TCanvas()
    h_tWZ.GetXaxis().SetTitle("tWZ MVA score")
    h_tWZ.GetXaxis().SetNdivisions(505)
    h_tWZ.SetLineColor(ROOT.kRed)
    h_tWZ.SetFillStyle(0)
    h_tWZ.SetLineWidth(2)
    h_ttZ.SetLineColor(ROOT.kAzure+4)
    h_ttZ.SetFillStyle(0)
    h_ttZ.SetLineWidth(2)
    h_tWZ.Draw("HIST")
    h_ttZ.Draw("HIST SAME")
    leg = ROOT.TLegend(0.2, 0.65, 0.5, 0.85)
    leg.SetFillStyle(0)
    leg.SetShadowColor(ROOT.kWhite)
    leg.SetBorderSize(0)
    leg.AddEntry(h_tWZ, "tWZ", "l")
    leg.AddEntry(h_ttZ, "ttZ", "l")
    leg.Draw()
    c2.SetTitle("")
    c2.RedrawAxis()
    c2.Print(os.path.join(plot_directory, model+"_score.pdf"))

Analysis.Tools.syncer.sync()
