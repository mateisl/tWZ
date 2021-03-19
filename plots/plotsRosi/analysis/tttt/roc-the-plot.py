from RootTools.core.standard import *
import Analysis.Tools.syncer
import os
import ROOT
import array

ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/TMB/Tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()

dirname = "/mnt/hephy/cms/robert.schoefbeck/www/tWZ/plots/analysisPlots/tttt_3l_v2_noData/RunII/all_log/trilep4tM-offZ1-minDLmass12-njet4p-btag2p/"
roc = {}
for name, filename in [
        ("lstm", os.path.join(dirname, "tttt_ttw_ttz_nonprompt_LSTM_TTTT.root")),
        ("flat", os.path.join(dirname, "tttt_ttw_ttz_nonprompt_TTTT.root")),
        ]:

    f = ROOT.TFile.Open(filename)
    canvas = f.Get(f.GetListOfKeys().At(0).GetName())
    sig = canvas.GetListOfPrimitives().At(3)
    sub = canvas.GetListOfPrimitives().At(4)
    sub.Scale(-1)
    sig.Add(sub)

    bkg = canvas.GetListOfPrimitives().At(1)
    sig.Scale(-1)
    bkg.Add(sig)
    sig.Scale(-1)

    print "sig",sig.GetName()
    print "bkg",bkg.GetName()

    sig.Scale(1./sig.Integral())
    bkg.Scale(1./bkg.Integral())

    sig_eff = []
    bkg_eff = []
    for i_bin in reversed(range(1,sig.GetNbinsX()+1)):
        sig_eff .append( sig.Integral(i_bin, sig.GetNbinsX()))
        bkg_eff .append( bkg.Integral(i_bin, sig.GetNbinsX()))
        #print i_bin, sig_eff, bkg_eff

    roc[name] = ROOT.TGraph(len(sig_eff), array.array('d',bkg_eff), array.array('d',sig_eff))

roc["lstm"].SetLineColor(ROOT.kRed)
roc["flat"].SetLineColor(ROOT.kBlue)
roc["lstm"].SetLineWidth(2)
roc["flat"].SetLineWidth(2)

l = ROOT.TLegend(0.4, 0.2, 0.9, 0.35)
l.SetFillStyle(0)
l.SetShadowColor(ROOT.kWhite)
l.SetBorderSize(0)

l.AddEntry( roc["flat"], "3l flat input" )
l.AddEntry( roc["lstm"], "3l including LSTM layer" )
c1 = ROOT.TCanvas()
roc['lstm'].Draw("AL")
roc['lstm'].SetTitle("")
roc['lstm'].GetXaxis().SetTitle("total background efficiency")
roc['lstm'].GetYaxis().SetTitle("t#bar{t}t#bar{t} signal efficiency")
d=0.4
roc['lstm'].GetXaxis().SetRangeUser(0,1-d)
roc['lstm'].GetYaxis().SetRangeUser(d,1)
roc['lstm'].SetMarkerStyle(0)
roc['flat'].SetMarkerStyle(0)

roc['flat'].Draw("L")

l.Draw()

ROOT.gStyle.SetOptStat(0)
c1.SetTitle("")
c1.RedrawAxis()
c1.Print(os.path.join(dirname, "roc.png"))
c1.Print(os.path.join(dirname, "roc.pdf"))
c1.Print(os.path.join(dirname, "roc.root"))
Analysis.Tools.syncer.sync()
