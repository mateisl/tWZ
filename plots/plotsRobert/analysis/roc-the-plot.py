from RootTools.core.standard import *
import Analysis.Tools.syncer
import os
import ROOT
import array


ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/TMB/Tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()

dirname = "/mnt/hephy/cms/robert.schoefbeck/www/tWZ/plots/analysisPlots/ttZ-flavor_noData/RunII/all_log/trilepT-minDLmass12-onZ1-njet3p-btag1p"
roc = {}
for name, filename in [
        ("3l-lstm-ttZISR", os.path.join(dirname, "ttZ_3l_flavor_LSTM_TTZ_ISR.root")),
        ("3l-ttZISR", os.path.join(dirname, "ttZ_3l_flavor_TTZ_ISR.root")),
        ]:

    f = ROOT.TFile.Open(filename)
    canvas = f.Get(f.GetListOfKeys().At(0).GetName())
    sig = canvas.GetListOfPrimitives().At(8)
    bkg = canvas.GetListOfPrimitives().At(2)

    bkg.Scale(-1)
    sig.Add(bkg)
    bkg.Scale(-1)

    bkg_bkg = canvas.GetListOfPrimitives().At(3)
    bkg_bkg.Scale(-1)
    bkg.Add(bkg_bkg)
    bkg_bkg.Scale(-1)

    print "sig",sig.GetName()
    print "bkg",bkg.GetName()
    print "bkg_bkg",bkg_bkg.GetName()

    sig.Scale(1./sig.Integral())
    bkg.Scale(1./bkg.Integral())

    sig_eff = []
    bkg_eff = []
    for i_bin in reversed(range(1,sig.GetNbinsX()+1)):
        sig_eff .append( sig.Integral(i_bin, sig.GetNbinsX()))
        bkg_eff .append( bkg.Integral(i_bin, sig.GetNbinsX()))
        #print i_bin, sig_eff, bkg_eff

    roc[name] = ROOT.TGraph(len(sig_eff), array.array('d',bkg_eff), array.array('d',sig_eff))

ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/TMB/Tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()

#dirname = "/mnt/hephy/cms/robert.schoefbeck/www/tWZ/plots/analysisPlots/tttt_3l_v2_noData/RunII/all_log/trilep4tM-offZ1-minDLmass12-njet4p-btag2p/"
#for name, filename in [
#        ("3l-lstm-ttZISR", os.path.join(dirname, "tttt_ttw_ttz_nonprompt_LSTM_TTTT.root")),
#        ("3l-ttZISR", os.path.join(dirname, "tttt_ttw_ttz_nonprompt_TTTT.root")),
#        ]:
#
#    f = ROOT.TFile.Open(filename)
#    canvas = f.Get(f.GetListOfKeys().At(0).GetName())
#    sig = canvas.GetListOfPrimitives().At(3)
#    sub = canvas.GetListOfPrimitives().At(4)
#    sub.Scale(-1)
#    sig.Add(sub)
#
#    bkg = canvas.GetListOfPrimitives().At(1)
#    sig.Scale(-1)
#    bkg.Add(sig)
#    sig.Scale(-1)
#
#    print "sig",sig.GetName()
#    print "bkg",bkg.GetName()
#
#    sig.Scale(1./sig.Integral())
#    bkg.Scale(1./bkg.Integral())
#
#    sig_eff = []
#    bkg_eff = []
#    for i_bin in reversed(range(1,sig.GetNbinsX()+1)):
#        sig_eff .append( sig.Integral(i_bin, sig.GetNbinsX()))
#        bkg_eff .append( bkg.Integral(i_bin, sig.GetNbinsX()))
#        #print i_bin, sig_eff, bkg_eff
#
#    roc[name] = ROOT.TGraph(len(sig_eff), array.array('d',bkg_eff), array.array('d',sig_eff))


roc["3l-lstm-ttZISR"].SetLineColor(ROOT.kBlue)
roc["3l-ttZISR"].SetLineColor(ROOT.kBlue)
roc["3l-lstm-ttZISR"].SetLineWidth(2)
roc["3l-ttZISR"].SetLineWidth(2)
roc["3l-ttZISR"].SetLineStyle(7)

l = ROOT.TLegend(0.33, 0.2, 0.9, 0.35)
l.SetFillStyle(0)
l.SetShadowColor(ROOT.kWhite)
l.SetBorderSize(0)

l.AddEntry( roc["3l-ttZISR"], "3l (flat input)" )
l.AddEntry( roc["3l-lstm-ttZISR"], "3l (flat+LSTM)" )
c1 = ROOT.TCanvas()

roc['3l-lstm-ttZISR'].Draw("AL")
roc['3l-lstm-ttZISR'].SetTitle("")
roc['3l-lstm-ttZISR'].GetXaxis().SetTitle("t#bar{t} + nonISR Z efficiency")
roc['3l-lstm-ttZISR'].GetYaxis().SetTitle("t#bar{t} + ISR Z efficiency")
d=0.0
roc['3l-lstm-ttZISR'].GetXaxis().SetRangeUser(0,1-d)
roc['3l-lstm-ttZISR'].GetYaxis().SetRangeUser(d,1)
roc['3l-lstm-ttZISR'].SetMarkerStyle(0)
roc['3l-ttZISR'].SetMarkerStyle(0)

roc['3l-lstm-ttZISR'].SetMarkerStyle(0)
roc['3l-ttZISR'].SetMarkerStyle(0)

roc['3l-lstm-ttZISR'].Draw("AL")
roc['3l-ttZISR'].Draw("L")

l.Draw()

ROOT.gStyle.SetOptStat(0)
c1.SetTitle("")
c1.RedrawAxis()
c1.Print(os.path.join(dirname, "roc.png"))
c1.Print(os.path.join(dirname, "roc.pdf"))
c1.Print(os.path.join(dirname, "roc.root"))
Analysis.Tools.syncer.sync()
