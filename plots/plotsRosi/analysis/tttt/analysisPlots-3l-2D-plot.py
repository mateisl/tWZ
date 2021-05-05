from RootTools.core.standard import *
import Analysis.Tools.syncer
import os
import ROOT
filename = "/mnt/hephy/cms/robert.schoefbeck/www//tWZ/plots/analysisPlots/tttt_3l_v2_noData/RunII/trilep4tM-offZ1-minDLmass12-njet4p-btag2p.pkl"
ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/TMB/Tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()

import pickle
hist_dict = pickle.load(file(filename))
mva  = 'tttt_ttw_ttz_nonprompt_LSTM'
classifier = 'TTTT_vs_TTZ'
#classifier = 'TTTT_vs_TTW'
bkgs = ['TTTT', 'TTZ', 'TTW'] 


niceName = {
    'TTTT':"t#bar{t}t#bar{t}",
    'TTZ':"t#bar{t}Z",
    'TTW':"t#bar{t}W",
}

ROOT.setTDRStyle()

c1 = ROOT.TCanvas()
h = {}
for bkg in bkgs:
    h[bkg] =  hist_dict[mva+'_'+bkg+'_'+classifier][0][0]
    h[bkg].Scale(1./h[bkg].Integral() )

h["TTW"].Scale(2)

h['TTTT']        .SetLineColor( ROOT.kOrange+1 )
h['TTTT']        .SetFillColor( ROOT.kOrange+1 )
h['TTTT']        .SetMarkerColor( ROOT.kOrange+1 )
h['TTTT']        .SetMarkerStyle(0)
h['TTZ']    .SetLineColor( ROOT.kAzure+4 )
h['TTZ']    .SetFillColor( ROOT.kAzure+4 )
h['TTZ']    .SetMarkerColor( ROOT.kAzure+4 )
h['TTZ']    .SetMarkerStyle(0)
h['TTW']    .SetLineColor( ROOT.kGreen+4 )
h['TTW']    .SetFillColor( ROOT.kGreen+4 )
h['TTW']    .SetMarkerColor( ROOT.kGreen+4 )
h['TTW']    .SetMarkerStyle(0)

l = ROOT.TLegend(0.65, 0.65, 0.95, 0.89)
l.SetFillStyle(0)
l.SetShadowColor(ROOT.kWhite)
l.SetBorderSize(0)

l.AddEntry( h['TTZ'], "t#bar{t}Z" )
l.AddEntry( h['TTW'], "t#bar{t}W" )
l.AddEntry( h['TTTT'], "t#bar{t}t#bar{t}" )
c1 = ROOT.TCanvas()
h['TTW'].Draw("BOX")
h['TTW'].SetTitle("")
h['TTW'].GetXaxis().SetTitle("p(%s)"%niceName[classifier.split('_vs_')[0]])
h['TTW'].GetYaxis().SetTitle("p(%s)"%niceName[classifier.split('_vs_')[1]])

h['TTZ'].Draw("BOXsame")
h['TTTT'].Draw("BOXsame")

c1.SetLogz(0)
l.Draw()
c1.RedrawAxis()
ROOT.gStyle.SetOptStat(0)
c1.SetTitle("")
c1.Print(filename.replace('.pkl', '_'+classifier+'_mva.png'))
c1.Print(filename.replace('.pkl', '_'+classifier+'_mva.pdf'))
c1.Print(filename.replace('.pkl', '_'+classifier+'_mva.root'))
Analysis.Tools.syncer.sync()
