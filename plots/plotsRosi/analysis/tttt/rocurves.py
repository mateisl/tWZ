from RootTools.core.standard import *
import Analysis.Tools.syncer
import os
import ROOT
import array

#import pkl file 
filename = {}
filename['4tM'] =  "/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/tttt_3l_v3_noData/RunII/trilep4tM-offZ1-minDLmass12-njet4p-btag2p.pkl" 
filename['4tL'] =  "/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/tttt_3l_v3_noData/RunII/trilep4tL-offZ1-minDLmass12-njet4p-btag2p.pkl" 
filename['4tVL'] =  "/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/tttt_3l_v3_noData/RunII/trilep4tVL-offZ1-minDLmass12-njet4p-btag2p.pkl" 
filename['4tT'] =  "/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/tttt_3l_v3_noData/RunII/trilep4tT-offZ1-minDLmass12-njet4p-btag2p.pkl" 
dirname = "/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/tttt_3l_v3_noData/RunII/"

files = filename.keys() 

ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/TMB/Tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()

roc = {}

import pickle
#hist_dict = pickle.load(file(filename))

histonames = {"flat":'tttt_ttw_ttz_nonprompt_TTTT', "lstm":'tttt_ttw_ttz_nonprompt_LSTM_TTTT'}
names = histonames.keys()

for f in files : 
    hist_dict = pickle.load(file(filename[f]))
    print f 
    for name in names: 
        sig = hist_dict[histonames[name]][0][2] # TTTT
        bkg = hist_dict[histonames[name]][0][0] # TTW
        bkg.Add(hist_dict[histonames[name]][0][1]) #TTZ 
        bkg.Add(hist_dict[histonames[name]][0][3]) #nonprompt 
        
        sig.Scale(1./sig.Integral())
        bkg.Scale(1./bkg.Integral())
        
        sig_eff = []
        bkg_eff = []
        for i_bin in reversed(range(1,sig.GetNbinsX()+1)):
            sig_eff .append( sig.Integral(i_bin, sig.GetNbinsX()))
            bkg_eff .append( bkg.Integral(i_bin, sig.GetNbinsX()))
            #print i_bin, sig_eff, bkg_eff
        
        roc[name + f] = ROOT.TGraph(len(sig_eff), array.array('d',bkg_eff), array.array('d',sig_eff))

        #adjust plots  
        
    roc["lstm"+f].SetLineColor(ROOT.kRed)
    roc["flat"+f].SetLineColor(ROOT.kBlue)
    roc["lstm"+f].SetLineWidth(2)
    roc["flat"+f].SetLineWidth(2)
    
    l = ROOT.TLegend(0.4, 0.2, 0.9, 0.35)
    l.SetFillStyle(0)
    l.SetShadowColor(ROOT.kWhite)
    l.SetBorderSize(0)
    
    l.AddEntry( roc["flat"+f],f+ "3l flat input" )
    l.AddEntry( roc["lstm"+f],f+ "3l including LSTM layer" )
    c1 = ROOT.TCanvas()
    roc['lstm'+f].Draw("AL")
    roc['lstm'+f].SetTitle("")
    roc['lstm'+f].GetXaxis().SetTitle("total background efficiency")
    roc['lstm'+f].GetYaxis().SetTitle("t#bar{t}t#bar{t} signal efficiency")
    d=0.4
    roc['lstm'+f].GetXaxis().SetRangeUser(0,1-d)
    roc['lstm'+f].GetYaxis().SetRangeUser(d,1)
    roc['lstm'+f].SetMarkerStyle(0)
    roc['flat'+f].SetMarkerStyle(0)
    
    roc['flat'+f].Draw("L")

##adjust plots  
#
#roc["lstm"].SetLineColor(ROOT.kRed)
#roc["flat"].SetLineColor(ROOT.kBlue)
#roc["lstm"].SetLineWidth(2)
#roc["flat"].SetLineWidth(2)
#
#l = ROOT.TLegend(0.4, 0.2, 0.9, 0.35)
#l.SetFillStyle(0)
#l.SetShadowColor(ROOT.kWhite)
#l.SetBorderSize(0)
#
#l.AddEntry( roc["flat"], "3l flat input" )
#l.AddEntry( roc["lstm"], "3l including LSTM layer" )
#c1 = ROOT.TCanvas()
#roc['lstm'].Draw("AL")
#roc['lstm'].SetTitle("")
#roc['lstm'].GetXaxis().SetTitle("total background efficiency")
#roc['lstm'].GetYaxis().SetTitle("t#bar{t}t#bar{t} signal efficiency")
#d=0.4
#roc['lstm'].GetXaxis().SetRangeUser(0,1-d)
#roc['lstm'].GetYaxis().SetRangeUser(d,1)
#roc['lstm'].SetMarkerStyle(0)
#roc['flat'].SetMarkerStyle(0)
#
#roc['flat'].Draw("L")

#c1.SetLogz(0)

    l.Draw()
    ROOT.gStyle.SetOptStat(0)
    c1.SetTitle("")
    c1.SetTitle("")
    c1.RedrawAxis()
    c1.Print(os.path.join(dirname,f+ "roc.png"))
    c1.Print(os.path.join(dirname,f+ "roc.pdf"))
    c1.Print(os.path.join(dirname,f+ "roc.root"))
Analysis.Tools.syncer.sync()
