from RootTools.core.standard import *
import Analysis.Tools.syncer
import os
import ROOT
import array

#import pkl file 
filename = {}
filename['ttZ_M'] =  "/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/ttZ_noData/RunII/dilepM-onZ1-minDLmass12-njet5p-btag1p.pkl" 
filename['ttZ_L'] =  "/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/ttZ_noData/RunII/dilepL-onZ1-minDLmass12-njet5p-btag1p.pkl" 
filename['ttZ_VL'] =  "/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/ttZ_noData/RunII/dilepVL-onZ1-minDLmass12-njet5p-btag1p.pkl" 
filename['ttZ_T'] =  "/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/ttZ_noData/RunII/dilepT-onZ1-minDLmass12-njet5p-btag1p.pkl" 
dirname = "/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/ttZ_noData/RunII/"

files = filename.keys() 

ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/TMB/Tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()

roc = {}

import pickle
#hist_dict = pickle.load(file(filename))

histonames = {"flat":'ttz_tt_dy_TTZ', "lstm":'ttz_tt_dy_LSTM_TTZ'}
names = histonames.keys()

for f in files : 
    hist_dict = pickle.load(file(filename[f]))
    print f 
    for name in names: 
        sig = hist_dict[histonames[name]][0][0] # TTZ
        bkg = hist_dict[histonames[name]][0][1] # DY
        bkg.Add(hist_dict[histonames[name]][0][2]) #TTLep 
        
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
    
    l.AddEntry( roc["flat"+f],f+ "flat input" )
    l.AddEntry( roc["lstm"+f],f+ "including LSTM layer" )
    c1 = ROOT.TCanvas()
    roc['lstm'+f].Draw("AL")
    roc['lstm'+f].SetTitle("")
    roc['lstm'+f].GetXaxis().SetTitle("total background efficiency")
    roc['lstm'+f].GetYaxis().SetTitle("ttZ signal efficiency")
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
