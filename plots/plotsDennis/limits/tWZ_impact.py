#!/usr/bin/env python

import ROOT, os
import array
import Analysis.Tools.syncer
from tWZ.Tools.user                      import plot_directory


def FindX(graph, y, xmin, xmax):
    x1=xmin
    x2=xmax
    stepsize = (xmax - xmin)/10000 
    foundx = False
    x = 0
    while not foundx:
        if x > xmax:
            break
        yval = graph.Eval(x)
        if yval > y:
            foundx = True
            x2 = x
        x += stepsize
    foundx = False
    x = 0
    while not foundx:
        if x < xmin:
            break
        yval = graph.Eval(x)
        if yval > y:
            foundx = True
            x1 = x
        x -= stepsize  
    return x1,x2  

WCnames = ["cHq1Re11", "cHq1Re22", "cHq1Re33", "cHq3Re11", "cHq3Re22", "cHq3Re33"]
channels = ["ZZ", "WZ", "ttZ", "combined"]


file = ROOT.TFile(plot_directory+"/Limits/Likelihoods.root")
file_noTWZ = ROOT.TFile(plot_directory+"/Limits_noTWZ/Likelihoods.root")

for WCname in WCnames:
    for channel in channels:
        graph = file.Get("Likelihood__"+channel+"__"+WCname)
        graph_noTWZ = file_noTWZ.Get("Likelihood__"+channel+"__"+WCname)
        c = ROOT.TCanvas("c", "c", 600, 600)
        ROOT.gStyle.SetLegendBorderSize(0)
        ROOT.gStyle.SetPadTickX(1)
        ROOT.gStyle.SetPadTickY(1)
        
        
        xlo = -10.0
        xhi = 10.0
        if "cHq3Re11" in WCname:
            xlo = -0.2
            xhi = 0.2
        xmin, xmax = FindX(graph, 10, xlo, xhi)
        graph.GetXaxis().SetLimits(xmin, xmax)    
        graph.SetLineColor(ROOT.kAzure+7)
        graph.SetLineWidth(2)
        graph.Draw('AL')
        graph.SetTitle('')
        graph.GetXaxis().SetTitle(WCname)
        graph.GetYaxis().SetTitle('-2 #Delta ln L')
        graph_noTWZ.Draw("L SAME")
        graph_noTWZ.SetLineColor(ROOT.kRed)
        graph_noTWZ.SetLineWidth(2)
        line68 =  ROOT.TLine(xmin, 1,    xmax, 1)
        line95 =  ROOT.TLine(xmin, 3.84, xmax, 3.84)
        line68.SetLineColor(13)
        line68.SetLineWidth(2)
        line68.SetLineStyle(2)
        line68.Draw("SAME")
        line95.SetLineColor(13)
        line95.SetLineWidth(2)
        line95.SetLineStyle(2)
        line95.Draw("SAME")
        graph.Draw("L SAME")
        graph_noTWZ.Draw("L SAME")
        leg = ROOT.TLegend(.4, .55, .6, .85)
        leg.AddEntry(graph, "tWZ EFT included", "l")
        leg.AddEntry(graph_noTWZ, "tWZ at SM", "l")
        leg.Draw()
        ROOT.gPad.RedrawAxis()
        outdir = plot_directory+"/Limits_tWZ_impact/"
        c.Print(outdir+"Limit__"+channel+"__"+WCname+".pdf")
        
        
