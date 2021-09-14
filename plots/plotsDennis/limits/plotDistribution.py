#!/usr/bin/env python

import ROOT,os
from math                                import sqrt
import Analysis.Tools.syncer
from tWZ.Tools.user                      import plot_directory

def getSysUncert(histup, histdown, histcentral):
    h_sys = histcentral.Clone()
    h_sys.Reset()
    Nbins = histcentral.GetSize()-2
    for i in range(Nbins):
        bin=i+1
        up      = histup.GetBinContent(bin)
        down    = histdown.GetBinContent(bin)
        central = histcentral.GetBinContent(bin)
        err = ( abs(central-up)+abs(central-down) )/2
        # print bin, abs(central-up), abs(central-down)
        h_sys.SetBinContent(bin, err)
    return h_sys

def getTotalUncert(central, syshists):
    h_total = central.Clone()
    Nbins = central.GetSize()-2
    for i in range(Nbins):
        bin=i+1
        staterr2 = h_total.GetBinError(bin)*h_total.GetBinError(bin)
        sum2 = 0
        for h in syshists:
            sum2 += h.GetBinContent(bin)*h.GetBinContent(bin)
        h_total.SetBinError(bin, sqrt(sum2+staterr2))
    return h_total

def getRatio(h1, h2):
    ratio = h1.Clone()
    Nbins = ratio.GetSize()-2
    for i in range(Nbins):
        bin=i+1
        if h2.GetBinContent(bin)==0:
            r=0
            e=0
        else:
            r = h1.GetBinContent(bin)/h2.GetBinContent(bin)
            e = h1.GetBinError(bin)/h2.GetBinContent(bin)
        ratio.SetBinContent(bin,r)
        ratio.SetBinError(bin,e)
    return ratio

colors = {
    'nonprompt':   ROOT.kMagenta-2,
    'singleTop':   40,
    'ttX' :        ROOT.kRed-10,
    'TTH' :        ROOT.kRed-10,
    'tWZ' :        ROOT.kRed,
    'ttZ' :        ROOT.kAzure+4,
    'TTW' :        ROOT.kGreen+2 ,
    'tZq' :        ROOT.kRed - 9,
    'TTTT' :       ROOT.kOrange+1,
    'WJetsToLNu' : ROOT.kRed-10,
    'WJets' :      ROOT.kRed-10,
    'diBoson' :    ROOT.kOrange,
    'multiBoson' : ROOT.kOrange,
    'ZZ' :         ROOT.kGreen+3,
    'WZ' :         ROOT.kAzure+6,
    'WW' :         ROOT.kOrange-7,
    'triBoson' :   ROOT.kOrange+1,
}

signal_colors = [ROOT.kAzure+7, ROOT.kRed-2, ROOT.kGreen+2, 798, 13, 15]

def plotDistribution(filename, region, histname, xtitle, backgrounds, WCname, SMpoint, EFTpoints, value_to_number, systs):
    file = ROOT.TFile(filename)
    leg = ROOT.TLegend(.55, .4, .9, .85)
    leg.SetBorderSize(0)

    sum = 0
    bkg_list = []
    h_stack = ROOT.THStack()
    for bkg in backgrounds:
        hist = file.Get(region+"__"+histname+"/"+bkg)
        hist.SetFillColor(colors[bkg])
        hist.SetLineColor(colors[bkg])
        integral = hist.Integral()
        sum+=integral
        bkg_list.append( (hist, integral, bkg) )
    h_totalbkg = hist.Clone()
    h_totalbkg.Reset()
    # Sort by integral
    def takeSecond(elem):
        return elem[1]
    bkg_list.sort(key=takeSecond)
    # Now fill stack
    for hist, integral, bkgname in bkg_list:
        h_stack.Add(hist)
        h_totalbkg.Add(hist)
    # Legend has to be filled in reverse order
    for hist, integral, bkgname in reversed(bkg_list):
        leg.AddEntry(hist, bkgname, "f")

    # Get SM signal
    h_SMsignal = file.Get(region+"__"+histname+"/EFT_"+WCname+"_"+str(value_to_number[WCname][SMpoint]))
    h_SMsignal.Add(h_totalbkg)
    h_SMsignal.SetLineWidth(2)
    h_SMsignal.SetLineColor(ROOT.kBlack)
    h_SMsignal.SetFillColor(0)
    leg.AddEntry(h_SMsignal, 'SM point', "l")

    # Copy h_SMsignal for nominal prediction
    h_totalMC = h_SMsignal.Clone()

    # Get Systematics
    h_sys_uncerts = [] # these will contian only the error ( (|central-up|+|central-down|) / 2)
    for sys in systs:
        # first get signal
        histup = file.Get(region+"__"+histname+"/EFT_"+WCname+"_"+str(value_to_number[WCname][SMpoint])+"__"+sys+"Up")
        histdown = file.Get(region+"__"+histname+"/EFT_"+WCname+"_"+str(value_to_number[WCname][SMpoint])+"__"+sys+"Down")
        for bkg in backgrounds:
            hup = file.Get(region+"__"+histname+"/"+bkg+"__"+sys+"Up")
            hdown = file.Get(region+"__"+histname+"/"+bkg+"__"+sys+"Down")
            histup.Add(hup)
            histdown.Add(hdown)
        h_sys_uncerts.append( getSysUncert(histup, histdown, h_totalMC) )
    h_sys_total = getTotalUncert(h_totalMC, h_sys_uncerts) # this contains total uncertainty (pure uncertainty value)
    h_sys_total.SetFillStyle(3245)
    h_sys_total.SetFillColor(13)
    h_sys_total.SetLineWidth(0)

    leg.AddEntry(h_sys_total, "Total uncertainty","f")



    # Get EFT Signals
    h_EFTsignals = []
    for i in range(len(EFTpoints)):
        hist = file.Get(region+"__"+histname+"/EFT_"+WCname+"_"+str(value_to_number[WCname][EFTpoints[i]]))
        hist.Add(h_totalbkg)
        hist.SetLineWidth(2)
        hist.SetLineColor(signal_colors[i])
        hist.SetFillColor(0)
        h_EFTsignals.append(hist)
        leg.AddEntry(hist, WCname+" = "+str(EFTpoints[i]),"l")

    # Draw plot
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c", "c", 600, 600)
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.31, 1, 1.0)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.19)
    pad1.Draw()
    pad1.cd()
    ymin = 0
    ymax = h_SMsignal.GetMaximum()*1.4
    xmin = h_SMsignal.GetXaxis().GetBinLowEdge(1)
    xmax = h_SMsignal.GetXaxis().GetBinUpEdge(h_SMsignal.GetSize()-2)
    binwidth = h_SMsignal.GetXaxis().GetBinWidth(1)

    h_SMsignal.GetXaxis().SetTitle('')
    h_SMsignal.GetYaxis().SetRangeUser(ymin, ymax)
    h_SMsignal.GetYaxis().SetNdivisions(505)

    h_SMsignal.GetYaxis().SetTitle('Events/'+str(int(binwidth))+' GeV')
    h_SMsignal.SetTitle('')
    h_SMsignal.Draw("HIST")

    h_sys_total.Draw("E2 SAME")
    for h in h_EFTsignals:
        h.Draw("HIST SAME")
    h_stack.Draw("HIST SAME")



    leg.Draw()

    # some style things
    titlefont = 43

    # This is for correct axis titles
    h_SMsignal.GetXaxis().SetLabelSize(0.)
    h_SMsignal.GetYaxis().SetLabelSize(0.)
    axis = ROOT.TGaxis( xmin, ymin, xmin, ymax, ymin, ymax, 505,"")
    axis.SetLabelOffset(0.01)
    axis.SetLabelFont(titlefont)
    axis.SetLabelSize(21)
    axis.SetNdivisions(505)

    axis.Draw()
    ROOT.gPad.RedrawAxis()


    # Second pad with ratio plots
    c.cd();
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetLeftMargin(0.19)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.38);
    pad2.Draw()
    pad2.cd()

    EFTratios = []
    isfirst = True
    for i in range(len(h_EFTsignals)):
        ratio = getRatio(h_EFTsignals[i], h_SMsignal)
        ratio.SetLineWidth(2)
        ratio.SetLineColor(signal_colors[i])
        ratio.SetTitle('')
        ratio.GetYaxis().SetRangeUser(0.5, 1.5)
        ratio.GetYaxis().SetNdivisions(505)
        ratio.GetYaxis().SetTitle('#frac{Data}{SM}')
        ratio.GetXaxis().SetTitle(xtitle)
        ratio.GetXaxis().SetTickLength(0.07)
        ratio.GetXaxis().SetTitleSize(25)
        ratio.GetXaxis().SetTitleFont(titlefont)
        ratio.GetXaxis().SetTitleOffset(4.0)
        ratio.GetXaxis().SetLabelFont(titlefont)
        ratio.GetXaxis().SetLabelSize(21)
        ratio.GetXaxis().SetLabelOffset(0.035)
        ratio.GetYaxis().CenterTitle()
        ratio.GetYaxis().SetTitleSize(22)
        ratio.GetYaxis().SetTitleFont(titlefont)
        ratio.GetYaxis().SetTitleOffset(2.2)
        ratio.GetYaxis().SetLabelFont(titlefont)
        ratio.GetYaxis().SetLabelSize(19)
        ratio.GetYaxis().SetLabelOffset(0.009)
        if isfirst: ratio.Draw("HIST")
        else :      ratio.Draw("HIST SAME")
        EFTratios.append(ratio)
        isfirst=False
    ratio_uncert = getRatio(h_sys_total,h_SMsignal)
    ratio_uncert.Draw("E2 SAME")
    for r in EFTratios: r.Draw("HIST SAME") #draw all rations again
    line = ROOT.TLine(xmin, 1, xmax, 1)
    line.SetLineWidth(2)
    line.SetLineColor(13)
    line.Draw()
    for ratio in EFTratios:
        ratio.Draw("HIST SAME")
    ROOT.gPad.RedrawAxis()
    c.Print(os.path.join(plot_directory, region+"__"+histname+"__"+WCname+".pdf"))
