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

def plotDistribution(dir, prefix, region, histname, xtitle, backgrounds, signals, WCname, SMpoint, EFTpoints, value_to_number, systs):
    filename_SM = "CombineInput_cHq1Re11_"+str(value_to_number[WCname][SMpoint])+".root"
    file = ROOT.TFile(dir+"/"+filename_SM)
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
        legname = bkgname
        if   bkgname == "ttZ": legname = "t#bar{t}Z"
        elif bkgname == "ttX": legname = "t#bar{t}X"
        leg.AddEntry(hist, legname, "f")

    # Get SM signal
    h_SMsignal = h_totalbkg.Clone() # Start with bkg and stack signal on top
    for signal in signals:
        hist = file.Get(region+"__"+histname+"/"+signal)
        h_SMsignal.Add(hist)
    h_SMsignal.SetLineWidth(2)
    h_SMsignal.SetLineColor(ROOT.kBlack)
    h_SMsignal.SetFillColor(0)

    # Copy h_SMsignal for nominal prediction
    h_totalMC = h_SMsignal.Clone()

    # Get Systematics
    h_sys_uncerts = [] # these will contian only the error ( (|central-up|+|central-down|) / 2)
    for sys in systs:
        # first get signals
        isFirst = True
        for signal in signals:
            hup = file.Get(region+"__"+histname+"/"+signal+"__"+sys+"Up")
            hdown = file.Get(region+"__"+histname+"/"+signal+"__"+sys+"Down")
            if isFirst:
                histup = hup.Clone()
                histdown = hdown.Clone()
                isFirst = False
            else:
                histup.Add(hup)
                histdown.Add(hdown)
        for bkg in backgrounds:
            hup = file.Get(region+"__"+histname+"/"+bkg+"__"+sys+"Up")
            hdown = file.Get(region+"__"+histname+"/"+bkg+"__"+sys+"Down")
            if isFirst:
                histup = hup.Clone()
                histdown = hdown.Clone()
                isFirst = False
            else:
                histup.Add(hup)
                histdown.Add(hdown)
        h_sys_uncerts.append( getSysUncert(histup, histdown, h_totalMC) )
    h_sys_total = getTotalUncert(h_totalMC, h_sys_uncerts) # this contains total uncertainty (pure uncertainty value)
    h_sys_total.SetFillStyle(3245)
    h_sys_total.SetFillColor(13)
    h_sys_total.SetLineWidth(0)

    leg.AddEntry(h_sys_total, "Total uncertainty","f")
    if len(signals): leg.AddEntry(h_SMsignal, 'SM point', "l")

    # Get EFT Signals
    h_EFTsignals = []
    filesEFT = []
    for i, EFTpoint in enumerate(EFTpoints):
        filename_EFT = "CombineInput_"+WCname+"_"+str(value_to_number[WCname][EFTpoint])+".root"
        filesEFT.append(ROOT.TFile(dir+"/"+filename_EFT))
        h_EFTsignals.append(h_totalbkg.Clone())
    for i, EFTpoint in enumerate(EFTpoints):
        for signal in signals:
            histeft = filesEFT[i].Get(region+"__"+histname+"/"+signal)
            h_EFTsignals[i].Add(histeft)
        h_EFTsignals[i].SetLineWidth(2)
        h_EFTsignals[i].SetLineColor(signal_colors[i])
        h_EFTsignals[i].SetFillColor(0)
        leg.AddEntry(h_EFTsignals[i], WCname+" = "+str(EFTpoint),"l")

    # Draw plot
    drawRatio = len(signals) > 0

    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c", "c", 600, 600)
    pady1 = 0.31 if drawRatio else 0.0
    pad1 = ROOT.TPad("pad1", "pad1", 0, pady1, 1, 1.0)
    if drawRatio:    pad1.SetBottomMargin(0.02)
    else:            pad1.SetBottomMargin(0.12)

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


    h_SMsignal.GetYaxis().SetTitle('Events')
    h_SMsignal.GetYaxis().SetTitleOffset(1.4)
    h_SMsignal.SetTitle('')
    h_SMsignal.SetMarkerStyle(0)

    h_SMsignal.Draw("HIST")

    for h in h_EFTsignals:
        h.Draw("HIST SAME")
    h_stack.Draw("HIST SAME")
    h_sys_total.SetMarkerStyle(0)
    h_sys_total.Draw("E2 SAME")

    leg.Draw()

    # some style things
    titlefont = 43

    # This is for correct axis titles
    if drawRatio:
        axis = ROOT.TGaxis( xmin, ymin, xmin, ymax, ymin, ymax, 505,"")
        h_SMsignal.GetXaxis().SetLabelSize(0.)
        h_SMsignal.GetYaxis().SetLabelSize(0.)
        axis.SetLabelOffset(0.01)
        axis.SetLabelFont(titlefont)
        axis.SetLabelSize(21)
        axis.SetNdivisions(505)
        axis.Draw()
    else:
        h_SMsignal.GetXaxis().SetTitle(xtitle)
        h_SMsignal.GetXaxis().SetNdivisions(505)
        h_SMsignal.GetXaxis().SetTickLength(0.07)
        h_SMsignal.GetXaxis().SetTitleSize(25)
        h_SMsignal.GetXaxis().SetTitleFont(titlefont)
        h_SMsignal.GetXaxis().SetTitleOffset(1.2)
        h_SMsignal.GetXaxis().SetLabelFont(titlefont)
        h_SMsignal.GetXaxis().SetLabelSize(21)
    ROOT.gPad.RedrawAxis()

    if drawRatio:
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
    plotname = os.path.join(plot_directory, "/CombineInput/", region+"__"+histname+"__"+WCname+".pdf")
    if prefix: plotname = os.path.join(plot_directory, "/CombineInput/", prefix+"__"+region+"__"+histname+".pdf")
    c.Print(plotname)
    file.Close()
    for f in filesEFT: f.Close()
