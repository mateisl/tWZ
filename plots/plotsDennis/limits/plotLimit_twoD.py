#!/usr/bin/env python

import ROOT, os
import ctypes
import array
import Analysis.Tools.syncer
from tWZ.Tools.user                      import plot_directory


################################################################################
################################################################################
################################################################################

dir = "/users/dennis.schwarz/CMSSW_10_6_0/src/tWZ/plots/plotsDennis/DataCards/data_twoD/"
prefix = "higgsCombine.part3E_"
suffix = ".MultiDimFit.mH120.root"

WC1 = "cHq1Re1122"
WC2 = "cHq1Re33"

channels = {
"1": "ZZ",
"2": "WZ",
"3": "ttZ",
"combined": "combined"
}


for ch in channels:
    Npoints1 = 21
    Npoints2 = 21
    minq_point1 = -1
    minq_point2 = -1
    minqval = 10000000000
    qvals = []
    minval1 = -4.0
    maxval1 = 4.0
    minval2 = -4.0
    maxval2 = 4.0
    WC1vals = []
    WC2vals = []
    for i in range(Npoints1):
        value1 = minval1 + ((maxval1-minval1)/(Npoints1-1))*i
        for j in range(Npoints2):
            value2 = minval2 + ((maxval2-minval2)/(Npoints2-1))*j
            file = ROOT.TFile(dir+prefix+str(i)+"_"+str(j)+"_"+ch+"_"+WC1+"_"+WC2+suffix)
            tree = file.Get("limit")
            if tree.GetEntry(0)<=0:
                print 'EMPTY TREE, leave q unchanged...'
                qvals.append(q)
                WC1vals.append(value1)
                WC2vals.append(value2)
            for event in file.limit:
                q = 2*(event.nll+event.nll0)
                qvals.append(q)
                WC1vals.append(value1)
                WC2vals.append(value2)
                if q<minqval:
                    minqval = q
                    minq_point1 = value1
                    minq_point2 = value2

    qdiff = []
    for q in qvals:
        qdiff.append(q-minqval)

    graph = ROOT.TGraph2D( Npoints1*Npoints2, array.array('d',WC1vals), array.array('d',WC2vals), array.array('d',qdiff) )

    ############################################################################
    # Set up Canvas
    c=ROOT.TCanvas("", "", 800, 600)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gPad.SetRightMargin(0.2)
    ROOT.gPad.SetLeftMargin(0.1)


    ############################################################################
    # Convert to 2D hist and
    nxbins = Npoints1*2
    nybins = Npoints2*2
    hist = graph.GetHistogram().Clone()
    graph.SetNpx( nxbins )
    graph.SetNpy( nybins )
    hist = graph.GetHistogram().Clone()

    #smoothing
    # if args.smooth: hist.Smooth()

    contours = [2.28, 5.99]# (68%, 95%) for 2D
    histsForCont = hist.Clone()
    c_contlist = ((ctypes.c_double)*(len(contours)))(*contours)
    histsForCont.SetContour(len(c_contlist),c_contlist)
    histsForCont.Draw("contzlist")
    c.Update()
    conts = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    cont_p1 = conts.At(0).Clone()
    cont_p2 = conts.At(1).Clone()


    ############################################################################
    # Start plotting
    hist.Draw("colz")
    hist.SetTitle("")
    hist.GetXaxis().SetTitle(WC1)
    hist.GetYaxis().SetTitle(WC2)
    hist.GetZaxis().SetTitle("-2 #Delta ln L")

    hist.GetYaxis().SetRangeUser(-4.0, 4.0)
    hist.GetZaxis().SetTitleOffset(1.3)


    #draw contour lines
    for conts in [cont_p2]:
        for cont in conts:
            cont.SetLineColor(ROOT.kOrange+7)
            cont.SetLineWidth(3)
            cont.Draw("same")
    for conts in [cont_p1]:
        for cont in conts:
            cont.SetLineColor(ROOT.kSpring-1)
            cont.SetLineWidth(3)
            cont.Draw("same")

    # Draw SM point
    SMpoint = ROOT.TGraph(1)
    SMpoint.SetName("SMpoint")
    SMpoint.SetPoint(0, 0, 0)
    SMpoint.SetMarkerStyle(20)
    SMpoint.SetMarkerSize(2)
    SMpoint.SetMarkerColor(ROOT.kOrange+1)
    SMpoint.Draw("p same")

    # Draw best fit point
    BFpoint = ROOT.TGraph(1)
    BFpoint.SetName("BFpoint")
    BFpoint.SetPoint(0, minq_point1, minq_point2)
    BFpoint.SetMarkerStyle(29)
    BFpoint.SetMarkerSize(3)
    BFpoint.SetMarkerColor(ROOT.kCyan-3)
    BFpoint.Draw("p same")


    ############################################################################
    # Legend
    leg = ROOT.TLegend(.5, .6, .75, .85)
    leg.AddEntry( BFpoint, "Best fit","p")
    leg.AddEntry( SMpoint, "SM","p")
    leg.AddEntry( cont_p1.At(0), "68%s CL"%"%", "l")
    leg.AddEntry( cont_p2.At(0), "95%s CL"%"%", "l")
    leg.Draw()

    ROOT.gPad.RedrawAxis()
    c.Print(os.path.join(plot_directory+"/Limits/", "Limit_2D_"+channels[ch]+"_"+WC1+"_"+WC2+".pdf"))


Analysis.Tools.syncer.sync()
