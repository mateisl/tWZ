#!/usr/bin/env python

import ROOT, os
import array
import Analysis.Tools.syncer
from tWZ.Tools.user                      import plot_directory

def CreateBands(graph, lower, upper):
    # If only one boundary exists, set other limit to end of the graph
    if lower==upper:
        if lower < 0:
            upper = max(graph.GetX())
        else:
            lower= min(graph.GetX())
    # First get x and y values according to graph
    stepsize = 0.01
    x = lower
    xvals = []
    yvals = []
    while x < upper:
        xvals.append(x)
        yvals.append(graph.Eval(x))
        x += stepsize
    # Store xvals and loop through the reverse again and write a point (x,0)
    # Has to be reverse to get the form
    # [(x1,y1), (x2,y2), ..., (xN,yN), (xN,0), ..., (x2,0), (x1,0)]
    # (0 is chosen since the plot is defined to have a minimum at 0)
    xvalsold=xvals
    for xval in reversed(xvalsold):
        xvals.append(xval)
        yvals.append(0.0)
    band = ROOT.TGraph(len(xvals), array.array('f',xvals), array.array('f',yvals))
    return band

def findXBoundaries(graph, minval, maxval, ymax):
    # If graph at minval/maxval is already greater than ymax, return minval/maxval
    # Otherwise return values where graph is for the first time < ymax
    lo = minval
    hi = maxval
    stepsize=0.01
    # Find point starting from minval
    if graph.Eval(minval) > ymax:
        x=minval
        while x < maxval:
            if graph.Eval(x)<ymax:
                lo = x
                break
            x+=stepsize
    # Find other point starting from maxval
    if graph.Eval(maxval) > ymax:
        x=maxval
        while x > minval:
            if graph.Eval(x)<ymax:
                hi = x
                break
            x-=stepsize
    return lo, hi

################################################################################
################################################################################
################################################################################

dir = "/users/dennis.schwarz/CMSSW_10_6_0/src/tWZ/plots/plotsDennis/DataCards/data/"
prefix = "higgsCombine.part3E_"
suffix = ".MultiDimFit.mH120.root"

WCnames = ["cHq1Re11", "cHq1Re22", "cHq1Re33", "cHq3Re11", "cHq3Re22", "cHq3Re33"]

channels = {
"1": "ZZ",
"2": "WZ",
"3": "ttZ",
"combined": "combined"
}


for WCname in WCnames:
    likelihoods = {}
    for ch in channels:
        Npoints = 51
        i_SMpoint = 25
        qval_SM = -1
        qvals = []
        qmin = 1000000000
        minval = -10.0
        maxval = 10.0
        if "cHq3Re11" in WCname:
            minval = -0.2
            maxval = 0.2
        WCvals = []
        for i in range(Npoints):
            value = minval + ((maxval-minval)/(Npoints-1))*i
            file = ROOT.TFile(dir+prefix+str(i)+"_"+ch+"_"+WCname+suffix)
            tree = file.Get("limit")
            if tree.GetEntry(0)<=0:
                print 'EMPTY TREE, leave q unchanged...'
                qvals.append(q)
                WCvals.append(value)
            for event in file.limit:
                # print "---------------"
                # if i==i_SMpoint: print "SM point"
                # print "nll  =", event.nll
                # print "nll0 =", event.nll0
                # Option 1
                q = 2*(event.nll+event.nll0)
                qvals.append(q)
                WCvals.append(value)
                if q<qmin:
                    qmin = q
                if i==i_SMpoint:
                    qval_SM = q
                # Option 2
                # if event.deltaNLL != 0:
                #     q = 2*event.deltaNLL
                #     qvals.append(q)
                #     if i==i_SMpoint:
                #         qval_SM = q

        qdiff = []
        for q in qvals:
            qdiff.append(q-qmin)

        ############################################################################
        # Start plotting
        graph = ROOT.TGraph( Npoints, array.array('d',WCvals), array.array('d',qdiff) )
        likelihoods[channels[ch]] = graph
        c=ROOT.TCanvas()
        ROOT.gStyle.SetLegendBorderSize(0)
        ROOT.gStyle.SetPadTickX(1)
        ROOT.gStyle.SetPadTickY(1)

        ymax = 10
        xlo, xhi = findXBoundaries(graph, minval, maxval, ymax)
        graph.GetXaxis().SetLimits(xlo, xhi)
        graph.GetHistogram().SetMinimum(0.)
        graph.GetHistogram().SetMaximum(ymax)
        graph.Draw('AL')
        graph.SetTitle('')
        graph.GetXaxis().SetTitle(WCname)
        graph.GetYaxis().SetTitle('-2 #Delta ln L')

        ############################################################################
        # Search for 65%CL and 95%CL points
        lines68=[]
        lines95=[]
        stepsize = 0.001
        accuracy = 0.01
        searching = True
        x = minval
        CL68vals = []
        CL95vals = []
        while searching:
            x += stepsize
            y = graph.Eval(x)
            # 68% CL
            if abs(y-1.00) < accuracy:
                lines68.append(ROOT.TLine(x,0.0, x, y))
                CL68vals.append(x)
                # Now go to point where if condition does not hold anymore
                # Otherwise one would have a huge number of lines that refer to the
                # same point
                samepoint=True
                while samepoint:
                    x+=stepsize
                    y = graph.Eval(x)
                    if abs(y-1.00) > 2*accuracy or x > maxval:
                        samepoint=False
            # 95% CL
            if abs(y-3.84) < accuracy:
                lines95.append(ROOT.TLine(x,0.0, x, y))
                CL95vals.append(x)
                samepoint=True
                while samepoint:
                    x+=stepsize
                    y = graph.Eval(x)
                    if abs(y-3.84) > 2*accuracy or x > maxval:
                        samepoint=False
            if x > maxval:
                searching=False


        ############################################################################
        # Draw 68% and 95% CL lines
        # for line in lines68:
        #     line.SetLineColor(ROOT.kRed)
        #     line.SetLineStyle(2)
        #     line.SetLineWidth(2)
        #     line.Draw("SAME")
        # for line in lines95:
        #     line.SetLineColor(ROOT.kRed)
        #     line.SetLineStyle(1)
        #     line.SetLineWidth(2)
        #     line.Draw("SAME")

        # Also draw bands
        if CL95vals: CL95band=CreateBands(graph, min(CL95vals), max(CL95vals))
        else:        CL95band=CreateBands(graph, minval, maxval)
        CL95band.SetFillColor(ROOT.kOrange);
        CL95band.Draw("f");
        if CL68vals: CL68band=CreateBands(graph, min(CL68vals), max(CL68vals))
        else:        CL68band=CreateBands(graph, minval, maxval)
        CL68band.SetFillColor(ROOT.kGreen+1);
        CL68band.Draw("f");
        ############################################################################
        # Horizontal lines at 1 and 3.84
        linesCL = [
            ROOT.TLine(xlo, 1,    xhi, 1),
            ROOT.TLine(xlo, 3.84, xhi, 3.84),
        ]
        for line in linesCL:
            line.SetLineColor(13)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw("SAME")

        ############################################################################
        # Legend
        leg = ROOT.TLegend(.5, .6, .85, .85)
        leg.AddEntry(graph, "Expected", "l")
        leg.AddEntry(CL68band, "68% CL", "f")
        leg.AddEntry(CL95band, "95% CL", "f")
        leg.Draw()

        ROOT.gPad.RedrawAxis()
        c.Print(os.path.join(plot_directory, "Limit_"+channels[ch]+"_"+WCname+".pdf"))
    ############################################################################
    # Compare Likelihoods of different channels
    if likelihoods.has_key('combined') and likelihoods.has_key('ttZ') and likelihoods.has_key('ZZ') and likelihoods.has_key('WZ'):
        c=ROOT.TCanvas()
        ROOT.gStyle.SetLegendBorderSize(0)
        ROOT.gStyle.SetPadTickX(1)
        ROOT.gStyle.SetPadTickY(1)

        likelihoods["combined"].GetXaxis().SetLimits(-10, 10)
        if "cHq3Re11" in WCname: likelihoods["combined"].GetXaxis().SetLimits(-0.2, 0.2)
        likelihoods["combined"].GetHistogram().SetMinimum(0.)
        likelihoods["combined"].GetHistogram().SetMaximum(20)
        likelihoods["combined"].SetLineColor(ROOT.kBlack)
        likelihoods["combined"].SetLineWidth(2)
        likelihoods["combined"].Draw('AL')
        likelihoods["combined"].SetTitle('')
        likelihoods["combined"].GetXaxis().SetTitle(WCname)
        likelihoods["combined"].GetYaxis().SetTitle('-2 #Delta ln L')
        likelihoods["ttZ"].SetLineWidth(2)
        likelihoods["ttZ"].SetLineColor(ROOT.kAzure+4)
        likelihoods["ttZ"].Draw('L SAME')
        likelihoods["ZZ"].SetLineWidth(2)
        likelihoods["ZZ"].SetLineColor(ROOT.kGreen+3)
        likelihoods["ZZ"].Draw('L SAME')
        likelihoods["WZ"].SetLineWidth(2)
        likelihoods["WZ"].SetLineColor(ROOT.kRed)
        likelihoods["WZ"].Draw('L SAME')

        linesCL = [
            ROOT.TLine(-10, 1,    10, 1),
            ROOT.TLine(-10, 3.84, 10, 3.84),
        ]
        if "cHq3Re11" in WCname:
            linesCL = [
                ROOT.TLine(-0.2, 1,    0.2, 1),
                ROOT.TLine(-0.2, 3.84, 0.2, 3.84),
            ]

        for line in linesCL:
            line.SetLineColor(13)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw("SAME")

        leg = ROOT.TLegend(.5, .6, .85, .85)
        leg.AddEntry(likelihoods["combined"], "Combination", "l")
        leg.AddEntry(likelihoods["ttZ"], "ttZ", "l")
        leg.AddEntry(likelihoods["ZZ"], "ZZ", "l")
        leg.AddEntry(likelihoods["WZ"], "WZ", "l")
        leg.Draw()

        ROOT.gPad.RedrawAxis()
        c.Print(os.path.join(plot_directory, "Limit_comparison_"+WCname+".pdf"))
    else:
        print 'Comparison plot is not created since not all regions are filled.'

Analysis.Tools.syncer.sync()
