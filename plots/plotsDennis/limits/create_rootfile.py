#!/usr/bin/env python

import ROOT
from math                                import sqrt
from plotDistribution import plotDistribution
import Analysis.Tools.syncer

def simulateFlatSystematic(hist, size):
    sys = hist.Clone()
    Nbins = hist.GetSize()-2
    for i in range(Nbins):
        bin = i+1
        content = hist.GetBinContent(bin)
        sys.SetBinContent(bin, (1+size)*content)
    return sys

def simulateSystematic(hist, size):
    sys = hist.Clone()
    Nbins = hist.GetSize()-2
    start = -1 * size
    stop = size
    stepsize = (stop-start)/Nbins
    for i in range(Nbins):
        bin = i+1
        value = start + i*stepsize
        content = hist.GetBinContent(bin)
        sys.SetBinContent(bin, (1+value)*content)
    return sys

def removeNegative(hist):
    Nbins = hist.GetSize()-2
    for i in range(Nbins):
        bin = i+1
        content = hist.GetBinContent(bin)
        if content < 0:
            hist.SetBinContent(bin, 0.0)
    return hist
################################################################################
### Setup
regions = ["WZ", "ZZ", "ttZ"]
# regions = ["WZ", "ttZ"]
print 'Reading regions:', regions

# histname
histname = "Z1_pt"
print 'Reading Histogram:', histname

# Directories
dirs = {
    "ZZ":  "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/ZZ_EFT_simple_v1_noData_SMEFTsim/Run2018/all/trilepVL-minDLmass12-onZ1-onZ2-nLeptons4/",
    "WZ":  "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/WZ_EFT_simple_v1_noData_SMEFTsim/Run2018/all/trilepT-minDLmass12-onZ1-met60/",
    "ttZ": "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/ttZ_EFT_simple_v1_noData_SMEFTsim/Run2018/all/trilepT-minDLmass12-onZ1-njet4p-deepjet1/",
}

# Define backgrounds
backgrounds = {
    "ZZ":  ["tWZ", "ttZ", "ttX", "tZq", "WZ", "triBoson",       "nonprompt"],
    "WZ":  ["tWZ", "ttZ", "ttX", "tZq",       "triBoson", "ZZ", "nonprompt"],
    "ttZ": ["tWZ",        "ttX", "tZq", "WZ", "triBoson", "ZZ", "nonprompt"],
}


# Define Signal points
signalnames = []
WCs = ["cHq1Re11", "cHq1Re22", "cHq1Re33"]
minval = -10.0
maxval = 10.0
Npoints = 51
goodnames = {}
value_to_number = {}
for WCname in WCs:
    ntv = {}
    for i in range(Npoints):
        value = round(minval + ((maxval-minval)/(Npoints-1))*i,2)
        signalname='%s=%3.2f'%(WCname, value)
        signalnames.append(signalname)
        goodnames[signalname]="EFT_%s_%s"%(WCname,i)
        ntv[value] = i
    value_to_number[WCname] = ntv


# Define Systematics
sysnames = [
    'BTag_b',
    'BTag_l',
    'Trigger',
    'PU',
    'JES',
]


################################################################################
### Read Histograms and write to outfile
print 'Write output file...'
outname = '/mnt/hephy/cms/dennis.schwarz/www/tWZ/limits/CombineInput.root'
outfile = ROOT.TFile(outname, 'recreate')
outfile.cd()

for region in regions:
    print 'Filling region', region
    file = ROOT.TFile(dirs[region]+'/Results.root')
    outfile.mkdir(region+"__"+histname)
    outfile.cd(region+"__"+histname)
    # Get backgrounds
    is_first_bkg=True
    for bkg in backgrounds[region]:
        hist = file.Get(histname+"__"+bkg)
        # hist.Rebin(5)
        hist = removeNegative(hist)
        hist.Write(bkg)
        # Add all SM samples for expected limits
        if is_first_bkg: histsum = hist
        else:            histsum.Add(hist)
        is_first_bkg = False
        # Mock Systematics
        # simulateSystematic(hist,  0.4).Write(bkg+"__testsysUp")
        # simulateSystematic(hist, -0.4).Write(bkg+"__testsysDown")

        # Systematics
        for sys in sysnames:
            if 'test' in sys: continue
            sysdirUP = dirs[region]
            sysdirUP = sysdirUP.replace('/Run', '_'+sys+'_UP/Run')
            sysdirDOWN = dirs[region]
            sysdirDOWN = sysdirDOWN.replace('/Run', '_'+sys+'_DOWN/Run')
            fileUP   = ROOT.TFile(sysdirUP+'/Results.root')
            fileDOWN = ROOT.TFile(sysdirDOWN+'/Results.root')
            outfile.cd(region+"__"+histname)
            histUP   = fileUP.Get(histname+"__"+bkg)
            histDOWN = fileDOWN.Get(histname+"__"+bkg)
            # histUP.Rebin(5)
            # histDOWN.Rebin(5)
            histUP   = removeNegative(histUP)
            histDOWN = removeNegative(histDOWN)
            histUP.Write(bkg+"__"+sys+"Up")
            histDOWN.Write(bkg+"__"+sys+"Down")
    # Get signals
    for signalname in signalnames:
        hist = file.Get(histname+"__"+region+"__"+signalname)
        # hist.Rebin(5)
        # The Signal is now including all bkg processes
        hist.Add(histsum,-1) # Subtract background sum (which was added in the plotting script)
        hist = removeNegative(hist)
        hist.Write(goodnames[signalname])
        # Mock Systematics
        # simulateSystematic(hist,  0.4).Write(goodnames[signalname]+"__testsysUp")
        # simulateSystematic(hist, -0.4).Write(goodnames[signalname]+"__testsysDown")

        # Systematics
        for sys in sysnames:
            if 'test' in sys: continue
            sysdirUP = dirs[region]
            sysdirUP = sysdirUP.replace('/Run', '_'+sys+'_UP/Run')
            sysdirDOWN = dirs[region]
            sysdirDOWN = sysdirDOWN.replace('/Run', '_'+sys+'_DOWN/Run')
            fileUP   = ROOT.TFile(sysdirUP+'/Results.root')
            fileDOWN = ROOT.TFile(sysdirDOWN+'/Results.root')
            outfile.cd(region+"__"+histname)
            histUP   = fileUP.Get(histname+"__"+region+"__"+signalname)
            histDOWN = fileDOWN.Get(histname+"__"+region+"__"+signalname)
            # histUP.Rebin(5)
            # histDOWN.Rebin(5)
            histUP.Add(histsum,-1)
            histDOWN.Add(histsum,-1)
            histUP   = removeNegative(histUP)
            histDOWN = removeNegative(histDOWN)
            histUP.Write(goodnames[signalname]+"__"+sys+"Up")
            histDOWN.Write(goodnames[signalname]+"__"+sys+"Down")
    # Write observed
    # Add signal at SM point to sum of backgrounds
    observed=histsum.Clone()
    signalSM=file.Get(histname+"__"+region)
    # signalSM.Rebin(5)
    signalSM=removeNegative(signalSM)
    observed.Add(signalSM)
    observed.Write("data_obs")
    outfile.cd()
outfile.Close()
print 'Written to ', outname

for region in ["WZ", "ZZ", "ttZ"]:
    for WCname in WCs:
        SMpoint = 0
        EFTpoints = [-2.0, 2.0]
        plotDistribution(outname, region, 'Z1_pt', 'Z #it{p}_{T}', backgrounds[region], WCname, SMpoint, EFTpoints, value_to_number, sysnames)
Analysis.Tools.syncer.sync()
