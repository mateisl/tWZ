#!/usr/bin/env python

import ROOT
from math                                import sqrt
import array
from plotDistribution import plotDistribution
import Analysis.Tools.syncer


import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--plotOnly',         action='store_true', default=False, help='only plot without re-creating root file?')
argParser.add_argument('--noPlots',          action='store_true', default=False, help='No plots?')
argParser.add_argument('--twoD',             action='store_true', default=False, help='2D limits?')
argParser.add_argument('--noTWZ',            action='store_true', default=False, help='Keep tWZ at SM point?')
args = argParser.parse_args()

################################################################################
### Functions
def setPseudoDataErrors(hist):
    newhist = hist.Clone()
    Nbins = hist.GetSize()-2
    for i in range(Nbins):
        bin = i+1
        content = hist.GetBinContent(bin)
        newhist.SetBinError(bin, sqrt(content))
    return newhist

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

def setupHist(hist, bins):
    hist = hist.Rebin(len(bins)-1, "h", array.array('d',bins))
    hist = removeNegative(hist)
    return hist

################################################################################
### Setup
regions = ["WZ", "ZZ", "ttZ"]
print 'Reading regions:', regions


# binning 
bins  = [0, 60, 120, 180, 240, 300, 400, 1000]

# histname
histname = "Z1_pt"
print 'Reading Histogram:', histname

# Directories
dirs = {
    "ZZ":  "/groups/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1_noData/Run2018/all/trilepVL-minDLmass12-onZ1-onZ2-nLeptons4/",
    "WZ":  "/groups/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1_noData/Run2018/all/trilepT-minDLmass12-onZ1-deepjet0-met60/",
    "ttZ": "/groups/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1_noData/Run2018/all/trilepT-minDLmass12-onZ1-njet3p-deepjet1p/",
}
outdir = "/groups/hephy/cms/dennis.schwarz/www/tWZ/CombineInput/"
if args.noTWZ:
    outdir = "/groups/hephy/cms/dennis.schwarz/www/tWZ/CombineInput_noTWZ"

# Define backgrounds
processes = ["ttZ", "WZ", "ZZ", "tWZ", "ttX", "tZq", "triBoson", "nonprompt"]
backgrounds = ["ttX", "tZq", "triBoson", "nonprompt"]
signals = ["ttZ", "tWZ", "WZ", "ZZ"]
if args.noTWZ:
    backgrounds = ["tWZ","ttX", "tZq", "triBoson", "nonprompt"]
    signals = ["ttZ", "WZ", "ZZ"]

# Define Signal points
signalnames = []
WCs = ["cHq1Re11", "cHq1Re22", "cHq1Re33", "cHq3Re11", "cHq3Re22", "cHq3Re33"]
Npoints = 51
SMpointName = ""
goodnames = {}
value_to_number = {}
if not args.twoD:
    for WCname in WCs:
        minval = -10.0
        maxval = 10.0
        if "cHq3Re11" in WCname:
            minval = -0.2
            maxval = 0.2
        vtn = {}
        for i in range(Npoints):
            value = round(minval + ((maxval-minval)/(Npoints-1))*i,4)
            signalname='%s=%3.4f'%(WCname, value)
            if abs(value) < 0.001:
                print 'Found SM point'
                SMpointName = signalname
            signalnames.append(signalname)
            goodnames[signalname]="%s_%s"%(WCname,i)
            vtn[value] = i
        value_to_number[WCname] = vtn
elif args.twoD:
    minval1  = -4.0
    maxval1  = 4.0
    minval2  = -4.0
    maxval2  = 4.0
    Npoints1 = 21
    Npoints2 = 21
    WC1 = 'cHq1Re1122'
    WC2 = 'cHq1Re33'
    for i in range(Npoints1):
        value1 = minval1 + ((maxval1-minval1)/(Npoints1-1))*i
        for j in range(Npoints2):
            value2 = minval2 + ((maxval2-minval2)/(Npoints2-1))*j
            signalname='%s=%3.4f, %s=%3.4f'%(WC1,value1,WC2,value2)
            if abs(value1)<0.001 and abs(value2)<0.001:
                print 'Found SM point'
                SMpointName = signalname
            signalnames.append(signalname)
            goodnames[signalname]="%s_%i_%s_%i"%(WC1,i,WC2,j)

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
if args.twoD: inname = 'Results_twoD.root'
else:         inname = 'Results.root'

if not args.plotOnly:
    for signalpoint in signalnames:
        print '--------------------------------------------------------'
        print 'Working on', signalpoint
        outname = outdir+'/CombineInput_'+goodnames[signalpoint]+'.root'
        if args.twoD:
            outname = outdir+'/CombineInput_2D_'+goodnames[signalpoint]+'.root'
        outfile = ROOT.TFile(outname, 'recreate')
        outfile.cd()
        for region in regions:
            print 'Filling region', region
            file = ROOT.TFile(dirs[region]+inname)
            outfile.mkdir(region+"__"+histname)
            outfile.cd(region+"__"+histname)
            for process in processes:
                if process in backgrounds:
                    name = histname+"__"+process
                elif process in signals:
                    name = histname+"__"+process+"__"+signalpoint
                hist = file.Get(name)
                hist = setupHist(hist, bins)
                hist.Write(process)
                # Systematics
                for sys in sysnames:
                    sysdirUP = dirs[region]
                    sysdirUP = sysdirUP.replace('/Run', '_'+sys+'_UP/Run')
                    sysdirDOWN = dirs[region]
                    sysdirDOWN = sysdirDOWN.replace('/Run', '_'+sys+'_DOWN/Run')
                    fileUP   = ROOT.TFile(sysdirUP+inname)
                    fileDOWN = ROOT.TFile(sysdirDOWN+inname)
                    outfile.cd(region+"__"+histname)
                    histUP   = fileUP.Get(name)
                    histDOWN = fileDOWN.Get(name)
                    histUP   = setupHist(histUP, bins)
                    histDOWN = setupHist(histDOWN, bins)
                    histUP.Write(process+"__"+sys+"Up")
                    histDOWN.Write(process+"__"+sys+"Down")

            # Write observed
            # Add all relevant samples
            is_first = True
            observed = ROOT.TH1F()
            for bkg in backgrounds:
                hist = file.Get(histname+"__"+bkg)
                if is_first:
                    observed = hist.Clone()
                    is_first = False
                else:
                    observed.Add(hist)
            # now add all signals at SM point
            for signal in signals:
                hist = file.Get(histname+"__"+signal+"__"+SMpointName)
                observed.Add(hist)
            observed = setupHist(observed, bins)
            observed = setPseudoDataErrors(observed)
            observed.Write("data_obs")
        outfile.cd()
        outfile.Close()
        print 'Written to ', outname


if not args.noPlots and not args.twoD:
    signals_at_SM = signals
    for s in signals_at_SM:
        s = s+"__"+SMpointName
    for region in ["WZ", "ZZ", "ttZ"]:
        for WCname in WCs:
            SMpoint = 0
            EFTpoints = [-2, 2]
            if 'cHq3Re11' in WCname: EFTpoints = [-0.200, 0.200]
            if args.noTWZ:
                plotDistribution(outdir,"noTWZ", region, 'Z1_pt', 'Z #it{p}_{T}', backgrounds, signals, WCname, SMpoint, EFTpoints, value_to_number, sysnames)
            else:
                plotDistribution(outdir,None, region, 'Z1_pt', 'Z #it{p}_{T}', backgrounds, signals, WCname, SMpoint, EFTpoints, value_to_number, sysnames)

        if args.noTWZ:
            plotDistribution(outdir,"noTWZ__SM", region, 'Z1_pt', 'Z #it{p}_{T}', backgrounds+signals_at_SM, [], WCname, SMpoint, [], value_to_number, sysnames)
        else:
            plotDistribution(outdir,"SM", region, 'Z1_pt', 'Z #it{p}_{T}', backgrounds+signals_at_SM, [], WCname, SMpoint, [], value_to_number, sysnames)
            
    Analysis.Tools.syncer.sync()
