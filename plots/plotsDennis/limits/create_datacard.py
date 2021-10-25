#!/usr/bin/env python

import ROOT, os, copy, time, sys

from tWZ.Tools.user                                      import combineReleaseLocation, limit_directory

from Analysis.Tools.cardFileWriter                       import cardFileWriter
from Analysis.Tools.cardFileWriter.CombineResults        import CombineResults


# Default Parameter
loggerChoices = ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "TRACE", "NOTSET"]

# Arguments
import argparse
argParser=argparse.ArgumentParser(description="Argument parser" )
argParser.add_argument( "--logLevel",           action="store",      default="INFO",            choices=loggerChoices,      help="Log level for logging" )
args=argParser.parse_args()

# Logging
import Analysis.Tools.logger as logger
logger = logger.get_logger(       args.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger( args.logLevel, logFile = None )

# Interpolate signal expectation
def calculateSignalPoint(newpoint,x1,x2,x3,y1,y2,y3):
    # y = ax^2 + bx + c
    a = (x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/((x1-x2)*(x1-x3)*(x3-x2))
    b = (x1*x1*(y2-y3)+x2*x2*(y3-y1)+x3*x3*(y1-y2))/((x1-x2)*(x1-x3)*(x2-x3))
    c = (x1*x1*(x2*y3-x3*y2)+x1*x1*(x3*x3*y2-x2*x2*y3)+x2*x3*y1*(x2-x3))/((x1-x2)*(x1-x3)*(x2-x3))
    y = a*newpoint*newpoint + b*newpoint + c
    return y

# Get binning from hist
def readBinning(hist, underflow, overflow):
    bins = []
    Nbins = hist.GetSize()-2
    firstbin = 1
    if underflow:
        firstbin=0
        Nbins+=1
    if overflow:
        Nbins+=1
    bins.append(hist.GetXaxis().GetBinLowEdge(firstbin))
    for i in range(Nbins):
        if i==(Nbins-1) and overflow:
            bins.append(-999)
        else:
            bins.append(hist.GetXaxis().GetBinUpEdge(i+firstbin))
    return bins

# Wrapper for dataCardWriter and combine
def writeCardFile(name):
    c = cardFileWriter.cardFileWriter()
    c.releaseLocation = combineReleaseLocation

    c.reset()
    c.setPrecision(3)
    # shapeString     = "lnN" if args.useTxt else "shape"
    shapeString     = "lnN"

    # Define regions and histname
    regions = ["WZ", "ZZ", "ttZ"]
    histname = "Z1_pt"

    dir_ZZ  = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/ZZ_EFT_v1_noData_SMEFTsim/Run2018/all/trilepVL-minDLmass12-onZ1-onZ2-nLeptons4/"
    dir_WZ  = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/WZ_EFT_v1_noData_SMEFTsim/Run2018/all/trilepT-minDLmass12-onZ1-met60/"
    dir_ttZ = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/ttZ_EFT_v1_noData_SMEFTsim/Run2018/all/trilepT-minDLmass12-onZ1-njet4p-deepjet1/"
    file_ZZ  = ROOT.TFile(dir_ZZ+'/Results.root')
    file_WZ  = ROOT.TFile(dir_WZ+'/Results.root')
    file_ttZ = ROOT.TFile(dir_ttZ+'/Results.root')

    # Define processes
    processes = ["tWZ", "ttZ", "ttX", "tZq", "WZ", "triBoson", "ZZ", "nonprompt"]
    signal = "cHq1Re11=2.00"

    hists_ZZ  = {}
    hists_WZ  = {}
    hists_ttZ = {}

    for p in processes:
        hists_ZZ[p]  = file_ZZ.Get(histname+"__"+p)
        hists_WZ[p]  = file_WZ.Get(histname+"__"+p)
        hists_ttZ[p] = file_ttZ.Get(histname+"__"+p)

    hists_ZZ["signal"]  = file_ZZ.Get(histname+"__ZZ__"+signal)
    hists_WZ["signal"]  = file_WZ.Get(histname+"__WZ__"+signal)
    hists_ttZ["signal"] = file_ttZ.Get(histname+"__ttZ__"+signal)

    allhists = {
        "ZZ":  hists_ZZ,
        "WZ":  hists_WZ,
        "ttZ": hists_ttZ
    }

    # Rebin
    for region in allhists:
        for process in allhists[region]:
            allhists[region][process].Rebin(10)

    # Define uncertainties
    c.addUncertainty( "MCstat", shapeString)
    c.addUncertainty( "test", shapeString)

    bincounter = 0
    for region in regions:
        bins = readBinning(allhists[region]["ttZ"],underflow=False,overflow=True)

        # Run over all bins and set values
        for bin, threshold in enumerate(bins):
            binname = "Bin%i"%bincounter
            c.addBin(binname, processes+["signal"])

            # Fill backgrounds and uncerts
            background_yield = 0
            for p in processes:
                # Fill bin content
                bincontent = allhists[region][p].GetBinContent(bin+1)
                if bincontent<0.0001:
                    bincontent = 0.1
                c.specifyExpectation( binname, p, bincontent)
                background_yield += bincontent
                # Fill MC stat uncert
                if bincontent < 0.001:
                    MCstat_unc = 0
                else:
                    allhists[region][p].GetBinError(bin+1)
                    MCstat_unc = allhists[region][p].GetBinError(bin+1)/bincontent
                c.specifyUncertainty( "MCstat", binname, p, 1 + MCstat_unc )
                # Test Uncertainty
                c.specifyUncertainty( "test", binname, p, 1.1 )

            # Fill Signal
            c.specifyExpectation( binname, "signal", allhists[region]["signal"].GetBinContent(bin+1) )
            # Fill Obervation (data or background sum)
            c.specifyObservation( binname, int( round( background_yield, 0 ) ) )
            bincounter+=1
            ####

    #### Write card files
    # Define Card File
    cardFileNameTxt   = os.path.join( limit_directory, name+".txt" )
    cardFileNameShape = cardFileNameTxt.replace( ".txt", "_shape.root" )
    cardFileNameTxt     = c.writeToFile( cardFileNameTxt )
    cardFileNameShape   = c.writeToShapeFile( cardFileNameShape )

writeCardFile("cardFile")
