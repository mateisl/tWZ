#!/usr/bin/env python

import ROOT
from math                                import sqrt


# import argparse
# argParser = argparse.ArgumentParser(description = "Argument parser")
# argParser.add_argument('--twoD',             action='store_true', default=False, help='No plots?')
# args = argParser.parse_args()

################################################################################
### Setup
regions = ["WZ", "ZZ"]
print 'Reading regions:', regions

# histname
histname = "Z1_pt"
print 'Reading Histogram:', histname

# Directories
dirs = {
    "ZZ":  "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1_noData/Run2018/all/trilepVL-minDLmass12-onZ1-onZ2-nLeptons4/",
    "WZ":  "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1_noData/Run2018/all/trilepT-minDLmass12-onZ1-deepjet0-met60/",
    "ttZ": "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1_noData/Run2018/all/trilepT-minDLmass12-onZ1-njet4p-deepjet1p/",
}

# Define backgrounds
processes = ["ttZ", "WZ", "ZZ", "tWZ", "ttX", "tZq", "triBoson", "nonprompt"]
signals = ["ttZ", "WZ", "ZZ"]

################################################################################
### Read Histograms
inname = 'Results.root'
for region in regions:
    integrals = {}
    integral_background = 0
    file = ROOT.TFile(dirs[region]+inname)
    for process in processes:
        if process in signals:
            h_c = file.Get(histname+"__"+process+"__cHq1Re11=0.0000")
        else:
            h_c = file.Get(histname+"__"+process)

        if process == region:
            integrals["Signal"] = h_c.Integral()
            h_v1 = file.Get(histname+"__"+process+"__cHq1Re11=2.0000")
            h_v2 = file.Get(histname+"__"+process+"__cHq1Re22=0.4000")
            integrals["Signal_Var1"] = h_v1.Integral()
            integrals["Signal_Var2"] = h_v2.Integral()
        else:
            integral_background += h_c.Integral()

    print "=========================="
    print region
    print "      - Signal fraction            =", integrals["Signal"]/(integrals["Signal"]+integral_background)
    print "      - Variation 1                =", abs(integrals["Signal"]-integrals["Signal_Var1"])/integrals["Signal"]
    print "      - Variation 1 (inkl bkg)     =", abs(integrals["Signal"]-integrals["Signal_Var1"])/(integrals["Signal"]+integral_background)
    print "      - Variation 2                =", abs(integrals["Signal"]-integrals["Signal_Var2"])/integrals["Signal"]
    print "      - Variation 2 (inkl bkg)     =", abs(integrals["Signal"]-integrals["Signal_Var2"])/(integrals["Signal"]+integral_background)
    # print integrals["Signal"], integrals["Signal_Var1"], integrals["Signal_Var2"], integral_background
