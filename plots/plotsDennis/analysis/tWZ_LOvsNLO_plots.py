#!/usr/bin/env python
import ROOT, os
import Analysis.Tools.syncer
from tWZ.Tools.user                              import plot_directory
from MyRootTools.plotter.Plotter                 import Plotter

################################################################################
# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
# argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
args = argParser.parse_args()

################################################################################

histnames = ["gen_Z1_pt", "gen_Z1_pt_bveto", "gen_Z1_pt_bveto_btag", "Z1_pt", "l1_pt"]
selections = [
    "minDLmass12",
    "onZ1",
    "trilepT-minDLmass12-onZ1",
    "trilepT-minDLmass12-onZ1-njet3p-deepjet1p",
    "trilepT-minDLmass12-onZ1-njet4p-deepjet1p",
    "trilepT-minDLmass12-onZ1-njet3p-deepjet1",
    "trilepT-minDLmass12-onZ1-njet4p-deepjet1",
]


xtitles = {
    "gen_Z1_pt": "Z p_{T}",
    "gen_Z1_pt_bveto": "Z p_{T}",
    "gen_Z1_pt_bveto_btag": "Z p_{T}",
    "Z1_pt":  "Z p_{T}",
    "l1_pt":  "Leading lepton p_{T}",
}

dir = plot_directory+"/tWZ_LOvsNLO/"

for selection in selections:
    file = ROOT.TFile("/groups/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/tWZ_LOvsNLO_v1/Run2018/all/"+selection+"/Results.root")
    for norm in [True, False]:
        for histname in histnames:
            LO = file.Get(histname+"__tWZ")
            NLO_DR = file.Get(histname+"__tWZ_NLO_DR")
            NLO_DS = file.Get(histname+"__tWZ_NLO_DS")
            if norm:
                LO.Scale(1/LO.Integral())
                NLO_DR.Scale(1/NLO_DR.Integral())
                NLO_DS.Scale(1/NLO_DS.Integral())
            # print "Factor =", LO.Integral()/NLO.Integral()
            plotname = "tWZ__"+histname+"__"+selection
            if norm:
                plotname = "tWZ__"+histname+"__NORM__"+selection
            p = Plotter(plotname)
            p.plot_dir = dir
            p.ytitle = "Events"
            if norm:
                p.ytitle = "a.u."                
            if histname in xtitles.keys():
                p.xtitle = xtitles[histname]
            p.drawRatio=True
            p.log = False
            p.addBackground(LO,      "tWZ LO",     13)
            p.addSignal(NLO_DR,      "tWZ NLO DR", ROOT.kRed)
            p.addSignal(NLO_DS,      "tWZ NLO DS", ROOT.kAzure+7)
            p.draw()
