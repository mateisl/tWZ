#!/usr/bin/env python

import ROOT, os
import Analysis.Tools.syncer
from tWZ.Tools.user                      import plot_directory


dir = "/users/dennis.schwarz/CMSSW_10_6_20/src/tWZ01j_LHE/"

files = [
    ("djr_q10_15.root", "DJR_xqcut10_qcut15","xqcut=10, qCut=15"),
    ("djr_q10_20.root", "DJR_xqcut10_qcut20","xqcut=10, qCut=20"),
    ("djr_q10_30.root", "DJR_xqcut10_qcut30","xqcut=10, qCut=30"),
    ("djr_q10_30_noEtaMax.root", "DJR_xqcut10_qcut30_noEtaMax","xqcut=20, setMad=on", "EtaMax = -1"),
    ("djr_q10_40.root", "DJR_xqcut10_qcut40","xqcut=10, qCut=40"),
    ("djr_q20_25.root", "DJR_xqcut20_qcut25","xqcut=20, qCut=25"),
    ("djr_q20_30.root", "DJR_xqcut20_qcut30","xqcut=20, qCut=30"),
    ("djr_q20_40.root", "DJR_xqcut20_qcut40","xqcut=20, qCut=40"),
    ("djr_q20_Mad.root", "DJR_xqcut20_setMad","xqcut=20, setMad=on"),
]


for (filename, name, legend) in files:
    file = ROOT.TFile(dir+filename)
    # canvas = ROOT.TCanvas()
    canvas = file.Get("djr01")
    infotext = ROOT.TLatex(3.5, 24, legend)
    infotext.SetNDC()
    infotext.SetTextAlign(13)
    infotext.SetTextFont(62)
    infotext.SetTextSize(0.05)
    infotext.SetX(0.18)
    infotext.SetY(0.84)
    infotext.Draw()
    canvas.Print(os.path.join(plot_directory, "DJR", name+".pdf"))
