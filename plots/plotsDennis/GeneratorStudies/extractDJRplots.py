#!/usr/bin/env python

import ROOT, os
import Analysis.Tools.syncer
from tWZ.Tools.user                      import plot_directory


import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--mode', action='store', default='01')
argParser.add_argument('--latex', action='store_true', default=False)
args = argParser.parse_args()

filedir = "/users/dennis.schwarz/CMSSW_10_6_20/src/tWZ01j_LHE/"

files = [
    ("djr_ttZ_q10_30.root", "DJR_ttZ_q10_30","xqcut=10, qCut=30"),
    ("djr_ttZ_q20_25.root", "DJR_ttZ_q20_25","xqcut=20, qCut=25"),
    ("djr_ttZ_q20_30.root", "DJR_ttZ_q20_30","xqcut=20, qCut=30"),
    ("djr_ttZ_q20_35.root", "DJR_ttZ_q20_35","xqcut=20, qCut=35"),
    ("djr_ttZ_q20_40.root", "DJR_ttZ_q20_40","xqcut=20, qCut=40"),
    ("djr_ttZ_q20_45.root", "DJR_ttZ_q20_45","xqcut=20, qCut=45"),
    ("djr_ttZ_q20_50.root", "DJR_ttZ_q20_50","xqcut=20, qCut=50"),
    ("djr_ttZ_q30_40.root", "DJR_ttZ_q30_40","xqcut=30, qCut=40"),
    ("djr_ttZ_q30_45.root", "DJR_ttZ_q30_45","xqcut=30, qCut=45"),
    ("djr_ttZ_q30_50.root", "DJR_ttZ_q30_50","xqcut=30, qCut=50"),
    ("djr_ttZ_q30_55.root", "DJR_ttZ_q30_55","xqcut=30, qCut=55"),
    ("djr_ttZ_q30_60.root", "DJR_ttZ_q30_60","xqcut=30, qCut=60"),
    ("djr_ttZ_q30_65.root", "DJR_ttZ_q30_65","xqcut=30, qCut=65"),
    # ("djr_q10_15.root", "DJR_xqcut10_qcut15","xqcut=10, qCut=15"),
    # ("djr_q10_20.root", "DJR_xqcut10_qcut20","xqcut=10, qCut=20"),
    # ("djr_q10_30.root", "DJR_xqcut10_qcut30","xqcut=10, qCut=30"),
    # ("djr_q10_30_noEtaMax.root", "DJR_xqcut10_qcut30_noEtaMax","xqcut=10, EtaMax = -1"),
    # ("djr_q10_40.root", "DJR_xqcut10_qcut40","xqcut=10, qCut=40"),
    # ("djr_q20_25.root", "DJR_xqcut20_qcut25","xqcut=20, qCut=25"),
    # ("djr_q20_30.root", "DJR_xqcut20_qcut30","xqcut=20, qCut=30"),
    # ("djr_q20_40.root", "DJR_xqcut20_qcut40","xqcut=20, qCut=40"),
    # ("djr_q20_Mad.root", "DJR_xqcut20_setMad","xqcut=20, setMad=on"),
]

dir = os.path.join(plot_directory, "DJR")
texname = dir+"/DJR"+args.mode+".tex"
outfile = open(texname,"w")
if args.latex:
    outfile.write("\\documentclass[aspectratio=169]{beamer}\n")
    outfile.write("\\usepackage[english]{babel}\n")
    outfile.write("\\usepackage{graphicx} \n")
    outfile.write("\\usepackage{epstopdf}\n")
    outfile.write("\\usetheme{Copenhagen}\n")
    outfile.write("\\usecolortheme{beaver}\n")
    outfile.write("\\begin{document}\n")

for (filename, name, legend) in files:
    file = ROOT.TFile(filedir+"/"+filename)
    # canvas = ROOT.TCanvas()
    canvas = file.Get("djr"+args.mode)
    infotext = ROOT.TLatex(3.5, 24, legend)
    infotext.SetNDC()
    infotext.SetTextAlign(13)
    infotext.SetTextFont(62)
    infotext.SetTextSize(0.05)
    infotext.SetX(0.18)
    infotext.SetY(0.84)
    infotext.Draw()
    outname = name.replace("DJR", "DJR"+args.mode)
    canvas.Print(dir+"/"+outname+".pdf")

    if  args.latex:
        outfile.write("%--------------------------------------------\n")
        outfile.write("\\begin{frame}\n")
        outfile.write("\\frametitle{"+outname.replace("_", " ")+"}\n")
        outfile.write("\\begin{figure}\n")
        outfile.write("\\centering\n")
        outfile.write("\\includegraphics[width=.5\\textwidth]{"+outname+".pdf} \n")
        outfile.write("\\end{figure}\n")
        outfile.write("\\end{frame}\n")

if args.latex:
    outfile.write("\\end{document}")
    outfile.close()
