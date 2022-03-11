#!/usr/bin/env python

import ROOT, os
import Analysis.Tools.syncer
from tWZ.Tools.user                      import plot_directory
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--mode', action='store', default='01')
argParser.add_argument('--process', action='store', default='ttZ')
argParser.add_argument('--latex', action='store_true', default=False)
args = argParser.parse_args()

filedir = "/users/dennis.schwarz/CMSSW_10_6_20/src/tWZ01j_LHE/"

files_ttZ = [
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
]

files_WZ = [
    ("djr_WZ_q10_20.root", "DJR_WZ_q10_20","xqcut=10, qCut=20"),
    ("djr_WZ_q10_25.root", "DJR_WZ_q10_25","xqcut=10, qCut=25"),
    ("djr_WZ_q10_30.root", "DJR_WZ_q10_30","xqcut=10, qCut=30"),
    ("djr_WZ_q20_25.root", "DJR_WZ_q20_25","xqcut=20, qCut=25"),
    ("djr_WZ_q20_30.root", "DJR_WZ_q20_30","xqcut=20, qCut=30"),
    # ("djr_WZ_q20_30_noEFT.root", "DJR_WZ_q20_30_noEFT","xqcut=20, qCut=30, noEFT"),
    ("djr_WZ_q20_35.root", "DJR_WZ_q20_35","xqcut=20, qCut=35"),
    ("djr_WZ_q20_40.root", "DJR_WZ_q20_40","xqcut=20, qCut=40"),
    ("djr_WZ_q20_45.root", "DJR_WZ_q20_45","xqcut=20, qCut=45"),
    ("djr_WZ_q20_50.root", "DJR_WZ_q20_50","xqcut=20, qCut=50"),
]

files_ZZ = [
    ("djr_ZZ_q10_20.root", "DJR_ZZ_q10_20","xqcut=10, qCut=20"),
    ("djr_ZZ_q10_25.root", "DJR_ZZ_q10_25","xqcut=10, qCut=25"),
    ("djr_ZZ_q10_30.root", "DJR_ZZ_q10_30","xqcut=10, qCut=30"),
    ("djr_ZZ_q20_25.root", "DJR_ZZ_q20_25","xqcut=20, qCut=25"),
    ("djr_ZZ_q20_30.root", "DJR_ZZ_q20_30","xqcut=20, qCut=30"),
    ("djr_ZZ_q20_35.root", "DJR_ZZ_q20_35","xqcut=20, qCut=35"),
    ("djr_ZZ_q20_40.root", "DJR_ZZ_q20_40","xqcut=20, qCut=40"),
    ("djr_ZZ_q20_45.root", "DJR_ZZ_q20_45","xqcut=20, qCut=45"),
    ("djr_ZZ_q20_50.root", "DJR_ZZ_q20_50","xqcut=20, qCut=50"),
]

files = files_ttZ
if args.process=="WZ": files = files_WZ
elif args.process=="ZZ": files = files_ZZ

dir = os.path.join(plot_directory, "DJR")
texname = dir+"/DJR"+args.mode+"_"+args.process+".tex"
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

Analysis.Tools.syncer.sync()
