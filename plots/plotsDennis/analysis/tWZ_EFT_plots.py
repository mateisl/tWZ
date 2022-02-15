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
argParser.add_argument('--latex',  action='store_true', default=False, help='Create slides with latex?', )
args = argParser.parse_args()

################################################################################

histnames = ["gen_Z1_pt", "gen_m_WZ", "gen_Z1_pt_bveto", "gen_m_WZ_bveto", "gen_Z1_pt_bveto_btag", "gen_m_WZ_bveto_btag"]
# file = ROOT.TFile("/groups/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_GEN_v1/Run2018/all/trilepVL-minDLmass12-onZ1/Results.root")
file = ROOT.TFile("/groups/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_GEN_v1/Run2018/all/onZ1/Results.root")
file2 = ROOT.TFile("/groups/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_GEN_v1/Run2018/all/offZ1/Results.root")

xtitles = {
    "gen_Z1_pt": "Z p_{T}",
    "gen_m_WZ": "m_{WZ}",
    "gen_Z1_pt_bveto": "Z p_{T}",
    "gen_m_WZ_bveto": "m_{WZ}",
    "gen_Z1_pt_bveto_btag": "Z p_{T}",
    "gen_m_WZ_bveto_btag": "m_{WZ}",
}

WCnames = [
    'cHq1Re11',
    'cHq1Re22',
    'cHq1Re33',
    'cHq3Re11',
    'cHq3Re22',
    'cHq3Re33',
    'cHuRe11',
    'cHuRe22',
    'cHuRe33',
    'cHdRe11',
    'cHdRe22',
    'cHdRe33',
    'cW'     ,
    'cWtil'  ,
]

processes = ["tWZ", "ttZ"]

dir = plot_directory+"/tWZ_EFT_studies/"

outfile = open(dir+"/AllPlots.tex","w")
if args.latex:
    outfile.write("\\documentclass[aspectratio=169]{beamer}\n")
    outfile.write("\\usepackage[english]{babel}\n")
    outfile.write("\\usepackage{graphicx} \n")
    outfile.write("\\usepackage{epstopdf}\n")
    outfile.write("\\usetheme{Copenhagen}\n")
    outfile.write("\\usecolortheme{beaver}\n")
    outfile.write("\\begin{document}\n")

for histname in histnames:
    for WCname in WCnames:
        if args.latex:
            nounderscorename = histname.replace("_","")
            outfile.write("%--------------------------------------------\n")
            outfile.write("\\begin{frame}\n")
            outfile.write("\\frametitle{"+nounderscorename+", "+WCname+"}\n")
            outfile.write("\\begin{columns}\n")
        for process in processes:
            if "m_WZ" in histname and process == "ttZ":
                continue
            # print histname+"__"+process+"__"+WCname+"=0.0000"
            central = file.Get(histname+"__"+process+"__"+WCname+"=0.0000")
            central2 = file2.Get(histname+"__"+process+"__"+WCname+"=0.0000")
            central.Add(central2)
            central.Scale(1/(60*central.GetBinWidth(1)))
            down = file.Get(histname+"__"+process+"__"+WCname+"=-1.0000")
            down2 = file2.Get(histname+"__"+process+"__"+WCname+"=-1.0000")
            down.Add(down2)
            down.Scale(1/(60*central.GetBinWidth(1)))
            up = file.Get(histname+"__"+process+"__"+WCname+"=1.0000")
            up2 = file2.Get(histname+"__"+process+"__"+WCname+"=1.0000")
            up.Add(up2)
            up.Scale(1/(60*central.GetBinWidth(1)))
            plotname = process+"__"+histname+"__"+WCname
            p = Plotter(plotname)
            p.plot_dir = dir
            p.ytitle = "#sigma [pb/GeV]"
            if histname in xtitles.keys():
                p.xtitle = xtitles[histname]
            p.drawRatio=True
            p.log = True
            p.addBackground(central, process+" SM",   13)
            p.addSignal(up,      process+" "+WCname+" =  1.0", ROOT.kRed)
            p.addSignal(down,    process+" "+WCname+" = -1.0", ROOT.kBlue)
            p.draw()
            if args.latex:
                outfile.write("\\begin{column}{.5\\textwidth}\n")
                outfile.write("\\begin{figure}\n")
                outfile.write("\\includegraphics[width=\\textwidth]{"+plotname+".pdf} \n")
                outfile.write("\\end{figure}\n")
                outfile.write("\\end{column}\n")
        if args.latex:
            outfile.write("\\end{columns}\n")
            outfile.write("\\end{frame}\n")

if args.latex:
    outfile.write("\\end{document}")
    outfile.close()
    # os.chdir(dir)
    # os.system("pdflatex AllPlots.tex")
