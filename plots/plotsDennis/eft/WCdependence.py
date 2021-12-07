import ROOT
from MyRootTools.plotter.PlotWCDependence import PlotWCDependence
import Analysis.Tools.syncer
from tWZ.Tools.user                      import plot_directory

filenames = {
    "WZ":  "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1_noData/Run2018/all/trilepT-minDLmass12-onZ1-deepjet0-met60/Results.root",
    "ZZ":  "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1_noData/Run2018/all/trilepVL-minDLmass12-onZ1-onZ2-nLeptons4/Results.root",
    "ttZ": "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1_noData/Run2018/all/trilepT-minDLmass12-onZ1-njet4p-deepjet1p/Results.root",
}

regions = ["WZ", "ZZ", "ttZ"]

WCnames = ["cHq1Re11","cHq1Re22", "cHq1Re33", "cHq3Re11", "cHq3Re22", "cHq3Re33"]
processes = ["ttZ", "WZ", "ZZ"]
colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen]
legends = ["t#bar{t}Z", "WZ", "ZZ"]

histname = "Z1_pt"
bins = [0, 1, 2, 3, 4, 5]

for region in regions:
    for WCname in WCnames:
        for bin in bins:
            p = PlotWCDependence(filenames[region], histname, bin, WCname, region)
            p.plot_dir = plot_directory+"/WCdependence/"
            p.setDescription("Z p_{T}")
            p.doRebin(5)
            for i,pname in enumerate(processes):
                p.addProcess(pname, colors[i], legends[i])
            p.draw()
