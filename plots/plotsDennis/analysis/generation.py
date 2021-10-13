import ROOT

regions = ["ZZ", "WZ", "ttZ"]

canvas_name = {
    "ZZ": "e9867383_b134_4df4_9ce4_e1de1dc9b0a5",
    "WZ": "93ed6b18_6d6a_413b_b984_fa2f191bd3db",
    "ttZ": "60b5f921_2f66_4c65_9cbb_d1b476eeaa02",
}

dirs = {
    "ZZ":  "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1_noData/Run2018/all/trilepVL-minDLmass12-onZ1-onZ2-nLeptons4/",
    "WZ":  "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1_noData/Run2018/all/trilepT-minDLmass12-onZ1-deepjet0-met60/",
    "ttZ": "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots/analysisPlots/EFT_SYS_v1_noData/Run2018/all/trilepT-minDLmass12-onZ1-njet4p-deepjet1p/",
}

filename = "ProductionMode.root"

Nbins = 8
bins = ["other", "1st", "2nd", "3rd", "1st+g", "2nd+g", "3rd+g", "gg"]

for region in regions:
    print region 
    file = ROOT.TFile(dirs[region]+filename)
    readObj = ROOT.TFileIter(file) # This gets you the first object in root file
    canvas = readObj

    hist = canvas.GetListOfPrimitives().At(1)
    h_other = canvas.GetListOfPrimitives().At(2)
    hist.Add(h_other,-1)
    for i in range(Nbins):
        bin = i+1
        print "    ", bins[i], "=", hist.GetBinContent(bin)/hist.Integral()
    print "-------------------------------------"
