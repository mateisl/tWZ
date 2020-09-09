#!/usr/bin/env python

import ROOT
from tWZ.Tools.user                      import HistogramsForFakeStudies


ffile = ROOT.TFile("{}/nVertexReweight.root".format(HistogramsForFakeStudies),"UPDATE")

data = ffile.Get("nVtxs_data")
QCD = ffile.Get("nVtxs_QCD")
Wjets = ffile.Get("nVtxs_WJetsToLNu")
TTbar = ffile.Get("nVtxs_TTbar")

bkg = QCD.Clone("SumOfBkgs")
bkg.Add(Wjets)
bkg.Add(TTbar)

print "data integral: {}".format(data.Integral())
print "bkg integral: {}".format(bkg.Integral())

bkg.Scale(data.Integral()/bkg.Integral())

print "normalize bkg"
print "data integral: {}".format(data.Integral())
print "bkg integral: {}".format(bkg.Integral())


bkg.Write()

nVtxWeight = data.Clone("nVtxWeightCalc")
nVtxWeight.Divide(bkg)

nVtxWeight.Write()


ffile.Close()
