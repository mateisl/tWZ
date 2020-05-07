from tWZ.samples.GENandSIM_private import *

ROOT.gROOT.ProcessLine(".L $CMSSW_BASE/src/tWZ/samples/scripts/plotdjr.C")

for sample in [tOrtbar_WZ01j_OLRLL_LO]:#[ tWZ01j_rwgt_filter ]: #tOrtbar_WZ01j_OLRLL_LO ]:
    c = ROOT.TChain("Events")
    for file_ in sample.files[:50]:
        c.Add(file_)

    print "Plotting..."
    ROOT.plotdjr( c, sample.name )
