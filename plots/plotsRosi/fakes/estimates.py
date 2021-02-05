#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT, os
ROOT.gROOT.SetBatch(True)
from math                                import sqrt, cos, sin, pi, atan2, cosh

# RootTools
from RootTools.core.standard             import *

# tWZ
from tWZ.Tools.user                      import plot_directory, cache_dir
from tWZ.Tools.helpers                   import getObjDict, getVarValue
from tWZ.Tools.cutInterpreter            import cutInterpreter
# Analysis
from Analysis.Tools.helpers              import deltaPhi, deltaR
from Analysis.Tools.metFilters           import getFilterCut
from Analysis.Tools.DirDB                import DirDB
import Analysis.Tools.syncer             as syncer

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
#argParser.add_argument('--noData',         action='store_true', default=False, help='also plot data?')
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
argParser.add_argument('--overwrite',                         action='store_true', help='Overwrite cache?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
#argParser.add_argument('--dataMCScaling',  action='store_true', help='Data MC scaling?', )
argParser.add_argument('--plot_directory', action='store', default='tWZ_fake_estimates') # -v2')
argParser.add_argument('--era',            action='store', type=str, default="Run2016")
argParser.add_argument('--mode',           action='store', type=str, default="mu", choices=["mu","ele"])
argParser.add_argument('--mT',             action='store', type=int, default=-1)
argParser.add_argument('--met',            action='store', type=int, default=-1)
args = argParser.parse_args()

# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"
#if args.noData:                       args.plot_directory += "_noData"

logger.info( "Working in era %s", args.era)
year = int(args.era[3:]) # this is designed to fail for era "RunII"

if args.era == "Run2016":
    import tWZ.samples.nanoTuples_fakes_2016_nanoAODv6_private_postProcessed as samples
    if args.mode=='mu':
        data_sample =  samples.DoubleMuon_Run2016
        triggers    = ["HLT_Mu3_PFJet40" ] #, "HLT_Mu8", "HLT_Mu17"]#, "HLT_Mu27"] HLT_Mu27 is actually in SingleMuon!
        mc = [ samples.QCD_mu, samples.WJetsToLNu_mu, samples.TTbar_mu]  #samples.QCD_pt_mu]
    elif args.mode=='ele':
        data_sample =  samples.DoubleEG_Run2016
        triggers    = ["HLT_Ele8_CaloIdM_TrackIdM_PFJet30" ]
        mc = [ samples.QCD_ele, samples.WJetsToLNu_ele, samples.TTbar_ele] # samples.QCD_pt_ele

bins = [ 
    (2.5,5, -1.2, 1.2),
    (5.,10.,  -1.2, 1.2),
    (10.,999., -1.2, 1.2),
    ]

triggerSelection = '('+"||".join(triggers)+')'
leptonSelection  = 'n%s_FOmvaTOPT==1'%args.mode
jetSelection     = 'Sum$(Jet_pt>40&&abs(Jet_eta)<2.4&&JetGood_cleaned_%s_FOmvaTOPT)>=1'%args.mode

selectionName    = "inclusive"
eventSelection   = "(1)"
if args.met>0:   
    selectionName = "metTo{met}".format(met=args.met)
    eventSelection= "met_pt<={met}".format(met=args.met)
if args.mT>0:
    mT_str = "Sum$({mode}_FOmvaTOPT_mT<{mT})==1".format(**args.__dict__) 
    if args.met>0:
        eventSelection+= "&&"+mT_str
        selectionName += "-mTTo{mT}".format(mT=args.mT)
    else: 
        eventSelection = mT_str 
        selectionName  = "mTTo{mT}".format(mT=args.mT)

logger.info("selection %s -> %s", selectionName, eventSelection)

data_selectionString = "&&".join([getFilterCut(isData=True, year=year), triggerSelection, leptonSelection, jetSelection, eventSelection])
data_sample.setSelectionString( data_selectionString )
mc_selectionString   = "&&".join([getFilterCut(isData=False, year=year), triggerSelection, leptonSelection, jetSelection, eventSelection])
for s in mc:
    s.setSelectionString( mc_selectionString )

max_events = -1
for sample in [data_sample] + mc:
    # remove 'tree_...' files until they are deleted
    sample.files = filter( lambda f:not f.split('/')[-1].startswith('tree_'), sample.files ) 
if args.small:
    for sample in [data_sample] + mc:
        sample.normalization = 1.
        #sample.reduceFiles( factor = 10 )
        sample.reduceFiles( factor=50 )
        #sample.scale /= sample.normalization
        #max_events = 100000

# Text on the plots
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)
tex.SetTextAlign(11) # align right

# fire up the cache
cache_dir_ = os.path.join(cache_dir, 'fake_pu_cache')
dirDB      = DirDB(cache_dir_)

pu_key = ( triggerSelection, leptonSelection, jetSelection, args.era, args.small)
if dirDB.contains( pu_key ) and not args.overwrite:
    reweight_histo = dirDB.get( pu_key )
    logger.info( "Found PU reweight in cache %s", cache_dir_ )
else:
    logger.info( "Didn't find PU reweight histo %r. Obtaining it now.", pu_key)

    data_nvtx_histo = data_sample.get1DHistoFromDraw( "PV_npvsGood", [100, 0, 100], weightString = "weight" )
    data_nvtx_histo.Scale(1./data_nvtx_histo.Integral())

    mc_histos  = [ s.get1DHistoFromDraw( "PV_npvsGood", [100, 0, 100], weightString = "weight*reweightBTag_SF") for s in mc]
    mc_histo     = mc_histos[0]
    for h in mc_histos[1:]:
        mc_histo.Add( h )

    mc_histo.Scale(1./mc_histo.Integral())

    reweight_histo = data_nvtx_histo.Clone()
    reweight_histo.Divide(mc_histo)
    
    dirDB.add( pu_key, reweight_histo ) 
    logger.info( "Added PU reweight to cache %s", cache_dir_ )

def nvtx_puRW( event, sample ):
    return reweight_histo.GetBinContent(reweight_histo.FindBin( event.PV_npvsGood ))

#lumi_scale                 = data_sample.lumi/1000
data_sample.scale   = 1.
for sample in mc:
    sample.weight   = nvtx_puRW

def drawObjects():
    lines = [
      (0.15, 0.95, 'CMS Preliminary'), 
      #(0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines] 

def drawPlots(plots):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, selectionName, args.mode, ("log" if log else "lin"))
    for plot in plots:
      if not max(l.GetMaximum() for l in sum(plot.histos,[])): continue # Empty plot

      _drawObjects = []

      if isinstance( plot, Plot):
          plotting.draw(plot,
            plot_directory = plot_directory_,
            ratio = {},
            logX = False, logY = log, sorting = True,
            yRange = (0.03, "auto") if log else (0.001, "auto"),
            #scaling = {0:1},
            legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
            drawObjects = drawObjects() + _drawObjects,
            copyIndexPHP = True, extensions = ["png"],
          )
            
# Read variables and sequences

read_variables = [
    "weight/F", "PV_npvsGood/I",
    "JetGood[pt/F,phi/F,eta/F]",
#    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
#    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
    ]

#read_variables += ["n%s_looseHybridIso/I"%args.mode, "%s_looseHybridIso[pt/F,eta/F,phi/F,mT/F,hybridIso/F]"%args.mode, "met_pt/F"]
read_variables += ["n%s_FOmvaTOPT/I"%args.mode, "%s_FOmvaTOPT[pt/F,eta/F,phi/F,mT/F]"%args.mode, "met_pt/F", "nmu_mvaTOPT/I", "nele_mvaTOPT/I"]

sequence       = []
def makeLeptons( event, sample ):
    collVars = ["pt","eta","phi","mT"]
    lep  = getObjDict(event, args.mode+'_FOmvaTOPT_', collVars, 0)
    for var in collVars:
        setattr( event, "lep_"+var, lep[var]  )
sequence.append( makeLeptons )

allPlots   = {}

data_sample.texName = "data"
data_sample.name    = "data"
data_sample.style   = styles.errorStyle(ROOT.kBlack)
#lumi_scale          = data_sample.lumi/1000


for sample in mc: sample.style = styles.fillStyle(sample.color)

stack = Stack(mc, [data_sample])

weight_ = lambda event, sample: event.weight 
# Use some defaults
Plot.setDefaults(stack = stack, weight = staticmethod(weight_))

plots = []

plots.append(Plot(
  name = 'nVtxs', texX = 'vertex multiplicity', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "PV_npvsGood/I" ),
  binning=[50,0,50],
  addOverFlowBin='upper',
))

plots.append(Plot(
  name = 'met_pt', texX = 'MET', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "met_pt/F" ),
  binning=[50,0,250],
  addOverFlowBin='upper',
))

plots.append(Plot(
  name = 'mT', texX = 'm_{T}', texY = 'Number of Events',
  attribute = lambda event, sample: event.lep_mT,
  binning=[40,0,200],
  addOverFlowBin='upper',
))

plots.append(Plot(
  name = 'pt', texX = 'p_{T}', texY = 'Number of Events',
  attribute = lambda event, sample: event.lep_pt,
  binning=[100,0,50],
  addOverFlowBin='upper',
))

plots.append(Plot(
  name = 'eta', texX = '#eta', texY = 'Number of Events',
  attribute = lambda event, sample: event.lep_eta,
  binning=[30,-3,3],
  addOverFlowBin='upper',
))

#tight-loose plots 
plots.append(Plot(
  name = 'LT_mu', texX = 'LT_mu', texY = 'Number of Events',
  attribute = lambda event, sample: event.nmu_mvaTOPT == 1,
  binning=[2,0,2],
  addOverFlowBin='upper',
))

plots.append(Plot(
  name = 'LT_ele', texX = 'LT_ele', texY = 'Number of Events',
  attribute = lambda event, sample: event.nele_mvaTOPT == 1,
  binning=[2,0,2],
  addOverFlowBin='upper',
))

plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events=max_events)

#drawPlots(plots)

#get histograms 
for plot in plots:
    print plot.name
    if plot.name == 'mT' :
        data_histo    = plot.histos[1][0]
        QCD_histo     = plot.histos[0][0].Clone()
        EWK_histos    = [plot.histos[0][pos].Clone()  for pos in range(len(mc)) if pos!=0]
        EWK_histo     = EWK_histos[0]
        EWK_histo.Add(EWK_histos[1])

EWK_weighted = EWK_histo.Clone()
QCD_weighted = QCD_histo.Clone()
#go to EWK dominated region to make some first scaling in order to make the template fit work
#if args.small : 
#   i_low = data_histo.GetXaxis().FindBin(80)
#   i_up  = data_histo.GetXaxis().FindBin(100)
#   
#   I_data = data_histo.Integral(i_low,i_up)
#   I_QCD = QCD_histo.Integral(i_low,i_up)
#   I_EWK = EWK_histo.Integral(i_low,i_up)
#   
#   I_scale = I_data / (I_QCD + I_EWK)
#   QCD_histo.Scale(I_scale)
#   EWK_histo.Scale(I_scale)

#tfractionfitting
tarray = ROOT.TObjArray(2)
tarray.Add( QCD_histo )
tarray.Add( EWK_histo )

data_mc_ratio = data_histo.Integral()/( QCD_histo.Integral() + EWK_histo.Integral())

fit = ROOT.TFractionFitter( data_histo, tarray ) #initialise

#fit.Constrain(0,0.0,1.0)
#fit.Constrain(1,0.0,1.0)

fit.SetRangeX(0,40)                               #random range; #fit->SetRangeX(1,15);  // use only the first 15 bins in the fit

print 'fittig'
status = fit.Fit()                                # perform the fit 
print status
if (int(status) != 0) :
      logger.info("The template fit failed - aborting")
      exit(0)

mTfit_histo = fit.GetPlot()
#QCD_histopostfit      = fit.GetMCPrediction(0)
#EWK_histopostfit      = fit.GetMCPrediction(1)

#getFitResults 
print 'FIT RESULTS XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx'
qcd   = ROOT.Double()
qcd_e = ROOT.Double()
fit.GetResult(0, qcd, qcd_e)
ewk   = ROOT.Double()
ewk_e = ROOT.Double()
fit.GetResult(1, ewk, ewk_e)

print qcd 
print qcd_e
print ewk
print ewk_e


#format histos 
data_histo.legendText = "data"
QCD_histo.legendText = "QCD"
EWK_histo.legendText = "EWK"
mTfit_histo.legendText = "mT_fit"
QCD_histo.SetLineColor(618)
QCD_histo.SetFillColor(0)
EWK_histo.SetLineColor(420)
EWK_histo.SetFillColor(0)
mTfit_histo.SetLineColor(862)
mTfit_histo.SetMarkerStyle(0)

EWK_weighted.legendText = "EWK_weighted"
QCD_weighted.legendText = "QCD_weighted"
EWK_weighted.SetLineColor(7)
EWK_weighted.SetFillColor(0)
QCD_weighted.SetLineColor(7)
QCD_weighted.SetFillColor(0)

for p in plots: 
    p.name += "_postFit"
    p.histos[0][0].Scale(qcd*data_histo.Integral()/p.histos[0][0].Integral())
    mc_ewk = p.histos[0][1].Integral()+p.histos[0][2].Integral()
    p.histos[0][1].Scale(ewk*data_histo.Integral()/mc_ewk)
    p.histos[0][2].Scale(ewk*data_histo.Integral()/mc_ewk)

drawPlots(plots)

#EWK_weighted.Scale(ewk)
#QCD_weighted.Scale(qcd)

histos = [[data_histo], [mTfit_histo]]
fitresultplot = Plot.fromHisto(
  name = 'fittemplate',
  histos = histos,
  texX = 'mT_templatefit',
  texY = 'events',
)

ewkhistos = [[EWK_histo],[EWK_weighted]]
ewkhisto = Plot.fromHisto(
  name = 'EWK',
  histos = ewkhistos,
  texX = 'EWK_templatefit',
  texY = 'events',
)

ewkoriginal = [[EWK_histo]]
EWKoriginal = Plot.fromHisto(
  name = 'EWKoriginal',
  histos = ewkoriginal,
  texX = 'EWK_original',
  texY = 'events',
)

qcdhistos = [[QCD_histo],[QCD_weighted]]
qcdhisto = Plot.fromHisto(
  name = 'QCD',
  histos = qcdhistos,
  texX = 'QCD_templatefit',
  texY = 'events',
)

plot_directory_ = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, selectionName, args.mode, 'log')
plotting.draw(fitresultplot,
              plot_directory = plot_directory_,
              logX    = False,
              logY    = True,
              sorting = False,
              #ratio   = None,
              ratio   = {},
              scaling =  {i+1:0 for i in range(len(histos)-1)},
)
plotting.draw(ewkhisto,
              plot_directory = plot_directory_,
              logX    = False,
              logY    = True,
              sorting = False,
              ratio   = None,
              #ratio   = {},
              #scaling =  {i+1:0 for i in range(len(histos)-1)},
)

plotting.draw(qcdhisto,
              plot_directory = plot_directory_,
              logX    = False,
              logY    = True,
              sorting = False,
              ratio   = None,
              #ratio   = {},
              #scaling =  {i+1:0 for i in range(len(histos)-1)},
)

#plotting.draw(EWKoriginal,
#              plot_directory = plot_directory_,
#              logX    = False,
#              logY    = True,
#              sorting = False,
#              ratio   = None,
#              #ratio   = {},
#              #scaling =  {i+1:0 for i in range(len(histos)-1)},
#)

