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

from array import *

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
argParser.add_argument('--plot_directory', action='store', default='tWZ_fakes_v3') # -v2')
argParser.add_argument('--era',            action='store', type=str, default="Run2016")
argParser.add_argument('--mode',           action='store', type=str, default="mu", choices=["mu","ele"])
argParser.add_argument('--mT',             action='store', type=int, default=-1)
argParser.add_argument('--met',            action='store', type=int, default=-1)
argParser.add_argument('--pt',             action='store', nargs='+', type=str, default=["-1", "-1"])
argParser.add_argument('--isoDiscriminator',            action='store', choices=["looseHybridIso", "FOmvaTOPT"],default="looseHybridIso")

args = argParser.parse_args()

if (args.isoDiscriminator == "looseHybridIso") :
  iso_disc = "hybridIso"
else :
  iso_disc = "mvaTOPT"


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


triggerSelection = '('+"||".join(triggers)+')'
leptonSelection  = 'n{}_{}==1'.format(args.mode,args.isoDiscriminator)
jetSelection     = 'Sum$(Jet_pt>40&&abs(Jet_eta)<2.4&&JetGood_cleaned_{}_{})>=1'.format(args.mode, args.isoDiscriminator)

max_events = -1
if args.small:
    for sample in [data_sample] + mc:
        sample.normalization = 1.
        # sample.reduceFiles( factor = 4 )
        sample.reduceFiles( to=4 )
        #sample.scale /= sample.normalization
        max_events = 30000

# Text on the plots
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)
tex.SetTextAlign(11) # align right


# fire up the cache
cache_dir_ = os.path.join(cache_dir, 'fake_cache')
dirDB      = DirDB(cache_dir_)

pu_key = ( triggerSelection, leptonSelection, jetSelection, args.era, args.small, "pu")
if dirDB.contains( pu_key ) and not args.overwrite:
    reweight_histo = dirDB.get( pu_key )
    logger.info( "Found PU reweight in cache %s", cache_dir_ )
else:
    logger.info( "Didn't find PU reweight histo %r. Obtaining it now.", pu_key)
    
    
    data_nvtx_histo = data_sample.get1DHistoFromDraw( "PV_npvsGood", [100/5, 0, 100], selectionString=data_preselectionString, weightString = "weight" )
    data_nvtx_histo.Scale(1./data_nvtx_histo.Integral())

    mc_histos  = [ s.get1DHistoFromDraw( "PV_npvsGood", [100/5, 0, 100], selectionString=mc_preselectionString, weightString = "weight*reweightBTag_SF") for s in mc]
    mc_histo     = mc_histos[0]
    for h in mc_histos[1:]:
        mc_histo.Add( h )

    mc_histo.Scale(1./mc_histo.Integral())

    reweight_histo = data_nvtx_histo.Clone()
    reweight_histo.Divide(mc_histo)
    
    dirDB.add( pu_key, reweight_histo ) 
    logger.info( "Added PU reweight to cache %s", cache_dir_ )

# define reweight
def nvtx_puRW( event, sample ):
    return reweight_histo.GetBinContent(reweight_histo.FindBin( event.PV_npvsGood ))




templateFit_QCD_key = ( triggerSelection, leptonSelection, jetSelection, args.era, args.small, "template_QCD")
templateFit_EWK_key = ( triggerSelection, leptonSelection, jetSelection, args.era, args.small, "template_EWK")
if dirDB.contains( templateFit_EWK_key ) and not args.overwrite :
    logger.info( "Found template fit in cache %s", cache_dir_ )
    EWK_sf = dirDB.get(templateFit_EWK_key)
    QCD_sf = dirDB.get(templateFit_QCD_key)
else :
    logger.info( "Didn't find template fit reweight histo %r. Obtaining it now.", templateFit_EWK_key)

    # ffile = ROOT.TFile("mT_mu_met5To200_smallFalse.root")
    ffile = ROOT.TFile("mT_mu_inclusive_smallFalse.root")

    for k in ffile.GetListOfKeys() :
      if "mT_data" in k.GetName() :
        data_mT_histo = ffile.Get(k.GetName())
      elif "mT_QCD" in k.GetName() :
        QCD_mT_histo = ffile.Get(k.GetName())
      elif "mT_WJets" in k.GetName() :
        WJ_mT_histo = ffile.Get(k.GetName())
      elif "mT_TTbar" in k.GetName() :
        TT_mT_histo = ffile.Get(k.GetName())
    

    EWK_mT_histo = WJ_mT_histo.Clone("EWK")
    EWK_mT_histo.Add(TT_mT_histo)

    # go to EWK dominated region to make some first scaling in order to make the template fit work
    i_low = data_mT_histo.GetXaxis().FindBin(80)
    i_up  = data_mT_histo.GetXaxis().FindBin(100)

    I_data = data_mT_histo.Integral(i_low,i_up)
    I_QCD = QCD_mT_histo.Integral(i_low,i_up)
    I_EWK = EWK_mT_histo.Integral(i_low,i_up)

    I_scale = I_data / (I_QCD + I_EWK)
    QCD_mT_histo.Scale(I_scale)
    EWK_mT_histo.Scale(I_scale)
    

    vcolor = [2,3]

    tarray = ROOT.TObjArray(len(vcolor))
    tarray.Add(QCD_mT_histo)
    tarray.Add(EWK_mT_histo)

    templateFit = ROOT.TFractionFitter(data_mT_histo,tarray)
    ROOT.SetOwnership( templateFit, False ) 
    for i in xrange(len(vcolor)) :
      templateFit.Constrain(i,0.0,1.0) # each mc template is allowed to be only between [0,1]

    #could exclude bins from fit range
    templateFit.SetRangeX(5,30)

    status = templateFit.Fit()
    
    if (int(status) != 0) :
      logger.info("The template fit failed - aborting")
      exit(0)

    template_hist = templateFit.GetPlot()
    EWK_sf = templateFit.GetMCPrediction(1)
    QCD_sf = templateFit.GetMCPrediction(0)
    

    c0 = ROOT.TCanvas()
    ROOT.gPad.SetLogy()
    QCD_mT_histo.SetLineColor(2)
    EWK_mT_histo.SetLineColor(3)
    template_hist.SetLineColor(4)
    EWK_sf.SetLineColor(6)
    QCD_mT_histo.SetFillStyle(3001)
    EWK_mT_histo.SetFillStyle(3001)
    template_hist.SetFillStyle(3001)
    EWK_sf.SetFillStyle(3001)
    
    
    QCD_mT_histo.Draw("HIST")
    EWK_mT_histo.Draw("HIST same")
    template_hist.Draw("HIST same")
    EWK_sf.Draw("HIST same")
    data_mT_histo.Draw("EP same")
    
    EWK_sf.Draw("EP same")
    c0.SaveAs("temp.png")

    c1 = ROOT.TCanvas()
    ROOT.gPad.SetLogy()
    EWK_sf.Divide(EWK_mT_histo)
    EWK_sf.Scale(I_scale)
    EWK_sf.Draw("HIST E same")
    c1.SaveAs("SF_EWK.png")

    c2 = ROOT.TCanvas()
    ROOT.gPad.SetLogy()
    QCD_sf.Divide(QCD_mT_histo)
    QCD_sf.Scale(I_scale)
    QCD_sf.Draw("HIST E same")
    c2.SaveAs("SF_QCD.png")

    
    dirDB.add( templateFit_EWK_key, EWK_sf ) 
    dirDB.add( templateFit_QCD_key, QCD_sf ) 
    logger.info( "Added template fit to cache %s", cache_dir_ )
    

def qcd_prescale_x_nvtx_puRW( event, sample ):
    return QCD_sf.GetBinContent(QCD_sf.FindBin( event.lep_mT )) * nvtx_puRW(event, sample)
def ewk_prescale_x_nvtx_puRW( event, sample ):
    return EWK_sf.GetBinContent(EWK_sf.FindBin( event.lep_mT )) * nvtx_puRW(event, sample)

    
    

selectionName    = "inclusive"
eventSelection   = "(1)"
if args.met>0:   
    selectionName = "metTo{met}".format(met=args.met)
    eventSelection= "met_pt<={met}".format(met=args.met)
if args.mT>0:
    mT_str = "Sum$({mode}_{isoDiscriminator}_mT<{mT})==1".format(**args.__dict__) 
    if args.met>0:
        eventSelection+= "&&"+mT_str
        selectionName += "-mTTo{mT}".format(mT=args.mT)
    else: 
        eventSelection = mT_str 
        selectionName  = "mTTo{mT}".format(mT=args.mT)
if args.pt != ["-1","-1"] :
  pt_lower = ""
  pt_upper = ""
  if args.pt[0] != "-1" :
    pt_lower = args.pt[0]
  if args.pt[1] != "-1" :
    pt_upper = args.pt[1]

  pt_str = "Sum$({mode}_{isoDiscriminator}_pt<{upper}) + Sum$({mode}_{isoDiscriminator}_pt>={lower})==2".format(mode=args.mode,isoDiscriminator=args.isoDiscriminator,upper=pt_upper,lower=pt_lower) 
    
  if eventSelection != "(1)":
      eventSelection+= "&&"+pt_str
      selectionName += "-pt{}To{}".format(pt_lower,pt_upper)
  else: 
      eventSelection = pt_str 
      selectionName  = "pt{}To{}".format(pt_lower,pt_upper)



logger.info("selection %s -> %s", selectionName, eventSelection)

data_selectionString = "&&".join([getFilterCut(isData=True, year=year), triggerSelection, leptonSelection, jetSelection, eventSelection])
data_sample.setSelectionString( data_selectionString )
mc_selectionString   = "&&".join([getFilterCut(isData=False, year=year), triggerSelection, leptonSelection, jetSelection, eventSelection])
for s in mc:
    s.setSelectionString( mc_selectionString )

#lumi_scale                 = data_sample.lumi/1000
data_sample.scale   = 1.
mc[0].weight = qcd_prescale_x_nvtx_puRW 
for sample in mc[1:]:
   sample.weight   = ewk_prescale_x_nvtx_puRW


def drawObjects():
    lines = [
      (0.15, 0.95, 'CMS Preliminary'), 
      #(0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines] 

def drawPlots(plots):
  for log in [True] : #[False, True]:
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
            # scaling = {0:1},
            legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
            drawObjects = drawObjects() + _drawObjects,
            # copyIndexPHP = True, 
            copyIndexPHP = False, 
            extensions = ["png"],
          )
            
# Read variables and sequences


read_variables = [
    "weight/F", "PV_npvsGood/I",
    "JetGood[pt/F,phi/F,eta/F]",
#    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
#    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
    ]

sequence       = []

read_variables += ["n{}_{}/I".format(args.mode,args.isoDiscriminator), "{}_{}[pt/F,eta/F,phi/F,mT/F]".format(args.mode,args.isoDiscriminator), "met_pt/F", "nmu_{}/I".format(iso_disc), "nele_{}/I".format(iso_disc)]

def makeLeptons( event, sample ):
    collVars = ["pt","eta","phi","mT",iso_disc]
    lep  = getObjDict(event, args.mode+'_{}_'.format(args.isoDiscriminator), collVars, 0)
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
  name = 'pt', texX = 'p_{T}', texY = 'Number of Events',
  attribute = lambda event, sample: event.lep_pt,
  binning=[100,0,50],
  addOverFlowBin='upper',
))

plots.append(Plot(
  name = 'pt_rebin', texX = 'p_{T}', texY = 'Number of Events',
  attribute = lambda event, sample: event.lep_pt,
  binning=Binning.fromThresholds([3.5,5,8,10,15,20,30,50]),
  addOverFlowBin='upper',
))

plots.append(Plot(
  name = 'eta', texX = '#eta', texY = 'Number of Events',
  attribute = lambda event, sample: event.lep_eta,
  binning=[30,-3,3],
  addOverFlowBin='upper',
))

plots.append(Plot(
  name = 'abseta', texX = '#eta', texY = 'Number of Events',
  attribute = lambda event, sample: abs(event.lep_eta),
  binning=[15,0,3],
  addOverFlowBin='upper',
))

plots.append(Plot(
  name = 'mT', texX = 'm_{T}', texY = 'Number of Events',
  attribute = lambda event, sample: event.lep_mT,
  binning=[40,0,200],
  addOverFlowBin='upper',
))

if args.mode == "mu" :
  plots.append(Plot(
    name = 'LT_mu', texX = 'LT_mu', texY = 'Number of Events',
    attribute = lambda event, sample: event.nmu_hybridIso == 1,
    binning=[2,0,2],
    addOverFlowBin='upper',
  ))
else :
  plots.append(Plot(
    name = 'LT_ele', texX = 'LT_ele', texY = 'Number of Events',
    attribute = lambda event, sample: event.nele_hybridIso == 1,
    binning=[2,0,2],
    addOverFlowBin='upper',
  ))

# plots.append(Plot(
#  name = 'LT', texX = 'LT_mu', texY = 'Number of Events',
#  attribute = lambda event, sample: event.lep_hybridIso,
#  binning=[2,0,2],
#  addOverFlowBin='upper',
# ))


#plots.append(Plot(
#  name = 'dR_jet0', texX = 'm_{T}', texY = 'Number of Events',
#  attribute = lambda event, sample: cos(event.lep_phi - event.JetGood_phi[0] ),
#  binning=[40,-1,1],
#  addOverFlowBin='upper',
#))




# plots = []
plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events=max_events)

drawPlots(plots)