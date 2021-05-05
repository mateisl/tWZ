#!/usr/bin/env python
''' analysis script for standard plots with systematic errors
'''

# Standard imports and batch mode
import ROOT
ROOT.gROOT.SetBatch(True)
import operator
import pickle, os, time, sys
from math                                import sqrt, cos, sin, pi, atan2

# RootTools
from RootTools.core.standard             import *

#Analysis / StopsDilepton / Samples
from StopsDilepton.tools.user            import plot_directory
from StopsDilepton.tools.helpers         import deltaPhi, add_histos
from Analysis.Tools.metFilters            import getFilterCut
from StopsDilepton.tools.cutInterpreter  import cutInterpreter
from StopsDilepton.tools.RecoilCorrector import RecoilCorrector
from StopsDilepton.tools.mt2Calculator   import mt2Calculator
from Analysis.Tools.puReweighting        import getReweightingFunction
from Analysis.Tools.DirDB                import DirDB


# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',          action='store',      default='INFO',     nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--signal',            action='store',      default=None,        nargs='?', choices=['None', "T2tt",'DM'], help="Add signal to plot")
argParser.add_argument('--plot_directory',    action='store',      default='v1')
argParser.add_argument('--selection',         action='store',            default='njet2p-btag1p-relIso0.12-looseLeptonVeto-mll20-met80-metSig5-dPhiJet0-dPhiJet1')
argParser.add_argument('--variation',         action='store',      default=None, help="Which systematic variation to run. Don't specify for producing plots.")
argParser.add_argument('--small',             action='store_true',     help='Run only on a small subset of the data?')
argParser.add_argument('--scaling',           action='store',      default=None, choices = [None, 'mc', 'top'],     help='Scale top to data in mt2ll<100?')
argParser.add_argument('--variation_scaling', action='store_true', help='Scale the variations individually to mimick bkg estimation?')
argParser.add_argument('--overwrite',         action='store_true',     help='Overwrite?')
argParser.add_argument('--mode',              action='store',      default = 'all', choices = ['mumu', 'ee', 'mue', 'all'],   help='Which mode?')
argParser.add_argument('--normalizeBinWidth', action='store_true', default=False,       help='normalize wider bins?')
argParser.add_argument('--reweightPU',         action='store', default='Central', choices=[ 'Central', 'VUp'] )
#argParser.add_argument('--recoil',             action='store', type=str,      default="Central", choices = ["nvtx", "VUp", "Central"])
argParser.add_argument('--era',                action='store', type=str,      default="Run2016")

args = argParser.parse_args()

# Logger
import StopsDilepton.tools.logger as logger
import RootTools.core.logger as logger_rt
import Analysis.Tools.logger as logger_an
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)
logger_an = logger_an.get_logger(args.logLevel, logFile = None)

# Year
if "2016" in args.era:
    year = 2016
elif "2017" in args.era:
    year = 2017
elif "2018" in args.era:
    year = 2018

logger.info( "Working in year %i", year )

if args.reweightPU == 'Central':
    nominalPuWeight, upPUWeight, downPUWeight = "reweightPU", "reweightPUUp", "reweightPUDown"
elif args.reweightPU == 'VUp':
    nominalPuWeight, upPUWeight, downPUWeight = "reweightPUVUp", "reweightPUVVUp", "reweightPUUp"

def jetSelectionModifier( sys, returntype = "func"):
    #Need to make sure all jet variations of the following observables are in the ntuple
    variiedJetObservables = ['nJetGood', 'nBTag', 'dl_mt2ll', 'dl_mt2blbl', 'MET_significance', 'met_pt', 'metSig']
    if returntype == "func":
        def changeCut_( string ):
            for s in variiedJetObservables:
                string = string.replace(s, s+'_'+sys)
            return string
        return changeCut_
    elif returntype == "list":
        return [ v+'_'+sys for v in variiedJetObservables ]

def metSelectionModifier( sys, returntype = 'func'):
    #Need to make sure all MET variations of the following observables are in the ntuple
    variiedMetObservables = ['dl_mt2ll', 'dl_mt2blbl', 'MET_significance', 'met_pt', 'metSig']
    if returntype == "func":
        def changeCut_( string ):
            for s in variiedMetObservables:
                string = string.replace(s, s+'_'+sys)
            return string
        return changeCut_
    elif returntype == "list":
        return [ v+'_'+sys for v in variiedMetObservables ]

# these are the nominal MC weights we always apply
nominalMCWeights = ["weight", "reweightLeptonSF", "reweightPU", "reweightDilepTrigger", "reweightBTag_SF"]

# weight the MC according to a variation
def MC_WEIGHT( variation, returntype = "string"):
    variiedMCWeights = list(nominalMCWeights)   # deep copy
    if variation.has_key('replaceWeight'):
        for i_w, w in enumerate(variiedMCWeights):
            if w == variation['replaceWeight'][0]:
                variiedMCWeights[i_w] = variation['replaceWeight'][1]
                break
        # Let's make sure we don't screw it up ... because we mostly do.
        if variiedMCWeights==nominalMCWeights:
            raise RuntimeError( "Tried to change weight %s to %s but didn't find it in list %r" % ( variation['replaceWeight'][0], variation['replaceWeight'][1], variiedMCWeights ))
    # multiply strings for ROOT weights
    if returntype == "string":
        return "*".join(variiedMCWeights)
    # create a function that multiplies the attributes of the event
    elif returntype == "func":
        getters = map( operator.attrgetter, variiedMCWeights)
        def weight_( event, sample):
            return reduce(operator.mul, [g(event) for g in getters], 1)
        return weight_
    elif returntype == "list":
        return variiedMCWeights

def data_weight( event, sample ):
    return event.weight

data_weight_string = "weight"

# Define all systematic variations
variations = {
    'central'           : {'read_variables': [ '%s/F'%v for v in nominalMCWeights ]},
    'jesTotalUp'        : {'selectionModifier':jetSelectionModifier('jesTotalUp'),               'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('jesTotalUp','list')]},
    'jesTotalDown'      : {'selectionModifier':jetSelectionModifier('jesTotalDown'),             'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('jesTotalDown','list')]},
    'unclustEnUp'       : {'selectionModifier':metSelectionModifier('unclustEnUp'),              'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('unclustEnUp','list')]},
    'unclustEnDown'     : {'selectionModifier':metSelectionModifier('unclustEnDown'),            'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('unclustEnDown','list')]},
    'PUUp'              : {'replaceWeight':(nominalPuWeight,upPUWeight),                         'read_variables' : [ '%s/F'%v for v in nominalMCWeights + [upPUWeight] ]},
    'PUDown'            : {'replaceWeight':(nominalPuWeight,downPUWeight),                       'read_variables' : [ '%s/F'%v for v in nominalMCWeights + [downPUWeight] ]},
    'BTag_SF_b_Down'    : {'replaceWeight':('reweightBTag_SF','reweightBTag_SF_b_Down'),         'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightBTag_SF_b_Down']]},  
    'BTag_SF_b_Up'      : {'replaceWeight':('reweightBTag_SF','reweightBTag_SF_b_Up'),           'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightBTag_SF_b_Up'] ]},
    'BTag_SF_l_Down'    : {'replaceWeight':('reweightBTag_SF','reweightBTag_SF_l_Down'),         'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightBTag_SF_l_Down']]},
    'BTag_SF_l_Up'      : {'replaceWeight':('reweightBTag_SF','reweightBTag_SF_l_Up'),           'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightBTag_SF_l_Up'] ]},
    'DilepTriggerDown'  : {'replaceWeight':('reweightDilepTrigger','reweightDilepTriggerDown'),  'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightDilepTriggerDown']]},
    'DilepTriggerUp'    : {'replaceWeight':('reweightDilepTrigger','reweightDilepTriggerUp'),    'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightDilepTriggerUp']]},
    'LeptonSFDown'      : {'replaceWeight':('reweightLeptonSF','reweightLeptonSFDown'),          'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightLeptonSFDown']]},
    'LeptonSFUp'        : {'replaceWeight':('reweightLeptonSF','reweightLeptonSFUp'),            'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightLeptonSFUp']]},
#    'TopPt':{},
#   'JERUp':{},
#   'JERDown':{},
}

# Add a default selection modifier that does nothing
for key, variation in variations.iteritems():
    if not variation.has_key('selectionModifier'):
        variation['selectionModifier'] = lambda string:string
    if not variation.has_key('read_variables'):
        variation['read_variables'] = [] 

# Check if we know the variation
if args.variation is not None and args.variation not in variations.keys():
    raise RuntimeError( "Variation %s not among the known: %s", args.variation, ",".join( variation.keys() ) )

# arguments & directory
plot_subdirectory = args.plot_directory
if args.signal == "DM":           plot_subdirectory += "_DM"
if args.signal == "T2tt":         plot_subdirectory += "_T2tt"
if args.small:                    plot_subdirectory += "_small"
if args.reweightPU:               plot_subdirectory += "_reweightPU%s"%args.reweightPU
#if args.recoil:                  plot_subdirectory  += '_recoil_'+args.recoil
    
if year == 2016:
    from StopsDilepton.samples.nanoTuples_Summer16_postProcessed import *
    from StopsDilepton.samples.nanoTuples_Run2016_17Jul2018_postProcessed import *
    Top_pow, TTXNoZ, TTZ_LO, multiBoson, DY_HT_LO = Top_pow_16, TTXNoZ_16, TTZ_16, multiBoson_16, DY_HT_LO_16
    mc              = [ Top_pow_16, TTXNoZ_16, TTZ_16, multiBoson_16, DY_HT_LO_16]
elif year == 2017:
    from StopsDilepton.samples.nanoTuples_Fall17_postProcessed import *
    from StopsDilepton.samples.nanoTuples_Run2017_31Mar2018_postProcessed import *
    Top_pow, TTXNoZ, TTZ_LO, multiBoson, DY_HT_LO = Top_pow_17, TTXNoZ_17, TTZ_17, multiBoson_17, DY_LO_17
    mc              = [ Top_pow_17, TTXNoZ_17, TTZ_17, multiBoson_17, DY_LO_17]
elif year == 2018:
    from StopsDilepton.samples.nanoTuples_Run2018_PromptReco_postProcessed import *
    from StopsDilepton.samples.nanoTuples_Autumn18_postProcessed import *
    Top_pow, TTXNoZ, TTZ_LO, multiBoson, DY_HT_LO = Top_pow_18, TTXNoZ_18, TTZ_18, multiBoson_18, DY_LO_18
    mc              = [ Top_pow_18, TTXNoZ_18, TTZ_18, multiBoson_18, DY_LO_18]

try:
  data_sample = eval(args.era)
except Exception as e:
  logger.error( "Didn't find %s", args.era )
  raise e

# postions of MC components in list
position = {s.name:i_s for i_s,s in enumerate(mc)}

#if args.recoil:
#    from Analysis.Tools.RecoilCorrector import RecoilCorrector
#    if args.recoil == "nvtx":
#        recoilCorrector = RecoilCorrector( os.path.join( "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/results/", "recoil_v4.3_fine_nvtx_loop", "%s_lepSel-njet1p-btag0-relIso0.12-looseLeptonVeto-mll20-onZ_recoil_fitResults_SF.pkl"%args.era ) )
#    elif args.recoil == "VUp":
#        recoilCorrector = RecoilCorrector( os.path.join( "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/results/", "recoil_v4.3_fine_VUp_loop", "%s_lepSel-njet1p-btag0-relIso0.12-looseLeptonVeto-mll20-onZ_recoil_fitResults_SF.pkl"%args.era ) )
#    elif args.recoil is "Central":
#        recoilCorrector = RecoilCorrector( os.path.join( "/afs/hephy.at/data/rschoefbeck01/StopsDilepton/results/", "recoil_v4.3_fine", "%s_lepSel-njet1p-btag0-relIso0.12-looseLeptonVeto-mll20-onZ_recoil_fitResults_SF.pkl"%args.era ) )
#

# Read variables and sequences
read_variables = ["weight/F", "l1_pt/F", "l2_pt/F", "l1_eta/F" , "l1_phi/F", "l2_eta/F", "l2_phi/F", "JetGood[pt/F,eta/F,phi/F]", "dl_mass/F", "dl_eta/F", "dl_mt2ll/F", "dl_mt2bb/F", "dl_mt2blbl/F",
                  "met_pt/F", "met_phi/F", "dl_pt/F", "dl_phi/F",
                  #"LepGood[pt/F,eta/F,miniRelIso/F]", "nGoodMuons/F", "nGoodElectrons/F", "l1_mIsoWP/F", "l2_mIsoWP/F",
                  "metSig/F", "ht/F", "nBTag/I", "nJetGood/I","run/I","event/l"]

sequence = []
#def corr_recoil( event, sample ):
#    mt2Calculator.reset()
#    if not sample.isData: 
#        # Parametrisation vector - # define qt as GenMET + leptons
#        qt_px = event.l1_pt*cos(event.l1_phi) + event.l2_pt*cos(event.l2_phi) + event.GenMET_pt*cos(event.GenMET_phi)
#        qt_py = event.l1_pt*sin(event.l1_phi) + event.l2_pt*sin(event.l2_phi) + event.GenMET_pt*sin(event.GenMET_phi)
#
#        qt = sqrt( qt_px**2 + qt_py**2 )
#        qt_phi = atan2( qt_py, qt_px )
#
#        #ref_phi = qt_phi
#        ref_phi = event.dl_phi
#
#        # compute fake MET 
#        fakeMET_x = event.met_pt*cos(event.met_phi) - event.GenMET_pt*cos(event.GenMET_phi)
#        fakeMET_y = event.met_pt*sin(event.met_phi) - event.GenMET_pt*sin(event.GenMET_phi)
#
#        fakeMET = sqrt( fakeMET_x**2 + fakeMET_y**2 )
#        fakeMET_phi = atan2( fakeMET_y, fakeMET_x )
#
#        # project fake MET on qT
#        fakeMET_para = fakeMET*cos( fakeMET_phi - ref_phi )
#        fakeMET_perp = fakeMET*cos( fakeMET_phi - ( ref_phi - pi/2) )
#
#        fakeMET_para_corr = - recoilCorrector.predict_para( ref_phi, qt, -fakeMET_para )
#        fakeMET_perp_corr = - recoilCorrector.predict_perp( ref_phi, qt, -fakeMET_perp )
#
#        # rebuild fake MET vector
#        fakeMET_px_corr = fakeMET_para_corr*cos(ref_phi) + fakeMET_perp_corr*cos(ref_phi - pi/2)
#        fakeMET_py_corr = fakeMET_para_corr*sin(ref_phi) + fakeMET_perp_corr*sin(ref_phi - pi/2)
#        #print "%s qt: %3.2f para %3.2f->%3.2f perp %3.2f->%3.2f fakeMET(%3.2f,%3.2f) -> (%3.2f,%3.2f)" % ( sample.name, qt, fakeMET_para, fakeMET_para_corr, fakeMET_perp, fakeMET_perp_corr, fakeMET, fakeMET_phi, sqrt( fakeMET_px_corr**2+fakeMET_py_corr**2), atan2( fakeMET_py_corr, fakeMET_px_corr) )
#   
#        for var in [""] + jme_systematics:
#            if var: var = "_"+var
#            met_px_corr = getattr(event, "met_pt"+var)*cos(getattr(event, "met_phi"+var)) - fakeMET_x + fakeMET_px_corr 
#            met_py_corr = getattr(event, "met_pt"+var)*sin(getattr(event, "met_phi"+var)) - fakeMET_y + fakeMET_py_corr
#    
#            setattr(event, "met_pt_corr"+var, sqrt( met_px_corr**2 + met_py_corr**2 ) )
#            setattr(event, "met_phi_corr"+var, atan2( met_py_corr, met_px_corr ) )
#            
#            mt2Calculator.setLeptons(event.l1_pt, event.l1_eta, event.l1_phi, event.l2_pt, event.l2_eta, event.l2_phi)
#            mt2Calculator.setMet(getattr(event,"met_pt_corr"+var), getattr(event,"met_phi_corr"+var))
#            setattr(event, "dl_mt2ll_corr"+var, mt2Calculator.mt2ll() )
#
#    else:
#        event.met_pt_corr  = event.met_pt 
#        event.met_phi_corr = event.met_phi
#
#        mt2Calculator.setLeptons(event.l1_pt, event.l1_eta, event.l1_phi, event.l2_pt, event.l2_eta, event.l2_phi)
#        mt2Calculator.setMet(event.met_pt_corr, event.met_phi_corr)
#        event.dl_mt2ll_corr =  mt2Calculator.mt2ll()
#
#    #print event.dl_mt2ll, event.dl_mt2ll_corr
#
#sequence.append( corr_recoil )

signals = []

# selection
offZ = "&&abs(dl_mass-91.1876)>15" if not (args.selection.count("onZ") or args.selection.count("allZ") or args.selection.count("offZ")) else ""
def getLeptonSelection( mode ):
  if   mode=="mumu": return "nGoodMuons==2&&nGoodElectrons==0&&isOS&&isMuMu" + offZ
  elif mode=="mue":  return "nGoodMuons==1&&nGoodElectrons==1&&isOS&&isEMu"
  elif mode=="ee":   return "nGoodMuons==0&&nGoodElectrons==2&&isOS&&isEE" + offZ

modes = ['mumu', 'mue', 'ee'] if args.mode=='all' else [ args.mode ]

allPlots   = {}

logger.info('Working on modes: %s', ','.join(modes))

if year == 2016:
    data_sample = Run2016
    data_sample.texName = "data (2016)"
elif year == 2017:
    data_sample = Run2017
    data_sample.texName = "data (2017)"
elif year == 2018:
    data_sample = Run2018
    data_sample.texName = "data (2018)"

# Define samples
data_sample.name           = "data"
data_sample.read_variables = ["event/I","run/I"]
data_sample.style          = styles.errorStyle(ROOT.kBlack)
data_sample.scale          = 1.
lumi_scale                 = data_sample.lumi/1000
logger.info('Lumi scale is ' + str(lumi_scale))
for sample in mc:
    sample.scale           = lumi_scale
    sample.style           = styles.fillStyle(sample.color, lineColor = sample.color)
    sample.read_variables  = ['Pileup_nTrueInt/F', 'GenMET_pt/F', 'GenMET_phi/F']
    # append variables for systematics
    if args.variation is not None:
        sample.read_variables+=list(set(variations[args.variation]['read_variables']))

# reduce if small
if args.small:
  data_sample.normalization = 1.
  data_sample.reduceFiles( to = 1 )
  data_sample.scale /= data_sample.normalization
  for sample in mc:
    sample.normalization = 1.
    #sample.reduceFiles( factor = 40 )
    sample.reduceFiles( to = 1 )
    sample.scale /= sample.normalization

# Fire up the cache
dirDB = DirDB(os.path.join(plot_directory, 'systematicPlots', plot_subdirectory, args.selection, 'cache'))

# loop over modes
for mode in modes:
    logger.info('Working on mode: %s', mode)

    # set selection
    data_sample.setSelectionString([getFilterCut(isData=True, year=year), getLeptonSelection(mode)])
    for sample in mc:
        sample.setSelectionString([getFilterCut(isData=False, year=year), getLeptonSelection(mode)])

    # Use some defaults
    Plot.setDefaults( selectionString = cutInterpreter.cutString(args.selection) )

    # if we're running a variation specify
    if args.variation is not None:
        selectionModifier = variations[args.variation]['selectionModifier']
        mc_weight         = MC_WEIGHT( variation = variations[args.variation], returntype='func')
    else:
        selectionModifier = None 
        mc_weight         = None 

    # Stack
    stack_mc   = Stack( mc )
    stack_data = Stack( data_sample )

    plots      = []

    mt2llBinning = [0,20,40,60,80,100,140,240,340]
    if args.variation == 'central':
        dl_mt2ll_data   = Plot(
            name        = "dl_mt2ll_data",
            texX        = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
            binning     = Binning.fromThresholds(mt2llBinning),
            stack       = stack_data,
            attribute   = TreeVariable.fromString( "dl_mt2ll/F" ),
            weight      = data_weight )
        plots.append( dl_mt2ll_data )

    dl_mt2ll_mc  = Plot(\
        name            = "dl_mt2ll_mc",
        texX            = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
        binning         = Binning.fromThresholds(mt2llBinning),
        stack           = stack_mc,
        attribute       = TreeVariable.fromString( selectionModifier("dl_mt2ll/F"))   if selectionModifier is not None else None,
        selectionString = selectionModifier(cutInterpreter.cutString(args.selection)) if selectionModifier is not None else None,
        weight          = mc_weight )
    plots.append( dl_mt2ll_mc )

    #  #mt2llCorrBinning = [0,20,40,60,80,100,140,240,340]
    #  #
    #  #dl_mt2llCorr_data  = Plot(
    #  #    name = "dl_mt2ll_corr_data",
    #  #    texX = 'M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
    #  #    binning=Binning.fromThresholds(mt2llBinning),
    #  #    stack = stack_data,
    #  #    attribute = TreeVariable.fromString( "dl_mt2ll/F" ),
    ## #?     attribute = lambda event, sample: event.dl_mt2ll_corr,
    #  #    weight = data_weight,
    #  #    )
    #  #plots.append( dl_mt2llCorr_data )
    #
    #  #dl_mt2llCorr_mc  = { sys:Plot(\
    #  #    name            = "dl_mt2ll_corr" if sys is None else "dl_mt2ll_corr_mc_%s" % sys,
    #  #    texX            = 'corrected M_{T2}(ll) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
    #  #    binning         = Binning.fromThresholds(mt2llBinning),
    #  #    stack           = sys_stacks[sys],
    #  #    attribute       = lambda event, sample: event.dl_mt2ll_corr if sys is None or sys in weight_systematics else getattr(event, "dl_mt2ll_corr_"+sys) ,
    ## #     attribute        = TreeVariable.fromString( "dl_mt2ll/F" ) if sys is None or sys in weight_systematics else TreeVariable.fromString( "dl_mt2ll_%s/F" % sys ),
    #  #    selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
    #  #    weight          = MC_WEIGHT( sys = sys )[0],
    #  #    ) for sys in all_systematics }
    #  #plots.extend( dl_mt2llCorr_mc.values() )
    #
    if args.selection.count('njet2'):
        if args.variation == 'central':
            dl_mt2blbl_data  = Plot( 
                name = "dl_mt2blbl_data",
                texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
                binning=Binning.fromThresholds([0,20,40,60,80,100,120,140,160,200,250,300,350]),
                stack = stack_data,
                attribute = TreeVariable.fromString( "dl_mt2blbl/F" ),
                weight = data_weight,
                ) 
            plots.append( dl_mt2blbl_data )
    
        dl_mt2blbl_mc  = Plot(
            name = "dl_mt2blbl_mc",
            texX = 'M_{T2}(blbl) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
            binning=Binning.fromThresholds([0,20,40,60,80,100,120,140,160,200,250,300,350]),
            stack = stack_mc,
            attribute = TreeVariable.fromString( selectionModifier("dl_mt2blbl/F") )      if selectionModifier is not None else None,
            selectionString = selectionModifier(cutInterpreter.cutString(args.selection)) if selectionModifier is not None else None,
            weight          = mc_weight )
        plots.append( dl_mt2blbl_mc )
    
    #  nBtagBinning = [5, 1, 6] if args.selection.count('btag1p') else [1,0,1]
    #
    #  nbtags_data  = Plot( 
    #      name = "nbtags_data",
    #      texX = 'number of b-tags (CSVM)', texY = 'Number of Events',
    #      stack = stack_data,
    #      attribute = TreeVariable.fromString('nBTag/I'),
    #      binning=nBtagBinning,
    #      weight = data_weight,
    #      ) 
    #  plots.append( nbtags_data )
    #
    #  nbtags_mc  = {sys: Plot(
    #      name = "nbtags" if sys is None else "nbtags_mc_%s" % sys,
    #      texX = 'number of b-tags (CSVM)', texY = 'Number of Events',
    #      stack = sys_stacks[sys],
    #      attribute = TreeVariable.fromString('nBTag/I') if sys is None or sys in weight_systematics or sys in met_systematics else TreeVariable.fromString( "nBTag_%s/I" % sys ),
    #      binning=nBtagBinning,
    #      selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
    #      weight = MC_WEIGHT( sys = sys )[0],
    #      ) for sys in all_systematics }
    #  plots.extend( nbtags_mc.values() )
    #
    #  jetBinning = [8,2,10] if args.selection.count('njet2') else [2,0,2]
    #
    #  njets_data  = Plot( 
    #      name = "njets_data",
    #      texX = 'number of jets', texY = 'Number of Events',
    #      stack = stack_data,
    #      attribute = TreeVariable.fromString('nJetGood/I'),
    #      binning=jetBinning,
    #      weight = data_weight,
    #      )
    #  plots.append( njets_data )
    #
    #  njets_mc  = {sys: Plot(
    #      name = "njets" if sys is None else "njets_mc_%s" % sys,
    #      texX = 'number of jets', texY = 'Number of Events',
    #      stack = sys_stacks[sys],
    #      attribute = TreeVariable.fromString('nJetGood/I') if sys is None or sys in weight_systematics or sys in met_systematics else TreeVariable.fromString( "nJetGood_%s/I" % sys ),
    #      binning= jetBinning,
    #      selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
    #      weight = MC_WEIGHT( sys = sys )[0],
    #      ) for sys in all_systematics }
    #  plots.extend( njets_mc.values() )
    #
    #  metBinning = [0,20,40,60,80] if args.selection.count('metInv') else [80,130,180,230,280,320,420,520,800] if args.selection.count('met80') else [0,80,130,180,230,280,320,420,520,800]
    #
    #  met_data  = Plot( 
    #      name = "met_data",
    #      texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 50 GeV' if args.normalizeBinWidth else "Number of Events",
    #      stack = stack_data, 
    #      attribute = TreeVariable.fromString( "met_pt/F" ),
    #      binning=Binning.fromThresholds( metBinning ),
    #      weight = data_weight,
    #      )
    #  plots.append( met_data )
    #
    #  met_mc  = {sys: Plot(
    #      name = "met_pt" if sys is None else "met_pt_mc_%s" % sys,
    #      texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 50 GeV' if args.normalizeBinWidth else "Number of Events",
    #      stack = sys_stacks[sys],
    #      attribute = TreeVariable.fromString('met_pt/F') if sys not in met_systematics else TreeVariable.fromString( "met_pt_%s/F" % sys ),
    #      binning=Binning.fromThresholds( metBinning ),
    #      selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
    #      weight = MC_WEIGHT( sys = sys )[0],
    #      ) for sys in all_systematics }
    #  plots.extend( met_mc.values() )
    #
    #  metBinning2 = [0,20,40,60,80] if args.selection.count('metInv') else [80,100,120,140,160,200,500] if args.selection.count('met80') else [0,80,100,120,140,160,200,500]
    #
    #  met2_data  = Plot(
    #      name = "met2_data",
    #      texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
    #      stack = stack_data,
    #      attribute = TreeVariable.fromString( "met_pt/F" ),
    #      binning=Binning.fromThresholds( metBinning2 ),
    #      weight = data_weight,
    #      )
    #  plots.append( met2_data )
    #
    #  met2_mc  = {sys: Plot(
    #      name = "met2_pt" if sys is None else "met2_pt_mc_%s" % sys,
    #      texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
    #      stack = sys_stacks[sys],
    #      attribute = TreeVariable.fromString('met_pt/F') if sys not in met_systematics else TreeVariable.fromString( "met_pt_%s/F" % sys ),
    #      binning=Binning.fromThresholds( metBinning2 ),
    #      selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
    #      weight = MC_WEIGHT( sys = sys )[0],
    #      ) for sys in all_systematics }
    #  plots.extend( met2_mc.values() )
    #
    #  metSigBinning = [0,2,4,6,8,10,12] if args.selection.count('POGMetSig0To12') else [12,16,20,24,28,32,36,40,45,50,55,60,65,70,75,80,85,90,95,100] if args.selection.count('POGMetSig12') else [0,4,8,12,16,20,24,28,32,36,40,45,50,55,60,65,70,75,80,85,90,95,100]
    #
    #  metSig_data  = Plot( 
    #      name = "MET_significance_data",
    #      texX = 'E_{T}^{miss} significance (GeV)', texY = 'Number of Events / 5 GeV' if args.normalizeBinWidth else "Number of Events",
    #      stack = stack_data, 
    #      attribute = TreeVariable.fromString( "MET_significance/F" ),
    #      binning=Binning.fromThresholds( metSigBinning ),
    #      weight = data_weight,
    #      )
    #  plots.append( metSig_data )
    #
    #  metSig_mc  = {sys: Plot(
    #      name = "MET_Significance" if sys is None else "MET_significance_mc_%s" % sys,
    #      texX = 'E_{T}^{miss} significance (GeV)', texY = 'Number of Events / 5 GeV' if args.normalizeBinWidth else "Number of Events",
    #      stack = sys_stacks[sys],
    #      attribute = TreeVariable.fromString('MET_significance/F') if sys not in jme_systematics else TreeVariable.fromString( "MET_significance_%s/F" % sys ),
    #      binning=Binning.fromThresholds( metSigBinning ),
    #      selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
    #      weight = MC_WEIGHT( sys = sys )[0],
    #      ) for sys in all_systematics }
    #  plots.extend( metSig_mc.values() )

    ## Make plot directory
    #plot_directory_ = os.path.join(plot_directory, 'systematicPlots', args.plot_directory, args.selection, args.era, mode)
    #try: 
    #    os.makedirs(plot_directory_)
    #except: 
    #    pass

    if args.variation is not None:
        key  = (args.era, mode, args.variation)
        if dirDB.contains(key) and not args.overwrite:
            normalisation_mc, normalisation_data, histos = dirDB.get( key )
            for i_p, h_s in enumerate(histos):
                plots[i_p].histos = h_s
            logger.info( "Loaded normalisations and histograms for %s in mode %s from cache.", args.era, mode)
        else:
            logger.info( "Obtain normalisations and histograms for %s in mode %s.", args.era, mode)
            # Calculate the normalisation yield for mt2ll<100
            normalization_selection_string = selectionModifier(cutInterpreter.cutString(args.selection + '-mt2llTo100'))
            mc_normalization_weight_string    = MC_WEIGHT(variations[args.variation], returntype='string')
            normalisation_mc = {s.name :s.scale*s.getYieldFromDraw(selectionString = normalization_selection_string, weightString = mc_normalization_weight_string)['val'] for s in mc}

            if args.variation == 'central':
                normalisation_data = data_sample.scale*data_sample.getYieldFromDraw( selectionString = normalization_selection_string, weightString = data_weight_string)['val']
            else:
                normalisation_data = -1

            logger.info( "Making plots.")
            plotting.fill(plots, read_variables = read_variables, sequence = sequence)

            # Delete lambda because we can't serialize it
            for plot in plots:
                del plot.weight

            # save
            dirDB.add( key, (normalisation_mc, normalisation_data, [plot.histos for plot in plots]), overwrite = args.overwrite)

            logger.info( "Done with %s in channel %s.", args.variation, mode)

if args.variation is not None:
    logger.info( "Done with modes %s and variation %s of selection %s. Quit now.", ",".join( modes ), args.variation, args.selection )
    sys.exit(0)

# Systematic pairs:( 'name', 'up', 'down' )
systematics = [\
    {'name':'JEC',         'pair':('jesTotalUp', 'jesTotalDown')},
    {'name':'Unclustered', 'pair':('unclustEnUp', 'unclustEnDown') },
    {'name':'PU',          'pair':('PUUp', 'PUDown')},
    {'name':'BTag_b',      'pair':('BTag_SF_b_Down', 'BTag_SF_b_Up' )},
    {'name':'BTag_l',      'pair':('BTag_SF_l_Down', 'BTag_SF_l_Up')},
    {'name':'trigger',     'pair':('DilepTriggerDown', 'DilepTriggerUp')},
    {'name':'leptonSF',    'pair':('LeptonSFDown', 'LeptonSFUp')},
    #{'name': 'TopPt',     'pair':(  'TopPt', 'central'),},
    #{'name': 'JER',       'pair':(  'JERUp', 'JERDown'),},
]

# loop over modes
missing_cmds   = []
variation_data = {}
for mode in modes:
    logger.info('Working on mode: %s', mode)
    logger.info('Now attempting to load all variations from dirDB %s', dirDB.directory)
   
    for variation in variations.keys():
        key  = (args.era, mode, variation)
        if dirDB.contains(key) and not args.overwrite:
            normalisation_mc, normalisation_data, histos = dirDB.get(key)
            variation_data[(mode, variation)] = {'histos':histos, 'normalisation_mc':normalisation_mc, 'normalisation_data':normalisation_data}
            logger.info( "Loaded normalisations and histograms for variation %s, era %s in mode %s from cache.", variation, args.era, mode)
        else:
            # prepare sub variation command
            cmd = ['python', 'systematicVariation.py']
            cmd.append('--logLevel=%s'%args.logLevel)
            if args.signal is not None: cmd.append( '--signal=%s'%args.signal )
            cmd.append('--plot_directory=%s'%args.plot_directory)
            cmd.append('--selection=%s'%args.selection)
            cmd.append('--variation=%s'%variation)
            if args.small: cmd.append('--small')
            cmd.append('--mode=%s'%args.mode)
            if args.normalizeBinWidth: cmd.append('--normalizeBinWidth')
            cmd.append('--reweightPU=%s'%args.reweightPU)
            cmd.append('--era=%s'%args.era)
            if args.overwrite: cmd.append('--overwrite')

            cmd_string = ' '.join( cmd )
            missing_cmds.append( cmd_string )
            logger.info("Missing variation %s, era %s in mode %s in cache. Need to run: \n%s", variation, args.era, mode, cmd_string)

# write missing cmds
missing_cmds = list(set(missing_cmds))
if len(missing_cmds)>0:
    with file( 'missing.sh', 'w' ) as f:
        f.write("#!/bin/sh\n")
        for cmd in missing_cmds:
            f.write( cmd + '\n')
    logger.info( "Written %i variation commands to ./missing.sh. Now I quit!", len(missing_cmds) )
    sys.exit(0)
    
# make 'all' and 'SF' from ee/mumu/mue
new_modes = []
all_modes = list(modes)
if 'mumu' in modes and 'ee' in modes:
    new_modes.append( ('SF', ('mumu', 'ee')) )
    all_modes.append( 'SF' )
if 'mumu' in modes and 'ee' in modes and 'mue' in modes:
    new_modes.append( ('all', ('mue', 'SF')) )
    all_modes.append( 'all' )
for variation in variations:
    for new_mode, old_modes in new_modes:
        new_key = ( new_mode, variation )
        variation_data[new_key] = {}
        # Adding up data_normalisation 
        if variation == 'central':
            variation_data[new_key]['normalisation_data'] = sum( variation_data[( old_mode, variation )]['normalisation_data'] for old_mode in old_modes )
        else:
            variation_data[new_key]['normalisation_data'] = -1 

        # Adding up mc normalisation
        sample_keys = variation_data[( old_modes[0], variation )]['normalisation_mc'].keys()
        variation_data[new_key]['normalisation_mc'] = {}
        for sample_key in sample_keys: 
            variation_data[new_key]['normalisation_mc'][sample_key] = variation_data[( old_modes[0], variation )]['normalisation_mc'][sample_key]
            for mode in old_modes[1:]:
                variation_data[new_key]['normalisation_mc'][sample_key] += variation_data[( mode, variation )]['normalisation_mc'][sample_key]

        # Adding up histos (clone old_modes[0] at 3rd level, then add)
        variation_data[new_key]['histos'] = [[[ h.Clone() for h in hs ] for hs in plot_histos ] for plot_histos in variation_data[( old_modes[0], variation )]['histos']]
        for mode in old_modes[1:]:
            for i_plot_histos, plot_histos in  enumerate(variation_data[( mode, variation )]['histos']):
                for i_hs, hs in enumerate(plot_histos):
                    for i_h, h in enumerate(hs):
                        variation_data[new_key]['histos'][i_plot_histos][i_hs][i_h].Add(h)
                    
#from RootTools.plot.Plot import addOverFlowBin1D
#for p in plots:
#  p.histos = allPlots[p.name]
#  for s in p.histos:
#    for h in s:
#      addOverFlowBin1D(h, "upper")
#      if h.Integral()==0: logger.warning( "Found empty histogram %s in results file %s", h.GetName(), result_file )


# SF for top central such that we get area normalisation 
dataMC_SF = {}
for mode in all_modes:
    # All SF to 1
    dataMC_SF[mode] = {variation:{s.name:1 for s in mc} for variation in variations} 
    yield_data = variation_data[(mode,'central')]['normalisation_data'] 
    if args.scaling == 'top': 
        # scale variations individually
        if args.variation_scaling:
            logger.info( "Scaling top yield to data for mt2ll<100 individually for all variations." )
            for variation in variations.keys():
                yield_non_top = sum( val for name, val in variation_data[(mode,variation)]['normalisation_mc'].iteritems() if name != Top_pow.name)
                yield_top     = variation_data[(mode,variation)]['normalisation_mc'][Top_pow.name]
                dataMC_SF[mode][variation][Top_pow.name] = (yield_data - yield_non_top)/yield_top
        # scale all variations with the central factor
        else:
            logger.info( "Scaling top yield to data for mt2ll<100 ( all variations are scaled by central SF)" )
            yield_non_top = sum( val for name, val in variation_data[(mode,'central')]['normalisation_mc'].iteritems() if name != Top_pow.name)
            yield_top     = variation_data[(mode,'central')]['normalisation_mc'][Top_pow.name]
            sf = (yield_data - yield_non_top)/yield_top
            for variation in variations.keys():
                dataMC_SF[mode][variation][Top_pow.name] = sf 
    elif args.scaling == 'mc':
        # scale variations individually
        if args.variation_scaling:
            logger.info( "Scaling MC yield to data for mt2ll<100 individually for all variations." )
            for variation in variations.keys():
                yield_mc = sum( val for name, val in variation_data[(mode,variation)]['normalisation_mc'].iteritems())
                for s in mc:
                    dataMC_SF[mode][variation][s.name] = yield_data/yield_mc
        # scale all variations with the central factor
        else:
            logger.info( "Scaling MC yield to data for mt2ll<100 ( all variations are scaled by central SF)" )
            yield_mc = sum( val for name, val in variation_data[(mode,'central')]['normalisation_mc'].iteritems())
            sf = yield_data/yield_mc
            for variation in variations.keys():
                for s in mc:
                    dataMC_SF[mode][variation][s.name] = sf 

def drawObjects( ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Preliminary'),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)'% ( lumi_scale ) ),
      ]
    return [tex.DrawLatex(*l) for l in lines]

# We plot now. 
for mode in all_modes:
    for i_plot, plot in enumerate(plots):
        
        # for central (=no variation), we store plot_data_1, plot_mc_1, plot_data_2, plot_mc_2, ...
        data_histo_list = variation_data[(mode, 'central')]['histos'][2*i_plot]
        mc_histo_list   = {'central': variation_data[(mode, 'central')]['histos'][2*i_plot+1] }
        # for the other variations, there is no data
        for variation in variations.keys():
            if variation=='central': continue
            mc_histo_list[variation] = variation_data[(mode, variation)]['histos'][i_plot]

        # copy styles and tex
        data_histo_list[0][0].style = data_sample.style
        data_histo_list[0][0].legendText = data_sample.texName
        for i_mc_hm, mc_h in enumerate( mc_histo_list['central'][0] ):
            mc_h.style = stack_mc[0][i_mc_hm].style
            mc_h.legendText = stack_mc[0][i_mc_hm].texName

        # perform the scaling
        for variation in variations.keys():
            for s in mc:
                mc_histo_list[variation][0][position[s.name]].Scale( dataMC_SF[mode][variation][s.name] ) 

        # Add histos, del the stack (which refers to MC only )
        plot.histos =  mc_histo_list['central'] + data_histo_list
        plot.stack  = Stack(  mc, [data_sample]) 
        
        # Make boxes and ratio boxes
        boxes           = []
        ratio_boxes     = []
        # Compute all variied MC sums
        total_mc_histo   = {variation:add_histos( mc_histo_list[variation][0]) for variation in variations.keys() }

        # loop over bins & compute shaded uncertainty boxes
        boxes   = []
        r_boxes = []
        for i_b in range(1, 1 + total_mc_histo['central'].GetNbinsX() ):
            # Only positive yields
            total_central_mc_yield = total_mc_histo['central'].GetBinContent(i_b)
            if total_central_mc_yield<=0: continue
            variance = 0.
            for systematic in systematics:
                # Use 'central-variation' (factor 1) and 0.5*(varUp-varDown)
                if 'central' in systematic['pair']: 
                    factor = 1
                else:
                    factor = 0.5
                # sum in quadrature
                variance += ( factor*(total_mc_histo[systematic['pair'][0]].GetBinContent(i_b) - total_mc_histo[systematic['pair'][1]].GetBinContent(i_b)) )**2

            sigma     = sqrt(variance)
            sigma_rel = sigma/total_central_mc_yield 

            box = ROOT.TBox( 
                    total_mc_histo['central'].GetXaxis().GetBinLowEdge(i_b),
                    max([0.03, (1-sigma_rel)*total_central_mc_yield]),
                    total_mc_histo['central'].GetXaxis().GetBinUpEdge(i_b), 
                    max([0.03, (1+sigma_rel)*total_central_mc_yield]) )
            box.SetLineColor(ROOT.kBlack)
            box.SetFillStyle(3444)
            box.SetFillColor(ROOT.kBlack)
            boxes.append(box)

            r_box = ROOT.TBox( 
                total_mc_histo['central'].GetXaxis().GetBinLowEdge(i_b),  
                max(0.1, 1-sigma_rel), 
                total_mc_histo['central'].GetXaxis().GetBinUpEdge(i_b), 
                min(1.9, 1+sigma_rel) )
            r_box.SetLineColor(ROOT.kBlack)
            r_box.SetFillStyle(3444)
            r_box.SetFillColor(ROOT.kBlack)
            ratio_boxes.append(r_box)

        for log in [False, True]:
            plot_directory_ = os.path.join(plot_directory, 'systematicPlots', plot_subdirectory, args.selection, args.era, mode + ("_log" if log else ""))
            #if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
            if    mode == "all": plot.histos[1][0].legendText = "Data (%s)"%args.era
            else:                plot.histos[1][0].legendText = "Data (%s, %s)"%(args.mode, args.era)

            _drawObjects = []

            plotting.draw(plot,
              plot_directory = plot_directory_,
              ratio = {'yRange':(0.1,1.9), 'drawObjects':ratio_boxes},
              logX = False, logY = log, sorting = False,
              yRange = (0.03, "auto") if log else (0.001, "auto"),
              #scaling = {0:1},
              legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
              drawObjects = drawObjects() + boxes,
              copyIndexPHP = True, extensions = ["png"],
            )         
