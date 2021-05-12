#!/usr/bin/env python
''' analysis script for standard plots with systematic errors
'''

# Standard imports and batch mode
import ROOT
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas() # do this to avoid version conflict in png.h with keras import ...
c1.Draw()
c1.Print('delete.png')
import operator
import pickle, os, time, sys
from math                                import sqrt, cos, sin, pi, atan2

# RootTools
from RootTools.core.standard             import *

# tWZ
from tWZ.Tools.user                      import plot_directory
from tWZ.Tools.cutInterpreter            import cutInterpreter
from tWZ.Tools.objectSelection           import cbEleIdFlagGetter, vidNestedWPBitMapNamingList
from tWZ.Tools.objectSelection           import lepString
from tWZ.Tools.helpers                   import add_histos

# Analysis
from Analysis.Tools.helpers              import deltaPhi, deltaR
from Analysis.Tools.puProfileCache       import *
from Analysis.Tools.puReweighting        import getReweightingFunction
from Analysis.Tools.DirDB                import DirDB
import Analysis.Tools.syncer
import numpy as np

################################################################################
# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
argParser.add_argument('--plot_directory', action='store', default='tWZ_sys_v1')
argParser.add_argument('--era',            action='store', type=str, default="Run2016")
argParser.add_argument('--selection',      action='store', default='trilepT-minDLmass12-onZ1-njet4p-btag1')
argParser.add_argument('--signal',            action='store',      default=None,        nargs='?', choices=['None', "T2tt",'DM'], help="Add signal to plot")
argParser.add_argument('--variation',         action='store',      default=None, help="Which systematic variation to run. Don't specify for producing plots.")
argParser.add_argument('--overwrite',         action='store_true',     help='Overwrite?')
argParser.add_argument('--normalizeBinWidth', action='store_true', default=False,       help='normalize wider bins?')
argParser.add_argument('--reweightPU',         action='store', default='Central', choices=[ 'Central', 'VUp'] )
args = argParser.parse_args()

################################################################################
# Logger
import tWZ.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

logger.info( "Working in era %s", args.era)

################################################################################
# Set names of PU weights
if args.reweightPU == 'Central':
    nominalPuWeight, upPUWeight, downPUWeight = "reweightPU", "reweightPUUp", "reweightPUDown"
elif args.reweightPU == 'VUp':
    nominalPuWeight, upPUWeight, downPUWeight = "reweightPUVUp", "reweightPUVVUp", "reweightPUUp"

################################################################################
# Selection modifier
def jetSelectionModifier( sys, returntype = "func"):
    #Need to make sure all jet variations of the following observables are in the ntuple
    variiedJetObservables = ['nJetGood', 'nBTag', 'met_pt']
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
    variiedMetObservables = ['met_pt']
    if returntype == "func":
        def changeCut_( string ):
            for s in variiedMetObservables:
                string = string.replace(s, s+'_'+sys)
            return string
        return changeCut_
    elif returntype == "list":
        return [ v+'_'+sys for v in variiedMetObservables ]

################################################################################
# Nominal MC weights that need to be applied
nominalMCWeights = ["weight", "reweightPU", "reweightTrigger", "reweightBTag_SF", "reweightL1Prefire"]

################################################################################
# Change the weight according to variation
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

################################################################################
# Data weight has no variations
def data_weight( event, sample ):
    return event.weight

data_weight_string = "weight"

################################################################################
# Define systematic variations
variations = {
    'central'           : {'read_variables': [ '%s/F'%v for v in nominalMCWeights ]},
    'jesTotalUp'        : {'selectionModifier':jetSelectionModifier('jesTotalUp'),               'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('jesTotalUp','list')]},
    'jesTotalDown'      : {'selectionModifier':jetSelectionModifier('jesTotalDown'),             'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('jesTotalDown','list')]},
    # 'unclustEnUp'       : {'selectionModifier':metSelectionModifier('unclustEnUp'),              'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('unclustEnUp','list')]},
    # 'unclustEnDown'     : {'selectionModifier':metSelectionModifier('unclustEnDown'),            'read_variables' : [ '%s/F'%v for v in nominalMCWeights + jetSelectionModifier('unclustEnDown','list')]},
    'PUUp'              : {'replaceWeight':(nominalPuWeight,upPUWeight),                         'read_variables' : [ '%s/F'%v for v in nominalMCWeights + [upPUWeight] ]},
    'PUDown'            : {'replaceWeight':(nominalPuWeight,downPUWeight),                       'read_variables' : [ '%s/F'%v for v in nominalMCWeights + [downPUWeight] ]},
    'BTag_SF_b_Down'    : {'replaceWeight':('reweightBTag_SF','reweightBTag_SF_b_Down'),         'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightBTag_SF_b_Down']]},
    'BTag_SF_b_Up'      : {'replaceWeight':('reweightBTag_SF','reweightBTag_SF_b_Up'),           'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightBTag_SF_b_Up'] ]},
    'BTag_SF_l_Down'    : {'replaceWeight':('reweightBTag_SF','reweightBTag_SF_l_Down'),         'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightBTag_SF_l_Down']]},
    'BTag_SF_l_Up'      : {'replaceWeight':('reweightBTag_SF','reweightBTag_SF_l_Up'),           'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightBTag_SF_l_Up'] ]},
    'TriggerDown'       : {'replaceWeight':('reweightTrigger','reweightTriggerDown'),            'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightTriggerDown']]},
    'TriggerUp'         : {'replaceWeight':('reweightTrigger','reweightTriggerUp'),              'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightTriggerUp']]},
    # 'LeptonSFDown'      : {'replaceWeight':('reweightLeptonSF','reweightLeptonSFDown'),          'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightLeptonSFDown']]},
    # 'LeptonSFUp'        : {'replaceWeight':('reweightLeptonSF','reweightLeptonSFUp'),            'read_variables' : [ '%s/F'%v for v in nominalMCWeights + ['reweightLeptonSFUp']]},
#    'TopPt':{},
#   'JERUp':{},
#   'JERDown':{},
}

# Define which systematics change the jet v4
jet_systematics = ['jesTotalUp', 'jesTotalDown', 'JERUp', 'JERDown']

################################################################################
# Add a default selection modifier that does nothing
for key, variation in variations.iteritems():
    if not variation.has_key('selectionModifier'):
        variation['selectionModifier'] = lambda string:string
    if not variation.has_key('read_variables'):
        variation['read_variables'] = []

################################################################################
# Check if we know the variation
if args.variation is not None and args.variation not in variations.keys():
    raise RuntimeError( "Variation %s not among the known: %s", args.variation, ",".join( variation.keys() ) )

################################################################################
# arguments & directory
plot_subdirectory = args.plot_directory
if args.small:                    plot_subdirectory += "_small"
if args.reweightPU:               plot_subdirectory += "_reweightPU%s"%args.reweightPU

################################################################################
# Define MC samples
from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed import *

if args.era == "Run2016":
    mc = [Summer16.TWZ_NLO_DR, Summer16.TTZ, Summer16.TTX_rare, Summer16.TZQ, Summer16.WZ, Summer16.triBoson, Summer16.ZZ, Summer16.nonprompt_3l]
elif args.era == "Run2017":
    mc = [Fall17.TWZ_NLO_DR, Fall17.TTZ, Fall17.TTX_rare, Fall17.TZQ, Fall17.WZ, Fall17.triBoson, Fall17.ZZ, Fall17.nonprompt_3l]
elif args.era == "Run2018":
    mc = [Autumn18.TWZ_NLO_DR, Autumn18.TTZ, Autumn18.TTX_rare, Autumn18.TZQ, Autumn18.WZ, Autumn18.triBoson, Autumn18.ZZ, Autumn18.nonprompt_3l]
elif args.era == "RunII":
    mc = [TWZ_NLO_DR, TTZ, TTX_rare, TZQ, WZ, triBoson, ZZ, nonprompt_3l]

################################################################################
# Define data sample
try:
  data_sample = eval(args.era)
except Exception as e:
  logger.error( "Didn't find %s", args.era )
  raise e

lumi_scale                 = data_sample.lumi/1000
data_sample.scale          = 1.
for sample in mc:
    sample.scale           = 1 # Scale MCs individually with lumi

if args.small:
    for sample in mc + [data_sample]:
        sample.normalization = 1.
        sample.reduceFiles( to = 1 )
        sample.scale /= sample.normalization

################################################################################
# postions of MC components in list
position = {s.name:i_s for i_s,s in enumerate(mc)}


################################################################################
# Read variables
read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I",
    "l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
    "l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
    "JetGood[pt/F,eta/F,phi/F,area/F,btagDeepB/F,index/I]",
    "Jet[pt/F,eta/F,phi/F,mass/F]",
    "lep[pt/F,eta/F,phi/F,pdgId/I,muIndex/I,eleIndex/I]",
    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I",
    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
]

################################################################################
# Sequences
sequence = []

################################################################################
# Define 3l selections
mu_string  = lepString('mu','VL')
ele_string = lepString('ele','VL')
def getLeptonSelection( mode ):
    if   mode=="mumumu": return "Sum$({mu_string})==3&&Sum$({ele_string})==0".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="mumue":  return "Sum$({mu_string})==2&&Sum$({ele_string})==1".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="muee":   return "Sum$({mu_string})==1&&Sum$({ele_string})==2".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=="eee":    return "Sum$({mu_string})==0&&Sum$({ele_string})==3".format(mu_string=mu_string,ele_string=ele_string)
    elif mode=='all':    return "Sum$({mu_string})+Sum$({ele_string})==3".format(mu_string=mu_string,ele_string=ele_string)

################################################################################
# Getter functor for lepton quantities
def lep_getter( branch, index, abs_pdg = None, functor = None, debug=False):
    if functor is not None:
        if abs_pdg == 13:
            def func_( event, sample ):
                if debug:
                    print "Returning", "Muon_%s"%branch, index, abs_pdg, "functor", functor, "result",
                    print functor(getattr( event, "Muon_%s"%branch )[event.lep_muIndex[index]]) if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
                return functor(getattr( event, "Muon_%s"%branch )[event.lep_muIndex[index]]) if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
        else:
            def func_( event, sample ):
                return functor(getattr( event, "Electron_%s"%branch )[event.lep_eleIndex[index]]) if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
    else:
        if abs_pdg == 13:
            def func_( event, sample ):
                if debug:
                    print "Returning", "Muon_%s"%branch, index, abs_pdg, "functor", functor, "result",
                    print getattr( event, "Muon_%s"%branch )[event.lep_muIndex[index]] if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
                return getattr( event, "Muon_%s"%branch )[event.lep_muIndex[index]] if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
        else:
            def func_( event, sample ):
                return getattr( event, "Electron_%s"%branch )[event.lep_eleIndex[index]] if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
    return func_


################################################################################
# Set up channels and values for plotting
yields     = {}
allPlots   = {}
modes      = ['mumumu','mumue','muee', 'eee']

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

################################################################################
# Fire up the cache
dirDB = DirDB(os.path.join(plot_directory, 'systematicPlots', plot_subdirectory, args.selection, 'cache'))

################################################################################
# loop over modes
for mode in modes:
    yields[mode] = {}
    data_sample.texName = "data"
    data_sample.setSelectionString([getLeptonSelection(mode)])
    data_sample.name           = "data"
    data_sample.style          = styles.errorStyle(ROOT.kBlack)
    lumi_scale                 = data_sample.lumi/1000

    for sample in mc: sample.style = styles.fillStyle(sample.color)

    for sample in mc:
      sample.setSelectionString([getLeptonSelection(mode)])

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


    ############################################################################
    # Set up plots
    plots      = []

    if args.variation == 'central':
        l1_pt_data   = Plot(
            name        = "l1_pt_data",
            texX        = 'p_{T}(leading lepton) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
            binning     = [300/20,0,300],
            stack       = stack_data,
            attribute   = TreeVariable.fromString( "l1_pt/F" ),
            weight      = data_weight )
        plots.append( l1_pt_data )

    l1_pt_mc  = Plot(\
        name            = "l1_pt_mc",
        texX            = 'p_{T}(leading lepton) (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
        binning         = [300/20,0,300],
        stack           = stack_mc,
        attribute       = TreeVariable.fromString( "l1_pt/F" ),
        weight          = mc_weight )
    plots.append( l1_pt_mc )

    if args.variation == 'central':
        met_data   = Plot(
            name        = "met_data",
            texX        = 'p_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
            binning     = [400/20,0,400],
            stack       = stack_data,
            attribute   = TreeVariable.fromString( "met_pt/F" ),
            weight      = data_weight )
        plots.append( met_data )

    met_mc   = Plot(
        name        = "met_mc",
        texX        = 'p_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV' if args.normalizeBinWidth else "Number of Events",
        binning     = [400/20,0,400],
        stack       = stack_mc,
        # attribute   = TreeVariable.fromString('met_pt/F'),
        attribute   = TreeVariable.fromString( "met_pt_%s/F" % args.variation ) if args.variation in jet_systematics else TreeVariable.fromString('met_pt/F'),
        selectionString = selectionModifier(cutInterpreter.cutString(args.selection)) if selectionModifier is not None else None,
        weight      = mc_weight )
    plots.append( met_mc )

    ############################################################################
    # Check DB for existing plots
    if args.variation is not None:
        key  = (args.era, mode, args.variation)
        if dirDB.contains(key) and not args.overwrite:
            normalisation_mc, normalisation_data, histos = dirDB.get( key )
            for i_p, h_s in enumerate(histos):
                plots[i_p].histos = h_s
            logger.info( "Loaded normalisations and histograms for %s in mode %s from cache.", args.era, mode)
        else:
            logger.info( "Obtain normalisations and histograms for %s in mode %s.", args.era, mode)
            normalization_selection_string = selectionModifier(cutInterpreter.cutString(args.selection))
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

################################################################################
# Systematic pairs:( 'name', 'up', 'down' )
systematics = [\
    {'name':'JEC',         'pair':('jesTotalUp', 'jesTotalDown')},
    # {'name':'Unclustered', 'pair':('unclustEnUp', 'unclustEnDown') },
    {'name':'PU',          'pair':('PUUp', 'PUDown')},
    {'name':'BTag_b',      'pair':('BTag_SF_b_Down', 'BTag_SF_b_Up' )},
    {'name':'BTag_l',      'pair':('BTag_SF_l_Down', 'BTag_SF_l_Up')},
    {'name':'trigger',     'pair':('TriggerDown', 'TriggerUp')},
    # {'name':'leptonSF',    'pair':('LeptonSFDown', 'LeptonSFUp')},
    # {'name': 'TopPt',     'pair':(  'TopPt', 'central'),},
    # {'name': 'JER',       'pair':(  'JERUp', 'JERDown'),},
]

################################################################################
# Check if all files are alread in DB, otherwise create cmd strings
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
            cmd = ['python', 'systematics.py']
            cmd.append('--logLevel=%s'%args.logLevel)
            if args.signal is not None: cmd.append( '--signal=%s'%args.signal )
            cmd.append('--plot_directory=%s'%args.plot_directory)
            cmd.append('--selection=%s'%args.selection)
            cmd.append('--variation=%s'%variation)
            if args.small: cmd.append('--small')
            if args.normalizeBinWidth: cmd.append('--normalizeBinWidth')
            cmd.append('--reweightPU=%s'%args.reweightPU)
            cmd.append('--era=%s'%args.era)
            if args.overwrite: cmd.append('--overwrite')

            cmd_string = ' '.join( cmd )
            missing_cmds.append( cmd_string )
            logger.info("Missing variation %s, era %s in mode %s in cache. Need to run: \n%s", variation, args.era, mode, cmd_string)

################################################################################
# write cmd strings in missing.sh
missing_cmds = list(set(missing_cmds))
if len(missing_cmds)>0:
    with file( 'missing.sh', 'w' ) as f:
        f.write("#!/bin/sh\n")
        for cmd in missing_cmds:
            f.write( cmd + '\n')
    logger.info( "Written %i variation commands to ./missing.sh. Now I quit!", len(missing_cmds) )
    sys.exit(0)


################################################################################
# make 'all' and 'SF' modes
# usage: new_modes.append( ('new name', (list modes to sum up)) )
new_modes = []
all_modes = list(modes)
if 'mumumu' in modes and 'mumue' in modes and 'muee' in modes and 'eee' in modes:
    new_modes.append( ('all', ('mumumu', 'mumue', 'muee', 'eee')) )
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



################################################################################
# SF for top central such that we get area normalisation
dataMC_SF = {}
for mode in all_modes:
    # All SF to 1
    # The dict will look like dataMC_SF[mode][variationname][samplename]
    dataMC_SF[mode] = {variation:{s.name:1 for s in mc} for variation in variations}
    yield_data = variation_data[(mode,'central')]['normalisation_data']

################################################################################
# Draw
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

################################################################################
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
                # print 'DATA/MC scaling:'
                # print dataMC_SF[mode][variation][s.name]
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
            else:                plot.histos[1][0].legendText = "Data (%s, %s)"%(mode, args.era)

            # for i in range(len[plot.histos]):
            #     print len(plot.histos[i])

            _drawObjects = []
            plotting.draw(plot,
              plot_directory = plot_directory_,
              ratio = {'yRange':(0.1,1.9), 'drawObjects':ratio_boxes},
              logX = False, logY = log, sorting = False,
              yRange = (0.03, "auto") if log else (0.001, "auto"),
              #scaling = {0:1}, # scaling is already preformed above
              legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
              drawObjects = drawObjects() + boxes,
              copyIndexPHP = True, extensions = ["png"],
            )
