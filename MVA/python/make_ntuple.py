#!/usr/bin/env python

# General
import os
import ROOT

# Analysis
#import Analysis.Tools.syncer
# RootTools
from RootTools.core.standard import *
# TopEFT
from tWZ.Tools.cutInterpreter    import cutInterpreter

# MVA configuration
from tWZ.MVA.MVA_TWZ_3l          import sequence, read_variables, mva_variables, all_mva_variables 

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--selection',          action='store', type=str,   default='trilepM-onZ1')
argParser.add_argument('--sample',             action='store', type=str,   default='TWZ_NLO_DR')
argParser.add_argument('--output_directory',   action='store', type=str,   default='.')
argParser.add_argument('--small',              action='store_true')

args = argParser.parse_args()

#Logger
import tWZ.Tools.logger as logger
logger = logger.get_logger("INFO", logFile = None )
import Analysis.Tools.logger as logger_an
logger_an = logger_an.get_logger("INFO", logFile = None )

# get samples
import tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed as samples
sample = getattr(samples, args.sample)

if args.small:
    sample.reduceFiles(to=1)

# selection
if args.selection == None:
    selectionString = "(1)"
else:
    selectionString = cutInterpreter.cutString( args.selection )
sample.setSelectionString( selectionString )

count  = int(sample.getYieldFromDraw( weightString="(1)" )["val"])
logger.info( "Found %i events for sample %s", count, sample.name )

# where the output goes
output_file  = os.path.join( args.output_directory, sample.name + ("_small" if args.small else "") + ".root" )

# reader
reader = sample.treeReader( \
    variables = map( TreeVariable.fromString, read_variables),
    )
reader.start()

#filler
def filler( event ):

    # fill extra variables
    #event.isTraining = isTraining
    #event.isSignal   = isSignal
    # write mva variables
    for name, func in mva_variables.iteritems():
#                setattr( event, name, func(reader.event) )
        setattr( event, name, func(reader.event, sample=None) )

# Create a maker. Maker class will be compiled. 
maker = TreeMaker(
    sequence  = [ filler ],
    variables = map(TreeVariable.fromString, 
#          ["isTraining/I", "isSignal/I"] + 
          ["%s/F"%var for var in mva_variables.keys()] 
        ),
    treeName = "Events"
    )
maker.start()

logger.info( "Starting event loop" )
counter=0
while reader.run():
    for func in sequence:
        func(reader.event)

    ## determine whether training or test
    #isTraining = self.samples[i_sample].training_test_list.pop(0)
    #isSignal   = (i_sample == 0)

    maker.run()
    counter += 1
    if counter%10000 == 0:
        logger.info("Written %i events.", counter)

nEventsTotal = maker.tree.GetEntries()

tmp_directory = ROOT.gDirectory
dirname = os.path.dirname(output_file)
if not os.path.exists(dirname):
    os.makedirs(dirname)

outputfile = ROOT.TFile.Open(output_file, 'recreate')
maker.tree.Write()
outputfile.Close()
tmp_directory.cd()
logger.info( "Written %s", output_file)
#
#      # Destroy the TTree
maker.clear()
logger.info( "Written %i events to %s",  nEventsTotal, output_file )
