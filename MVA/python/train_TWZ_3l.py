#!/usr/bin/env python

# Analysis
from Analysis.TMVA.Trainer       import Trainer
from Analysis.TMVA.Reader        import Reader
from Analysis.TMVA.defaults      import default_methods, default_factory_settings 
import Analysis.Tools.syncer

# TopEFT
from tWZ.Tools.user              import plot_directory, mva_directory
from tWZ.Tools.cutInterpreter    import cutInterpreter

# MVA configuration
from tWZ.MVA.MVA_TWZ_3l          import sequence, read_variables, mva_variables, all_mva_variables 
from tWZ.MVA.MVA_TWZ_3l          import * 

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--plot_directory',     action='store',             default=None)
argParser.add_argument('--selection',          action='store', type=str,   default='trilepM-onZ1')
argParser.add_argument('--trainingFraction',   action='store', type=float, default=0.5)
argParser.add_argument('--small',              action='store_true')
argParser.add_argument('--overwrite',          action='store_true')

args = argParser.parse_args()

#Logger
import tWZ.Tools.logger as logger
logger = logger.get_logger("INFO", logFile = None )
import Analysis.Tools.logger as logger_an
logger_an = logger_an.get_logger("INFO", logFile = None )

if args.plot_directory == None:
    args.plot_directory = plot_directory

if args.selection == None:
    selectionString = "(1)"
else:
    selectionString = cutInterpreter.cutString( args.selection )

# Samples
#from tWZ.samples.nanoTuples_RunII_nanoAODv6_private_postProcessed    import *
from tWZ.samples.nanoTuples_Summer16_nanoAODv6_private_postProcessed import *

signal = TWZ_NLO_DR 
#signal.reduceFiles(factor=20)
#TTZ.reduceFiles(factor=3)

# TTZ
backgrounds = [ TTZ ]

samples = backgrounds + [signal]
for sample in samples:
    sample.setSelectionString( selectionString )
    if args.small:
        sample.reduceFiles(to = 1)

#mvas = [mlp_np7,mlp_np10, mlp_np20, mlp_np30, mlp_np40]
#mvas = [all_bdt_nTr1000_maxD1_mNS5, all_bdt_nTr1000_maxD1_mNS10, all_bdt_nTr1000_maxD3_mNS10, all_bdt_nTr1000_maxD4_mNS20,
#        all_mlp_np5s0c5e0c8, all_mlp_np5s0c3e0c5, all_mlp_np5s0c5e1, all_mlp_np7s0c3e0c5, all_mlp_np7s0c5e1, all_mlp_np7s0c5e0c8, all_mlp_ncnc1s0c3e0c5, all_mlp_ncnc1s0c5e1, all_mlp_ncnc1s0c5e0c8,
#        all_mlp_ncnp5c1s0c3e0c5, all_mlp_ncnp5c1s0c5e1, all_mlp_ncnp5c1s0c5e0c8, all_mlp_np20, all_mlp_np30, all_mlp_np40, all_mlp_oldconfig_ncnp5c1s0c3e0c5, all_mlp_oldconfig_ncnp5c1s0c3e0c3, all_mlp_oldconfig_np7s0c3e0c8, all_mlp_oldconfig_np7c1s0c5e0c5 ]
mvas = [ all_mlp_oldconfig_ncnp5c1s0c3e0c3 ] #, all_mlp_oldconfig_np7s0c3e0c8, all_mlp_oldconfig_np7c1s0c5e0c5 ] 

## TMVA Trainer instance
trainer = Trainer( 
    signal = signal, 
    backgrounds = backgrounds, 
    output_directory = mva_directory, 
    mva_variables    = all_mva_variables,
    label            = "TWZ_3l", 
    fractionTraining = args.trainingFraction, 
    )

weightString = "(1)"
trainer.createTestAndTrainingSample( 
    read_variables   = read_variables,   
    sequence         = sequence,
    weightString     = weightString,
    overwrite        = args.overwrite,
    mva_variables    = all_mva_variables 
    )

#trainer.addMethod(method = default_methods["BDT"])
#trainer.addMethod(method = default_methods["MLP"])

for mva in mvas:
    trainer.addMethod(method = mva)

trainer.trainMVA( factory_settings = default_factory_settings )
trainer.plotEvaluation(plot_directory = os.path.join( plot_directory, "MVA", "TWZ_3l_all") )

#reader.addMethod(method = bdt1)
#reader.addMethod(method = default_methods["MLP"])

#print reader.evaluate("BDT", met_pt=100, ht=-210, Z1_pt_4l=100, lnonZ1_pt=100, lnonZ1_eta=0)
#prinMt reader.evaluate("BDT", met_pt=120, ht=-210)
