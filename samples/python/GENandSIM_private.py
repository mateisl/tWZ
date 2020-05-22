'''
GENSIM private production
'''

import copy, os, sys
from RootTools.fwlite.FWLiteSample import FWLiteSample
import ROOT

def get_parser():
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser for samples file")
    argParser.add_argument('--overwrite',   action='store_true',    help="Overwrite current entry in db?")
    return argParser

# Logging
if __name__=="__main__":
    import Samples.Tools.logger as logger
    logger = logger.get_logger("INFO", logFile = None )
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger("INFO", logFile = None )
    options = get_parser().parse_args()
    ov = options.overwrite

else:
    import logging
    logger = logging.getLogger(__name__)
    ov = False

from Samples.Tools.config import dbDir, redirector, redirector_global
dbFile = dbDir+"DB_tWZ_GENSIM.sql"

logger.info("Using db file: %s", dbFile)

tWZ01j_rwgt_filter = FWLiteSample.fromDAS("tWZ01j_rwgt_filter", "/tWZ01j_rwgt_filter_2/ttschida-tWZ01j_rwgt_filter_2_RAWSIMoutput-d21de898ad3d5a9f1861cf5edf14f0e3/USER", instance="phys03", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
#tOrtbar_WZ01j_OLRLL_LO = FWLiteSample.fromDAS("tOrtbar_WZ01j_OLRLL_LO", "/tOrtbar_WZ01j_OLRLL_LO/schoef-tOrtbar_WZ01j_OLRLL_LO-129aa8aa5b17c813bcab74208422fa28/USER", instance="phys03", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)

# GENSIM
#tOrtbar_WZ01j_OLRLL_LO = FWLiteSample.fromDAS("tOrtbar_WZ01j_OLRLL_LO", "/gensim_tOrtbar_WZ01j_OLRLL_LO/schoef-gensim_tOrtbar_WZ01j_OLRLL_LO_RAWSIMoutput-f6a5011659c881b9d217cfc35ad761a3/USER", instance="phys03", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
#tOrtbar_WZ01j_OLRLL_LO = FWLiteSample.fromDAS("tOrtbar_WZ01j_OLRLL_LO", "/tOrtbar_WZ01j_OLRLL_LO_ext4/schoef-tOrtbar_WZ01j_OLRLL_LO_ext4-c88179374d61ba4483b39a09a725f437/USER", instance="phys03", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)

tOrtbar_WZ01j_OLRLL_LO = FWLiteSample.fromDAS("tOrtbar_WZ01j_OLRLL_LO", "/tOrtbar_WZ01j_OLRLL_LO_ext5/schoef-tOrtbar_WZ01j_OLRLL_LO_ext5-c88179374d61ba4483b39a09a725f437/USER", instance="phys03", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)

tWZ_NLO = FWLiteSample.fromDAS("tWZ_NLO", "/tWZ_NLO_v16_ext3/schoef-tWZ_NLO_v16_ext3-a1e791adef4fc0a2f32451f12fdbd583/USER", instance="phys03", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)

#tWZ_LO_central = FWLiteSample.fromDAS("tWZ_LO_central", "/ST_tWll_5f_LO_13TeV-MadGraph-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
tWZ_LO_central = FWLiteSample.fromDAS("tWZ_LO_central", "/ST_tWll_5f_LO_13TeV-MadGraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)

allSamples = [tWZ01j_rwgt_filter, tOrtbar_WZ01j_OLRLL_LO, tWZ_NLO, tWZ_LO_central] 

for sample in allSamples:
    sample.isData = False

from Samples.Tools.AutoClass import AutoClass
samples = AutoClass( allSamples )
