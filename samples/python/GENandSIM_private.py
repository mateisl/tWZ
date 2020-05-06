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

tOrtbar_WZ01j_OLRLL_LO = FWLiteSample.fromDAS("tOrtbar_WZ01j_OLRLL_LO", "/gensim_tOrtbar_WZ01j_OLRLL_LO/schoef-gensim_tOrtbar_WZ01j_OLRLL_LO_RAWSIMoutput-f6a5011659c881b9d217cfc35ad761a3/USER", instance="phys03", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)

allSamples = [tWZ01j_rwgt_filter, tOrtbar_WZ01j_OLRLL_LO] 

for sample in allSamples:
    sample.isData = False

from Samples.Tools.AutoClass import AutoClass
samples = AutoClass( allSamples )
