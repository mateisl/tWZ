import copy, os, sys
from RootTools.core.Sample import Sample
import ROOT

def get_parser():
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser for samples file")
    argParser.add_argument('--overwrite',   action='store_true',    help="Overwrite current entry in db?")
    argParser.add_argument('--update',      action='store_true',    help="Update current entry in db?")
    argParser.add_argument('--check_completeness', action='store_true',    help="Check competeness?")
    return argParser

# Logging
if __name__=="__main__":
    import Samples.Tools.logger as logger
    logger = logger.get_logger("INFO", logFile = None )
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger("INFO", logFile = None )
    options = get_parser().parse_args()
    ov = options.overwrite
    if options.update:
        ov = 'update'
else:
    import logging
    logger = logging.getLogger(__name__)
    ov = False

# Redirector
try:
    redirector = sys.modules['__main__'].redirector
except:
    if "clip" in os.getenv("HOSTNAME").lower():
        from Samples.Tools.config import redirector_clip_local as redirector
    else:
        from Samples.Tools.config import redirector as redirector

#from Samples.Tools.config import  redirector_global


# DB
from Samples.Tools.config import dbDir
dbFile = dbDir+"/DB_tWZ_private_legacy.sql"

logger.info("Using db file: %s", dbFile)

yt_tWZ01j = Sample.nanoAODfromDAS("yt_tWZ01j", "/tWZ01j_rwgt/schoef-crab_ttschida-Summer16-mAOD949-bd3e7bcff6c9bcad356ea4ed7e4f08b4_legacy_nano_v4-b9659cf3bef5e21efe24288a402778f7/USER", dbFile=dbFile, redirector=redirector, instance="phys03", overwrite=ov, xSection=0.2279)
yt_tWZ01j.reweight_pkl = '/afs/hephy.at/data/rschoefbeck01/gridpacks/Yt/tZZ1j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_16_tarball.pkl'

tWZ01j_filter_efficiency = (601438./(30*10**6))*(10**6/363784.)

yt_tWZ01j_filter = Sample.nanoAODfromDAS("yt_tWZ01j_filter", "/tWZ01j_rwgt_filter_2/schoef-crab_ttschida-Summer16-mAOD949-bd3e7bcff6c9bcad356ea4ed7e4f08b4_legacy_nano_v4-b9659cf3bef5e21efe24288a402778f7/USER", dbFile=dbFile, redirector=redirector, instance="phys03", overwrite=ov, xSection=0.2279*tWZ01j_filter_efficiency)
yt_tWZ01j_filter.reweight_pkl = '/afs/hephy.at/data/rschoefbeck01/gridpacks/Yt/tZZ1j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_16_tarball.pkl'#same

tWZ01j_private = [yt_tWZ01j, yt_tWZ01j_filter]

allSamples = tWZ01j_private 

for s in allSamples:
    s.isData = False

from Samples.Tools.AutoClass import AutoClass
samples = AutoClass( allSamples )
if __name__=="__main__":
    if options.check_completeness:
        samples.check_completeness( cores=20 )

