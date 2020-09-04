import copy, os, sys
from RootTools.core.Sample import Sample
import ROOT

# Logging
import logging
logger = logging.getLogger(__name__)

from tWZ.samples.color import color

# Data directory
try:
    directory_ = sys.modules['__main__'].directory_
except:
    import tWZ.samples.nanoAODv6_locations as locations
    directory_ = locations.mc_2016

logger.info("Loading MC samples from directory %s", directory_)

def make_dirs( dirs ):
    return [ os.path.join( directory_, dir_ ) for dir_ in dirs ]

dirs = {}
# tWZ

tWZ_LO   = Sample.fromDirectory(name="tWZ_LO", treeName="Events", isData=False, color=ROOT.kMagenta, texName="tWZ_LO", directory= make_dirs(['tWZ_LO']))
tWZ_NLO  = Sample.fromDirectory(name="tWZ_NLO", treeName="Events", isData=False, color=ROOT.kGreen, texName="tWZ_NLO", directory= make_dirs(['tWZ_NLO'])) 

#dirs['tWZ_DR']       = ["tWll_thad_Wlept_DR","tWll_tlept_Whad_DR","tWll_tlept_Wlept_DR"]
tWZ_DR               = Sample.fromDirectory(name="tWZ_DR", treeName="Events", isData=False, color=ROOT.kAzure-1, texName="tWZ_DR", directory= make_dirs(["tWll_thad_Wlept_DR","tWll_tlept_Whad_DR","tWll_tlept_Wlept_DR"]))

#dirs['tWZ_DS']       = ["tWll_thad_Wlept_DS","tWll_tlept_Whad_DS","tWll_tlept_Wlept_DS"]
tWZ_DS               = Sample.fromDirectory(name="tWZ_DS", treeName="Events", isData=False, color=ROOT.kSpring+3, texName="tWZ_DS", directory= make_dirs(["tWll_thad_Wlept_DS","tWll_tlept_Whad_DS","tWll_tlept_Wlept_DS"]))

