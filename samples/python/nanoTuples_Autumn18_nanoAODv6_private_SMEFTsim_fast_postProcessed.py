import copy, os, sys
from RootTools.core.Sample import Sample
import ROOT

# Logging
import logging
logger = logging.getLogger(__name__)

from tWZ.samples.color import color

#directory_ = "/scratch-cbe/users/robert.schoefbeck/tWZ/nanoTuples/tWZ_nAODv6_private_v6/2018/trilep/"
directory_ = "/groups/hephy/cms/robert.schoefbeck/tWZ/nanoTuples/tWZ_nAODv6_private_v6/2018/trilep/"
directory_lepWfilter_ = "/scratch-cbe/users/dennis.schwarz/tWZ/nanoTuples/tWZ_nAODv6_private_v6/2018/trilep/"

logger.info("Loading FastSim SMEFTSim MC samples from directory %s", directory_)

def make_dirs( dirs ):
    return [ os.path.join( directory_, dir_ ) for dir_ in dirs ]

def make_dirs_lepWfilter( dirs ):
    return [ os.path.join( directory_lepWfilter_, dir_ ) for dir_ in dirs ]

dirs = {}

dirs['WZ']               = ["WZTo3L1Nu_fast"]
WZ = Sample.fromDirectory(name="WZ", treeName="Events", isData=False, color=color.WZ, texName="WZ", directory=make_dirs( dirs['WZ']))
WZ.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/flavor/vec/WZTo3L1Nu-vec_reweight_card.pkl"

dirs['ZZ']               = ["ZZ_fast"]
ZZ = Sample.fromDirectory(name="ZZ", treeName="Events", isData=False, color=color.ZZ, texName="ZZ", directory=make_dirs( dirs['ZZ']))
ZZ.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/flavor/vec/ZZ-vec_reweight_card.pkl"

dirs['ttZ01j']             = ["ttZ01j_fast"]
ttZ01j         = Sample.fromDirectory(name="ttZ01j",      treeName="Events", isData=False, color=color.TTZ, texName="t#bar{t}Z", directory=make_dirs( dirs['ttZ01j']))
ttZ01j.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/flavor/vec/ttZ01j-vec_reweight_card.pkl"

dirs['ttZ01j_lepWFilter']             = ["ttZ01j_lepWFilter"]
ttZ01j_lepWFilter = Sample.fromDirectory(name="ttZ01j_lepWFilter",      treeName="Events", isData=False, color=color.TTZ, texName="t#bar{t}Z", directory=make_dirs_lepWfilter( dirs['ttZ01j_lepWFilter']))
ttZ01j_lepWFilter.reweight_pkl = "/eos/vbc/group/cms/dennis.schwarz/gridpacks/ttZ01j_ZToLL/ttZ01j-vec_reweight_card.pkl"

dirs['tWZToLL01j_lepWFilter']             = ["tWZToLL01j_lepWFilter"]
tWZToLL01j_lepWFilter = Sample.fromDirectory(name="tWZToLL01j_lepWFilter",      treeName="Events", isData=False, color=color.TWZ, texName="tWZ", directory=make_dirs_lepWfilter( dirs['tWZToLL01j_lepWFilter']))
tWZToLL01j_lepWFilter.reweight_pkl = "/eos/vbc/group/cms/dennis.schwarz/gridpacks/tWZ01j_ZtoLL/tWZ-vec_reweight_card.pkl"
