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
    import tWZ.samples.nanoAODv6_private_locations as locations
    directory_ = locations.mc_2018

logger.info("Loading MC samples from directory %s", directory_)

def make_dirs( dirs ):
    return [ os.path.join( directory_, dir_ ) for dir_ in dirs ]

dirs = {}

# 2016 analysis
## tWZ_sample, TTZtoLLNuNu, TTX_rare_TWZ, TZQ, WZ_amcatnlo, rare, ZZ, nonprompt_TWZ_3l

# diboson
#dirs['diBosonInclusive'] = ["WW", "WZ", "ZZ"]

# diboson exclusive
#dirs['WZ']               = ["WZTo1L1Nu2Q",  "WZTo1L3Nu", "WZTo2L2Q", "WZTo3LNu"]
dirs['WZ']               = ["WZTo1L3Nu", "WZTo2L2Q", "WZTo3LNu_amcatnlo"]
WZ = Sample.fromDirectory(name="WZ", treeName="Events", isData=False, color=color.WZ, texName="WZ", directory=make_dirs( dirs['WZ']))

dirs['ZZ']               = ["ZZTo2L2Q_redBy5", "ZZTo2Q2Nu", "ZZTo4L"]
ZZ = Sample.fromDirectory(name="ZZ", treeName="Events", isData=False, color=color.ZZ, texName="ZZ", directory=make_dirs( dirs['ZZ']))

dirs['WW']               = ["WWToLNuQQ"]
dirs['VVTo2L2Nu']        = ["VVTo2L2Nu"]

# covered by VVTo2L2Nu:
#dirs['WWTo2L2Nu']       = ['WWTo2L2Nu']
#dirs['ZZTo2L2Nu']       = ['ZZTo2L2Nu', "ZZTo2L2Q"]
#dirs['ZZTo2L2Nu_2']     = ['ZZTo2L2Nu']

dirs['triBoson']         = ["WWW_4F","WWZ","WZZ","ZZZ"]
triBoson = Sample.fromDirectory(name="triBoson", treeName="Events", isData=False, color=color.triBoson, texName="VVV", directory=make_dirs( dirs['triBoson']))

# TWZ (nominal)
dirs['TWZ']             = ['tWnunu', 'tWll']
TWZ  = Sample.fromDirectory(name="TWZ", treeName="Events", isData=False, color=color.TWZ, texName="tWZ", directory=make_dirs( dirs['TWZ']))

dirs['TWZ_NLO_DR']             = ['tWll_thad_Wlept_DR', 'tWll_tlept_Wlept_DR', 'tWll_tlept_Whad_DR']
TWZ_NLO_DR  = Sample.fromDirectory(name="TWZ_NLO_DR", treeName="Events", isData=False, color=color.TWZ, texName="tWZ(NLO)", directory=make_dirs( dirs['TWZ_NLO_DR']))

dirs['TWZ_NLO_DS']             = ['tWll_thad_Wlept_DS', 'tWll_tlept_Wlept_DS', 'tWll_tlept_Whad_DS']
TWZ_NLO_DS  = Sample.fromDirectory(name="TWZ_NLO_DS", treeName="Events", isData=False, color=color.TWZ, texName="tWZ(NLO)", directory=make_dirs( dirs['TWZ_NLO_DS']))

# TTZ
dirs['TTZToLLNuNu']     = ['TTZToLLNuNu', 'TTZToLLNuNu_m1to10'] 
dirs['TTZToQQ']         = ['TTZToQQ']
dirs['TTZ']             = ['TTZToLLNuNu', 'TTZToLLNuNu_m1to10', "TTZToQQ"]
TTZToLLNuNu = Sample.fromDirectory(name="ToLLNuNu", treeName="Events", isData=False, color=color.TTZ, texName="t#bar{t}Z #rightarrow ll#nu#nu", directory=make_dirs( dirs['TTZToLLNuNu']))
TTZ         = Sample.fromDirectory(name="TTZ",      treeName="Events", isData=False, color=color.TTZ, texName="t#bar{t}Z", directory=make_dirs( dirs['TTZ']))

# TTX
dirs['TZQ']             = ['tZq_ll']#, 'tZq_nunu'] 
TZQ = Sample.fromDirectory(name="TZQ", treeName="Events", isData=False, color=color.TZQ, texName="tZQ", directory=make_dirs( dirs['TZQ']))

#dirs['TTH']             = ['TTHbbLep', 'TTHnobb_pow']
#dirs['THX']             = ['THW', 'THQ']
#dirs['TTTT']            = ['TTTT']
dirs['TTW']             = ['TTWToLNu', 'TTWToQQ']
dirs['TTVV']            = ['TTWW', 'TTWZ','TTZZ']
dirs['TTX_rare']        = ["TTTT", "TTHbbLep", "TTHnobb_pow", "THW", "THQ"] + dirs['TTW'] + dirs['TTVV'] # same as TTX_rare but without tZq_ll
TTX_rare = Sample.fromDirectory(name="TTX_rare", treeName="Events", isData=False, color=color.TTX_rare, texName="t/t#bar{t}+(t#bar{t}/H/W/VV)", directory=make_dirs( dirs['TTX_rare']))

TTW  = Sample.fromDirectory(name="TTW", treeName="Events", isData=False, color=color.TTW, texName="t#bar{t}W", directory=make_dirs( dirs['TTW']))
dirs['TTTT']            = ['TTTT']
TTTT = Sample.fromDirectory(name="TTTT", treeName="Events", isData=False, color=color.TTTT, texName="t#bar{t}W", directory=make_dirs( dirs['TTTT']))

# TT
dirs['TTLep']           = ['TTLep_pow']
TTLep = Sample.fromDirectory(name="TTLep",      treeName="Events", isData=False, color=color.TTZ, texName="t#bar{t}Lep", directory=make_dirs( dirs['TTLep']))
dirs['singleTop']       = ['TBar_tch_pow_redBy15', 'TBar_tWch', 'T_tch_pow_redBy15', 'T_tWch']
dirs['Top']             = dirs['TTLep'] + dirs['singleTop']
Top  = Sample.fromDirectory(name="Top", treeName="Events", isData=False, color=color.TTJets, texName="t/t#bar{t}", directory=make_dirs(dirs['Top']))

# DY
dirs['DY_M4']           = ['DYJetsToLL_M4to50_HT400to600', 'DYJetsToLL_M4to50_HT200to400', 'DYJetsToLL_M4to50_HT100to200', 'DYJetsToLL_M4to50_HT600toInf', 'DYJetsToLL_M4to50_HT70to100', 'DYJetsToLL_M50_HT100to200']
dirs['DY_HT_LO']        = ['DYJetsToLL_M50_HT70to100', 'DYJetsToLL_M10to50_LO_lheHT70' , 'DYJetsToLL_M50_HT600to800', 'DYJetsToLL_M50_HT800to1200', 'DYJetsToLL_M50_HT1200to2500', 'DYJetsToLL_M50_HT2500toInf', 'DYJetsToLL_M50_HT200to400', 'DYJetsToLL_M50_HT400to600']
dirs['DY_LO']           = [ 'DYJetsToLL_M10to50_LO', 'DYJetsToLL_M50_LO','DYJetsToLL_M50_LO_lheHT70' ]
dirs['DY']              = dirs['DY_M4'] + dirs['DY_HT_LO'] + dirs['DY_LO']
DY  = Sample.fromDirectory(name="DY", treeName="Events", isData=False, color=color.DY, texName="DY", directory=make_dirs(dirs['DY_LO']))

dirs['nonprompt_3l']    = dirs['WW'] + dirs['VVTo2L2Nu'] + dirs['singleTop'] + dirs['DY_LO'] + dirs['TTLep']
dirs['nonprompt_4l']    = dirs['WW'] + dirs['VVTo2L2Nu'] + dirs['singleTop'] + dirs['TZQ'] + dirs["WZ"] + dirs['DY_LO'] + dirs['TTLep'] 
nonprompt_3l = Sample.fromDirectory(name="nonprompt_3l", treeName="Events", isData=False, color=color.nonprompt, texName="nonprompt_3l", directory=make_dirs( dirs['nonprompt_3l']))
nonprompt_4l = Sample.fromDirectory(name="nonprompt_4l", treeName="Events", isData=False, color=color.nonprompt, texName="nonprompt_4l", directory=make_dirs( dirs['nonprompt_4l']))

# Let's see if we need these:
#dirs['pseudoData']      = dirs['TTZtoLLNuNu'] + dirs["WZ_powheg"] + dirs['TTW'] + dirs['TTX'] + dirs['rare'] + dirs['nonprompt_TWZ_3l'] +dirs['ZZ']
#dirs['pseudoDataPriv']  = dirs['ewkDM_ttZ_ll_noH'] + dirs["WZ"] + dirs['TTW'] + dirs['TTX'] + dirs['TZQ'] + dirs['rare'] + dirs['nonprompt']
#dirs['background']      = dirs["WZ_amcatnlo"] + dirs['TTW'] + dirs['TTX'] + dirs['TZQ'] + dirs['rare']


#yt_TWZ_filter   = Sample.fromDirectory(name="yt_TWZ_filter",    treeName="Events", isData=False, color=ROOT.kRed,               texName="TWZ",                                  directory="/afs/hephy.at/data/rschoefbeck01/cmgTuples/TopEFT_PP_2016_mva_v21/dilepTiny/yt_tWZ_filter/")
#yt_TWZ          = Sample.fromDirectory(name="yt_TWZ",           treeName="Events", isData=False, color=ROOT.kRed,               texName="TWZ",                                  directory="/afs/hephy.at/data/rschoefbeck01/cmgTuples/TopEFT_PP_2016_mva_v21/dilepTiny/yt_tWZ/")
#yt_TZZ          = Sample.fromDirectory(name="yt_TZZ",           treeName="Events", isData=False, color=ROOT.kRed,               texName="TZZ",                                  directory="/afs/hephy.at/data/rschoefbeck01/cmgTuples/TopEFT_PP_2016_mva_v22/dilepTiny/yt_tZZ/")
##yt_TWW          = Sample.fromDirectory(name="yt_TWW",           treeName="Events", isData=False, color=ROOT.kRed,               texName="TWW",                                  directory="/afs/hephy.at/data/rschoefbeck01/cmgTuples/TopEFT_PP_2016_mva_v22/dilepTiny/yt_tWW_0/")
#yt_TWW          = Sample.fromDirectory(name="yt_TWW",           treeName="Events", isData=False, color=ROOT.kRed,               texName="TWW",                                  directory="/afs/hephy.at/data/cms04/ttschida/cmgTuples/yt_tWW")
#dim6top_TTW     = Sample.fromDirectory(name="dim6top_TTW",      treeName="Events", isData=False, color=ROOT.kYellow,            texName="TTW",                                  directory="/afs/hephy.at/data/rschoefbeck01/cmgTuples/TopEFT_PP_2016_mva_v21/dilepTiny/dim6top_ttW/")
#yt_TWZ_filter.reweight_pkl = "/afs/hephy.at/data/cms04/ttschida/gridpacks/tWZ01j_new/tWZ01j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_16_tarball.pkl"

## set sample selection strings for the nonprompt and Zgamma sample
#nonpromptMC.setSelectionString('nLeptons_FO_3l_genPrompt<=2')
#ZGTo2LG.setSelectionString('nLeptons_FO_3l_genPrompt>2')
