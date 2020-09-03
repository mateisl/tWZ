from RootTools.core.standard import *

from tWZ.samples.nanoTuples_Run2016_nanoAODv6_private_postProcessed import Run2016
from tWZ.samples.nanoTuples_Run2017_nanoAODv6_private_postProcessed import Run2017
from tWZ.samples.nanoTuples_Run2018_nanoAODv6_private_postProcessed import Run2018

RunII      = Sample.combine( "RunII", [Run2016, Run2017, Run2018] ) 
RunII.lumi = Run2016.lumi + Run2017.lumi + Run2018.lumi

lumi_year  = {2016:Run2016.lumi, 2017:Run2017.lumi, 2018:Run2018.lumi}

import tWZ.samples.nanoTuples_Summer16_nanoAODv6_private_postProcessed as Summer16
import tWZ.samples.nanoTuples_Fall17_nanoAODv6_private_postProcessed   as Fall17
import tWZ.samples.nanoTuples_Autumn18_nanoAODv6_private_postProcessed as Autumn18

TWZ          = Sample.combine( "TWZ", [Summer16.TWZ, Fall17.TWZ, Autumn18.TWZ] )
TWZ_NLO_DR   = Sample.combine( "TWZ_NLO_DR", [Summer16.TWZ_NLO_DR, Fall17.TWZ_NLO_DR, Autumn18.TWZ_NLO_DR] )
TWZ_NLO_DS   = Sample.combine( "TWZ_NLO_DS", [Summer16.TWZ_NLO_DS, Fall17.TWZ_NLO_DS, Autumn18.TWZ_NLO_DS] )
TTZ          = Sample.combine( "TTZ", [Summer16.TTZ, Fall17.TTZ, Autumn18.TTZ] )
TTX_rare     = Sample.combine( "TTX_rare", [Summer16.TTX_rare, Fall17.TTX_rare, Autumn18.TTX_rare] )
TZQ          = Sample.combine( "TZQ", [Summer16.TZQ, Fall17.TZQ, Autumn18.TZQ] )
WZ           = Sample.combine( "WZ", [Summer16.WZ, Fall17.WZ, Autumn18.WZ] )
triBoson     = Sample.combine( "triBoson", [Summer16.triBoson, Fall17.triBoson, Autumn18.triBoson] )
ZZ           = Sample.combine( "ZZ", [Summer16.ZZ, Fall17.ZZ, Autumn18.ZZ] )
nonprompt_3l = Sample.combine( "nonprompt_3l", [Summer16.nonprompt_3l, Fall17.nonprompt_3l, Autumn18.nonprompt_3l] )
Xgamma       = Sample.combine( "Xgamma", [Summer16.TWZ, Fall17.TWZ, Autumn18.TWZ])
