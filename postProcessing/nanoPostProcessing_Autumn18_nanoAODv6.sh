## DY
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --reduceSizeBy 2 --sample DYJetsToLL_M50_LO #SPLIT40
##python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M50 #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --reduceSizeBy 2 --sample DYJetsToLL_M10to50_LO #SPLIT40

## full stats
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M50_LO #SPLIT40
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M10to50_LO #SPLIT29
#
## HT binned samples ##
## high mass #
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --LHEHTCut 70 --sample DYJetsToLL_M50_LO #SPLIT40
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M50_HT70to100 #SPLIT16
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M50_HT100to200 #SPLIT12
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M50_HT200to400 #SPLIT12
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M50_HT400to600 DYJetsToLL_M50_HT400to600_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M50_HT600to800 #SPLIT13
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M50_HT800to1200 #SPLIT6
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M50_HT1200to2500 #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M50_HT2500toInf #SPLIT2
# low mass #
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --LHEHTCut 70 --sample DYJetsToLL_M10to50_LO #SPLIT28
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M4to50_HT70to100 #SPLIT16
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M4to50_HT100to200 #SPLIT11
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M4to50_HT200to400 #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M4to50_HT400to600 #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample DYJetsToLL_M4to50_HT600toInf #SPLIT5
#
# top
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --reduceSizeBy 3 --sample TTLep_pow #SPLIT50
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --flagTTBar --sample TTLep_pow  #SPLIT59
#python nanoPostProcessing.py  --overwrite --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --flagTTBar --sample TTSingleLep_pow #SPLIT80
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --flagTTGamma --sample TTGLep #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --reduceSizeBy 5 --sample TToLeptons_sch_amcatnlo #SPLIT5
python nanoPostProcessing.py --overwrite --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample T_tWch #SPLIT3
python nanoPostProcessing.py --overwrite --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TBar_tWch #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --reduceSizeBy 15 --sample T_tch_pow #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --reduceSizeBy 15 --sample TBar_tch_pow #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample tZq_ll #SPLIT17
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample tWll 
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TTWToLNu #SPLIT7
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TTWToQQ #SPLIT3
##python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TTW_LO #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TTZToLLNuNu #SPLIT15
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TTZToLLNuNu_m1to10 #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TTZToQQ #SPLIT2
##python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TTZ_LO #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TTHbbLep #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TTHnobb_pow #SPLIT12
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample THQ #SPLIT4
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample THW #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TTTT #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample tWnunu 

##
##
## di/multi boson

python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample VVTo2L2Nu #SPLIT13
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample WZTo3LNu_amcatnlo #SPLIT16
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample WZTo2L2Q #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample WZTo1L3Nu #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample WWTo1L1Nu2Q #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample WWToLNuQQ #SPLIT18
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample ZZTo4L #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample ZZTo2Q2Nu #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --reduceSizeBy 5 --sample ZZTo2L2Q #SPLIT5
# additional samples for ARC studies
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample WWTo2L2Nu #SPLIT11
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample ZZTo2L2Nu #SPLIT11
##
###python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample WW #SPLIT5
###python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample ZZ #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample WWW_4F #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample WWZ 
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample WZZ #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample ZZZ 
##
### rare
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TTWZ 
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TTWW #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample TTZZ 
#
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample tWll_thad_Wlept_DR #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample tWll_thad_Wlept_DS #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample tWll_tlept_Whad_DR #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample tWll_tlept_Whad_DS #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample tWll_tlept_Wlept_DR #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v2 --sample tWll_tlept_Wlept_DS #SPLIT10
