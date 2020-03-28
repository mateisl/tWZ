## DY
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --reduceSizeBy 2 --sample DYJetsToLL_M50_LO #SPLIT40
##python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50 #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --reduceSizeBy 2 --sample DYJetsToLL_M10to50_LO #SPLIT40

## full stats
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_LO #SPLIT40
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M10to50_LO #SPLIT29
#
## HT binned samples ##
## high mass #
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --LHEHTCut 70 --sample DYJetsToLL_M50_LO #SPLIT40
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT70to100 #SPLIT16
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT100to200 #SPLIT12
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT200to400 #SPLIT12
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT400to600 DYJetsToLL_M50_HT400to600_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT600to800 #SPLIT13
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT800to1200 #SPLIT6
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT1200to2500 #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT2500toInf #SPLIT2
# low mass #
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --LHEHTCut 70 --sample DYJetsToLL_M10to50_LO #SPLIT28
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M4to50_HT70to100 #SPLIT16
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M4to50_HT100to200 #SPLIT11
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M4to50_HT200to400 #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M4to50_HT400to600 #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M4to50_HT600toInf #SPLIT5
#
# top
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --reduceSizeBy 3 --sample TTLep_pow #SPLIT50
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --flagTTBar --sample TTLep_pow  #SPLIT59
#python nanoPostProcessing.py  --overwrite --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --flagTTBar --sample TTSingleLep_pow #SPLIT80
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --flagTTGamma --sample TTGLep #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --reduceSizeBy 5 --sample TToLeptons_sch_amcatnlo #SPLIT5
python nanoPostProcessing.py --overwrite --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample T_tWch #SPLIT3
python nanoPostProcessing.py --overwrite --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TBar_tWch #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --reduceSizeBy 15 --sample T_tch_pow #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --reduceSizeBy 15 --sample TBar_tch_pow #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample tZq_ll #SPLIT17
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample tWll 
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TTWToLNu #SPLIT7
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TTWToQQ #SPLIT3
##python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TTW_LO #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TTZToLLNuNu #SPLIT15
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TTZToLLNuNu_m1to10 #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TTZToQQ #SPLIT2
##python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TTZ_LO #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TTHbbLep #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TTHnobb_pow #SPLIT12
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample THQ #SPLIT4
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample THW #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TTTT #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample tWnunu 

##
##
## di/multi boson

python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample VVTo2L2Nu #SPLIT13
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample WZTo3LNu_amcatnlo #SPLIT16
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample WZTo2L2Q #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample WZTo1L3Nu #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample WWTo1L1Nu2Q #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample WWToLNuQQ #SPLIT18
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample ZZTo4L #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample ZZTo2Q2Nu #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --reduceSizeBy 5 --sample ZZTo2L2Q #SPLIT5
# additional samples for ARC studies
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample WWTo2L2Nu #SPLIT11
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample ZZTo2L2Nu #SPLIT11
##
###python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample WW #SPLIT5
###python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample ZZ #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample WWW_4F #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample WWZ 
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample WZZ #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample ZZZ 
##
### rare
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TTWZ 
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TTWW #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6 --sample TTZZ 
#
