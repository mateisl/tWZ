## DY
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --reduceSizeBy 2 --sample DYJetsToLL_M50_LO #SPLIT40
##python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M50 #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --reduceSizeBy 2 --sample DYJetsToLL_M10to50_LO #SPLIT40

## full stats
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M50_LO #SPLIT40
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M10to50_LO #SPLIT40
#
## HT binned samples ##
## high mass #
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --LHEHTCut 70 --sample DYJetsToLL_M50_LO #SPLIT40
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M50_HT70to100 #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M50_HT100to200 #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M50_HT200to400 #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M50_HT400to600 DYJetsToLL_M50_HT400to600_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M50_HT600to800 #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M50_HT800to1200 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M50_HT1200to2500 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M50_HT2500toInf #SPLIT10
# low mass #
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --LHEHTCut 70 --sample DYJetsToLL_M10to50_LO #SPLIT40
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M4to50_HT70to100 #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M4to50_HT100to200 #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M4to50_HT200to400 #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M4to50_HT400to600 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample DYJetsToLL_M4to50_HT600toInf #SPLIT10
#
# top
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --reduceSizeBy 3 --sample TTLep_pow #SPLIT80
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --flagTTBar --sample TTLep_pow  #SPLIT80
#python nanoPostProcessing.py  --overwrite --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --flagTTBar --sample TTSingleLep_pow #SPLIT80
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --flagTTGamma --sample TTGLep #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --reduceSizeBy 5 --sample TToLeptons_sch_amcatnlo #SPLIT20
python nanoPostProcessing.py --overwrite --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample T_tWch #SPLIT20
python nanoPostProcessing.py --overwrite --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TBar_tWch #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --reduceSizeBy 15 --sample T_tch_pow #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --reduceSizeBy 15 --sample TBar_tch_pow #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample tZq_ll #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample tWll #SPLIT11
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TTWToLNu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TTWToQQ #SPLIT14
##python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TTW_LO #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TTZToLLNuNu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TTZToLLNuNu_m1to10 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TTZToQQ #SPLIT20
##python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TTZ_LO #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TTHbbLep #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TTHnobb_pow #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample THQ #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample THW #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TTTT #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample tWnunu #SPLIT5

##
##
## di/multi boson



python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample VVTo2L2Nu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample WZTo3LNu_amcatnlo #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample WZTo2L2Q #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample WZTo1L3Nu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample WWTo1L1Nu2Q #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample WWToLNuQQ #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample ZZTo4L #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample ZZTo2Q2Nu #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --reduceSizeBy 5 --sample ZZTo2L2Q #SPLIT20
# additional samples for ARC studies
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample WWTo2L2Nu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample ZZTo2L2Nu #SPLIT20
##
###python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample WW #SPLIT5
###python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample ZZ #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample WWW_4F #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample WWZ #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample WZZ #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample ZZZ #SPLIT5
##
### rare
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TTWZ #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TTWW #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample TTZZ #SPLIT5
#
#python nanoPostProcessing.py  --forceProxy --skim trilep  --year 2018 --processingEra tWZ_v3 --susySignal --fastSim --sample SMS_T2tt_mStop_150to250 #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim trilep  --year 2018 --processingEra tWZ_v3 --susySignal --fastSim --sample SMS_T2tt_mStop_250to350 #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim trilep  --year 2018 --processingEra tWZ_v3 --susySignal --fastSim --sample SMS_T2tt_mStop_350to400 #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim trilep  --year 2018 --processingEra tWZ_v3 --susySignal --fastSim --sample SMS_T2tt_mStop_400to1200 #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim trilep  --year 2018 --processingEra tWZ_v3 --susySignal --fastSim --sample SMS_T2tt_mStop_1200to2000 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim trilep  --year 2018 --processingEra tWZ_v3 --susySignal --fastSim --sample SMS_T2bW #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim trilep  --year 2018 --processingEra tWZ_v3 --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p05 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim trilep  --year 2018 --processingEra tWZ_v3 --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p05_mN1_700_1000 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim trilep  --year 2018 --processingEra tWZ_v3 --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p5 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim trilep  --year 2018 --processingEra tWZ_v3 --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p5_mN1_700_1300 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim trilep  --year 2018 --processingEra tWZ_v3 --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p95_mN1_0_650 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim trilep  --year 2018 --processingEra tWZ_v3 --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p95_mN1_700_1600 #SPLIT45
#
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2018 --processingEra tWZ_v3 --sample ttH_HToInvisible # SPLIT20

#python nanoPostProcessing.py  --forceProxy --skim inclusive  --year 2018 --processingEra tWZ_v3 --susySignal --sample SMS_T2tt_mStop_175_mLSP_1 #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 175,0  --year 2018 --processingEra stops_2018_nano_v3p20_1 --susySignal --fastSim --sample SMS_T2tt_mStop_150to250
#python nanoPostProcessing.py  --forceProxy --skim inclusive  --year 2018 --processingEra tWZ_v3 --susySignal --sample SMS_T2tt_mStop_250_mLSP_50 #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 250,50  --year 2018 --processingEra stops_2018_nano_v3p20_1 --susySignal --fastSim --sample SMS_T2tt_mStop_150to250 #SPLIT1
#python nanoPostProcessing.py  --forceProxy --skim inclusive  --year 2018 --processingEra tWZ_v3 --susySignal --sample SMS_T2tt_mStop_250_mLSP_75 #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 250,75  --year 2018 --processingEra stops_2018_nano_v3p20_1 --susySignal --fastSim --sample SMS_T2tt_mStop_150to250 #SPLIT1
#python nanoPostProcessing.py  --forceProxy --skim inclusive  --year 2018 --processingEra tWZ_v3 --susySignal --sample SMS_T2tt_mStop_250_mLSP_100 #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 250,100  --year 2018 --processingEra stops_2018_nano_v3p20_1 --susySignal --fastSim --sample SMS_T2tt_mStop_150to250 #SPLIT1
#
#python nanoPostProcessing.py  --forceProxy --skim trilep  --year 2018 --processingEra tWZ_v3 --susySignal --sample SMS_T2tt_mStop_650_mLSP_350 #SPLIT5
