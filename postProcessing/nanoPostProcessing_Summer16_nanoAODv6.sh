python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --flagTTBar --sample TTLep_pow #SPLIT50
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --flagTTBar --sample TTSingleLep_pow #SPLIT50
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --flagTTBar --sample TTGJets TTGJets_ext #SPLIT50
#
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TTZToQQ 
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TTZToLLNuNu_ext2 TTZToLLNuNu_ext3 #SPLIT9
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TTZToLLNuNu_m1to10 

python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_ext2 #SPLIT40
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M10to50 #SPLIT30

python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_LO_ext1 DYJetsToLL_M50_LO_ext2 #SPLIT50
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_LO_ext1 DYJetsToLL_M50_LO_ext2 --LHEHTCut 70 #SPLIT50
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT70to100 #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT100to200_ext #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT200to400 DYJetsToLL_M50_HT200to400_ext #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT400to600 DYJetsToLL_M50_HT400to600_ext #SPLIT11
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT600to800 #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT800to1200 #SPLIT4
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT1200to2500 
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M50_HT2500toInf 

python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M10to50_LO #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M10to50_LO --LHEHTCut 70 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M5to50_HT70to100 #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M5to50_HT100to200 DYJetsToLL_M5to50_HT100to200_ext #SPLIT9
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M5to50_HT200to400 DYJetsToLL_M5to50_HT200to400_ext #SPLIT4
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M5to50_HT400to600 DYJetsToLL_M5to50_HT400to600_ext #SPLIT4
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample DYJetsToLL_M5to50_HT600toInf #SPLIT2

python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TToLeptons_sch_amcatnlo #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample T_tch_pow #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TBar_tch_pow #SPLIT20
python nanoPostProcessing.py --overwrite --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample T_tWch_ext #SPLIT6
python nanoPostProcessing.py --overwrite --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TBar_tWch_ext #SPLIT5

python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TTHbb #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TTHnobb_pow #SPLIT7
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample THQ #SPLIT6
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample THW #SPLIT2
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TGJets TGJets_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample tZq_ll_ext #SPLIT15
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample tZq_nunu #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample tWll #SPLIT1
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample tWnunu #SPLIT1
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TTWToLNu_ext2 #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TTWToQQ #SPLIT2
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TTGJets TTGJets_ext #SPLIT20
#
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample VVTo2L2Nu VVTo2L2Nu_ext #SPLIT6
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample WWToLNuQQ #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample ZZTo2Q2Nu #SPLIT18
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample ZZTo2L2Q 
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample ZZTo2L2Nu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample ZZTo4L #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample WZTo1L1Nu2Q #SPLIT18
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample WZTo1L3Nu #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample WZTo2L2Q #SPLIT18
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample WZTo3LNu_ext #SPLIT14
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample WWTo2L2Nu #SPLIT2

python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TTTT #SPLIT1
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TTWW 
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TTWZ #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample TTZZ 

python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample WWW_4F #SPLIT1
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample WWZ #SPLIT1
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample WZG #SPLIT1
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample WZZ 
python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample ZZZ 

##python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample T2tt_mStop_850_mLSP_100 #SPLIT2
##python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --sample T2tt_mStop_500_mLSP_325 #SPLIT2
##
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T2tt_mStop_150to250 #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T2tt_mStop_250to350 #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T2tt_mStop_350to400 #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T2tt_mStop_400to1200 #SPLIT40
#
#python nanoPostProcessing.py  --forceProxy --overwrite --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T2bW #SPLIT30
#
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p05 #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p05_mN1_700_1000 #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p5 #SPLIT32
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p5_mN1_700_1300 #SPLIT32
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p95 #SPLIT32
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p95_mN1_700_1600 #SPLIT32
#
#python nanoPostProcessing.py  --forceProxy --skim inclusive --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_175_mLSP_1 #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 175,0 --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T2tt_mStop_150to250 #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim inclusive --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_250_mLSP_50 #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 250,50 --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T2tt_mStop_150to250 #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim inclusive --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_250_mLSP_75 #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 250,75 --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T2tt_mStop_150to250 #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim inclusive --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_250_mLSP_100 #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 250,100 --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T2tt_mStop_150to250 #SPLIT5
#
#
## FullSim corridor
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_150_mLSP_50 #SPLIT1
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_175_mLSP_1 #SPLIT1
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_200_mLSP_50 #SPLIT1
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_225_mLSP_50 #SPLIT1
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_250_mLSP_50 #SPLIT1
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_250_mLSP_75 #SPLIT1
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_250_mLSP_150 #SPLIT1
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_300_mLSP_150 #SPLIT1
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_325_mLSP_150 #SPLIT1
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_350_mLSP_150 #SPLIT1
#
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_650_mLSP_350 #SPLIT1
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --sample SMS_T2tt_mStop_850_mLSP_100 #SPLIT1
#
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T8bbstausnu_XCha0p5_mStop_200to1800_XStau0p25 #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T8bbstausnu_mStop_200to1800_XCha0p5_XStau0p5 #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim trilep --year 2016 --processingEra tWZ_nAODv6 --susySignal --fastSim --sample SMS_T8bbstausnu_mStop_200to1800_XCha0p5_XStau0p75 #SPLIT5
