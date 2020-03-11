###
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TTZToQQ #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TTZToLLNuNu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TTZToLLNuNu_m1to10 #SPLIT4
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --flagTTBar --sample TTLep_pow #SPLIT77
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --flagTTBar --sample TTSingleLep_pow #SPLIT80
##
#python nanoPostProcessing.py/  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M50_ext1 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M50_LO DYJetsToLL_M50_LO_ext1 #SPLIT50
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M10to50_LO #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M50_LO DYJetsToLL_M50_LO_ext1 --LHEHTCut 70 #SPLIT32
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M50_HT70to100 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M50_HT100to200 #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M50_HT200to400 DYJetsToLL_M50_HT200to400_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M50_HT400to600 DYJetsToLL_M50_HT400to600_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M50_HT600to800 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M50_HT800to1200 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M50_HT1200to2500 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M50_HT2500toInf #SPLIT10
#
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M10to50_LO --LHEHTCut 100 --overwrite #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M4to50_HT100to200 DYJetsToLL_M4to50_HT100to200_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M4to50_HT200to400 DYJetsToLL_M4to50_HT200to400_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M4to50_HT400to600 DYJetsToLL_M4to50_HT400to600_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample DYJetsToLL_M4to50_HT600toInf DYJetsToLL_M4to50_HT600toInf_ext #SPLIT10

python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TToLeptons_sch_amcatnlo #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample T_tch_pow #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TBar_tch_pow #SPLIT20
python nanoPostProcessing.py --overwrite --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample T_tWch_ext #SPLIT20
python nanoPostProcessing.py --overwrite --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TBar_tWch_ext #SPLIT20

python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TTHbb #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TTHnobb_pow #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample THQ #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample THW #SPLIT19
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TGJets #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample tZq_ll #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample tZq_nunu #SPLIT15
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample tWll #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample tWnunu #SPLIT6
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TTWW #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TTWZ #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TTZZ #SPLIT1
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TTWToLNu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TTWToQQ #SPLIT9
python nanoPostProcessing.py   --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TTW_LO #SPLIT20
python nanoPostProcessing.py   --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TTGJets #SPLIT20
#
#
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample VVTo2L2Nu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample WWToLNuQQ #SPLIT20
#python nanoPostProcessing.py   --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample ZZTo2Q2Nu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample ZZTo2L2Q #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample ZZTo4L #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample WZTo1L1Nu2Q #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample WZTo2L2Q #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample WZTo1L3Nu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample WZTo3LNu_amcatnlo #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample TTTT #SPLIT12
#
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample WWTo2L2Nu #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample ZZTo2L2Nu #SPLIT20

python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample WWW_4F #SPLIT11
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample WWZ_4F #SPLIT12
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample WZG #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample WZZ #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample ZZZ #SPLIT11

python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample WW #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample WZ #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample ZZ #SPLIT5

##python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23   --skipGenLepMatching --sample ttH_HToInvisible #SPLIT20
#
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T2tt_mStop_150to250 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T2tt_mStop_250to350 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T2tt_mStop_350to400 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T2tt_mStop_400to1200 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T2tt_mStop_1200to2000 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T2bW #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p05 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p05_mN1_700_1000 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p5 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p5_mN1_700_1300 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p95 #SPLIT45
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T8bbllnunu_XCha0p5_XSlep0p95_mN1_700_1600 #SPLIT45

## FullSim corridor
#python nanoPostProcessing.py  --forceProxy --skim inclusive  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --sample SMS_T2tt_mStop_175_mLSP_1
#python nanoPostProcessing.py  --forceProxy --skim inclusive  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --sample SMS_T2tt_mStop_250_mLSP_50
#python nanoPostProcessing.py  --forceProxy --skim inclusive  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --sample SMS_T2tt_mStop_250_mLSP_75
#python nanoPostProcessing.py  --forceProxy --skim inclusive  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --sample SMS_T2tt_mStop_250_mLSP_100
#python nanoPostProcessing.py  --forceProxy --skim inclusive  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --sample SMS_T2tt_mStop_650_mLSP_350
#python nanoPostProcessing.py  --forceProxy --skim inclusive  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --sample SMS_T2tt_mStop_850_mLSP_100
#
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 175,1 --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T2tt_mStop_150to250
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 250,50 --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T2tt_mStop_150to250
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 250,75 --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T2tt_mStop_150to250
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 250,100 --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T2tt_mStop_150to250
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 650,350 --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T2tt_mStop_400to1200
#python nanoPostProcessing.py  --forceProxy --skim inclusive --massSkim 850,100 --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T2tt_mStop_400to1200

#2017 T8bbstausnu
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T8bbstausnu_mStop_200to1800_XCha0p5_XStau0p5 #SPLIT30
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T8bbstausnu_mStop_200to1800_XCha0p5_XStau0p5 #SPLIT30
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T8bbstausnu_XCha0p5_mStop_200to1800_XStau0p25 #SPLIT30
#
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T8bbstausnu_mStop_200to1800_XCha0p5_XStau0p5 #SPLIT30
#python nanoPostProcessing.py  --forceProxy --skim dilep  --year 2017 --processingEra stops_2017_nano_v0p23 --skipGenLepMatching --susySignal --fastSim --sample SMS_T8bbstausnu_mStop_200to1800_XCha0p5_XStau0p75 #SPLIT30
#
