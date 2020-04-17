###
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TTZToQQ #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TTZToLLNuNu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TTZToLLNuNu_m1to10 #SPLIT4
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --flagTTBar --sample TTLep_pow #SPLIT77
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --flagTTBar --sample TTSingleLep_pow #SPLIT80
##
#python nanoPostProcessing.py/  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M50_ext1 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M50_LO DYJetsToLL_M50_LO_ext1 #SPLIT50
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M10to50_LO #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M50_LO DYJetsToLL_M50_LO_ext1 --LHEHTCut 70 #SPLIT32
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M50_HT70to100 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M50_HT100to200 #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M50_HT200to400 DYJetsToLL_M50_HT200to400_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M50_HT400to600 DYJetsToLL_M50_HT400to600_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M50_HT600to800 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M50_HT800to1200 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M50_HT1200to2500 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M50_HT2500toInf #SPLIT10
#
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M10to50_LO --LHEHTCut 100 --overwrite #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M4to50_HT100to200 DYJetsToLL_M4to50_HT100to200_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M4to50_HT200to400 DYJetsToLL_M4to50_HT200to400_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M4to50_HT400to600 DYJetsToLL_M4to50_HT400to600_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample DYJetsToLL_M4to50_HT600toInf DYJetsToLL_M4to50_HT600toInf_ext #SPLIT10

python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TToLeptons_sch_amcatnlo #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample T_tch_pow #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TBar_tch_pow #SPLIT20
python nanoPostProcessing.py --overwrite --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample T_tWch_ext #SPLIT20
python nanoPostProcessing.py --overwrite --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TBar_tWch_ext #SPLIT20

python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TTHbb #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TTHnobb_pow #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample THQ #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample THW #SPLIT19
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TGJets #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample tZq_ll #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample tZq_nunu #SPLIT15
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample tWll #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample tWnunu #SPLIT6
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TTWW #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TTWZ #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TTZZ #SPLIT1
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TTWToLNu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TTWToQQ #SPLIT9
python nanoPostProcessing.py   --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TTW_LO #SPLIT20
python nanoPostProcessing.py   --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TTGJets #SPLIT20
#
#
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample VVTo2L2Nu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample WWToLNuQQ #SPLIT20
#python nanoPostProcessing.py   --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample ZZTo2Q2Nu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample ZZTo2L2Q #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample ZZTo4L #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample WZTo1L1Nu2Q #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample WZTo2L2Q #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample WZTo1L3Nu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample WZTo3LNu_amcatnlo #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample WZTo3LNu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample TTTT #SPLIT12
#
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample WWTo2L2Nu #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample ZZTo2L2Nu #SPLIT20

python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample WWW_4F #SPLIT11
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample WWZ_4F #SPLIT12
#python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample WZG #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample WZZ #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample ZZZ #SPLIT11

python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample WW #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample WZ #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4_v1 --sample ZZ #SPLIT5
