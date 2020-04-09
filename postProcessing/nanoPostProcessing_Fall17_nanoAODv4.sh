###
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTZToQQ #SPLIT8
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTZToLLNuNu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTZToLLNuNu_m1to10 #SPLIT4
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --flagTTBar --sample TTLep_pow #SPLIT77
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --flagTTBar --sample TTSingleLep_pow #SPLIT80
##
#python nanoPostProcessing.py/  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M50_ext1 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M50_LO DYJetsToLL_M50_LO_ext1 #SPLIT50
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M10to50_LO #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M50_LO DYJetsToLL_M50_LO_ext1 --LHEHTCut 70 #SPLIT32
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M50_HT70to100 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M50_HT100to200 #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M50_HT200to400 DYJetsToLL_M50_HT200to400_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M50_HT400to600 DYJetsToLL_M50_HT400to600_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M50_HT600to800 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M50_HT800to1200 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M50_HT1200to2500 #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M50_HT2500toInf #SPLIT10
#
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M10to50_LO --LHEHTCut 100 --overwrite #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M4to50_HT100to200 DYJetsToLL_M4to50_HT100to200_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M4to50_HT200to400 DYJetsToLL_M4to50_HT200to400_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M4to50_HT400to600 DYJetsToLL_M4to50_HT400to600_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample DYJetsToLL_M4to50_HT600toInf DYJetsToLL_M4to50_HT600toInf_ext #SPLIT10

python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TToLeptons_sch_amcatnlo #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample T_tch_pow #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TBar_tch_pow #SPLIT20
python nanoPostProcessing.py --overwrite --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample T_tWch_ext #SPLIT20
python nanoPostProcessing.py --overwrite --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TBar_tWch_ext #SPLIT20

python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTHbb #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTHnobb_pow #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample THQ #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample THW #SPLIT19
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TGJets #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample tZq_ll #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample tZq_nunu #SPLIT15
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample tWll #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample tWnunu #SPLIT6
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTWW #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTWZ #SPLIT3
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTZZ #SPLIT1
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTWToLNu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTWToQQ #SPLIT9
python nanoPostProcessing.py   --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTW_LO #SPLIT20
python nanoPostProcessing.py   --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTGJets #SPLIT20
#
#
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample VVTo2L2Nu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample WWToLNuQQ #SPLIT20
#python nanoPostProcessing.py   --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample ZZTo2Q2Nu #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample ZZTo2L2Q #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample ZZTo4L #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample WZTo1L1Nu2Q #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample WZTo2L2Q #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample WZTo1L3Nu #SPLIT20
#python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample WZTo3LNu_amcatnlo #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample WZTo3LNu_ext #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTTT #SPLIT12
#
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample WWTo2L2Nu #SPLIT10
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample ZZTo2L2Nu #SPLIT20

python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample WWW_4F #SPLIT11
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample WWZ_4F #SPLIT12
#python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample WZG #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample WZZ #SPLIT2
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample ZZZ #SPLIT11

python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample WW #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample WZ #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample ZZ #SPLIT5
