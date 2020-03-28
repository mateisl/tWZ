
#python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4 --year 2016 --processingEra tWZ_nAODv4 --sample TToLeptons_sch_amcatnlo --nJobs 8 --job 2
#python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4 --year 2016 --processingEra tWZ_nAODv4 --sample tZq_nunu --nJobs 10 --job 9
#python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4 --year 2016 --processingEra tWZ_nAODv4 --sample WZTo1L3Nu --nJobs 6 --job 4
#
#
#python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TTZToLLNuNu_m1to10 --nJobs 4 --job 1
#python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TToLeptons_sch_amcatnlo --nJobs 20 --job 5
#python nanoPostProcessing.py --overwrite --forceProxy --skim trilep --nanoAODv4  --year 2017 --processingEra tWZ_nAODv4 --sample TBar_tWch_ext --nJobs 20 --job 5

python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4 --triggerSelection --year 2018 --processingEra tWZ_nAODv4  --overwrite --sample EGamma_Run2018D_17Sep2018 #SPLIT30
python nanoPostProcessing.py  --forceProxy --skim trilep --nanoAODv4 --triggerSelection --year 2018 --processingEra tWZ_nAODv4  --overwrite --sample SingleMuon_Run2018D_17Sep2018 #SPLIT30


