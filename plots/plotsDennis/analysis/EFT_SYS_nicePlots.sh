# ZZ
python EFT_SYS.py  --noData --selection=trilepVL-minDLmass12-onZ1-onZ2-nLeptons4 --nicePlots

# ttZ 1
python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-njet3-nBTag1p --nicePlots
# ttZ 2
python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-njet4p-nBTag1 --nicePlots
# python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-njet4p-nBTag1p --nicePlots --doTTbarReco


# WZ
python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-nBTag0-met60 --nicePlots
# python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-nBTag0 --nicePlots

# WZ control regions
# python EFT_SYS.py  --selection=trilepVetoT-minDLmass12-onZ1-nBTag0-met60 --nicePlots
# python EFT_SYS.py  --selection=trilepVetoT-minDLmass12-offZ1-nBTag0-met60 --nicePlots
# python EFT_SYS.py  --selection=trilepT-minDLmass12-offZ1-nBTag0-met60 --nicePlots
# 
# python EFT_SYS.py  --selection=trilepVetoT-minDLmass12-onZ1-nBTag0 --nicePlots
