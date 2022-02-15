# ZZ
python EFT_SYS.py  --noData --selection=trilepVL-minDLmass12-onZ1-onZ2-nLeptons4 --nicePlots

# ttZ
python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-njet3-deepjet1 --nicePlots
python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-njet3-deepjet1p --nicePlots

python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-njet3p-deepjet1 --nicePlots
python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-njet3p-deepjet1p --nicePlots
python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-njet3p-deepjet2p --nicePlots

python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-njet4p-deepjet1 --nicePlots
python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-njet4p-deepjet1p --nicePlots
python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-njet4p-deepjet2p --nicePlots

# WZ
python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-deepjet0-met60 --nicePlots
python EFT_SYS.py  --noData --selection=trilepT-minDLmass12-onZ1-deepjet0 --nicePlots

# WZ control regions
python EFT_SYS.py  --selection=trilepVetoT-minDLmass12-onZ1-deepjet0-met60 --nicePlots
python EFT_SYS.py  --selection=trilepVetoT-minDLmass12-offZ1-deepjet0-met60 --nicePlots
python EFT_SYS.py  --selection=trilepT-minDLmass12-offZ1-deepjet0-met60 --nicePlots

python EFT_SYS.py  --selection=trilepVetoT-minDLmass12-onZ1-deepjet0 --nicePlots
