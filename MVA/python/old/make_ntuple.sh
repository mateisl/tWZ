#!/bin/sh
python make_ntuple.py  --output /scratch/$USER/tWZ/ --sample TWZ_NLO_DR 
python make_ntuple.py  --output /scratch/$USER/tWZ/ --sample WZ 
python make_ntuple.py  --output /scratch/$USER/tWZ/ --sample TTZ 
python make_ntuple.py  --output /scratch/$USER/tWZ/ --sample nonprompt_3l
