#!/bin/bash

combineTool.py -M Impacts -d $1.root -m 125 --doInitialFit --robustFit 1 --rMax 10 --rMin -10
combineTool.py -M Impacts -d $1.root -m 125 --robustFit 1 --doFits --rMax 10 --rMin -10
combineTool.py -M Impacts -d $1.root -m 125 -o $1.json
plotImpacts.py -i  $1.json -o $1
