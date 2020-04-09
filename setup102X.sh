# 
eval `scram runtime -sh`
cd $CMSSW_BASE/src

# nanoAOD tools (for MET Significance, JEC/JER...)
git clone -b stopsDilepton https://github.com/HephyAnalysisSW/nanoAOD-tools.git PhysicsTools/NanoAODTools
cd $CMSSW_BASE/src

# RootTools (for plotting, sample handling, processing)
git clone https://github.com/HephyAnalysisSW/RootTools.git
cd $CMSSW_BASE/src

# Shared samples (miniAOD/nanoAOD)
git clone https://github.com/HephyAnalysisSW/Samples.git
cd $CMSSW_BASE/src

# Shared analysis tools and data
git clone https://github.com/HephyAnalysisSW/Analysis.git
cd $CMSSW_BASE/src

#scram b -j9

## Get combine
## Latest recommendations at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#setting-up-the-environment-and-installation
#cd $CMSSW_BASE/src
#git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
#cd HiggsAnalysis/CombinedLimit
#git fetch origin
#git checkout v8.0.1
#scramv1 b clean; scramv1 b # always make a clean build
#
## for combineTools
#cd $CMSSW_BASE/src
#wget https://raw.githubusercontent.com/cms-analysis/CombineHarvester/master/CombineTools/scripts/sparse-checkout-https.sh; source sparse-checkout-https.sh
#scram b -j 8

cd $CMSSW_BASE/src/tWZ
#git fetch origin
#git checkout -b 101X origin/101X

#compile
cd $CMSSW_BASE/src && scram b -j 8 
