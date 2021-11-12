#!/bin/bash
npoints=50
echo 'Running over' ${npoints}+1 'points'

cd /users/dennis.schwarz/CMSSW_10_6_0/src/tWZ/plots/plotsDennis/DataCards/data_statonly

for WCname in cHq1Re11 cHq1Re22 cHq1Re33
do
  # First combine regions to a final card and name region 'combined'
  for i in {0..${npoints}}
  do
    combineCards.py topEFT_${WCname}_1_13TeV_${i}.txt topEFT_${WCname}_2_13TeV_${i}.txt topEFT_${WCname}_3_13TeV_${i}.txt > topEFT_${WCname}_combined_13TeV_${i}.txt
  done
  # Now Run over all txt files, convert to root and run combine
  for region in 1 2 3 combined
  do
    for i in {0..${npoints}}
    do
      text2workspace.py topEFT_${WCname}_${region}_13TeV_${i}.txt -m ${i} -o topEFT_${WCname}_${region}_13TeV_${i}.root
      combine -M MultiDimFit topEFT_${WCname}_${region}_13TeV_${i}.root -n .part3E_${i}_${region}_${WCname} --saveNLL --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --expectSignal=1 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.0001 --freezeParameters r  --setParameterRanges r=-10,10
      echo '--------------------------------------------------------------------'
      echo '     done with' ${WCname} ', region' ${region}, 'point' ${i}
      echo '--------------------------------------------------------------------'
    done
  done
done

cd ..
