#!/bin/bash
npoints1=20 # put one less here because of bash loop logic
npoints2=20

echo 'Running over' ${npoints1}+1 'and' ${npoints2}+1 'points'

cd /users/dennis.schwarz/CMSSW_10_6_0/src/tWZ/plots/plotsDennis/DataCards/data_twoD

WC1=cHq1Re11
WC2=cHq3Re11

for i in {0..${npoints1}}
do
  # First combine regions to a final card and name region 'combined'
  for j in {0..${npoints2}}
  do
    combineCards.py topEFT_2D_${WC1}_${WC2}_1_13TeV_${i}_${j}.txt topEFT_2D_${WC1}_${WC2}_2_13TeV_${i}_${j}.txt topEFT_2D_${WC1}_${WC2}_3_13TeV_${i}_${j}.txt > topEFT_2D_${WC1}_${WC2}_combined_13TeV_${i}_${j}.txt
  done
done
# Now Run over all txt files, convert to root and run combine
for region in 1 2 3 combined
do
  for i in {0..${npoints1}}
  do
    for j in {0..${npoints2}}
    do
      text2workspace.py topEFT_2D_${WC1}_${WC2}_${region}_13TeV_${i}_${j}.txt -m ${i}${j} -o topEFT_2D_${WC1}_${WC2}_${region}_13TeV_${i}_${j}.root
      combine -M MultiDimFit topEFT_2D_${WC1}_${WC2}_${region}_13TeV_${i}_${j}.root -n .part3E_${i}_${j}_${region}_${WC1}_${WC2} --saveNLL --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --expectSignal=1 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.0001 --freezeParameters r  --setParameterRanges r=-10,10
      echo '--------------------------------------------------------------------'
      echo '     done with' ${WC1}_${WC2} ', region' ${region}, 'point' ${i} ${j}
      echo '--------------------------------------------------------------------'
    done
  done
done

cd ..
