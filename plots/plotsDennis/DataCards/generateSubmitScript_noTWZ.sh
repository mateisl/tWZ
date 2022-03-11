#!/bin/bash
npoints=50
dir=/users/dennis.schwarz/CMSSW_10_6_0/src/tWZ/plots/plotsDennis/DataCards/data_noTWZ/
submitdir=/users/dennis.schwarz/CMSSW_10_6_0/src/tWZ/plots/plotsDennis/DataCards/submitScripts

output=submitScript_noTWZ.sh

rm -f ${output}

# First combine regions to a final card and name region 'combined'
echo 'Combining combine cards...'
cd ${dir}
for WCname in cHq1Re11 cHq1Re22 cHq1Re33 cHq3Re11 cHq3Re22 cHq3Re33
do
  echo "   - ${WCname}"
  for i in {0..${npoints}}
  do
    combineCards.py topEFT_${WCname}_1_13TeV_${i}.txt topEFT_${WCname}_2_13TeV_${i}.txt topEFT_${WCname}_3_13TeV_${i}.txt topEFT_${WCname}_4_13TeV_${i}.txt > topEFT_${WCname}_combined_13TeV_${i}.txt
  done
done
cd ..

# Now write a python file for each region and WCname
cd ${submitdir}
echo "Write python files in ${submitdir}"
for WCname in cHq1Re11 cHq1Re22 cHq1Re33 cHq3Re11 cHq3Re22 cHq3Re33
do
  for region in 1 2 3 4 combined
  do
    pythonoutput=noTWZ__${WCname}__${region}.py
    echo "   - ${pythonoutput}"
    echo "python ${pythonoutput}" >> ${output}
    rm -f ${pythonoutput}
    echo "#!/usr/bin/env python" >> ${pythonoutput}
    echo "import os" >> ${pythonoutput}
    echo "os.chdir('${dir}')" >> ${pythonoutput}
    for i in {0..${npoints}}
    do
      command1="text2workspace.py topEFT_${WCname}_${region}_13TeV_${i}.txt -m ${i} -o topEFT_${WCname}_${region}_13TeV_${i}.root"
      command2="combine -M MultiDimFit topEFT_${WCname}_${region}_13TeV_${i}.root -n .part3E_${i}_${region}_${WCname} --saveNLL --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --expectSignal=1 --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.0001 --freezeParameters r  --setParameterRanges r=-10,10"
      echo "os.system('${command1}')" >> ${pythonoutput}
      echo "os.system('${command2}')" >> ${pythonoutput}
    done
  done
done
cd ..
