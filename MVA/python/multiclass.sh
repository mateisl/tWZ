# python $CMSSW_BASE/src/Analysis/MVA/python/multiclass.py --name tWZ_3l_ttz --input_directory /scratch-cbe/users/$USER/tWZ/training-ntuples-tWZ-v1/MVA-training --config_module tWZ.MVA.configs --config tWZ_3l --output_directory /mnt/hephy/cms/dennis.schwarz/tWZ/

# python $CMSSW_BASE/src/Analysis/MVA/python/multiclass.py --name tWZ_3l_ttz_topReco --input_directory /scratch-cbe/users/$USER/tWZ/training-ntuples-tWZ-topReco-v1/MVA-training --config_module tWZ.MVA.configs --config tWZ_3l_topReco --output_directory /mnt/hephy/cms/dennis.schwarz/tWZ/

# python $CMSSW_BASE/src/Analysis/MVA/python/multiclass_generator.py --add_LSTM --name tWZ_3l_ttz_topReco --input_directory /scratch-cbe/users/$USER/tWZ/training-ntuples-tWZ-topReco-v1/MVA-training --config_module tWZ.MVA.configs --config tWZ_3l_topReco --output_directory /mnt/hephy/cms/$USER/tWZ/

# python $CMSSW_BASE/src/Analysis/MVA/python/multiclass.py --name tWZ_3l_ttz_deltaEta --input_directory /scratch-cbe/users/$USER/tWZ/training-ntuples-tWZ-deltaEta-v1/MVA-training --config_module tWZ.MVA.configs --config tWZ_3l_deltaEta --output_directory /mnt/hephy/cms/dennis.schwarz/tWZ/

python $CMSSW_BASE/src/Analysis/MVA/python/multiclass.py --name tWZ_3l_ttz_LSTM --add_LSTM --input_directory /scratch-cbe/users/$USER/tWZ/training-ntuples-tWZ-v1/MVA-training --config_module tWZ.MVA.configs --config tWZ_3l --output_directory /mnt/hephy/cms/dennis.schwarz/tWZ/


# python $CMSSW_BASE/src/Analysis/MVA/python/multiclass.py --name tWZ_3l_ttz_dummy --input_directory /scratch-cbe/users/$USER/tWZ/training-ntuples-tWZ-dummy-v1/MVA-training --config_module tWZ.MVA.configs --config tWZ_3l_dummy --output_directory /mnt/hephy/cms/dennis.schwarz/tWZ/
