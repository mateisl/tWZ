import os

if os.environ["USER"] in ["schoef", "rschoefbeck", "schoefbeck"]:

    postprocessing_output_directory = "/afs/hephy.at/data/cms06/nanoTuples/"
    plot_directory                  = "/afs/hephy.at/user/r/rschoefbeck/www/tWZ"
    cache_dir                       = "/afs/hephy.at/data/cms01/tWZ/signals/caches/"
    # Analysis result files
    analysis_results                = "/afs/hephy.at/data/cms05/tWZ/results/v1/"
    dpm_directory                   = "/dpm/oeaw.ac.at/home/cms/store/user/schoef/"
    cern_proxy_certificate          = "/afs/cern.ch/user/s/schoef/private/.proxy"
    # directory with veto lists

if os.environ["USER"] in ["robert.schoefbeck"]:
    #postprocessing_output_directory = "/mnt/hephy/cms/robert.schoefbeck/tWZ/nanoTuples"
    postprocessing_output_directory = "/scratch-cbe/users/robert.schoefbeck/tWZ/nanoTuples"
    postprocessing_tmp_directory    = "/scratch/hephy/cms/robert.schoefbeck/tWZ/tmp/"
    plot_directory                  = "/mnt/hephy/cms/robert.schoefbeck/www/tWZ/plots"
    cache_dir                       = "/groups/hephy/cms/robert.schoefbeck/tWZ/caches"
    # Analysis result files
    analysis_results                = "/groups/hephy/cms/robert.schoefbeck/tWZ/results/v1"
    mva_directory                   = "/groups/hephy/cms/robert.schoefbeck/tWZ/MVA"
    cern_proxy_certificate          = "/users/robert.schoefbeck/.private/.proxy"

if os.environ["USER"] in ["rosmarie.schoefbeck"]:
    #postprocessing_output_directory = "/mnt/hephy/cms/rosmarie.schoefbeck/tWZ/nanoTuples"
    postprocessing_output_directory = "/scratch-cbe/users/rosmarie.schoefbeck/tWZ/nanoTuples"
    postprocessing_tmp_directory    = "/scratch/hephy/cms/rosmarie.schoefbeck/tWZ/tmp/"
    plot_directory                  = "/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots"
    cache_dir                       = "/mnt/hephy/cms/rosmarie.schoefbeck/tWZ/caches"
    # Analysis result files
    analysis_results                = "/mnt/hephy/cms/rosmarie.schoefbeck/tWZ/results/v1"
    mva_directory                   = "/mnt/hephy/cms/rosmarie.schoefbeck/tWZ/MVA"
    mva_directory_twz_wz            = "/mnt/hephy/cms/rosmarie.schoefbeck/tWZ_WZ/MVA"
    mva_directory_wz                = "/mnt/hephy/cms/rosmarie.schoefbeck/WZ/MVA"
    cern_proxy_certificate          = "/users/rosmarie.schoefbeck/.private/.proxy"

if os.environ['USER'] in ['janik.andrejkovic']:
    postprocessing_output_directory = "/mnt/hephy/cms/janik.andrejkovic/tWZ/nanoTuples"
    postprocessing_tmp_directory    = "/scratch/hephy/cms/janik.andrejkovic/tWZ/tmp/"
    plot_directory                  = "/mnt/hephy/cms/janik.andrejkovic/www/tWZ/plots"
    cache_dir                       = "/mnt/hephy/cms/janik.andrejkovic/tWZ/caches"
    # Analysis result files
    analysis_results                = "/mnt/hephy/cms/janik.andrejkovic/tWZ/results/v1"
    mva_directory                   = "/mnt/hephy/cms/janik.andrejkovic/tWZ/MVA"
    cern_proxy_certificate          = "/users/janik.andrejkovic/.private/.proxy"
    HistogramsForFakeStudies        = "/users/janik.andrejkovic/public/test/CMSSW_10_2_18/src/tWZ/fake/helperHistograms"

if os.environ["USER"] in ["dennis.schwarz"]:
    postprocessing_output_directory = "/scratch-cbe/users/dennis.schwarz/tWZ/nanoTuples"
    postprocessing_tmp_directory    = "/scratch/hephy/cms/dennis.schwarz/tWZ/tmp/"
    plot_directory                  = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/plots"
    cache_dir                       = "/mnt/hephy/cms/dennis.schwarz/tWZ/caches"
    # Analysis result files
    analysis_results                = "/mnt/hephy/cms/dennis.schwarz/tWZ/results/v1"
    mva_directory                   = "/mnt/hephy/cms/dennis.schwarz/tWZ/MVA"
    cern_proxy_certificate          = "/users/dennis.schwarz/.private/.proxy"
    combineReleaseLocation          = "/users/dennis.schwarz/CMSSW_10_6_0/src"
    limit_directory                 = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/limits"
