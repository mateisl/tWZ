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
    postprocessing_output_directory = "/mnt/hephy/cms/robert.schoefbeck/tWZ/nanoTuples"
    postprocessing_tmp_directory    = "/scratch/hephy/cms/robert.schoefbeck/tWZ/tmp/"
    plot_directory                  = "/mnt/hephy/cms/robert.schoefbeck/www/tWZ/plots"
    cache_dir                       = "/mnt/hephy/cms/robert.schoefbeck/tWZ/caches"
    # Analysis result files
    analysis_results                = "/mnt/hephy/cms/robert.schoefbeck/tWZ/results/v1"
    cern_proxy_certificate          = "/users/robert.schoefbeck/.private/.proxy"
