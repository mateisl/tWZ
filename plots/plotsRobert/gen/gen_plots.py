''' FWLiteReader example: Loop over a sample and write some data to a histogram.
'''
# Standard imports
import os
import logging
import ROOT

#RootTools
from RootTools.core.standard import *

# argParser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel', 
      action='store',
      nargs='?',
      choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],
      default='INFO',
      help="Log level for logging"
)

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)


import tWZ.samples.GENandSIM_private as samples
## from DAS
s = tOrtbar_WZ01j_OLRLL_LO

s.reduceFiles( to = 1 )

products = {
#    'slimmedJets':{'type':'vector<pat::Jet>', 'label':("slimmedJets", "", "reRECO")} 
      'genJets': {'type':'vector<reco::GenJet>', 'label':( "ak4GenJets" ) } ,
    }

r = s.fwliteReader( products = products )

r.start()

while r.run():
    print r.event.evt, r.event.lumi, r.event.run, "Number of genJets", r.event.genJets.size()
    break


