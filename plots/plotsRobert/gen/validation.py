''' FWLiteReader example: Loop over a sample and write some data to a histogram.
'''
# Standard imports
import os
import logging
import ROOT
import array

#RootTools
from RootTools.core.standard import *

#Helper
from math import pi

# argParser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel', action='store', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], default='INFO', help="Log level for logging")

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

max_events  = 10000
max_files   = 10

LO  = FWLiteSample.fromDAS("LO",  "/tOrtbar_WZ01j_OLRLL_LO_ext5/schoef-tOrtbar_WZ01j_OLRLL_LO_ext5-c88179374d61ba4483b39a09a725f437/USER", instance="phys03", maxN=max_files)
NLO = FWLiteSample.fromDAS("NLO", "/tWZ_NLO_v16_ext3/schoef-tWZ_NLO_v16_ext3-a1e791adef4fc0a2f32451f12fdbd583/USER", instance="phys03", maxN=max_files)

LO.color = ROOT.kBlue
NLO.color= ROOT.kRed

samples = [LO, NLO]

def jetID(j):
    if abs(j.eta())<3.0:
        jetId = True if abs(j.eta())>2.4 else j.chargedHadronEnergyFraction()>0 and j.chargedMultiplicity()>0 and j.chargedEmEnergyFraction()<0.99
        return jetId and j.neutralHadronEnergyFraction()<0.99 and j.neutralEmEnergyFraction()<0.99 and j.chargedMultiplicity()+j.neutralMultiplicity()>1
    else:
        return j.neutralEmEnergyFraction()<0.9 and j.neutralMultiplicity()>10

products = {
  'jets': {'type': 'vector<reco::GenJet>',      'label':"ak4GenJets"},
  #'gp':   {'type': 'vector<reco::GenParticle>', 'label':"genParticles"},
}

plots = {
    'j0_pt': {'name': 'j0_pt', 'binning':[20,0,2000]},
    'j1_pt': {'name': 'j1_pt', 'binning':[20,0,2000]},
    'nJet':  {'name': 'nJet',  'binning':[10,0,10]},
}

# Make histos
for name, plot in plots.iteritems():
    plot['histos'] = [ROOT.TH1F("_".join([plot['name'],sample.name]),"_".join([plot['name'],sample.name]), *plot['binning']) for sample in samples]
    for i_sample, sample in enumerate(samples):
        plot['histos'][i_sample].style = styles.lineStyle( sample.color )
        plot['histos'][i_sample].legendText = sample.name 

for i_sample, sample in enumerate([LO, NLO]):
    r = sample.fwliteReader( products = products )
    r.start()
    count=0
    while r.run( ):
        if max_events is not None and max_events>0 and count>=max_events:break
        count += 1
        if count % 100 == 0: print "At %i" % count

        jets = [ j for j in r.products['jets'] if j.pt()>30 ] # if jetID( j )]
        #print jets
        #print len(jets)
        plots['nJet']['histos'][i_sample].Fill( len(jets) )
        if len(jets)>0:
            plots['j0_pt']['histos'][i_sample].Fill( jets[0].pt() )
        if len(jets)>1:
            plots['j1_pt']['histos'][i_sample].Fill( jets[1].pt() )


for name, plot in plots.iteritems():
    for h in  plot['histos']:
        if h.Integral()>0: h.Scale(1./h.Integral())
    p = Plot.fromHisto(name = name, histos = [ [h] for h in plot['histos']], texX = name )
    plotting.draw(p, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = {0:1}, yRange = (10**-3, 1), logY = True, logX = False)

