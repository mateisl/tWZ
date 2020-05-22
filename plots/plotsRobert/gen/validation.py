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

import tWZ.Tools.user as user
import tWZ.samples.GENandSIM_private as samples
import tWZ.Tools.helpers as helpers

LO      = samples.tOrtbar_WZ01j_OLRLL_LO 
NLO     = samples.tWZ_NLO
central = samples.tWZ_LO_central 

LO.color = ROOT.kBlue
NLO.color= ROOT.kRed
central.color = ROOT.kGreen

def jetID(j):
    if abs(j.eta())<3.0:
        jetId = True if abs(j.eta())>2.4 else j.chargedHadronEnergyFraction()>0 and j.chargedMultiplicity()>0 and j.chargedEmEnergyFraction()<0.99
        return jetId and j.neutralHadronEnergyFraction()<0.99 and j.neutralEmEnergyFraction()<0.99 and j.chargedMultiplicity()+j.neutralMultiplicity()>1
    else:
        return j.neutralEmEnergyFraction()<0.9 and j.neutralMultiplicity()>10

def genLep( p, mother_pdgIds = None):
    return abs(p.pdgId()) in [11, 13] and p.numberOfMothers()==1 and ( mother_pdgIds is None or abs(p.mother(0).pdgId()) in mother_pdgIds)

GEN_products = {
  'jets': {'type':'vector<reco::GenJet>',       'label':"ak4GenJets"},
  'gp':   {'type':'vector<reco::GenParticle>',  'label':"genParticles"},
}
MAOD_products = {
  'jets': {'type': 'vector<reco::GenJet>',          'label':"slimmedGenJets"},
  'gp':   {'type':'vector<pat::PackedGenParticle>', 'label':"packedGenParticles"},
}

LO.products  = GEN_products
NLO.products = GEN_products
central.products = MAOD_products

samples = [central, LO, NLO]

plots = {
    'j0_pt': {'name': 'j0_pt', 'binning':[20,0,2000]},
    'j0_nd': {'name': 'j0_nd', 'binning':[20,0,50]},
    'j1_pt': {'name': 'j1_pt', 'binning':[20,0,2000]},
    'nJet':  {'name': 'nJet',  'binning':[12,0,12]},
    'nJet_central':  {'name': 'nJet_central',  'binning':[12,0,12]},
}

# Make histos
for name, plot in plots.iteritems():
    plot['histos'] = [ROOT.TH1F("_".join([plot['name'],sample.name]),"_".join([plot['name'],sample.name]), *plot['binning']) for sample in samples]
    for i_sample, sample in enumerate(samples):
        plot['histos'][i_sample].style = styles.lineStyle( sample.color )
        plot['histos'][i_sample].legendText = sample.name 

for i_sample, sample in enumerate(samples):
    sample.reduceFiles(to=max_files)
    r = sample.fwliteReader( products = sample.products )
    r.start()
    count=0
    while r.run( ):
        if max_events is not None and max_events>0 and count>=max_events:break
        count += 1
        if count % 10 == 0: print "At %i" % count

        lep_from_W = filter( lambda p: genLep(p, mother_pdgIds=[24]), r.products['gp'] )
        lep_from_Z = filter( lambda p: genLep(p, mother_pdgIds=[23]), r.products['gp'] )

        if not (len(lep_from_Z)==2 and len(lep_from_W)==1): continue

        jets = [ j for j in r.products['jets'] if j.pt()>30 and min([helpers.deltaR({'phi':j.phi(), 'eta':j.eta()}, {'phi':l.phi(),'eta':l.eta()}) for l in lep_from_W+lep_from_Z],999.)>0.4 ] # if jetID( j )]
        #print jets
        #print len(jets)
        plots['nJet']['histos'][i_sample].Fill( len(jets) )
        plots['nJet_central']['histos'][i_sample].Fill( len([j for j in jets if abs(j.eta())<2.4]) )
        if len(jets)>0:
            plots['j0_pt']['histos'][i_sample].Fill( jets[0].pt() )
            plots['j0_nd']['histos'][i_sample].Fill( jets[0].numberOfDaughters() )
        if len(jets)>1:
            plots['j1_pt']['histos'][i_sample].Fill( jets[1].pt() )



for name, plot in plots.iteritems():
    for h in  plot['histos']:
        if h.Integral()>0: h.Scale(1./h.Integral())
    p = Plot.fromHisto(name = name, histos = [ [h] for h in plot['histos']], texX = name )
    plotting.draw(p, plot_directory = os.path.join(user.plot_directory, "validation"), ratio = {0:1}, yRange = (10**-3, 1), logY = True, logX = False)

