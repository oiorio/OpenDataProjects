import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

from mclist import slist
from mclistqcd import *
from mclisttt import slisttt
from mclistttdl import slistttdl
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/0019AB30-B9B7-E311-9E28-003048FF86CA.root",
#'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/DYJetsToLL_M-10To50_TuneZ2_7TeV-pythia6/AODSIM/PU_S13_START53_LV6-v1/00000/001988CE-B20B-E411-942E-0025905A6118.root'
#        'file:myfile.root'
    )
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("/tmp/oiorio/TestFileSkimQCDBig2.root")) 
process.source.fileNames=slist_v1
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(68224) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(262570) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1024) )
#process.source.fileNames=slist_v11
#process.TFileService.fileName=cms.string("/tmp/oiorio/TestFileSkimQCD_11.root")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.source.fileNames=slist
process.TFileService = cms.Service("TFileService", fileName = cms.string("/tmp/oiorio/TestFileSkimZData.root")) 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(550) )
process.source.fileNames=slisttt
process.TFileService.fileName=cms.string("/tmp/oiorio/TestFileSkimTTSmol.root")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
##process.source.fileNames=slistttdl
#process.TFileService.fileName=cms.string("/tmp/oiorio/TestFileSkimTTBig.root")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
#process.source.fileNames=slistttdl
#process.TFileService.fileName=cms.string("/tmp/oiorio/TestFileSkimTTDL.root")
                                                                                                                     
process.trees = cms.EDAnalyzer('AnalyzerTT',
muons=cms.InputTag("muons"),
electrons=cms.InputTag("gsfElectrons"),
globalMuons=cms.InputTag("globalMuons"),
jets=cms.InputTag("ak5PFJets"),
met=cms.InputTag("pfMet"),
vertexLabel=cms.InputTag("offlinePrimaryVertices"),
lhes = cms.InputTag('source'),
#lhes = cms.InputTag('externalLHEProducer'),
doMCMatch=cms.untracked.bool(False),
#doMCMatch=cms.untracked.bool(True),
)


process.p = cms.Path(process.trees)
