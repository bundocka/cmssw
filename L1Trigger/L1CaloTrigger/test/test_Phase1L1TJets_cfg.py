import FWCore.ParameterSet.Config as cms
from math import pi
import FWCore.Utilities.FileUtils as FileUtils # ADDED
import FWCore.ParameterSet.VarParsing as VarParsing # ADDED


process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

options = VarParsing.VarParsing ('analysis')
# get and parse the command line arguments

options.register('skipEvents',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to skip")
options.register('outFile',
                 'L1Ntuple.root',
                  VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file')

options.parseArguments()

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

fileList = FileUtils.loadListFromFile('ttqcd.list')
readFiles = cms.untracked.vstring(*fileList)

process.source = process.source = cms.Source("PoolSource",
   skipEvents = cms.untracked.uint32(options.skipEvents), #added
   fileNames = readFiles,
  #fileNames = cms.untracked.vstring(
  #  "file:pf500.root",
  #)
)

process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_cff')

process.load('L1Trigger.Phase2L1ParticleFlow.l1pfJetMet_cff')

#process.out = cms.OutputModule("PoolOutputModule",
#  fileName = cms.untracked.string('myOutputFile.root'),
#  outputCommands = cms.untracked.vstring(
#    "drop *",
#    "keep *_Phase1L1TJetProducer_*_*",
#    "keep *_ak4GenJetsNoNu_*_*",
#    "keep *_Phase1L1TJetCalibrator_*_*",
#    "keep *_*_*_*",
#  ),
#)

process.load("L1Trigger.L1TNtuples.l1PhaseIPFJetTreeProducer_cfi")


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('L1Ntuple.root')
)


process.p = cms.Path(process.Phase1L1TJetsSequence+process.l1PFJets)

#process.e = cms.EndPath(process.out)
