import FWCore.ParameterSet.Config as cms

l1Phase2CaloTree = cms.EDAnalyzer(
    "L1Phase2CaloTreeProducer",
    sumToken = cms.untracked.InputTag("simCaloPhase2Digis",""),
    maxL1Phase2Calo = cms.uint32(60)
)



