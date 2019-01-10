import FWCore.ParameterSet.Config as cms

l1Phase2CaloTree = cms.EDAnalyzer(
    "L1Phase2CaloTreeProducer",
    sumToken = cms.untracked.InputTag("caloPhase2Digis","EtSum"),
    maxL1Phase2Calo = cms.uint32(60)
)



