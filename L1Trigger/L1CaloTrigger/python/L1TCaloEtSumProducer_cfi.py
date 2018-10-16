import FWCore.ParameterSet.Config as cms

L1TCaloEtSumProducer = cms.EDProducer('L1TCaloEtSumProducer',
  ecalColl = cms.InputTag("l1EGammaCrystalsProducer","L1EGXtalClusterNoCuts"),
  ecalEtMin = cms.double(0.5),
  hcalDigis = cms.VInputTag(cms.InputTag('simHcalTriggerPrimitiveDigis')),
  hgCalColl = cms.InputTag("hgcalTriggerPrimitiveDigiProducer", "tower"),
  metEtaMax = cms.int32(28),
  metEtaMaxHF = cms.int32(41),
  ettEtaMax = cms.int32(28),
  ettEtaMaxHF = cms.int32(41),
  nTowThresholdHw = cms.int32(0),
  nTowEtaMax = cms.int32(41)
)
