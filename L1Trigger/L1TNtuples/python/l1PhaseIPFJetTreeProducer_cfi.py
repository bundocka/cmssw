import FWCore.ParameterSet.Config as cms

l1PhaseIPFJetTree = cms.EDAnalyzer("L1PhaseIPFJetTreeProducer",
   l1PhaseIPFJets = cms.InputTag("Phase1L1TJetCalibrator", "Phase1L1TJetFromPfCandidates"),
   maxL1Extra = cms.uint32(20)
)

#### Gen level tree

from L1Trigger.L1TNtuples.l1GeneratorTree_cfi  import l1GeneratorTree
genTree=l1GeneratorTree.clone()

runmenutree=cms.Path(l1PhaseIPFJetTree*genTree)




