#include "L1Trigger/L1TNtuples/interface/L1AnalysisPhaseIPFJet.h"

L1Analysis::L1AnalysisPhaseIPFJet::L1AnalysisPhaseIPFJet()
{
}

L1Analysis::L1AnalysisPhaseIPFJet::~L1AnalysisPhaseIPFJet()
{

}



void L1Analysis::L1AnalysisPhaseIPFJet::SetPhaseIPFJet(const edm::Handle< vector<reco::CaloJet> > phaseIPFJets, unsigned maxL1Extra)
{
  for (unsigned int i=0; i<phaseIPFJets->size() && l1extra_.nPhaseIPFJets<maxL1Extra; i++){
    if (phaseIPFJets->at(i).pt()){
      l1extra_.phaseIPFJetEt .push_back(phaseIPFJets->at(i).pt());
      l1extra_.phaseIPFJetEta.push_back(phaseIPFJets->at(i).eta());
      l1extra_.phaseIPFJetPhi.push_back(phaseIPFJets->at(i).phi());
      l1extra_.nPhaseIPFJets++;
    }
  }
}
