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

void L1Analysis::L1AnalysisPhaseIPFJet::SetPFJet(const edm::Handle< vector<l1t::PFJet> > ak4PFJets, unsigned maxL1Extra)
{
  for (unsigned int i=0; i<ak4PFJets->size() && l1extra_.nPhaseIPFJets<maxL1Extra; i++){
    if (ak4PFJets->at(i).pt()){
      l1extra_.ak4PFJetEt .push_back(ak4PFJets->at(i).pt());
      l1extra_.ak4PFJetEta.push_back(ak4PFJets->at(i).eta());
      l1extra_.ak4PFJetPhi.push_back(ak4PFJets->at(i).phi());
      l1extra_.nPhaseIPFJets++;
    }
  }
}
