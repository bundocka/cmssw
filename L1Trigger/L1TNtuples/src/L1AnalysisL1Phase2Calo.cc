#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1Phase2Calo.h"

L1Analysis::L1AnalysisL1Phase2Calo::L1AnalysisL1Phase2Calo()
{
}

L1Analysis::L1AnalysisL1Phase2Calo::~L1AnalysisL1Phase2Calo()
{

}

void L1Analysis::L1AnalysisL1Phase2Calo::SetSum(const edm::Handle<l1t::EtSumBxCollection> sums, unsigned maxL1Phase2Calo)
{
  for (int ibx = sums->getFirstBX(); ibx <= sums->getLastBX(); ++ibx) {
    for (l1t::EtSumBxCollection::const_iterator it=sums->begin(ibx); it!=sums->end(ibx) && l1Phase2Calo_.nSums<maxL1Phase2Calo; it++) {
      int type = static_cast<int>( it->getType() ); 
      l1Phase2Calo_.sumType. push_back( type ); 
      l1Phase2Calo_.sumEt. push_back( it->et() ); 
      l1Phase2Calo_.sumPhi.push_back( it->phi() );
      l1Phase2Calo_.sumIEt. push_back( it->hwPt() ); 
      l1Phase2Calo_.sumIPhi.push_back( it->hwPhi() );
      l1Phase2Calo_.sumBx. push_back( ibx );
      l1Phase2Calo_.nSums++;
    }
  }
}


