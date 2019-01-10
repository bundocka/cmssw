#ifndef __L1Analysis_L1AnalysisL1Phase2Calo_H__
#define __L1Analysis_L1AnalysisL1Phase2Calo_H__

//-------------------------------------------------------------------------------
// Created 02/08/2018 - A.Bundock
// 
// 
// Original code : L1Trigger/L1Ntuples/L1AnalysisL1Upgrade - Jim Brooke
//-------------------------------------------------------------------------------

#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "L1AnalysisL1Phase2CaloDataFormat.h"

namespace L1Analysis
{
  class L1AnalysisL1Phase2Calo 
  {
  public:
    enum {TEST=0};
    L1AnalysisL1Phase2Calo();
    ~L1AnalysisL1Phase2Calo();
    void Reset() {l1Phase2Calo_.Reset();}
    void SetSum  (const edm::Handle<l1t::EtSumBxCollection>  sums, unsigned maxL1Phase2Calo);
    L1AnalysisL1Phase2CaloDataFormat * getData() {return &l1Phase2Calo_;}

  private :
    L1AnalysisL1Phase2CaloDataFormat l1Phase2Calo_;
  }; 
}
#endif


