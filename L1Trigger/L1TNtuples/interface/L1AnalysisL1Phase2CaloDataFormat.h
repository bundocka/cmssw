#ifndef __L1Analysis_L1AnalysisL1Phase2CaloDataFormat_H__
#define __L1Analysis_L1AnalysisL1Phase2CaloDataFormat_H__

//-------------------------------------------------------------------------------
// Created 20/04/2010 - E. Conte, A.C. Le Bihan
// 
// 
// Original code : L1TriggerDPG/L1Ntuples/L1Phase2CaloTreeProducer - Jim Brooke
//-------------------------------------------------------------------------------


#include <vector>

namespace L1Analysis
{

  // copied from DataFormats/L1Trigger/interface/EtSum.h, for use in standalone ROOT macros which use this class.
  struct L1AnalysisL1Phase2CaloDataFormat
  {
  
    L1AnalysisL1Phase2CaloDataFormat(){ Reset();};
    ~L1AnalysisL1Phase2CaloDataFormat(){};
    
    void Reset()
    {
      nSums = 0;
      sumType.clear();
      sumEt.clear();
      sumPhi.clear();
      sumIEt.clear();
      sumIPhi.clear();
      sumBx.clear();

    }
   
    unsigned short int nSums;
    std::vector<short int> sumType;
    std::vector<float> sumEt;
    std::vector<float> sumPhi;
    std::vector<short int> sumIEt;
    std::vector<short int> sumIPhi;
    std::vector<float> sumBx;

  }; 
}
#endif


