#ifndef __L1Analysis_L1AnalysisPhaseIPFJetDataFormat_H__
#define __L1Analysis_L1AnalysisPhaseIPFJetDataFormat_H__

//-------------------------------------------------------------------------------
// Created 20/04/2010 - E. Conte, A.C. Le Bihan
// 
// 
// Original code : UserCode/L1TriggerDPG/L1ExtraTreeProducer - Jim Brooke
//-------------------------------------------------------------------------------


#include <vector>

namespace L1Analysis
{
  struct L1AnalysisPhaseIPFJetDataFormat
  {
    L1AnalysisPhaseIPFJetDataFormat(){Reset();};
    ~L1AnalysisPhaseIPFJetDataFormat(){};
    
    void Reset()
    {

      nPhaseIPFJets = 0;
      phaseIPFJetEt.clear();
      phaseIPFJetEta.clear();
      phaseIPFJetPhi.clear();

    }
 
    unsigned short int nPhaseIPFJets;
    std::vector<double> phaseIPFJetEt;
    std::vector<double> phaseIPFJetEta;
    std::vector<double> phaseIPFJetPhi;


  }; 
}
#endif


