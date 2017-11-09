///
/// Description: Firmware headers
///
/// Implementation:
///    Concrete firmware implementations
///
/// \author: Jim Brooke - University of Bristol
///

//
//

#ifndef Stage2Layer2EtSumAlgorithmFirmware_H
#define Stage2Layer2EtSumAlgorithmFirmware_H

#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2EtSumAlgorithm.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"

namespace l1t {

  // Imp1 is for v1 and v2
  class Stage2Layer2EtSumAlgorithmFirmwareImp1 : public Stage2Layer2EtSumAlgorithm {
  public:
    Stage2Layer2EtSumAlgorithmFirmwareImp1(CaloParamsHelper* params);
    virtual ~Stage2Layer2EtSumAlgorithmFirmwareImp1();
    virtual void processEvent(const std::vector<l1t::CaloTower> & towers,
			      std::vector<l1t::EtSum> & sums);
    virtual void resetEnergySums();

  private:
    CaloParamsHelper* params_;
    int towEtMetThresh_;
    int towEtSumEtThresh_;
    int towEtEcalSumThresh_;
    int metEtaMax_;
    int metEtaMaxHF_;
    int ettEtaMax_;
    int ettEtaMaxHF_;
    int nTowThresholdHw_;
    int nTowEtaMax_;
    int ex_, ey_, et_;	
    int exHF_, eyHF_, etHF_;
    int etem_;		
    unsigned int mb0_, mb1_;	
    unsigned int ntowers_;		
  };
}

#endif
