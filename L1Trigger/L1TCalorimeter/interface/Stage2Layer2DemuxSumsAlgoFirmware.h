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

#ifndef Stage2Layer2DemuxSumsAlgoFirmware_H
#define Stage2Layer2DemuxSumsAlgoFirmware_H

#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2DemuxSumsAlgo.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"
#include "L1Trigger/L1TCalorimeter/interface/Cordic.h"

namespace l1t {

  // Imp1 is for v1 and v2
  class Stage2Layer2DemuxSumsAlgoFirmwareImp1 : public Stage2Layer2DemuxSumsAlgo {
  public:
    Stage2Layer2DemuxSumsAlgoFirmwareImp1(CaloParamsHelper* params);
    virtual ~Stage2Layer2DemuxSumsAlgoFirmwareImp1();
    virtual void processEvent(const std::vector<l1t::EtSum> & inputSums,
			      std::vector<l1t::EtSum> & outputSums);
    
    virtual void calibrateEnergySums();
    virtual void resetEnergySums();
    
  private:

    CaloParamsHelper* params_;

    Cordic cordic_;

    unsigned int et_;             			     
    unsigned int etHF_;	    			     
    unsigned int etem_;	    			     
    int  metx_;	    			     
    int  mety_;	    			     
    int  metxHF_;	    			     
    int  metyHF_;	    			     
    unsigned int met_;	    			     
    unsigned int metHF_;	    			     
    unsigned int calibMet_;	    			     
    unsigned int calibMetHF_;   			     
    unsigned int ht_;	    			     
    unsigned int mht_;	    			     
    int  mhtx_;	    			     
    int  mhty_;	    			     
    int  mhtxHF_;	    			     
    int  mhtyHF_;	    			     
    unsigned int mhtHF_;	    			     
    int metPhi_;	    			     
    int metPhiHF_;	    			     
    int calibMetPhi_;   			     
    int calibMetHFPhi_; 			     
    int mhtPhi_;	    			     
    int mhtPhiHF_;	    			     
    unsigned int mbp0_;	     
    unsigned int mbm0_;	     
    unsigned int mbp1_;	     
    unsigned int mbm1_;	     
    unsigned int ntow_;

  };

}

#endif
