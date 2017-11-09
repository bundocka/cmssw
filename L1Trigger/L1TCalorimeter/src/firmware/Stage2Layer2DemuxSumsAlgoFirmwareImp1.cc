///
/// \class l1t::Stage2Layer2SumsAlgorithmFirmwareImp1
///
/// \author:
///
/// Description:

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2DemuxSumsAlgoFirmware.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"

#include <vector>
#include <algorithm>


l1t::Stage2Layer2DemuxSumsAlgoFirmwareImp1::Stage2Layer2DemuxSumsAlgoFirmwareImp1(CaloParamsHelper* params) :
  params_(params), cordic_(Cordic(144*16,17,8))  // These are the settings in the hardware - should probably make this configurable
{
}


l1t::Stage2Layer2DemuxSumsAlgoFirmwareImp1::~Stage2Layer2DemuxSumsAlgoFirmwareImp1() {


}


void l1t::Stage2Layer2DemuxSumsAlgoFirmwareImp1::processEvent(const std::vector<l1t::EtSum> & inputSums,
                                                              std::vector<l1t::EtSum> & outputSums) {
  // set all energy sums back to zero!
  resetEnergySums();

  bool metSat(0), metHFSat(0), mhtSat(0), mhtHFSat(0);

  // Add up the x, y and scalar components
  for (auto&& eSum : inputSums)
  {
      switch (eSum.getType()) {

      case l1t::EtSum::EtSumType::kTotalEt:
        et_ += eSum.hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalEtEm:
        etem_ += eSum.hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalEtx:
	if(eSum.hwPt()==0x7fffffff) metSat=true;
        else metx_ += eSum.hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalEty:
	if(eSum.hwPt()==0x7fffffff) metSat=true;
        else mety_ += eSum.hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalHt:
        ht_ += eSum.hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalHtx:
	if(eSum.hwPt()==0x7fffffff) mhtSat=true;
        else mhtx_ += eSum.hwPt();
	break;

      case l1t::EtSum::EtSumType::kTotalHty:
	if(eSum.hwPt()==0x7fffffff) mhtSat=true;
	else mhty_ += eSum.hwPt();
        break;
	
      case l1t::EtSum::EtSumType::kTotalEtxHF:
	if(eSum.hwPt()==0x7fffffff) metHFSat=true;
        else metxHF_ += eSum.hwPt();
        break;
	
      case l1t::EtSum::EtSumType::kTotalEtyHF:
	if(eSum.hwPt()==0x7fffffff) metHFSat=true;
        else metyHF_ += eSum.hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalHtxHF:
	if(eSum.hwPt()==0x7fffffff) mhtHFSat=true;
	else mhtxHF_ += eSum.hwPt();
        break;
	
      case l1t::EtSum::EtSumType::kTotalHtyHF:
	if(eSum.hwPt()==0x7fffffff) mhtHFSat=true;
	else mhtyHF_ += eSum.hwPt();
        break;

      case l1t::EtSum::EtSumType::kMinBiasHFP0:
	mbp0_ = eSum.hwPt();
	break;

      case l1t::EtSum::EtSumType::kMinBiasHFM0:
	mbm0_ = eSum.hwPt();
	break;

      case l1t::EtSum::EtSumType::kMinBiasHFP1:
	mbp1_ = eSum.hwPt();
	break;

      case l1t::EtSum::EtSumType::kMinBiasHFM1:
	mbm1_ = eSum.hwPt();
	break;

      case l1t::EtSum::EtSumType::kTowerCount:
	ntow_ = eSum.hwPt();
	break;

      default:
        continue; // Should throw an exception or something?
      }
    }

  
  // Final MET calculation
  if ( (metx_ != 0 || mety_ != 0) && !metSat ) cordic_( metx_ , mety_ , metPhi_ , met_ );
  // sets the met scale back to the original range for output into GT, this corresponds to
  // the previous scaling of sin/cos factors in calculation of metx and mety by 2^10 = 1024
  met_ >>= 10; 

  // Final METHF calculation
  if ( (metxHF_ != 0 || metyHF_ != 0) && !metHFSat ) cordic_( metxHF_ , metyHF_ , metPhiHF_ , metHF_ );
  metHF_ >>= 10;


  // Final MHT calculation
  if ( (mhtx_ != 0 || mhty_ != 0) && !mhtSat ) cordic_( mhtx_ , mhty_ , mhtPhi_ , mht_ );
  // sets the mht scale back to the original range for output into GT, the other 4
  // bits are brought back just before the accumulation of ring sum in MP jet sum algorithm
  mht_ >>= 6; 

  if ( (mhtxHF_ != 0 || mhtyHF_ != 0) && !mhtHFSat ) cordic_( mhtxHF_ , mhtyHF_ , mhtPhiHF_ , mhtHF_ );
  mhtHF_ >>= 6; 


  if (et_>0xFFF)   et_   = 0xFFF;
  if (etem_>0xFFF) etem_ = 0xFFF;
  if (ht_>0xFFF)   ht_   = 0xFFF;

  if(metSat) met_=0xFFF;
  if(metHFSat) metHF_=0xFFF;
  if(mhtSat) mht_=0xFFF;
  if(mhtHFSat) mhtHF_=0xFFF;

  //if (mhtx_>0xFFF) mhtx_ = 0xFFF;
  //if (mhty_>0xFFF) mhty_ = 0xFFF;
  //mhtPhi_ = (111 << 4);
  //mhtPhiHF_ = (111 << 4); // to match hw value if undefined
  
  // calibrate energy sums
  calibrateEnergySums();

  // Make final collection
  math::XYZTLorentzVector p4;

  l1t::EtSum etSumTotalEt(p4,l1t::EtSum::EtSumType::kTotalEt,et_,0,0,0);
  l1t::EtSum etSumTotalEtEm(p4,l1t::EtSum::EtSumType::kTotalEtEm,etem_,0,0,0);
  l1t::EtSum etSumMissingEt(p4,l1t::EtSum::EtSumType::kMissingEt,met_,0,metPhi_>>4,0);
  l1t::EtSum etSumMissingEtHF(p4,l1t::EtSum::EtSumType::kMissingEtHF,metHF_,0,metPhiHF_>>4,0);
  l1t::EtSum htSumht(p4,l1t::EtSum::EtSumType::kTotalHt,ht_,0,0,0);
  l1t::EtSum htSumMissingHt(p4,l1t::EtSum::EtSumType::kMissingHt,mht_,0,mhtPhi_>>4,0);
  l1t::EtSum htSumMissingHtHF(p4,l1t::EtSum::EtSumType::kMissingHtHF,mhtHF_,0,mhtPhiHF_>>4,0);
  l1t::EtSum etSumMinBiasHFP0(p4,l1t::EtSum::EtSumType::kMinBiasHFP0,mbp0_,0,0,0);
  l1t::EtSum etSumMinBiasHFM0(p4,l1t::EtSum::EtSumType::kMinBiasHFM0,mbm0_,0,0,0);
  l1t::EtSum etSumMinBiasHFP1(p4,l1t::EtSum::EtSumType::kMinBiasHFP1,mbp1_,0,0,0);
  l1t::EtSum etSumMinBiasHFM1(p4,l1t::EtSum::EtSumType::kMinBiasHFM1,mbm1_,0,0,0);
  l1t::EtSum etSumTowCount(p4,l1t::EtSum::EtSumType::kTowerCount,ntow_,0,0,0);

  outputSums.push_back(etSumTotalEt);
  outputSums.push_back(etSumTotalEtEm);
  outputSums.push_back(etSumMinBiasHFP0);
  outputSums.push_back(htSumht);
  outputSums.push_back(etSumMinBiasHFM0);
  outputSums.push_back(etSumMissingEt);
  outputSums.push_back(etSumMinBiasHFP1);
  outputSums.push_back(htSumMissingHt);
  outputSums.push_back(etSumMinBiasHFM1);
  outputSums.push_back(etSumMissingEtHF);
  outputSums.push_back(htSumMissingHtHF);
  outputSums.push_back(etSumTowCount);

}


void l1t::Stage2Layer2DemuxSumsAlgoFirmwareImp1::calibrateEnergySums()
{

  if(params_->etSumCalibrationType() == "LUT"){

    //calibrate Metx & Mety
    unsigned int metxShiftSat = ( metx_ >> 10 ) & 0x7FF; // saturates at 2 TeV
    unsigned int metyShiftSat = ( mety_ >> 10 ) & 0x7FF;
    bool metxIsNeg = ( metx_ >> 31 ) & 0x1;
    bool metyIsNeg = ( mety_ >> 31 ) & 0x1;
    
    unsigned int metxComp = params_->etSumCompressionLUT()->data( metxShiftSat | ( metxIsNeg << 11 ) );
    unsigned int metyComp = params_->etSumCompressionLUT()->data( metyShiftSat | ( metyIsNeg << 11 ) );
    
    unsigned int metLUTAddr = ( metyComp << 6 ) | metxComp;

    unsigned int metAddMult = params_->etSumCalibrationLUT()->data(metLUTAddr);
    unsigned int metAdd = ( metAddMult >> 10 ) & 0xFF;
    unsigned int metMult = metAddMult & 0x3FF;

    unsigned int calibMet =  ( ( met_ ) * ( metMult / 512 ) ) + ( metAdd );
    if (calibMet > 0x1FFFFFF) calibMet = 0x1FFFFFF;
    

  } else {
    if(params_->etSumCalibrationType() != "None" && params_->etSumCalibrationType() != "none") 
      edm::LogError("l1t|stage 2") << "Invalid etSum calibration type in calo params. Not applying etSum calibration to Stage 2 MET" << std::endl;
    return;
  } 

  if(params_->etSumPhiCalibrationType() == "LUT"){


  } else {
    if(params_->etSumPhiCalibrationType() != "None" && params_->etSumPhiCalibrationType() != "none")
      edm::LogError("l1t|stage 2") << "Invalid etSum Phi calibration type in calo params. Not applying etSum Phi calibration to Stage 2 MET" << std::endl;
    return;
  }
  
  if(params_->etSumHFCalibrationType() == "LUT"){
    
    //calibrate MetxHF and MetyHF
    unsigned int metxHFShiftSat = ( metxHF_ >> 10 ) & 0x7FF;
    unsigned int metyHFShiftSat = ( metyHF_ >> 10 ) & 0x7FF;
    bool metxHFIsNeg = ( metxHF_ >> 31 ) & 0x1;
    bool metyHFIsNeg = ( metyHF_ >> 31 ) & 0x1;
    
    unsigned int metxHFComp = params_->etSumCompressionLUT()->data( metxHFShiftSat | ( metxHFIsNeg << 11 ) );
    unsigned int metyHFComp = params_->etSumCompressionLUT()->data( metyHFShiftSat | ( metyHFIsNeg << 11 ) );

    unsigned int metHFLUTAddr = ( metyHFComp << 6 ) | metxHFComp;
    
    unsigned int metHFAddMult = params_->etSumCalibrationLUT()->data(metHFLUTAddr);
    unsigned int metHFAdd = ( metHFAddMult >> 10 ) & 0xFF;
    unsigned int metHFMult = metHFAddMult & 0x3FF;

    unsigned int calibMetHF =  ( ( metHF_ ) * ( metHFMult / 512 ) ) + ( metHFAdd );
    if (calibMetHF > 0x1FFFFFF) calibMetHF = 0x1FFFFFF;


  } else {
    if(params_->etSumHFCalibrationType() != "None" && params_->etSumHFCalibrationType() != "none")
      edm::LogError("l1t|stage 2") << "Invalid etSumHF calibration type in calo params. Not applying etSumHF calibration to Stage 2 MET" << std::endl;
    return;
  }
   

  // calibrate metHF Phi
  if(params_->etSumHFPhiCalibrationType() == "LUT"){
    



    
  } else {
    if(params_->etSumHFPhiCalibrationType() != "None" && params_->etSumHFPhiCalibrationType() != "none")
      edm::LogError("l1t|stage 2") << "Invalid etSumHF Phi calibration type in calo params. Not applying etSumHF Phi calibration to Stage 2 MET" << std::endl;
    return;
  }
  


  
  if(params_->etSumEttCalibrationType() == "LUT"){
    //calibrate Et
    unsigned int etCalibLUTAddr = et_;
    unsigned int etAdd  = ( ( params_->etSumEttCalibrationLUT()->data(etCalibLUTAddr) >> 10 ) & 0xFF ); 
    unsigned int etMult = ( params_->etSumEttCalibrationLUT()->data(etCalibLUTAddr) & 0x3FF );
    unsigned int etCorr = ( et_ * (etMult/512) ) + etAdd;
    et_ = etCorr;
 
    //calibrate EtHF
    unsigned int etHFCalibLUTAddr = etHF_;
    unsigned int etHFAdd  = ( ( params_->etSumEttCalibrationLUT()->data(etHFCalibLUTAddr) >> 10 ) & 0xFF ); 
    unsigned int etHFMult = ( params_->etSumEttCalibrationLUT()->data(etHFCalibLUTAddr) & 0x3FF );
    unsigned int etHFCorr = ( etHF_ * (etHFMult/512) ) + etHFAdd;
    etHF_ = etHFCorr;
 
  } else {
    if(params_->etSumEttCalibrationType() != "None" && params_->etSumEttCalibrationType() != "none") 
      edm::LogError("l1t|stage 2") << "Invalid etSumEtt calibration type in calo params. Not applying etSumEtt calibration to Stage 2 ETT" << std::endl;
    return;
  }
 
  if(params_->etSumEcalSumCalibrationType() == "LUT"){
    //calibrate Etem
    unsigned int etemCalibLUTAddr = etem_;
    unsigned int etemAdd  = ( ( params_->etSumEcalSumCalibrationLUT()->data(etemCalibLUTAddr) >> 10 ) & 0xFF ); 
    unsigned int etemMult = ( params_->etSumEcalSumCalibrationLUT()->data(etemCalibLUTAddr) & 0x3FF );
    unsigned int etemCorr = ( etem_ * (etemMult/512) ) + etemAdd;
    etem_ = etemCorr;
 
  } else {
    if(params_->etSumEcalSumCalibrationType() != "None" && params_->etSumEcalSumCalibrationType() != "none") 
      edm::LogError("l1t|stage 2") << "Invalid etSumEcalSum calibration type in calo params. Not applying etSumEcalSum calibration to Stage 2 ECal Et" << std::endl;
    return;
  }
  
   
}

void l1t::Stage2Layer2DemuxSumsAlgoFirmwareImp1::resetEnergySums()
{
  
  et_= 0;
  etHF_ = 0;
  etem_= 0;    
  metx_= 0;    
  mety_= 0;    
  metxHF_= 0;  
  metyHF_= 0;  
  met_= 0;     
  metHF_= 0;
  calibMet_ = 0;
  calibMetHF_ = 0;
  ht_= 0;      
  mht_= 0;     
  mhtx_= 0;    
  mhty_= 0;    
  mhtxHF_= 0;  
  mhtyHF_= 0;  
  mhtHF_= 0;   
  metPhi_= 0;  
  metPhiHF_= 0;
  calibMetPhi_ = 0;
  calibMetHFPhi_ = 0;
  mhtPhi_= 0;  
  mhtPhiHF_= 0;
  mbp0_= 0;    
  mbm0_= 0;    
  mbp1_= 0;    
  mbm1_= 0;    
  ntow_= 0;                  

}
