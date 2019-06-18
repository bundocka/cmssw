#ifndef RecoLocalCalo_HGCalRecAlgos_HGCalRecHitAbsAlgo_HH
#define RecoLocalCalo_HGCalRecAlgos_HGCalRecHitAbsAlgo_HH

/** \class HGCalRecHitAbsAlgo
  *  Template algorithm to make rechits from uncalibrated rechits
  *
  *
  *
  *  \author Valeri Andreev
  */

#include <vector>
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCUncalibratedRecHit.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

class HGCalRecHitAbsAlgo {
public:
  /// Constructor
  //HGCalRecHitAbsAlgo() { };

  /// Destructor
  virtual ~HGCalRecHitAbsAlgo(){};

  inline void getEventSetup(const edm::EventSetup& es) { rhtools_.getEventSetup(es); }

  /// make rechits from dataframes
  virtual void setLayerWeights(const std::vector<float>& weights){};

  virtual void setADCToGeVConstant(const float value) = 0;
  virtual HGCRecHit makeRecHit(const HGCUncalibratedRecHit& uncalibRH, const uint32_t& flags) const = 0;

protected:
  hgcal::RecHitTools rhtools_;
};
#endif
