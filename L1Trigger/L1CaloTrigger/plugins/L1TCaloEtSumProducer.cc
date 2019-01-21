// -*- C++ -*-
//
// Package:    L1Trigger/L1TCaloEtSumProducer
// Class:      L1TCaloEtSumProducer
// 
/**\class L1TCaloEtSumProducer L1TCaloEtSumProducer.cc L1Trigger/L1TCaloEtSumProducer/plugins/L1TCaloEtSumProducer.cc

 Description: Phase 2 L1T EtSum (MET/ETT) producer

 Implementation:
     Ported from Phase 1 Stage 2
*/
//
// Original Author:  Aaron Bundock
//         Created:  Wed, 25 Jul 2018 16:25:51 GMT


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/Phase2L1CaloTrig/interface/L1EGCrystalCluster.h"
#include "DataFormats/L1THGCal/interface/HGCalTower.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"
#include "L1Trigger/L1TCalorimeter/interface/Cordic.h"


// class declaration
class L1TCaloEtSumProducer : public edm::stream::EDProducer<> {
public:
  explicit L1TCaloEtSumProducer(const edm::ParameterSet&);
  ~L1TCaloEtSumProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;

  void buildTowers(edm::Event &event, const edm::EventSetup&);
  
  Cordic cordic_;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------

  edm::EDGetTokenT<l1slhc::L1EGCrystalClusterCollection> ecalColl_;
  double ecalEtMin_;

  std::vector<edm::EDGetTokenT<HcalTrigPrimDigiCollection>> hcalDigis_;
  edm::ESHandle<CaloTPGTranscoder> decoder_;

  edm::EDGetTokenT<l1t::HGCalTowerBxCollection> hgCalColl_;

  std::vector<l1t::CaloTower> towers_;
  std::vector<l1t::EtSum> etsums_;

  int metEtaMax_;
  int metEtaMaxHF_;
  int ettEtaMax_;
  int ettEtaMaxHF_;
  int nTowThresholdHw_;
  int nTowEtaMax_;
  
};

// constants, enums and typedefs

// static data member definitions

// constructors and destructor


L1TCaloEtSumProducer::L1TCaloEtSumProducer(const edm::ParameterSet& iConfig) :
  cordic_(Cordic(144*16,17,8)),
  ecalColl_(consumes<l1slhc::L1EGCrystalClusterCollection>(iConfig.getParameter<edm::InputTag>("ecalColl"))),
  ecalEtMin_(iConfig.getParameter<double>("ecalEtMin")),
  hgCalColl_(consumes<l1t::HGCalTowerBxCollection>(iConfig.getParameter<edm::InputTag>("hgCalColl"))),
  metEtaMax_(iConfig.getParameter<int>("metEtaMax")),
  metEtaMaxHF_(iConfig.getParameter<int>("metEtaMaxHF")),
  ettEtaMax_(iConfig.getParameter<int>("ettEtaMax")),
  ettEtaMaxHF_(iConfig.getParameter<int>("ettEtaMaxHF")),
  nTowThresholdHw_(iConfig.getParameter<int>("nTowThresholdHw")),
  nTowEtaMax_(iConfig.getParameter<int>("nTowEtaMax"))
{

  //register your products
  produces<l1t::CaloTowerBxCollection> ();
  produces<l1t::EtSumBxCollection>();

   //if do put with a label
  //produces<ExampleData2>("label");
 
   //if you want to put into the Run
   //produces<ExampleData2,InRun>();

  //now do what ever other initialization is needed
  for (auto & tag : iConfig.getParameter<std::vector<edm::InputTag>>("hcalDigis")) {
    hcalDigis_.push_back(consumes<HcalTrigPrimDigiCollection>(tag));
  }

}


L1TCaloEtSumProducer::~L1TCaloEtSumProducer()
{
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
}


// member functions

// ------------ method called to produce the data  ------------
void
L1TCaloEtSumProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  etsums_.clear();
  towers_.clear();
  
  buildTowers(iEvent, iSetup);

  math::XYZTLorentzVector p4;
  int ntowers = 0;

  int ex(0), ey(0), et(0);
  int exHF(0), eyHF(0), etHF(0);
  int etem(0);
  bool ettSat(false), ettHFSat(false), ecalEtSat(false), metSat(false), metHFSat(false);

  // etaSide=1 is positive eta, etaSide=-1 is negative eta
  for (int etaSide=1; etaSide>=-1; etaSide-=2) {

    unsigned int mb0(0), mb1(0);
    int towEtMetThresh_ = 0;
    int towEtSumEtThresh_ = 0;
    int towEtEcalSumThresh_ = 0;

    for (unsigned absieta=1; absieta<=(unsigned int)l1t::CaloTools::mpEta(l1t::CaloTools::kHFEnd); absieta++) {

      int ieta = etaSide * absieta;
      
      int ringEx(0), ringEy(0), ringEt(0);
      int ringExHF(0), ringEyHF(0), ringEtHF(0);
      int ringEtEm(0);
      unsigned int ringNtowers(0);
      unsigned int ringMB0(0), ringMB1(0);

      for (int iphi=1; iphi<=l1t::CaloTools::kHBHENrPhi; iphi++) {
	
	l1t::CaloTower tower = l1t::CaloTools::getTower(towers_, l1t::CaloTools::caloEta(ieta), iphi);


	// MET without HF

	if (tower.hwPt()>towEtMetThresh_ && l1t::CaloTools::mpEta(abs(tower.hwEta()))<=l1t::CaloTools::mpEta(metEtaMax_) && !metSat) {

	  // x- and -y coefficients are truncated by after multiplication of Et by trig coefficient.
	  // The trig coefficients themselves take values [-1023,1023] and so were scaled by
	  // 2^10 = 1024, which requires bitwise shift to the right of the final value by 10 bits.
	  // This is accounted for at ouput of demux (see Stage2Layer2DemuxSumsAlgoFirmwareImp1.cc)
	  if(tower.hwPt() == l1t::CaloTools::kSatHcal || tower.hwPt() == l1t::CaloTools::kSatEcal || tower.hwPt() == l1t::CaloTools::kSatTower) metSat=true;
	  ringEx += (int) (tower.hwPt() * l1t::CaloTools::cos_coeff[iphi - 1] );
	  ringEy += (int) (tower.hwPt() * l1t::CaloTools::sin_coeff[iphi - 1] );	    
	}

	// MET *with* HF
	if (tower.hwPt()>towEtMetThresh_ && l1t::CaloTools::mpEta(abs(tower.hwEta()))<=l1t::CaloTools::mpEta(metEtaMaxHF_) && !metHFSat) {
	  if(tower.hwPt() == l1t::CaloTools::kSatHcal || tower.hwPt() == l1t::CaloTools::kSatEcal || tower.hwPt() == l1t::CaloTools::kSatTower) metHFSat=true;
	  ringExHF += (int) (tower.hwPt() * l1t::CaloTools::cos_coeff[iphi - 1] );
	  ringEyHF += (int) (tower.hwPt() * l1t::CaloTools::sin_coeff[iphi - 1] );	    
	}

	// scalar sum
	if (tower.hwPt()>towEtSumEtThresh_ && l1t::CaloTools::mpEta(abs(tower.hwEta()))<=l1t::CaloTools::mpEta(ettEtaMax_) && !ettSat){
	  if(tower.hwPt() == l1t::CaloTools::kSatHcal || tower.hwPt() == l1t::CaloTools::kSatEcal || tower.hwPt() == l1t::CaloTools::kSatTower) ettSat=true;
	  ringEt += tower.hwPt();
	}
  
	// scalar sum including HF
	if (tower.hwPt()>towEtSumEtThresh_ && l1t::CaloTools::mpEta(abs(tower.hwEta()))<=l1t::CaloTools::mpEta(ettEtaMaxHF_) && !ettHFSat) {
	  if(tower.hwPt() == l1t::CaloTools::kSatHcal || tower.hwPt() == l1t::CaloTools::kSatEcal || tower.hwPt() == l1t::CaloTools::kSatTower) ettHFSat=true;
	  ringEtHF += tower.hwPt();
	}
	
        // scalar sum (EM)
        if (tower.hwPt()>towEtEcalSumThresh_ && l1t::CaloTools::mpEta(abs(tower.hwEta()))<=l1t::CaloTools::mpEta(ettEtaMax_) && !ecalEtSat){
	  if(tower.hwPt() == l1t::CaloTools::kSatEcal || tower.hwPt() == l1t::CaloTools::kSatTower) ecalEtSat=true;
          ringEtEm += tower.hwEtEm();
	}

	// count HF tower HCAL flags
	if (l1t::CaloTools::mpEta(abs(tower.hwEta()))>=l1t::CaloTools::mpEta(l1t::CaloTools::kHFBegin) &&
	    l1t::CaloTools::mpEta(abs(tower.hwEta()))<=l1t::CaloTools::mpEta(l1t::CaloTools::kHFEnd) &&
	    (tower.hwQual() & 0x4) > 0) 
	  ringMB0 += 1;
	  
        // tower counting 
	if (tower.hwPt()>nTowThresholdHw_ && l1t::CaloTools::mpEta(abs(tower.hwEta()))<=nTowEtaMax_) 
	  ringNtowers += 1;
      }    
      
      ex += ringEx;
      ey += ringEy;
      et += ringEt;
      etHF += ringEtHF;
      exHF += ringExHF;
      eyHF += ringEyHF;

      etem  += ringEtEm;

      mb0 += ringMB0;
      mb1 += ringMB1;

      ntowers += ringNtowers;
    }

    if (mb0>0xf) mb0 = 0xf;
    if (mb1>0xf) mb1 = 0xf;


    l1t::EtSum::EtSumType type0 = l1t::EtSum::EtSumType::kMinBiasHFP0;
    l1t::EtSum::EtSumType type1 = l1t::EtSum::EtSumType::kMinBiasHFP1;
    if (etaSide<0) {
      type0 = l1t::EtSum::EtSumType::kMinBiasHFM0;
      type1 = l1t::EtSum::EtSumType::kMinBiasHFM1;
    } 
    l1t::EtSum etSumMinBias0(p4,type0,mb0,0,0,0);
    l1t::EtSum etSumMinBias1(p4,type1,mb1,0,0,0);
    
    etsums_.push_back(etSumMinBias0);
    etsums_.push_back(etSumMinBias1);
    
  }

  // saturate energy sums if saturated TP/tower

  if(ecalEtSat || etem > 0xFFF) etem = 0xFFF;
  if(ettSat || et > 0xFFF) et = 0xFFF;
  if(ettHFSat || etHF > 0xFFF) etHF = 0xFFF;
  
  unsigned int met(0), metHF(0);
  int metPhi(0), metPhiHF(0);

  // Final MET calculation
  if ( (ex != 0 || ey != 0) && !metSat ) cordic_( ex , ey , metPhi , met );
  // sets the met scale back to the original range for output into GT, this corresponds to
  // the previous scaling of sin/cos factors in calculation of metx and mety by 2^10 = 1024
  met >>= 10; 

  // Final METHF calculation
  if ( (exHF != 0 || eyHF != 0) && !metHFSat ) cordic_( exHF , eyHF , metPhiHF , metHF );
  metHF >>= 10;

  if(metSat || met > 0xFFF) met=0xFFF;
  if(metHFSat || metHF > 0xFFF) metHF=0xFFF;

  // Make final collection

  l1t::EtSum etSumTotalEt(p4,l1t::EtSum::EtSumType::kTotalEt,et,0,0,0);
  l1t::EtSum etSumTotalEtEm(p4,l1t::EtSum::EtSumType::kTotalEtEm,etem,0,0,0);
  l1t::EtSum etSumMissingEt(p4,l1t::EtSum::EtSumType::kMissingEt,met,0,metPhi>>4,0);
  l1t::EtSum etSumMissingEtHF(p4,l1t::EtSum::EtSumType::kMissingEtHF,metHF,0,metPhiHF>>4,0);
  l1t::EtSum etSumTowCount(p4,l1t::EtSum::EtSumType::kTowerCount,ntowers,0,0,0);

  etsums_.push_back(etSumTotalEt);
  etsums_.push_back(etSumTotalEtEm);
  etsums_.push_back(etSumMissingEt);
  etsums_.push_back(etSumMissingEtHF);
  etsums_.push_back(etSumTowCount);


  std::unique_ptr<l1t::EtSumBxCollection> etsums (new l1t::EtSumBxCollection(0, 0, 0));
  for( auto etsum = etsums_.begin(); etsum != etsums_.end(); ++etsum)
    etsums->push_back(0, l1t::CaloTools::etSumP4Demux(*etsum));
  
  iEvent.put(std::move(etsums));

 
}


void
L1TCaloEtSumProducer::buildTowers(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // ECAL 5x5
  edm::Handle<l1slhc::L1EGCrystalClusterCollection> ecalClusters;
  iEvent.getByToken(ecalColl_, ecalClusters);

  for(auto it = ecalClusters->begin(), ed = ecalClusters->end(); it != ed; ++it) {
    if (it->e5x5() <= ecalEtMin_) continue;
    //if (it->e5x5() > 2) std::cout << "e5x5 = " << it->e5x5() << std::endl;
    l1t::CaloTower tow = l1t::CaloTools::getTower(towers_, l1t::CaloTools::caloEta(it->hwEta()), it->hwPhi());
    if(tow.hwPt()==0){
      tow.setHwPt(it->e5x5());
      tow.setHwEta(it->hwEta());
      tow.setHwPhi(it->hwPhi());
      towers_.push_back(tow);
    } else {
      tow.setHwPt(it->e5x5()+tow.hwPt());
    }
  }
  
  // HCAL TPs
  if (!hcalDigis_.empty()) {
    iSetup.get<CaloTPGRecord>().get(decoder_);
    edm::Handle<HcalTrigPrimDigiCollection> hcalTPs;
    for (const auto & token : hcalDigis_) {
      iEvent.getByToken(token, hcalTPs);
      for (const auto & itr : *hcalTPs) {
	HcalTrigTowerDetId id = itr.id();
	double et = decoder_->hcaletValue(itr.id(), itr.t0());
	//if (et > 2) std::cout << "HB Et = " << et << std::endl;
	if (et <= 0) continue;
	l1t::CaloTower tow = l1t::CaloTools::getTower(towers_, l1t::CaloTools::caloEta(id.ieta()), id.iphi());
	if(tow.hwPt()==0){
	  tow.setHwPt(et);
	  tow.setHwEta(id.ieta());
	  tow.setHwPhi(id.iphi());
	  towers_.push_back(tow);
	} else {
	  tow.setHwPt(et+tow.hwPt());
	}
      }
    }
  } 

  // HGC TPs
  edm::Handle<l1t::HGCalTowerBxCollection> hgCalTowers;
  iEvent.getByToken(hgCalColl_, hgCalTowers);
  for(auto it = hgCalTowers->begin(), end = hgCalTowers->end(); it != end; ++it) {
    //if (it->etHad() > 2) std::cout << "HGC Had Et = " << it->etHad() << std::endl;
    //if (it->etEm() > 2) std::cout << "HGC Em Et = " << it->etEm() << std::endl;
    l1t::CaloTower tow = l1t::CaloTools::getTower(towers_, l1t::CaloTools::caloEta(it->hwEta()), it->hwPhi());
    if(tow.hwPt()==0){
      tow.setHwPt(it->etHad()+it->etEm());
      tow.setHwEta(it->hwEta());
      tow.setHwPhi(it->hwPhi());
      towers_.push_back(tow);
    } else {
      tow.setHwPt(it->etHad()+it->etEm()+tow.hwPt());
    }
  }

  

}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
L1TCaloEtSumProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
L1TCaloEtSumProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
L1TCaloEtSumProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
L1TCaloEtSumProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
L1TCaloEtSumProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
L1TCaloEtSumProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TCaloEtSumProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TCaloEtSumProducer);
