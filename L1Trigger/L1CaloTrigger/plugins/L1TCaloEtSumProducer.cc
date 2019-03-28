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
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/CaloClusterer.h"


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

  edm::EDGetTokenT<l1slhc::L1EGCrystalClusterCollection> ecalColl_;
  double ecalEtMin_;

  std::vector<edm::EDGetTokenT<HcalTrigPrimDigiCollection>> hcalDigis_;
  edm::ESHandle<CaloTPGTranscoder> decoder_;

  edm::EDGetTokenT<l1t::HGCalTowerBxCollection> hgCalColl_;
  edm::EDGetTokenT<l1t::HGCalMulticlusterBxCollection> hgCalMClusts_;

  std::vector<l1t::EtSum> etsums_;

  int metEtaMax_;
  int metEtaMaxHF_;
  int ettEtaMax_;
  int ettEtaMaxHF_;

};

L1TCaloEtSumProducer::L1TCaloEtSumProducer(const edm::ParameterSet& iConfig) :
  ecalColl_(consumes<l1slhc::L1EGCrystalClusterCollection>(iConfig.getParameter<edm::InputTag>("ecalColl"))),
  ecalEtMin_(iConfig.getParameter<double>("ecalEtMin")),
  hgCalColl_(consumes<l1t::HGCalTowerBxCollection>(iConfig.getParameter<edm::InputTag>("hgCalColl"))),
  hgCalMClusts_(consumes<l1t::HGCalMulticlusterBxCollection>(iConfig.getParameter<edm::InputTag>("hgCalMClusts"))),
  metEtaMax_(iConfig.getParameter<int>("metEtaMax")),
  metEtaMaxHF_(iConfig.getParameter<int>("metEtaMaxHF")),
  ettEtaMax_(iConfig.getParameter<int>("ettEtaMax")),
  ettEtaMaxHF_(iConfig.getParameter<int>("ettEtaMaxHF"))
{

  //register your products
  produces<l1t::EtSumBxCollection>();

  for (auto & tag : iConfig.getParameter<std::vector<edm::InputTag>>("hcalDigis")) {
    hcalDigis_.push_back(consumes<HcalTrigPrimDigiCollection>(tag));
  }

}


L1TCaloEtSumProducer::~L1TCaloEtSumProducer()
{
}


void
L1TCaloEtSumProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  etsums_.clear();

  math::XYZTLorentzVector p4;

  double exHF(0), eyHF(0), etHF(0);
  double ex(0), ey(0), et(0), etEm(0);
  unsigned int mbP0(0), mbP1(0), mbM0(0), mbM1(0);
  
  //bool ettSat(false), ettHFSat(false), ecalEtSat(false), metSat(false), metHFSat(false);

  
  // ECAL clusters
  edm::Handle<l1slhc::L1EGCrystalClusterCollection> ecalClusters;
  iEvent.getByToken(ecalColl_, ecalClusters);
  for(auto it = ecalClusters->begin(), end = ecalClusters->end(); it != end; ++it) {
    //if (it->calibratedPt() < ecalEtMin_) continue;
    //if (it->calibratedPt() <= 0.5) continue;
    //if(it->calibratedPt() > 2) std::cout << "calibratedPt = " << it->calibratedPt() << ", eta = " << it->eta() << ", phi = " << it->phi() << std::endl;
    ex += it->calibratedPt()*cos(it->phi());
    ey += it->calibratedPt()*sin(it->phi());
    et += it->calibratedPt();
    etEm += it->calibratedPt();
    exHF += it->calibratedPt()*cos(it->phi());
    eyHF += it->calibratedPt()*sin(it->phi());
    etHF += it->calibratedPt();
  }
  
  // HCAL TPs
  if (!hcalDigis_.empty()) {
    iSetup.get<CaloTPGRecord>().get(decoder_);
    edm::Handle<HcalTrigPrimDigiCollection> hcalTPs;
    for (const auto & token : hcalDigis_) {
      iEvent.getByToken(token, hcalTPs);
      for (const auto & itr : *hcalTPs) {
	HcalTrigTowerDetId id = itr.id();
	double tpEt = decoder_->hcaletValue(itr.id(), itr.t0());
	double tpEta = l1t::CaloTools::towerEta(id.ieta());
	double tpPhi = l1t::CaloTools::towerPhi(id.ieta(), id.iphi());
	if (tpEt <= 0) continue;
	//if(tpEt>2) std::cout << "HCAL TP Et = " << tpEt << ", eta = " << id.ieta() << ", phi = " << id.iphi() << std::endl;
	if(abs(tpEta) < 3.0){
	  ex += tpEt*cos(tpPhi);
	  ey += tpEt*sin(tpPhi);
	  et += tpEt;
	  exHF += tpEt*cos(tpPhi);
	  eyHF += tpEt*sin(tpPhi);
	  etHF += tpEt;
	  
	  //if(et>2) std::cout << "Tow EtHad = " << tow.etHad() << ", eta = " << tow.hwEta() << ", phi = " << tow.hwPhi() << std::endl;
	}else{
	  exHF += tpEt*cos(tpPhi);
	  eyHF += tpEt*sin(tpPhi);
	  etHF += tpEt;
	  if(itr.t0().fineGrain(1)) {
	    if(tpEta > 3.0) mbP0 += 1;
	    if(tpEta <-3.0) mbM0 += 1;
	  }
	}
      }
    }
  } 
  
  // HGC TPs

  // HGCAL towers
  //edm::Handle<l1t::HGCalTowerBxCollection> hgCalTowers;
  //iEvent.getByToken(hgCalColl_, hgCalTowers);
  //for(auto it = hgCalTowers->begin(), end = hgCalTowers->end(); it != end; ++it) {
    //if(it->etEm() > 2) std::cout << "HGC Had Et = " << it->etHad() << ", Em Et = " << it->etEm() << ", eta = " << it->eta() << ", phi = " << it->phi() << std::endl;
    //if(it->etHad() > 2) std::cout << "HGC Had Et = " << it->etHad() << ", Em Et = " << it->etEm() << ", eta = " << it->eta() << ", phi = " << it->phi() << std::endl;
    //if(it->etEm() < 0.5 && it->etHad() < 0.5) continue;
    //double tpEt = it->etHad() + it->etEm();
    //ex   += tpEt*cos(it->phi());
    //ey   += tpEt*sin(it->phi());
    //et   += tpEt;
    //exHF   += tpEt*cos(it->phi());
    //eyHF   += tpEt*sin(it->phi());
    //etHF   += tpEt;
    //etEm += it->etEm();
    //if(it->etEm() > 2)  std::cout << "Tow EtHad = " << tow.etHad() << ", tow EtEm = " << tow.etEm() << ", eta = " << tow.hwEta() << ", phi = " << tow.hwPhi() << std::endl;
    //if(it->etHad() > 2) std::cout << "Tow EtHad = " << tow.etHad() << ", tow EtEm = " << tow.etEm() << ", eta = " << tow.hwEta() << ", phi = " << tow.hwPhi() << std::endl;
  //}
 

  // HGCAL clusters
  edm::Handle<l1t::HGCalMulticlusterBxCollection> multiclusters;
  iEvent.getByToken(hgCalMClusts_, multiclusters);

  for(auto it = multiclusters->begin(0), ed = multiclusters->end(0); it != ed; ++it) {
    float pt  = it->pt();
  //  float hoe = it->hOverE();
    ex   += pt*cos(it->phi());
    ey   += pt*sin(it->phi());
    et   += pt;
    exHF   += pt*cos(it->phi());
    eyHF   += pt*sin(it->phi());
    etHF   += pt;
    //  etEm += hoe == -1 ? 0 : pt /= 1 + hoe;
  }
   


  if (mbP0>0xf) mbP0 = 0xf; 
  if (mbP1>0xf) mbP1 = 0xf;
  if (mbM0>0xf) mbM0 = 0xf; 
  if (mbM1>0xf) mbM1 = 0xf;
  l1t::EtSum etSumMBP0(p4,l1t::EtSum::EtSumType::kMinBiasHFP0,mbP0,0,0,0);
  l1t::EtSum etSumMBP1(p4,l1t::EtSum::EtSumType::kMinBiasHFP1,mbP1,0,0,0);
  l1t::EtSum etSumMBM0(p4,l1t::EtSum::EtSumType::kMinBiasHFM0,mbM0,0,0,0);
  l1t::EtSum etSumMBM1(p4,l1t::EtSum::EtSumType::kMinBiasHFM1,mbM1,0,0,0);
  etsums_.push_back(etSumMBP0); 
  etsums_.push_back(etSumMBP1);
  etsums_.push_back(etSumMBM0); 
  etsums_.push_back(etSumMBM1);
   
   
  double met(0), metHF(0);
  double metPhi(0), metHFPhi(0);
   
  // Final MET calculation
  met = sqrt(ex*ex + ey*ey);
  metPhi = atan(ey/ex);
   
  // Final METHF calculation
  metHF = sqrt(exHF*exHF + eyHF*eyHF);
  metHFPhi = atan(eyHF/exHF);
   
  // Make final collection
  l1t::EtSum etSumTotalEt(p4,l1t::EtSum::EtSumType::kTotalEt,et,0,0,0);
  l1t::EtSum etSumTotalEtEm(p4,l1t::EtSum::EtSumType::kTotalEtEm,etEm,0,0,0);
  l1t::EtSum etSumMissingEt(p4,l1t::EtSum::EtSumType::kMissingEt,met,0,metPhi,0);
  l1t::EtSum etSumMissingEtHF(p4,l1t::EtSum::EtSumType::kMissingEtHF,metHF,0,metHFPhi,0);
   
  etsums_.push_back(etSumTotalEt);
  etsums_.push_back(etSumTotalEtEm);
  etsums_.push_back(etSumMissingEt);
  etsums_.push_back(etSumMissingEtHF);


  std::unique_ptr<l1t::EtSumBxCollection> etsums (new l1t::EtSumBxCollection(0, 0, 0));
  for( auto etsum = etsums_.begin(); etsum != etsums_.end(); ++etsum)
    etsums->push_back(0, l1t::CaloTools::etSumP4Demux(*etsum));
  
  iEvent.put(std::move(etsums));

  //std::cout << "Total ET = " << et << ", METHF = " << metHF << std::endl;
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
