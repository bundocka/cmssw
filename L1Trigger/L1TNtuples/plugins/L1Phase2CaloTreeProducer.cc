// -*- C++ -*-
//
// Package:    L1TriggerDPG/L1Ntuples
// Class:      L1Phase2CaloTreeProducer
// 
/**\class L1Phase2CaloTreeProducer L1Phase2CaloTreeProducer.cc L1TriggerDPG/L1Ntuples/src/L1Phase2CaloTreeProducer.cc

Description: Produce L1 Extra tree

Implementation:
     
*/
//
// Original Author:  
//         Created:  
// $Id: L1Phase2CaloTreeProducer.cc,v 1.8 2012/08/29 12:44:03 jbrooke Exp $
//
//


// system include files
#include <memory>

// framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// data formats
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

// ROOT output stuff
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"

#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1Phase2Calo.h"

//
// class declaration
//

class L1Phase2CaloTreeProducer : public edm::EDAnalyzer {
public:
  explicit L1Phase2CaloTreeProducer(const edm::ParameterSet&);
  ~L1Phase2CaloTreeProducer() override;
  
  
private:
  void beginJob(void) override ;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

public:
  
  L1Analysis::L1AnalysisL1Phase2Calo* l1Phase2Calo;
  L1Analysis::L1AnalysisL1Phase2CaloDataFormat * l1Phase2CaloData;

private:

  unsigned maxL1Phase2Calo_;

  // output file
  edm::Service<TFileService> fs_;
  
  // tree
  TTree * tree_;
 
  // EDM input tags
  edm::EDGetTokenT<l1t::EtSumBxCollection> sumToken_;
  
};



L1Phase2CaloTreeProducer::L1Phase2CaloTreeProducer(const edm::ParameterSet& iConfig)
{

  sumToken_ = consumes<l1t::EtSumBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("sumToken"));
  
  maxL1Phase2Calo_ = iConfig.getParameter<unsigned int>("maxL1Phase2Calo");
 
  l1Phase2Calo     = new L1Analysis::L1AnalysisL1Phase2Calo();
  l1Phase2CaloData = l1Phase2Calo->getData();
  
  // set up output
  tree_=fs_->make<TTree>("L1Phase2CaloTree", "L1Phase2CaloTree");
  tree_->Branch("L1Phase2Calo", "L1Analysis::L1AnalysisL1Phase2CaloDataFormat", &l1Phase2CaloData, 32000, 3);

}


L1Phase2CaloTreeProducer::~L1Phase2CaloTreeProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
L1Phase2CaloTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  l1Phase2Calo->Reset();

  edm::Handle<l1t::EtSumBxCollection> sums;
  
  iEvent.getByToken(sumToken_, sums);
  
  if (sums.isValid()){ 
    l1Phase2Calo->SetSum(sums, maxL1Phase2Calo_);  
  } else {
    edm::LogWarning("MissingProduct") << "L1Phase2Calo EtSums not found. Branch will not be filled" << std::endl;
  }

  tree_->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void 
L1Phase2CaloTreeProducer::beginJob(void)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1Phase2CaloTreeProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1Phase2CaloTreeProducer);
