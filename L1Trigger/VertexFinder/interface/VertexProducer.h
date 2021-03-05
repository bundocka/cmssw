#ifndef __L1Trigger_VertexFinder_VertexProducer_h__
#define __L1Trigger_VertexFinder_VertexProducer_h__


#include <map>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"

#include "L1Trigger/VertexFinder/interface/AlgoSettings.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"


namespace l1tVertexFinder {
class AlgoSettings;
}

class VertexProducer : public edm::EDProducer {

public:
  explicit VertexProducer(const edm::ParameterSet&);
  ~VertexProducer() {}

private:
  typedef edm::View<TTTrack<Ref_Phase2TrackerDigi_>> TTTrackCollectionView;

  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

private:
  const edm::EDGetTokenT<TTTrackCollectionView> l1TracksToken_;
  const edm::EDGetTokenT<TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > ttTrackMCTruthToken_;

  tensorflow::GraphDef* cnnAssGraph_;
  tensorflow::Session* cnnAssSesh_;
  tensorflow::GraphDef* cnnTrkGraph_;
  tensorflow::Session* cnnTrkSesh_;
  tensorflow::GraphDef* cnnPVZ0Graph_;
  tensorflow::Session* cnnPVZ0Sesh_;

  l1tVertexFinder::AlgoSettings settings_;
};

#endif
