#include <L1Trigger/VertexFinder/interface/VertexProducer.h>

#include <iostream>
#include <set>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/L1TVertex/interface/Vertex.h"

#include "L1Trigger/VertexFinder/interface/VertexFinder.h"

#include "L1Trigger/VertexFinder/interface/RecoVertexWithTP.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

using namespace l1tVertexFinder;
using namespace std;


VertexProducer::VertexProducer(const edm::ParameterSet& iConfig) :
  l1TracksToken_(consumes<TTTrackCollectionView>(iConfig.getParameter<edm::InputTag>("l1TracksInputTag"))),
  ttTrackMCTruthToken_(consumes< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >(iConfig.getParameter<edm::InputTag>("mcTruthTrackInputTag"))),
  settings_(AlgoSettings(iConfig))
{
  // Get configuration parameters

  switch (settings_.vx_algo()) {
    case Algorithm::GapClustering:
      cout << "L1T vertex producer: Finding vertices using a gap clustering algorithm " << endl;
      break;
    case Algorithm::AgglomerativeHierarchical:
      cout << "L1T vertex producer: Finding vertices using a Simple Merge Clustering algorithm " << endl;
      break;
    case Algorithm::DBSCAN:
      cout << "L1T vertex producer: Finding vertices using a DBSCAN algorithm " << endl;
      break;
    case Algorithm::PVR:
      cout << "L1T vertex producer: Finding vertices using a PVR algorithm " << endl;
      break;
    case Algorithm::AdaptiveVertexReconstruction:
      cout << "L1T vertex producer: Finding vertices using an AdaptiveVertexReconstruction algorithm " << endl;
      break;
    case Algorithm::HPV:
      cout << "L1T vertex producer: Finding vertices using an Highest Pt Vertex algorithm " << endl;
      break;
    case Algorithm::Kmeans:
      cout << "L1T vertex producer: Finding vertices using a kmeans algorithm" << endl;
      break;
    case Algorithm::Generator:
      cout << "\n\n ** L1T vertex producer: Using ** GENERATOR ** vertex (average of TP z0s) ** \n\n" << endl;
      break;
  }

  // Tame debug printout.
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(4);

  //--- Define EDM output to be written to file (if required)
  produces<l1t::VertexCollection>("l1vertices");
  produces<l1t::VertexCollection>("l1vertextdr");

  if(settings_.vx_cnn_trk_assoc()){
    std::cout << "loading graph from " << settings_.vx_cnn_graph() << std::endl;
    // load the graph
    cnnGraph_ = tensorflow::loadGraphDef(settings_.vx_cnn_graph());
    // create a new session and add the graphDef
    cnnSesh_ = tensorflow::createSession(cnnGraph_);
  }
}

void VertexProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}

void VertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<TTTrackCollectionView> l1TracksHandle;
  iEvent.getByToken(l1TracksToken_, l1TracksHandle);

  std::vector<L1Track> l1Tracks;
  l1Tracks.reserve(l1TracksHandle->size());

  for (const auto& track : l1TracksHandle->ptrs())
    l1Tracks.push_back(L1Track(track));

  std::vector<const L1Track*> l1TrackPtrs;
  l1TrackPtrs.reserve(l1Tracks.size());
  for (const auto& track : l1Tracks) {
  //if (track.pt() < 100 or track.getNumStubs() > 5)
    if ((track.pt() > 2 && track.pt() < 500 &&
      abs(track.z0()) < 15 && track.chi2dof() < 100 && track.getNumStubs() > 4))
     l1TrackPtrs.push_back(&track);
  }
  //}

  // FIXME: Check with Davide if the tracks should be filtered using the following cuts
  //   fittedTracks[i].second.accepted() and fittedTracks[i].second.chi2dof()< settings_->chi2OverNdfCut()
  VertexFinder vf(l1TrackPtrs, settings_);

  switch (settings_.vx_algo()) {
    case Algorithm::GapClustering:
      vf.GapClustering();
      break;
    case Algorithm::AgglomerativeHierarchical:
      vf.AgglomerativeHierarchicalClustering();
      break;
    case Algorithm::DBSCAN:
      vf.DBSCAN();
      break;
    case Algorithm::PVR:
      vf.PVR();
      break;
    case Algorithm::AdaptiveVertexReconstruction:
      vf.AdaptiveVertexReconstruction();
      break;
    case Algorithm::HPV:
      vf.HPV();
      break;
    case Algorithm::Kmeans:
      vf.Kmeans();
      break;
    case Algorithm::Generator:
      std::vector<const L1Track*> pvTracks;
      for (const auto& track : l1Tracks) {
        edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTTrackHandle;
        iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);
        edm::Ptr< TrackingParticle > tpMatch = MCTruthTTTrackHandle->findTrackingParticlePtr(track.getTTTrackPtr());
        if(tpMatch.isNull())
          continue;
        if(tpMatch->eventId().event() == 0)
          pvTracks.push_back(&track);
      }
      vf.Generator(pvTracks);
      break;
  }

  vf.TDRalgorithm();
  vf.SortVerticesInZ0();
  vf.FindPrimaryVertex();

  // //=== Store output EDM track and hardware stub collections.
  std::unique_ptr<l1t::VertexCollection> lProduct(new std::vector<l1t::Vertex>());

  for (const auto& vtx : vf.Vertices()) {
    std::vector<edm::Ptr<l1t::Vertex::Track_t>> lVtxTracks;
    lVtxTracks.reserve(vtx.tracks().size());
    for (const auto& t : vtx.tracks())
      lVtxTracks.push_back(t->getTTTrackPtr());
    lProduct->emplace_back(l1t::Vertex(vtx.z0(), lVtxTracks));
  }
  iEvent.put(std::move(lProduct), "l1vertices");

  // //=== Store output EDM track and hardware stub collections.
  std::unique_ptr<l1t::VertexCollection> lProductTDR(new std::vector<l1t::Vertex>());
  std::vector<edm::Ptr<l1t::Vertex::Track_t>> lVtxTracksTDR;
  lVtxTracksTDR.reserve(vf.TDRPrimaryVertex().tracks().size());
  // use normal tracks or cnn tracks
  if(settings_.vx_cnn_trk_assoc()){
    std::vector<const L1Track*> cnnPVTracks;
    vf.cnnTrkAssociation(vf.TDRPrimaryVertex().z0(), cnnPVTracks, cnnSesh_);
    for (const auto& t : cnnPVTracks)
      lVtxTracksTDR.emplace_back(t->getTTTrackPtr());
  } else {
    for (const auto& t : vf.TDRPrimaryVertex().tracks())
      lVtxTracksTDR.emplace_back(t->getTTTrackPtr());
  }
  lProductTDR->emplace_back(l1t::Vertex(vf.TDRPrimaryVertex().z0(), lVtxTracksTDR));
  iEvent.put(std::move(lProductTDR), "l1vertextdr");
}

void VertexProducer::endJob()
{

  if(settings_.vx_cnn_trk_assoc()){
    // close the session
    tensorflow::closeSession(cnnSesh_);
    cnnSesh_ = nullptr;

    // delete the graph
    delete cnnGraph_;
    cnnGraph_ = nullptr;
  }
}

DEFINE_FWK_MODULE(VertexProducer);
