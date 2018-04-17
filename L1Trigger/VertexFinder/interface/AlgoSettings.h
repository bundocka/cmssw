#ifndef __L1Trigger_VertexFinder_AlgoSettings_h__
#define __L1Trigger_VertexFinder_AlgoSettings_h__


#include <vector>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"



namespace l1tVertexFinder {

enum class Algorithm {
  GapClustering,
  AgglomerativeHierarchical,
  DBSCAN,
  PVR,
  AdaptiveVertexReconstruction,
  HPV,
  Kmeans
};


class AlgoSettings {

public:
  AlgoSettings(const edm::ParameterSet& iConfig);
  ~AlgoSettings() {}

  //=== Vertex Reconstruction configuration
  // Vertex Reconstruction algo
  Algorithm vx_algo() const { return vx_algo_; }
  /// For Agglomerative cluster algorithm, select a definition of distance between clusters
  unsigned int vx_distanceType() const { return vx_distanceType_; }
  // Assumed Vertex Distance
  float vx_distance() const { return vx_distance_; }
  // Assumed Vertex Resolution
  float vx_resolution() const { return vx_resolution_; }
  // Minimum number of tracks to accept vertex
  unsigned int vx_minTracks() const { return vx_minTracks_; }
  // Compute the z0 position of the vertex with a mean weighted with track momenta
  bool vx_weightedmean() const { return vx_weightedmean_; }
  /// Chi2 cut for the Adaptive Vertex Recostruction Algorithm
  float vx_chi2cut() const { return vx_chi2cut_; }
  /// TDR assumed vertex width
  float tdr_vx_width() const { return tdr_vx_width_; }
  float vx_dbscan_pt() const { return vx_dbscan_pt_; }
  unsigned int vx_dbscan_mintracks() const { return vx_dbscan_mintracks_; }

  unsigned int vx_kmeans_iterations() const { return vx_kmeans_iterations_; }
  unsigned int vx_kmeans_nclusters() const { return vx_kmeans_nclusters_; }
  float vx_TrackMinPt() const { return vx_TrackMinPt_; }



  //=== Debug printout
  unsigned int debug() const { return debug_; }


  //=== Hard-wired constants
  // EJC Check this.  Found stub at r = 109.504 with flat geometry in 81X, so increased tracker radius for now.
  double trackerOuterRadius() const { return 120.2; } // max. occuring stub radius.
  // EJC Check this.  Found stub at r = 20.664 with flat geometry in 81X, so decreased tracker radius for now.
  double trackerInnerRadius() const { return 20; }   // min. occuring stub radius.
  double trackerHalfLength() const { return 270.; }  // half-length  of tracker.
  double layerIDfromRadiusBin() const { return 6.; } // When counting stubs in layers, actually histogram stubs in distance from beam-line with this bin size.

private:
  static const std::map<std::string, Algorithm> algoNameMap;

  // Parameter sets for differents types of configuration parameter.
  edm::ParameterSet vertex_;

  // Vertex Reconstruction configuration
  Algorithm vx_algo_;
  float vx_distance_;
  float vx_resolution_;
  unsigned int vx_distanceType_;
  unsigned int vx_minTracks_;
  bool vx_weightedmean_;
  float vx_chi2cut_;
  float tdr_vx_width_;
  float vx_TrackMinPt_;
  float vx_dbscan_pt_;
  float vx_dbscan_mintracks_;
  unsigned int vx_kmeans_iterations_;
  unsigned int vx_kmeans_nclusters_;

  // Debug printout
  unsigned int debug_;
};

} // end namespace l1tVertexFinder

#endif
