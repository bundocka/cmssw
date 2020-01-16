#ifndef L1Trigger_Phase2L1ParticleFlow_PUAlgoBase_h
#define L1Trigger_Phase2L1ParticleFlow_PUAlgoBase_h

#include "L1Trigger/Phase2L1ParticleFlow/interface/Region.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1TVertex/interface/Vertex.h"

namespace l1tpf_impl { 

  class PUAlgoBase {
    public:
        PUAlgoBase( const edm::ParameterSet& ) ;
        virtual ~PUAlgoBase() ;

        /// global operations
        enum VertexAlgo { OldVtxAlgo, TPVtxAlgo, ExternalVtxAlgo, ExternalCNNVtxAlgo };
        virtual void doVertexing(std::vector<Region> &rs, VertexAlgo algo, float& pvdz, l1t::Vertex& pv) const ; // region is not const since it sets the fromPV bit of the tracks

        virtual void runChargedPV(Region &r, float z0) const ;
        
        virtual const std::vector<std::string> & puGlobalNames() const ;
        virtual void doPUGlobals(const std::vector<Region> &rs, float npu, std::vector<float> & globals) const = 0;
        virtual void runNeutralsPU(Region &r, float npu, const std::vector<float> & globals) const = 0;

    protected:
        int debug_;
        float etaCharged_, vtxRes_;
        bool vtxAdaptiveCut_; 
  };

} // end namespace

#endif
