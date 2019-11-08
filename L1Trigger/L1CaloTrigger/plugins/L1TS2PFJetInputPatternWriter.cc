// -*- C++ -*-
//
// Package:    L1Trigger/L1CaloTrigger
// Class:      L1TS2PFJetInputPatternWriter
// 
/**\class L1TS2PFJetInputPatternWriter L1TS2PFJetInputPatternWriter.cc L1Trigger/L1TCalorimeter/plugins/L1TS2PFJetInputPatternWriter.cc

   Description: 

   Implementation:

*/
//
// Original Author:  Aaron Bundock
//         Created:  Fri, 26 Jul 2018 14:20:25 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

//
// class declaration
//

class L1TS2PFJetInputPatternWriter : public edm::EDAnalyzer {
public:
  explicit L1TS2PFJetInputPatternWriter(const edm::ParameterSet&);
  ~L1TS2PFJetInputPatternWriter() override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  uint64_t convertEtaOrPhi(double etaOrPhi, double centre);
  
private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<l1t::PFCandidate> > pfToken_;
  std::string filename_;
  std::string outDir_;

  // constants
  unsigned nChan_;  // number of channels per quad
  unsigned nQuad_;
  unsigned nLink_;
  unsigned nHeaderFrames_;
  unsigned nPayloadFrames_;
  unsigned nClearFrames_;
  unsigned nFrame_;
  unsigned nFrameFile_;
  unsigned nEvents_;
  unsigned nFramesPerFile_;
  float    ptLSB_;
  float    etaPhiLSB_;
  
  // data arranged by link and frame
  std::vector< std::vector<uint64_t> > data_;

  // data valid flags (just one per frame for now)
  std::vector<int> dataValid_;

  // map of towers onto links/frames
  std::map< int, int > map_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1TS2PFJetInputPatternWriter::L1TS2PFJetInputPatternWriter(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed

  // register what you consume and keep token for later access:
  pfToken_ = consumes<std::vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pfTag"));  
  filename_ = iConfig.getUntrackedParameter<std::string>("filename");
  outDir_ = iConfig.getUntrackedParameter<std::string>("outDir");

  nChan_ = 4;
  nQuad_ = 18;
  ptLSB_ = 0.25;
  etaPhiLSB_ = 0.0043633231;

  nHeaderFrames_ = iConfig.getUntrackedParameter<unsigned>("nHeaderFrames");
  nPayloadFrames_ = iConfig.getUntrackedParameter<unsigned>("nPayloadFrames");
  nClearFrames_ = iConfig.getUntrackedParameter<unsigned>("nClearFrames");
  nFrame_ = 0;
  nFrameFile_ = 0;
  nEvents_ = 0;
  nFramesPerFile_ = 977;

  nLink_ = nChan_ * nQuad_;
  data_.resize(nLink_);
  LogDebug("L1TDebug") << "Preparing for " << nLink_ << " links" << std::endl;

}


L1TS2PFJetInputPatternWriter::~L1TS2PFJetInputPatternWriter()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1TS2PFJetInputPatternWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  //count events
  nEvents_++;

  edm::Handle<std::vector<l1t::PFCandidate>> pfHandle;
  iEvent.getByToken(pfToken_, pfHandle);

  std::vector<l1t::PFCandidate> pfPartsA;
  std::vector<l1t::PFCandidate> pfPartsB;
  
  for (std::vector<l1t::PFCandidate>::const_iterator pfIt = pfHandle->begin(); pfIt != pfHandle->end(); pfIt++){
    // select first two "small" regions for current fw
    if(pfIt->eta() >= 0 && pfIt->eta() < 0.75 && pfIt->phi() >= 0 && pfIt->phi() < 0.7)
      pfPartsA.push_back(*pfIt);
    if(pfIt->eta() >= 0.75 && pfIt->eta() < 1.5 && pfIt->phi() >= 0 && pfIt->phi() < 0.7)
      pfPartsB.push_back(*pfIt);
  }

  //if(pfPartsA.size()==0 && pfPartsB.size()==0)
  //  return;
  

  if(nFrame_ == 0 || nFrameFile_ == 0){
    //first empty frames
    while(nFrameFile_ < 8){
      if(nFrameFile_ > 4)
	dataValid_.push_back( 1 );
      else
	dataValid_.push_back( 0 );
      for ( unsigned iQuad=0; iQuad<nQuad_; ++iQuad ) {
	for ( unsigned iChan=0; iChan<nChan_; ++iChan ) {
	  uint iLink = (iQuad*nChan_)+iChan;
	  if(iLink==40 and nFrameFile_ == 5)
	    data_.at(iLink).push_back(5839227683695833246);
	  else
	    data_.at(iLink).push_back(0);
	  continue;
	}
      }
      nFrame_++;
      nFrameFile_++;
    }    
  }  

  // loop over frames
  for ( unsigned iFrame=0; iFrame<nPayloadFrames_; ++iFrame ) {
    dataValid_.push_back( 1 );
    // loop over links
    for ( unsigned iQuad=0; iQuad<nQuad_; ++iQuad ) {
      for ( unsigned iChan=0; iChan<nChan_; ++iChan ) {

        // get tower ieta, iphi for link
	uint iLink = (iQuad*nChan_)+iChan;

	uint64_t data=0;     

	if((nFrameFile_%17) == 8){
	  if(iLink > 39 && pfPartsA.size() > (iLink-40)){
	    data |= ((uint64_t)floor(pfPartsA.at(iLink-40).pt()  / ptLSB_ )     & 0xffff) << 32;
	    data |= convertEtaOrPhi(pfPartsA.at(iLink-40).eta(), 0.375);
	    data |= convertEtaOrPhi(pfPartsA.at(iLink-40).phi(), 0.35)  << 10;
	    //std::cout << std::fixed << std::setprecision(2) << pfPartsA.at(iLink).pt() << "\t" <<  
	    // pfPartsA.at(iLink).eta() << "\t" << pfPartsA.at(iLink).phi() << std::endl;

	  }
	}
	if((nFrameFile_%17) == 10){
	  if(iLink > 39 && pfPartsB.size() > (iLink-40)){
	    data |= ((uint64_t)floor(pfPartsB.at(iLink-40).pt()  / ptLSB_ )     & 0xffff) << 32;
	    data |= convertEtaOrPhi(pfPartsB.at(iLink-40).eta(), 1.125);
	    data |= convertEtaOrPhi(pfPartsB.at(iLink-40).phi(), 0.35)  << 10;
	    //std::cout << std::fixed << std::setprecision(2) << pfPartsB.at(iLink).pt() << "\t" <<  
	    //  pfPartsB.at(iLink).eta() << "\t" << pfPartsB.at(iLink).phi() << std::endl;
	  }
	}
	// add data to output
	data_.at(iLink).push_back( data );
      }
    }
    nFrame_+=1;
    nFrameFile_+=1;
    if(nFrame_%nFramesPerFile_ == 0) nFrameFile_ = 0;
  }
}

uint64_t
  L1TS2PFJetInputPatternWriter::convertEtaOrPhi(double etaOrPhi, double centre){

  int sign = (etaOrPhi-centre)/abs(etaOrPhi-centre);
  uint64_t ietaOrPhi = (uint64_t)floor((abs(etaOrPhi-centre)/etaPhiLSB_));
  //double etaOrPhiRe = ((ietaOrPhi*etaPhiLSB_)*sign)+centre;


  if(sign<0){
    ietaOrPhi = (~ietaOrPhi & 0x3ff) +1;//1024-ietaOrPhi;//(~ietaOrPhi & 0x3ff) +1;
    //etaOrPhiRe = (((~(ietaOrPhi-1) & 0x3ff) * etaPhiLSB_)*sign)+centre;//(((1024-ietaOrPhi)*etaPhiLSB_)*sign)+centre;//(((~(ietaOrPhi-1) & 0x1ff) * etaPhiLSB_)*sign)+centre;
  }

  //std::cout << "etaOrPhi = " << etaOrPhi << ", rel etaOrPhi = " << etaOrPhi - centre << std::endl;
  //std::cout << "etaOrPhiRe = " << etaOrPhiRe << ", rel etaOrPhiRe = " << etaOrPhiRe - centre << std::endl;


  return ietaOrPhi;

}


// ------------ method called once each job just before starting event loop  ------------
void 
L1TS2PFJetInputPatternWriter::beginJob()
{


}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1TS2PFJetInputPatternWriter::endJob() 
{

  //frames per event
  unsigned int framesPerEv = nHeaderFrames_ + nPayloadFrames_ + nClearFrames_;

  //events per file
  unsigned int evPerFile = floor(nFramesPerFile_/framesPerEv);

  //number of output files
  unsigned int nOutFiles = ceil((float)nEvents_/(float)evPerFile);

  LogDebug("L1TDebug") << "Read " << nFrame_ << " frames" << std::endl;
  LogDebug("L1TDebug") << "Read " << nEvents_ << " events" << std::endl;
  LogDebug("L1TDebug") << "Writing " << nOutFiles << " files" << std::endl;
  LogDebug("L1TDebug") << "Output directory: ./" << outDir_ << "/" << std::endl;

  //files
  std::vector< std::ofstream > outFiles(nOutFiles);
    
  //make output files and write to them
  for(uint itFile=0; itFile<nOutFiles; ++itFile){
    std::stringstream outFilename;
    outFilename << outDir_ << "/" << filename_ << "_" << itFile << ".txt";
    outFiles[itFile] = std::ofstream(outFilename.str());
    LogDebug("L1TDebug") << "Writing to file: ./" << outFilename.str() << std::endl;
    std::cout << "Writing to file: ./" << outFilename.str() << std::endl;

    outFiles[itFile] << "Board SRNTY_TEST" << std::endl;
    
    // quad/chan numbers
    outFiles[itFile] << " Quad/Chan :      ";
    for ( unsigned i=0; i<nQuad_; ++i ) {
      for ( unsigned j=0; j<nChan_; ++j ) {
	outFiles[itFile] << "  q" << setfill('0') << setw(2) << i << "c" << j << "            ";
      }
    }
    outFiles[itFile] << std::endl;

    // link numbers
    outFiles[itFile] << "      Link :     ";
    for ( unsigned i=0; i<nQuad_; ++i ) {
      for ( unsigned j=0; j<nChan_; ++j ) {
	outFiles[itFile] << "    " << setfill('0') << setw(2) << (i*nChan_)+j << "             ";
      }
    }

    outFiles[itFile] << std::endl;

    // then the data
    unsigned iFileFrame=0;
    for ( unsigned iFrame=itFile*nFramesPerFile_; iFrame<(itFile*nFramesPerFile_+  nFramesPerFile_); ++iFrame ) {
      if( iFrame <= nFrame_  && iFrame < (framesPerEv*nEvents_)){
	outFiles[itFile] << "Frame " << std::dec << std::setw(4) << std::setfill('0') << iFileFrame << " : ";
	for ( unsigned iQuad=0; iQuad<nQuad_; ++iQuad ) {
	  for ( unsigned iChan=0; iChan<nChan_; ++iChan ) {
	    unsigned iLink = (iQuad*nChan_)+iChan;
	    if (iLink<data_.size() && iFrame<data_.at(iLink).size()) {
	      outFiles[itFile] << std::hex << ::std::setw(1) << dataValid_.at(iFrame) << "v" << std::hex << std::setw(16) << std::setfill('0') << data_.at(iLink).at(iFrame) << " ";
	    }
	    else {
	      //std::cerr << "Out of range : " << iLink << ", " << iFrame << std::endl;
	      outFiles[itFile] << std::hex << ::std::setw(1) << 0 << "v" << std::hex << std::setw(16) << std::setfill('0') << 0 << " ";
	    }
	  }
	}
      }
      outFiles[itFile] << std::endl;
      iFileFrame++;
    }
    outFiles[itFile].close();
  }
  
}



// ------------ method called when starting to processes a run  ------------
/*
  void 
  L1TS2PFJetInputPatternWriter::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  L1TS2PFJetInputPatternWriter::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  L1TS2PFJetInputPatternWriter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  L1TS2PFJetInputPatternWriter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TS2PFJetInputPatternWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TS2PFJetInputPatternWriter);
