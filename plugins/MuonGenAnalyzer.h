#ifndef MuonAnalysis_MuonAnalyzer_plugins_MuonGenAnalyzer
#define MuonAnalysis_MuonAnalyzer_plugins_MuonGenAnalyzer

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "NtupleContent.h"
#include "TLorentzVector.h"
#include "helper.h"

class MuonGenAnalyzer {
 public:
  MuonGenAnalyzer();
  virtual ~MuonGenAnalyzer();

  void SetInputs(const edm::Event &,
                 const edm::EDGetTokenT<reco::GenParticleCollection> &,
                 const int &);
  void FillNtuple(NtupleContent &);
  std::vector<reco::GenParticleRef> gmuons;

 private:
  edm::Handle<reco::GenParticleCollection> gens;
};

#endif
