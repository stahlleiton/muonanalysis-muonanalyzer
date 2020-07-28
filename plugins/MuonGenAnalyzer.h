#ifndef MuonAnalysis_MuonAnalyzer_plugins_MuonGenAnalyzer
#define MuonAnalysis_MuonAnalyzer_plugins_MuonGenAnalyzer

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "NtupleContent.h"
#include "helper.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


class MuonGenAnalyzer{
  public:
    MuonGenAnalyzer();
    virtual ~MuonGenAnalyzer();

    void SetInputs(const edm::Event &,const edm::EDGetTokenT<edm::View<reco::GenParticle> > & , const int & );
    void FillNtuple(NtupleContent &);

  private:  
    edm::Handle<edm::View<reco::GenParticle> > gens;
    TLorentzVector gmuon1,gmuon2;
    int gcharge1,gcharge2;
//    unsigned reco_idx1,reco_idx2;
    bool success=true;
};

#endif
