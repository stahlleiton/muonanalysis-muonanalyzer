// Package:    MuonAnalyzer for Run 3
//             version 2.0
//
/**\class 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:
//                george karathanasis
//         Created:  Thu, 20 feb 2020 17:40:23 GMT
//
//

// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include <vector>
#include "TTree.h"
#include "DataFormats/Common/interface/Ref.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TLorentzVector.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "NtupleContent.h"
#include "helper.h"
#include "KlFitter.h"
#include "MuonBranches.h"
#include "MuonGenAnalyzer.h"

using namespace std;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
class MuonFullAODAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  typedef std::vector<std::pair<reco::Muon, reco::TransientTrack>> RecoTrkAndTransientTrkCollection;
  explicit MuonFullAODAnalyzer(const edm::ParameterSet&);
  ~MuonFullAODAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  bool HLTaccept(const edm::Event&, NtupleContent&, std::vector<std::string>&);
  void HLTmuon(const edm::Event&,
               std::vector<float>&,
               std::vector<float>&,
               std::vector<float>&,
               std::vector<std::string>&,
               const int&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  edm::EDGetToken muonsToken_;
  edm::EDGetToken tracksToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genToken_;
  std::vector<std::string> HLTPaths_;
  std::vector<std::string> HLTFilters_;
  std::vector<std::string> ProbePaths_;
  std::vector<std::string> ProbeFilters_;

  const double trgDRwindow_;
  const unsigned int tagQual_;
  const StringCutObjectSelector<reco::Muon> tagSelection_;  //kinematic cuts for tag
  const bool HighPurity_;
  const StringCutObjectSelector<reco::Track> probeSelection_;  //kinematic cuts for probe
  const double pairMassMin_;
  const double pairMassMax_;
  const double pairDz_;
  const bool RequireVtxCreation_;  //if true skip pairs that do not create gthat do not have a vertex
  const double minSVtxProb_;       //min probability of a vertex to be kept. If <0 inactive
  const double maxdz_trk_mu_;
  const double maxpt_relative_dif_trk_mu_;
  const double maxdr_trk_mu_;
  const unsigned momPdgId_;
  const double genRecoDrMatch_;
  const int debug_;

  edm::Service<TFileService> fs;
  TTree* t1;
  NtupleContent nt;
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
MuonFullAODAnalyzer::MuonFullAODAnalyzer(const edm::ParameterSet& iConfig)
    :  // inputs
      beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
      vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      muonsToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      trgresultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      trigobjectsToken_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
      genToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("gen"))),
      HLTPaths_(iConfig.getParameter<std::vector<std::string>>("triggerPaths")),
      HLTFilters_(iConfig.getParameter<std::vector<std::string>>("triggerFilters")),
      ProbePaths_(iConfig.getParameter<std::vector<std::string>>("ProbePaths")),
      ProbeFilters_(iConfig.getParameter<std::vector<std::string>>("ProbeFilters")),
      trgDRwindow_(iConfig.getParameter<double>("trgDRwindow")),
      tagQual_(iConfig.getParameter<unsigned>("tagQuality")),
      tagSelection_(iConfig.getParameter<std::string>("tagSelection")),
      HighPurity_(iConfig.getParameter<bool>("probeHPyrity")),
      probeSelection_(iConfig.getParameter<std::string>("probeSelection")),
      pairMassMin_(iConfig.getParameter<double>("pairMassMin")),
      pairMassMax_(iConfig.getParameter<double>("pairMassMax")),
      pairDz_(iConfig.getParameter<double>("pairDz")),
      RequireVtxCreation_(iConfig.getParameter<bool>("RequireVtxCreation")),
      minSVtxProb_(iConfig.getParameter<double>("minSVtxProb")),
      maxdz_trk_mu_(iConfig.getParameter<double>("maxDzProbeTrkMuon")),
      maxpt_relative_dif_trk_mu_(iConfig.getParameter<double>("maxRelPtProbeTrkMuon")),
      maxdr_trk_mu_(iConfig.getParameter<double>("maxDRProbeTrkMuon")),
      momPdgId_(iConfig.getParameter<unsigned>("momPdgId")),
      genRecoDrMatch_(iConfig.getParameter<double>("genRecoDrMatch")),
      debug_(iConfig.getParameter<int>("debug"))

{
  //  edm::ParameterSet runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");
}

MuonFullAODAnalyzer::~MuonFullAODAnalyzer() {
  // cout<<"total "<<trg_counter<<" fires "<<fire_counter<<" l3 "<<l3_counter<<endl;
  // do anything here that needs to be done at desctruction time
}

bool MuonFullAODAnalyzer::HLTaccept(const edm::Event& iEvent, NtupleContent& nt, std::vector<std::string>& HLTPaths) {
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  edm::TriggerNames trigName;
  trigName = iEvent.triggerNames(*trigResults);
  bool EvtFire = false;
  unsigned int ipath = 0;
  for (auto path : HLTPaths) {
    bool TrgFire = false;
    for (unsigned int itrg = 0; itrg < trigResults->size(); ++itrg) {
      TString TrigPath = trigName.triggerName(itrg);
      if (!trigResults->accept(itrg))
        continue;
      if (!TrigPath.Contains(path))
        continue;
      EvtFire = true;
      TrgFire = true;
    }
    nt.trigger[ipath] = TrgFire;
    ipath++;
  }
  return EvtFire;
}

void MuonFullAODAnalyzer::HLTmuon(const edm::Event& iEvent,
                                  std::vector<float>& trg_pt,
                                  std::vector<float>& trg_eta,
                                  std::vector<float>& trg_phi,
                                  std::vector<std::string>& HLTFilters,
                                  const int& debug_) {
  edm::Handle<trigger::TriggerEvent> triggerObjects;
  iEvent.getByToken(trigobjectsToken_, triggerObjects);
  trigger::TriggerObjectCollection allTriggerObjects = triggerObjects->getObjects();
  for (auto ifilter : HLTFilters) {
    size_t filterIndex = (*triggerObjects).filterIndex(edm::InputTag(ifilter, "", "HLT"));
    if (filterIndex < (*triggerObjects).sizeFilters()) {
      const trigger::Keys& keys = (*triggerObjects).filterKeys(filterIndex);
      for (size_t j = 0; j < keys.size(); j++) {
        trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
        if (fabs(foundObject.id()) != 13)
          continue;
        trg_pt.push_back(foundObject.pt());
        trg_eta.push_back(foundObject.eta());
        trg_phi.push_back(foundObject.phi());
        if (debug_ > 0)
          std::cout << "trg muon " << foundObject.pt() << std::endl;
      }
    }
  }
}

//
void MuonFullAODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  //Get data
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_, theBeamSpot);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  //skip evts if there are no vertices
  if (vertices->size() == 0)
    return;
  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);
  edm::Handle<std::vector<reco::Track>> tracks;
  iEvent.getByToken(tracksToken_, tracks);
  edm::ESHandle<MagneticField> bField;
  iSetup.get<IdealMagneticFieldRecord>().get(bField);

  //information about run
  nt.ClearBranches();
  nt.run = iEvent.id().run();
  nt.ls = iEvent.luminosityBlock();
  nt.fromFullAOD = true;
  nt.BSpot_x = theBeamSpot->x0();
  nt.BSpot_y = theBeamSpot->y0();
  nt.BSpot_z = theBeamSpot->z0();

  if (debug_ > 0)
    std::cout << "new Evt " << nt.run << std::endl;

  reco::TrackBase::Point vertex_point;
  bool goodVtx = false;
  for (const reco::Vertex& vtx : *vertices) {
    if (vtx.isFake() || !vtx.isValid())
      continue;
    nt.pv_x = vtx.x();
    nt.pv_y = vtx.y();
    nt.pv_z = vtx.z();
    goodVtx = true;
    break;
  }
  if (!goodVtx)
    return;  //skipping in absence of good vertex
  vertex_point.SetCoordinates(nt.pv_x, nt.pv_y, nt.pv_z);

  //check if path fired, if so save hlt muons
  if (!HLTaccept(iEvent, nt, HLTPaths_))
    return;
  HLTmuon(iEvent, nt.trg_pt, nt.trg_eta, nt.trg_phi, HLTFilters_, debug_);
  HLTmuon(iEvent, nt.prb_pt, nt.prb_eta, nt.prb_phi, ProbeFilters_, debug_);

  // gen information
  MuonGenAnalyzer genmu;
  std::vector<unsigned> matched_muon_idx;
  std::vector<unsigned> matched_track_idx;
  if (!iEvent.isRealData()) {
    genmu.SetInputs(iEvent, genToken_, momPdgId_);
    genmu.FillNtuple(nt);
    auto reco_match_genmu1 = MatchReco<reco::Muon>(*muons, nt.genmu1_eta, nt.genmu1_phi, genRecoDrMatch_);
    auto reco_match_genmu2 = MatchReco<reco::Muon>(*muons, nt.genmu2_eta, nt.genmu2_phi, genRecoDrMatch_);
    if (reco_match_genmu1.first)
      matched_muon_idx.push_back(reco_match_genmu1.second);

    if (reco_match_genmu2.first)
      matched_muon_idx.push_back(reco_match_genmu2.second);

    reco_match_genmu1 = MatchReco<reco::Track>(*tracks, nt.genmu1_eta, nt.genmu1_phi, genRecoDrMatch_);
    reco_match_genmu2 = MatchReco<reco::Track>(*tracks, nt.genmu2_eta, nt.genmu2_phi, genRecoDrMatch_);
    if (reco_match_genmu1.first)
      matched_track_idx.push_back(reco_match_genmu1.second);
    if (reco_match_genmu2.first)
      matched_track_idx.push_back(reco_match_genmu2.second);
  }

  // match hlt with offline muon
  std::vector<unsigned> trg_idx;
  for (unsigned itrg = 0; itrg < nt.trg_pt.size(); ++itrg) {
    float minDR = 1000;
    unsigned idx = 0;
    for (auto& mu : *muons) {
      if (minDR < deltaR(nt.trg_eta[itrg], nt.trg_phi[itrg], mu.eta(), mu.phi()))
        continue;
      minDR = deltaR(nt.trg_eta[itrg], nt.trg_phi[itrg], mu.eta(), mu.phi());
      idx = &mu - &muons->at(0);
    }
    if (debug_ > 0)
      std::cout << "Trg " << itrg << " min DR " << minDR << std::endl;
    if (minDR < trgDRwindow_)
      trg_idx.push_back(idx);
    if (minDR < trgDRwindow_ && debug_ > 0)
      std::cout << "matched ! " << std::endl;
  }
  nt.nmuons = muons->size();
  nt.ntag = trg_idx.size();

  //select tags
  RecoTrkAndTransientTrkCollection tag_trkttrk;
  std::vector<bool> genmatched_tag;
  for (const reco::Muon& mu : *muons) {
    if (!mu.passed(pow(2, tagQual_)))
      continue;
    if (!tagSelection_(mu))
      continue;
    if (std::find(trg_idx.begin(), trg_idx.end(), &mu - &muons->at(0)) == trg_idx.end())
      continue;
    tag_trkttrk.emplace_back(std::make_pair(mu, reco::TransientTrack(*mu.bestTrack(), &(*bField))));
    if (std::find(matched_muon_idx.begin(), matched_muon_idx.end(), &mu - &muons->at(0)) != matched_muon_idx.end())
      genmatched_tag.push_back(true);
    else
      genmatched_tag.push_back(false);
  }

  if (debug_ > 0)
    std::cout << "tag muons " << tag_trkttrk.size() << std::endl;

  std::pair<std::vector<unsigned>, std::vector<unsigned>> trk_muon_map;
  for (const reco::Muon& mu : *muons) {
    float minDR = 1000;
    unsigned int idx_trk;
    if (debug_ > 1)
      std::cout << "new trk-muon map entry pt " << mu.pt() << " eta " << mu.eta() << " phi " << mu.phi() << std::endl;
    for (const reco::Track& trk : *tracks) {
      if (mu.charge() != trk.charge())
        continue;
      if (fabs(mu.vz() - trk.vz()) > maxdz_trk_mu_ && maxdz_trk_mu_ > 0)
        continue;
      if (fabs(mu.pt() - trk.pt()) / mu.pt() > maxpt_relative_dif_trk_mu_ && maxpt_relative_dif_trk_mu_ > 0)
        continue;
      float DR = deltaR(mu.eta(), mu.phi(), trk.eta(), trk.phi());
      if (debug_ > 1)
        std::cout << "   DR " << DR << "  " << mu.eta() << "  " << mu.phi() << "  " << trk.eta() << "  " << trk.phi()
                  << std::endl;
      if (minDR < DR)
        continue;
      minDR = DR;
      idx_trk = &trk - &tracks->at(0);
    }
    if (minDR > maxdr_trk_mu_)
      continue;
    trk_muon_map.first.push_back(idx_trk);
    trk_muon_map.second.push_back(&mu - &muons->at(0));
  }

  if (debug_ > 0)
    std::cout << "matched trk-mu " << trk_muon_map.first.size() << std::endl;

  //select probes
  for (auto& tag : tag_trkttrk) {
    if (debug_ > 0)
      std::cout << " new tag pt " << tag.first.pt() << " eta " << tag.first.eta() << " phi " << tag.first.phi()
                << std::endl;
    for (const reco::Track& probe : *tracks) {
      if (debug_ > 1)
        std::cout << "    probe pt " << probe.pt() << " eta " << probe.eta() << " phi " << probe.phi() << "  charge "
                  << probe.charge() << std::endl;
      if (HighPurity_ && !probe.quality(Track::highPurity))
        continue;
      if (tag.first.charge() == probe.charge())
        continue;

      //apply cuts on pairs selected will be saved
      if (!probeSelection_(probe))
        continue;
      if (fabs(tag.first.vz() - probe.vz()) > pairDz_ && pairDz_ > 0)
        continue;
      float mass = DimuonMass(tag.first.pt(), tag.first.eta(), tag.first.phi(), probe.pt(), probe.eta(), probe.phi());

      if (mass < pairMassMin_ || mass > pairMassMax_)
        continue;
      std::vector<reco::TransientTrack> trk_pair = {tag.second, reco::TransientTrack(probe, &(*bField))};
      KlFitter vtx(trk_pair);
      if (RequireVtxCreation_ && !vtx.status())
        continue;
      if (minSVtxProb_ > 0 && vtx.prob() < minSVtxProb_)
        continue;
      math::PtEtaPhiMLorentzVector mu1(tag.first.pt(), tag.first.eta(), tag.first.phi(), MU_MASS);
      math::PtEtaPhiMLorentzVector mu2(probe.pt(), probe.eta(), probe.phi(), MU_MASS);

      FillTagBranches<reco::Muon, reco::Track>(tag.first, *tracks, nt);
      nt.tag_isMatchedGen = genmatched_tag[&tag - &tag_trkttrk[0]];

      std::vector<unsigned>::iterator it =
          std::find(trk_muon_map.first.begin(), trk_muon_map.first.end(), &probe - &tracks->at(0));

      if (it == trk_muon_map.first.end()) {
        reco::Muon fakeMuon;
        fakeMuon.setP4(mu2);
        FillProbeBranches<reco::Muon, reco::Track>(fakeMuon, *tracks, nt, false);
        if (debug_ > 0)
          std::cout << "  unsuccessful probe " << std::endl;
      } else {
        unsigned idx = std::distance(trk_muon_map.first.begin(), it);
        if (debug_ > 0)
          std::cout << "  successful probe pt " << muons->at(trk_muon_map.second[idx]).pt() << " eta "
                    << muons->at(trk_muon_map.second[idx]).eta() << " phi " << muons->at(trk_muon_map.second[idx]).phi()
                    << std::endl;
        FillProbeBranches<reco::Muon, reco::Track>(muons->at(trk_muon_map.second[idx]), *tracks, nt, true);
        /*for ( const std::string path: ProbePaths_){
           if ( deltaR(muons->at(trk_muon_map.second[idx]).eta(),muons->at(trk_muon_map.second[idx]).phi(),nt.prb_eta[&path-&ProbePaths_[0]],nt.prb_phi[&path-&ProbePaths_[0]])<trgDRwindow_)
             nt.probe_trg[&path-&ProbePaths_[0]]=true;
           else nt.probe_trg[&path-&ProbePaths_[0]]=false;       
        }*/
      }

      if (std::find(matched_track_idx.begin(), matched_track_idx.end(), &probe - &tracks->at(0)) !=
          matched_track_idx.end())
        nt.probe_isMatchedGen = true;
      else
        nt.probe_isMatchedGen = false;

      FillPairBranches<reco::Muon, reco::Track>(tag.first, probe, nt);
      nt.iprobe++;
      nt.probe_isHighPurity = probe.quality(Track::highPurity);
      vtx.fillNtuple(nt);
      t1->Fill();
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void MuonFullAODAnalyzer::beginJob() {
  t1 = fs->make<TTree>("Events", "Events");
  nt.SetTree(t1);
  nt.CreateBranches(HLTPaths_);
  // nt.CreateExtraTrgBranches(ProbePaths_);
}

// ------------ method called once each job just after ending the event loop  ------------
void MuonFullAODAnalyzer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonFullAODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

///////////////////////

//define this as a plug-in
DEFINE_FWK_MODULE(MuonFullAODAnalyzer);
