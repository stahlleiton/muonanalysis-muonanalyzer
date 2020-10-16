// Package:    MuonAnalyzer for Run 3
//             version 2.0
//
/**\class

 Description: Ntuplizer class for full AOD files
*/
//
// Original Author:
//                george karathanasis
//         Created:  Thu, 20 feb 2020 17:40:23 GMT
//
// Modified:
//                Andre Frankenthal (Sept. 2020)
//
//

// system include files
#include <iostream>
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <vector>
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TRegexp.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "KlFitter.h"
#include "MuonBranches.h"
#include "MuonGenAnalyzer.h"
#include "NtupleContent.h"
#include "helper.h"

using namespace std;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
class MuonFullAODAnalyzer : public edm::one::EDAnalyzer<> {
 public:
  typedef std::vector<std::pair<reco::Muon, reco::TransientTrack>>
      RecoTrkAndTransientTrkCollection;
  explicit MuonFullAODAnalyzer(const edm::ParameterSet&);
  ~MuonFullAODAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  void beginJob() override;
  bool HLTaccept(const edm::Event&, NtupleContent&, std::vector<std::string>&);
  typedef std::map<std::string, std::vector<std::tuple<float, float, float, std::set<size_t>, std::set<size_t> > > > TRGInfo;
  void HLTmuon(const edm::Event&, TRGInfo&, std::vector<std::string>&, const int&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupSummaryToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  edm::EDGetToken muonsToken_;
  edm::EDGetToken tracksToken_;
  edm::EDGetToken dSAToken_;
  edm::EDGetToken dglToken_;
  edm::EDGetToken cosmicToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genToken_;
  std::vector<std::string> HLTPaths_;
  std::vector<std::string> HLTFilters_;
  std::vector<std::string> ProbePaths_;
  std::vector<std::string> ProbeFilters_;

  const bool isMC_;
  const double trgDRwindow_;
  const unsigned int tagQual_;
  const StringCutObjectSelector<reco::Muon>
      tagSelection_;  // kinematic cuts for tag
  const bool HighPurity_;
  const StringCutObjectSelector<reco::Track>
      probeSelection_;  // kinematic cuts for probe
  const double pairMassMin_;
  const double pairMassMax_;
  const double pairDz_;
  const bool RequireVtxCreation_;  // if true skip pairs that do not create
                                   // that do not have a vertex
  const double
      minSVtxProb_;  // min probability of a vertex to be kept. If <0 inactive
  const double maxdz_trk_mu_;
  const double maxpt_relative_dif_trk_mu_;
  const double maxdr_trk_mu_;
  const double maxdr_trk_dsa_;
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
      rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("Rho"))),
      pileupSummaryToken_(consumes<std::vector<PileupSummaryInfo>>(
          iConfig.getParameter<edm::InputTag>("pileupInfo"))),
      beamSpotToken_(consumes<reco::BeamSpot>(
          iConfig.getParameter<edm::InputTag>("beamSpot"))),
      vtxToken_(consumes<std::vector<reco::Vertex>>(
          iConfig.getParameter<edm::InputTag>("vertices"))),
      muonsToken_(consumes<std::vector<reco::Muon>>(
          iConfig.getParameter<edm::InputTag>("muons"))),
      tracksToken_(consumes<std::vector<reco::Track>>(
          iConfig.getParameter<edm::InputTag>("tracks"))),
      dSAToken_(consumes<std::vector<reco::Track>>(
          iConfig.getParameter<edm::InputTag>("dSAmuons"))),
      dglToken_(consumes<std::vector<reco::Track>>(
          iConfig.getParameter<edm::InputTag>("dGlmuons"))),
      cosmicToken_(consumes<std::vector<reco::Track>>(
          iConfig.getParameter<edm::InputTag>("staCosmic"))),
      trgresultsToken_(consumes<edm::TriggerResults>(
          iConfig.getParameter<edm::InputTag>("triggerResults"))),
      trigobjectsToken_(consumes<trigger::TriggerEvent>(
          iConfig.getParameter<edm::InputTag>("triggerObjects"))),
      genToken_(consumes<edm::View<reco::GenParticle>>(
          iConfig.getParameter<edm::InputTag>("gen"))),
      HLTPaths_(iConfig.getParameter<std::vector<std::string>>("triggerPaths")),
      HLTFilters_(
          iConfig.getParameter<std::vector<std::string>>("triggerFilters")),
      ProbePaths_(iConfig.getParameter<std::vector<std::string>>("ProbePaths")),
      ProbeFilters_(
          iConfig.getParameter<std::vector<std::string>>("ProbeFilters")),
      isMC_(iConfig.getParameter<bool>("isMC")),
      trgDRwindow_(iConfig.getParameter<double>("trgDRwindow")),
      tagQual_(iConfig.getParameter<unsigned>("tagQuality")),
      tagSelection_(iConfig.getParameter<std::string>("tagSelection")),
      HighPurity_(iConfig.getParameter<bool>("probeHPurity")),
      probeSelection_(iConfig.getParameter<std::string>("probeSelection")),
      pairMassMin_(iConfig.getParameter<double>("pairMassMin")),
      pairMassMax_(iConfig.getParameter<double>("pairMassMax")),
      pairDz_(iConfig.getParameter<double>("pairDz")),
      RequireVtxCreation_(iConfig.getParameter<bool>("RequireVtxCreation")),
      minSVtxProb_(iConfig.getParameter<double>("minSVtxProb")),
      maxdz_trk_mu_(iConfig.getParameter<double>("maxDzProbeTrkMuon")),
      maxpt_relative_dif_trk_mu_(
          iConfig.getParameter<double>("maxRelPtProbeTrkMuon")),
      maxdr_trk_mu_(iConfig.getParameter<double>("maxDRProbeTrkMuon")),
      maxdr_trk_dsa_(iConfig.getParameter<double>("maxDRProbeTrkDSA")),
      momPdgId_(iConfig.getParameter<unsigned>("momPdgId")),
      genRecoDrMatch_(iConfig.getParameter<double>("genRecoDrMatch")),
      debug_(iConfig.getParameter<int>("debug"))

{
  //  edm::ParameterSet
  //  runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");
}

MuonFullAODAnalyzer::~MuonFullAODAnalyzer() {
  // cout << "total " << trg_counter << " fires " << fire_counter << " l3"
  // << l3_counter << endl; do anything here that needs to be done at desctruction
  // time
}

bool MuonFullAODAnalyzer::HLTaccept(const edm::Event& iEvent, NtupleContent& nt,
                                    std::vector<std::string>& HLTPaths) {
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
      if (!trigResults->accept(itrg)) continue;
      if (!TrigPath.Contains(TRegexp(TString(path)))) continue;
      EvtFire = true;
      TrgFire = true;
      break;
    }
    nt.trigger[ipath] = TrgFire;
    ipath++;
  }
  return EvtFire;
}

void MuonFullAODAnalyzer::HLTmuon(const edm::Event& iEvent,
                                  TRGInfo& trgInfo,
                                  std::vector<std::string>& HLTFilters,
                                  const int& debug_) {
  edm::Handle<trigger::TriggerEvent> triggerObjects;
  iEvent.getByToken(trigobjectsToken_, triggerObjects);
  trigger::TriggerObjectCollection allTriggerObjects = triggerObjects->getObjects();
  // get filter information per trigger object
  std::map<size_t, std::set<size_t> > tagFilters, probeFilters;
  for (size_t filterIndex=0; filterIndex < triggerObjects->sizeFilters(); filterIndex++) {
    TString TrigFilter = triggerObjects->filterLabel(filterIndex);
    const trigger::Keys& keys = (*triggerObjects).filterKeys(filterIndex);
    for (size_t ipath=0; ipath<HLTFilters.size(); ipath++) {
      if (!TrigFilter.Contains(TRegexp(TString(HLTFilters[ipath])))) continue;
      for (size_t j = 0; j < keys.size(); j++) tagFilters[keys[j]].insert(ipath);
    }
    for (size_t ipath=0; ipath<ProbeFilters_.size(); ipath++) {
      if (!TrigFilter.Contains(TRegexp(TString(ProbeFilters_[ipath])))) continue;
      for (size_t j = 0; j < keys.size(); j++) probeFilters[keys[j]].insert(ipath);
    }
  }
  if (tagFilters.empty() && probeFilters.empty()) return;
  // add trigger objects per collection
  const auto& cK = triggerObjects->collectionKeys();
  for (size_t i=1; i<cK.size(); i++) {
    const auto& coll = triggerObjects->collectionTag(i).encode();
    for (size_t j=cK[i-1]; j<cK[i]; j++) {
      trigger::TriggerObject foundObject = (allTriggerObjects)[j];
      if (fabs(foundObject.id()) != 13) continue;
      trgInfo[coll].emplace_back(std::make_tuple(foundObject.pt(), foundObject.eta(), foundObject.phi(), tagFilters[j], probeFilters[j]));
      if (debug_ > 0) std::cout << "Trg muon " << foundObject.pt() << std::endl;
    }
  }
}

//
void MuonFullAODAnalyzer::analyze(const edm::Event& iEvent,
                                  const edm::EventSetup& iSetup) {
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  // Get data
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_, theBeamSpot);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  // Skip evts if there are no vertices
  if (vertices->empty()) return;
  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);
  edm::Handle<std::vector<reco::Track>> tracks;
  iEvent.getByToken(tracksToken_, tracks);
  edm::Handle<std::vector<reco::Track>> dSAmuons;
  iEvent.getByToken(dSAToken_, dSAmuons);
  edm::Handle<std::vector<reco::Track>> dGlmuons;
  iEvent.getByToken(dglToken_, dGlmuons);
  edm::Handle<std::vector<reco::Track>> staCosmic;
  iEvent.getByToken(cosmicToken_, staCosmic);
  edm::ESHandle<MagneticField> bField;
  iSetup.get<IdealMagneticFieldRecord>().get(bField);

  // Information about run
  nt.ClearBranches();
  nt.run = iEvent.id().run();
  nt.ls = iEvent.luminosityBlock();
  nt.fromFullAOD = true;
  nt.BSpot_x = theBeamSpot->x0();
  nt.BSpot_y = theBeamSpot->y0();
  nt.BSpot_z = theBeamSpot->z0();
  nt.nvertices = vertices->size();

  // Pileup information
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);
  nt.Rho = *rhoHandle;

  float trueNumInteractions = -1;
  int puNumInteractions = -1;
  if (isMC_) {
    edm::Handle<std::vector<PileupSummaryInfo>> PupInfo;
    iEvent.getByToken(pileupSummaryToken_, PupInfo);

    for (auto PVI : *PupInfo) {
      int BX = PVI.getBunchCrossing();
      if (BX == 0) {
        trueNumInteractions = PVI.getTrueNumInteractions();
        puNumInteractions = PVI.getPU_NumInteractions();
        continue;
      }
    }
  }

  nt.trueNumInteractions = trueNumInteractions;
  nt.puNumInteractions = puNumInteractions;

  if (debug_ > 0) std::cout << "New Evt " << nt.run << std::endl;

  reco::TrackBase::Point vertex_point;
  bool goodVtx = false;
  for (const reco::Vertex& vtx : *vertices) {
    if (vtx.isFake() || !vtx.isValid()) continue;
    nt.pv_x = vtx.x();
    nt.pv_y = vtx.y();
    nt.pv_z = vtx.z();
    goodVtx = true;
    break;
  }
  if (!goodVtx) return;  // skipping in absence of good vertex
  vertex_point.SetCoordinates(nt.pv_x, nt.pv_y, nt.pv_z);

  // check if path fired, if so save hlt muons
  if (!HLTaccept(iEvent, nt, HLTPaths_)) return;
  TRGInfo trgInfo;
  auto filters = HLTFilters_;
  filters.insert(filters.begin(), ProbeFilters_.begin(), ProbeFilters_.end());
  HLTmuon(iEvent, trgInfo, filters, debug_);

  // gen information
  MuonGenAnalyzer genmu;
  std::vector<unsigned> matched_muon_idx;
  std::vector<unsigned> matched_track_idx;
  if (!iEvent.isRealData()) {
    genmu.SetInputs(iEvent, genToken_, momPdgId_);
    genmu.FillNtuple(nt);
    auto reco_match_genmu1 = MatchReco<reco::Muon>(
        *muons, nt.genmu1_eta, nt.genmu1_phi, genRecoDrMatch_);
    auto reco_match_genmu2 = MatchReco<reco::Muon>(
        *muons, nt.genmu2_eta, nt.genmu2_phi, genRecoDrMatch_);
    if (reco_match_genmu1.first)
      matched_muon_idx.push_back(reco_match_genmu1.second);

    if (reco_match_genmu2.first)
      matched_muon_idx.push_back(reco_match_genmu2.second);

    reco_match_genmu1 = MatchReco<reco::Track>(*tracks, nt.genmu1_eta,
                                               nt.genmu1_phi, genRecoDrMatch_);
    reco_match_genmu2 = MatchReco<reco::Track>(*tracks, nt.genmu2_eta,
                                               nt.genmu2_phi, genRecoDrMatch_);
    if (reco_match_genmu1.first)
      matched_track_idx.push_back(reco_match_genmu1.second);
    if (reco_match_genmu2.first)
      matched_track_idx.push_back(reco_match_genmu2.second);
  }

  // match hlt with offline muon for tags
  std::vector<std::vector<size_t> > trg_tag_idx(HLTFilters_.size());
  for (size_t imu=0; imu<muons->size(); imu++) {
    const auto& mu = muons->at(imu);
    for (const auto& coll : trgInfo) {
      float minDR = 1000;
      size_t idx = 0;
      for (size_t iobj=0; iobj<coll.second.size(); iobj++) {
        const auto& obj = coll.second[iobj];
        const auto& dR = deltaR(std::get<1>(obj), std::get<2>(obj), mu.eta(), mu.phi());
        if (minDR < dR) continue;
        minDR = dR;
        idx = iobj;
      }
      if (minDR >= trgDRwindow_) continue;
      const auto& obj = coll.second[idx];
      for (size_t ipath=0; ipath<HLTFilters_.size(); ipath++) {
        if (std::get<3>(obj).find(ipath)==std::get<3>(obj).end()) continue;
        trg_tag_idx[ipath].push_back(imu);
      }
    }
  }
  nt.nmuons = muons->size();
  nt.ntag = trg_tag_idx.size();

  // match hlt with offline muon for probes
  std::vector<std::vector<size_t> > trg_prb_idx(ProbeFilters_.size());
  for (size_t imu=0; imu<tracks->size(); imu++) {
    const auto& mu = tracks->at(imu);
    for (const auto& coll : trgInfo) {
      float minDR = 1000;
      size_t idx = 0;
      for (size_t iobj=0; iobj<coll.second.size(); iobj++) {
        const auto& obj = coll.second[iobj];
        const auto& dR = deltaR(std::get<1>(obj), std::get<2>(obj), mu.eta(), mu.phi());
        if (minDR < dR) continue;
        minDR = dR;
        idx = iobj;
      }
      if (minDR >= trgDRwindow_) continue;
      const auto& obj = coll.second[idx];
      for (size_t ipath=0; ipath<ProbeFilters_.size(); ipath++) {
        if (std::get<4>(obj).find(ipath)==std::get<4>(obj).end()) continue;
        trg_prb_idx[ipath].push_back(imu);
      }
    }
  }

  // select tags
  RecoTrkAndTransientTrkCollection tag_trkttrk;
  std::vector<bool> genmatched_tag;
  for (const reco::Muon& mu : *muons) {
    if (!mu.passed(pow(2, tagQual_))) continue;
    if (!tagSelection_(mu)) continue;
    bool isTrgMatch = false;
    for (const auto& trg : trg_tag_idx) {
      if (std::find(trg.begin(), trg.end(), &mu - &muons->at(0)) != trg.end()) {
        isTrgMatch = true;
        break;
      }
    }
    if (!isTrgMatch) continue;
    tag_trkttrk.emplace_back(
        std::make_pair(mu, reco::TransientTrack(*mu.bestTrack(), &(*bField))));
    if (std::find(matched_muon_idx.begin(), matched_muon_idx.end(),
                  &mu - &muons->at(0)) != matched_muon_idx.end())
      genmatched_tag.push_back(true);
    else
      genmatched_tag.push_back(false);
  }

  if (debug_ > 0) std::cout << "Tag muons " << tag_trkttrk.size() << std::endl;

  // mapping with muon object
  std::pair<std::vector<unsigned>, std::vector<unsigned>> trk_muon_map;
  for (const reco::Muon& mu : *muons) {
    float minDR = 1000;
    unsigned int idx_trk;
    if (debug_ > 1)
      std::cout << "New trk-muon map entry pt " << mu.pt() << " eta "
                << mu.eta() << " phi " << mu.phi() << std::endl;
    for (const reco::Track& trk : *tracks) {
      if (mu.charge() != trk.charge()) continue;
      if (mu.innerTrack().isNull()) continue;
      const auto& muTrk = *mu.innerTrack();
      if (fabs(muTrk.vz() - trk.vz()) > maxdz_trk_mu_ && maxdz_trk_mu_ > 0)
        continue;
      if (fabs(muTrk.pt() - trk.pt()) / muTrk.pt() > maxpt_relative_dif_trk_mu_ &&
          maxpt_relative_dif_trk_mu_ > 0)
        continue;
      float DR = deltaR(muTrk.eta(), muTrk.phi(), trk.eta(), trk.phi());
      if (debug_ > 1)
        std::cout << "   DR " << DR << "  " << mu.eta() << "  " << mu.phi()
                  << "  " << trk.eta() << "  " << trk.phi() << std::endl;
      if (minDR < DR) continue;
      minDR = DR;
      idx_trk = &trk - &tracks->at(0);
    }
    if (minDR > maxdr_trk_mu_) continue;
    trk_muon_map.first.push_back(idx_trk);
    trk_muon_map.second.push_back(&mu - &muons->at(0));
  }

  if (debug_ > 0)
    std::cout << "Matched trk-mu " << trk_muon_map.first.size() << std::endl;

  // mapping with displaced standalone muon
  std::pair<std::vector<unsigned>, std::vector<unsigned>> trk_dsA_map;
  for (const reco::Track& dsA : *dSAmuons) {
    float minDeltaR = 1000;
    unsigned int idx_dSAtrk;
    for (const reco::Track& trks : *tracks) {
      float deltaDR = deltaR(dsA.eta(), dsA.phi(), trks.eta(), trks.phi());
      if (minDeltaR < deltaDR) continue;
      minDeltaR = deltaDR;
      idx_dSAtrk = &trks - &tracks->at(0);
    }
    if (minDeltaR > maxdr_trk_dsa_) continue;
    trk_dsA_map.first.push_back(idx_dSAtrk);
    trk_dsA_map.second.push_back(&dsA - &dSAmuons->at(0));
  }
  if (debug_ > 0)
    std::cout << "Matched trk-dSA " << trk_dsA_map.first.size() << std::endl;

  // mapping with global displaced muon
  std::pair<std::vector<unsigned>, std::vector<unsigned>> trk_dglobal_map;
  for (const reco::Track& dgl : *dGlmuons) {
    float minDltR = 1000;
    unsigned int idx_dgltrk;
    for (const reco::Track& tk : *tracks) {
      float dltDR = deltaR(dgl.eta(), dgl.phi(), tk.eta(), tk.phi());
      if (minDltR < dltDR) continue;
      minDltR = dltDR;
      idx_dgltrk = &tk - &tracks->at(0);
    }
    if (minDltR > maxdr_trk_dsa_) continue;
    trk_dglobal_map.first.push_back(idx_dgltrk);
    trk_dglobal_map.second.push_back(&dgl - &dGlmuons->at(0));
  }
  if (debug_ > 0)
    std::cout << "Matched trk-dgl " << trk_dglobal_map.first.size()
              << std::endl;

  // mapping with standalone cosmic muon
  std::pair<std::vector<unsigned>, std::vector<unsigned>> trk_cosmic_map;
  for (const reco::Track& cosmicMu : *staCosmic) {
    float minDtR = 1000;
    unsigned int idx_cosmictrk;
    for (const reco::Track& tks : *tracks) {
      float dtDR = deltaR(cosmicMu.eta(), cosmicMu.phi(), tks.eta(), tks.phi());
      if (minDtR < dtDR) continue;
      minDtR = dtDR;
      idx_cosmictrk = &tks - &tracks->at(0);
    }
    if (minDtR > maxdr_trk_dsa_) continue;
    trk_cosmic_map.first.push_back(idx_cosmictrk);
    trk_cosmic_map.second.push_back(&cosmicMu - &staCosmic->at(0));
  }
  if (debug_ > 0)
    std::cout << "Matched trk-cosmic " << trk_cosmic_map.first.size()
              << std::endl;

  // select probes
  for (auto& tag : tag_trkttrk) {
    if (debug_ > 0)
      std::cout << "New tag pt " << tag.first.pt() << " eta "
                << tag.first.eta() << " phi " << tag.first.phi() << std::endl;
    for (const reco::Track& probe : *tracks) {
      if (debug_ > 1)
        std::cout << "    Probe pt " << probe.pt() << " eta " << probe.eta()
                  << " phi " << probe.phi() << "  charge " << probe.charge()
                  << std::endl;
      if (HighPurity_ && !probe.quality(Track::highPurity)) continue;
      if (tag.first.charge() == probe.charge()) continue;

      // apply cuts on pairs selected will be saved
      if (!probeSelection_(probe)) continue;
      if (fabs(tag.first.vz() - probe.vz()) > pairDz_ && pairDz_ > 0) continue;
      float mass = DimuonMass(tag.first.pt(), tag.first.eta(), tag.first.phi(),
                              probe.pt(), probe.eta(), probe.phi());

      if (mass < pairMassMin_ || mass > pairMassMax_) continue;
      std::vector<reco::TransientTrack> trk_pair = {
          tag.second, reco::TransientTrack(probe, &(*bField))};
      KlFitter vtx(trk_pair);
      if (RequireVtxCreation_ && !vtx.status()) continue;
      if (minSVtxProb_ > 0 && vtx.prob() < minSVtxProb_) continue;
      math::PtEtaPhiMLorentzVector mu1(tag.first.pt(), tag.first.eta(),
                                       tag.first.phi(), MU_MASS);
      math::PtEtaPhiMLorentzVector mu2(probe.pt(), probe.eta(), probe.phi(),
                                       MU_MASS);

      nt.tag_isMatchedGen = genmatched_tag[&tag - &tag_trkttrk[0]];

      FillTagBranches<reco::Muon, reco::Track>(tag.first, *tracks, nt);

      if (std::find(matched_track_idx.begin(), matched_track_idx.end(),
                    &probe - &tracks->at(0)) != matched_track_idx.end())
        nt.probe_isMatchedGen = true;
      else
        nt.probe_isMatchedGen = false;

      std::vector<unsigned>::iterator it =
          std::find(trk_muon_map.first.begin(), trk_muon_map.first.end(),
                    &probe - &tracks->at(0));
      std::vector<unsigned>::iterator itdsa =
          std::find(trk_dsA_map.first.begin(), trk_dsA_map.first.end(),
                    &probe - &tracks->at(0));
      std::vector<unsigned>::iterator itdgl =
          std::find(trk_dglobal_map.first.begin(), trk_dglobal_map.first.end(),
                    &probe - &tracks->at(0));
      std::vector<unsigned>::iterator itcosmic =
          std::find(trk_cosmic_map.first.begin(), trk_cosmic_map.first.end(),
                    &probe - &tracks->at(0));

      // probe trigger matching
      for (size_t ipath=0; ipath<trg_prb_idx.size(); ipath++) {
        const auto& trg = trg_prb_idx[ipath];
        nt.probe_trg[ipath] = (std::find(trg.begin(), trg.end(), &probe - &tracks->at(0)) != trg.end());
      }

      if (it == trk_muon_map.first.end()) {
        reco::Muon fakeMuon;
        fakeMuon.setP4(mu2);
        FillProbeBranches<reco::Muon, reco::Track>(fakeMuon, *tracks, nt,
                                                   false);
        if (debug_ > 0) std::cout << "  Unsuccessful probe " << std::endl;
      } else {
        unsigned idx = std::distance(trk_muon_map.first.begin(), it);
        if (debug_ > 0)
          std::cout << "  Successful probe pt "
                    << muons->at(trk_muon_map.second[idx]).pt() << " eta "
                    << muons->at(trk_muon_map.second[idx]).eta() << " phi "
                    << muons->at(trk_muon_map.second[idx]).phi() << std::endl;
        FillProbeBranches<reco::Muon, reco::Track>(
            muons->at(trk_muon_map.second[idx]), *tracks, nt, true);
        /*for ( const std::string path: ProbePaths_){
           if (
        deltaR(muons->at(trk_muon_map.second[idx]).eta(),muons->at(
                trk_muon_map.second[idx]).phi(),nt.prb_eta[&path-&ProbePaths_[0]],
                nt.prb_phi[&path-&ProbePaths_[0]])<trgDRwindow_)
             nt.probe_trg[&path-&ProbePaths_[0]]=true;
           else nt.probe_trg[&path-&ProbePaths_[0]]=false;
        }*/
      }

      if (itdsa == trk_dsA_map.first.end()) {
        FillProbeBranchesdSA<reco::Track>(probe, nt, false);
      } else {
        unsigned itdsax = std::distance(trk_dsA_map.first.begin(), itdsa);
        if (debug_ > 0)
          std::cout << "Successful probe dsA "
                    << dSAmuons->at(trk_dsA_map.second[itdsax]).pt() << " eta "
                    << dSAmuons->at(trk_dsA_map.second[itdsax]).eta() << " phi "
                    << dSAmuons->at(trk_dsA_map.second[itdsax]).phi()
                    << std::endl;
        FillProbeBranchesdSA<reco::Track>(
            dSAmuons->at(trk_dsA_map.second[itdsax]), nt, true);
      }

      if (itdgl == trk_dglobal_map.first.end()) {
        FillProbeBranchesdgl<reco::Track>(probe, nt, false);
      } else {
        unsigned itdglx = std::distance(trk_dglobal_map.first.begin(), itdgl);
        if (debug_ > 0)
          std::cout << "Successful probe displaced global "
                    << dGlmuons->at(trk_dglobal_map.second[itdglx]).pt()
                    << " eta "
                    << dGlmuons->at(trk_dglobal_map.second[itdglx]).eta()
                    << " phi "
                    << dGlmuons->at(trk_dglobal_map.second[itdglx]).phi()
                    << std::endl;
        FillProbeBranchesdgl<reco::Track>(
            dGlmuons->at(trk_dglobal_map.second[itdglx]), nt, true);
      }

      if (itcosmic == trk_cosmic_map.first.end()) {
        FillProbeBranchesCosmic<reco::Track>(probe, nt, false);
      } else {
        unsigned itcosmicx =
            std::distance(trk_cosmic_map.first.begin(), itcosmic);
        if (debug_ > 0)
          std::cout << "Successful probe cosmic "
                    << staCosmic->at(trk_cosmic_map.second[itcosmicx]).pt()
                    << " eta "
                    << staCosmic->at(trk_cosmic_map.second[itcosmicx]).eta()
                    << " phi "
                    << staCosmic->at(trk_cosmic_map.second[itcosmicx]).phi()
                    << std::endl;
        FillProbeBranchesCosmic<reco::Track>(
            staCosmic->at(trk_cosmic_map.second[itcosmicx]), nt, true);
      }

      FillPairBranches<reco::Muon, reco::Track>(tag.first, probe, nt);
      nt.iprobe++;
      nt.probe_isHighPurity = probe.quality(Track::highPurity);
      vtx.fillNtuple(nt);
      t1->Fill();
    }
  }
}

// ------------ method called once each job just before starting event loop
// ------------
void MuonFullAODAnalyzer::beginJob() {
  t1 = fs->make<TTree>("Events", "Events");
  nt.SetTree(t1);
  nt.CreateBranches(HLTPaths_);
  nt.CreateExtraTrgBranches(ProbePaths_);
}

// ------------ method called once each job just after ending the event loop
// ------------
void MuonFullAODAnalyzer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void MuonFullAODAnalyzer::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

///////////////////////

// define this as a plug-in
DEFINE_FWK_MODULE(MuonFullAODAnalyzer);
