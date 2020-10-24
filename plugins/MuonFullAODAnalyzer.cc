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
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
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
#include "TrackingTools/IPTools/interface/IPTools.h"

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
  typedef std::map<std::string, pat::TriggerObjectStandAloneCollection> TRGInfo;
  explicit MuonFullAODAnalyzer(const edm::ParameterSet&);
  ~MuonFullAODAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  void beginJob() override;
  bool HLTaccept(const edm::Event&, NtupleContent&);
  void HLTmuon(const edm::Event&, TRGInfo&, const int&);
  pat::Muon MakePatMuon(const reco::Muon&, const bool&, const reco::Vertex&, const edm::ESHandle<MagneticField>&);
  math::XYZVector MomentumAt2ndMuonStation(const reco::Muon& muon, const edm::ESHandle<GlobalTrackingGeometry>& geometry);
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
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::EDGetTokenT<int> centToken_;
  std::map<std::string, std::vector<std::string>> HLTPaths_;
  std::map<std::string, std::vector<std::string>> HLTFilters_;

  const edm::ParameterSet probeFlags_;
  const double trgDRwindow_, trgRelDPtwindow_;
  const double trgL2DRwindow_, trgL2RelDPtwindow_;
  const double trgL1DRwindow_, trgL1DEtawindow_;
  const StringCutObjectSelector<pat::Muon> tagSelection_;  // kinematic cuts for tag
  const StringCutObjectSelector<pat::Muon> probeSelection_;  // kinematic cuts for probe
  const StringCutObjectSelector<reco::CompositeCandidate> pairSelection_; // kinematic cuts for pair
  const bool RequireVtxCreation_;  // if true skip pairs that do not create
                                   // that do not have a vertex
  const double minSVtxProb_;  // min probability of a vertex to be kept. If <0 inactive
  const double maxdr_trk_dsa_;
  const unsigned momPdgId_;
  const double genRecoDrMatch_, genRecoRelDPtMatch_;
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
      genToken_(consumes<reco::GenParticleCollection>(
          iConfig.getParameter<edm::InputTag>("gen"))),
      centToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("centrality"))),
      probeFlags_(iConfig.getParameter<edm::ParameterSet>("probeFlags")),
      trgDRwindow_(iConfig.getParameter<double>("trgDRwindow")),
      trgRelDPtwindow_(iConfig.getParameter<double>("trgRelDPtwindow")),
      trgL2DRwindow_(iConfig.getParameter<double>("trgL2DRwindow")),
      trgL2RelDPtwindow_(iConfig.getParameter<double>("trgL2RelDPtwindow")),
      trgL1DRwindow_(iConfig.getParameter<double>("trgL1DRwindow")),
      trgL1DEtawindow_(iConfig.getParameter<double>("trgL1DEtawindow")),
      tagSelection_(iConfig.getParameter<std::string>("tagSelection")),
      probeSelection_(iConfig.getParameter<std::string>("probeSelection")),
      pairSelection_(iConfig.getParameter<std::string>("pairSelection")),
      RequireVtxCreation_(iConfig.getParameter<bool>("RequireVtxCreation")),
      minSVtxProb_(iConfig.getParameter<double>("minSVtxProb")),
      maxdr_trk_dsa_(iConfig.getParameter<double>("maxDRProbeTrkDSA")),
      momPdgId_(iConfig.getParameter<unsigned>("momPdgId")),
      genRecoDrMatch_(iConfig.getParameter<double>("genRecoDrMatch")),
      genRecoRelDPtMatch_(iConfig.getParameter<double>("genRecoRelDPtMatch")),
      debug_(iConfig.getParameter<int>("debug"))

{
  //  edm::ParameterSet
  //  runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");
  HLTPaths_["tag"] = iConfig.getParameter<std::vector<std::string>>("triggerPaths");
  HLTFilters_["tag"] = iConfig.getParameter<std::vector<std::string>>("triggerFilters");
  HLTPaths_["probe"] = iConfig.getParameter<std::vector<std::string>>("ProbePaths");
  HLTFilters_["probe"] = iConfig.getParameter<std::vector<std::string>>("ProbeFilters");
}

MuonFullAODAnalyzer::~MuonFullAODAnalyzer() {
  // cout << "total " << trg_counter << " fires " << fire_counter << " l3"
  // << l3_counter << endl; do anything here that needs to be done at desctruction
  // time
}

bool MuonFullAODAnalyzer::HLTaccept(const edm::Event& iEvent, NtupleContent& nt) {
  edm::Handle<edm::TriggerResults> trgR;
  iEvent.getByToken(trgresultsToken_, trgR);
  const auto& trgN = iEvent.triggerNames(*trgR);
  // fill probe trigger results
  for (const auto& path : HLTPaths_["probe"]) {
    bool trgFire = false;
    for (size_t i = 0; i < trgR->size(); i++)
      if (trgR->accept(i) && TString(trgN.triggerName(i)).Contains(TRegexp(path))) {
        trgFire = true;
        break;
      }
    nt.trigger.at(path) = trgFire;
  }
  // return true if at least one tag trigger fired
  for (const auto& path : HLTPaths_["tag"])
    for (size_t i = 0; i < trgR->size(); i++)
      if (trgR->accept(i) && TString(trgN.triggerName(i)).Contains(TRegexp(path)))
        return true;
  return false;
}

void MuonFullAODAnalyzer::HLTmuon(const edm::Event& iEvent,
                                  TRGInfo& trgInfo,
                                  const int& debug_) {
  edm::Handle<trigger::TriggerEvent> triggerObjects;
  iEvent.getByToken(trigobjectsToken_, triggerObjects);
  if (!triggerObjects.isValid()) return;
  // define set of filters
  std::set<std::string> HLTFilters(HLTFilters_["tag"].begin(), HLTFilters_["tag"].end());
  std::copy(HLTFilters_["probe"].begin(), HLTFilters_["probe"].end(), std::inserter(HLTFilters, HLTFilters.end()));
  // get filter information per trigger object
  std::map<size_t, std::set<std::string>> objFilters;
  for (size_t filterIndex=0; filterIndex < triggerObjects->sizeFilters(); filterIndex++) {
    const auto& name = triggerObjects->filterLabel(filterIndex);
    const auto& keys = triggerObjects->filterKeys(filterIndex);
    for (const auto& filter : HLTFilters)
      if (TString(name).Contains(TRegexp(filter)))
        for (const auto& j : keys) objFilters[j].insert(filter);
  }
  if (objFilters.empty()) return;
  // add trigger objects per collection
  const auto& colKeys = triggerObjects->collectionKeys();
  const auto& objects = triggerObjects->getObjects();
  for (const auto& c : {"L3Muon", "L2Muon", "Stage2Digis:Muon"})
    for (size_t i=0; i<colKeys.size(); i++) {
      const auto& col = triggerObjects->collectionTag(i).encode();
      if (col.find(c)==std::string::npos) continue;
      for (size_t j=(i<1 ? 0 : colKeys[i-1]); j<colKeys[i]; j++) {
        trgInfo[col].emplace_back(objects[j]);
        trgInfo[col].back().setCollection(col);
        for (const auto& name : objFilters[j])
          trgInfo[col].back().addFilterLabel(name);
        if (debug_ > 0) std::cout << c << " pt: " << trgInfo[col].back().pt() << std::endl;
      }
    }
}

pat::Muon MuonFullAODAnalyzer::MakePatMuon(const reco::Muon& src, const bool& isMuon,
                                           const reco::Vertex& vertex, const edm::ESHandle<MagneticField>& bField) {
  pat::Muon muon(src);
  muon.addUserInt("isMuon", isMuon);
  // add track information
  const auto& gTrack = muon.globalTrack();
  if (gTrack.isNonnull() && gTrack.isAvailable())
    muon.setNormChi2(gTrack->chi2() / gTrack->ndof());
  auto bTrack = muon.muonBestTrack();
  if (bTrack.isNull() || !bTrack.isAvailable())
    bTrack = muon.innerTrack();
  if (bTrack.isNonnull() && bTrack.isAvailable())
    muon.setNumberOfValidHits(bTrack->numberOfValidHits());
  // add vertex information
  const auto& iTrack = muon.innerTrack();
  if (iTrack.isNonnull() && iTrack.isAvailable()) {
    muon.addUserFloat("innerTrack_dxy", iTrack->dxy(vertex.position()));
    muon.addUserFloat("innerTrack_dz", iTrack->dz(vertex.position()));
  }
  if (bTrack.isNonnull() && bTrack.isAvailable()) {
    muon.addUserFloat("muonBestTrack_dxy", bTrack->dxy(vertex.position()));
    muon.addUserFloat("muonBestTrack_dz", bTrack->dz(vertex.position()));
    const auto& mom = GlobalVector(bTrack->px(), bTrack->py(), bTrack->pz());
    const auto& tt = reco::TransientTrack(*bTrack, &(*bField));
    const auto& ip = IPTools::signedTransverseImpactParameter(tt, mom, vertex);
    muon.setDB(ip.second.value(), ip.second.error(), pat::Muon::PV2D);
  }
  return muon;
}

math::XYZVector MuonFullAODAnalyzer::MomentumAt2ndMuonStation(const reco::Muon& muon,
                                                              const edm::ESHandle<GlobalTrackingGeometry>& geometry) {
  GlobalPoint pos;
  for (const auto& m : muon.matches()) {
    int station = 999;
    if (m.id.subdetId() == MuonSubdetId::DT)
      station = DTChamberId(m.id.rawId()).station();
    else if (m.id.subdetId()==MuonSubdetId::CSC)
      station = CSCDetId(m.id.rawId()).station();
    if (station > 3) continue;
    const auto& geo = geometry->idToDet(m.id.rawId());
    if (geo) pos = geo->toGlobal(LocalPoint(m.x, m.y, 0));
    if (station == 2) break;
  }
  if (pos==GlobalPoint()) return {};
  return math::PtEtaPhiMLorentzVector(muon.pt(), pos.eta(), pos.phi()-(M_PI/144.), 0).Vect(); // L1 phi offset: 1.25*pi/180
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
  edm::ESHandle<GlobalTrackingGeometry> geometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(geometry);

  // Information about run
  nt.ClearBranches();
  nt.run = iEvent.id().run();
  nt.ls = iEvent.luminosityBlock();
  nt.event = iEvent.id().event();
  nt.fromFullAOD = true;
  nt.BSpot_x = theBeamSpot->x0();
  nt.BSpot_y = theBeamSpot->y0();
  nt.BSpot_z = theBeamSpot->z0();
  nt.nvertices = vertices->size();

  edm::Handle<int> cent;
  iEvent.getByToken(centToken_, cent);
  if (cent.isValid()) nt.cent = *cent;

  // Pileup information
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);
  nt.Rho = *rhoHandle;

  if (!iEvent.isRealData()) {
    edm::Handle<std::vector<PileupSummaryInfo>> PupInfo;
    iEvent.getByToken(pileupSummaryToken_, PupInfo);
    for (const auto& PVI : *PupInfo) {
      if (PVI.getBunchCrossing() != 0) continue;
      nt.trueNumInteractions = PVI.getTrueNumInteractions();
      nt.puNumInteractions = PVI.getPU_NumInteractions();
      break;
    }
  }

  if (debug_ > 0) std::cout << "New Evt " << nt.run << std::endl;

  reco::Vertex vertex;
  for (const reco::Vertex& vtx : *vertices) {
    if (vtx.isFake() || !vtx.isValid()) continue;
    nt.pv_x = vtx.x();
    nt.pv_y = vtx.y();
    nt.pv_z = vtx.z();
    vertex = vtx;
    break;
  }
  if (!vertex.isValid()) return;  // skipping in absence of good vertex

  // fill HLT trigger decisions
  if (!HLTaccept(iEvent, nt)) return;

  // select probes
  std::map<std::string, pat::MuonCollection> source;
  for (size_t idx=0; idx<tracks->size(); idx++) {
    const auto& trk = reco::TrackRef(tracks, idx);
    pat::Muon probe;
    size_t imu;
    if (MatchByTrackRef(imu, trk, *muons))
      probe = MakePatMuon(muons->at(imu), true, vertex, bField);
    else {
      math::XYZTLorentzVector p4(trk->px(), trk->py(), trk->pz(), std::sqrt(trk->p()*trk->p() + MU_MASS*MU_MASS));
      reco::Muon mu(trk->charge(), p4, trk->vertex());
      mu.setInnerTrack(trk);
      mu.setBestTrack(reco::Muon::InnerTrack);
      probe = MakePatMuon(mu, false, vertex, bField);
    }
    if (probeSelection_(probe)) source["probe"].push_back(probe);
  }
  if (source["probe"].empty()) return;

  // select tags
  for (const auto& mu : *muons) {
    const auto& tag = MakePatMuon(mu, true, vertex, bField);
    if (tagSelection_(tag)) source["tag"].push_back(tag);
  }
  if (source["tag"].empty()) return;
  nt.nmuons = muons->size();

  // perform trigger-reco matching
  TRGInfo trgInfo;
  HLTmuon(iEvent, trgInfo, debug_);
  for (const auto& c : trgInfo)
    for (auto& s : source)
      for (auto& mu : s.second) {
        size_t iobj;
        if (c.first.find("L3Muon")!=std::string::npos) {
          if (MatchByDeltaR(iobj, mu, c.second, trgDRwindow_, trgRelDPtwindow_))
            mu.addTriggerObjectMatch(c.second[iobj]);
        }
        else if (c.first.find("L2Muon")!=std::string::npos) {
          if (MatchByDeltaR(iobj, mu, c.second, trgL2DRwindow_, trgL2RelDPtwindow_))
            mu.addTriggerObjectMatch(c.second[iobj]);
        }
        else if (c.first.find("Stage2Digis:Muon")!=std::string::npos) {
          const auto& mom = MomentumAt2ndMuonStation(mu, geometry);
          if (mom.r()>0. && MatchByDeltaR(iobj, mom, c.second, trgL1DRwindow_, -1., trgL1DEtawindow_))
            mu.addTriggerObjectMatch(c.second[iobj]);
        }
      }

  // remove tags that did not trigger
  auto r = std::remove_if(source["tag"].begin(), source["tag"].end(), [&](auto& tag) -> bool {
    for (const auto& f : HLTFilters_["tag"])
      if (!tag.triggerObjectMatchesByFilter(f).empty())
        return false;
    return true;
  });
  source["tag"].erase(r, source["tag"].end());
  if (source["tag"].empty()) return;
  nt.ntag = source["tag"].size();

  // extract gen information
  MuonGenAnalyzer genmu;
  if (!iEvent.isRealData()) {
    genmu.SetInputs(iEvent, genToken_, momPdgId_);
    genmu.FillNtuple(nt);
    // perform gen-track matching
    size_t idx;
    std::map<reco::TrackRef, reco::GenParticleRef> genRecoMap;
    for (const auto& gmu : genmu.gmuons)
      if (MatchByDeltaR(idx, *gmu, *tracks, genRecoDrMatch_, genRecoRelDPtMatch_))
        genRecoMap.emplace(reco::TrackRef(tracks, idx), gmu);
    // add gen particle
    for (auto& s : source)
      for (auto& mu : s.second)
        if (genRecoMap[mu.track()].isNonnull())
          mu.addGenParticleRef(genRecoMap[mu.track()]);
  }

  // perform general - muon track matching
  std::map<std::string, std::map<reco::TrackRef, reco::TrackRef>> trk_map;
  std::map<std::string, edm::Handle<std::vector<reco::Track>>> muTrkMap = {
    {"dSTA", dSAmuons}, // displaced standalone muon
    {"dGlb", dGlmuons}, // global displaced muon
    {"cSTA", staCosmic} // standalone cosmic muon
  };
  for (const auto& c : muTrkMap) {
    for (const auto& mu : *c.second) {
      size_t idx;
      if (MatchByDeltaR(idx, mu, *tracks, maxdr_trk_dsa_))
        trk_map[c.first].emplace(reco::TrackRef(tracks, idx), reco::TrackRef(c.second, &mu - &c.second->at(0)));
    }
    if (debug_ > 0)
      std::cout << "Matched trk-" << c.first << " " << trk_map[c.first].size() << std::endl;
  }

  // select tag-probe pairs
  for (const auto& tag : source["tag"]) {
    if (debug_ > 0)
      std::cout << "New tag pt " << tag.pt() << " eta "
                << tag.eta() << " phi " << tag.phi() << std::endl;
    for (const auto& probe : source["probe"]) {
      if (debug_ > 1)
        std::cout << "    Probe pt " << probe.pt() << " eta " << probe.eta()
                  << " phi " << probe.phi() << "  charge " << probe.charge()
                  << std::endl;
      if (tag.charge() == probe.charge()) continue;

      const auto& tagP4 = math::PtEtaPhiMLorentzVector(tag.pt(), tag.eta(), tag.phi(), MU_MASS);
      const auto& prbP4 = math::PtEtaPhiMLorentzVector(probe.pt(), probe.eta(), probe.phi(), MU_MASS);
      reco::CompositeCandidate pair(0, tagP4+prbP4);
      pair.addDaughter(tag, "tag");
      pair.addDaughter(probe, "probe");

      // apply cuts on pairs selected will be saved
      if (!pairSelection_(pair)) continue;
      KlFitter vtx({reco::TransientTrack(*tag.bestTrack(), &(*bField)),
                    reco::TransientTrack(*probe.bestTrack(), &(*bField))});
      if (RequireVtxCreation_ && !vtx.status()) continue;
      if (minSVtxProb_ > 0 && vtx.prob() < minSVtxProb_) continue;

      // fill tag information
      nt.tag_isMatchedGen = tag.genParticleRef().isNonnull();;
      FillTagBranches(tag, *tracks, nt);

      // fill probe information
      nt.probe_isMatchedGen = probe.genParticleRef().isNonnull();
      for (const auto& f : HLTFilters_["probe"])
        if (!probe.triggerObjectMatchesByFilter(f).empty())
          nt.probe_trg.at(f) = true;
      FillProbeBranches(probe, *tracks, nt, probe.userInt("isMuon"));
      if (debug_ > 0)
        std::cout << "  " << (probe.userInt("isMuon")?"S":"Uns") << "uccessful"
                  << " probe pt " << probe.pt() << " eta " << probe.eta() << " phi "
                  << probe.phi() << std::endl;

      // fill probe flags
      for (const auto& name : probeFlags_.getParameterNames()) {
        const auto& flag = probeFlags_.getParameter<std::string>(name);
        nt.probe_flag.at(name) = (*std::make_unique<StringCutObjectSelector<pat::Muon>>(flag))(probe);
      }

      // fill muon track information
      for (auto& c : trk_map) {
        const auto& trk = c.second[probe.track()];
        if (c.first=="dSTA")
          FillProbeBranchesdSA(trk.isNonnull() ? *trk : *probe.track(), nt, trk.isNonnull());
        else if (c.first=="dGlb")
          FillProbeBranchesdgl(trk.isNonnull() ? *trk : *probe.track(), nt, trk.isNonnull());
        else if (c.first=="cSTA")
          FillProbeBranchesCosmic(trk.isNonnull() ? *trk : *probe.track(), nt, trk.isNonnull());
        if (debug_ > 0 && trk.isNonnull())
          std::cout << "Successful probe " << c.first <<  " pt "
                    << trk->pt() << " eta " << trk->eta() << " phi " << trk->phi()
                    << std::endl;
      }

      // fill tag-probe pair information
      FillPairBranches(tag, probe, nt);
      nt.iprobe++;
      nt.probe_isHighPurity = probe.track()->quality(Track::highPurity);
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
  nt.CreateBranches(HLTPaths_["probe"]);
  nt.CreateExtraTrgBranches(HLTFilters_["probe"]);
  nt.CreateProbeFlagBranches(probeFlags_.getParameterNames());
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
