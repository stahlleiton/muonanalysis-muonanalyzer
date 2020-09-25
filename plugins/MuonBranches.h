//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// filling functions for aod and miniaod tag/probe

#ifndef MuonAnalysis_MuonAnalyzer_MuonBranches
#define MuonAnalysis_MuonAnalyzer_MuonBranches

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include <type_traits>
#include "NtupleContent.h"
#include "helper.h"

template <typename MUON, typename TRK>
inline void FillTagBranches(const MUON &muon, const std::vector<TRK> &tracks,
                            NtupleContent &nt) {
  nt.tag_pt = muon.pt();
  nt.tag_eta = muon.eta();
  nt.tag_phi = muon.phi();
  nt.tag_isLoose = muon.passed(reco::Muon::CutBasedIdLoose);
  nt.tag_isMedium = muon.passed(reco::Muon::CutBasedIdMedium);
  nt.tag_isTight = muon.passed(reco::Muon::CutBasedIdTight);
  nt.tag_isSoft = muon.passed(reco::Muon::SoftCutBasedId);
  nt.tag_isHighPt = muon.passed(reco::Muon::CutBasedIdTrkHighPt);
  float iso04 =
      (TrackerEnergy04<TRK>(muon.eta(), muon.phi(), tracks) - muon.pt()) /
      muon.pt();
  nt.tag_relIso04 = (iso04 > 0) ? iso04 : 0;
}

template <typename MUON, typename TRK>
inline void FillProbeBranches(const MUON &mu, const std::vector<TRK> &tracks,
                              NtupleContent &nt, bool success) {
  nt.probe_pt = mu.pt();
  nt.probe_eta = mu.eta();
  nt.probe_phi = mu.phi();
  float iso04 =
      (TrackerEnergy04<TRK>(mu.eta(), mu.phi(), tracks) - mu.pt()) / mu.pt();
  nt.probe_relIso04 = (iso04 > 0) ? iso04 : 0;
  // success --> muon obj and track match in dR
  if (success) {
    nt.probe_isLoose = mu.passed(reco::Muon::CutBasedIdLoose);
    nt.probe_isMedium = mu.passed(reco::Muon::CutBasedIdMedium);
    nt.probe_isTight = mu.passed(reco::Muon::CutBasedIdTight);
    nt.probe_isSoft = mu.passed(reco::Muon::SoftCutBasedId);
    nt.probe_isHighPt = mu.passed(reco::Muon::CutBasedIdTrkHighPt);
    nt.probe_isPF = mu.isPFMuon();
    nt.probe_isSA = mu.isStandAloneMuon();
    nt.probe_isTracker = mu.isTrackerMuon();
    nt.probe_isGlobal = mu.isGlobalMuon();
    if (mu.globalTrack().isNonnull())
      nt.probe_trkChi2 = mu.globalTrack()->normalizedChi2();
    else
      nt.probe_trkChi2 = -99;
    if (mu.innerTrack().isNonnull() && mu.innerTrack().isAvailable()) {
      nt.probe_validFraction = mu.innerTrack()->validFraction();
      nt.probe_trackerLayers =
          mu.innerTrack()->hitPattern().trackerLayersWithMeasurement();
      nt.probe_pixelLayers =
          mu.innerTrack()->hitPattern().pixelLayersWithMeasurement();
      nt.probe_dxy = mu.innerTrack()->dxy(
          reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
      nt.probe_dz = mu.innerTrack()->dz(
          reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
    } else {
      nt.probe_validFraction = -99;
      nt.probe_trackerLayers = -99;
      nt.probe_pixelLayers = -99;
      nt.probe_dxy = -99;
      nt.probe_dz = -99;
    }
    if (mu.outerTrack().isNonnull() && mu.outerTrack().isAvailable()) {
      nt.probe_muonStations =
          mu.outerTrack()->hitPattern().muonStationsWithValidHits();
      nt.probe_muonHits = mu.outerTrack()->hitPattern().numberOfValidMuonHits();
      nt.probe_DTHits = mu.outerTrack()->hitPattern().numberOfValidMuonDTHits();
      nt.probe_CSCHits =
          mu.outerTrack()->hitPattern().numberOfValidMuonCSCHits();
    } else {
      nt.probe_muonStations = -99;
      nt.probe_muonHits = -99;
      nt.probe_DTHits = -99;
      nt.probe_CSCHits = -99;
    }
    nt.probe_positionChi2 = mu.combinedQuality().chi2LocalPosition;
    nt.probe_trkKink = mu.combinedQuality().trkKink;
    //     nt.probe_segmentCompatibility=
    nt.probe_isMuMatched = true;
  }
  // no successs (no match)
  else {
    nt.probe_isLoose = false;
    nt.probe_isMedium = false;
    nt.probe_isTight = false;
    nt.probe_isSoft = false;
    nt.probe_isHighPt = false;
    nt.probe_isMuMatched = false;
    nt.probe_isPF = false;
    nt.probe_isSA = false;
    nt.probe_isTracker = false;
    nt.probe_isGlobal = false;
    nt.probe_validFraction = -99;
    nt.probe_trkChi2 = -99;
    nt.probe_positionChi2 = -99;
    nt.probe_trkKink = -99;
    nt.probe_trackerLayers = -99;
    nt.probe_pixelLayers = -99;
    nt.probe_dxy = -99;
    nt.probe_dz = -99;
    nt.probe_muonStations = -99;
    nt.probe_muonHits = -99;
    nt.probe_DTHits = -99;
    nt.probe_CSCHits = -99;
    nt.probe_pterr = -99;
  }
}

template <typename TRK>
inline void FillProbeBranchesdSA(const TRK &trk, NtupleContent &nt,
                                 bool passdSA) {
  nt.probe_isdSA = passdSA;

  nt.probe_dsa_pt = trk.pt();
  nt.probe_dsa_eta = trk.eta();
  nt.probe_dsa_phi = trk.phi();

  if (passdSA) {
    nt.probe_dsa_dxy =
        trk.dxy(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
    nt.probe_dsa_dz = trk.dz(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
    nt.probe_dsa_muonStations = trk.hitPattern().muonStationsWithValidHits();
    nt.probe_dsa_muonHits = trk.hitPattern().numberOfValidMuonHits();
    nt.probe_dsa_DTHits = trk.hitPattern().numberOfValidMuonDTHits();
    nt.probe_dsa_CSCHits = trk.hitPattern().numberOfValidMuonCSCHits();
    nt.probe_dsa_pterr = trk.ptError() / trk.pt();
    nt.probe_dsa_trkChi2 = trk.normalizedChi2();
  } else {
    nt.probe_dsa_dxy = -99;
    nt.probe_dsa_dz = -99;
    nt.probe_dsa_muonStations = -99;
    nt.probe_dsa_muonHits = -99;
    nt.probe_dsa_DTHits = -99;
    nt.probe_dsa_CSCHits = -99;
    nt.probe_dsa_pterr = -99;
    nt.probe_dsa_trkChi2 = -99;
  }
}

template <typename TRK>
inline void FillProbeBranchesdgl(const TRK &trk, NtupleContent &nt,
                                 bool passdgl) {
  nt.probe_isdGlobal = passdgl;
}

template <typename TRK>
inline void FillProbeBranchesCosmic(const TRK &trk, NtupleContent &nt,
                                    bool passcosmic) {
  nt.probe_isCosmic = passcosmic;
}

template <typename MUO, typename TRK>
inline void FillPairBranches(const MUO &muon, const TRK &trk,
                             NtupleContent &nt) {
  math::PtEtaPhiMLorentzVector mu1(muon.pt(), muon.eta(), muon.phi(), MU_MASS);
  math::PtEtaPhiMLorentzVector mu2(trk.pt(), trk.eta(), trk.phi(), MU_MASS);
  nt.pair_pt = (mu1 + mu2).pt();
  nt.pair_mass = (mu1 + mu2).mass();
  nt.pair_eta = (mu1 + mu2).eta();
  nt.pair_phi = (mu1 + mu2).phi();
  nt.pair_dz = muon.vz() - trk.vz();
}

#endif
