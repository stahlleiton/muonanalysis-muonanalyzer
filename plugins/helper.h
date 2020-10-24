//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// helper functions

#ifndef MuonAnalysis_MuonAnalyzer_plugins_helper
#define MuonAnalysis_MuonAnalyzer_plugins_helper
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"

const float MU_MASS = 0.10565837;
template <class T, class C>
inline double DimuonMass(const T& mu1, const C& mu2) {
  math::PtEtaPhiMLorentzVector mu1P4(mu1.pt(), mu1.eta(), mu1.phi(), MU_MASS);
  math::PtEtaPhiMLorentzVector mu2P4(mu2.pt(), mu2.eta(), mu2.phi(), MU_MASS);
  return (mu1P4 + mu2P4).mass();
}

template <typename TRK>
inline float TrackerEnergy04(const float eta_muon, const float phi_muon,
                             const std::vector<TRK> tracks) {
  float energy = 0;
  for (const auto& trk : tracks) {
    if (deltaR(eta_muon, phi_muon, trk.eta(), trk.phi()) > 0.4) continue;
    energy += trk.pt();
  }
  return energy;
}

template <class C>
bool MatchByDeltaR(size_t& index,
                   const math::XYZVector& s, const std::vector<C>& mV,
                   const double& maxDR=1000., const double& maxRelDPt=-1., const double& maxDEta=-1.) {
  std::pair<double, size_t> match({maxDR, 0});
  for (const auto& m : mV) {
    if (maxRelDPt > 0 && std::abs(s.rho() - m.pt())/m.pt() >= maxRelDPt) continue;
    if (maxDEta > 0 && std::abs(s.eta() - m.eta()) >= maxDEta) continue;
    const auto& dR = deltaR(s.eta(), s.phi(), m.eta(), m.phi());
    if (dR < match.first)
      match = std::make_pair(dR, &m - &mV.at(0));
  }
  if (match.first >= maxDR) return false;
  index = match.second;
  return true;
}

template <class T, class C>
bool MatchByDeltaR(size_t& index,
                   const T& s, const std::vector<C>& mV,
                   const double& maxDR=1000., const double& maxRelDPt=-1., const double& maxDEta=-1.) {
  math::XYZVector mom(s.px(), s.py(), s.pz());
  return MatchByDeltaR(index, mom, mV, maxDR, maxRelDPt, maxDEta);
}

template <class C>
bool MatchByTrackRef(size_t& index,
                     const reco::TrackRef& s, const std::vector<C>& mV) {
  for (const auto& m : mV) {
    if (m.track().isNull() || m.track()!=s) continue;
    index = &m - &mV.at(0);
    return true;
  }
  return false;
}

#endif
