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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "helper.h"
#include "NtupleContent.h"
#include <type_traits>

template<typename MUON,typename TRK>
inline void FillTagBranches( const MUON & muon, const std::vector<TRK> &tracks, NtupleContent & nt){
  nt.tag_pt=muon.pt();
  nt.tag_eta=muon.eta();
  nt.tag_phi=muon.phi();
  nt.tag_isLoose=muon.passed(reco::Muon::CutBasedIdLoose);
  nt.tag_isMedium=muon.passed(reco::Muon::CutBasedIdMedium);
  nt.tag_isTight=muon.passed(reco::Muon::CutBasedIdTight);
  nt.tag_isSoft=muon.passed(reco::Muon::SoftCutBasedId);
  nt.tag_isHighPt=muon.passed(reco::Muon::CutBasedIdTrkHighPt);
  float iso04= (TrackerEnergy04<TRK> (muon.eta(),muon.phi(), tracks) -muon.pt() )/ muon.pt();
  nt.tag_relIso04=(iso04>0) ? iso04 : 0;

}


template<typename MUON,typename TRK>
inline void FillProbeBranches( const MUON & mu, const std::vector<TRK> &tracks, NtupleContent & nt, bool success){
  nt.probe_pt=mu.pt();
  nt.probe_eta=mu.eta();
  nt.probe_phi=mu.phi();
  float iso04= (TrackerEnergy04<TRK> (mu.eta(),mu.phi(), tracks) -mu.pt() )/ mu.pt();
  nt.probe_relIso04=(iso04>0) ? iso04 : 0;
  if (success){
     nt.probe_isLoose=mu.passed(reco::Muon::CutBasedIdLoose);
     nt.probe_isMedium=mu.passed(reco::Muon::CutBasedIdMedium);
     nt.probe_isTight=mu.passed(reco::Muon::CutBasedIdTight);
     nt.probe_isSoft=mu.passed(reco::Muon::SoftCutBasedId);
     nt.probe_isHighPt=mu.passed(reco::Muon::CutBasedIdTrkHighPt);
     nt.probe_isPF=mu.isPFMuon();
     nt.probe_isGlobal=mu.isGlobalMuon();
     if ( mu.globalTrack().isNonnull() )
        nt.probe_trkChi2=mu.globalTrack()->normalizedChi2();
     else
        nt.probe_trkChi2=-99;
     if ( mu.innerTrack().isNonnull() &&  mu.innerTrack().isAvailable() ){
        nt.probe_validFraction=mu.innerTrack()->validFraction();
        nt.probe_trackerLayers=mu.innerTrack()->hitPattern().trackerLayersWithMeasurement();
        nt.probe_pixelLayers=mu.innerTrack()->hitPattern().pixelLayersWithMeasurement();
        nt.probe_dxy=mu.innerTrack()->dxy(reco::TrackBase::Point(nt.pv_x,nt.pv_y,nt.pv_z));
        nt.probe_dz=mu.innerTrack()->dz(reco::TrackBase::Point(nt.pv_x,nt.pv_y,nt.pv_z));
    } else{
        nt.probe_validFraction=-99;   nt.probe_trackerLayers=-99;   
        nt.probe_pixelLayers=-99;     nt.probe_dxy=-99;  
        nt.probe_dz=-99;
    }
     nt.probe_positionChi2=mu.combinedQuality().chi2LocalPosition;
     nt.probe_trkKink=mu.combinedQuality().trkKink;
//     nt.probe_segmentCompatibility= 
     nt.probe_isMuMatched=true;
  }
   else{
     
     nt.probe_isLoose=false;        nt.probe_isMedium=false;
     nt.probe_isTight=false;        nt.probe_isSoft=false;
     nt.probe_isHighPt=false;       nt.probe_isMuMatched=false;
     nt.probe_isPF=false;           nt.probe_isGlobal=false;
     nt.probe_validFraction=-99;    nt.probe_trkChi2=-99;
     nt.probe_positionChi2=-99;     nt.probe_trkKink=-99;
     nt. probe_trackerLayers=-99;   nt. probe_pixelLayers=-99;
     nt.probe_dxy=-99;              nt.probe_dz=-99;
  }
}

template<typename MUO, typename TRK>
inline void FillPairBranches( const MUO & muon, const TRK & trk, NtupleContent & nt){
   math::PtEtaPhiMLorentzVector mu1(muon.pt(), muon.eta(), muon.phi(), MU_MASS);
   math::PtEtaPhiMLorentzVector mu2(trk.pt(), trk.eta(), trk.phi(), MU_MASS); 
   nt.pair_pt=(mu1+mu2).pt();   
   nt.pair_mass=(mu1+mu2).mass();        
   nt.pair_eta=(mu1+mu2).eta();
   nt.pair_phi=(mu1+mu2).phi();   
   nt.pair_dz=muon.vz()-trk.vz();
}


#endif
