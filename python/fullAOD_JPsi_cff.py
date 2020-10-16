'''author gkarathanasis
option for AOD run'''

import FWCore.ParameterSet.Config as cms

#paths and corresponding l3 filters
#Path=["HLT_Mu8_v","HLT_Mu17_v","HLT_Mu19_v","HLT_Mu20_v","HLT_IsoMu20_v","HLT_IsoMu24_v","HLT_Mu50"]  #WARNING lower than 10 path!!!!
#Filter=["hltL3fL1sMu5L1f0L2f5L3Filtered8","hltL3fL1sMu15DQlqL1f0L2f10L3Filtered17","hltL3fL1sMu15DQlqL1f0L2f10L3Filtered19","hltL3fL1sMu18L1f0L2f10QL3Filtered20Q","hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"]

#PbPb18 J/psi
Path=["HLT_HIL3Mu5_NHitQ10_v","HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v","HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v","HLT_HIL3Mu3_v","HLT_HIL3Mu5_v", "HLT_HIL3Mu7_v", "HLT_HIL3Mu12_v"]  #WARNING lower than 10 path!!!!
Filter=["hltL3fL1sL1SingleMuOpenL1f0L2f0L3Filtered5NHitQ10","hltL3f0L3Mu0L2Mu0DR3p5FilteredNHitQ10M1to5","hltL3f0L3Mu2p5NHitQ10L2Mu2FilteredM7toinf","hltL3fL1sSingleMu3L1f0L2f0L3Filtered3", "hltL3fL1sSingleMu3OR5L1f0L2f0L3Filtered5","hltL3fL1sSingleMu3OR5L1f0L2f0L3Filtered7","hltL3fL1sSingleMu7L1f0L2f0L3Filtered12"]

#PbPb18 J/psi
Path_tag=["HLT_HIL3Mu5_NHitQ10_v*",]  #WARNING lower than 10 path!!!!
Filter_tag=["hltL3fL1sL1SingleMuOpenL1f0L2f0L3Filtered5NHitQ10"]

if len(Path)>10:
  print ("WARNING either put less than 11 paths, or increase the path quota from NtupleContent.h/.cc")
  exit()

muon = cms.EDAnalyzer('MuonFullAODAnalyzer',
           isMC=cms.bool(False),
           pileupInfo=cms.InputTag('addPileupInfo'),
           Rho=cms.InputTag('fixedGridRhoFastjetAll'),
           beamSpot=cms.InputTag('offlineBeamSpot'),
           vertices=cms.InputTag("offlinePrimaryVertices"),
           muons=cms.InputTag("muons"),
           tracks=cms.InputTag("generalTracks"),
           dSAmuons=cms.InputTag("displacedStandAloneMuons"),
           dGlmuons=cms.InputTag("displacedGlobalMuons"),
           staCosmic=cms.InputTag("cosmicMuons"),
           triggerResults=cms.InputTag("TriggerResults::HLT"),
           triggerObjects=cms.InputTag('hltTriggerSummaryAOD'),
           triggerPaths=cms.vstring(Path_tag),
           triggerFilters=cms.vstring(Filter_tag),
           gen = cms.InputTag("genParticles"),
           ProbePaths=cms.vstring(Path),
           ProbeFilters=cms.vstring(Filter),
           trgDRwindow= cms.double(0.3), # dr winwow hlt mu/offline
           tagQuality = cms.uint32(3),
           tagSelection = cms.string("((abs(eta)<1.2 && pt>=3.5) || (1.2<=abs(eta) && abs(eta)<2.1 && pt>=5.47-1.89*abs(eta)) || (2.1<=abs(eta) && abs(eta)<2.4 && pt>=1.5)) && abs(eta)<2.4"),
           #tagSelection = cms.string("pt()>0"),
           probeHPurity = cms.bool(False),
           probeSelection = cms.string("((abs(eta)<1.2 && pt>=3.5) || (1.2<=abs(eta) && abs(eta)<2.1 && pt>=5.47-1.89*abs(eta)) || (2.1<=abs(eta) && abs(eta)<2.4 && pt>=1.5)) && abs(eta)<2.4"),
           #probeSelection = cms.string("pt()>0"),
           pairMassMin = cms.double(2),
           pairMassMax = cms.double(4),
           pairDz = cms.double(10.1),
           RequireVtxCreation = cms.bool(False),
           minSVtxProb = cms.double(-0.01),
           maxDzProbeTrkMuon = cms.double(10.01), # max Dz(mu1,mu2)
           maxRelPtProbeTrkMuon = cms.double(20.0),# max [pt(mu)-pt(trk)]/pt(trk) for probe/offline
           maxDRProbeTrkMuon =  cms.double(0.03), # max DR for probe/offline
           maxDRProbeTrkDSA =  cms.double(0.2), # max DR for general track and dSA
           momPdgId= cms.uint32(443),
           genRecoDrMatch = cms.double(0.03),
           debug = cms.int32(0)

)

fullAODSequence=cms.Sequence(muon)
