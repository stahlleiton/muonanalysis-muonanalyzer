'''author gkarathanasis
option for AOD run'''

import FWCore.ParameterSet.Config as cms

#paths and corresponding l3 filters
Path=["HLT_Mu8_v","HLT_Mu17_v","HLT_Mu19_v","HLT_Mu20_v","HLT_IsoMu20_v","HLT_IsoMu24_v","HLT_Mu50"]  #WARNING lower than 10 path!!!!
Filter=["hltL3fL1sMu5L1f0L2f5L3Filtered8","hltL3fL1sMu15DQlqL1f0L2f10L3Filtered17","hltL3fL1sMu15DQlqL1f0L2f10L3Filtered19","hltL3fL1sMu18L1f0L2f10QL3Filtered20Q","hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"]

if len(Path)>10:
  print "WARNING either put less than 11 paths, or increase the path quota from NtupleContent.h/.cc"
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
           triggerPaths=cms.vstring(Path),
           triggerFilters=cms.vstring(Filter),
           gen = cms.InputTag("genParticles"),
           ProbePaths=cms.vstring(Path),
           ProbeFilters=cms.vstring(Filter),
           trgDRwindow= cms.double(0.05), # dr winwow hlt mu/offline
           tagQuality = cms.uint32(0),
           tagSelection = cms.string("pt()>10"),
           probeHPurity = cms.bool(False),
           probeSelection = cms.string("pt()>5"),
           pairMassMin = cms.double(60.0),
           pairMassMax = cms.double(140.0),
           pairDz = cms.double(-1),
           RequireVtxCreation = cms.bool(False),
           minSVtxProb = cms.double(0.01),
           maxDzProbeTrkMuon = cms.double(0.01), # max Dz(mu1,mu2)
           maxRelPtProbeTrkMuon = cms.double(1.0),# max [pt(mu)-pt(trk)]/pt(trk) for probe/offline
           maxDRProbeTrkMuon =  cms.double(0.03), # max DR for probe/offline
           maxDRProbeTrkDSA =  cms.double(0.4), # max DR for general track and dSA
           momPdgId= cms.uint32(23),
           genRecoDrMatch = cms.double(0.03),
           debug = cms.int32(0)
           
)

fullAODSequence=cms.Sequence(muon)
