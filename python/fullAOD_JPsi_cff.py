'''author gkarathanasis
option for AOD run'''

import FWCore.ParameterSet.Config as cms

#paths and corresponding l3 filters
Path=["HLT_Mu7p5_Track7_Jpsi", "HLT_Mu7p5_Track3p5_Jpsi", "HLT_Mu7p5_Track2_Jpsi"]  #WARNING lower than 10 path!!!!
Filter=["hltL3fLMu7p5TrackL3Filtered7p5", "hltL3fLMu7p5TrackL3Filtered7p5", "hltL3fLMu7p5TrackL3Filtered7p5"]

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
           #tagSelection = cms.string("pt()>10"),
           tagSelection = cms.string("pt()>0"),
           probeHPurity = cms.bool(False),
           #probeSelection = cms.string("pt()>5"),
           probeSelection = cms.string("pt()>0"),
           #pairMassMin = cms.double(60.0),
           #pairMassMax = cms.double(140.0),
           pairMassMin = cms.double(2.8),
           pairMassMax = cms.double(3.4),
           #pairDz = cms.double(-1),
           pairDz = cms.double(10.1),
           #RequireVtxCreation = cms.bool(False),
           RequireVtxCreation = cms.bool(True),
           minSVtxProb = cms.double(0.01),
           maxDzProbeTrkMuon = cms.double(0.01), # max Dz(mu1,mu2)
           maxRelPtProbeTrkMuon = cms.double(1.0),# max [pt(mu)-pt(trk)]/pt(trk) for probe/offline
           maxDRProbeTrkMuon =  cms.double(0.03), # max DR for probe/offline
           maxDRProbeTrkDSA =  cms.double(0.2), # max DR for general track and dSA
           momPdgId= cms.uint32(443),
           genRecoDrMatch = cms.double(0.03),
           debug = cms.int32(0)
           
)

fullAODSequence=cms.Sequence(muon)
