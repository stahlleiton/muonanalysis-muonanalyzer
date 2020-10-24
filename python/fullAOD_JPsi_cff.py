'''author gkarathanasis
option for AOD run'''

import FWCore.ParameterSet.Config as cms

#PbPb18 J/psi settings
Path_probe=["HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v"]
Filter_probe=["hltL3f0L3Mu0L2Mu0DR3p5FilteredNHitQ10M1to5", "hltL3f0DR3p5L3FilteredNHitQ10", "hltL2fDoubleMuOpenL2DR3p5PreFiltered0", "hltL1fL1sL1DoubleMuOpenMAXdR3p5L1Filtered0"]

Path_tag=["HLT_HIL3Mu5_NHitQ10_v",]
Filter_tag=["hltL3fL1sL1SingleMu*OpenL1f0L2f0L3Filtered5NHitQ10"]

InAcceptance = '((abs(eta)<1.2 && pt>=3.5) || (1.2<=abs(eta) && abs(eta)<2.1 && pt>=5.47-1.89*abs(eta)) || (2.1<=abs(eta) && abs(eta)<2.4 && pt>=1.5))'
TightId = "passed('CutBasedIdTight')"
HybridSoftId = "(isGlobalMuon && isTrackerMuon && "\
               "innerTrack.hitPattern.trackerLayersWithMeasurement > 5 && "\
               "innerTrack.hitPattern.pixelLayersWithMeasurement > 0 && "\
               "abs(userFloat('innerTrack_dxy')) < 0.3 && abs(userFloat('innerTrack_dz')) < 20.)"
LooseId = "(isTrackerMuon || isGlobalMuon || isPFMuon)"

muon = cms.EDAnalyzer('MuonFullAODAnalyzer',
           isMC=cms.bool(False),
           pileupInfo=cms.InputTag('addPileupInfo'),
           Rho=cms.InputTag('fixedGridRhoFastjetAll'),
           beamSpot=cms.InputTag('offlineBeamSpot'),
           vertices=cms.InputTag("offlinePrimaryVerticesRecovery"),
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
           centrality = cms.InputTag("centralityBin","HFtowers"),
           ProbePaths=cms.vstring(Path_probe),
           ProbeFilters=cms.vstring(Filter_probe),
           probeFlags = cms.PSet(
             isHybridSoft = cms.string(HybridSoftId),
             InAcceptance = cms.string(InAcceptance),
           ),
           trgDRwindow = cms.double(0.1), # dr winwow hlt L3 mu/offline
           trgRelDPtwindow = cms.double(10.0), # rel dpt winwow hlt L3 mu/offline
           trgL2DRwindow = cms.double(0.3), # dr winwow hlt L2 mu/offline
           trgL2RelDPtwindow = cms.double(10.0), # rel dpt winwow hlt L2 mu/offline
           trgL1DRwindow = cms.double(0.3), # dr winwow hlt L1 mu/offline
           trgL1DEtawindow = cms.double(0.2), # deta winwow hlt L1 mu/offline
           tagSelection = cms.string(InAcceptance+" && "+HybridSoftId),
           probeSelection = cms.string(InAcceptance+" && "+LooseId),
           pairSelection = cms.string("2.0 < mass < 4.0"),
           RequireVtxCreation = cms.bool(False),
           minSVtxProb = cms.double(-0.01),
           maxDRProbeTrkDSA =  cms.double(0.1), # max DR for general track and dSA
           momPdgId= cms.uint32(443),
           genRecoDrMatch = cms.double(0.1),
           genRecoRelDPtMatch = cms.double(10.0),
           debug = cms.int32(0)
)

from MuonAnalysis.MuonAnalyzer.OfflinePrimaryVerticesRecovery_cfi import *
from MuonAnalysis.MuonAnalyzer.collisionEventSelection_cff import *
from RecoHI.HiCentralityAlgos.CentralityBin_cfi import *
centralityBin.Centrality = cms.InputTag("hiCentrality")
centralityBin.centralityVariable = cms.string("HFtowers")
centralityBin.nonDefaultGlauberModel = cms.string("")
primaryVertexFilter.src = cms.InputTag("offlinePrimaryVerticesRecovery")

fullAODSequence=cms.Sequence(offlinePrimaryVerticesRecovery + collisionEventSelectionAODv2 + centralityBin + muon)
