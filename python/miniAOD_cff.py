'''author: g. karathanasis
parameters for miniAOD'''

import FWCore.ParameterSet.Config as cms


Path=["HLT_Mu8_v","HLT_Mu17_v","HLT_Mu19_v","HLT_Mu20_v","HLT_IsoMu20_v","HLT_IsoMu24_v","HLT_Mu50"]  #paths for tag muon


muon = cms.EDAnalyzer('MuonMiniAODAnalyzer',
           beamSpot=cms.InputTag('offlineBeamSpot'),
           vertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
           muons=cms.InputTag("slimmedMuons"),
           triggerResults=cms.InputTag("TriggerResults::HLT"),
           PFCands=cms.InputTag("packedPFCandidates"),
           lostTracks=cms.InputTag("lostTracks"),
           gen = cms.InputTag("prunedGenParticles"),
           HLTPaths=cms.vstring(Path),
           ProbePaths=cms.vstring(Path),
           tagQuality = cms.uint32(0), # quality of tag muon following muonSelector convention
           tagSelection = cms.string("pt()>0"), # string to pass cuts on tag
           ProbeHPyrity = cms.bool(True), # skips non High purity probes
           probeSelection = cms.string("pt()>0"), #string for probe
           pairMassMin = cms.double(2.9), # min mass of mu pair
           pairMassMax = cms.double(3.3), # max mss of mu pair
           pairDz = cms.double(10.1), #max Dz of mu1,mu2
           RequireVtxCreation = cms.bool(False), # if true kills pairs w/o vtx
           minSVtxProb = cms.double(-0.01), # min prob of mu pair
           maxDzProbeTrkMuon = cms.double(0.01), # max Dz(mu1,mu2)
           maxRelPtProbeTrkMuon = cms.double(1.0),# max [pt(mu)-pt(trk)]/pt(trk) for probe/offline
           maxDRProbeTrkMuon =  cms.double(0.03), # max DR for probe/offline
           momPdgId = cms.uint32(443),
           genRecoDrMatch= cms.double(0.03)
           
)

miniAODSequence=cms.Sequence(muon)
