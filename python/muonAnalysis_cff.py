'''author: g karathanasis
load the AOD and miniAOD selections
'''
import FWCore.ParameterSet.Config as cms


def  muonAnalysis_customizeFullAOD(process):
   process.load("MuonAnalysis.MuonAnalyzer.fullAOD_cff")
   process.muSequence = cms.Sequence(process.fullAODSequence)
   return process

def  muonAnalysis_customizeMiniAOD(process):
   process.load("MuonAnalysis.MuonAnalyzer.miniAOD_cff")
   process.muSequence = cms.Sequence( process.miniAODSequence)
   return process
