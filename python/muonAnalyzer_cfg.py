'''Author: g. karathanasis. georgios.karathanasis@cern.ch
cfg to run tag and probe ntuple for muon POG. It runs both on AOD and miniAOD
Modified by Andre Frankenthal (as2872@cornell.edu) -- September 2020
usage: cmsRun muonAnalizer_cfg.py option1=value1 option2=value2
'''

from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

# defaults
options.maxEvents = -1

options.register('isFullAOD', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Set to True for AOD datatier"
)

options.register('isMC', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Set to True for MC"
)

options.register('era', 'Run2018',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Era (Run2018/Run2017/Run2016)"
)

options.register('isRun2018D', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Is Run2018D (different global tag)"
)

options.register('globalTag', 'NOTSET',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)

options.register('numThreads', 8,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Number of threads (for CRAB vs non-CRAB execution)"
)

options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report frequency"
)

options.parseArguments()

if '2018' in options.era:
    if options.isRun2018D:
        globaltag = '102X_dataRun2_Prompt_v15'
    else:
        globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
elif '2017' in options.era:
    globaltag = '94X_dataRun2_v11' if not options.isMC else '94X_mc2017_realistic_v17'
elif '2016' in options.era:
    globaltag = '80X_dataRun2_2016SeptRepro_v7' if not options.isMC else '80X_mcRun2_asymptotic_2016_TrancheIV_v8'

if options._beenSet['globalTag']:
    globaltag = options.globalTag


if options.outputFile == "":
    options.outputFile = "output"
    if options.isMC:
        options.outputFile += "_mc"
    else:
        options.outputFile += "_data" 
    if options.isFullAOD:
        options.outputFile += "_full"
    else:
        options.outputFile += "_mini"
    options.outputFile += ".root"


process = cms.Process("MuonAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,globaltag, '')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop *_ctppsPixelClusters_*_*'
    )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(options.numThreads)
)

from MuonAnalysis.MuonAnalyzer.muonAnalysis_cff import *
#process.load("MuonAnalysis.MuonAnalyzer.muonAnalysis_cff")
if options.isFullAOD:
    process = muonAnalysis_customizeFullAOD(process)
else:
    process = muonAnalysis_customizeMiniAOD(process)

process.analysis_step = cms.Path(process.muSequence)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outputFile)
)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.fevt = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(),
    fileName = cms.untracked.string("edm_output.root")
)


from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
