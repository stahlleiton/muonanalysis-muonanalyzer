'''Author: g. karathanasis. georgios.karathanasis@cern.ch
cfg to run tag and probe ntuple for muon POG. It runs both on AOD and miniAOD
Modified by Andre Frankenthal (as2872@cornell.edu) -- September 2020
usage: cmsRun run_muonAnalyzer_cfg.py option1=value1 option2=value2
'''

from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

# defaults
options.maxEvents = -1

options.register('isFullAOD', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Set to False for MiniAOD datatier"
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

options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report frequency"
)

options.register('numThreads', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Number of CMSSW threads" 
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


if options.isFullAOD:
    if options.isMC and len(options.inputFiles)==0:
        options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/80002/FF9AF238-78B6-CF48-BC7C-05025D85A45C.root')
    elif not options.isMC and len(options.inputFiles)==0:
        options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/60004/FB123080-071C-F64D-BAFD-F2F292F7FC64.root')
else:
    if options.isMC and len(options.inputFiles)==0:
        options.inputFiles.append('/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/FAACE9E0-1D0E-204E-9960-078F095EA34C.root')
    elif not options.isMC and len(options.inputFiles)==0:
        options.inputFiles.append('/store/data/Run2018D/SingleMuon/MINIAOD/PromptReco-v2/000/325/159/00000/BE9BB28C-8AEC-1B4B-A7BD-AD1C9A0D67A8.root')


if options.outputFile=="":
    options.outputFile="output"
    if options.isMC:
        options.outputFile+="_mc"
    else:
        options.outputFile+="_data" 
    if options.isFullAOD:
        options.outputFile+="_full"
    else:
        options.outputFile+="_mini"
    options.outputFile+=".root"


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
    secondaryFileNames=cms.untracked.vstring(),
    inputCommands=cms.untracked.vstring(
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
   
