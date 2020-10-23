from CRABClient.UserUtilities import config, getUsernameFromCRIC
config = config()

config.section_('General')
config.General.requestName = 'TnP_HISingleMuon_ReReco_PbPb5TeV_2018_JPsi_Data_TestForPP_PPCodeRuns327000-327100_v3'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_muonAnalyzer_cfg.py'
config.JobType.pyCfgParams = ['isFullAOD=True', 'isMC=False','globalTag=103X_dataRun2_Prompt_v3','reportEvery=100000', 'maxEvents=-1']
config.JobType.maxMemoryMB = 2000
config.JobType.maxJobRuntimeMin = 500
#config.JobType.outputFiles = ['TnPtree_PbPb_with_pp_AOD.root']
config.JobType.allowUndistributedCMSSW = True #allow use of SL7

config.section_('Data')
config.Data.inputDataset ='/HISingleMuon/HIRun2018A-04Apr2019-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 10
config.Data.splitting = 'LumiBased'
#config.Data.totalUnits = 10
#config.Data.runRange = '326381-327564'
config.Data.runRange = '327000-327100'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON_HF_and_MuonPhys.txt'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/TagAndProbe/PbPbTest/TnP_PbPb/%s' % (getUsernameFromCRIC(), config.General.requestName)
config.Data.publication = False


config.section_('Site')
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN']
config.Site.storageSite = 'T2_CH_CERN'
