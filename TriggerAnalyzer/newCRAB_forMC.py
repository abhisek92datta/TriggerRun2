from WMCore.Configuration import Configuration
config = Configuration() 
config.section_("General") 
config.General.requestName = 'triggerTree_MC_80X_tthbb'
#config.General.requestName = 'triggerTree_MC_80X_tthnonbb'
#config.General.requestName = 'triggerTree_MC_80X_ttjets_semilep'
#config.General.requestName = 'triggerTree_MC_80X_ttjets_dilep'
#config.General.requestName = 'triggerTree_MC_80X_drell_yan'

config.section_("JobType") 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'test/trigger_analyzer_mc_25ns_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['data']

config.section_("Data") 
config.Data.inputDataset = '/ttHTobb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-premix_withHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'
#config.Data.inputDataset = '/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-premix_withHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-premix_withHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-premix_withHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'

config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/' 
config.Data.allowNonValidInputDataset = True 
config.Data.splitting = 'FileBased' 
config.Data.unitsPerJob = 10
##config.Data.totalUnits = 
config.Data.publication = True 
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' 
config.Data.outputDatasetTag = 'triggerTree_MC_80X_tthbb' 
#config.Data.outputDatasetTag = 'triggerTree_MC_80X_tthnonbb' 
#config.Data.outputDatasetTag = 'triggerTree_MC_80X_ttjets_semilep'
#config.Data.outputDatasetTag = 'triggerTree_MC_80X_ttjets_dilep'
#config.Data.outputDatasetTag = 'triggerTree_MC_80X_drell_yan'
config.Data.outLFNDirBase = '/store/user/abdatta/Trigger_Analysis/Ntuples_without_event_sel/' 
config.Data.ignoreLocality = True

config.section_("Site") 
config.Site.storageSite = 'T3_US_FNALLPC'
