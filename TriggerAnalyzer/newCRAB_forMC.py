from WMCore.Configuration import Configuration
config = Configuration() 
config.section_("General") 
config.General.requestName = 'triggerTree_MC_80X_tth'
#config.General.requestName = 'triggerTree_MC_80X_ttjets'

config.section_("JobType") 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'test/trigger_analyzer_mc_25ns_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data") 
config.Data.inputDataset = '/ttHTobb_M125_13TeV_powheg_pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'
#config.Data.inputDataset = ' /TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/MINIAODSIM'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/' 
config.Data.allowNonValidInputDataset = True 
config.Data.splitting = 'FileBased' 
#config.Data.lumiMask = 'data/Cert_271036-274421_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.unitsPerJob = 10
##config.Data.totalUnits = 
config.Data.publication = True 
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' 
config.Data.outputDatasetTag = 'triggerTree_MC_80X_tth' 
#config.Data.outputDatasetTag = 'triggerTree_MC_80X_ttjets'
#config.Data.outLFNDirBase = '/store/user/sflowers/treeMaker/Summer16_June13_Trigger/'
config.Data.outLFNDirBase = '/store/user/abdatta/Trigger_Analysis/Ntuples_without_event_sel/' 
config.Data.ignoreLocality = True

config.section_("Site") 
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.storageSite = 'T3_US_FNALLPC'
