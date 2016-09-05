from WMCore.Configuration import Configuration
config = Configuration() 
config.section_("General") 
config.General.requestName = 'triggerTree_SingleElectron_Run2016B'

config.section_("JobType") 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'test/trigger_analyzer_data_80x_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data") 
config.Data.inputDataset = '/SingleElectron/Run2016B-PromptReco-v2/MINIAOD'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/' 
config.Data.allowNonValidInputDataset = True 
config.Data.splitting = 'LumiBased' 
config.Data.lumiMask = 'data/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON_prescale.txt'
config.Data.unitsPerJob = 20
##config.Data.totalUnits = 
config.Data.publication = True 
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' 
config.Data.outputDatasetTag = 'triggerTree_SingleElectron_Run2016B' 
config.Data.outLFNDirBase = '/store/user/abdatta/Trigger_Analysis/Ntuples_without_event_sel/' 
config.Data.ignoreLocality = True

config.section_("Site") 
config.Site.storageSite = 'T3_US_FNALLPC'
