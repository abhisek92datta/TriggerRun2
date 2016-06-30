from WMCore.Configuration import Configuration
config = Configuration() 
config.section_("General") 
config.General.requestName = 'triggerTree_SingleElectron_Run2016B_PromptReco-Try2'

config.section_("JobType") 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'test/trigger_analyzer_data_76x_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data") 
config.Data.inputDataset = '/SingleElectron/Run2016B-PromptReco-v2/MINIAOD'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/' 
config.Data.allowNonValidInputDataset = True 
config.Data.splitting = 'LumiBased' 
config.Data.lumiMask = 'data/Cert_271036-274421_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.unitsPerJob = 20
##config.Data.totalUnits = 
config.Data.publication = True 
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' 
config.Data.outputDatasetTag = 'Run2015D-16Run2016B-v1_SilverJSON_2016_01_31_triggerTree_v1_211' 
#config.Data.outLFNDirBase = '/store/user/sflowers/treeMaker/Summer16_June13_Trigger/'
config.Data.outLFNDirBase = '/store/user/abdatta/' 
config.Data.ignoreLocality = True

config.section_("Site") 
config.Site.storageSite = 'T3_US_Cornell'
#config.Site.storageSite = 'T3_US_FNALLPC'
