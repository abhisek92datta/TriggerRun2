from WMCore.Configuration import Configuration
config = Configuration() 
config.section_("General") 
config.General.requestName = 'triggerTree_SingleElectron_Run2016B'

config.section_("JobType") 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'test/trigger_analyzer_data_80x_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['data']
config.JobType.maxJobRuntimeMin = 3000

config.section_("Data") 
#Single Electron
#for B
config.Data.inputDataset = '/SingleElectron/Run2016B-03Feb2017_ver2-v2/MINIAOD'
#for C,D,E,F,G
#config.Data.inputDataset = '/SingleElectron/Run2016C-03Feb2017-v1/MINIAOD'
#for H
#config.Data.inputDataset = '/SingleElectron/Run2016H-03Feb2017_ver2-v1/MINIAOD'
#config.Data.inputDataset = '/SingleElectron/Run2016H-03Feb2017_ver3-v1/MINIAOD'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/' 
config.Data.allowNonValidInputDataset = True 
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'data/JSON/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_unprescaled.txt'
config.Data.unitsPerJob = 20
##config.Data.totalUnits = 
config.Data.publication = True 
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' 
config.Data.outputDatasetTag = 'triggerTree_SingleElectron_Run2016B'
config.Data.outLFNDirBase = '/store/user/abdatta/Trigger_Analysis/' 
config.Data.ignoreLocality = True

config.section_("Site") 
config.Site.storageSite = 'T3_US_FNALLPC'
