# Trigger Analysis for Single Electron HLT

source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH="slc6_amd64_gcc530"

cmsrel CMSSW_8_0_26_patch1

cd CMSSW_8_0_26_patch1/src/

cmsenv

cd ../../

git clone https://github.com/abhisek92datta/TriggerRun2.git

cp TriggerRun2/TriggerRun_run_recipe.sh .

chmod 755 TriggerRun_run_recipe.sh

./TriggerRun_run_recipe.sh

rm -rf TriggerRun2.git

rm -rf TriggerRun_run_recipe.sh

cd CMSSW_8_0_26_patch1/src/TriggerRun2/TriggerAnalyzer/

RUN :

For MC :

To Run Locally:

1. put desired MC filename in : test/ trigger_analyzer_mc_25ns_cfg.py
2. cmsRun test/trigger_analyzer_mc_25ns_cfg.py

To Run on CRAB :

1. remove filename from test/trigger_analyzer_mc_25ns_cfg.py
2. put desired MC Dataset in : newCRAB_forMC.py
3. in newCRAB_forMC.py : config.Data.splitting = 'FileBased'
    
                       : config.Data.unitsPerJob = 10  (put desired number)

                       : config.Data.outLFNDirBase = '/store/user/abdatta/' (your own storage area)

                       : config.Site.storageSite = 'T3_US_FNALLPC' (your own storage site)

4. crab submit -c newCRAB_forMC.py

For DATA :

To Run Locally:

1. put desired DATA filename in : test/trigger_analyzer_data_80x_cfg.py
2. put latest LUMI filename in: test/trigger_analyzer_data_80x_cfg.py
3. cmsRun test/trigger_analyzer_data_80x_cfg.py

To Run on CRAB :

1. remove filename from test/trigger_analyzer_data_80x_cfg.py
2. put desired DATA dataset filename in : newCRAB_forData.py
3. put LUMI file name in : newCRAB_forData.py
4. in newCRAB_forData.py : config.Data.splitting = 'LumiBased'

                         : config.Data.unitsPerJob = 20  (put desired number)

                         : config.Data.outLFNDirBase = '/store/user/abdatta/' (your own storage area)

                         : config.Site.storageSite = 'T3_US_FNALLPC' (your own storage site)

5. crab submit -c newCRAB_forData.py

To check CRAB status : crab status -d <crab_output_directory_name>

To resubmit : crab resubmit -d <crab_output_directory_name>


After Ntuples created and stored in storage area, 
Plotting efficiencies using  : macros/Trigger_Express_June.C

cd macros/

1. Put desired file name with complete path in  : Trigger_Express_June.C 
2. Put the names of the triggers whose efficiency desired in : Trigger_Express_June.C 
3. scram b -j 32
4. root -l -q Trigger_Express_June.C
  
  or

  root -l -q Trigger_Express_June.C+'(-1,5,1)' 
 
  root -l -q Trigger_Express_June.C+'(-1,5,2)' 
 
  root -l -q Trigger_Express_June.C+'(-1,5,3)' 
 
  root -l -q Trigger_Express_June.C+'(-1,5,4)' 
 
  root -l -q Trigger_Express_June.C+'(-1,5,5)' 

  hadd Express_Trigger_June_Data.root Express_Trigger_June_Data_*


To put together the histograms in different root files in 1 root file :

root -l -q combine_plots.C
















