voms-proxy-init -voms cms

cd CMSSW_8_0_26_patch1/src/

git cms-init

git cms-merge-topic riga:badGlobalMuonTagger_fix

git cms-merge-topic rafaellopesdesa:RegressionCheckNegEnergy

git cms-merge-topic shervin86:Moriond2017_JEC_energyScales

cd EgammaAnalysis/ElectronTools/data

git clone git@github.com:ECALELFS/ScalesSmearings.git

cd ../../..

mv ../../TriggerRun2.git .

git clone https://github.com/abhisek92datta/MiniAOD.git

cd MiniAOD

git checkout CMSSW_8_0_8_ICHEP_Leptons

cd ..

scramv1 b -j 32

scramv1 b -j 32

scramv1 b -j 32

scramv1 b -j 32
