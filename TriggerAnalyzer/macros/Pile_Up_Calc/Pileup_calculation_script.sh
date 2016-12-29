scram b -j 32

python pileupCalc.py -i '../../data/JSON/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON_unprescaled.txt' --inputLumiJSON '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt' --calcMode 'true' --minBiasXsec 69200 --maxPileupBin 75 --numPileupBins 75  output.root

root -l -q PU_data_hist_prod.C

g++ PU_hist_calc.cc -o PU_hist_calc

./PU_hist_calc

cp PU_weights.txt ../../data/PU_weight/


