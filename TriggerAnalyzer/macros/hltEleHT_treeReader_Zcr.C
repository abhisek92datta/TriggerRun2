#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "Math/Interpolator.h"


#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif


#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "TriggerRun2/TriggerAnalyzer/interface/TriggerStudyEventVars.h"

#endif


//*****************************************************************************
typedef std::vector< TLorentzVector >          vecTLorentzVector;
typedef std::vector<int>                       vint;
typedef std::vector<double>                    vdouble;
typedef std::vector<std::vector<double> >      vvdouble;

double reweightPU( int nPU );
float DeltaR(float eta1,float phi1,float eta2,float phi2);

//*****************************************************************************

void hltEleHT_treeReader_Zcr( int insample=1, int maxNentries=-1, int Njobs=1, int jobN=1, double intLumi=-1 ) {

  ////
  std::cout << "   ===> load the root files! " << std::endl;

  std::string sampleType = ( insample>=0 ) ? "mc" : "data";
  std::string str_jobN;
  std::stringstream stream;
  stream << jobN;
  str_jobN = stream.str();

  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
  double mySample_xSec_ = 1.;
  double mySample_nGen_ = 1.;
  std::string mySample_sampleName_ = "delete";
  std::string mySample_inputDir_ = "";
  std::vector<std::string> mySample_inputDirs_;
  if( insample==2500 ){
    mySample_xSec_ = 831.76;//https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
    mySample_nGen_ = 25407710;//25357774;//25446993;
    mySample_sampleName_ = "TTJets";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1_yggdrasilTree_v3/150105_181019/0000/");
    //mySample_inputDir_ = "/uscms_data/d2/dpuigh/TTH/miniAOD/CMSSW_7_2_3/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/");
  }
  else if( insample==2325 ){
    mySample_xSec_ = 3*2008.4; 
    mySample_nGen_ = 19310834;//59000;//28825132;
    mySample_sampleName_ = "DYJets_M50_13TeV_Spring15_Asympt25ns";
    // mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3_triggerTree_v1/150925_213046/0000/");
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_74X_mcRun2_asymptotic_v2-v3_triggerTree_v1/151009_135127/0000/");
  }
  else if( insample==2300 ){
    mySample_xSec_ = 3*2008.4; 
    mySample_nGen_ = 13351018;//19925500;//59000;//28825132;
    mySample_sampleName_ = "DYJets_M50_13TeV_Spring15";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2_yggdrasilTree_v1/150716_231755/0000/");
  }
  else if( insample==2359 ){
    mySample_xSec_ = 3*2008.4; 
    mySample_nGen_ = 200173;//299269;//13351018;//19925500;//59000;//28825132;
    mySample_sampleName_ = "DYJets_M50_13TeV_Spring15_StartupFlat10to50bx50Raw";
    mySample_inputDirs_.push_back("/uscms_data/d2/dpuigh/TTH/triggerRun2/CMSSW_7_4_7/src/TreeMaker_DYJetsToLL_M50_StartupFlat10to50bx50Raw_MCRUN2_74_V8_StartupFlat10to50bx50Raw_MCRUN2_74_V8/");
  }
  else if( insample==2358 ){
    mySample_xSec_ = 3*2008.4; 
    mySample_nGen_ = 200173;//299269;//13351018;//19925500;//59000;//28825132;
    mySample_sampleName_ = "DYJets_M50_13TeV_Spring15_StartupFlat10to50bx50Raw_OriginalRCTcalib";
    mySample_inputDirs_.push_back("/uscms_data/d2/dpuigh/TTH/triggerRun2/CMSSW_7_4_7/src/TreeMaker_DYJetsToLL_M50_StartupFlat10to50bx50Raw_MCRUN2_74_V8_StartupFlat10to50bx50Raw_MCRUN2_74_V8_OriginalRCTcalib/");
  }
  else if( insample==2400 ){
    mySample_xSec_ = 20508.9;  
    mySample_nGen_ = 10017462;
    mySample_sampleName_ = "WJetsToLNu";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WJetsToLNu_13TeV-madgraph-pythia8-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1_yggdrasilTree_v1/150217_010312/0000/");
  }
  else if( insample==2524 ){
    mySample_xSec_ = 1.152;  
    mySample_nGen_ = 246521;
    mySample_sampleName_ = "TTWJets";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTWJets_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1_yggdrasilTree_v1/150217_005352/0000/");
  }
  else if( insample==2523 ){
    mySample_xSec_ = 2.232;  
    mySample_nGen_ = 249275;
    mySample_sampleName_ = "TTZJets";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTZJets_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1_yggdrasilTree_v1/150217_005607/0000/");
  }
  else if( insample==2510 ){
    mySample_xSec_ = 2.232;  
    mySample_nGen_ = 500000;
    mySample_sampleName_ = "TToLeptons_s";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1_yggdrasilTree_v1/150217_005853/0000/");
  }
  else if( insample==2511 ){
    mySample_xSec_ = 2.232;  
    mySample_nGen_ = 250000;
    mySample_sampleName_ = "TBarToLeptons_s";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1_yggdrasilTree_v1/150217_004555/0000/");
  }
  else if( insample==2512 ){
    mySample_xSec_ = 2.232;  
    mySample_nGen_ = 3991000;
    mySample_sampleName_ = "TToLeptons_t";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1_yggdrasilTree_v1/150217_005929/0000/");
  }
  else if( insample==2513 ){
    mySample_xSec_ = 2.232;  
    mySample_nGen_ = 1999800;
    mySample_sampleName_ = "TBarToLeptons_t";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1_yggdrasilTree_v1/150217_004732/0000/");
  }
  else if( insample==2514 ){
    mySample_xSec_ = 35.6;  
    mySample_nGen_ = 986100;
    mySample_sampleName_ = "T_tW_DR";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1_yggdrasilTree_v1/150217_010006/0000/");
  }
  else if( insample==2515 ){
    mySample_xSec_ = 35.6;  
    mySample_nGen_ = 971800;
    mySample_sampleName_ = "Tbar_tW_DR";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1_yggdrasilTree_v1/150217_010035/0000/");
  }
  else if( insample==9125 ){
    mySample_xSec_ = 0.5085 * 1.0;// YR3 * BR(all)  
    mySample_nGen_ = 199700;//199000;
    mySample_sampleName_ = "TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola_PU20bx25_tsg_PHYS14_25_V1-v2";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v2_v1_yggdrasilTree_v1/150217_004834/0000/");
  }
  else if( insample==-13 ){
    mySample_sampleName_ = "SingleElectron_Run2015D_PromptReco_254231_258158";
    // //mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleElectron/Run2015B-PromptReco-v1_DCSONLY_251244_251562_yggdrasilTree_v1/150713_201455/0000/");
    // //mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleElectron/Run2015C-PromptReco-v1_GoldenJSON_254231_256869_triggerTree_v1/150925_212953/0000/");
    // mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleElectron/Run2015D-PromptReco-v3_GoldenJSON_254231_256869_triggerTree_v1/150928_210759/0000/");
    // mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleElectron/Run2015D-PromptReco-v3_GoldenJSON_254231_257599_triggerTree_v1/151005_145059/0000/");
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleElectron/Run2015D-PromptReco-v3_GoldenJSON_254231_258158_triggerTree_v1/151012_183037/0000/");
  }




  //std::string treefilename = mySample_inputDir_ + "trigger_analyzer*.root";
  //std::string treefilename = mySample_inputDir_ + "yggdrasil_treeMaker_99.root";
  //std::string treefilename = "/uscms_data/d2/dpuigh/TTH/triggerRun2/CMSSW_7_4_12/src/TriggerRun2/TriggerAnalyzer/trigger_analyzer.root";

  std::string s_end = "histo_" + str_jobN + ".root";
  if( Njobs==1 ) s_end = "histo.root";

  std::string histofilename = Form("HistoFiles/hltEleHT_treeReader_Zcr_%s_%s", mySample_sampleName_.c_str(), s_end.c_str());

  TChain *chain = new TChain("triggeranalzyer/triggerTree");
  for( int iFile=0; iFile<int(mySample_inputDirs_.size()); iFile++ ){
    std::string treefilename = mySample_inputDirs_[iFile] + "trigger_analyzer*.root";
    std::cout << "  treefilename " << iFile << ": " << treefilename.c_str() << std::endl;
    chain->Add(treefilename.c_str());
  }

  std::cout << "  histofilename = " << histofilename.c_str() << std::endl;

  //////////////////////////////////////////////////////////////////////////
  ///  Tree branches/leaves
  //////////////////////////////////////////////////////////////////////////

  triggerStudyEventVars *eve=0;
  chain->SetBranchAddress("eve.", &eve );

  //////////////////////////////////////////////////////////////////////////
  ///  Histogram making
  //////////////////////////////////////////////////////////////////////////


  TFile histofile(histofilename.c_str(),"recreate");
  histofile.cd();


  bool verbose = false;

  //////////////////////////////////////////////////////////////////////////
  ///  Histograms
  //////////////////////////////////////////////////////////////////////////


  std::vector<int> runsInData;
  runsInData.push_back(254231);
  runsInData.push_back(254232);
  runsInData.push_back(254790);
  runsInData.push_back(254852);
  runsInData.push_back(254879);
  runsInData.push_back(254906);
  runsInData.push_back(254907);
  runsInData.push_back(254914);

  runsInData.push_back(256630);
  runsInData.push_back(256673);
  runsInData.push_back(256674);
  runsInData.push_back(256675);
  runsInData.push_back(256677);
  runsInData.push_back(256729);
  runsInData.push_back(256734);
  runsInData.push_back(256801);
  runsInData.push_back(256842);
  runsInData.push_back(256843);

  int Nruns = int(runsInData.size());
  TH2D* h_run_passL1SingleEG25 = new TH2D("h_run_passL1SingleEG25",";Run;Pass L1SingleEG25", Nruns, 0, Nruns, 120, -0.1, 1.1 );
  for( unsigned int iRun=0; iRun<Nruns; iRun++ ){
    TString runNum = Form("%d",runsInData[iRun]);
    h_run_passL1SingleEG25->GetXaxis()->SetBinLabel(iRun+1,runNum);
  }

  double window = 10.;
  double pdgZmass = 91.1876;
  double MinMass  = pdgZmass - window;
  double MaxMass  = pdgZmass + window;

  TH1::SetDefaultSumw2();

  TH1D* h_numEvents = new TH1D("h_numEvents",";Number of events", 2, 0, 2 );
  TH1D* h_numEvents_wgt = new TH1D("h_numEvents_wgt",";Predicted number of events", 2, 0, 2 );

  TH1D* h_numPVs = new TH1D("h_numPVs",";Number of PVs", 50, 0-0.5, 50-0.5 );
  TH1D* h_numPVs_PUwgt = new TH1D("h_numPVs_PUwgt",";Number of PVs", 50, 0-0.5, 50-0.5 );


  int NumCuts = 8;
  TH1D* h_event_selection  = new TH1D("h_event_selection",";cut", NumCuts, 0, NumCuts );
  h_event_selection->GetXaxis()->SetBinLabel(1,"All");
  h_event_selection->GetXaxis()->SetBinLabel(2,"HLT Ele27WP");
  h_event_selection->GetXaxis()->SetBinLabel(3,">=1 tag electron");
  h_event_selection->GetXaxis()->SetBinLabel(4,">=1 tag-probe pair");
  h_event_selection->GetXaxis()->SetBinLabel(5,"Opposite charge");
  h_event_selection->GetXaxis()->SetBinLabel(6,"|M(ee) - MZ|<10");
  h_event_selection->GetXaxis()->SetBinLabel(7,"L1 HTT > 125");
  h_event_selection->GetXaxis()->SetBinLabel(8,"HLT HT200");


  TH1D* h_numTag = new TH1D("h_numTag",";number of tag electrons", 5, 0, 5 );
  TH1D* h_numProbe = new TH1D("h_numProbe",";number of probe electrons", 5, 0, 5 );
  TH1D* h_numTagProbe = new TH1D("h_numTagProbe",";number of electron tag-probe pairs", 5, 0, 5 );

  TH1D* h_diele_mass   = new TH1D("h_diele_mass",";tag-probe mass", 200, 0., 200. );
  TH1D* h_diele_deltaR = new TH1D("h_diele_deltaR",";tag-probe #DeltaR", 61, 0., 6.1 );

  TH1D* h_deltaR_jet_ele = new TH1D("h_deltaR_jet_ele",";#DeltaR(jet,ele)", 61, 0., 6.1 );

  TH1D* h_diele_mass_closestZmass   = new TH1D("h_diele_mass_closestZmass",";tag-probe mass", 200, 0., 200. );


  TH1D* h_L1HTT = new TH1D("h_L1HTT",";L1 HTT", 300, 0., 300. );

  TH1D* h_HT30 = new TH1D("h_HT30",";H_{T} (p_{T}>30)", 500, 0., 500. );
  TH1D* h_HT30_L1HTT125 = new TH1D("h_HT30_L1HTT125",";H_{T} (p_{T}>30)", 500, 0., 500. );
  TH1D* h_HT30_L1HTT125_passHLTEle27HT200 = new TH1D("h_HT30_L1HTT125_passHLTEle27HT200",";H_{T} (p_{T}>30)", 500, 0., 500. );

  TH1D* h_HT30_v2 = new TH1D("h_HT30_v2",";H_{T} (p_{T}>30)", 500, 0., 500. );
  TH1D* h_HT30_v2_L1HTT = new TH1D("h_HT30_v2_L1HTT",";H_{T} (p_{T}>30)", 500, 0., 500. );
  TH1D* h_HT30_v2_L1HTT_HLTHT = new TH1D("h_HT30_v2_L1HTT_HLTHT",";H_{T} (p_{T}>30)", 500, 0., 500. );


  TH1D* h_met_pt = new TH1D("h_met_pt",";MET", 300, 0., 300. );

  TH1D* h_jet_pt = new TH1D("h_jet_pt",";jet p_{T}", 300, 0., 300. );
  TH1D* h_jet_1_pt = new TH1D("h_jet_1_pt",";jet 1 p_{T}", 300, 0., 300. );
  TH1D* h_jet_2_pt = new TH1D("h_jet_2_pt",";jet 2 p_{T}", 300, 0., 300. );
  TH1D* h_jet_3_pt = new TH1D("h_jet_3_pt",";jet 3 p_{T}", 300, 0., 300. );
  TH1D* h_jet_4_pt = new TH1D("h_jet_4_pt",";jet 4 p_{T}", 300, 0., 300. );

  TH1D* h_jet_eta = new TH1D("h_jet_eta",";jet #eta", 52, -2.5, 2.5 );
  TH1D* h_jet_phi = new TH1D("h_jet_phi",";jet #phi", 64, -3.2, 3.2 );

  TH1D* h_jet_csv = new TH1D("h_jet_csv",";jet CSV", 152, -0.5, 1.02 );

  TH1D* h_numJet = new TH1D("h_numJet",";Number of Jets", 6, 0-0.5, 6-0.5 );
  TH1D* h_numBtag = new TH1D("h_numBtag",";Number of b-tagged Jets", 4, 0-0.5, 4-0.5 );

  TH1D* h_numL1EG25 = new TH1D("h_numL1EG25",";Number of L1EG25", 5, 0-0.5, 5-0.5 );



  TH1D* h_l1EG25_pt = new TH1D("h_l1EG25_pt",";L1 EG p_{T}", 80, 0., 80. );
  TH1D* h_l1EG25_eta = new TH1D("h_l1EG25_eta",";L1 EG #eta", 52, -2.5, 2.5 );
  TH1D* h_l1EG25_phi = new TH1D("h_l1EG25_phi",";L1 EG #phi", 64, -3.2, 3.2 );

  TH1D* h_l1EG25_1_pt = new TH1D("h_l1EG25_1_pt",";L1 EG 1 p_{T}", 80, 0., 80. );
  TH1D* h_l1EG25_2_pt = new TH1D("h_l1EG25_2_pt",";L1 EG 2 p_{T}", 80, 0., 80. );

  TH1D* h_l1EG25_1_eta = new TH1D("h_l1EG25_1_eta",";L1 EG 1 #eta", 52, -2.5, 2.5 );
  TH1D* h_l1EG25_2_eta = new TH1D("h_l1EG25_2_eta",";L1 EG 2 #eta", 52, -2.5, 2.5 );


  TH1D* h_numhltL1EG25 = new TH1D("h_numhltL1EG25",";Number of hltL1EG25", 5, 0-0.5, 5-0.5 );

  TH1D* h_hltL1EG25_pt  = new TH1D("h_hltL1EG25_pt",";hltL1EG25 EG p_{T}", 80, 0., 80. );
  TH1D* h_hltL1EG25_eta = new TH1D("h_hltL1EG25_eta",";hltL1EG25 EG #eta", 52, -2.5, 2.5 );
  TH1D* h_hltL1EG25_phi = new TH1D("h_hltL1EG25_phi",";hltL1EG25 EG #phi", 64, -3.2, 3.2 );

  TH1D* h_numhltL1IsoEG22OrEG25 = new TH1D("h_numhltL1IsoEG22OrEG25",";Number of hltL1IsoEG22OrEG25", 5, 0-0.5, 5-0.5 );

  TH1D* h_hltL1IsoEG22OrEG25_pt  = new TH1D("h_hltL1IsoEG22OrEG25_pt",";hltL1IsoEG22OrEG25 EG p_{T}", 80, 0., 80. );
  TH1D* h_hltL1IsoEG22OrEG25_eta = new TH1D("h_hltL1IsoEG22OrEG25_eta",";hltL1IsoEG22OrEG25 EG #eta", 52, -2.5, 2.5 );
  TH1D* h_hltL1IsoEG22OrEG25_phi = new TH1D("h_hltL1IsoEG22OrEG25_phi",";hltL1IsoEG22OrEG25 EG #phi", 64, -3.2, 3.2 );



  // int Nptbins = 17;
  // double ptbins[] = { 10, 15, 20, 30, 40, 50, 60, 70, 80, 100, 125, 150, 200, 250, 300, 400, 500, 600 };
  int Nptbins = 17;
  double ptbins[] = { 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 125, 150, 200, 300, 500 };
  double ptMax = 200.;

  int Netabins = 12;
  double etabins[] = { -2.5, -2.1, -1.65, -1.2, -0.9, -0.45, 0, 0.45, 0.9, 1.2, 1.65, 2.1, 2.5 };

  int Nphibins = 16;
  double phibins[] = { -3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2 };

  int Npvbins = 14;
  double pvbins[] = { 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17, 35 };

  TH1D* h_probe_pt_pL1T_limRun = new TH1D("h_probe_pt_pL1T_limRun",";probe p_{T}", Nptbins, ptbins );
  TH1D* h_probe_pt_all_limRun  = new TH1D("h_probe_pt_all_limRun", ";probe p_{T}", Nptbins, ptbins );

  TH1D* h_probe_pt_pHLT_fL1T = new TH1D("h_probe_pt_pHLT_fL1T",";probe p_{T}", Nptbins, ptbins );
  TH1D* h_probe_pt_pL1TOR = new TH1D("h_probe_pt_pL1TOR",";probe p_{T}", Nptbins, ptbins );

  TH1D* h_probe_pt_pHLTall = new TH1D("h_probe_pt_pHLTall",";probe p_{T}", Nptbins, ptbins );

  TH1D* h_probe_pt_pL1T = new TH1D("h_probe_pt_pL1T",";probe p_{T}", Nptbins, ptbins );
  TH1D* h_probe_pt_fL1T = new TH1D("h_probe_pt_fL1T",";probe p_{T}", Nptbins, ptbins );
  TH1D* h_probe_pt_pHLT = new TH1D("h_probe_pt_pHLT",";probe p_{T}", Nptbins, ptbins );
  TH1D* h_probe_pt_fHLT = new TH1D("h_probe_pt_fHLT",";probe p_{T}", Nptbins, ptbins );
  TH1D* h_probe_pt_all  = new TH1D("h_probe_pt_all", ";probe p_{T}", Nptbins, ptbins );

  TH1D* h_probe_eta_pL1TOR = new TH1D("h_probe_eta_pL1TOR",";probe #eta", Netabins, etabins );

  TH1D* h_probe_eta_pL1T = new TH1D("h_probe_eta_pL1T",";probe #eta", Netabins, etabins );
  TH1D* h_probe_eta_fL1T = new TH1D("h_probe_eta_fL1T",";probe #eta", Netabins, etabins );
  TH1D* h_probe_eta_pHLT = new TH1D("h_probe_eta_pHLT",";probe #eta", Netabins, etabins );
  TH1D* h_probe_eta_fHLT = new TH1D("h_probe_eta_fHLT",";probe #eta", Netabins, etabins );
  TH1D* h_probe_eta_all  = new TH1D("h_probe_eta_all", ";probe #eta", Netabins, etabins );

  TH1D* h_probe_phi_pL1TOR = new TH1D("h_probe_phi_pL1TOR",";probe #phi", Nphibins, phibins );

  TH1D* h_probe_phi_pL1T = new TH1D("h_probe_phi_pL1T",";probe #phi", Nphibins, phibins );
  TH1D* h_probe_phi_fL1T = new TH1D("h_probe_phi_fL1T",";probe #phi", Nphibins, phibins );
  TH1D* h_probe_phi_pHLT = new TH1D("h_probe_phi_pHLT",";probe #phi", Nphibins, phibins );
  TH1D* h_probe_phi_fHLT = new TH1D("h_probe_phi_fHLT",";probe #phi", Nphibins, phibins );
  TH1D* h_probe_phi_all  = new TH1D("h_probe_phi_all", ";probe #phi", Nphibins, phibins );

  TH1D* h_probe_numPVs_pL1TOR = new TH1D("h_probe_numPVs_pL1TOR",";Number of PVs", Npvbins, pvbins );

  TH1D* h_probe_numPVs_pL1T = new TH1D("h_probe_numPVs_pL1T",";Number of PVs", Npvbins, pvbins );
  TH1D* h_probe_numPVs_fL1T = new TH1D("h_probe_numPVs_fL1T",";Number of PVs", Npvbins, pvbins );
  TH1D* h_probe_numPVs_pHLT = new TH1D("h_probe_numPVs_pHLT",";Number of PVs", Npvbins, pvbins );
  TH1D* h_probe_numPVs_fHLT = new TH1D("h_probe_numPVs_fHLT",";Number of PVs", Npvbins, pvbins );
  TH1D* h_probe_numPVs_all  = new TH1D("h_probe_numPVs_all", ";Number of PVs", Npvbins, pvbins );


  //////////////////////////////////////////////////////////////////////////
  /////
  //////////////////////////////////////////////////////////////////////////

  int numEvents_all = 0;
  int numEvents_wgt_gen = 0;
  double numEvents_wgt_lumi = 0;
  double numEvents_wgt = 0;

  int nentries = chain->GetEntries();
  std::cout << "\n\t Number of entries = " << nentries << std::endl;
  std::cout << "\t Max number of entries = " << maxNentries << std::endl;
  std::cout << "\n" << std::endl;

  int use_nentries = std::max( maxNentries, nentries);

  int NeventsPerJob = int( double(use_nentries)/double(Njobs) + 0.000001 ) + 1;

  int firstEvent = (jobN-1)*NeventsPerJob + 1;
  int lastEvent  = firstEvent + NeventsPerJob;
  if( jobN==Njobs ) lastEvent = -1;
  if( jobN==1 ) firstEvent = 0;

  int cnt = 0;
  std::cout << "========  Starting Event Loop  ========" << std::endl;
  for (Long64_t ievt=0; ievt<chain->GetEntries();ievt++) {    //Long64_t
    cnt++;
    if( ievt<firstEvent ) continue;
    if( ievt==lastEvent ) break;

    if( ievt==1 )        std::cout << "     Event " << ievt << std::endl;
    if( ievt%10000==0 && ievt!=1 ) std::cout << "           " << ievt << "\t" 
					     << int(double(ievt-firstEvent)/double(NeventsPerJob)*100) << "% done" << std::endl;

    //if( ievt==(maxNentries+1) ) break;
    if( ievt==(maxNentries+1) && ievt!=0 ) break;

    chain->GetEntry(ievt);

    // vint hlt_accept = eve->hlt_accept_;
    // vstring hlt_name = eve->hlt_name_;

    // int Nhlt = int( hlt_accept.size() );
    // for( int iHLT=0; iHLT<Nhlt; iHLT++ ){
    //   std::string name = hlt_name[iHLT];
    //   int accept = hlt_accept[iHLT];
    //   printf("\t %90s,\t %5d \n", name.c_str(), accept);

    //   if( accept ){
    // 	TAxis * axis = h_hlt->GetXaxis();
    // 	if( !axis ) continue;
    // 	int bin_num = axis->FindBin(name.c_str());
    // 	int bn = bin_num - 1;
    // 	h_hlt->Fill(bn, 1);
    //   }
    // }

    int run = eve->run_;

    int numPVs = eve->numPVs_;

    double scalePU = ( insample < 0 ) ? 1. : reweightPU(numPVs);

    //// --------- various weights: PU, topPt, triggerSF, leptonSF...
    // double  wgt_topPtSF = eve->wgt_topPt_;
    double Xsec = mySample_xSec_;//eve->wgt_xs_;
    double nGen = ( maxNentries>0 ) ? maxNentries : mySample_nGen_;//eve->wgt_nGen_;
    if( maxNentries==1000000 && insample==2300 ) nGen = 670121;
    double lumi = ( intLumi > 0 ) ? intLumi : 10000 ;

    double wgt_gen = ( insample > 0 ) ? eve->wgt_generator_ : 1;
    wgt_gen = ( wgt_gen > 0 ) ? 1. : -1.;
    double wgt_lumi = ( insample > 0 ) ? lumi * (Xsec/nGen) : 1;//"weight_PU*topPtWgt*osTriggerSF*lepIDAndIsoSF*"; // various weights

    double wgt = wgt_gen * wgt_lumi;

    h_numEvents->Fill(0.5,1.);
    h_numEvents->Fill(1.5,wgt_gen);

    h_numEvents_wgt->Fill(0.5,wgt_lumi);
    h_numEvents_wgt->Fill(1.5,wgt);


    numEvents_all += 1.;
    numEvents_wgt_gen += wgt_gen;
    numEvents_wgt_lumi += wgt_lumi;
    numEvents_wgt += wgt;


    ///////////////////
    ////// selections
    ///////////////////

    // Don't use data runs before 256630 -> L1_SingleEG25 is prescaled
    if( run<256630 && insample<0 ) continue;


    int ncut=0;
    h_event_selection->Fill(0.5+ncut++, wgt); // all

    bool pass_trigger = false;
    if( insample<0 ) pass_trigger = ( eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_==1 );
    else             pass_trigger = ( eve->pass_HLT_Ele27_WP85_Gsf_v_==1 );
    //else             pass_trigger = ( eve->pass_HLT_Ele27_eta2p1_WP75_Gsf_v_==1 );
    

    if( !pass_trigger ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // trigger


    h_numPVs->Fill(numPVs,wgt);


    // Remove PU reweighting for now
    //scalePU = 1;

    wgt *= scalePU;

    h_numPVs_PUwgt->Fill(numPVs,wgt);

    double HT30=0;
    for( int iJet=0; iJet<int(eve->jet_nocc_pt_.size()); iJet++ ){
      double pt  = eve->jet_nocc_pt_[iJet];
      double eta = eve->jet_nocc_eta_[iJet];

      if( !(pt>30. && fabs(eta)<3.0) ) continue;

      HT30 += pt;
    }

    double MET = eve->pfMET_pt_;



    ////////////


    double mass_leplep = eve->mass_leplep_;
    double dR_leplep = eve->dR_leplep_;
    int oppositeLepCharge = eve->oppositeLepCharge_;


    std::vector<int> ind_tag_ele;
    std::vector<int> ind_probe_ele;
    for( int iLep=0; iLep<int(eve->lepton_pt_.size()); iLep++ ){
      bool isElectron = ( eve->lepton_isMuon_[iLep]==0 );
      if( !isElectron ) continue;

      bool isEB = ( eve->lepton_isEB_[iLep]==1 );

      bool isSpring15M = ( eve->lepton_isSpring15M_[iLep]==1 );

      if( !isSpring15M ) continue;
      if( eve->lepton_inCrack_[iLep]==1 ) continue;


      bool matchHLT = false;
      if( insample<0 ) matchHLT = ( eve->lepton_ele_matchHLT_hltEle27WPLooseGsfTrackIsoFilter_[iLep]==1 );
      else             matchHLT = ( eve->lepton_ele_matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter_[iLep]==1 );

      double pt  = eve->lepton_pt_[iLep];
      double eta = eve->lepton_eta_[iLep];

      if( pt > 30 && abs(eta)<2.1 && matchHLT ) ind_tag_ele.push_back(iLep);

      if( pt > 20 && abs(eta)<2.5 ) ind_probe_ele.push_back(iLep);
    }


    int numTag = int(ind_tag_ele.size());
    int numProbe = int(ind_probe_ele.size());

    h_numTag->Fill(numTag,wgt);
    h_numProbe->Fill(numProbe,wgt);

    int numTagProbe=0;
    bool hasPairInZWindow=false;
    bool OppositeSign=false;
    double diele_mass_closestZmass = -99;
    for( int iTag=0; iTag<numTag; iTag++ ){
      int ind1 = ind_tag_ele[iTag];

      TLorentzVector myLep1;
      myLep1.SetPtEtaPhiE( eve->lepton_pt_[ind1], eve->lepton_eta_[ind1], eve->lepton_phi_[ind1], eve->lepton_energy_[ind1] );

      int charge1 = eve->lepton_trkCharge_[ind1];

      for( int iProbe=0; iProbe<numProbe; iProbe++ ){
	int ind2 = ind_probe_ele[iProbe];

	if( ind2==ind1 ) continue;

	TLorentzVector myLep2;
	myLep2.SetPtEtaPhiE( eve->lepton_pt_[ind2], eve->lepton_eta_[ind2], eve->lepton_phi_[ind2], eve->lepton_energy_[ind2] );

	int charge2 = eve->lepton_trkCharge_[ind2];

	if( charge1==charge2 ) continue;

	OppositeSign = true;

	TLorentzVector diele = myLep1 + myLep2;
	double mass = diele.M();
	double deltaR = myLep1.DeltaR(myLep2);

	h_diele_mass->Fill(mass,wgt);
	h_diele_deltaR->Fill(deltaR,wgt);

	if( fabs(mass - pdgZmass) < fabs(diele_mass_closestZmass - pdgZmass) ) diele_mass_closestZmass = mass;
	
	if( !(mass > MinMass && mass < MaxMass) ) continue;
	hasPairInZWindow = true;

	numTagProbe++;

	double pt = myLep2.Pt();
	double eta = eve->lepton_scEta_[ind2];//myLep2.Eta();
	double phi = myLep2.Phi();


	bool passHLT=false;
	bool passL1T=false;
	bool passL1TOR=false;
	for( int iHLT=0; iHLT<int(eve->hltEle27WP85Gsf_pt_.size()); iHLT++ ){
	  std::string hltname = eve->hltEle27WP85Gsf_filter_[iHLT];

	  double hlteta = eve->hltEle27WP85Gsf_eta_[iHLT];
	  double hltphi = eve->hltEle27WP85Gsf_phi_[iHLT];

	  double deltaR = DeltaR( eta, phi, hlteta, hltphi );
	  if( deltaR<0.5 ){
	    if( hltname=="hltL1sL1SingleEG25" ) passL1T = true;
	    if( hltname=="hltL1sL1SingleIsoEG22erOrSingleEG25" ) passL1TOR = true;
	    if( hltname=="hltEle27WPLooseGsfTrackIsoFilter" ||
		hltname=="hltL1EG25Ele27WP85GsfTrackIsoFilter" ) passHLT = true;
	  }
	}

	if( fabs(eta) < 2.1 && pt > 30 ){
	  double does_passL1T = (passL1T) ? 1.0 : 0;
	  if( h_run_passL1SingleEG25 ){
	    TString useRunNum = Form("%d",run);
	    TAxis * axis = h_run_passL1SingleEG25->GetXaxis();
	    if( !axis ) continue;
	    int bin_num = axis->FindBin(useRunNum);
	    int bn = bin_num - 1;
	    h_run_passL1SingleEG25->Fill(bn, does_passL1T, wgt);
	  }
	}

	if( fabs(eta) < 2.1 ){
	  h_probe_pt_all->Fill(pt, wgt);
	  if( passHLT ) h_probe_pt_pHLTall->Fill(pt, wgt);
	  if( passHLT && !passL1T ) h_probe_pt_pHLT_fL1T->Fill(pt, wgt);
	  if( passL1TOR ) h_probe_pt_pL1TOR->Fill(pt, wgt);
	  if( passL1T ){
	    h_probe_pt_pL1T->Fill(pt, wgt);
	    if( passHLT ) h_probe_pt_pHLT->Fill(pt, wgt);
	    else          h_probe_pt_fHLT->Fill(pt, wgt);
	  }
	  else{
	    h_probe_pt_fL1T->Fill(pt, wgt);
	  }

	  if( run>=256630 ){
	    h_probe_pt_all_limRun->Fill(pt, wgt);
	    if( passL1T ) h_probe_pt_pL1T_limRun->Fill(pt, wgt);
	  }
	}

	if( pt > 30 ){
	  h_probe_eta_all->Fill(eta, wgt);
	  if( passL1TOR ) h_probe_eta_pL1TOR->Fill(eta, wgt);
	  if( passL1T ){
	    h_probe_eta_pL1T->Fill(eta, wgt);
	    if( passHLT ) h_probe_eta_pHLT->Fill(eta, wgt);
	    else          h_probe_eta_fHLT->Fill(eta, wgt);
	  }
	  else{
	    h_probe_eta_fL1T->Fill(eta, wgt);
	  }
	}

	if( pt > 30 && fabs(eta) < 2.1 ){
	  h_probe_phi_all->Fill(phi, wgt);
	  if( passL1TOR ) h_probe_phi_pL1TOR->Fill(phi, wgt);
	  if( passL1T ){
	    h_probe_phi_pL1T->Fill(phi, wgt);
	    if( passHLT ) h_probe_phi_pHLT->Fill(phi, wgt);
	    else          h_probe_phi_fHLT->Fill(phi, wgt);
	  }
	  else{
	    h_probe_phi_fL1T->Fill(phi, wgt);
	  }
	}
	if( pt > 30 && fabs(eta) < 2.1 ){
	  h_probe_numPVs_all->Fill(numPVs, wgt);
	  if( passL1TOR ) h_probe_numPVs_pL1TOR->Fill(numPVs, wgt);
	  if( passL1T ){
	    h_probe_numPVs_pL1T->Fill(numPVs, wgt);
	    if( passHLT ) h_probe_numPVs_pHLT->Fill(numPVs, wgt);
	    else          h_probe_numPVs_fHLT->Fill(numPVs, wgt);
	  }
	  else{
	    h_probe_numPVs_fL1T->Fill(numPVs, wgt);
	  }
	}
      }
    }

    h_numTagProbe->Fill(numTagProbe,wgt);

    h_diele_mass_closestZmass->Fill(diele_mass_closestZmass,wgt);


    if( !(numTag >= 1) ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // >= 1 tag

    if( !(numProbe >= 2) ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // >= 1 tag-probe pair

    if( !(OppositeSign) ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // opposite sign

    if( !(hasPairInZWindow) ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // in Z window



    vdouble jetPts;
    int numJet = 0;
    int numBtag = 0;
    for( int iJet=0; iJet<int(eve->jet_pt_.size()); iJet++ ){
      TLorentzVector myJet;
      myJet.SetPtEtaPhiE( eve->jet_pt_[iJet], eve->jet_eta_[iJet], eve->jet_phi_[iJet], eve->jet_energy_[iJet] );

      double pt  = eve->jet_pt_[iJet];
      double eta = eve->jet_eta_[iJet];
      double phi = eve->jet_phi_[iJet];

      if( !(pt>30. && fabs(eta)<2.4) ) continue;

      bool hasEleOverlap = false;
      for( int iEle=0; iEle<int(ind_probe_ele.size()); iEle++ ){
	int ind = ind_probe_ele[iEle];

	TLorentzVector myLep;
	myLep.SetPtEtaPhiE( eve->lepton_pt_[ind], eve->lepton_eta_[ind], eve->lepton_phi_[ind], eve->lepton_energy_[ind] );

	double dR = myJet.DeltaR(myLep);
	h_deltaR_jet_ele->Fill(dR,wgt);

	if( dR < 0.3 ) hasEleOverlap = true;
      }

      if( hasEleOverlap ) continue;

      double csv =  eve->jet_csv_[iJet];
      if( csv<0 && csv>-9 ) csv = -0.2;
      else if( csv < -5 )   csv = -0.4;

      h_jet_pt->Fill(pt,wgt);
      h_jet_eta->Fill(eta,wgt);
      h_jet_phi->Fill(phi,wgt);
      h_jet_csv->Fill(csv,wgt);

      jetPts.push_back(pt);

      numJet++;
      if( csv > 0.890 ) numBtag++;
    }

    h_met_pt->Fill(MET,wgt);

    std::sort(jetPts.begin(), jetPts.end());
    std::reverse(jetPts.begin(), jetPts.end());

    if( jetPts.size()>=1 ) h_jet_1_pt->Fill(jetPts[0],wgt);
    if( jetPts.size()>=2 ) h_jet_2_pt->Fill(jetPts[1],wgt);
    if( jetPts.size()>=3 ) h_jet_3_pt->Fill(jetPts[2],wgt);
    if( jetPts.size()>=4 ) h_jet_4_pt->Fill(jetPts[3],wgt);

    h_numJet->Fill(numJet,wgt);
    h_numBtag->Fill(numBtag,wgt);

    double leadL1EG_pt = -99;
    double leadL1EG_eta = -99;

    double subleadL1EG_pt = -99;
    double subleadL1EG_eta = -99;

    int numL1EG=0;
    for( int iEG=0; iEG<int(eve->hltL1SingleEG25_pt_.size()); iEG++ ){
      double pt = eve->hltL1SingleEG25_pt_[iEG];
      double eta = eve->hltL1SingleEG25_eta_[iEG];
      double phi = eve->hltL1SingleEG25_phi_[iEG];

      h_l1EG25_pt->Fill(pt,wgt);
      h_l1EG25_eta->Fill(eta,wgt);
      h_l1EG25_phi->Fill(phi,wgt);

      if( pt > leadL1EG_pt ){
	subleadL1EG_pt = leadL1EG_pt;
	subleadL1EG_eta = leadL1EG_eta;

	leadL1EG_pt = pt;
	leadL1EG_eta = eta;
      }
      else if( pt > subleadL1EG_pt ){
	subleadL1EG_pt = pt;
	subleadL1EG_eta = eta;
      }
      numL1EG++;
    }

    h_l1EG25_1_pt->Fill(leadL1EG_pt,wgt);
    h_l1EG25_2_pt->Fill(subleadL1EG_pt,wgt);

    h_l1EG25_1_eta->Fill(leadL1EG_eta,wgt);
    h_l1EG25_2_eta->Fill(subleadL1EG_eta,wgt);

    h_numL1EG25->Fill(numL1EG,wgt);


    int numhltL1EG25=0;
    int numhltL1IsoEG22OrEG25=0;
    for( int iHLT=0; iHLT<int(eve->hltEle27WP85Gsf_pt_.size()); iHLT++ ){
      std::string name = eve->hltEle27WP85Gsf_filter_[iHLT];

      double pt  = eve->hltEle27WP85Gsf_pt_[iHLT];
      double eta = eve->hltEle27WP85Gsf_eta_[iHLT];
      double phi = eve->hltEle27WP85Gsf_phi_[iHLT];

      if( name=="hltL1sL1SingleEG25" ){
	numhltL1EG25++;

	h_hltL1EG25_pt->Fill(pt,wgt);
	h_hltL1EG25_eta->Fill(eta,wgt);
	h_hltL1EG25_phi->Fill(phi,wgt);
      }
      else if( name=="hltL1sL1SingleIsoEG22erOrSingleEG25" ){
	numhltL1IsoEG22OrEG25++;

	h_hltL1IsoEG22OrEG25_pt->Fill(pt,wgt);
	h_hltL1IsoEG22OrEG25_eta->Fill(eta,wgt);
	h_hltL1IsoEG22OrEG25_phi->Fill(phi,wgt);
      }
    }

    h_numhltL1EG25->Fill(numhltL1EG25,wgt);
    h_numhltL1IsoEG22OrEG25->Fill(numhltL1IsoEG22OrEG25,wgt);

    h_L1HTT->Fill(eve->L1HTT_,wgt);

    h_HT30->Fill(HT30,wgt);



    bool passHLTEle27HT200 = false;
    if( insample < 0 ) passHLTEle27HT200 = ( eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_==1 );
    else               passHLTEle27HT200 = ( eve->pass_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v_==1 );


    h_HT30_v2->Fill(HT30,wgt);
    double minL1HTT = ( insample < 0 ) ? 100. : 125.;
    if( eve->L1HTT_ > minL1HTT ){
      h_HT30_v2_L1HTT->Fill(HT30,wgt);
      if( passHLTEle27HT200 ) h_HT30_v2_L1HTT_HLTHT->Fill(HT30,wgt);
    }

    if( !(eve->L1HTT_ > 125.) ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // L1 HTT > 125

    h_HT30_L1HTT125->Fill(HT30,wgt);

    if( !(passHLTEle27HT200) ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // Ele + HT pass

    h_HT30_L1HTT125_passHLTEle27HT200->Fill(HT30,wgt);

  } // end loop over events

  std::cout << "**************************************************************" << std::endl;
  std::cout << "\t Number of raw events = " << numEvents_all << std::endl;
  std::cout << "\t Number of gen weighted events = " << numEvents_wgt_gen << std::endl;
  std::cout << "\t Number of lumi weighted events = " << numEvents_wgt_lumi << std::endl;
  std::cout << "\t Number of gen * lumi weighted events = " << numEvents_wgt << std::endl;
  std::cout << "**************************************************************" << std::endl;



  //TEfficiency* h_eff_probe_pt  = new TEfficiency(*h_probe_pt_pass,*h_probe_pt_all);
  //TEfficiency* h_eff_probe_eta = new TEfficiency(*h_probe_eta_pass,*h_probe_eta_all);
  //TEfficiency* h_eff_probe_phi = new TEfficiency(*h_probe_phi_pass,*h_probe_phi_all);

  //h_eff_probe_eta->Write("h_eff_probe_eta");
  //h_eff_probe_phi->Write("h_eff_probe_phi");

  std::cout << " Done! " << std::endl;

  histofile.Write();
  histofile.Close();

}


double reweightPU( int nPU ){

  double PUscale[50];

  PUscale[0] = 6.83938;
  PUscale[1] = 2.30668;
  PUscale[2] = 1.91805;
  PUscale[3] = 1.83689;
  PUscale[4] = 1.89282;
  PUscale[5] = 1.94686;
  PUscale[6] = 1.97328;
  PUscale[7] = 1.93306;
  PUscale[8] = 1.82852;
  PUscale[9] = 1.67631;
  PUscale[10] = 1.48487;
  PUscale[11] = 1.27255;
  PUscale[12] = 1.0549;
  PUscale[13] = 0.845621;
  PUscale[14] = 0.665044;
  PUscale[15] = 0.508729;
  PUscale[16] = 0.389761;
  PUscale[17] = 0.289674;
  PUscale[18] = 0.217065;
  PUscale[19] = 0.161309;
  PUscale[20] = 0.118954;
  PUscale[21] = 0.0885512;
  PUscale[22] = 0.0675094;
  PUscale[23] = 0.0491225;
  PUscale[24] = 0.0386613;
  PUscale[25] = 0.0296305;
  PUscale[26] = 0.0252724;
  PUscale[27] = 0.0175392;
  PUscale[28] = 0.0135296;
  PUscale[29] = 0.011069;
  PUscale[30] = 0.00914227;
  PUscale[31] = 0.00730368;
  PUscale[32] = 0.00265843;
  PUscale[33] = 0.0021888;
  PUscale[34] = 0;
  PUscale[35] = 0;
  PUscale[36] = 0;
  PUscale[37] = 0;
  PUscale[38] = 0;
  PUscale[39] = 0;
  PUscale[40] = 0;
  PUscale[41] = 0;
  PUscale[42] = 0;
  PUscale[43] = 0;
  PUscale[44] = 0;
  PUscale[45] = 0;
  PUscale[46] = 0;
  PUscale[47] = 0;
  PUscale[48] = 0;
  PUscale[49] = 0;


  if( nPU>49 ) nPU = 49;

  return PUscale[nPU];
}

float DeltaR(float eta1,float phi1,float eta2,float phi2){
  float deltaPhi = TMath::Abs(phi1-phi2);
  float deltaEta = eta1-eta2;
  if(deltaPhi > TMath::Pi())
    deltaPhi = TMath::TwoPi() - deltaPhi;
  return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}
/*

To do:
1. HT efficiency of our path
2. Process more MC
3. Look at W and top control regions


x. Fix b-tag in treeMaker
x. CMS preliminary
x. Amount of data

TEfficiency* eff_pt_hlt = new TEfficiency(*h_probe_pt_pHLTall,*h_probe_pt_all);


root -l HistoFiles/hltEleHT_treeReader_SingleElectron_Run2015B_PromptReco_251244_251562_histo.root


TH1D* h_probe_pt_pL1T = (TH1D*)_file0->Get("h_probe_pt_pL1T");
TH1D* h_probe_pt_pL1TOR = (TH1D*)_file0->Get("h_probe_pt_pL1TOR");

TH1D* h_probe_eta_pL1T = (TH1D*)_file0->Get("h_probe_eta_pL1T");
TH1D* h_probe_eta_pL1TOR = (TH1D*)_file0->Get("h_probe_eta_pL1TOR");

TH1D* h_probe_phi_pL1T = (TH1D*)_file0->Get("h_probe_phi_pL1T");
TH1D* h_probe_phi_pL1TOR = (TH1D*)_file0->Get("h_probe_phi_pL1TOR");



TH1D* h_pv_all = new TH1D("h_pv_all",";Number of PVs", 50, 0-0.5, 50-0.5 );
chain->Draw("eve.numPVs_","eve.pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_==1");


h_numPVs->SetMarkerColor(kBlack);
h_numPVs_251244->SetMarkerColor(kRed);
h_numPVs_251251->SetMarkerColor(kBlue);
h_numPVs_251252->SetMarkerColor(kOrange);
h_numPVs_251561->SetMarkerColor(kMagenta);
h_numPVs_251562->SetMarkerColor(kGreen+1);

h_numPVs->SetMarkerStyle(20);
h_numPVs_251244->SetMarkerStyle(20);
h_numPVs_251251->SetMarkerStyle(20);
h_numPVs_251252->SetMarkerStyle(20);
h_numPVs_251561->SetMarkerStyle(20);
h_numPVs_251562->SetMarkerStyle(20);


TLegend *legend = new TLegend(0.6,0.53,0.89,0.89);

legend->SetFillColor(kWhite);
legend->SetLineColor(kWhite);
legend->SetShadowColor(kWhite);
legend->SetTextFont(42);
legend->SetTextSize(0.035);

legend->AddEntry(h_numPVs,"All runs","p");
legend->AddEntry(h_numPVs_251244,"251244","p");
legend->AddEntry(h_numPVs_251251,"251251","p");
legend->AddEntry(h_numPVs_251252,"251252","p");
legend->AddEntry(h_numPVs_251561,"251561","p");
legend->AddEntry(h_numPVs_251562,"251562","p");

h_numPVs->SetStats(0);

h_numPVs->Draw();
h_numPVs_251244->Draw("same");
h_numPVs_251251->Draw("same");
h_numPVs_251252->Draw("same");
h_numPVs_251561->Draw("same");
h_numPVs_251562->Draw("same");

legend->Draw();

c1->Print("numPVs_by_run.png");


h_numPVs->Scale( 1 / h_numPVs->Integral() );
h_numPVs_251244->Scale( 1 / h_numPVs_251244->Integral() );
h_numPVs_251251->Scale( 1 / h_numPVs_251251->Integral() );
h_numPVs_251252->Scale( 1 / h_numPVs_251252->Integral() );
h_numPVs_251561->Scale( 1 / h_numPVs_251561->Integral() );
h_numPVs_251562->Scale( 1 / h_numPVs_251562->Integral() );


h_numPVs->SetMaximum(0.11);

h_numPVs->Draw();
h_numPVs_251244->Draw("same");
h_numPVs_251251->Draw("same");
h_numPVs_251252->Draw("same");
h_numPVs_251561->Draw("same");
h_numPVs_251562->Draw("same");

legend->Draw();


c1->Print("numPVs_by_run_norm.png");










root -l HistoFiles/hltEleHT_treeReader_SingleElectron_Run2015B_PromptReco_251244_251562_histo.root HistoFiles/hltEleHT_treeReader_DYJets_M50_13TeV_Spring15_StartupFlat10to50bx50Raw_histo.root HistoFiles/hltEleHT_treeReader_DYJets_M50_13TeV_Spring15_StartupFlat10to50bx50Raw_OriginalRCTcalib_histo.root

TH1D* h_data_L1HTT = (TH1D*)_file0->Get("h_L1HTT");
TH1D* h_mc_L1HTT = (TH1D*)_file1->Get("h_L1HTT");
TH1D* h_mc_L1HTT_newRCTcalib = (TH1D*)_file1->Get("h_L1HTT_newRCTcalib");
TH1D* h_mc_L1HTT_oldRCTcalib = (TH1D*)_file2->Get("h_L1HTT_newRCTcalib");


h_data_L1HTT->Rebin(4);
h_mc_L1HTT->Rebin(4);
h_mc_L1HTT_newRCTcalib->Rebin(4);
h_mc_L1HTT_oldRCTcalib->Rebin(4);

h_mc_L1HTT->Scale(9.583464e-01);
h_mc_L1HTT_newRCTcalib->Scale(9.583464e-01);
h_mc_L1HTT_oldRCTcalib->Scale(9.583464e-01);


h_data_L1HTT->SetMarkerColor(kBlack);
h_mc_L1HTT->SetMarkerColor(kRed);
h_mc_L1HTT_newRCTcalib->SetMarkerColor(kGreen+1);
h_mc_L1HTT_oldRCTcalib->SetMarkerColor(kBlue);

h_data_L1HTT->SetLineColor(kBlack);
h_mc_L1HTT->SetLineColor(kRed);
h_mc_L1HTT_newRCTcalib->SetLineColor(kGreen+1);
h_mc_L1HTT_oldRCTcalib->SetLineColor(kBlue);

h_data_L1HTT->SetLineWidth(2);
h_mc_L1HTT->SetLineWidth(2);
h_mc_L1HTT_newRCTcalib->SetLineWidth(2);
h_mc_L1HTT_oldRCTcalib->SetLineWidth(2);

h_data_L1HTT->SetMarkerStyle(20);
h_mc_L1HTT->SetMarkerStyle(20);
h_mc_L1HTT_newRCTcalib->SetMarkerStyle(20);
h_mc_L1HTT_oldRCTcalib->SetMarkerStyle(20);

TLegend *legend = new TLegend(0.05,0.92,0.95,0.99);

legend->SetFillColor(kWhite);
legend->SetLineColor(kWhite);
legend->SetShadowColor(kWhite);
legend->SetTextFont(42);
legend->SetTextSize(0.032);

legend->SetNColumns(4);

legend->AddEntry(h_data_L1HTT,"Data","lp");
legend->AddEntry(h_mc_L1HTT,"MC standard","l");
legend->AddEntry(h_mc_L1HTT_newRCTcalib,"MC new RCT calib","l");
legend->AddEntry(h_mc_L1HTT_oldRCTcalib,"MC old RCT calib","l");

h_data_L1HTT->SetStats(0);
h_data_L1HTT->GetYaxis()->SetTitle("Number of events");
h_data_L1HTT->GetYaxis()->SetTitleOffset(1.2);


TCanvas* c1 = new TCanvas("c1","",700,600);

h_data_L1HTT->Draw("pe1");
h_mc_L1HTT->Draw("histsame");
h_mc_L1HTT_newRCTcalib->Draw("histsame");
h_mc_L1HTT_oldRCTcalib->Draw("histsame");

legend->Draw();


c1->Print("L1HTT_RCTcalib_data_mc_Zee_lin.png");
c1->Print("L1HTT_RCTcalib_data_mc_Zee_lin.pdf");

c1->SetLogy(1);

c1->Print("L1HTT_RCTcalib_data_mc_Zee_log.png");
c1->Print("L1HTT_RCTcalib_data_mc_Zee_log.pdf");



TEfficiency* h_eff_HT_L1_HLT = new TEfficiency(*h_HT30_L1HTT125_passHLTEle27HT200,*h_HT30);
TEfficiency* h_eff_HT_HLT = new TEfficiency(*h_HT30_L1HTT125_passHLTEle27HT200,*h_HT30_L1HTT125);



root -l HistoFiles/backup_v2/hltEleHT_treeReader_Zcr_DYJets_M50_13TeV_Spring15_Asympt25ns_histo.root HistoFiles/hltEleHT_treeReader_Zcr_DYJets_M50_13TeV_Spring15_Asympt25ns_histo.root

root -l HistoFiles/backup_v2/hltEleHT_treeReader_Zcr_SingleElectron_Run2015D_PromptReco_254231_257599_histo.root HistoFiles/hltEleHT_treeReader_SingleElectron_Run2015CD_PromptReco_254231_256869_histo.root

TH1D* h_es_new = (TH1D*)_file0->Get("h_event_selection")->Clone("h_es_new");
TH1D* h_es_old = (TH1D*)_file1->Get("h_event_selection")->Clone("h_es_old");

TH1D* h_ratio_new_over_old = (TH1D*)h_es_new->Clone("h_ratio_new_over_old");
h_ratio_new_over_old->Divide(h_es_old);

h_ratio_new_over_old->Draw("pe1");



TEfficiency* eff_numPVs_hlt = new TEfficiency(*h_probe_numPVs_pHLT,*h_probe_numPVs_all);

TEfficiency* eff_numPVs_hlt_l1t = new TEfficiency(*h_probe_numPVs_pHLT,*h_probe_numPVs_pL1T);

TEfficiency* eff_numPVs_l1t = new TEfficiency(*h_probe_numPVs_pL1T,*h_probe_numPVs_all);


h_HT30_L1HTT125_passHLTEle27HT200->Rebin(10);
h_HT30_L1HTT125->Rebin(10);
h_HT30->Rebin(10);


TEfficiency* eff_ht_hlt_all = new TEfficiency(*h_HT30_L1HTT125_passHLTEle27HT200,*h_HT30);
TEfficiency* eff_ht_hlt_l1t = new TEfficiency(*h_HT30_L1HTT125_passHLTEle27HT200,*h_HT30_L1HTT125);
TEfficiency* eff_ht_l1t_all = new TEfficiency(*h_HT30_L1HTT125,*h_HT30);

eff_ht_hlt_all->SetLineColor(kBlue);
eff_ht_hlt_l1t->SetLineColor(kRed);
eff_ht_l1t_all->SetLineColor(kGreen);

eff_ht_hlt_all->Draw();
eff_ht_hlt_l1t->Draw("same");
eff_ht_l1t_all->Draw("same");



*/
