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

void hltEleHT_treeReader_TTcr( int insample=1, int maxNentries=-1, int Njobs=1, int jobN=1, double intLumi=-1, double iSys_=0 ) {

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
    mySample_nGen_ = 45420978;//25357774;//25446993;
    mySample_sampleName_ = "TT_13TeV_Spring15_Asympt25ns";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_74X_mcRun2_asymptotic_v2_ext3-v1_triggerTree_v1/151009_135600/0000/");
    //mySample_inputDir_ = "/uscms_data/d2/dpuigh/TTH/miniAOD/CMSSW_7_2_3/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/");
  }
  else if( insample==2325 ){
    mySample_xSec_ = 3*2008.4; 
    mySample_nGen_ = 19310834;//59000;//28825132;
    mySample_sampleName_ = "DYJets_M50_13TeV_Spring15_Asympt25ns";
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
    mySample_xSec_ = 0.5085 * 0.577;// YR3 * BR(all)  
    mySample_nGen_ = 3933404;//199000;
    mySample_sampleName_ = "ttHTobb_M125_13TeV_powheg_pythia8_Spring15_Asympt25ns";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ttHTobb_M125_13TeV_powheg_pythia8/RunIISpring15DR74-Asympt25ns_74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151009_135215/0000/");
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

  std::string histofilename = Form("HistoFiles/hltEleHT_treeReader_TTcr_%s_%s", mySample_sampleName_.c_str(), s_end.c_str());

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

  double window = 10.;
  double pdgZmass = 91.1876;
  double MinMass  = pdgZmass - window;
  double MaxMass  = pdgZmass + window;

  TH1::SetDefaultSumw2();

  TH1D* h_numEvents = new TH1D("h_numEvents",";Number of events", 2, 0, 2 );
  TH1D* h_numEvents_wgt = new TH1D("h_numEvents_wgt",";Predicted number of events", 2, 0, 2 );

  TH1D* h_numPVs = new TH1D("h_numPVs",";Number of PVs", 50, 0-0.5, 50-0.5 );
  TH1D* h_numPVs_PUwgt = new TH1D("h_numPVs_PUwgt",";Number of PVs", 50, 0-0.5, 50-0.5 );


  int NumCuts = 9;
  TH1D* h_event_selection  = new TH1D("h_event_selection",";cut", NumCuts, 0, NumCuts );
  h_event_selection->GetXaxis()->SetBinLabel(1,"All");
  h_event_selection->GetXaxis()->SetBinLabel(2,"HLT Ele27WP");
  h_event_selection->GetXaxis()->SetBinLabel(3,">=1 electron");
  h_event_selection->GetXaxis()->SetBinLabel(4,"==1 electron");
  h_event_selection->GetXaxis()->SetBinLabel(5,">=4 jets");
  h_event_selection->GetXaxis()->SetBinLabel(6,">=2 b-jets");
  h_event_selection->GetXaxis()->SetBinLabel(7,"L1 HTT > 100");
  h_event_selection->GetXaxis()->SetBinLabel(8,"L1 HTT > 125");
  h_event_selection->GetXaxis()->SetBinLabel(9,"HLT HT200");


  TH1D* h_numEle = new TH1D("h_numEle",";number of electrons", 5, 0, 5 );

  TH1D* h_deltaR_jet_ele = new TH1D("h_deltaR_jet_ele",";#DeltaR(jet,ele)", 61, 0., 6.1 );
  TH1D* h_minDR_ele_reco_hlt = new TH1D("h_minDR_ele_reco_hlt",";#DeltaR(reco ele,hlt ele)", 122, 0., 6.1 );

  TH1D* h_met_pt = new TH1D("h_met_pt",";MET", 300, 0., 300. );

  TH1D* h_jet_pt = new TH1D("h_jet_pt",";jet p_{T}", 300, 0., 300. );
  TH1D* h_jet_1_pt = new TH1D("h_jet_1_pt",";jet 1 p_{T}", 300, 0., 300. );
  TH1D* h_jet_2_pt = new TH1D("h_jet_2_pt",";jet 2 p_{T}", 300, 0., 300. );
  TH1D* h_jet_3_pt = new TH1D("h_jet_3_pt",";jet 3 p_{T}", 300, 0., 300. );
  TH1D* h_jet_4_pt = new TH1D("h_jet_4_pt",";jet 4 p_{T}", 300, 0., 300. );

  TH1D* h_jet_eta = new TH1D("h_jet_eta",";jet #eta", 52, -2.5, 2.5 );
  TH1D* h_jet_phi = new TH1D("h_jet_phi",";jet #phi", 64, -3.2, 3.2 );

  TH1D* h_jet_csv = new TH1D("h_jet_csv",";jet CSV", 152, -0.5, 1.02 );

  TH1D* h_numJet = new TH1D("h_numJet",";Number of Jets", 8, 0-0.5, 8-0.5 );
  TH1D* h_numBtag = new TH1D("h_numBtag",";Number of b-tagged Jets", 5, 0-0.5, 5-0.5 );

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
  int Nptbins = 14;
  double ptbins[] = { 10, 15, 20, 30, 40, 50, 60, 70, 80, 100, 125, 150, 200, 300, 500 };
  double ptMax = 200.;

  int Netabins = 12;
  double etabins[] = { -2.5, -2.1, -1.65, -1.2, -0.9, -0.45, 0, 0.45, 0.9, 1.2, 1.65, 2.1, 2.5 };

  int Nphibins = 16;
  double phibins[] = { -3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2 };

  TH1D* h_electron_pt  = new TH1D("h_electron_pt", ";electron p_{T}", Nptbins, ptbins );
  TH1D* h_electron_eta = new TH1D("h_electron_eta",";electron #eta", Netabins, etabins );
  TH1D* h_electron_phi = new TH1D("h_electron_phi",";electron #phi", Nphibins, phibins );

  TH1D* h_hltEle_pt  = new TH1D("h_hltEle_pt", ";HLT electron p_{T}", Nptbins, ptbins );
  TH1D* h_hltEle_eta = new TH1D("h_hltEle_eta",";HLT electron #eta", Netabins, etabins );
  TH1D* h_hltEle_phi = new TH1D("h_hltEle_phi",";HLT electron #phi", Nphibins, phibins );


  double HTmax = 1000.;
  int numHTbins = 1000;
  double L1HTmax = HTmax;
  int numL1HTTbins = numHTbins;

  TH2D* h_HT30_HT30er = new TH2D("h_HT30_HT30er",";reco H_{T};reco H_{T}", numHTbins, 0, HTmax, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT = new TH2D("h_HT30_L1HTT",";reco H_{T};reco H_{T}", numHTbins, 0, HTmax, numL1HTTbins, 0, L1HTmax );

  TH1D* h_L1HTT = new TH1D("h_L1HTT",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_HT30 = new TH1D("h_HT30",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125 = new TH1D("h_HT30_L1HTT125",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_passHLTEle27HT200 = new TH1D("h_HT30_L1HTT125_passHLTEle27HT200",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100 = new TH1D("h_HT30_L1HTT100",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_passHLTEle27HT200 = new TH1D("h_HT30_L1HTT100_passHLTEle27HT200",";reco H_{T}", numHTbins, 0, HTmax );

  TH1D* h_HT30er = new TH1D("h_HT30er",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30er_L1HTT125 = new TH1D("h_HT30er_L1HTT125",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30er_L1HTT125_passHLTEle27HT200 = new TH1D("h_HT30er_L1HTT125_passHLTEle27HT200",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30er_L1HTT100 = new TH1D("h_HT30er_L1HTT100",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30er_L1HTT100_passHLTEle27HT200 = new TH1D("h_HT30er_L1HTT100_passHLTEle27HT200",";reco H_{T}", numHTbins, 0, HTmax );

  // reco elePt
  TH2D* h_L1HTT_elePt = new TH2D("h_L1HTT_elePt",";electron p_{T};L1 HTT", Nptbins, ptbins, numL1HTTbins, 0, L1HTmax );
  TH2D* h_HT30_elePt = new TH2D("h_HT30_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT125_elePt = new TH2D("h_HT30_L1HTT125_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT125_passHLTEle27HT200_elePt = new TH2D("h_HT30_L1HTT125_passHLTEle27HT200_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT100_elePt = new TH2D("h_HT30_L1HTT100_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT100_passHLTEle27HT200_elePt = new TH2D("h_HT30_L1HTT100_passHLTEle27HT200_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );

  TH2D* h_HT30_4j_elePt = new TH2D("h_HT30_4j_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_4j_L1HTT125_elePt = new TH2D("h_HT30_4j_L1HTT125_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_4j_L1HTT125_passHLTEle27HT200_elePt = new TH2D("h_HT30_4j_L1HTT125_passHLTEle27HT200_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_4j_L1HTT100_elePt = new TH2D("h_HT30_4j_L1HTT100_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_4j_L1HTT100_passHLTEle27HT200_elePt = new TH2D("h_HT30_4j_L1HTT100_passHLTEle27HT200_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );

  TH2D* h_HT30er_elePt = new TH2D("h_HT30er_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30er_L1HTT125_elePt = new TH2D("h_HT30er_L1HTT125_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30er_L1HTT125_passHLTEle27HT200_elePt = new TH2D("h_HT30er_L1HTT125_passHLTEle27HT200_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30er_L1HTT100_elePt = new TH2D("h_HT30er_L1HTT100_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30er_L1HTT100_passHLTEle27HT200_elePt = new TH2D("h_HT30er_L1HTT100_passHLTEle27HT200_elePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );

  TH2D* h_HT30_numJet = new TH2D("h_HT30_numJet",";electron p_{T};reco H_{T}", 8, 0-0.5, 8-0.5, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT125_numJet = new TH2D("h_HT30_L1HTT125_numJet",";electron p_{T};reco H_{T}", 8, 0-0.5, 8-0.5, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT125_passHLTEle27HT200_numJet = new TH2D("h_HT30_L1HTT125_passHLTEle27HT200_numJet",";electron p_{T};reco H_{T}", 8, 0-0.5, 8-0.5, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT100_numJet = new TH2D("h_HT30_L1HTT100_numJet",";electron p_{T};reco H_{T}", 8, 0-0.5, 8-0.5, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT100_passHLTEle27HT200_numJet = new TH2D("h_HT30_L1HTT100_passHLTEle27HT200_numJet",";electron p_{T};reco H_{T}", 8, 0-0.5, 8-0.5, numHTbins, 0, HTmax );

  // eb
  TH2D* h_L1HTT_eleEBPt = new TH2D("h_L1HTT_eleEBPt",";electron p_{T};L1 HTT", Nptbins, ptbins, numL1HTTbins, 0, L1HTmax );
  TH2D* h_HT30_eleEBPt = new TH2D("h_HT30_eleEBPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT125_eleEBPt = new TH2D("h_HT30_L1HTT125_eleEBPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT125_passHLTEle27HT200_eleEBPt = new TH2D("h_HT30_L1HTT125_passHLTEle27HT200_eleEBPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT100_eleEBPt = new TH2D("h_HT30_L1HTT100_eleEBPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT100_passHLTEle27HT200_eleEBPt = new TH2D("h_HT30_L1HTT100_passHLTEle27HT200_eleEBPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );

  TH2D* h_HT30er_eleEBPt = new TH2D("h_HT30er_eleEBPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30er_L1HTT125_eleEBPt = new TH2D("h_HT30er_L1HTT125_eleEBPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30er_L1HTT125_passHLTEle27HT200_eleEBPt = new TH2D("h_HT30er_L1HTT125_passHLTEle27HT200_eleEBPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30er_L1HTT100_eleEBPt = new TH2D("h_HT30er_L1HTT100_eleEBPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30er_L1HTT100_passHLTEle27HT200_eleEBPt = new TH2D("h_HT30er_L1HTT100_passHLTEle27HT200_eleEBPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );

  // ee
  TH2D* h_L1HTT_eleEEPt = new TH2D("h_L1HTT_eleEEPt",";electron p_{T};L1 HTT", Nptbins, ptbins, numL1HTTbins, 0, L1HTmax );
  TH2D* h_HT30_eleEEPt = new TH2D("h_HT30_eleEEPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT125_eleEEPt = new TH2D("h_HT30_L1HTT125_eleEEPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT125_passHLTEle27HT200_eleEEPt = new TH2D("h_HT30_L1HTT125_passHLTEle27HT200_eleEEPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT100_eleEEPt = new TH2D("h_HT30_L1HTT100_eleEEPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT100_passHLTEle27HT200_eleEEPt = new TH2D("h_HT30_L1HTT100_passHLTEle27HT200_eleEEPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );

  TH2D* h_HT30er_eleEEPt = new TH2D("h_HT30er_eleEEPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30er_L1HTT125_eleEEPt = new TH2D("h_HT30er_L1HTT125_eleEEPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30er_L1HTT125_passHLTEle27HT200_eleEEPt = new TH2D("h_HT30er_L1HTT125_passHLTEle27HT200_eleEEPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30er_L1HTT100_eleEEPt = new TH2D("h_HT30er_L1HTT100_eleEEPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30er_L1HTT100_passHLTEle27HT200_eleEEPt = new TH2D("h_HT30er_L1HTT100_passHLTEle27HT200_eleEEPt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );

  ////


  TH1D* h_L1HTT_elePt30to50 = new TH1D("h_L1HTT_elePt30to50",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_HT30_elePt30to50 = new TH1D("h_HT30_elePt30to50",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_elePt30to50 = new TH1D("h_HT30_L1HTT125_elePt30to50",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_passHLTEle27HT200_elePt30to50 = new TH1D("h_HT30_L1HTT125_passHLTEle27HT200_elePt30to50",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_elePt30to50 = new TH1D("h_HT30_L1HTT100_elePt30to50",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_passHLTEle27HT200_elePt30to50 = new TH1D("h_HT30_L1HTT100_passHLTEle27HT200_elePt30to50",";reco H_{T}", numHTbins, 0, HTmax );

  TH1D* h_L1HTT_elePt50to80 = new TH1D("h_L1HTT_elePt50to80",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_HT30_elePt50to80 = new TH1D("h_HT30_elePt50to80",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_elePt50to80 = new TH1D("h_HT30_L1HTT125_elePt50to80",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_passHLTEle27HT200_elePt50to80 = new TH1D("h_HT30_L1HTT125_passHLTEle27HT200_elePt50to80",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_elePt50to80 = new TH1D("h_HT30_L1HTT100_elePt50to80",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_passHLTEle27HT200_elePt50to80 = new TH1D("h_HT30_L1HTT100_passHLTEle27HT200_elePt50to80",";reco H_{T}", numHTbins, 0, HTmax );

  TH1D* h_L1HTT_elePt80toInf = new TH1D("h_L1HTT_elePt80toInf",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_HT30_elePt80toInf = new TH1D("h_HT30_elePt80toInf",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_elePt80toInf = new TH1D("h_HT30_L1HTT125_elePt80toInf",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_passHLTEle27HT200_elePt80toInf = new TH1D("h_HT30_L1HTT125_passHLTEle27HT200_elePt80toInf",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_elePt80toInf = new TH1D("h_HT30_L1HTT100_elePt80toInf",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_passHLTEle27HT200_elePt80toInf = new TH1D("h_HT30_L1HTT100_passHLTEle27HT200_elePt80toInf",";reco H_{T}", numHTbins, 0, HTmax );

  // hlt elePt
  TH2D* h_L1HTT_hltElePt = new TH2D("h_L1HTT_hltElePt",";electron p_{T};L1 HTT", Nptbins, ptbins, numL1HTTbins, 0, L1HTmax );
  TH2D* h_HT30_hltElePt = new TH2D("h_HT30_hltElePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT125_hltElePt = new TH2D("h_HT30_L1HTT125_hltElePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT125_passHLTEle27HT200_hltElePt = new TH2D("h_HT30_L1HTT125_passHLTEle27HT200_hltElePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT100_hltElePt = new TH2D("h_HT30_L1HTT100_hltElePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );
  TH2D* h_HT30_L1HTT100_passHLTEle27HT200_hltElePt = new TH2D("h_HT30_L1HTT100_passHLTEle27HT200_hltElePt",";electron p_{T};reco H_{T}", Nptbins, ptbins, numHTbins, 0, HTmax );

  TH1D* h_L1HTT_hltElePt30to50 = new TH1D("h_L1HTT_hltElePt30to50",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_HT30_hltElePt30to50 = new TH1D("h_HT30_hltElePt30to50",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_hltElePt30to50 = new TH1D("h_HT30_L1HTT125_hltElePt30to50",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_passHLTEle27HT200_hltElePt30to50 = new TH1D("h_HT30_L1HTT125_passHLTEle27HT200_hltElePt30to50",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_hltElePt30to50 = new TH1D("h_HT30_L1HTT100_hltElePt30to50",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_passHLTEle27HT200_hltElePt30to50 = new TH1D("h_HT30_L1HTT100_passHLTEle27HT200_hltElePt30to50",";reco H_{T}", numHTbins, 0, HTmax );

  TH1D* h_L1HTT_hltElePt50to80 = new TH1D("h_L1HTT_hltElePt50to80",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_HT30_hltElePt50to80 = new TH1D("h_HT30_hltElePt50to80",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_hltElePt50to80 = new TH1D("h_HT30_L1HTT125_hltElePt50to80",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_passHLTEle27HT200_hltElePt50to80 = new TH1D("h_HT30_L1HTT125_passHLTEle27HT200_hltElePt50to80",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_hltElePt50to80 = new TH1D("h_HT30_L1HTT100_hltElePt50to80",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_passHLTEle27HT200_hltElePt50to80 = new TH1D("h_HT30_L1HTT100_passHLTEle27HT200_hltElePt50to80",";reco H_{T}", numHTbins, 0, HTmax );

  TH1D* h_L1HTT_hltElePt80toInf = new TH1D("h_L1HTT_hltElePt80toInf",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_HT30_hltElePt80toInf = new TH1D("h_HT30_hltElePt80toInf",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_hltElePt80toInf = new TH1D("h_HT30_L1HTT125_hltElePt80toInf",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_passHLTEle27HT200_hltElePt80toInf = new TH1D("h_HT30_L1HTT125_passHLTEle27HT200_hltElePt80toInf",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_hltElePt80toInf = new TH1D("h_HT30_L1HTT100_hltElePt80toInf",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_passHLTEle27HT200_hltElePt80toInf = new TH1D("h_HT30_L1HTT100_passHLTEle27HT200_hltElePt80toInf",";reco H_{T}", numHTbins, 0, HTmax );


  // >=4 jets
  TH1D* h_L1HTT_4j = new TH1D("h_L1HTT_4j",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_HT30_4j = new TH1D("h_HT30_4j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_4j = new TH1D("h_HT30_L1HTT125_4j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_passHLTEle27HT200_4j = new TH1D("h_HT30_L1HTT125_passHLTEle27HT200_4j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_4j = new TH1D("h_HT30_L1HTT100_4j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_passHLTEle27HT200_4j = new TH1D("h_HT30_L1HTT100_passHLTEle27HT200_4j",";reco H_{T}", numHTbins, 0, HTmax );

  TH1D* h_L1HTT_eq4j = new TH1D("h_L1HTT_eq4j",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_HT30_eq4j = new TH1D("h_HT30_eq4j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_eq4j = new TH1D("h_HT30_L1HTT125_eq4j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_passHLTEle27HT200_eq4j = new TH1D("h_HT30_L1HTT125_passHLTEle27HT200_eq4j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_eq4j = new TH1D("h_HT30_L1HTT100_eq4j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_passHLTEle27HT200_eq4j = new TH1D("h_HT30_L1HTT100_passHLTEle27HT200_eq4j",";reco H_{T}", numHTbins, 0, HTmax );

  TH1D* h_L1HTT_eq5j = new TH1D("h_L1HTT_eq5j",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_HT30_eq5j = new TH1D("h_HT30_eq5j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_eq5j = new TH1D("h_HT30_L1HTT125_eq5j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_passHLTEle27HT200_eq5j = new TH1D("h_HT30_L1HTT125_passHLTEle27HT200_eq5j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_eq5j = new TH1D("h_HT30_L1HTT100_eq5j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_passHLTEle27HT200_eq5j = new TH1D("h_HT30_L1HTT100_passHLTEle27HT200_eq5j",";reco H_{T}", numHTbins, 0, HTmax );

  TH1D* h_L1HTT_ge6j = new TH1D("h_L1HTT_ge6j",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_HT30_ge6j = new TH1D("h_HT30_ge6j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_ge6j = new TH1D("h_HT30_L1HTT125_ge6j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_passHLTEle27HT200_ge6j = new TH1D("h_HT30_L1HTT125_passHLTEle27HT200_ge6j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_ge6j = new TH1D("h_HT30_L1HTT100_ge6j",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_passHLTEle27HT200_ge6j = new TH1D("h_HT30_L1HTT100_passHLTEle27HT200_ge6j",";reco H_{T}", numHTbins, 0, HTmax );

  // >=4 jets  >=2 b-tags
  TH1D* h_L1HTT_4j2t = new TH1D("h_L1HTT_4j2t",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_HT30_4j2t = new TH1D("h_HT30_4j2t",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_4j2t = new TH1D("h_HT30_L1HTT125_4j2t",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT125_passHLTEle27HT200_4j2t = new TH1D("h_HT30_L1HTT125_passHLTEle27HT200_4j2t",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_4j2t = new TH1D("h_HT30_L1HTT100_4j2t",";reco H_{T}", numHTbins, 0, HTmax );
  TH1D* h_HT30_L1HTT100_passHLTEle27HT200_4j2t = new TH1D("h_HT30_L1HTT100_passHLTEle27HT200_4j2t",";reco H_{T}", numHTbins, 0, HTmax );



  std::vector<std::string> cat_labels;
  cat_labels.push_back("incl4j2t");
  cat_labels.push_back("4j2t");
  cat_labels.push_back("5j2t");
  cat_labels.push_back("6j2t");
  cat_labels.push_back("4j3t");
  cat_labels.push_back("5j3t");
  cat_labels.push_back("6j3t");
  cat_labels.push_back("4j4t");
  cat_labels.push_back("5j4t");
  cat_labels.push_back("6j4t");

  int NumCat = int(cat_labels.size());

  TH1D* h_category_yield = new TH1D("h_category_yield", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_recoHT300 = new TH1D("h_category_yield_recoHT300", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_recoHT300_L1HTT125 = new TH1D("h_category_yield_recoHT300_L1HTT125", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_recoHT300_L1HTT125_passHLTEle27HT200 = new TH1D("h_category_yield_recoHT300_L1HTT125_passHLTEle27HT200", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_recoHT300_L1HTT100 = new TH1D("h_category_yield_recoHT300_L1HTT100", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_recoHT300_L1HTT100_passHLTEle27HT200 = new TH1D("h_category_yield_recoHT300_L1HTT100_passHLTEle27HT200", ";category", NumCat, 0, NumCat );

  TH1D* h_category_yield_recoHTer300 = new TH1D("h_category_yield_recoHTer300", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_recoHTer300_L1HTT125 = new TH1D("h_category_yield_recoHTer300_L1HTT125", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_recoHTer300_L1HTT125_passHLTEle27HT200 = new TH1D("h_category_yield_recoHTer300_L1HTT125_passHLTEle27HT200", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_recoHTer300_L1HTT100 = new TH1D("h_category_yield_recoHTer300_L1HTT100", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_recoHTer300_L1HTT100_passHLTEle27HT200 = new TH1D("h_category_yield_recoHTer300_L1HTT100_passHLTEle27HT200", ";category", NumCat, 0, NumCat );

  for( int iCat=0; iCat<NumCat; iCat++ ){
    h_category_yield->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());
    h_category_yield_recoHT300->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());
    h_category_yield_recoHT300_L1HTT125->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());
    h_category_yield_recoHT300_L1HTT125_passHLTEle27HT200->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());
    h_category_yield_recoHT300_L1HTT100->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());
    h_category_yield_recoHT300_L1HTT100_passHLTEle27HT200->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());

    h_category_yield_recoHTer300->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());
    h_category_yield_recoHTer300_L1HTT125->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());
    h_category_yield_recoHTer300_L1HTT125_passHLTEle27HT200->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());
    h_category_yield_recoHTer300_L1HTT100->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());
    h_category_yield_recoHTer300_L1HTT100_passHLTEle27HT200->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat].c_str());
  }

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


    int ncut=0;
    h_event_selection->Fill(0.5+ncut++, wgt); // all

    bool pass_trigger = false;
    if( insample<0 ) pass_trigger = ( eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_==1 );
    else             pass_trigger = ( eve->pass_HLT_Ele27_WP85_Gsf_v_==1 );
    
    if( !pass_trigger ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // trigger

    h_numPVs->Fill(numPVs,wgt);

    // Remove PU reweighting for now
    //scalePU = 1;

    wgt *= scalePU;

    h_numPVs_PUwgt->Fill(numPVs,wgt);

    //////

    double HT30=0;
    double HT30er=0;
    for( int iJet=0; iJet<int(eve->jet_nocc_pt_.size()); iJet++ ){
      double pt  = eve->jet_nocc_pt_[iJet];
      double eta = eve->jet_nocc_eta_[iJet];

      if( !(pt>30. && fabs(eta)<3.0) ) continue;
      HT30 += pt;

      if( !(pt>30. && fabs(eta)<2.4) ) continue;
      HT30er += pt;
    }

    double MET = eve->pfMET_pt_;

    ////////////

    std::vector<int> ind_ele;
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

      if( pt > 30 && abs(eta)<2.1 && matchHLT ) ind_ele.push_back(iLep);
    }


    int numEle = int(ind_ele.size());

    h_numEle->Fill(numEle,wgt);


    if( !(numEle >= 1) ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // >= 1 electron

    if( !(numEle == 1) ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // == 1 electron


    int ind = ind_ele[0];

    TLorentzVector myEle;
    myEle.SetPtEtaPhiE( eve->lepton_pt_[ind], eve->lepton_eta_[ind], eve->lepton_phi_[ind], eve->lepton_energy_[ind] );

    double elePt = myEle.Pt();
    if( elePt>ptMax ) elePt = ptMax-0.0001;
    double eleEta = myEle.Eta();
    double elePhi = myEle.Phi();

    bool eleEB = ( eve->lepton_isEB_[ind]==1 );


    h_electron_pt->Fill(elePt, wgt);
    h_electron_eta->Fill(eleEta, wgt);
    h_electron_phi->Fill(elePhi, wgt);


    double hltElePt = -99;
    double hltEleEta = -99;
    double hltElePhi = -99;

    double minDR = 999;
    for( int iHLT=0; iHLT<int(eve->hltEle27WP85Gsf_pt_.size()); iHLT++ ){
      std::string name = eve->hltEle27WP85Gsf_filter_[iHLT];

      //if( name!="hltL1EG25Ele27WP85GsfTrackIsoFilter" ) continue;

      if( !((name=="hltL1EG25Ele27WP85GsfTrackIsoFilter" && insample>=0) ||
	    (name=="hltEle27WPLooseGsfTrackIsoFilter"    && insample< 0)) ) continue;

      double hltpt  = eve->hltEle27WP85Gsf_pt_[iHLT];
      double hlteta = eve->hltEle27WP85Gsf_eta_[iHLT];
      double hltphi = eve->hltEle27WP85Gsf_phi_[iHLT];

      double deltaR = DeltaR( eleEta, elePhi, hlteta, hltphi );

      if( deltaR < minDR ){
	minDR = deltaR;
	hltElePt = hltpt;
	hltEleEta = hlteta;
	hltElePhi = hltphi;
      }
    }

    if( minDR > 0.3 ) hltElePt = -99;

    h_minDR_ele_reco_hlt->Fill(minDR, wgt);



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
      for( int iEle=0; iEle<int(ind_ele.size()); iEle++ ){
	int ind = ind_ele[iEle];

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

      if( ((name=="hltL1EG25Ele27WP85GsfTrackIsoFilter" && insample>=0) ||
	   (name=="hltEle27WPLooseGsfTrackIsoFilter"    && insample< 0)) ){
	if( pt > hltElePt ){
	  hltElePt = pt;
	  hltEleEta = eta;
	  hltElePhi = phi;
	}
      }

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



    h_hltEle_pt->Fill(hltElePt, wgt);
    h_hltEle_eta->Fill(hltEleEta, wgt);
    h_hltEle_phi->Fill(hltElePhi, wgt);




    h_numhltL1EG25->Fill(numhltL1EG25,wgt);
    h_numhltL1IsoEG22OrEG25->Fill(numhltL1IsoEG22OrEG25,wgt);

    bool passHLTEle27HT200 = false;
    if( insample < 0 ) passHLTEle27HT200 = ( eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_==1 );
    else               passHLTEle27HT200 = ( eve->pass_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v_==1 );


    h_HT30_HT30er->Fill(HT30,HT30er,wgt);
    h_HT30_L1HTT->Fill(HT30,eve->L1HTT_,wgt);

    h_L1HTT->Fill(eve->L1HTT_,wgt);
    h_HT30->Fill(HT30,wgt);
    h_HT30er->Fill(HT30er,wgt);
    if( (eve->L1HTT_ > 125.) ){
      h_HT30_L1HTT125->Fill(HT30,wgt);
      if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200->Fill(HT30,wgt);

      h_HT30er_L1HTT125->Fill(HT30er,wgt);
      if( passHLTEle27HT200 ) h_HT30er_L1HTT125_passHLTEle27HT200->Fill(HT30er,wgt);
    }
    if( (eve->L1HTT_ > 100.) ){
      h_HT30_L1HTT100->Fill(HT30,wgt);
      if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200->Fill(HT30,wgt);

      h_HT30er_L1HTT100->Fill(HT30er,wgt);
      if( passHLTEle27HT200 ) h_HT30er_L1HTT100_passHLTEle27HT200->Fill(HT30er,wgt);
    }


    // 2D
    h_L1HTT_elePt->Fill(elePt,eve->L1HTT_,wgt);
    h_HT30_elePt->Fill(elePt,HT30,wgt);
    h_HT30er_elePt->Fill(elePt,HT30er,wgt);
    if( (eve->L1HTT_ > 125.) ){
      h_HT30_L1HTT125_elePt->Fill(elePt,HT30,wgt);
      if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_elePt->Fill(elePt,HT30,wgt);

      h_HT30er_L1HTT125_elePt->Fill(elePt,HT30er,wgt);
      if( passHLTEle27HT200 ) h_HT30er_L1HTT125_passHLTEle27HT200_elePt->Fill(elePt,HT30er,wgt);
    }
    if( (eve->L1HTT_ > 100.) ){
      h_HT30_L1HTT100_elePt->Fill(elePt,HT30,wgt);
      if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_elePt->Fill(elePt,HT30,wgt);

      h_HT30er_L1HTT100_elePt->Fill(elePt,HT30er,wgt);
      if( passHLTEle27HT200 ) h_HT30er_L1HTT100_passHLTEle27HT200_elePt->Fill(elePt,HT30er,wgt);
    }

    h_HT30_numJet->Fill(numJet,HT30,wgt);
    if( (eve->L1HTT_ > 125.) ){
      h_HT30_L1HTT125_numJet->Fill(numJet,HT30,wgt);
      if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_numJet->Fill(numJet,HT30,wgt);
    }
    if( (eve->L1HTT_ > 100.) ){
      h_HT30_L1HTT100_numJet->Fill(numJet,HT30,wgt);
      if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_numJet->Fill(numJet,HT30,wgt);
    }

    if( numJet>=4 ){
      h_HT30_4j_elePt->Fill(elePt,HT30,wgt);
      if( (eve->L1HTT_ > 125.) ){
	h_HT30_4j_L1HTT125_elePt->Fill(elePt,HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_4j_L1HTT125_passHLTEle27HT200_elePt->Fill(elePt,HT30,wgt);
      }
      if( (eve->L1HTT_ > 100.) ){
	h_HT30_4j_L1HTT100_elePt->Fill(elePt,HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_4j_L1HTT100_passHLTEle27HT200_elePt->Fill(elePt,HT30,wgt);
      }
    }

    //// EB and EE
    if( eleEB ){
      h_L1HTT_eleEBPt->Fill(elePt,eve->L1HTT_,wgt);
      h_HT30_eleEBPt->Fill(elePt,HT30,wgt);
      h_HT30er_eleEBPt->Fill(elePt,HT30er,wgt);
      if( (eve->L1HTT_ > 125.) ){
	h_HT30_L1HTT125_eleEBPt->Fill(elePt,HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_eleEBPt->Fill(elePt,HT30,wgt);

	h_HT30er_L1HTT125_eleEBPt->Fill(elePt,HT30er,wgt);
	if( passHLTEle27HT200 ) h_HT30er_L1HTT125_passHLTEle27HT200_eleEBPt->Fill(elePt,HT30er,wgt);
      }
      if( (eve->L1HTT_ > 100.) ){
	h_HT30_L1HTT100_eleEBPt->Fill(elePt,HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_eleEBPt->Fill(elePt,HT30,wgt);

	h_HT30er_L1HTT100_eleEBPt->Fill(elePt,HT30er,wgt);
	if( passHLTEle27HT200 ) h_HT30er_L1HTT100_passHLTEle27HT200_eleEBPt->Fill(elePt,HT30er,wgt);
      }
    }
    else {
      h_L1HTT_eleEEPt->Fill(elePt,eve->L1HTT_,wgt);
      h_HT30_eleEEPt->Fill(elePt,HT30,wgt);
      h_HT30er_eleEEPt->Fill(elePt,HT30er,wgt);
      if( (eve->L1HTT_ > 125.) ){
	h_HT30_L1HTT125_eleEEPt->Fill(elePt,HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_eleEEPt->Fill(elePt,HT30,wgt);

	h_HT30er_L1HTT125_eleEEPt->Fill(elePt,HT30er,wgt);
	if( passHLTEle27HT200 ) h_HT30er_L1HTT125_passHLTEle27HT200_eleEEPt->Fill(elePt,HT30er,wgt);
      }
      if( (eve->L1HTT_ > 100.) ){
	h_HT30_L1HTT100_eleEEPt->Fill(elePt,HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_eleEEPt->Fill(elePt,HT30,wgt);

	h_HT30er_L1HTT100_eleEEPt->Fill(elePt,HT30er,wgt);
	if( passHLTEle27HT200 ) h_HT30er_L1HTT100_passHLTEle27HT200_eleEEPt->Fill(elePt,HT30er,wgt);
      }
    }


    if( elePt>=30. && elePt<50. ){
      h_L1HTT_elePt30to50->Fill(eve->L1HTT_,wgt);
      h_HT30_elePt30to50->Fill(HT30,wgt);
      if( (eve->L1HTT_ > 125.) ){
	h_HT30_L1HTT125_elePt30to50->Fill(HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_elePt30to50->Fill(HT30,wgt);
      }
      if( (eve->L1HTT_ > 100.) ){
	h_HT30_L1HTT100_elePt30to50->Fill(HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_elePt30to50->Fill(HT30,wgt);
      }
    }
    else if( elePt>=50. && elePt<80. ){
      h_L1HTT_elePt50to80->Fill(eve->L1HTT_,wgt);
      h_HT30_elePt50to80->Fill(HT30,wgt);
      if( (eve->L1HTT_ > 125.) ){
	h_HT30_L1HTT125_elePt50to80->Fill(HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_elePt50to80->Fill(HT30,wgt);
      }
      if( (eve->L1HTT_ > 100.) ){
	h_HT30_L1HTT100_elePt50to80->Fill(HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_elePt50to80->Fill(HT30,wgt);
      }
    }
    else if( elePt>=80. ){
      h_L1HTT_elePt80toInf->Fill(eve->L1HTT_,wgt);
      h_HT30_elePt80toInf->Fill(HT30,wgt);
      if( (eve->L1HTT_ > 125.) ){
	h_HT30_L1HTT125_elePt80toInf->Fill(HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_elePt80toInf->Fill(HT30,wgt);
      }
      if( (eve->L1HTT_ > 100.) ){
	h_HT30_L1HTT100_elePt80toInf->Fill(HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_elePt80toInf->Fill(HT30,wgt);
      }
    }


    // hlt Ele
    if( hltElePt>0 ){
      // 2D
      h_L1HTT_hltElePt->Fill(hltElePt,eve->L1HTT_,wgt);
      h_HT30_hltElePt->Fill(hltElePt,HT30,wgt);
      if( (eve->L1HTT_ > 125.) ){
	h_HT30_L1HTT125_hltElePt->Fill(hltElePt,HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_hltElePt->Fill(hltElePt,HT30,wgt);
      }
      if( (eve->L1HTT_ > 100.) ){
	h_HT30_L1HTT100_hltElePt->Fill(hltElePt,HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_hltElePt->Fill(hltElePt,HT30,wgt);
      }

      if( hltElePt>=30. && hltElePt<50. ){
	h_L1HTT_hltElePt30to50->Fill(eve->L1HTT_,wgt);
	h_HT30_hltElePt30to50->Fill(HT30,wgt);
	if( (eve->L1HTT_ > 125.) ){
	  h_HT30_L1HTT125_hltElePt30to50->Fill(HT30,wgt);
	  if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_hltElePt30to50->Fill(HT30,wgt);
	}
	if( (eve->L1HTT_ > 100.) ){
	  h_HT30_L1HTT100_hltElePt30to50->Fill(HT30,wgt);
	  if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_hltElePt30to50->Fill(HT30,wgt);
	}
      }
      else if( hltElePt>=50. && hltElePt<80. ){
	h_L1HTT_hltElePt50to80->Fill(eve->L1HTT_,wgt);
	h_HT30_hltElePt50to80->Fill(HT30,wgt);
	if( (eve->L1HTT_ > 125.) ){
	  h_HT30_L1HTT125_hltElePt50to80->Fill(HT30,wgt);
	  if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_hltElePt50to80->Fill(HT30,wgt);
	}
	if( (eve->L1HTT_ > 100.) ){
	  h_HT30_L1HTT100_hltElePt50to80->Fill(HT30,wgt);
	  if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_hltElePt50to80->Fill(HT30,wgt);
	}
      }
      else if( hltElePt>=80. ){
	h_L1HTT_hltElePt80toInf->Fill(eve->L1HTT_,wgt);
	h_HT30_hltElePt80toInf->Fill(HT30,wgt);
	if( (eve->L1HTT_ > 125.) ){
	  h_HT30_L1HTT125_hltElePt80toInf->Fill(HT30,wgt);
	  if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_hltElePt80toInf->Fill(HT30,wgt);
	}
	if( (eve->L1HTT_ > 100.) ){
	  h_HT30_L1HTT100_hltElePt80toInf->Fill(HT30,wgt);
	  if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_hltElePt80toInf->Fill(HT30,wgt);
	}
      }

    }




    if( !(numJet>=4) ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // >=4 jets


    h_L1HTT_4j->Fill(eve->L1HTT_,wgt);
    h_HT30_4j->Fill(HT30,wgt);
    if( (eve->L1HTT_ > 125.) ){
      h_HT30_L1HTT125_4j->Fill(HT30,wgt);
      if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_4j->Fill(HT30,wgt);
    }
    if( (eve->L1HTT_ > 100.) ){
      h_HT30_L1HTT100_4j->Fill(HT30,wgt);
      if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_4j->Fill(HT30,wgt);
    }

    if( numJet==4 ){
      h_L1HTT_eq4j->Fill(eve->L1HTT_,wgt);
      h_HT30_eq4j->Fill(HT30,wgt);
      if( (eve->L1HTT_ > 125.) ){
	h_HT30_L1HTT125_eq4j->Fill(HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_eq4j->Fill(HT30,wgt);
      }
      if( (eve->L1HTT_ > 100.) ){
	h_HT30_L1HTT100_eq4j->Fill(HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_eq4j->Fill(HT30,wgt);
      }
    }
    else if( numJet==5 ){
      h_L1HTT_eq5j->Fill(eve->L1HTT_,wgt);
      h_HT30_eq5j->Fill(HT30,wgt);
      if( (eve->L1HTT_ > 125.) ){
	h_HT30_L1HTT125_eq5j->Fill(HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_eq5j->Fill(HT30,wgt);
      }
      if( (eve->L1HTT_ > 100.) ){
	h_HT30_L1HTT100_eq5j->Fill(HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_eq5j->Fill(HT30,wgt);
      }
    }
    else if( numJet>=6 ){
      h_L1HTT_ge6j->Fill(eve->L1HTT_,wgt);
      h_HT30_ge6j->Fill(HT30,wgt);
      if( (eve->L1HTT_ > 125.) ){
	h_HT30_L1HTT125_ge6j->Fill(HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_ge6j->Fill(HT30,wgt);
      }
      if( (eve->L1HTT_ > 100.) ){
	h_HT30_L1HTT100_ge6j->Fill(HT30,wgt);
	if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_ge6j->Fill(HT30,wgt);
      }
    }


    if( !(numBtag>=2) ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // >=2 b-tags


    h_L1HTT_4j2t->Fill(eve->L1HTT_,wgt);
    h_HT30_4j2t->Fill(HT30,wgt);
    if( (eve->L1HTT_ > 125.) ){
      h_HT30_L1HTT125_4j2t->Fill(HT30,wgt);
      if( passHLTEle27HT200 ) h_HT30_L1HTT125_passHLTEle27HT200_4j2t->Fill(HT30,wgt);
    }
    if( (eve->L1HTT_ > 100.) ){
      h_HT30_L1HTT100_4j2t->Fill(HT30,wgt);
      if( passHLTEle27HT200 ) h_HT30_L1HTT100_passHLTEle27HT200_4j2t->Fill(HT30,wgt);
    }



    int this_category = -1;
    if( numJet==4 && numBtag==2) this_category=1;
    if( numJet==5 && numBtag==2) this_category=2;
    if( numJet>=6 && numBtag==2) this_category=3;	
    if( numJet==4 && numBtag==3) this_category=4;
    if( numJet==5 && numBtag==3) this_category=5;
    if( numJet>=6 && numBtag==3) this_category=6;
    if( numJet==4 && numBtag>=4) this_category=7;
    if( numJet==5 && numBtag>=4) this_category=8;
    if( numJet>=6 && numBtag>=4) this_category=9;


    h_category_yield->Fill(0.5,wgt);
    h_category_yield->Fill(this_category,wgt);

    if( HT30>300 ){
      h_category_yield_recoHT300->Fill(0.5,wgt);
      h_category_yield_recoHT300->Fill(this_category,wgt);

      if( eve->L1HTT_ > 125. ){
	h_category_yield_recoHT300_L1HTT125->Fill(0.5,wgt);
	h_category_yield_recoHT300_L1HTT125->Fill(this_category,wgt);

	if( passHLTEle27HT200 ){
	  h_category_yield_recoHT300_L1HTT125_passHLTEle27HT200->Fill(0.5,wgt);
	  h_category_yield_recoHT300_L1HTT125_passHLTEle27HT200->Fill(this_category,wgt);
	}
      }
      if( eve->L1HTT_ > 100. ){
	h_category_yield_recoHT300_L1HTT100->Fill(0.5,wgt);
	h_category_yield_recoHT300_L1HTT100->Fill(this_category,wgt);

	if( passHLTEle27HT200 ){
	  h_category_yield_recoHT300_L1HTT100_passHLTEle27HT200->Fill(0.5,wgt);
	  h_category_yield_recoHT300_L1HTT100_passHLTEle27HT200->Fill(this_category,wgt);
	}
      }
    }

    if( HT30er>300 ){
      h_category_yield_recoHTer300->Fill(0.5,wgt);
      h_category_yield_recoHTer300->Fill(this_category,wgt);

      if( eve->L1HTT_ > 125. ){
	h_category_yield_recoHTer300_L1HTT125->Fill(0.5,wgt);
	h_category_yield_recoHTer300_L1HTT125->Fill(this_category,wgt);

	if( passHLTEle27HT200 ){
	  h_category_yield_recoHTer300_L1HTT125_passHLTEle27HT200->Fill(0.5,wgt);
	  h_category_yield_recoHTer300_L1HTT125_passHLTEle27HT200->Fill(this_category,wgt);
	}
      }
      if( eve->L1HTT_ > 100. ){
	h_category_yield_recoHTer300_L1HTT100->Fill(0.5,wgt);
	h_category_yield_recoHTer300_L1HTT100->Fill(this_category,wgt);

	if( passHLTEle27HT200 ){
	  h_category_yield_recoHTer300_L1HTT100_passHLTEle27HT200->Fill(0.5,wgt);
	  h_category_yield_recoHTer300_L1HTT100_passHLTEle27HT200->Fill(this_category,wgt);
	}
      }
    }


    if( !(eve->L1HTT_ > 100.) ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // Ele + HT pass

    if( !(eve->L1HTT_ > 125.) ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // Ele + HT pass

    if( !(passHLTEle27HT200) ) continue;
    h_event_selection->Fill(0.5+ncut++, wgt); // Ele + HT pass



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

  PUscale[0] = 7.37587;
  PUscale[1] = 4.64965;
  PUscale[2] = 3.60649;
  PUscale[3] = 3.1161;
  PUscale[4] = 2.92148;
  PUscale[5] = 2.71865;
  PUscale[6] = 2.51325;
  PUscale[7] = 2.26984;
  PUscale[8] = 1.98533;
  PUscale[9] = 1.69326;
  PUscale[10] = 1.4068;
  PUscale[11] = 1.14087;
  PUscale[12] = 0.895973;
  PUscale[13] = 0.68638;
  PUscale[14] = 0.519752;
  PUscale[15] = 0.385682;
  PUscale[16] = 0.284182;
  PUscale[17] = 0.20762;
  PUscale[18] = 0.151705;
  PUscale[19] = 0.111675;
  PUscale[20] = 0.0786714;
  PUscale[21] = 0.0580111;
  PUscale[22] = 0.0428043;
  PUscale[23] = 0.0308221;
  PUscale[24] = 0.0252111;
  PUscale[25] = 0.0165618;
  PUscale[26] = 0.0140687;
  PUscale[27] = 0.0102577;
  PUscale[28] = 0.00680251;
  PUscale[29] = 0.00851538;
  PUscale[30] = 0.00127876;
  PUscale[31] = 0;
  PUscale[32] = 0.00176625;
  PUscale[33] = 0;
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


h_HT30->Rebin(10);
h_HT30_L1HTT125->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200->Divide(h_HT30);

h_HT30_L1HTT125->Divide(h_HT30);

h_HT30_L1HTT125_passHLTEle27HT200->Draw();
h_HT30_L1HTT125->Draw();



h_HT30_4j->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_4j->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_4j->Divide(h_HT30_4j);
h_HT30_L1HTT125_passHLTEle27HT200_4j->SetMarkerColor(kMagenta+2);
h_HT30_L1HTT125_passHLTEle27HT200_4j->SetMarkerStyle(20);

h_HT30_4j2t->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_4j2t->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_4j2t->Divide(h_HT30_4j2t);
h_HT30_L1HTT125_passHLTEle27HT200_4j2t->SetMarkerColor(kOrange+2);
h_HT30_L1HTT125_passHLTEle27HT200_4j2t->SetMarkerStyle(20);

h_HT30->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200->Divide(h_HT30);
h_HT30_L1HTT125_passHLTEle27HT200->SetMarkerColor(kBlack);
h_HT30_L1HTT125_passHLTEle27HT200->SetMarkerStyle(20);

h_HT30_elePt30to50->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_elePt30to50->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_elePt30to50->Divide(h_HT30_elePt30to50);
h_HT30_L1HTT125_passHLTEle27HT200_elePt30to50->SetMarkerColor(kBlue);
h_HT30_L1HTT125_passHLTEle27HT200_elePt30to50->SetMarkerStyle(20);

h_HT30_elePt50to80->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_elePt50to80->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_elePt50to80->Divide(h_HT30_elePt50to80);
h_HT30_L1HTT125_passHLTEle27HT200_elePt50to80->SetMarkerColor(kRed);
h_HT30_L1HTT125_passHLTEle27HT200_elePt50to80->SetMarkerStyle(20);

h_HT30_elePt80toInf->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_elePt80toInf->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_elePt80toInf->Divide(h_HT30_elePt80toInf);
h_HT30_L1HTT125_passHLTEle27HT200_elePt80toInf->SetMarkerColor(kGreen+1);
h_HT30_L1HTT125_passHLTEle27HT200_elePt80toInf->SetMarkerStyle(20);


h_HT30_L1HTT125_passHLTEle27HT200->SetStats(0);
h_HT30_L1HTT125_passHLTEle27HT200->GetYaxis()->SetRangeUser(0.0,1.2);

h_HT30_L1HTT125_passHLTEle27HT200->Draw("pe1");
h_HT30_L1HTT125_passHLTEle27HT200_elePt30to50->Draw("pe1same");
h_HT30_L1HTT125_passHLTEle27HT200_elePt50to80->Draw("pe1same");
h_HT30_L1HTT125_passHLTEle27HT200_elePt80toInf->Draw("pe1same");







//////


h_HT30->Rebin(10);
h_HT30_L1HTT125->Rebin(10);
h_HT30_L1HTT125->Divide(h_HT30);
h_HT30_L1HTT125->SetMarkerColor(kBlack);
h_HT30_L1HTT125->SetMarkerStyle(20);

h_HT30_elePt30to50->Rebin(10);
h_HT30_L1HTT125_elePt30to50->Rebin(10);
h_HT30_L1HTT125_elePt30to50->Divide(h_HT30_elePt30to50);
h_HT30_L1HTT125_elePt30to50->SetMarkerColor(kBlue);
h_HT30_L1HTT125_elePt30to50->SetMarkerStyle(20);

h_HT30_elePt50to80->Rebin(10);
h_HT30_L1HTT125_elePt50to80->Rebin(10);
h_HT30_L1HTT125_elePt50to80->Divide(h_HT30_elePt50to80);
h_HT30_L1HTT125_elePt50to80->SetMarkerColor(kRed);
h_HT30_L1HTT125_elePt50to80->SetMarkerStyle(20);

h_HT30_elePt80toInf->Rebin(10);
h_HT30_L1HTT125_elePt80toInf->Rebin(10);
h_HT30_L1HTT125_elePt80toInf->Divide(h_HT30_elePt80toInf);
h_HT30_L1HTT125_elePt80toInf->SetMarkerColor(kGreen+1);
h_HT30_L1HTT125_elePt80toInf->SetMarkerStyle(20);


h_HT30_L1HTT125->SetStats(0);
h_HT30_L1HTT125->GetYaxis()->SetRangeUser(0.0,1.2);

h_HT30_L1HTT125->Draw("pe1");
h_HT30_L1HTT125_elePt30to50->Draw("pe1same");
h_HT30_L1HTT125_elePt50to80->Draw("pe1same");
h_HT30_L1HTT125_elePt80toInf->Draw("pe1same");





////// hlt ele


h_HT30->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200->Divide(h_HT30);
h_HT30_L1HTT125_passHLTEle27HT200->SetMarkerColor(kBlack);
h_HT30_L1HTT125_passHLTEle27HT200->SetMarkerStyle(20);

h_HT30_hltElePt30to50->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt30to50->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt30to50->Divide(h_HT30_hltElePt30to50);
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt30to50->SetMarkerColor(kBlue);
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt30to50->SetMarkerStyle(20);

h_HT30_hltElePt50to80->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt50to80->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt50to80->Divide(h_HT30_hltElePt50to80);
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt50to80->SetMarkerColor(kRed);
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt50to80->SetMarkerStyle(20);

h_HT30_hltElePt80toInf->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt80toInf->Rebin(10);
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt80toInf->Divide(h_HT30_hltElePt80toInf);
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt80toInf->SetMarkerColor(kGreen+1);
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt80toInf->SetMarkerStyle(20);


h_HT30_L1HTT125_passHLTEle27HT200->SetStats(0);
h_HT30_L1HTT125_passHLTEle27HT200->GetYaxis()->SetRangeUser(0.0,1.2);

h_HT30_L1HTT125_passHLTEle27HT200->Draw("pe1");
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt30to50->Draw("pe1same");
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt50to80->Draw("pe1same");
h_HT30_L1HTT125_passHLTEle27HT200_hltElePt80toInf->Draw("pe1same");



h_category_yield->SetMarkerColor(kBlack);
h_category_yield->SetMarkerStyle(20);

h_category_yield_recoHT300->SetMarkerColor(kBlue);
h_category_yield_recoHT300->SetMarkerStyle(20);

h_category_yield_recoHT300_L1HTT125->SetMarkerColor(kRed);
h_category_yield_recoHT300_L1HTT125->SetMarkerStyle(20);

h_category_yield_recoHT300_L1HTT125_passHLTEle27HT200->SetMarkerColor(kGreen+1);
h_category_yield_recoHT300_L1HTT125_passHLTEle27HT200->SetMarkerStyle(20);

h_category_yield->Draw("pe1");
h_category_yield_recoHT300->Draw("pe1same");
h_category_yield_recoHT300_L1HTT125->Draw("pe1same");
h_category_yield_recoHT300_L1HTT125_passHLTEle27HT200->Draw("pe1same");



h_HT30_elePt->RebinY(10);
h_HT30_L1HTT125_passHLTEle27HT200_elePt->RebinY(10);
h_HT30_L1HTT125_passHLTEle27HT200_elePt->Divide(h_HT30_elePt);
h_HT30_L1HTT125_passHLTEle27HT200_elePt->Draw("colz");
h_HT30_L1HTT125_passHLTEle27HT200_elePt->ProjectionY("h_HT_elePt30to40",4,4);


h_HT30er_elePt->RebinY(10);
h_HT30er_L1HTT125_passHLTEle27HT200_elePt->RebinY(10);
h_HT30er_L1HTT125_passHLTEle27HT200_elePt->Divide(h_HT30er_elePt);
h_HT30er_L1HTT125_passHLTEle27HT200_elePt->Draw("colz");
h_HT30er_L1HTT125_passHLTEle27HT200_elePt->ProjectionY("h_HTer_elePt30to40",4,4);


h_HT30_eleEBPt->RebinY(10);
h_HT30_L1HTT125_passHLTEle27HT200_eleEBPt->RebinY(10);
h_HT30_L1HTT125_passHLTEle27HT200_eleEBPt->Divide(h_HT30_eleEBPt);
h_HT30_L1HTT125_passHLTEle27HT200_eleEBPt->Draw("colz");
h_HT30_L1HTT125_passHLTEle27HT200_eleEBPt->ProjectionY("h_HT_eleEBPt30to40",4,4);

h_HT30_eleEEPt->RebinY(10);
h_HT30_L1HTT125_passHLTEle27HT200_eleEEPt->RebinY(10);
h_HT30_L1HTT125_passHLTEle27HT200_eleEEPt->Divide(h_HT30_eleEEPt);
h_HT30_L1HTT125_passHLTEle27HT200_eleEEPt->Draw("colz");
h_HT30_L1HTT125_passHLTEle27HT200_eleEEPt->ProjectionY("h_HT_eleEEPt30to40",4,4);


h_HT_elePt30to40->SetMarkerColor(kBlack);
h_HTer_elePt30to40->SetMarkerColor(kBlue);
h_HT_eleEBPt30to40->SetMarkerColor(kRed);
h_HT_eleEEPt30to40->SetMarkerColor(kGreen+1);

h_HT_elePt30to40->SetMarkerStyle(20);
h_HTer_elePt30to40->SetMarkerStyle(20);
h_HT_eleEBPt30to40->SetMarkerStyle(20);
h_HT_eleEEPt30to40->SetMarkerStyle(20);

h_HT_elePt30to40->Draw("pe1");
h_HTer_elePt30to40->Draw("pe1same");
h_HT_eleEBPt30to40->Draw("pe1same");
h_HT_eleEEPt30to40->Draw("pe1same");



root -b -q macros/head13TeV.C macros/hltEleHT_treeReader_ttbarCR.C+'(2500,-1,10,1,10000,0)' >&! out_hltEleHT_treeReader_ttbarCR_2500_1_2015_08_11_1400.log &
root -b -q macros/head13TeV.C macros/hltEleHT_treeReader_ttbarCR.C+'(2500,-1,10,2,10000,0)' >&! out_hltEleHT_treeReader_ttbarCR_2500_2_2015_08_11_1400.log &
root -b -q macros/head13TeV.C macros/hltEleHT_treeReader_ttbarCR.C+'(2500,-1,10,3,10000,0)' >&! out_hltEleHT_treeReader_ttbarCR_2500_3_2015_08_11_1400.log &
root -b -q macros/head13TeV.C macros/hltEleHT_treeReader_ttbarCR.C+'(2500,-1,10,4,10000,0)' >&! out_hltEleHT_treeReader_ttbarCR_2500_4_2015_08_11_1400.log &
root -b -q macros/head13TeV.C macros/hltEleHT_treeReader_ttbarCR.C+'(2500,-1,10,5,10000,0)' >&! out_hltEleHT_treeReader_ttbarCR_2500_5_2015_08_11_1400.log &
root -b -q macros/head13TeV.C macros/hltEleHT_treeReader_ttbarCR.C+'(2500,-1,10,6,10000,0)' >&! out_hltEleHT_treeReader_ttbarCR_2500_6_2015_08_11_1400.log &
root -b -q macros/head13TeV.C macros/hltEleHT_treeReader_ttbarCR.C+'(2500,-1,10,7,10000,0)' >&! out_hltEleHT_treeReader_ttbarCR_2500_7_2015_08_11_1400.log &
root -b -q macros/head13TeV.C macros/hltEleHT_treeReader_ttbarCR.C+'(2500,-1,10,8,10000,0)' >&! out_hltEleHT_treeReader_ttbarCR_2500_8_2015_08_11_1400.log &
root -b -q macros/head13TeV.C macros/hltEleHT_treeReader_ttbarCR.C+'(2500,-1,10,9,10000,0)' >&! out_hltEleHT_treeReader_ttbarCR_2500_9_2015_08_11_1400.log &
root -b -q macros/head13TeV.C macros/hltEleHT_treeReader_ttbarCR.C+'(2500,-1,10,10,10000,0)' >&! out_hltEleHT_treeReader_ttbarCR_2500_10_2015_08_11_1400.log &


hadd -f HistoFiles/hltEleHT_treeReader_ttbarCR_TTJets_powheg_Spring15_Asympt25ns_histo.root HistoFiles/hltEleHT_treeReader_ttbarCR_TTJets_powheg_Spring15_Asympt25ns_histo_*.root
rm HistoFiles/hltEleHT_treeReader_ttbarCR_TTJets_powheg_Spring15_Asympt25ns_histo_*.root



*/
