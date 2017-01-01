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
float DeltaPhi(float phi1,float phi2);

//*****************************************************************************

void dileptonCR_optimize_treeReader( int insample=1, int maxNentries=-1, int Njobs=1, int jobN=1, double intLumi=-1, int ttCat_=-1 ) {

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
    mySample_nGen_ = 115091972;//25357774;//25446993;
    mySample_sampleName_ = "ttbar";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2_triggerTree_v1/151015_183717/0000/");
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1_triggerTree_v1/151015_183747/0000/");
    if( ttCat_>=0 ){
      if( ttCat_==0 ) mySample_sampleName_ = "ttlf";
      if( ttCat_==1 ) mySample_sampleName_ = "ttcc";
      if( ttCat_==2 ) mySample_sampleName_ = "ttb";
      if( ttCat_==3 ) mySample_sampleName_ = "tt2b";
      if( ttCat_==4 ) mySample_sampleName_ = "ttbb";
    }
  }
  else if( insample==2300 ){
    mySample_xSec_ = 6025.2; 
    mySample_nGen_ = 19310834;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M50";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3_triggerTree_v1/151015_183818/0000/");
  }
  else if( insample==2310 ){
    mySample_xSec_ = 22635.09; 
    mySample_nGen_ = 22217467;//19925500;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M10to50";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_183844/0000/");
  }
  else if( insample==2400 ){
    mySample_xSec_ = 61526.7;  
    mySample_nGen_ = 16518218;
    mySample_sampleName_ = "WJets";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184003/0000/");
  }
  else if( insample==2524 ){
    mySample_xSec_ = 0.435;  
    mySample_nGen_ = 430330;
    mySample_sampleName_ = "ttW_had";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184416/0000/");
  }
  else if( insample==2525 ){
    mySample_xSec_ = 0.21;  
    mySample_nGen_ = 129850;
    mySample_sampleName_ = "ttW_lep";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184445/0000/");
  }
  else if( insample==2523 ){
    mySample_xSec_ = 0.611;  
    mySample_nGen_ = 351398;
    mySample_sampleName_ = "ttZ_had";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184313/0000/");
  }
  else if( insample==2522 ){
    mySample_xSec_ = 0.263;  
    mySample_nGen_ = 184990;
    mySample_sampleName_ = "ttZ_lep";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184349/0000/");
  }
  else if( insample==2510 ){
    mySample_xSec_ = 136.02;  
    mySample_nGen_ = 3299800;
    mySample_sampleName_ = "st_tchan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184110/0000/");
  }
  else if( insample==2511 ){
    mySample_xSec_ = 80.95;  
    mySample_nGen_ = 1695400;
    mySample_sampleName_ = "stbar_tchan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184034/0000/");
  }
  else if( insample==2512 ){
    mySample_xSec_ = 35.9;  
    mySample_nGen_ = 995600;
    mySample_sampleName_ = "st_tWchan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184134/0000/");
  }
  else if( insample==2513 ){
    mySample_xSec_ = 35.9;  
    mySample_nGen_ = 1000000;
    mySample_sampleName_ = "stbar_tWchan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184202/0000/");
  }
  else if( insample==2514 ){
    mySample_xSec_ = 10.32;  
    mySample_nGen_ = 613384;
    mySample_sampleName_ = "st_schan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184243/0000/");
  }
  else if( insample==2424 ){
    mySample_xSec_ = 118.7;  
    mySample_nGen_ = 994416;
    mySample_sampleName_ = "WW";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WW_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184512/0000/");
  }
  else if( insample==2423 ){
    mySample_xSec_ = 44.9;  
    mySample_nGen_ = 991232;
    mySample_sampleName_ = "WZ";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_184539/0000/");
  }
  else if( insample==2323 ){
    mySample_xSec_ = 15.4;  
    mySample_nGen_ = 996168;
    mySample_sampleName_ = "ZZ";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3_triggerTree_v1/151015_184614/0000/");
  }
  else if( insample==9125 ){
    mySample_xSec_ = 0.2934;// YR3 * BR(all)  
    mySample_nGen_ = 3933404;//199000;
    mySample_sampleName_ = "ttHTobb";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ttHTobb_M125_13TeV_powheg_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_triggerTree_v1/151015_183910/0000/");
  }
  else if( insample==8125 ){
    mySample_xSec_ = 0.2151;// YR3 * BR(all)  
    mySample_nGen_ = 3796398;//199000;
    mySample_sampleName_ = "ttHnonbb";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ttHToNonbb_M125_13TeV_powheg_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2_triggerTree_v1/151015_183937/0000/");
  }
  else if( insample==-13 ){
    mySample_sampleName_ = "SingleElectron_Run2015D_254231_258158";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleElectron/Run2015D-PromptReco-v3_GoldenJSON_254231_258158_triggerTree_v1/151012_183037/0000/");
  }
  else if( insample==-11 ){
    mySample_sampleName_ = "SingleMuon_Run2015D_254231_258158";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleMuon/Run2015D-PromptReco-v3_GoldenJSON_254231_258158_triggerTree_v1/151015_184732/0000/");
  }




  std::string s_end = "histo_" + str_jobN + ".root";
  if( Njobs==1 ) s_end = "histo.root";

  std::string histofilename = Form("HistoFiles/dileptonCR_optimize_treeReader_%s_%s", mySample_sampleName_.c_str(), s_end.c_str());

  TChain *chain = new TChain("triggeranalzyer/triggerTree");
  for( int iFile=0; iFile<int(mySample_inputDirs_.size()); iFile++ ){
    std::string treefilename = mySample_inputDirs_[iFile] + "trigger_analyzer*.root";
    std::cout << "  treefilename " << iFile << ": " << treefilename.c_str() << std::endl;
    chain->Add(treefilename.c_str());
  }
  //if( insample==2500 ) chain->Add("/eos/uscms/store/user/puigh/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1_triggerTree_v1/151015_183747/0000/trigger_analyzer_22*.root");
  //if( insample==2300 ) chain->Add("/eos/uscms/store/user/puigh/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3_triggerTree_v1/151015_183818/0000/trigger_analyzer_2*.root");

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


  TH1D* h_additionalJetEventId = new TH1D("h_additionalJetEventId",";additionalJetEventId", 201, -100-0.5, 101-0.5 );



  int NumCuts = 15;

  TH1D* h_ee_event_selection  = new TH1D("h_ee_event_selection",";cut", NumCuts, 0, NumCuts );
  h_ee_event_selection->GetXaxis()->SetBinLabel(1,"All");
  h_ee_event_selection->GetXaxis()->SetBinLabel(2,"HLT");
  h_ee_event_selection->GetXaxis()->SetBinLabel(3,"==2 l (20,20)");
  h_ee_event_selection->GetXaxis()->SetBinLabel(4,"==2 l (20,15)");
  h_ee_event_selection->GetXaxis()->SetBinLabel(5,">=2 j (20)");
  h_ee_event_selection->GetXaxis()->SetBinLabel(6,">=1 b (WPM)");
  // Look into Z mass veto here
  h_ee_event_selection->GetXaxis()->SetBinLabel(7,"mass veto");
  h_ee_event_selection->GetXaxis()->SetBinLabel(8,"MET>20");
  h_ee_event_selection->GetXaxis()->SetBinLabel(9,"MET>30");
  // Look into Z veto here
  h_ee_event_selection->GetXaxis()->SetBinLabel(10,"Zveto");
  // Look into MET cut here
  h_ee_event_selection->GetXaxis()->SetBinLabel(11,"MET>50");
  // Look into jet cuts here
  h_ee_event_selection->GetXaxis()->SetBinLabel(12,"==2 j (20)");
  h_ee_event_selection->GetXaxis()->SetBinLabel(13,"==2 j (30,20)");
  h_ee_event_selection->GetXaxis()->SetBinLabel(14,">=2 j (20,20)");
  h_ee_event_selection->GetXaxis()->SetBinLabel(15,">=2 j (30,20)");

  TH1D* h_mm_event_selection  = new TH1D("h_mm_event_selection",";cut", NumCuts, 0, NumCuts );
  h_mm_event_selection->GetXaxis()->SetBinLabel(1,"All");
  h_mm_event_selection->GetXaxis()->SetBinLabel(2,"HLT");
  h_mm_event_selection->GetXaxis()->SetBinLabel(3,"==2 l (20,20)");
  h_mm_event_selection->GetXaxis()->SetBinLabel(4,"==2 l (20,15)");
  h_mm_event_selection->GetXaxis()->SetBinLabel(5,">=2 j (20)");
  h_mm_event_selection->GetXaxis()->SetBinLabel(6,">=1 b (WPM)");
  // Look into Z mass veto here
  h_mm_event_selection->GetXaxis()->SetBinLabel(7,"mass veto");
  h_mm_event_selection->GetXaxis()->SetBinLabel(8,"MET>20");
  h_mm_event_selection->GetXaxis()->SetBinLabel(9,"MET>30");
  // Look into Z veto here
  h_mm_event_selection->GetXaxis()->SetBinLabel(10,"Zveto");
  // Look into MET cut here
  h_mm_event_selection->GetXaxis()->SetBinLabel(11,"MET>50");
  // Look into jet cuts here
  h_mm_event_selection->GetXaxis()->SetBinLabel(12,"==2 j (20)");
  h_mm_event_selection->GetXaxis()->SetBinLabel(13,"==2 j (30,20)");
  h_mm_event_selection->GetXaxis()->SetBinLabel(14,">=2 j (20,20)");
  h_mm_event_selection->GetXaxis()->SetBinLabel(15,">=2 j (30,20)");

  TH1D* h_em_event_selection  = new TH1D("h_em_event_selection",";cut", NumCuts, 0, NumCuts );
  h_em_event_selection->GetXaxis()->SetBinLabel(1,"All");
  h_em_event_selection->GetXaxis()->SetBinLabel(2,"HLT");
  h_em_event_selection->GetXaxis()->SetBinLabel(3,"==2 l (20,20)");
  h_em_event_selection->GetXaxis()->SetBinLabel(4,"==2 l (20,15)");
  h_em_event_selection->GetXaxis()->SetBinLabel(5,">=2 j (20)");
  h_em_event_selection->GetXaxis()->SetBinLabel(6,">=1 b (WPM)");
  // Look into Z mass veto here
  h_em_event_selection->GetXaxis()->SetBinLabel(7,"mass veto");
  h_em_event_selection->GetXaxis()->SetBinLabel(8,"MET>20");
  h_em_event_selection->GetXaxis()->SetBinLabel(9,"MET>30");
  // Look into Z veto here
  h_em_event_selection->GetXaxis()->SetBinLabel(10,"Zveto");
  // Look into MET cut here
  h_em_event_selection->GetXaxis()->SetBinLabel(11,"MET>50");
  // Look into jet cuts here
  h_em_event_selection->GetXaxis()->SetBinLabel(12,"==2 j (20)");
  h_em_event_selection->GetXaxis()->SetBinLabel(13,"==2 j (30,20)");
  h_em_event_selection->GetXaxis()->SetBinLabel(14,">=2 j (20,20)");
  h_em_event_selection->GetXaxis()->SetBinLabel(15,">=2 j (30,20)");

  TH1D* h_ll_event_selection  = new TH1D("h_ll_event_selection",";cut", NumCuts, 0, NumCuts );
  h_ll_event_selection->GetXaxis()->SetBinLabel(1,"All");
  h_ll_event_selection->GetXaxis()->SetBinLabel(2,"HLT");
  h_ll_event_selection->GetXaxis()->SetBinLabel(3,"==2 l (20,20)");
  h_ll_event_selection->GetXaxis()->SetBinLabel(4,"==2 l (20,15)");
  h_ll_event_selection->GetXaxis()->SetBinLabel(5,">=2 j (20)");
  h_ll_event_selection->GetXaxis()->SetBinLabel(6,">=1 b (WPM)");
  // Look into Z mass veto here
  h_ll_event_selection->GetXaxis()->SetBinLabel(7,"mass veto");
  h_ll_event_selection->GetXaxis()->SetBinLabel(8,"MET>20");
  h_ll_event_selection->GetXaxis()->SetBinLabel(9,"MET>30");
  // Look into Z veto here
  h_ll_event_selection->GetXaxis()->SetBinLabel(10,"Zveto");
  // Look into MET cut here
  h_ll_event_selection->GetXaxis()->SetBinLabel(11,"MET>50");
  // Look into jet cuts here
  h_ll_event_selection->GetXaxis()->SetBinLabel(12,"==2 j (20)");
  h_ll_event_selection->GetXaxis()->SetBinLabel(13,"==2 j (30,20)");
  h_ll_event_selection->GetXaxis()->SetBinLabel(14,">=2 j (20,20)");
  h_ll_event_selection->GetXaxis()->SetBinLabel(15,">=2 j (30,20)");



  TH1D* h_deltaR_jet_lep1 = new TH1D("h_deltaR_jet_lep1",";#DeltaR(jet,lep1)", 61, 0., 6.1 );
  TH1D* h_deltaR_jet_lep2 = new TH1D("h_deltaR_jet_lep2",";#DeltaR(jet,lep2)", 61, 0., 6.1 );




  //// Histograming 
  double metmax   = 500.;
  int NmetBins   = int( metmax/10. + 0.0001 );


  //
  // TwoLep
  //

  // After >=2 jet && >=1 btag 
  TH1D* h_ll_mht20_pt = new TH1D("h_ll_mht20_pt",";MHT (p_{T}>20)", NmetBins, 0, metmax );
  TH1D* h_ll_mht30_pt = new TH1D("h_ll_mht30_pt",";MHT (p_{T}>30)", NmetBins, 0, metmax );

  TH1D* h_ll_diLepMass = new TH1D("h_ll_diLepMass",";M(lep,lep)", 150, 0, 150 );

  TH2D* h_ll_mht20_diLepMass = new TH2D("h_ll_mht20_diLepMass",";MHT (p_{T}>20);M(lep,lep)", NmetBins, 0, metmax, 150, 0, 150 );
  TH2D* h_ll_mht30_diLepMass = new TH2D("h_ll_mht30_diLepMass",";MHT (p_{T}>30);M(lep,lep)", NmetBins, 0, metmax, 150, 0, 150 );

  // After Z veto
  TH1D* h_ll_met_pt = new TH1D("h_ll_met_pt",";MET", NmetBins, 0, metmax );
  TH1D* h_ll_mht20_pt_passZveto = new TH1D("h_ll_mht20_pt_passZveto",";MHT (p_{T}>20)", NmetBins, 0, metmax );
  TH1D* h_ll_mht30_pt_passZveto = new TH1D("h_ll_mht30_pt_passZveto",";MHT (p_{T}>30)", NmetBins, 0, metmax );

  // After MET
  TH1D* h_ll_numJet20 = new TH1D("h_ll_numJet20",";Number of Jets", 7, 1-0.5, 8-0.5 );
  TH1D* h_ll_numBtag20 = new TH1D("h_ll_numBtag20",";Number of b-tagged Jets", 5, 0-0.5, 5-0.5 );

  TH1D* h_ll_numJet30 = new TH1D("h_ll_numJet30",";Number of Jets", 7, 1-0.5, 8-0.5 );
  TH1D* h_ll_numBtag30 = new TH1D("h_ll_numBtag30",";Number of b-tagged Jets", 5, 0-0.5, 5-0.5 );

  TH2D* h_ll_numJet20_numJet30 = new TH2D("h_ll_numJet20_numJet30",";Number of Jets (p_{T}>20);Number of Jets (p_{T}>30)", 7, 1-0.5, 8-0.5, 7, 1-0.5, 8-0.5 );
  TH2D* h_ll_numBtag20_numBtag30 = new TH2D("h_ll_numBtag20_numBtag30",";Number of b-tagged Jets (p_{T}>20);Number of b-tagged Jets (p_{T}>30)", 5, 0-0.5, 5-0.5, 5, 0-0.5, 5-0.5 );



  //
  // TwoEle
  //

  // After >=2 jet && >=1 btag 
  TH1D* h_ee_mht20_pt = new TH1D("h_ee_mht20_pt",";MHT (p_{T}>20)", NmetBins, 0, metmax );
  TH1D* h_ee_mht30_pt = new TH1D("h_ee_mht30_pt",";MHT (p_{T}>30)", NmetBins, 0, metmax );

  TH1D* h_ee_diLepMass = new TH1D("h_ee_diLepMass",";M(lep,lep)", 150, 0, 150 );

  TH2D* h_ee_mht20_diLepMass = new TH2D("h_ee_mht20_diLepMass",";MHT (p_{T}>20);M(lep,lep)", NmetBins, 0, metmax, 150, 0, 150 );
  TH2D* h_ee_mht30_diLepMass = new TH2D("h_ee_mht30_diLepMass",";MHT (p_{T}>30);M(lep,lep)", NmetBins, 0, metmax, 150, 0, 150 );

  // After Z veto
  TH1D* h_ee_met_pt = new TH1D("h_ee_met_pt",";MET", NmetBins, 0, metmax );
  TH1D* h_ee_mht20_pt_passZveto = new TH1D("h_ee_mht20_pt_passZveto",";MHT (p_{T}>20)", NmetBins, 0, metmax );
  TH1D* h_ee_mht30_pt_passZveto = new TH1D("h_ee_mht30_pt_passZveto",";MHT (p_{T}>30)", NmetBins, 0, metmax );

  // After MET
  TH1D* h_ee_numJet20 = new TH1D("h_ee_numJet20",";Number of Jets", 7, 1-0.5, 8-0.5 );
  TH1D* h_ee_numBtag20 = new TH1D("h_ee_numBtag20",";Number of b-tagged Jets", 5, 0-0.5, 5-0.5 );

  TH1D* h_ee_numJet30 = new TH1D("h_ee_numJet30",";Number of Jets", 7, 1-0.5, 8-0.5 );
  TH1D* h_ee_numBtag30 = new TH1D("h_ee_numBtag30",";Number of b-tagged Jets", 5, 0-0.5, 5-0.5 );

  TH2D* h_ee_numJet20_numJet30 = new TH2D("h_ee_numJet20_numJet30",";Number of Jets (p_{T}>20);Number of Jets (p_{T}>30)", 7, 1-0.5, 8-0.5, 7, 1-0.5, 8-0.5 );
  TH2D* h_ee_numBtag20_numBtag30 = new TH2D("h_ee_numBtag20_numBtag30",";Number of b-tagged Jets (p_{T}>20);Number of b-tagged Jets (p_{T}>30)", 5, 0-0.5, 5-0.5, 5, 0-0.5, 5-0.5 );


  //
  // TwoMu
  //

  // After >=2 jet && >=1 btag 
  TH1D* h_mm_mht20_pt = new TH1D("h_mm_mht20_pt",";MHT (p_{T}>20)", NmetBins, 0, metmax );
  TH1D* h_mm_mht30_pt = new TH1D("h_mm_mht30_pt",";MHT (p_{T}>30)", NmetBins, 0, metmax );

  TH1D* h_mm_diLepMass = new TH1D("h_mm_diLepMass",";M(lep,lep)", 150, 0, 150 );

  TH2D* h_mm_mht20_diLepMass = new TH2D("h_mm_mht20_diLepMass",";MHT (p_{T}>20);M(lep,lep)", NmetBins, 0, metmax, 150, 0, 150 );
  TH2D* h_mm_mht30_diLepMass = new TH2D("h_mm_mht30_diLepMass",";MHT (p_{T}>30);M(lep,lep)", NmetBins, 0, metmax, 150, 0, 150 );

  // After Z veto
  TH1D* h_mm_met_pt = new TH1D("h_mm_met_pt",";MET", NmetBins, 0, metmax );
  TH1D* h_mm_mht20_pt_passZveto = new TH1D("h_mm_mht20_pt_passZveto",";MHT (p_{T}>20)", NmetBins, 0, metmax );
  TH1D* h_mm_mht30_pt_passZveto = new TH1D("h_mm_mht30_pt_passZveto",";MHT (p_{T}>30)", NmetBins, 0, metmax );

  // After MET
  TH1D* h_mm_numJet20 = new TH1D("h_mm_numJet20",";Number of Jets", 7, 1-0.5, 8-0.5 );
  TH1D* h_mm_numBtag20 = new TH1D("h_mm_numBtag20",";Number of b-tagged Jets", 5, 0-0.5, 5-0.5 );

  TH1D* h_mm_numJet30 = new TH1D("h_mm_numJet30",";Number of Jets", 7, 1-0.5, 8-0.5 );
  TH1D* h_mm_numBtag30 = new TH1D("h_mm_numBtag30",";Number of b-tagged Jets", 5, 0-0.5, 5-0.5 );

  TH2D* h_mm_numJet20_numJet30 = new TH2D("h_mm_numJet20_numJet30",";Number of Jets (p_{T}>20);Number of Jets (p_{T}>30)", 7, 1-0.5, 8-0.5, 7, 1-0.5, 8-0.5 );
  TH2D* h_mm_numBtag20_numBtag30 = new TH2D("h_mm_numBtag20_numBtag30",";Number of b-tagged Jets (p_{T}>20);Number of b-tagged Jets (p_{T}>30)", 5, 0-0.5, 5-0.5, 5, 0-0.5, 5-0.5 );

  //
  // EleMu
  //

  // After >=2 jet && >=1 btag 
  TH1D* h_em_mht20_pt = new TH1D("h_em_mht20_pt",";MHT (p_{T}>20)", NmetBins, 0, metmax );
  TH1D* h_em_mht30_pt = new TH1D("h_em_mht30_pt",";MHT (p_{T}>30)", NmetBins, 0, metmax );

  TH1D* h_em_diLepMass = new TH1D("h_em_diLepMass",";M(lep,lep)", 150, 0, 150 );

  TH2D* h_em_mht20_diLepMass = new TH2D("h_em_mht20_diLepMass",";MHT (p_{T}>20);M(lep,lep)", NmetBins, 0, metmax, 150, 0, 150 );
  TH2D* h_em_mht30_diLepMass = new TH2D("h_em_mht30_diLepMass",";MHT (p_{T}>30);M(lep,lep)", NmetBins, 0, metmax, 150, 0, 150 );

  // After Z veto
  TH1D* h_em_met_pt = new TH1D("h_em_met_pt",";MET", NmetBins, 0, metmax );
  TH1D* h_em_mht20_pt_passZveto = new TH1D("h_em_mht20_pt_passZveto",";MHT (p_{T}>20)", NmetBins, 0, metmax );
  TH1D* h_em_mht30_pt_passZveto = new TH1D("h_em_mht30_pt_passZveto",";MHT (p_{T}>30)", NmetBins, 0, metmax );

  // After MET
  TH1D* h_em_numJet20 = new TH1D("h_em_numJet20",";Number of Jets", 7, 1-0.5, 8-0.5 );
  TH1D* h_em_numBtag20 = new TH1D("h_em_numBtag20",";Number of b-tagged Jets", 5, 0-0.5, 5-0.5 );

  TH1D* h_em_numJet30 = new TH1D("h_em_numJet30",";Number of Jets", 7, 1-0.5, 8-0.5 );
  TH1D* h_em_numBtag30 = new TH1D("h_em_numBtag30",";Number of b-tagged Jets", 5, 0-0.5, 5-0.5 );

  TH2D* h_em_numJet20_numJet30 = new TH2D("h_em_numJet20_numJet30",";Number of Jets (p_{T}>20);Number of Jets (p_{T}>30)", 7, 1-0.5, 8-0.5, 7, 1-0.5, 8-0.5 );
  TH2D* h_em_numBtag20_numBtag30 = new TH2D("h_em_numBtag20_numBtag30",";Number of b-tagged Jets (p_{T}>20);Number of b-tagged Jets (p_{T}>30)", 5, 0-0.5, 5-0.5, 5, 0-0.5, 5-0.5 );



  // TwoLep

  int NcsvBins = 106;
  int NumPtBins = 6;


  TH1D* h_ge2j20_probeJet_jet1_flavor = new TH1D("h_ge2j20_probeJet_jet1_flavor",";flavor", 17, -6-0.5, 11-0.5 );

  TH1D* h_e2j20_tagJet_flavor = new TH1D("h_e2j20_tagJet_flavor",";flavor", 17, -6-0.5, 11-0.5 );
  TH1D* h_j30j20_tagJet_flavor = new TH1D("h_j30j20_tagJet_flavor",";flavor", 17, -6-0.5, 11-0.5 );
  TH1D* h_ge2j20_tagJet_flavor = new TH1D("h_ge2j20_tagJet_flavor",";flavor", 17, -6-0.5, 11-0.5 );

  TH1D* h_e2j20_numProbeJet20 = new TH1D("h_e2j20_numProbeJet20",";Number probe jets (p_{T}>20)", 6, 0-0.5, 6-0.5 );
  TH1D* h_j30j20_numProbeJet20 = new TH1D("h_j30j20_numProbeJet20",";Number probe jets (p_{T}>20)", 6, 0-0.5, 6-0.5 );
  TH1D* h_ge2j20_numProbeJet20 = new TH1D("h_ge2j20_numProbeJet20",";Number probe jets (p_{T}>20)", 6, 0-0.5, 6-0.5 );

  TH2D* h_e2j20_numProbeJet20_numBProbeJet20  = new TH2D("h_e2j20_numProbeJet20_numBProbeJet20",";Number probe jets (p_{T}>20);Number b probe jets (p_{T}>20)", 6, 0-0.5, 6-0.5, 6, 0-0.5, 6-0.5 );
  TH2D* h_j30j20_numProbeJet20_numBProbeJet20 = new TH2D("h_j30j20_numProbeJet20_numBProbeJet20",";Number probe jets (p_{T}>20);Number b probe jets (p_{T}>20)", 6, 0-0.5, 6-0.5, 6, 0-0.5, 6-0.5 );
  TH2D* h_ge2j20_numProbeJet20_numBProbeJet20 = new TH2D("h_ge2j20_numProbeJet20_numBProbeJet20",";Number probe jets (p_{T}>20);Number b probe jets (p_{T}>20)", 6, 0-0.5, 6-0.5, 6, 0-0.5, 6-0.5 );


  TH1D* h_e2j20_b_probeJet_csv[NumPtBins];
  TH1D* h_e2j20_nonb_probeJet_csv[NumPtBins];
  TH1D* h_e2j20_probeJet_flavor[NumPtBins];

  TH1D* h_e2j20_b_deltaR_probeJet_tagJet[NumPtBins];
  TH1D* h_e2j20_b_deltaEta_probeJet_tagJet[NumPtBins];
  TH1D* h_e2j20_b_deltaPhi_probeJet_tagJet[NumPtBins];
  TH1D* h_e2j20_nonb_deltaR_probeJet_tagJet[NumPtBins];
  TH1D* h_e2j20_nonb_deltaEta_probeJet_tagJet[NumPtBins];
  TH1D* h_e2j20_nonb_deltaPhi_probeJet_tagJet[NumPtBins];


  TH1D* h_j30j20_b_probeJet_csv[NumPtBins];
  TH1D* h_j30j20_nonb_probeJet_csv[NumPtBins];
  TH1D* h_j30j20_probeJet_flavor[NumPtBins];

  TH1D* h_j30j20_b_deltaR_probeJet_tagJet[NumPtBins];
  TH1D* h_j30j20_b_deltaEta_probeJet_tagJet[NumPtBins];
  TH1D* h_j30j20_b_deltaPhi_probeJet_tagJet[NumPtBins];
  TH1D* h_j30j20_nonb_deltaR_probeJet_tagJet[NumPtBins];
  TH1D* h_j30j20_nonb_deltaEta_probeJet_tagJet[NumPtBins];
  TH1D* h_j30j20_nonb_deltaPhi_probeJet_tagJet[NumPtBins];


  

  TH1D* h_ge2j20_b_probeJet_firstJet_csv[NumPtBins];
  TH1D* h_ge2j20_nonb_probeJet_firstJet_csv[NumPtBins];
  TH1D* h_ge2j20_probeJet_firstJet_flavor[NumPtBins];

  TH1D* h_ge2j20_b_probeJet_csv[NumPtBins];
  TH1D* h_ge2j20_nonb_probeJet_csv[NumPtBins];
  TH1D* h_ge2j20_probeJet_flavor[NumPtBins];

  TH1D* h_ge2j20_b_deltaR_probeJet_tagJet[NumPtBins];
  TH1D* h_ge2j20_b_deltaEta_probeJet_tagJet[NumPtBins];
  TH1D* h_ge2j20_b_deltaPhi_probeJet_tagJet[NumPtBins];
  TH1D* h_ge2j20_nonb_deltaR_probeJet_tagJet[NumPtBins];
  TH1D* h_ge2j20_nonb_deltaEta_probeJet_tagJet[NumPtBins];
  TH1D* h_ge2j20_nonb_deltaPhi_probeJet_tagJet[NumPtBins];


  for( int iPt=0; iPt<NumPtBins; iPt++ ){
    h_e2j20_b_probeJet_csv[iPt] = new TH1D(Form("h_e2j20_b_probeJet_csv_Pt%d",iPt),";jet CSV", NcsvBins, -0.05, 1.01 );
    h_e2j20_nonb_probeJet_csv[iPt] = new TH1D(Form("h_e2j20_nonb_probeJet_csv_Pt%d",iPt),";jet CSV", NcsvBins, -0.05, 1.01 );
    h_e2j20_probeJet_flavor[iPt] = new TH1D(Form("h_e2j20_probeJet_flavor_Pt%d",iPt),";flavor", 17, -6-0.5, 11-0.5 );

    h_e2j20_b_deltaR_probeJet_tagJet[iPt] = new TH1D(Form("h_e2j20_b_deltaR_probeJet_tagJet_Pt%d",iPt),";#DeltaR(tag,probe)", 65, -0.01, 6.49 );
    h_e2j20_b_deltaEta_probeJet_tagJet[iPt] = new TH1D(Form("h_e2j20_b_deltaEta_probeJet_tagJet_Pt%d",iPt),";#Delta#eta(tag,probe)", 50, -0.01, 4.99 );
    h_e2j20_b_deltaPhi_probeJet_tagJet[iPt] = new TH1D(Form("h_e2j20_b_deltaPhi_probeJet_tagJet_Pt%d",iPt),";#Delta#phi(tag,probe)", 32, -0.01, 3.19 );
    h_e2j20_nonb_deltaR_probeJet_tagJet[iPt] = new TH1D(Form("h_e2j20_nonb_deltaR_probeJet_tagJet_Pt%d",iPt),";#DeltaR(tag,probe)", 65, -0.01, 6.49 );
    h_e2j20_nonb_deltaEta_probeJet_tagJet[iPt] = new TH1D(Form("h_e2j20_nonb_deltaEta_probeJet_tagJet_Pt%d",iPt),";#Delta#eta(tag,probe)", 50, -0.01, 4.99 );
    h_e2j20_nonb_deltaPhi_probeJet_tagJet[iPt] = new TH1D(Form("h_e2j20_nonb_deltaPhi_probeJet_tagJet_Pt%d",iPt),";#Delta#phi(tag,probe)", 32, -0.01, 3.19 );

    h_j30j20_b_probeJet_csv[iPt] = new TH1D(Form("h_j30j20_b_probeJet_csv_Pt%d",iPt),";jet CSV", NcsvBins, -0.05, 1.01 );
    h_j30j20_nonb_probeJet_csv[iPt] = new TH1D(Form("h_j30j20_nonb_probeJet_csv_Pt%d",iPt),";jet CSV", NcsvBins, -0.05, 1.01 );
    h_j30j20_probeJet_flavor[iPt] = new TH1D(Form("h_j30j20_probeJet_flavor_Pt%d",iPt),";flavor", 17, -6-0.5, 11-0.5 );

    h_j30j20_b_deltaR_probeJet_tagJet[iPt] = new TH1D(Form("h_j30j20_b_deltaR_probeJet_tagJet_Pt%d",iPt),";#DeltaR(tag,probe)", 65, -0.01, 6.49 );
    h_j30j20_b_deltaEta_probeJet_tagJet[iPt] = new TH1D(Form("h_j30j20_b_deltaEta_probeJet_tagJet_Pt%d",iPt),";#Delta#eta(tag,probe)", 50, -0.01, 4.99 );
    h_j30j20_b_deltaPhi_probeJet_tagJet[iPt] = new TH1D(Form("h_j30j20_b_deltaPhi_probeJet_tagJet_Pt%d",iPt),";#Delta#phi(tag,probe)", 32, -0.01, 3.19 );
    h_j30j20_nonb_deltaR_probeJet_tagJet[iPt] = new TH1D(Form("h_j30j20_nonb_deltaR_probeJet_tagJet_Pt%d",iPt),";#DeltaR(tag,probe)", 65, -0.01, 6.49 );
    h_j30j20_nonb_deltaEta_probeJet_tagJet[iPt] = new TH1D(Form("h_j30j20_nonb_deltaEta_probeJet_tagJet_Pt%d",iPt),";#Delta#eta(tag,probe)", 50, -0.01, 4.99 );
    h_j30j20_nonb_deltaPhi_probeJet_tagJet[iPt] = new TH1D(Form("h_j30j20_nonb_deltaPhi_probeJet_tagJet_Pt%d",iPt),";#Delta#phi(tag,probe)", 32, -0.01, 3.19 );

    h_ge2j20_b_probeJet_firstJet_csv[iPt] = new TH1D(Form("h_ge2j20_b_probeJet_firstJet_csv_Pt%d",iPt),";jet CSV", NcsvBins, -0.05, 1.01 );
    h_ge2j20_nonb_probeJet_firstJet_csv[iPt] = new TH1D(Form("h_ge2j20_nonb_probeJet_firstJet_csv_Pt%d",iPt),";jet CSV", NcsvBins, -0.05, 1.01 );
    h_ge2j20_probeJet_firstJet_flavor[iPt] = new TH1D(Form("h_ge2j20_probeJet_firstJet_flavor_Pt%d",iPt),";flavor", 17, -6-0.5, 11-0.5 );

    h_ge2j20_b_probeJet_csv[iPt] = new TH1D(Form("h_ge2j20_b_probeJet_csv_Pt%d",iPt),";jet CSV", NcsvBins, -0.05, 1.01 );
    h_ge2j20_nonb_probeJet_csv[iPt] = new TH1D(Form("h_ge2j20_nonb_probeJet_csv_Pt%d",iPt),";jet CSV", NcsvBins, -0.05, 1.01 );
    h_ge2j20_probeJet_flavor[iPt] = new TH1D(Form("h_ge2j20_probeJet_flavor_Pt%d",iPt),";flavor", 17, -6-0.5, 11-0.5 );

    h_ge2j20_b_deltaR_probeJet_tagJet[iPt] = new TH1D(Form("h_ge2j20_b_deltaR_probeJet_tagJet_Pt%d",iPt),";#DeltaR(tag,probe)", 65, -0.01, 6.49 );
    h_ge2j20_b_deltaEta_probeJet_tagJet[iPt] = new TH1D(Form("h_ge2j20_b_deltaEta_probeJet_tagJet_Pt%d",iPt),";#Delta#eta(tag,probe)", 50, -0.01, 4.99 );
    h_ge2j20_b_deltaPhi_probeJet_tagJet[iPt] = new TH1D(Form("h_ge2j20_b_deltaPhi_probeJet_tagJet_Pt%d",iPt),";#Delta#phi(tag,probe)", 32, -0.01, 3.19 );
    h_ge2j20_nonb_deltaR_probeJet_tagJet[iPt] = new TH1D(Form("h_ge2j20_nonb_deltaR_probeJet_tagJet_Pt%d",iPt),";#DeltaR(tag,probe)", 65, -0.01, 6.49 );
    h_ge2j20_nonb_deltaEta_probeJet_tagJet[iPt] = new TH1D(Form("h_ge2j20_nonb_deltaEta_probeJet_tagJet_Pt%d",iPt),";#Delta#eta(tag,probe)", 50, -0.01, 4.99 );
    h_ge2j20_nonb_deltaPhi_probeJet_tagJet[iPt] = new TH1D(Form("h_ge2j20_nonb_deltaPhi_probeJet_tagJet_Pt%d",iPt),";#Delta#phi(tag,probe)", 32, -0.01, 3.19 );

  }


  //////////////////////////////////////////////////////////////////////////
  /////
  //////////////////////////////////////////////////////////////////////////

  std::string ee_path_name  = "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  std::string em_path_name1 = "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v";
  std::string em_path_name2 = "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v";
  std::string mm_path_name1 = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
  std::string mm_path_name2 = "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";


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


    int additionalJetEventId = -99;
    if( insample==2500 ){
      additionalJetEventId = eve->additionalJetEventId_;

      bool keepTTbarEvent = true;
      if( ttCat_>=0 ){
	if( ttCat_==0 && additionalJetEventId==0 ) keepTTbarEvent = true;
	else if( ttCat_==1 && (additionalJetEventId>=41 && additionalJetEventId<=45) ) keepTTbarEvent = true;
	else if( ttCat_==2 && additionalJetEventId==51 ) keepTTbarEvent = true;
	else if( ttCat_==3 && additionalJetEventId==52 ) keepTTbarEvent = true;
	else if( ttCat_==4 && (additionalJetEventId>=53 && additionalJetEventId<=55) ) keepTTbarEvent = true;
	else keepTTbarEvent = false;
      }

      // if( keepTTbarEvent==false && ttCat_>=0 ){
      // 	std::cout << " ERROR!! keepTTbarEvent == false! additionalJetEventId = " << additionalJetEventId << " and ttCat_ = " << ttCat_ << std::endl;
      // }

      if( !keepTTbarEvent ) continue;
      //std::cout << " I have made it past the gate! " << std::endl;
    }

    h_additionalJetEventId->Fill(additionalJetEventId);

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


    h_numPVs->Fill(numPVs,wgt);


    wgt *= scalePU;

    h_numPVs_PUwgt->Fill(numPVs,wgt);

    ///////////////////
    ////// selections
    ///////////////////


    h_ll_event_selection->Fill(0.5, wgt); // all
    h_ee_event_selection->Fill(0.5, wgt); // all
    h_mm_event_selection->Fill(0.5, wgt); // all
    h_em_event_selection->Fill(0.5, wgt); // all

    bool pass_trigger_ee = false;
    bool pass_trigger_mm1 = false;
    bool pass_trigger_mm2 = false;
    bool pass_trigger_em1 = false;
    bool pass_trigger_em2 = false;

    vint hlt_accept = eve->hlt_accept_;
    vstring hlt_name = eve->hlt_name_;

    int Nhlt = int( hlt_accept.size() );
    for( int iHLT=0; iHLT<Nhlt; iHLT++ ){
      std::string name = hlt_name[iHLT];
      int accept = hlt_accept[iHLT];

      if( name.find(ee_path_name)!=std::string::npos ){
	if( accept ) pass_trigger_ee = true;
	else         pass_trigger_ee = false;
      }
      else if( name.find(mm_path_name1)!=std::string::npos ){
	if( accept ) pass_trigger_mm1 = true;
	else         pass_trigger_mm1 = false;
      }
      else if( name.find(mm_path_name2)!=std::string::npos ){
	if( accept ) pass_trigger_mm2 = true;
	else         pass_trigger_mm2 = false;
      }
      else if( name.find(em_path_name1)!=std::string::npos ){
	if( accept ) pass_trigger_em1 = true;
	else         pass_trigger_em1 = false;
      }
      else if( name.find(em_path_name2)!=std::string::npos ){
	if( accept ) pass_trigger_em2 = true;
	else         pass_trigger_em2 = false;
      }
    }

    bool pass_trigger_mm = (pass_trigger_mm1 || pass_trigger_mm2);
    bool pass_trigger_em = (pass_trigger_em1 || pass_trigger_em2);

    bool pass_trigger_ll = (pass_trigger_ee || pass_trigger_mm || pass_trigger_em);

    // put in stuff for data
    // if( insample<0 ){
    //   pass_trigger_ele = ( insample==-13 && pass_trigger_ele );
    //   pass_trigger_mu  = ( insample==-11 && pass_trigger_mu  );
    // }

    if( pass_trigger_ll ) h_ll_event_selection->Fill(0.5+1, wgt); // all
    if( pass_trigger_ee ) h_ee_event_selection->Fill(0.5+1, wgt); // all
    if( pass_trigger_mm ) h_mm_event_selection->Fill(0.5+1, wgt); // all
    if( pass_trigger_em ) h_em_event_selection->Fill(0.5+1, wgt); // all


    // Pass Trigger selection
    if( !(pass_trigger_ll) ) continue;


    //
    // Lepton selection
    //

    std::vector<int> ind_ele;
    std::vector<int> ind_mu;
    std::vector<int> ind_ele_loose;
    std::vector<int> ind_mu_loose;
    for( int iLep=0; iLep<int(eve->lepton_pt_.size()); iLep++ ){
      bool isMuon = ( eve->lepton_isMuon_[iLep]==1 );

      bool isSpring15M = ( eve->lepton_isSpring15M_[iLep]==1 );
      bool isTight = ( eve->lepton_isTight_[iLep]==1 );
      bool isLoose = ( eve->lepton_isLoose_[iLep]==1 );
      bool isCrack = ( eve->lepton_inCrack_[iLep]==1 );

      double pt  = eve->lepton_pt_[iLep];
      double eta = eve->lepton_eta_[iLep];

      double relIso = eve->lepton_relIsoR04_[iLep];

      if( isMuon ){
	if( pt > 20 && abs(eta)<2.4 && isLoose && relIso < 0.15 ) ind_mu.push_back(iLep);
	if( pt > 15 && abs(eta)<2.4 && isLoose && relIso < 0.15 ) ind_mu_loose.push_back(iLep);
      }

      if( !isMuon ){
	if( pt > 20 && abs(eta)<2.4 && isSpring15M && !isCrack ) ind_ele.push_back(iLep);
	if( pt > 15 && abs(eta)<2.4 && isSpring15M && !isCrack ) ind_ele_loose.push_back(iLep);
      }
    }

    int numEle = int( ind_ele.size() );
    int numMu  = int( ind_mu.size() );

    int numEle_loose = int( ind_ele_loose.size() );
    int numMu_loose  = int( ind_mu_loose.size() );

    bool TwoTightEle = ( numEle==2 && numEle_loose==2 && numMu_loose==0 );
    bool TwoTightMu  = ( numMu==2 && numMu_loose==2 && numEle_loose==0 );


    bool TwoEle = ( numEle>=1 && numEle_loose==2 && numMu_loose==0 );
    bool TwoMu  = ( numMu>=1 && numMu_loose==2 && numEle_loose==0 );

    bool TightEleMu = ( numEle==1 && numMu==1 && numMu_loose==1 && numEle_loose==1 );

    bool EleMu = ( (numEle==1 || numMu==1) && numMu_loose==1 && numEle_loose==1 );


    // Require trigger fired for each type
    TwoTightEle = ( TwoTightEle && pass_trigger_ee );
    TwoTightMu  = ( TwoTightMu  && pass_trigger_mm );
    TightEleMu  = ( TightEleMu  && pass_trigger_em );

    TwoEle = ( TwoEle && pass_trigger_ee );
    TwoMu  = ( TwoMu  && pass_trigger_mm );
    EleMu  = ( EleMu  && pass_trigger_em );


    bool TwoTightLep = ( TwoTightEle || TwoTightMu || TightEleMu );
    bool TwoLep = ( TwoEle || TwoMu || EleMu );


    if( TwoTightLep ) h_ll_event_selection->Fill(0.5+2, wgt); // all
    if( TwoTightEle ) h_ee_event_selection->Fill(0.5+2, wgt); // all
    if( TwoTightMu  ) h_mm_event_selection->Fill(0.5+2, wgt); // all
    if( TightEleMu  ) h_em_event_selection->Fill(0.5+2, wgt); // all


    if( TwoLep ) h_ll_event_selection->Fill(0.5+3, wgt); // all
    if( TwoEle ) h_ee_event_selection->Fill(0.5+3, wgt); // all
    if( TwoMu  ) h_mm_event_selection->Fill(0.5+3, wgt); // all
    if( EleMu  ) h_em_event_selection->Fill(0.5+3, wgt); // all


    // Require exactly one lepton
    if( !TwoLep ) continue;


    int lepInd1 = -1, lepInd2 = -1;
    if( TwoEle ){
      lepInd1 = ind_ele_loose[0];
      lepInd2 = ind_ele_loose[1];
    }
    else if( TwoMu ){
      lepInd1 = ind_mu_loose[0];
      lepInd2 = ind_mu_loose[1];
    }
    else if( EleMu ){
      lepInd1 = ind_ele_loose[0];
      lepInd2 = ind_mu_loose[0];
    }

    if( lepInd1<0 || lepInd2<0 ) std::cout << " OOOHHHH SSHHIIITTTT" << std::endl;

    TLorentzVector myLep1;
    myLep1.SetPtEtaPhiE( eve->lepton_pt_[lepInd1], eve->lepton_eta_[lepInd1], eve->lepton_phi_[lepInd1], eve->lepton_energy_[lepInd1] );

    TLorentzVector myLep2;
    myLep2.SetPtEtaPhiE( eve->lepton_pt_[lepInd2], eve->lepton_eta_[lepInd2], eve->lepton_phi_[lepInd2], eve->lepton_energy_[lepInd2] );

    int charge1 = eve->lepton_trkCharge_[lepInd1];
    int charge2 = eve->lepton_trkCharge_[lepInd2];

    
    //
    // REQUIRE OPPOSITE CHARGE
    //

    if( !(charge1*charge2<0) ) continue;



    vdouble jetPts30;
    int numJet30 = 0;
    int numBtag30 = 0;
    double HT30=0;
    TLorentzVector sumJet30;
    bool firstJet30 = false;

    vdouble jetPts20;
    int numJet20 = 0;
    int numBtag20 = 0;
    double HT20=0;
    TLorentzVector sumJet20;
    bool firstJet20 = false;
    std::vector<int> ind_jet20;

    for( int iJet=0; iJet<int(eve->jet_pt_.size()); iJet++ ){
      double pt  = eve->jet_pt_[iJet];
      double eta = eve->jet_eta_[iJet];
      double phi = eve->jet_phi_[iJet];
      double csv =  eve->jet_csv_[iJet];

      if( csv<0 && csv>-9 ) csv = -0.2;
      else if( csv < -5 )   csv = -0.4;

      if( csv > 1.0 ) csv = 1.0;

      if( !(pt>20. && fabs(eta)<2.4) ) continue;

      TLorentzVector myJet;
      myJet.SetPtEtaPhiE( eve->jet_pt_[iJet], eve->jet_eta_[iJet], eve->jet_phi_[iJet], eve->jet_energy_[iJet] );

      double dR1 = myJet.DeltaR(myLep1);
      h_deltaR_jet_lep1->Fill(dR1,wgt);

      if( dR1 < 0.4 ) continue;

      double dR2 = myJet.DeltaR(myLep2);
      h_deltaR_jet_lep2->Fill(dR2,wgt);

      if( dR2 < 0.4 ) continue;

      ind_jet20.push_back(iJet);

      if( !firstJet20 ){
	firstJet20 = true;
	sumJet20 = myJet;
      }
      else sumJet20 += myJet;

      HT20 += pt;
      jetPts20.push_back(pt);

      numJet20++;
      if( csv > 0.890 ) numBtag20++;

      if( !(pt>30. && fabs(eta)<2.4) ) continue;


      if( !firstJet30 ){
	firstJet30 = true;
	sumJet30 = myJet;
      }
      else sumJet30 += myJet;

      HT30 += pt;
      jetPts30.push_back(pt);

      numJet30++;
      if( csv > 0.890 ) numBtag30++;
    }


    if( !(numJet20>=2) ) continue;

    if( TwoLep ) h_ll_event_selection->Fill(0.5+4, wgt); // all
    if( TwoEle ) h_ee_event_selection->Fill(0.5+4, wgt); // all
    if( TwoMu  ) h_mm_event_selection->Fill(0.5+4, wgt); // all
    if( EleMu  ) h_em_event_selection->Fill(0.5+4, wgt); // all

    if( !(numBtag20>=1) ) continue;

    if( TwoLep ) h_ll_event_selection->Fill(0.5+5, wgt); // all
    if( TwoEle ) h_ee_event_selection->Fill(0.5+5, wgt); // all
    if( TwoMu  ) h_mm_event_selection->Fill(0.5+5, wgt); // all
    if( EleMu  ) h_em_event_selection->Fill(0.5+5, wgt); // all


    TLorentzVector diLep = myLep1 + myLep2;

    TLorentzVector sum30 = diLep + sumJet30;
    TLorentzVector sum20 = diLep + sumJet20;

    double met = eve->pfMET_pt_;
    double met_phi = eve->pfMET_phi_;
    double mht30 = sum30.Pt();
    double mht30_phi = sum30.Phi();
    double mht20 = sum20.Pt();
    double mht20_phi = sum20.Phi();


    double diLepMass = diLep.M();

    bool passZveto = ( EleMu || 
		       (diLepMass < (65.5 + 3*mht20/8)) || 
		       (diLepMass > (108 - mht20/4)) || 
		       (diLepMass < (79 - 3*mht20/4)) || 
		       (diLepMass > (99 + mht20/2)) 
		       );


    if( TwoLep ){
      h_ll_mht20_pt->Fill(mht20,wgt);
      h_ll_mht30_pt->Fill(mht30,wgt);

      h_ll_diLepMass->Fill(diLepMass,wgt);

      h_ll_mht20_diLepMass->Fill(mht20,diLepMass,wgt);
      h_ll_mht30_diLepMass->Fill(mht30,diLepMass,wgt);
    }
    if( TwoEle ){
      h_ee_mht20_pt->Fill(mht20,wgt);
      h_ee_mht30_pt->Fill(mht30,wgt);

      h_ee_diLepMass->Fill(diLepMass,wgt);

      h_ee_mht20_diLepMass->Fill(mht20,diLepMass,wgt);
      h_ee_mht30_diLepMass->Fill(mht30,diLepMass,wgt);
    }
    if( TwoMu  ){
      h_mm_mht20_pt->Fill(mht20,wgt);
      h_mm_mht30_pt->Fill(mht30,wgt);

      h_mm_diLepMass->Fill(diLepMass,wgt);

      h_mm_mht20_diLepMass->Fill(mht20,diLepMass,wgt);
      h_mm_mht30_diLepMass->Fill(mht30,diLepMass,wgt);
    }
    if( EleMu  ){
      h_em_mht20_pt->Fill(mht20,wgt);
      h_em_mht30_pt->Fill(mht30,wgt);

      h_em_diLepMass->Fill(diLepMass,wgt);

      h_em_mht20_diLepMass->Fill(mht20,diLepMass,wgt);
      h_em_mht30_diLepMass->Fill(mht30,diLepMass,wgt);
    }


    if( EleMu || fabs(diLepMass-91)>10 ){
      if( TwoLep ) h_ll_event_selection->Fill(0.5+6, wgt); // all
      if( TwoEle ) h_ee_event_selection->Fill(0.5+6, wgt); // all
      if( TwoMu  ) h_mm_event_selection->Fill(0.5+6, wgt); // all
      if( EleMu  ) h_em_event_selection->Fill(0.5+6, wgt); // all

      if( met>20 ){
	if( TwoLep ) h_ll_event_selection->Fill(0.5+7, wgt); // all
	if( TwoEle ) h_ee_event_selection->Fill(0.5+7, wgt); // all
	if( TwoMu  ) h_mm_event_selection->Fill(0.5+7, wgt); // all
	if( EleMu  ) h_em_event_selection->Fill(0.5+7, wgt); // all

	if( met>30 ){
	  if( TwoLep ) h_ll_event_selection->Fill(0.5+8, wgt); // all
	  if( TwoEle ) h_ee_event_selection->Fill(0.5+8, wgt); // all
	  if( TwoMu  ) h_mm_event_selection->Fill(0.5+8, wgt); // all
	  if( EleMu  ) h_em_event_selection->Fill(0.5+8, wgt); // all
	}
      }
    }

    //
    // REQUIRE PASS VETO
    //
    if( !passZveto ) continue;


    if( TwoLep ) h_ll_event_selection->Fill(0.5+9, wgt); // all
    if( TwoEle ) h_ee_event_selection->Fill(0.5+9, wgt); // all
    if( TwoMu  ) h_mm_event_selection->Fill(0.5+9, wgt); // all
    if( EleMu  ) h_em_event_selection->Fill(0.5+9, wgt); // all


    if( TwoLep ){
      h_ll_met_pt->Fill(met,wgt);
      h_ll_mht20_pt_passZveto->Fill(mht20,wgt);
      h_ll_mht30_pt_passZveto->Fill(mht30,wgt);
    }
    if( TwoEle ){
      h_ee_met_pt->Fill(met,wgt);
      h_ee_mht20_pt_passZveto->Fill(mht20,wgt);
      h_ee_mht30_pt_passZveto->Fill(mht30,wgt);
    }
    if( TwoMu  ){
      h_mm_met_pt->Fill(met,wgt);
      h_mm_mht20_pt_passZveto->Fill(mht20,wgt);
      h_mm_mht30_pt_passZveto->Fill(mht30,wgt);
    }
    if( EleMu  ){
      h_em_met_pt->Fill(met,wgt);
      h_em_mht20_pt_passZveto->Fill(mht20,wgt);
      h_em_mht30_pt_passZveto->Fill(mht30,wgt);
    }


    //
    // REQUIRE PASS MET CUT
    //
    if( !(met>50) ) continue;

    if( TwoLep ) h_ll_event_selection->Fill(0.5+10, wgt); // all
    if( TwoEle ) h_ee_event_selection->Fill(0.5+10, wgt); // all
    if( TwoMu  ) h_mm_event_selection->Fill(0.5+10, wgt); // all
    if( EleMu  ) h_em_event_selection->Fill(0.5+10, wgt); // all


    if( TwoLep ){
      h_ll_numJet20->Fill(numJet20,wgt);
      h_ll_numBtag20->Fill(numBtag20,wgt);

      h_ll_numJet30->Fill(numJet30,wgt);
      h_ll_numBtag30->Fill(numBtag30,wgt);

      h_ll_numJet20_numJet30->Fill(numJet20,numJet30,wgt);
      h_ll_numBtag20_numBtag30->Fill(numBtag20,numBtag30,wgt);
    }
    if( TwoEle ){
      h_ee_numJet20->Fill(numJet20,wgt);
      h_ee_numBtag20->Fill(numBtag20,wgt);

      h_ee_numJet30->Fill(numJet30,wgt);
      h_ee_numBtag30->Fill(numBtag30,wgt);

      h_ee_numJet20_numJet30->Fill(numJet20,numJet30,wgt);
      h_ee_numBtag20_numBtag30->Fill(numBtag20,numBtag30,wgt);
    }
    if( TwoMu  ){
      h_mm_numJet20->Fill(numJet20,wgt);
      h_mm_numBtag20->Fill(numBtag20,wgt);

      h_mm_numJet30->Fill(numJet30,wgt);
      h_mm_numBtag30->Fill(numBtag30,wgt);

      h_mm_numJet20_numJet30->Fill(numJet20,numJet30,wgt);
      h_mm_numBtag20_numBtag30->Fill(numBtag20,numBtag30,wgt);
    }
    if( EleMu  ){
      h_em_numJet20->Fill(numJet20,wgt);
      h_em_numBtag20->Fill(numBtag20,wgt);

      h_em_numJet30->Fill(numJet30,wgt);
      h_em_numBtag30->Fill(numBtag30,wgt);

      h_em_numJet20_numJet30->Fill(numJet20,numJet30,wgt);
      h_em_numBtag20_numBtag30->Fill(numBtag20,numBtag30,wgt);
    }


    if( numJet20==2 ){
      if( TwoLep ) h_ll_event_selection->Fill(0.5+11, wgt); // all
      if( TwoEle ) h_ee_event_selection->Fill(0.5+11, wgt); // all
      if( TwoMu  ) h_mm_event_selection->Fill(0.5+11, wgt); // all
      if( EleMu  ) h_em_event_selection->Fill(0.5+11, wgt); // all

      for( int iJet=0; iJet<int(ind_jet20.size()); iJet++ ){
	double pt  = eve->jet_pt_[iJet];
	double eta = eve->jet_eta_[iJet];
	double phi = eve->jet_phi_[iJet];
	double csv =  eve->jet_csv_[iJet];
	int flavor = eve->jet_flavor_[iJet];

	if( csv < 0.0 ) csv = -0.01;
	if( csv > 1.0 ) csv = 1.0;

	if( flavor==21 ) flavor=7;

	// require tag jet
	if( !(csv>0.89 && pt>20 && fabs(eta)<2.4) ) continue;

	h_e2j20_tagJet_flavor->Fill(flavor,wgt);
	if( abs(flavor)==5 ) h_e2j20_tagJet_flavor->Fill(9,wgt);
	else                 h_e2j20_tagJet_flavor->Fill(10,wgt);

	int numProbeJet20=0;
	int numBProbeJet20=0;
	for( int kJet=0; kJet<int(ind_jet20.size()); kJet++ ){
	  if( kJet==iJet ) continue;

	  double probePt  = eve->jet_pt_[kJet];
	  double probeEta = eve->jet_eta_[kJet];
	  double probePhi = eve->jet_phi_[kJet];
	  double probeCSV =  eve->jet_csv_[kJet];
	  int probeFlavor = eve->jet_flavor_[kJet];

	  if( probeCSV < 0.0 ) probeCSV = -0.01;
	  if( probeCSV > 1.0 ) probeCSV = 1.0;

	  int iPt = -1;
	  if (probePt>=19.99    && probePt<30)  iPt = 0;
	  else if (probePt>=30  && probePt<40)  iPt = 1;
	  else if (probePt>=40  && probePt<60)  iPt = 2;
	  else if (probePt>=60  && probePt<100) iPt = 3;
	  else if (probePt>=100 && probePt<160) iPt = 4;
	  else if (probePt>=160 )               iPt = 5;

	  if( iPt==0 ){
	    numProbeJet20++;
	    if( abs(probeFlavor)==5 ) numBProbeJet20++;
	  }

	  if( probeFlavor==21 ) probeFlavor=7;

	  h_e2j20_probeJet_flavor[iPt]->Fill(probeFlavor,wgt);
	  if( abs(probeFlavor)==5 ) h_e2j20_probeJet_flavor[iPt]->Fill(9,wgt);
	  else                      h_e2j20_probeJet_flavor[iPt]->Fill(10,wgt);

	  double deltaR = DeltaR(eta, phi, probeEta, probePhi);
	  double deltaEta = fabs(eta - probeEta);
	  double deltaPhi = DeltaPhi(phi, probePhi);

	  if( abs(probeFlavor)==5 ){
	    h_e2j20_b_probeJet_csv[iPt]->Fill(probeCSV,wgt);

	    h_e2j20_b_deltaR_probeJet_tagJet[iPt]->Fill(deltaR,wgt);
	    h_e2j20_b_deltaEta_probeJet_tagJet[iPt]->Fill(deltaEta,wgt);
	    h_e2j20_b_deltaPhi_probeJet_tagJet[iPt]->Fill(deltaPhi,wgt);
	  }
	  else {
	    h_e2j20_nonb_probeJet_csv[iPt]->Fill(probeCSV,wgt);

	    h_e2j20_nonb_deltaR_probeJet_tagJet[iPt]->Fill(deltaR,wgt);
	    h_e2j20_nonb_deltaEta_probeJet_tagJet[iPt]->Fill(deltaEta,wgt);
	    h_e2j20_nonb_deltaPhi_probeJet_tagJet[iPt]->Fill(deltaPhi,wgt);
	  }
	}
	h_e2j20_numProbeJet20->Fill(numProbeJet20,wgt);
	h_e2j20_numProbeJet20_numBProbeJet20->Fill(numProbeJet20,numBProbeJet20,wgt);
      }
    }

    if( numJet20>=2 && numJet30>=1 && numJet30<=2 ){
      if( TwoLep ) h_ll_event_selection->Fill(0.5+12, wgt); // all
      if( TwoEle ) h_ee_event_selection->Fill(0.5+12, wgt); // all
      if( TwoMu  ) h_mm_event_selection->Fill(0.5+12, wgt); // all
      if( EleMu  ) h_em_event_selection->Fill(0.5+12, wgt); // all

      for( int iJet=0; iJet<int(ind_jet20.size()); iJet++ ){
	double pt  = eve->jet_pt_[iJet];
	double eta = eve->jet_eta_[iJet];
	double phi = eve->jet_phi_[iJet];
	double csv =  eve->jet_csv_[iJet];
	int flavor = eve->jet_flavor_[iJet];

	if( csv < 0.0 ) csv = -0.01;
	if( csv > 1.0 ) csv = 1.0;

	if( flavor==21 ) flavor=7;

	// require tag jet
	if( !(csv>0.89 && pt>30 && fabs(eta)<2.4) ) continue;

	h_j30j20_tagJet_flavor->Fill(flavor,wgt);
	if( abs(flavor)==5 ) h_j30j20_tagJet_flavor->Fill(9,wgt);
	else                 h_j30j20_tagJet_flavor->Fill(10,wgt);

	int numProbeJet20=0;
	int numBProbeJet20=0;
	for( int kJet=0; kJet<int(ind_jet20.size()); kJet++ ){
	  if( kJet==iJet ) continue;

	  double probePt  = eve->jet_pt_[kJet];
	  double probeEta = eve->jet_eta_[kJet];
	  double probePhi = eve->jet_phi_[kJet];
	  double probeCSV =  eve->jet_csv_[kJet];
	  int probeFlavor = eve->jet_flavor_[kJet];

	  if( probeCSV < 0.0 ) probeCSV = -0.01;
	  if( probeCSV > 1.0 ) probeCSV = 1.0;

	  int iPt = -1;
	  if (probePt>=19.99    && probePt<30)  iPt = 0;
	  else if (probePt>=30  && probePt<40)  iPt = 1;
	  else if (probePt>=40  && probePt<60)  iPt = 2;
	  else if (probePt>=60  && probePt<100) iPt = 3;
	  else if (probePt>=100 && probePt<160) iPt = 4;
	  else if (probePt>=160 )               iPt = 5;

	  if( iPt==0 ){
	    numProbeJet20++;
	    if( abs(probeFlavor)==5 ) numBProbeJet20++;
	  }

	  if( probeFlavor==21 ) probeFlavor=7;

	  h_j30j20_probeJet_flavor[iPt]->Fill(probeFlavor,wgt);
	  if( abs(probeFlavor)==5 ) h_j30j20_probeJet_flavor[iPt]->Fill(9,wgt);
	  else                      h_j30j20_probeJet_flavor[iPt]->Fill(10,wgt);

	  double deltaR = DeltaR(eta, phi, probeEta, probePhi);
	  double deltaEta = fabs(eta - probeEta);
	  double deltaPhi = DeltaPhi(phi, probePhi);

	  if( abs(probeFlavor)==5 ){
	    h_j30j20_b_probeJet_csv[iPt]->Fill(probeCSV,wgt);

	    h_j30j20_b_deltaR_probeJet_tagJet[iPt]->Fill(deltaR,wgt);
	    h_j30j20_b_deltaEta_probeJet_tagJet[iPt]->Fill(deltaEta,wgt);
	    h_j30j20_b_deltaPhi_probeJet_tagJet[iPt]->Fill(deltaPhi,wgt);
	  }
	  else {
	    h_j30j20_nonb_probeJet_csv[iPt]->Fill(probeCSV,wgt);

	    h_j30j20_nonb_deltaR_probeJet_tagJet[iPt]->Fill(deltaR,wgt);
	    h_j30j20_nonb_deltaEta_probeJet_tagJet[iPt]->Fill(deltaEta,wgt);
	    h_j30j20_nonb_deltaPhi_probeJet_tagJet[iPt]->Fill(deltaPhi,wgt);
	  }
	}
	h_j30j20_numProbeJet20->Fill(numProbeJet20,wgt);
	h_j30j20_numProbeJet20_numBProbeJet20->Fill(numProbeJet20,numBProbeJet20,wgt);
      }
    }


    if( numJet20>=2 ){
      if( TwoLep ) h_ll_event_selection->Fill(0.5+13, wgt); // all
      if( TwoEle ) h_ee_event_selection->Fill(0.5+13, wgt); // all
      if( TwoMu  ) h_mm_event_selection->Fill(0.5+13, wgt); // all
      if( EleMu  ) h_em_event_selection->Fill(0.5+13, wgt); // all

      for( int iJet=0; iJet<int(ind_jet20.size()); iJet++ ){
	double pt  = eve->jet_pt_[iJet];
	double eta = eve->jet_eta_[iJet];
	double phi = eve->jet_phi_[iJet];
	double csv =  eve->jet_csv_[iJet];
	int flavor = eve->jet_flavor_[iJet];

	if( csv < 0.0 ) csv = -0.01;
	if( csv > 1.0 ) csv = 1.0;

	if( flavor==21 ) flavor=7;

	// require tag jet
	if( !(csv>0.89 && pt>20 && fabs(eta)<2.4) ) continue;

	h_ge2j20_tagJet_flavor->Fill(flavor,wgt);
	if( abs(flavor)==5 ) h_ge2j20_tagJet_flavor->Fill(9,wgt);
	else                 h_ge2j20_tagJet_flavor->Fill(10,wgt);

	int numProbeJet20=0;
	int numBProbeJet20=0;
	std::vector<double> probe_jet_pt;
	std::vector<int> probe_jet_ipt;
	std::vector<int> probe_jet_flavor;
	std::vector<double> probe_jet_csv;
	for( int kJet=0; kJet<int(ind_jet20.size()); kJet++ ){
	  if( kJet==iJet ) continue;

	  double probePt  = eve->jet_pt_[kJet];
	  double probeEta = eve->jet_eta_[kJet];
	  double probePhi = eve->jet_phi_[kJet];
	  double probeCSV =  eve->jet_csv_[kJet];
	  int probeFlavor = eve->jet_flavor_[kJet];

	  if( probeCSV < 0.0 ) probeCSV = -0.01;
	  if( probeCSV > 1.0 ) probeCSV = 1.0;

	  int iPt = -1;
	  if (probePt>=19.99    && probePt<30)  iPt = 0;
	  else if (probePt>=30  && probePt<40)  iPt = 1;
	  else if (probePt>=40  && probePt<60)  iPt = 2;
	  else if (probePt>=60  && probePt<100) iPt = 3;
	  else if (probePt>=100 && probePt<160) iPt = 4;
	  else if (probePt>=160 )               iPt = 5;

	  probe_jet_pt.push_back(probePt);
	  probe_jet_ipt.push_back(iPt);
	  probe_jet_flavor.push_back(probeFlavor);
	  probe_jet_csv.push_back(probeCSV);

	  if( iPt==0 ){
	    numProbeJet20++;
	    if( abs(probeFlavor)==5 ) numBProbeJet20++;
	  }

	  if( probeFlavor==21 ) probeFlavor=7;

	  h_ge2j20_probeJet_flavor[iPt]->Fill(probeFlavor,wgt);
	  if( abs(probeFlavor)==5 ) h_ge2j20_probeJet_flavor[iPt]->Fill(9,wgt);
	  else                      h_ge2j20_probeJet_flavor[iPt]->Fill(10,wgt);

	  double deltaR = DeltaR(eta, phi, probeEta, probePhi);
	  double deltaEta = fabs(eta - probeEta);
	  double deltaPhi = DeltaPhi(phi, probePhi);

	  if( abs(probeFlavor)==5 ){
	    h_ge2j20_b_probeJet_csv[iPt]->Fill(probeCSV,wgt);

	    h_ge2j20_b_deltaR_probeJet_tagJet[iPt]->Fill(deltaR,wgt);
	    h_ge2j20_b_deltaEta_probeJet_tagJet[iPt]->Fill(deltaEta,wgt);
	    h_ge2j20_b_deltaPhi_probeJet_tagJet[iPt]->Fill(deltaPhi,wgt);
	  }
	  else {
	    h_ge2j20_nonb_probeJet_csv[iPt]->Fill(probeCSV,wgt);

	    h_ge2j20_nonb_deltaR_probeJet_tagJet[iPt]->Fill(deltaR,wgt);
	    h_ge2j20_nonb_deltaEta_probeJet_tagJet[iPt]->Fill(deltaEta,wgt);
	    h_ge2j20_nonb_deltaPhi_probeJet_tagJet[iPt]->Fill(deltaPhi,wgt);
	  }
	}
	h_ge2j20_numProbeJet20->Fill(numProbeJet20,wgt);
	h_ge2j20_numProbeJet20_numBProbeJet20->Fill(numProbeJet20,numBProbeJet20,wgt);

	int n = int(probe_jet_pt.size());

	Int_t idx[n];
	double pt_sorted[n];
	for( int j=0; j<n; j++ ) pt_sorted[j] = probe_jet_pt[j];

	TMath::Sort(n,pt_sorted,idx,true);

	if( n>0 ){
	  double pt1  = probe_jet_pt[idx[0]];
	  int ipt1 = probe_jet_ipt[idx[0]];
	  int flavor1 = probe_jet_flavor[idx[0]];
	  double csv1 = probe_jet_csv[idx[0]];

	  if( !(ipt1<0) ){
	    h_ge2j20_probeJet_firstJet_flavor[ipt1]->Fill(flavor1,wgt);
	    if( abs(flavor1)==5 ) h_ge2j20_probeJet_firstJet_flavor[ipt1]->Fill(9,wgt);
	    else                  h_ge2j20_probeJet_firstJet_flavor[ipt1]->Fill(10,wgt);

	    if( abs(flavor1)==5 ) h_ge2j20_b_probeJet_firstJet_csv[ipt1]->Fill(csv1,wgt);
	    else                  h_ge2j20_nonb_probeJet_firstJet_csv[ipt1]->Fill(csv1,wgt);
	  }

	  if( flavor1==21 ) flavor1=7;
	  if( pt1>19.999 && pt1<30 ){
	    h_ge2j20_probeJet_jet1_flavor->Fill(flavor1,wgt);
	    if( abs(flavor1)==5 ) h_ge2j20_probeJet_jet1_flavor->Fill(9,wgt);
	    else                  h_ge2j20_probeJet_jet1_flavor->Fill(10,wgt);
	  }
	}
      }
    }

    if( numJet20>=1 && numJet30>=1 && (numJet20>=2 || numJet30>=2) ){
      if( TwoLep ) h_ll_event_selection->Fill(0.5+14, wgt); // all
      if( TwoEle ) h_ee_event_selection->Fill(0.5+14, wgt); // all
      if( TwoMu  ) h_mm_event_selection->Fill(0.5+14, wgt); // all
      if( EleMu  ) h_em_event_selection->Fill(0.5+14, wgt); // all
    }


  } // end loop over events

  std::cout << "**************************************************************" << std::endl;
  std::cout << "\t Number of raw events = " << numEvents_all << std::endl;
  std::cout << "\t Number of gen weighted events = " << numEvents_wgt_gen << std::endl;
  std::cout << "\t Number of lumi weighted events = " << numEvents_wgt_lumi << std::endl;
  std::cout << "\t Number of gen * lumi weighted events = " << numEvents_wgt << std::endl;
  std::cout << "**************************************************************" << std::endl;



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

float DeltaPhi(float phi1,float phi2){
  float deltaPhi = TMath::Abs(phi1-phi2);
  if(deltaPhi > TMath::Pi())
    deltaPhi = TMath::TwoPi() - deltaPhi;
  return deltaPhi;
}
