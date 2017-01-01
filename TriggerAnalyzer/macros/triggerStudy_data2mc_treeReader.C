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
//#include "TriggerRun2/TriggerAnalyzer/interface/BTagCalibrationStandalone.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"

#endif


//*****************************************************************************
typedef std::vector< TLorentzVector >          vecTLorentzVector;
typedef std::vector<int>                       vint;
typedef std::vector<double>                    vdouble;
typedef std::vector<std::vector<double> >      vvdouble;

double reweightPU( int nPU, int iSys );
float DeltaR(float eta1,float phi1,float eta2,float phi2);

// ------------ csv applying functions -------------
void fillCSVhistos(TFile *fileHF, TFile *fileLF);
double get_csv_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCSV, vint jetFlavor, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );

double get_btv_csv_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCSV, vint jetFlavor, TString sys_name, bool use_csv, bool verbose );

void fillCSVEffhistos(TFile *file);
double get_csv_efficiency( double jetPt, double jetEta, int jetFlavor );


int PtBinsHF_ = 5;

// CSV reweighting
TH1D* h_csv_wgt_hf[9][5];
TH1D* c_csv_wgt_hf[9][5];
TH1D* h_csv_wgt_lf[9][4][3];

// BTV efficiency
TH2D* h_a_jet_pt_eta_eff_ = NULL;
TH2D* h_b_jet_pt_eta_eff_ = NULL;
TH2D* h_c_jet_pt_eta_eff_ = NULL;
TH2D* h_l_jet_pt_eta_eff_ = NULL;

TH2D* h_a_jet_pt_eta_all_ = NULL;
TH2D* h_b_jet_pt_eta_all_ = NULL;
TH2D* h_c_jet_pt_eta_all_ = NULL;
TH2D* h_l_jet_pt_eta_all_ = NULL;


// setup calibration readers
BTagCalibration calib_btv_csvv2("csvv2", "CSVv2.csv");
BTagCalibrationReader reader_btv_mujets(&calib_btv_csvv2, // calibration instance
					BTagEntry::OP_MEDIUM, // operating point
					"mujets", // measurement type
					"central"); // systematics type

BTagCalibrationReader reader_btv_mujets_up(&calib_btv_csvv2, // calibration instance
					BTagEntry::OP_MEDIUM, // operating point
					"mujets", // measurement type
					"up"); // systematics type

BTagCalibrationReader reader_btv_mujets_down(&calib_btv_csvv2, // calibration instance
					     BTagEntry::OP_MEDIUM, // operating point
					     "mujets", // measurement type
					     "down"); // systematics type


BTagCalibrationReader reader_btv_comb(&calib_btv_csvv2, // calibration instance
					BTagEntry::OP_MEDIUM, // operating point
					"comb", // measurement type
					"central"); // systematics type

BTagCalibrationReader reader_btv_comb_up(&calib_btv_csvv2, // calibration instance
					BTagEntry::OP_MEDIUM, // operating point
					"comb", // measurement type
					"up"); // systematics type

BTagCalibrationReader reader_btv_comb_down(&calib_btv_csvv2, // calibration instance
					   BTagEntry::OP_MEDIUM, // operating point
					   "comb", // measurement type
					   "down"); // systematics type


//*****************************************************************************

void triggerStudy_data2mc_treeReader( int insample=1, int maxNentries=-1, int Njobs=1, int jobN=1, double intLumi=-1, int ttCat_=-1, bool useHTbins_=false, bool useCondor_=false ) {

  // set seed for reproducibility
  TRandom3 r(12345);

  std::string inputFileHF = "data/csv_rwt_hf_IT_FlatSF_2015_07_27.root";
  std::string inputFileLF = "data/csv_rwt_lf_IT_FlatSF_2015_07_27.root";

  TFile* f_CSVwgt_HF = new TFile("/uscms_data/d2/dpuigh/TTH/triggerRun2/CMSSW_7_4_15/src/csvReweightingRun2/csvTreeMaker/data/csv_rwt_fit_hf_2015_12_14.root");
  TFile* f_CSVwgt_LF = new TFile("/uscms_data/d2/dpuigh/TTH/triggerRun2/CMSSW_7_4_15/src/csvReweightingRun2/csvTreeMaker/data/csv_rwt_fit_lf_2015_12_14.root");

  fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);


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
    mySample_nGen_ = 19757190;//19757190//116591749//25357774;//25446993;
    mySample_sampleName_ = "ttbar";
    // standard
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164524/0000/");
    // ext
    //mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2_ext3-v1_triggerTree_v1/151115_164541/0000/");
    if( ttCat_>=0 ){
      if( ttCat_==0 ) mySample_sampleName_ = "ttlf";
      if( ttCat_==1 ) mySample_sampleName_ = "ttcc";
      if( ttCat_==2 ) mySample_sampleName_ = "ttb";
      if( ttCat_==3 ) mySample_sampleName_ = "tt2b";
      if( ttCat_==4 ) mySample_sampleName_ = "ttbb";
    }
  }
  else if( insample==2501 ){
    mySample_xSec_ = 831.76;//https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
    mySample_nGen_ = 11158755;//25357774;//25446993;
    mySample_sampleName_ = "ttbar_MGMLM";
    // standard
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151118_182129/0000/");
    if( ttCat_>=0 ){
      if( ttCat_==0 ) mySample_sampleName_ = "ttlf_MGMLM";
      if( ttCat_==1 ) mySample_sampleName_ = "ttcc_MGMLM";
      if( ttCat_==2 ) mySample_sampleName_ = "ttb_MGMLM";
      if( ttCat_==3 ) mySample_sampleName_ = "tt2b_MGMLM";
      if( ttCat_==4 ) mySample_sampleName_ = "ttbb_MGMLM";
    }
  }
  else if( insample==2502 ){
    mySample_xSec_ = 831.76;//https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
    mySample_nGen_ = 14188545;//42784971//25357774;//25446993;
    mySample_sampleName_ = "ttbar_aMC";
    // standard
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151118_181902/0000/");
    if( ttCat_>=0 ){
      if( ttCat_==0 ) mySample_sampleName_ = "ttlf_aMC";
      if( ttCat_==1 ) mySample_sampleName_ = "ttcc_aMC";
      if( ttCat_==2 ) mySample_sampleName_ = "ttb_aMC";
      if( ttCat_==3 ) mySample_sampleName_ = "tt2b_aMC";
      if( ttCat_==4 ) mySample_sampleName_ = "ttbb_aMC";
    }
  }
  else if( insample==2300 ){
    mySample_xSec_ = 6025.2; 
    mySample_nGen_ = 19310834;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M50";
    if( useHTbins_ ){
      mySample_xSec_ = 6025.2 * 9.59938e-01; // FIXME
      mySample_nGen_ = 19310834 * 9.59938e-01;// FIXME
      mySample_sampleName_ = "ZJets_M50_HT0to100";
    }
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164558/0000/");
  }
  else if( insample==2301 ){
    mySample_xSec_ = 6025.2 * 3.02711e-02;//139.4 * 1.23; 
    mySample_nGen_ = 2725655;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M50_HT100to200";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_165424/0000/");
  }
  else if( insample==2302 ){
    mySample_xSec_ = 6025.2 * 8.26415e-03;//42.75 * 1.23; 
    mySample_nGen_ = 973937;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M50_HT200to400";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_165441/0000/");
  }
  else if( insample==2303 ){
    mySample_xSec_ = 6025.2 * 1.12585e-03;//5.497 * 1.23; 
    mySample_nGen_ = 1067758;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M50_HT400to600";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_165459/0000/");
  }
  else if( insample==2304 ){
    mySample_xSec_ = 6025.2 * 4.01406e-04;//2.21 * 1.23; 
    mySample_nGen_ = 998912;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M50_HT600toInf";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_165518/0000/");
  }
  else if( insample==2305 ){
    mySample_xSec_ = 6025.2 * 9.59938e-01; 
    mySample_nGen_ = 9042031;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M50_MGMLM";
    if( useHTbins_ ){
      mySample_xSec_ = 6025.2 * 9.59938e-01; // FIXME
      mySample_nGen_ = 9042031 * 9.59938e-01;// FIXME
      mySample_sampleName_ = "ZJets_M50_MGMLM_HT0to100";
    }
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151119_132217/0000/");
  }
  else if( insample==2310 ){
    mySample_xSec_ = 18610;//22635.09; 
    mySample_nGen_ = 22217467;//19925500;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M10to50";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164616/0000/");
  }
  else if( insample==2400 ){
    mySample_xSec_ = 61526.7;  
    mySample_nGen_ = 16518218;
    mySample_sampleName_ = "WJets";
    if( useHTbins_ ){
      mySample_xSec_ = 61526.7 * 9.648036e-01; // FIXME
      mySample_nGen_ = 16518218 * 9.648036e-01;// FIXME
      mySample_sampleName_ = "WJets_HT0To100";
    }
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164715/0000/");
  }
  else if( insample==2401 ){
    mySample_xSec_ = 1345 * 1.21;  
    mySample_nGen_ = 10152718;
    mySample_sampleName_ = "WJets_HT100To200";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_165300/0000/");
  }
  else if( insample==2402 ){
    mySample_xSec_ = 359.7 * 1.21;  
    mySample_nGen_ = 5221599;
    mySample_sampleName_ = "WJets_HT200To400";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_165317/0000/");
  }
  else if( insample==2403 ){
    mySample_xSec_ = 48.91 * 1.21;  
    mySample_nGen_ = 1745914;
    mySample_sampleName_ = "WJets_HT400To600";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_165336/0000/");
  }
  else if( insample==2404 ){
    mySample_xSec_ = 18.77 * 1.21;  
    mySample_nGen_ = 1039152;
    mySample_sampleName_ = "WJets_HT600ToInf";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_165404/0000/");
  }
  else if( insample==2405 ){
    mySample_xSec_ = 61526.7;  
    mySample_nGen_ = 72207128;
    mySample_sampleName_ = "WJets_MGMLM";
    if( useHTbins_ ){
      mySample_xSec_ = 61526.7 * 9.648036e-01; // FIXME
      mySample_nGen_ = 72207128 * 9.648036e-01;// FIXME
      mySample_sampleName_ = "WJets_MGMLM_HT0To100";
    }
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151119_132138/0000/");
  }
  else if( insample==2524 ){
    mySample_xSec_ = 0.435;  
    mySample_nGen_ = 430330;
    mySample_sampleName_ = "ttW_had";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_165130/0000/");
  }
  else if( insample==2525 ){
    mySample_xSec_ = 0.21;  
    mySample_nGen_ = 129850;
    mySample_sampleName_ = "ttW_lep";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_165147/0000/");
  }
  else if( insample==2523 ){
    mySample_xSec_ = 0.611;  
    mySample_nGen_ = 351398;
    mySample_sampleName_ = "ttZ_had";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164913/0000/");
  }
  else if( insample==2522 ){
    mySample_xSec_ = 0.2529;  
    mySample_nGen_ = 184990;
    mySample_sampleName_ = "ttZ_lep";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164931/0000/");
  }
  else if( insample==2510 ){
    mySample_xSec_ = 136.02 * 1./3.; 
    mySample_nGen_ = 3299800;
    mySample_sampleName_ = "st_tchan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164755/0000/");
  }
  else if( insample==2511 ){
    mySample_xSec_ = 80.95 * 1./3.;  
    mySample_nGen_ = 1695400;
    mySample_sampleName_ = "stbar_tchan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164733/0000/");
  }
  else if( insample==2512 ){
    mySample_xSec_ = 35.85;  
    mySample_nGen_ = 995600;
    mySample_sampleName_ = "st_tWchan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164817/0000/");
  }
  else if( insample==2513 ){
    mySample_xSec_ = 35.85;  
    mySample_nGen_ = 1000000;
    mySample_sampleName_ = "stbar_tWchan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164835/0000/");
  }
  else if( insample==2514 ){
    mySample_xSec_ = 10.32 * 1./3.;  
    mySample_nGen_ = 613384;
    mySample_sampleName_ = "st_schan";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164856/0000/");
  }
  else if( insample==2424 ){
    mySample_xSec_ = 118.7;  
    mySample_nGen_ = 994416;
    mySample_sampleName_ = "WW";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WW_TuneCUETP8M1_13TeV-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_165204/0000/");
  }
  else if( insample==2423 ){
    mySample_xSec_ = 44.9;  
    mySample_nGen_ = 991232;
    mySample_sampleName_ = "WZ";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_165221/0000/");
  }
  else if( insample==2323 ){
    mySample_xSec_ = 15.4;  
    mySample_nGen_ = 996168;
    mySample_sampleName_ = "ZZ";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_165240/0000/");
  }
  else if( insample==1001 ){
    mySample_xSec_ = 27850000;  
    mySample_nGen_ = 81637494;
    mySample_sampleName_ = "QCD_HT100to200";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151210_174712/0000/");
  }
  else if( insample==1002 ){
    mySample_xSec_ = 1717000;  
    mySample_nGen_ = 18718905;
    mySample_sampleName_ = "QCD_HT200to300";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151210_175438/0000/");
  }
  else if( insample==1003 ){
    mySample_xSec_ = 351300;  
    mySample_nGen_ = 19826197;
    mySample_sampleName_ = "QCD_HT300to500";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151210_175811/0000/");
  }
  else if( insample==1004 ){
    mySample_xSec_ = 31630;  
    mySample_nGen_ = 19664159;
    mySample_sampleName_ = "QCD_HT500to700";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151210_182545/0000/");
  }
  else if( insample==1005 ){
    mySample_xSec_ = 6802;  
    mySample_nGen_ = 15356448;
    mySample_sampleName_ = "QCD_HT700to1000";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151210_182934/0000/");
  }
  else if( insample==1006 ){
    mySample_xSec_ = 1206;  
    mySample_nGen_ = 4963895;
    mySample_sampleName_ = "QCD_HT1000to1500";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151210_183402/0000/");
  }
  else if( insample==1007 ){
    mySample_xSec_ = 120.4;  
    mySample_nGen_ = 3868886;
    mySample_sampleName_ = "QCD_HT1500to2000";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151210_185729/0000/");
  }
  else if( insample==1008 ){
    mySample_xSec_ = 25.24;  
    mySample_nGen_ = 1912529;
    mySample_sampleName_ = "QCD_HT2000toInf";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151210_191142/0000/");
  }
  else if( insample==9125 ){
    mySample_xSec_ = 0.2934;// YR3 * BR(all)  
    mySample_nGen_ = 3933404;//199000;
    mySample_sampleName_ = "ttHTobb";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ttHTobb_M125_13TeV_powheg_pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164638/0000/");
  }
  else if( insample==8125 ){
    mySample_xSec_ = 0.2151;// YR3 * BR(all)  
    mySample_nGen_ = 3796398;//199000;
    mySample_sampleName_ = "ttHnonbb";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/ttHToNonbb_M125_13TeV_powheg_pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164658/0000/");
  }
  else if( insample==-13 ){
    mySample_sampleName_ = "SingleElectron_Run2015D_05Oct2015_PromptRecov4";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleElectron/Run2015D-05Oct2015-v1_SilverJSON_2015_11_13_triggerTree_v1/151115_164410/0000/");
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleElectron/Run2015D-PromptReco-v4_SilverJSON_2015_11_13_triggerTree_v1/151115_164446/0000/");
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleElectron/Run2015D-PromptReco-v4_SilverJSON_2015_11_13_triggerTree_v1/151115_164446/0001/");
  }
  else if( insample==-11 ){
    mySample_sampleName_ = "SingleMuon_Run2015D_05Oct2015_PromptRecov4";
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleMuon/Run2015D-05Oct2015-v1_SilverJSON_2015_11_13_triggerTree_v1/151115_164428/0000/");
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleMuon/Run2015D-PromptReco-v4_SilverJSON_2015_11_13_triggerTree_v1/151115_164504/0000/");
    mySample_inputDirs_.push_back("/eos/uscms/store/user/puigh/SingleMuon/Run2015D-PromptReco-v4_SilverJSON_2015_11_13_triggerTree_v1/151115_164504/0001/");
  }



  std::string s_end = "histo_" + str_jobN + ".root";
  if( Njobs==1 ) s_end = "histo.root";

  std::string histofilename = Form("HistoFiles/triggerStudy_data2mc_treeReader_%s_%s", mySample_sampleName_.c_str(), s_end.c_str());
  if( useCondor_ ) histofilename = Form("triggerStudy_data2mc_treeReader_%s_%s", mySample_sampleName_.c_str(), s_end.c_str());


  TChain *chain = new TChain("triggeranalzyer/triggerTree");
  for( int iFile=0; iFile<int(mySample_inputDirs_.size()); iFile++ ){
    std::string treefilename = mySample_inputDirs_[iFile] + "trigger_analyzer*.root";
    std::cout << "  treefilename " << iFile << ": " << treefilename.c_str() << std::endl;
    chain->Add(treefilename.c_str());
  }

  //if( insample==2500 ) chain->Add("/eos/uscms/store/user/puigh/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2_ext3-v1_triggerTree_v1/151115_164541/0000/trigger_analyzer_22*.root");
  //if( insample==2300 ) chain->Add("/eos/uscms/store/user/puigh/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164558/0000/trigger_analyzer_22*.root");

  //chain->Add("/eos/uscms/store/user/puigh/SingleElectron/Run2015D-PromptReco-v4_SilverJSON_2015_11_13_triggerTree_v1/151115_164446/0000/trigger_analyzer_22*.root");

   //chain->Add("/uscms_data/d2/dpuigh/TTH/triggerRun2/CMSSW_7_4_15/src/TriggerRun2/TriggerAnalyzer/trigger_analyzer.root");

  std::string csvefffilename = Form("/uscms_data/d2/dpuigh/TTH/triggerRun2/CMSSW_7_4_15/src/TriggerRun2/TriggerAnalyzer/backup_2015_12_12_HistoFiles/ttHbb_data2mc_treeReader_%s_histo.root", mySample_sampleName_.c_str());

  // if( insample>=1001 && insample<=1008 ) csvefffilename = "/uscms_data/d2/dpuigh/TTH/triggerRun2/CMSSW_7_4_15/src/TriggerRun2/TriggerAnalyzer/backup_2015_12_05_HistoFiles/ttHbb_data2mc_treeReader_ttbar_histo_eff.root";

  std::cout << "  csv eff filename = " << csvefffilename.c_str() << std::endl;

  TFile* f_CSVeff = new TFile(csvefffilename.c_str());

  fillCSVEffhistos(f_CSVeff);


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

  TH1::SetDefaultSumw2();

  bool verbose_ = !true;

  //////////////////////////////////////////////////////////////////////////
  ///  Histograms
  //////////////////////////////////////////////////////////////////////////

  double window = 10.;
  double pdgZmass = 91.1876;
  double MinMass  = pdgZmass - window;
  double MaxMass  = pdgZmass + window;

  double MaxLHEHT = 2500;

  int MaxNjet = 12;
  int MaxNbtag = 5;


  /// General Plots
  TH1D* h_lheHT = new TH1D("h_lheHT",";LHE H_{T}", 2500, 0, MaxLHEHT );
  TH1D* h_lheHT_wgt = new TH1D("h_lheHT_wgt",";LHE H_{T}", 2500, 0, MaxLHEHT );
  TH1D* h_lheHT_wgt_afterHT = new TH1D("h_lheHT_wgt_afterHT",";LHE H_{T}", 2500, 0, MaxLHEHT );


  TH1D* h_numEvents = new TH1D("h_numEvents",";Number of events", 2, 0, 2 );
  TH1D* h_numEvents_wgt = new TH1D("h_numEvents_wgt",";Predicted number of events", 6, 0, 6 );

  TH1D* h_numTruePVs = new TH1D("h_numTruePVs",";Number of True PVs", 50, 0, 50 );

  TH1D* h_numPVs = new TH1D("h_numPVs",";Number of PVs", 50, 0-0.5, 50-0.5 );

  TH1D* h_additionalJetEventId = new TH1D("h_additionalJetEventId",";additionalJetEventId", 201, -100-0.5, 101-0.5 );

  TH1D* h_numPVs_wgt = new TH1D("h_numPVs_wgt",";Number of PVs", 50, 0-0.5, 50-0.5 );
  TH1D* h_numPVs_noPUwgt = new TH1D("h_numPVs_noPUwgt",";Number of PVs", 50, 0-0.5, 50-0.5 );


  /// Dilepton control region

  // double HTmax = 1000.;
  // int numHTbins = 1000;
  double L1HTmax = 1000.;
  int numL1HTTbins = 1000;


  double HTmax = 2000.;
  std::vector<double> ht_xbins;
  for( int iBin=0; iBin<500; iBin+=10 ) ht_xbins.push_back(double(iBin));
  for( int iBin=500; iBin<1000; iBin+=50 ) ht_xbins.push_back(double(iBin));
  for( int iBin=1000; iBin<2000; iBin+=250 ) ht_xbins.push_back(double(iBin));
  ht_xbins.push_back(2000);
  int Nhtbins = int(ht_xbins.size());

  double htbins[65];
  for( int iBin=0; iBin<Nhtbins; iBin++ ) htbins[iBin] = ht_xbins[iBin];

  Nhtbins -= 1;


  // int Nwidehtbins = Nhtbins + 1;
  // double widehtbins[66];
  // for( int iBin=0; iBin<Nhtbins; iBin++ ) widehtbins[iBin] = htbins[iBin];
  // widehtbins[65] = 5000;

  // Nwidehtbins -= 1;


  std::vector<double> ele_ht_xbins;
  for( int iBin=0; iBin<500; iBin+=50 ) ele_ht_xbins.push_back(double(iBin));
  for( int iBin=500; iBin<1000; iBin+=250 ) ele_ht_xbins.push_back(double(iBin));
  for( int iBin=1000; iBin<2000; iBin+=500 ) ele_ht_xbins.push_back(double(iBin));
  ele_ht_xbins.push_back(2000);
  int Nelehtbins = int(ele_ht_xbins.size());

  double elehtbins[15];
  for( int iBin=0; iBin<Nelehtbins; iBin++ ) elehtbins[iBin] = ele_ht_xbins[iBin];

  Nelehtbins -= 1;


  TH1D* h_ll_diLepMass = new TH1D("h_ll_diLepMass",";M(lep,lep)", 100, 50, 150 );
  TH1D* h_ee_diLepMass = new TH1D("h_ee_diLepMass",";M(e,e)", 100, 50, 150 );
  TH1D* h_mm_diLepMass = new TH1D("h_mm_diLepMass",";M(#mu,#mu)", 100, 50, 150 );


  // int Nptbins = 17;
  // double ptbins[] = { 10, 15, 20, 30, 40, 50, 60, 70, 80, 100, 125, 150, 200, 250, 300, 400, 500, 600 };
  // int Nptbins = 17;
  // double ptbins[] = { 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 125, 150, 200, 300, 500 };
  int Nptbins = 13;
  double ptbins[] = { 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 500 };

  // int Netabins = 12;
  // double etabins[] = { -2.5, -2.1, -1.65, -1.2, -0.9, -0.45, 0, 0.45, 0.9, 1.2, 1.65, 2.1, 2.5 };
  int Netabins = 10;
  double etabins[] = { -2.5, -2.1, -1.566, -1.4442, -0.8, 0, 0.8, 1.4442, 1.566, 2.1, 2.5 };

  int Nphibins = 16;
  double phibins[] = { -3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2 };

  int Npvbins = 14;
  double pvbins[] = { 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17, 35 };

  TH1D* h_probe_pt_pHLT_fL1T = new TH1D("h_probe_pt_pHLT_fL1T",";probe p_{T}", Nptbins, ptbins );

  TH1D* h_probe_pt_pHLTall = new TH1D("h_probe_pt_pHLTall",";probe p_{T}", Nptbins, ptbins );

  TH1D* h_probe_pt_pL1T = new TH1D("h_probe_pt_pL1T",";probe p_{T}", Nptbins, ptbins );
  TH1D* h_probe_pt_fL1T = new TH1D("h_probe_pt_fL1T",";probe p_{T}", Nptbins, ptbins );
  TH1D* h_probe_pt_pHLT = new TH1D("h_probe_pt_pHLT",";probe p_{T}", Nptbins, ptbins );
  TH1D* h_probe_pt_fHLT = new TH1D("h_probe_pt_fHLT",";probe p_{T}", Nptbins, ptbins );
  TH1D* h_probe_pt_all  = new TH1D("h_probe_pt_all", ";probe p_{T}", Nptbins, ptbins );

  TH1D* h_probe_eta_pL1T = new TH1D("h_probe_eta_pL1T",";probe #eta", Netabins, etabins );
  TH1D* h_probe_eta_fL1T = new TH1D("h_probe_eta_fL1T",";probe #eta", Netabins, etabins );
  TH1D* h_probe_eta_pHLT = new TH1D("h_probe_eta_pHLT",";probe #eta", Netabins, etabins );
  TH1D* h_probe_eta_fHLT = new TH1D("h_probe_eta_fHLT",";probe #eta", Netabins, etabins );
  TH1D* h_probe_eta_all  = new TH1D("h_probe_eta_all", ";probe #eta", Netabins, etabins );

  TH1D* h_probe_phi_pL1T = new TH1D("h_probe_phi_pL1T",";probe #phi", Nphibins, phibins );
  TH1D* h_probe_phi_fL1T = new TH1D("h_probe_phi_fL1T",";probe #phi", Nphibins, phibins );
  TH1D* h_probe_phi_pHLT = new TH1D("h_probe_phi_pHLT",";probe #phi", Nphibins, phibins );
  TH1D* h_probe_phi_fHLT = new TH1D("h_probe_phi_fHLT",";probe #phi", Nphibins, phibins );
  TH1D* h_probe_phi_all  = new TH1D("h_probe_phi_all", ";probe #phi", Nphibins, phibins );

  TH1D* h_probe_numPVs_pL1T = new TH1D("h_probe_numPVs_pL1T",";Number of PVs", Npvbins, pvbins );
  TH1D* h_probe_numPVs_fL1T = new TH1D("h_probe_numPVs_fL1T",";Number of PVs", Npvbins, pvbins );
  TH1D* h_probe_numPVs_pHLT = new TH1D("h_probe_numPVs_pHLT",";Number of PVs", Npvbins, pvbins );
  TH1D* h_probe_numPVs_fHLT = new TH1D("h_probe_numPVs_fHLT",";Number of PVs", Npvbins, pvbins );
  TH1D* h_probe_numPVs_all  = new TH1D("h_probe_numPVs_all", ";Number of PVs", Npvbins, pvbins );

  TH1D* h_probe_HT30_pL1T = new TH1D("h_probe_HT30_pL1T",";reco H_{T}", Nelehtbins, elehtbins );
  TH1D* h_probe_HT30_fL1T = new TH1D("h_probe_HT30_fL1T",";reco H_{T}", Nelehtbins, elehtbins );
  TH1D* h_probe_HT30_pHLT = new TH1D("h_probe_HT30_pHLT",";reco H_{T}", Nelehtbins, elehtbins );
  TH1D* h_probe_HT30_fHLT = new TH1D("h_probe_HT30_fHLT",";reco H_{T}", Nelehtbins, elehtbins );
  TH1D* h_probe_HT30_all  = new TH1D("h_probe_HT30_all", ";reco H_{T}", Nelehtbins, elehtbins );



  TH2D* h_probe_pt_eta_pHLTall = new TH2D("h_probe_pt_eta_pHLTall",";probe p_{T};probe #eta", Nptbins, ptbins, Netabins, etabins );
  TH2D* h_probe_pt_eta_pHLT_fL1T = new TH2D("h_probe_pt_eta_pHLT_fL1T",";probe p_{T};probe #eta", Nptbins, ptbins, Netabins, etabins );

  TH2D* h_probe_pt_eta_pL1T = new TH2D("h_probe_pt_eta_pL1T",";probe p_{T};probe #eta", Nptbins, ptbins, Netabins, etabins );
  TH2D* h_probe_pt_eta_fL1T = new TH2D("h_probe_pt_eta_fL1T",";probe p_{T};probe #eta", Nptbins, ptbins, Netabins, etabins );
  TH2D* h_probe_pt_eta_pHLT = new TH2D("h_probe_pt_eta_pHLT",";probe p_{T};probe #eta", Nptbins, ptbins, Netabins, etabins );
  TH2D* h_probe_pt_eta_fHLT = new TH2D("h_probe_pt_eta_fHLT",";probe p_{T};probe #eta", Nptbins, ptbins, Netabins, etabins );
  TH2D* h_probe_pt_eta_all  = new TH2D("h_probe_pt_eta_all", ";probe p_{T};probe #eta", Nptbins, ptbins, Netabins, etabins );


  TH1D* h_2e_L1HTT = new TH1D("h_2e_L1HTT",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_2e_HT30 = new TH1D("h_2e_HT30",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_2e_HT30_L1HTT125 = new TH1D("h_2e_HT30_L1HTT125",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_2e_HT30_L1HTT125_passHLTEle27HT200 = new TH1D("h_2e_HT30_L1HTT125_passHLTEle27HT200",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_2e_HT30_L1HTT100 = new TH1D("h_2e_HT30_L1HTT100",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_2e_HT30_L1HTT100_passHLTEle27HT200 = new TH1D("h_2e_HT30_L1HTT100_passHLTEle27HT200",";reco H_{T}", Nhtbins, htbins );

  TH1D* h_2e_HT30nocc = new TH1D("h_2e_HT30nocc",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_2e_HT30nocc_L1HTT125 = new TH1D("h_2e_HT30nocc_L1HTT125",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_2e_HT30nocc_L1HTT125_passHLTEle27HT200 = new TH1D("h_2e_HT30nocc_L1HTT125_passHLTEle27HT200",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_2e_HT30nocc_L1HTT100 = new TH1D("h_2e_HT30nocc_L1HTT100",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_2e_HT30nocc_L1HTT100_passHLTEle27HT200 = new TH1D("h_2e_HT30nocc_L1HTT100_passHLTEle27HT200",";reco H_{T}", Nhtbins, htbins );

  TH2D* h_2e_HT30nocc_HT30 = new TH2D("h_2e_HT30nocc_HT30",";reco H_{T} (no cc);reco H_{T}", Nhtbins, htbins, Nhtbins, htbins );


  int MaxNjet_2e = 7;
  int MinNjet_2e = 0;
  TH1D* h_2e_numJet = new TH1D("h_2e_numJet",";Number of Jets", MaxNjet_2e - MinNjet_2e, MinNjet_2e - 0.5, MaxNjet_2e - 0.5 );

  int MaxNbtag_2e = 5;
  int MinNbtag_2e = 0;
  TH1D* h_2e_numBtag = new TH1D("h_2e_numBtag",";Number of Btags", MaxNbtag_2e - MinNbtag_2e, MinNbtag_2e - 0.5, MaxNbtag_2e - 0.5 );


  /// lepton + jets control region

  std::vector<TString> cat_labels;
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

  TH1D* h_category_yield    = new TH1D("h_category_yield", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_1e = new TH1D("h_category_yield_1e", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_1m = new TH1D("h_category_yield_1m", ";category", NumCat, 0, NumCat );

  TH1D* h_category_yield_1e_HLT_Ele27 = new TH1D("h_category_yield_1e_HLT_Ele27", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_1e_HLT_Ele23 = new TH1D("h_category_yield_1e_HLT_Ele23", ";category", NumCat, 0, NumCat );

  TH1D* h_category_yield_1e_recoHT300 = new TH1D("h_category_yield_1e_recoHT300", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_1e_recoHT300_L1HTT125 = new TH1D("h_category_yield_1e_recoHT300_L1HTT125", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_1e_recoHT300_L1HTT125_passHLTEle27HT200 = new TH1D("h_category_yield_1e_recoHT300_L1HTT125_passHLTEle27HT200", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_1e_recoHT300_L1HTT100 = new TH1D("h_category_yield_1e_recoHT300_L1HTT100", ";category", NumCat, 0, NumCat );
  TH1D* h_category_yield_1e_recoHT300_L1HTT100_passHLTEle27HT200 = new TH1D("h_category_yield_1e_recoHT300_L1HTT100_passHLTEle27HT200", ";category", NumCat, 0, NumCat );

  for( int iCat=0; iCat<NumCat; iCat++ ){
    h_category_yield->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
    h_category_yield_1e->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
    h_category_yield_1m->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);

    h_category_yield_1e_HLT_Ele27->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
    h_category_yield_1e_HLT_Ele23->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);

    h_category_yield_1e_recoHT300->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
    h_category_yield_1e_recoHT300_L1HTT125->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
    h_category_yield_1e_recoHT300_L1HTT125_passHLTEle27HT200->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
    h_category_yield_1e_recoHT300_L1HTT100->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
    h_category_yield_1e_recoHT300_L1HTT100_passHLTEle27HT200->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
  }


  TH1D* h_1e4j2t_L1HTT = new TH1D("h_1e4j2t_L1HTT",";L1 HTT", numL1HTTbins, 0, L1HTmax );
  TH1D* h_1e4j2t_HT30 = new TH1D("h_1e4j2t_HT30",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_1e4j2t_HT30_L1HTT125 = new TH1D("h_1e4j2t_HT30_L1HTT125",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_1e4j2t_HT30_L1HTT125_passHLTEle27HT200 = new TH1D("h_1e4j2t_HT30_L1HTT125_passHLTEle27HT200",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_1e4j2t_HT30_L1HTT100 = new TH1D("h_1e4j2t_HT30_L1HTT100",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_1e4j2t_HT30_L1HTT100_passHLTEle27HT200 = new TH1D("h_1e4j2t_HT30_L1HTT100_passHLTEle27HT200",";reco H_{T}", Nhtbins, htbins );

  TH1D* h_1e4j2t_HT30nocc = new TH1D("h_1e4j2t_HT30nocc",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_1e4j2t_HT30nocc_L1HTT125 = new TH1D("h_1e4j2t_HT30nocc_L1HTT125",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_1e4j2t_HT30nocc_L1HTT125_passHLTEle27HT200 = new TH1D("h_1e4j2t_HT30nocc_L1HTT125_passHLTEle27HT200",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_1e4j2t_HT30nocc_L1HTT100 = new TH1D("h_1e4j2t_HT30nocc_L1HTT100",";reco H_{T}", Nhtbins, htbins );
  TH1D* h_1e4j2t_HT30nocc_L1HTT100_passHLTEle27HT200 = new TH1D("h_1e4j2t_HT30nocc_L1HTT100_passHLTEle27HT200",";reco H_{T}", Nhtbins, htbins );

  TH2D* h_1e4j2t_HT30nocc_HT30 = new TH2D("h_1e4j2t_HT30nocc_HT30",";reco H_{T} (no cc);reco H_{T}", Nhtbins, htbins, Nhtbins, htbins );

  int MaxNjet_1l4j2t = 11;
  int MinNjet_1l4j2t = 4;
  TH1D* h_1l4j2t_numJet = new TH1D("h_1l4j2t_numJet",";Number of Jets", MaxNjet_1l4j2t - MinNjet_1l4j2t, MinNjet_1l4j2t - 0.5, MaxNjet_1l4j2t - 0.5 );

  int MaxNbtag_1l4j2t = 7;
  int MinNbtag_1l4j2t = 2;
  TH1D* h_1l4j2t_numBtag = new TH1D("h_1l4j2t_numBtag",";Number of Btags", MaxNbtag_1l4j2t - MinNbtag_1l4j2t, MinNbtag_1l4j2t - 0.5, MaxNbtag_1l4j2t - 0.5 );





  //////////////////////////////////////////////////////////////////////////
  /////
  //////////////////////////////////////////////////////////////////////////

  std::string ele_path_name = ( insample<0 ) ? "HLT_Ele27_eta2p1_WPLoose_Gsf_v" : "HLT_Ele27_WP85_Gsf_v";
  std::string mu_path_name  = ( insample<0 ) ? "HLT_IsoMu20_v" : "HLT_IsoMu20_v";

  int numEvents_all = 0;
  int numEvents_wgt_gen = 0;
  double numEvents_wgt_lumi = 0;
  double numEvents_wgt_gen_lumi = 0;
  double numEvents_wgt_gen_lumi_pu = 0;
  double numEvents_wgt_gen_lumi_pu_csv = 0;

  double numEvents_wgt_gen_lumi_pu_4j2t = 0;
  double numEvents_wgt_gen_lumi_pu_csv_4j2t = 0;


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

    // std::cout << "     Event " << ievt << std::endl;

    chain->GetEntry(ievt);

    if( verbose_ ) std::cout << " ===> test 0 " << std::endl;
    // If using HT bins, only take low HT bins of inclusive samples
    double lheHT = eve->lheHT_;
    lheHT = std::min( lheHT, MaxLHEHT-0.00001 );
    
    h_lheHT->Fill(lheHT);


    int additionalJetEventId = -99;
    if( insample>=2500 && insample<=2504 ){
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

      if( !keepTTbarEvent ) continue;
    }

    h_additionalJetEventId->Fill(additionalJetEventId);

    int numTruePVs = eve->numTruePVs_;


    int numPVs = eve->numPVs_;

    //double wgt_pu = ( insample < 0 ) ? 1. : reweightPU(numPVs);
    double wgt_pu = ( insample < 0 ) ? 1. : reweightPU(numTruePVs,0);

    int run = eve->run_;

    //// --------- various weights: PU, topPt, triggerSF, leptonSF...
    // double  wgt_topPtSF = eve->wgt_topPt_;
    double Xsec = mySample_xSec_;//eve->wgt_xs_;
    double nGen = ( maxNentries>0 ) ? maxNentries : mySample_nGen_;//eve->wgt_nGen_;
    if( maxNentries==1000000 && insample==2300 ) nGen = 670121;
    double lumi = ( intLumi > 0 ) ? intLumi : 10000 ;

    double wgt_gen = ( insample > 0 ) ? eve->wgt_generator_ : 1;
    wgt_gen = ( wgt_gen > 0 ) ? 1. : -1.;
    double wgt_lumi = ( insample > 0 ) ? lumi * (Xsec/nGen) : 1;//"weight_PU*topPtWgt*osTriggerSF*lepIDAndIsoSF*"; // various weights

    h_numEvents->Fill(0.5,1.);
    h_numEvents->Fill(1.5,wgt_gen);

    h_numEvents_wgt->Fill(0.5,1.);
    h_numEvents_wgt->Fill(1.5,wgt_gen);
    h_numEvents_wgt->Fill(2.5,wgt_lumi);
    h_numEvents_wgt->Fill(3.5,wgt_gen * wgt_lumi);
    h_numEvents_wgt->Fill(4.5,wgt_gen * wgt_lumi * wgt_pu);


    double wgt = wgt_gen * wgt_lumi * wgt_pu;

    numEvents_all += 1.;
    numEvents_wgt_gen += wgt_gen;
    numEvents_wgt_lumi += wgt_lumi;
    numEvents_wgt_gen_lumi += wgt_gen * wgt_lumi;
    numEvents_wgt_gen_lumi_pu += wgt_gen * wgt_lumi * wgt_pu;


    h_numTruePVs->Fill(numTruePVs,wgt_gen);

    h_lheHT_wgt->Fill(lheHT,wgt_gen * wgt_lumi);


    // Skip HT range, if necessary
    if( useHTbins_ && (insample==2300 || insample==2305 || insample==2400 || insample==2405) ){
      if( lheHT > 100 ) continue;
    }

    h_lheHT_wgt_afterHT->Fill(lheHT,wgt_gen * wgt_lumi);


    ///////////////////
    ////// selections
    ///////////////////

    bool pass_trigger_ele = false;
    if( insample<0 ) pass_trigger_ele = ( eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_==1 );
    else             pass_trigger_ele = ( eve->pass_HLT_Ele27_WP85_Gsf_v_==1 );

    bool pass_trigger_mu  = false;
    if( insample<0 ) pass_trigger_mu = ( eve->pass_HLT_IsoMu20_v_==1 );
    else             pass_trigger_mu = ( eve->pass_HLT_IsoMu20_v_==1 );



    if( insample<0 ){
      pass_trigger_ele = ( insample==-13 && pass_trigger_ele );
      pass_trigger_mu  = ( insample==-11 && pass_trigger_mu  );
    }


    ///// REMOVE trigger requirement
    // Pass Trigger selection
    //if( !(pass_trigger_mu || pass_trigger_ele) ) continue;

    // require at least one PV
    if( numPVs<1 ) continue;


    h_numPVs->Fill(numPVs,wgt);



    h_numPVs_wgt->Fill(numPVs,wgt_gen * wgt_lumi * wgt_pu);

    h_numPVs_noPUwgt->Fill(numPVs,wgt_gen * wgt_lumi);


    //
    // Lepton selection
    //

    std::vector<int> ind_ele;
    std::vector<int> ind_mu;
    std::vector<int> ind_ele_loose;
    std::vector<int> ind_mu_loose;
    for( int iLep=0; iLep<int(eve->lepton_pt_.size()); iLep++ ){
      bool isMuon = ( eve->lepton_isMuon_[iLep]==1 );

      bool isTrigMVAM = ( eve->lepton_isTrigMVAM_[iLep]==1 );
      bool isTight = ( eve->lepton_isTight_[iLep]==1 );
      bool isLoose = ( eve->lepton_isLoose_[iLep]==1 );
      bool isCrack = ( eve->lepton_inCrack_[iLep]==1 );

      double pt  = eve->lepton_pt_[iLep];
      double eta = eve->lepton_eta_[iLep];

      double relIso = eve->lepton_relIso_[iLep];
      double relIsoR04 = eve->lepton_relIsoR04_[iLep];

      ////
      if( isMuon ){
	if( pt > 30 && abs(eta)<2.1 && isTight ) ind_mu.push_back(iLep);
	// FIXME should use 0.25, but need to fix MiniAODHelper
	if( pt > 15 && abs(eta)<2.4 && isLoose && relIsoR04 < 0.20 ) ind_mu_loose.push_back(iLep);
      }

      if( !isMuon ){
	if( pt > 30 && abs(eta)<2.1 && isTrigMVAM && !isCrack && relIso < 0.10 ) ind_ele.push_back(iLep);
	if( pt > 15 && abs(eta)<2.4 && isTrigMVAM && !isCrack && relIso < 0.15 ) ind_ele_loose.push_back(iLep);
      }
    }

    int numEle = int( ind_ele.size() );
    int numMu  = int( ind_mu.size() );

    int numEle_loose = int( ind_ele_loose.size() );
    int numMu_loose  = int( ind_mu_loose.size() );

    bool oneEle = ( numEle==1 && numEle_loose==1 && numMu_loose==0 );
    bool oneMu  = ( numMu==1 && numMu_loose==1 && numEle_loose==0 );

    // oneEle = ( oneEle && pass_trigger_ele );
    // oneMu  = ( oneMu  && pass_trigger_mu );

    bool oneLep = ( oneEle || oneMu );


    bool TwoEle = ( numEle>=1 && numEle_loose==2 && numMu_loose==0 );
    bool TwoMu  = ( numMu>=1 && numMu_loose==2 && numEle_loose==0 );

    TwoEle = ( TwoEle && pass_trigger_ele );
    TwoMu  = ( TwoMu  && pass_trigger_mu );

    bool TwoLep = ( TwoEle || TwoMu );

    // Require at least two leptons
    if( TwoLep ){

      int lepInd1 = -1, lepInd2 = -1;
      if( TwoEle ){
	lepInd1 = ind_ele_loose[0];
	lepInd2 = ind_ele_loose[1];
      }
      else if( TwoMu ){
	lepInd1 = ind_mu_loose[0];
	lepInd2 = ind_mu_loose[1];
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

      TLorentzVector diLep = myLep1 + myLep2;

      double diLepMass = diLep.M();

      //
      // REQUIRE LOOSE DILEPTON MASS
      //

      if( !(diLepMass > 50) ) continue;

      h_ll_diLepMass->Fill(diLepMass,wgt);
      if( TwoEle )     h_ee_diLepMass->Fill(diLepMass,wgt);
      else if( TwoMu ) h_mm_diLepMass->Fill(diLepMass,wgt);


      //
      // REQUIRE TIGHT DILEPTON MASS
      //

      if( !(fabs(diLepMass - 91.2)<10) ) continue;

      // Require two electrons
      if( !TwoEle ) continue;

      bool matchHLT1 = false;
      if( insample<0 ) matchHLT1 = ( eve->lepton_ele_matchHLT_hltEle27WPLooseGsfTrackIsoFilter_[lepInd1]==1 );
      else             matchHLT1 = ( eve->lepton_ele_matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter_[lepInd1]==1 );

      bool matchHLT2 = false;
      if( insample<0 ) matchHLT2 = ( eve->lepton_ele_matchHLT_hltEle27WPLooseGsfTrackIsoFilter_[lepInd2]==1 );
      else             matchHLT2 = ( eve->lepton_ele_matchHLT_hltL1EG25Ele27WP85GsfTrackIsoFilter_[lepInd2]==1 );


      // Must have at least one tag
      if( !(matchHLT1 || matchHLT2 ) ) continue;

      // If more than one tag, pick probe randomly

      bool useFirstLep = ( r.Binomial(1,0.5) < 0.5 );

      int iTag = -1;
      int iProbe = -1;
      if( useFirstLep ){
	iTag = lepInd1;
	iProbe = lepInd2;
      }
      else{
	iTag = lepInd2;
	iProbe = lepInd1;
      }


      /// Jet specific

      /////
      /////  Loop over systematics
      /////

      vdouble use_jet_pt = eve->jet_pt_;
      vdouble use_jet_eta = eve->jet_eta_;
      vdouble use_jet_phi = eve->jet_phi_;
      vdouble use_jet_energy = eve->jet_energy_;

      vdouble use_jet_csv = eve->jet_csv_;
      vdouble use_jet_pileupJetId_fullDiscriminant = eve->jet_pileupJetId_fullDiscriminant_;

      vint use_jet_hadronFlavour = eve->jet_hadronFlavour_;



      // Set lepton scale factor equal to 1 for now
      double wgt_lepIdSF = 1.;

      // PU systematic
      double use_wgt_pu = wgt_pu;

      // Factorization and normalization scale uncertainty
      double wgt_LHEscale = 1.;
      if( insample>-1 ){
	double originalXWGTUP = eve->originalXWGTUP_;
	vdouble lhe_event_weights = eve->LHEEvent_weights_;

	if( lhe_event_weights.size()>6 && originalXWGTUP!=0 ){
	  int use_weight_index = 0;

	  wgt_LHEscale = lhe_event_weights[use_weight_index] / originalXWGTUP;
	}
      }



      vdouble jetPts;
      vdouble jetEtas;
      vdouble jetCSVs;
      vint    jetFlavors;

      vint ind_jet;

      int numJet = 0;
      int numBtag = 0;

      int numJet_b = 0;
      int numJet_c = 0;
      int numJet_l = 0;

      double HT30 = 0;
      TLorentzVector sumJet;
      bool firstJet = false;

      for( int iJet=0; iJet<int(use_jet_pt.size()); iJet++ ){
	TLorentzVector myJet;
	myJet.SetPtEtaPhiE( use_jet_pt[iJet], use_jet_eta[iJet], use_jet_phi[iJet], use_jet_energy[iJet] );

	if( !firstJet ) sumJet = myJet;
	else            sumJet += myJet;

	double pt  = use_jet_pt[iJet];
	double eta = use_jet_eta[iJet];

	if( !(pt>30. && fabs(eta)<2.4) ) continue;

	double dR1 = myJet.DeltaR(myLep1);
	double dR2 = myJet.DeltaR(myLep2);

	if( dR1 < 0.4 ) continue;
	if( dR2 < 0.4 ) continue;

	double csv = use_jet_csv[iJet];
	if( csv < 0.0 ) csv = -0.05;
	if( csv > 1.0 ) csv = 1.0;

	int flavor = use_jet_hadronFlavour[iJet];

	ind_jet.push_back(iJet);

	HT30 += pt;

	jetPts.push_back(pt);
	jetEtas.push_back(eta);
	jetCSVs.push_back(csv);
	jetFlavors.push_back(flavor);

	numJet++;
	if( csv > 0.890 ){
	  numBtag++;
	}

	if( abs(flavor)==5 )      numJet_b++;
	else if( abs(flavor)==4 ) numJet_c++;
	else                      numJet_l++;
      }


      double HT30nocc = 0;
      for( int iJet=0; iJet<int(eve->jet_nocc_pt_.size()); iJet++ ){
	double pt  = eve->jet_nocc_pt_[iJet];
	double eta = eve->jet_nocc_eta_[iJet];

	if( !(pt>30. && fabs(eta)<3.0) ) continue;

	HT30nocc += pt;
      }

      double fill_HT30 = std::min( HTmax - 0.0001, HT30 );
      double fill_HT30nocc = std::min( HTmax - 0.0001, HT30nocc );

      int iSys = 0;

      // Calculate CSV weight
      double wgt_csv_hf, wgt_csv_lf, wgt_csv_cf;
      double wgt_csv = ( insample<0 ) ? 1 : get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);


      wgt = wgt_gen * wgt_lumi * use_wgt_pu * wgt_lepIdSF * wgt_LHEscale * wgt_csv;


      int fill_numJet  = std::min( MaxNjet_2e-1, numJet );
      int fill_numBtag = std::min( MaxNbtag_2e-1, numBtag );

      h_2e_numJet->Fill(fill_numJet,wgt);

      h_2e_numBtag->Fill(fill_numBtag,wgt);



      bool passHLTEle27HT200 = false;
      if( insample < 0 ) passHLTEle27HT200 = ( eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_==1 );
      else               passHLTEle27HT200 = ( eve->pass_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v_==1 );

      h_2e_L1HTT->Fill(eve->L1HTT_,wgt);

      h_2e_HT30->Fill(fill_HT30,wgt);
      if( (eve->L1HTT_ > 125.) ){
	h_2e_HT30_L1HTT125->Fill(fill_HT30,wgt);
	if( passHLTEle27HT200 ) h_2e_HT30_L1HTT125_passHLTEle27HT200->Fill(fill_HT30,wgt);
      }
      if( (eve->L1HTT_ > 100.) ){
	h_2e_HT30_L1HTT100->Fill(fill_HT30,wgt);
	if( passHLTEle27HT200 ) h_2e_HT30_L1HTT100_passHLTEle27HT200->Fill(fill_HT30,wgt);
      }

      h_2e_HT30nocc->Fill(fill_HT30nocc,wgt);
      if( (eve->L1HTT_ > 125.) ){
	h_2e_HT30nocc_L1HTT125->Fill(fill_HT30nocc,wgt);
	if( passHLTEle27HT200 ) h_2e_HT30nocc_L1HTT125_passHLTEle27HT200->Fill(fill_HT30nocc,wgt);
      }
      if( (eve->L1HTT_ > 100.) ){
	h_2e_HT30nocc_L1HTT100->Fill(fill_HT30nocc,wgt);
	if( passHLTEle27HT200 ) h_2e_HT30nocc_L1HTT100_passHLTEle27HT200->Fill(fill_HT30nocc,wgt);
      }

      h_2e_HT30nocc_HT30->Fill(fill_HT30nocc,fill_HT30,wgt);

      /// Lepton specific

      TLorentzVector myTag;
      myTag.SetPtEtaPhiE( eve->lepton_pt_[iTag], eve->lepton_eta_[iTag], eve->lepton_phi_[iTag], eve->lepton_energy_[iTag] );

      TLorentzVector myProbe;
      myProbe.SetPtEtaPhiE( eve->lepton_pt_[iProbe], eve->lepton_eta_[iProbe], eve->lepton_phi_[iProbe], eve->lepton_energy_[iProbe] );


      double pt = myProbe.Pt();
      double eta = eve->lepton_scEta_[iProbe];//myLep2.Eta();
      double phi = myProbe.Phi();


      bool passHLT=false;
      bool passL1T=false;
      for( int iHLT=0; iHLT<int(eve->hltEle27WP85Gsf_pt_.size()); iHLT++ ){
	std::string hltname = eve->hltEle27WP85Gsf_filter_[iHLT];

	double hlteta = eve->hltEle27WP85Gsf_eta_[iHLT];
	double hltphi = eve->hltEle27WP85Gsf_phi_[iHLT];

	double deltaR = DeltaR( eta, phi, hlteta, hltphi );
	if( deltaR<0.5 ){
	  if( hltname=="hltL1sL1SingleEG25" ) passL1T = true;
	  if( hltname=="hltEle27WPLooseGsfTrackIsoFilter" ||
	      hltname=="hltL1EG25Ele27WP85GsfTrackIsoFilter" ) passHLT = true;
	}
      }


      h_probe_pt_eta_all->Fill(pt, eta, wgt);
      if( passHLT ) h_probe_pt_eta_pHLTall->Fill(pt, eta, wgt);
      if( passHLT && !passL1T ) h_probe_pt_eta_pHLT_fL1T->Fill(pt, eta, wgt);
      if( passL1T ){
	h_probe_pt_eta_pL1T->Fill(pt, eta, wgt);
	if( passHLT ) h_probe_pt_eta_pHLT->Fill(pt, eta, wgt);
	else          h_probe_pt_eta_fHLT->Fill(pt, eta, wgt);
      }
      else{
	h_probe_pt_eta_fL1T->Fill(pt, eta, wgt);
      }


      if( fabs(eta) < 2.1 ){
	h_probe_pt_all->Fill(pt, wgt);
	if( passHLT ) h_probe_pt_pHLTall->Fill(pt, wgt);
	if( passHLT && !passL1T ) h_probe_pt_pHLT_fL1T->Fill(pt, wgt);
	if( passL1T ){
	  h_probe_pt_pL1T->Fill(pt, wgt);
	  if( passHLT ) h_probe_pt_pHLT->Fill(pt, wgt);
	  else          h_probe_pt_fHLT->Fill(pt, wgt);
	}
	else{
	  h_probe_pt_fL1T->Fill(pt, wgt);
	}
      }

      if( pt > 30 ){
	h_probe_eta_all->Fill(eta, wgt);
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
	if( passL1T ){
	  h_probe_numPVs_pL1T->Fill(numPVs, wgt);
	  if( passHLT ) h_probe_numPVs_pHLT->Fill(numPVs, wgt);
	  else          h_probe_numPVs_fHLT->Fill(numPVs, wgt);
	}
	else{
	  h_probe_numPVs_fL1T->Fill(numPVs, wgt);
	}
      }


      if( pt > 30 && fabs(eta) < 2.1 ){
	h_probe_HT30_all->Fill(HT30, wgt);
	if( passL1T ){
	  h_probe_HT30_pL1T->Fill(HT30, wgt);
	  if( passHLT ) h_probe_HT30_pHLT->Fill(HT30, wgt);
	  else          h_probe_HT30_fHLT->Fill(HT30, wgt);
	}
	else{
	  h_probe_HT30_fL1T->Fill(HT30, wgt);
	}
      }
    } // end require exactly two leptons TwoLep



    // Require exactly one lepton
    if( oneLep ){

      int lepInd = ( numEle>0 ) ? ind_ele[0] : ind_mu[0];

      TLorentzVector myLep;
      myLep.SetPtEtaPhiE( eve->lepton_pt_[lepInd], eve->lepton_eta_[lepInd], eve->lepton_phi_[lepInd], eve->lepton_energy_[lepInd] );


      vdouble use_jet_pt = eve->jet_pt_;
      vdouble use_jet_eta = eve->jet_eta_;
      vdouble use_jet_phi = eve->jet_phi_;
      vdouble use_jet_energy = eve->jet_energy_;

      vdouble use_jet_csv = eve->jet_csv_;
      vdouble use_jet_pileupJetId_fullDiscriminant = eve->jet_pileupJetId_fullDiscriminant_;

      vint use_jet_hadronFlavour = eve->jet_hadronFlavour_;



      // Set lepton scale factor equal to 1 for now
      double wgt_lepIdSF = 1.;

      // PU systematic
      double use_wgt_pu = wgt_pu;



      // Factorization and normalization scale uncertainty
      double wgt_LHEscale = 1.;
      if( insample>-1 ){
	double originalXWGTUP = eve->originalXWGTUP_;
	vdouble lhe_event_weights = eve->LHEEvent_weights_;

	if( lhe_event_weights.size()>6 && originalXWGTUP!=0 ){
	  int use_weight_index = 0;
	  wgt_LHEscale = lhe_event_weights[use_weight_index] / originalXWGTUP;
	}
      }



      vdouble jetPts;
      vdouble jetEtas;
      vdouble jetCSVs;
      vint    jetFlavors;

      vint ind_jet;

      int numJet = 0;
      int numBtag = 0;

      double HT30 = 0;
      TLorentzVector sumJet;
      bool firstJet = false;

      for( int iJet=0; iJet<int(use_jet_pt.size()); iJet++ ){
	TLorentzVector myJet;
	myJet.SetPtEtaPhiE( use_jet_pt[iJet], use_jet_eta[iJet], use_jet_phi[iJet], use_jet_energy[iJet] );

	if( !firstJet ) sumJet = myJet;
	else            sumJet += myJet;

	double pt  = use_jet_pt[iJet];
	double eta = use_jet_eta[iJet];

	if( !(pt>30. && fabs(eta)<2.4) ) continue;

	double dR = myJet.DeltaR(myLep);

	if( dR < 0.4 ) continue;

	double csv = use_jet_csv[iJet];
	if( csv < 0.0 ) csv = -0.05;
	if( csv > 1.0 ) csv = 1.0;

	int flavor = use_jet_hadronFlavour[iJet];

	ind_jet.push_back(iJet);

	HT30 += pt;


	jetPts.push_back(pt);
	jetEtas.push_back(eta);
	jetCSVs.push_back(csv);
	jetFlavors.push_back(flavor);

	numJet++;
	if( csv > 0.890 ){
	  numBtag++;
	}
      }


      double HT30nocc = 0;
      for( int iJet=0; iJet<int(eve->jet_nocc_pt_.size()); iJet++ ){
	double pt  = eve->jet_nocc_pt_[iJet];
	double eta = eve->jet_nocc_eta_[iJet];

	if( !(pt>30. && fabs(eta)<3.0) ) continue;

	HT30nocc += pt;
      }


      double fill_HT30 = std::min( HTmax - 0.0001, HT30 );
      double fill_HT30nocc = std::min( HTmax - 0.0001, HT30nocc );

      int iSys = 0;

      // Calculate CSV weight
      double wgt_csv_hf, wgt_csv_lf, wgt_csv_cf;
      double wgt_csv = ( insample<0 ) ? 1 : get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);

      wgt = wgt_gen * wgt_lumi * use_wgt_pu * wgt_lepIdSF * wgt_LHEscale * wgt_csv;




      // Require at least four jets
      if( !(numJet>=4) ) continue;



      // Require at least two b-tags
      if( !(numBtag>=2) ) continue;

      numEvents_wgt_gen_lumi_pu_4j2t += wgt_gen * wgt_lumi * wgt_pu;


      int fill_numJet  = std::min( MaxNjet_1l4j2t-1, numJet );
      int fill_numBtag = std::min( MaxNbtag_1l4j2t-1, numBtag );

      h_1l4j2t_numJet->Fill(fill_numJet,wgt);
      h_1l4j2t_numBtag->Fill(fill_numBtag,wgt);



      ////////////


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

      if( oneEle ){
	h_category_yield_1e->Fill(0.5,wgt);
	h_category_yield_1e->Fill(this_category,wgt);

	if( eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_==1 || eve->pass_HLT_Ele27_WP85_Gsf_v_==1 ){
	  h_category_yield_1e_HLT_Ele27->Fill(0.5,wgt);
	  h_category_yield_1e_HLT_Ele27->Fill(this_category,wgt);
	}
	if( eve->pass_HLT_Ele23_WPLoose_Gsf_v_==1 || eve->pass_HLT_Ele22_eta2p1_WP75_Gsf_v_==1 ){
	  h_category_yield_1e_HLT_Ele23->Fill(0.5,wgt);
	  h_category_yield_1e_HLT_Ele23->Fill(this_category,wgt);
	}
      }
      else if( oneMu ){
	h_category_yield_1m->Fill(0.5,wgt);
	h_category_yield_1m->Fill(this_category,wgt);
      }


      bool passHLTEle27HT200 = false;
      if( insample < 0 ) passHLTEle27HT200 = ( eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_==1 );
      else               passHLTEle27HT200 = ( eve->pass_HLT_Ele27_eta2p1_WP85_Gsf_HT200_v_==1 );

      oneEle = ( oneEle && pass_trigger_ele );

      if( oneEle && HT30>300 ){
	h_category_yield_1e_recoHT300->Fill(0.5,wgt);
	h_category_yield_1e_recoHT300->Fill(this_category,wgt);

	if( eve->L1HTT_ > 125. ){
	  h_category_yield_1e_recoHT300_L1HTT125->Fill(0.5,wgt);
	  h_category_yield_1e_recoHT300_L1HTT125->Fill(this_category,wgt);

	  if( passHLTEle27HT200 ){
	    h_category_yield_1e_recoHT300_L1HTT125_passHLTEle27HT200->Fill(0.5,wgt);
	    h_category_yield_1e_recoHT300_L1HTT125_passHLTEle27HT200->Fill(this_category,wgt);
	  }
	}
	if( eve->L1HTT_ > 100. ){
	  h_category_yield_1e_recoHT300_L1HTT100->Fill(0.5,wgt);
	  h_category_yield_1e_recoHT300_L1HTT100->Fill(this_category,wgt);

	  if( passHLTEle27HT200 ){
	    h_category_yield_1e_recoHT300_L1HTT100_passHLTEle27HT200->Fill(0.5,wgt);
	    h_category_yield_1e_recoHT300_L1HTT100_passHLTEle27HT200->Fill(this_category,wgt);
	  }
	}
      }


      if( oneEle ){
	h_1e4j2t_L1HTT->Fill(eve->L1HTT_,wgt);

	h_1e4j2t_HT30->Fill(fill_HT30,wgt);
	if( (eve->L1HTT_ > 125.) ){
	  h_1e4j2t_HT30_L1HTT125->Fill(fill_HT30,wgt);
	  if( passHLTEle27HT200 ) h_1e4j2t_HT30_L1HTT125_passHLTEle27HT200->Fill(fill_HT30,wgt);
	}
	if( (eve->L1HTT_ > 100.) ){
	  h_1e4j2t_HT30_L1HTT100->Fill(fill_HT30,wgt);
	  if( passHLTEle27HT200 ) h_1e4j2t_HT30_L1HTT100_passHLTEle27HT200->Fill(fill_HT30,wgt);
	}


	h_1e4j2t_HT30nocc->Fill(fill_HT30nocc,wgt);
	if( (eve->L1HTT_ > 125.) ){
	  h_1e4j2t_HT30nocc_L1HTT125->Fill(fill_HT30nocc,wgt);
	  if( passHLTEle27HT200 ) h_1e4j2t_HT30nocc_L1HTT125_passHLTEle27HT200->Fill(fill_HT30nocc,wgt);
	}
	if( (eve->L1HTT_ > 100.) ){
	  h_1e4j2t_HT30nocc_L1HTT100->Fill(fill_HT30nocc,wgt);
	  if( passHLTEle27HT200 ) h_1e4j2t_HT30nocc_L1HTT100_passHLTEle27HT200->Fill(fill_HT30nocc,wgt);
	}

	h_1e4j2t_HT30nocc_HT30->Fill(fill_HT30nocc,fill_HT30,wgt);

      }

    } // end if oneLep
  } // end loop over events

  std::cout << "**************************************************************" << std::endl;
  std::cout << "\t Number of all events = " << numEvents_all << std::endl;
  std::cout << "\t Number of gen weighted events = " << numEvents_wgt_gen << std::endl;
  std::cout << "\t Number of lumi weighted events = " << numEvents_wgt_lumi << std::endl;
  std::cout << "\t Number of gen * lumi weighted events = " << numEvents_wgt_gen_lumi << std::endl;
  std::cout << "\t Number of gen * lumi * PU weighted events = " << numEvents_wgt_gen_lumi_pu << std::endl;
  std::cout << "\t Number of gen * lumi * PU * CSV weighted events = " << numEvents_wgt_gen_lumi_pu_csv << std::endl;
  std::cout << "\t Number of gen * lumi * PU weighted events (>=4j,>=2b) = " << numEvents_wgt_gen_lumi_pu_4j2t << std::endl;
  std::cout << "\t Number of gen * lumi * PU * CSV weighted events (>=4j,>=2b) = " << numEvents_wgt_gen_lumi_pu_csv_4j2t  << std::endl;
  std::cout << "**************************************************************" << std::endl;



  std::cout << " Done! " << std::endl;

  histofile.Write();
  histofile.Close();

}


double reweightPU( int nPU, int iSys ){

  double PUscale[50];

  if( iSys==0 ){
  PUscale[0] = 110.875;
  PUscale[1] = 140.546;
  PUscale[2] = 98.7647;
  PUscale[3] = 32.3564;
  PUscale[4] = 18.0196;
  PUscale[5] = 3.41092;
  PUscale[6] = 2.02401;
  PUscale[7] = 2.63715;
  PUscale[8] = 3.51755;
  PUscale[9] = 3.42482;
  PUscale[10] = 3.09016;
  PUscale[11] = 2.72066;
  PUscale[12] = 2.13645;
  PUscale[13] = 1.42201;
  PUscale[14] = 0.791452;
  PUscale[15] = 0.369213;
  PUscale[16] = 0.152774;
  PUscale[17] = 0.0649893;
  PUscale[18] = 0.033305;
  PUscale[19] = 0.0189498;
  PUscale[20] = 0.00993347;
  PUscale[21] = 0.0043885;
  PUscale[22] = 0.0016274;
  PUscale[23] = 0.000543008;
  PUscale[24] = 0.000188179;
  PUscale[25] = 8.05937e-05;
  PUscale[26] = 4.40758e-05;
  PUscale[27] = 2.86902e-05;
  PUscale[28] = 2.07412e-05;
  PUscale[29] = 1.56476e-05;
  PUscale[30] = 1.11086e-05;
  PUscale[31] = 6.70872e-06;
  PUscale[32] = 3.27221e-06;
  PUscale[33] = 1.36206e-06;
  PUscale[34] = 4.99768e-07;
  PUscale[35] = 1.70217e-07;
  PUscale[36] = 5.35254e-08;
  PUscale[37] = 1.65078e-08;
  PUscale[38] = 4.90562e-09;
  PUscale[39] = 1.3586e-09;
  PUscale[40] = 3.61389e-10;
  PUscale[41] = 9.30517e-11;
  PUscale[42] = 2.34651e-11;
  PUscale[43] = 5.51348e-12;
  PUscale[44] = 1.30187e-12;
  PUscale[45] = 2.87977e-13;
  PUscale[46] = 1.49461e-13;
  PUscale[47] = 4.44305e-14;
  PUscale[48] = 3.883e-14;
  PUscale[49] = 0;
  }
  else if( iSys==1 ){
  PUscale[0] = 93.4113;
  PUscale[1] = 128.951;
  PUscale[2] = 90.6777;
  PUscale[3] = 28.8355;
  PUscale[4] = 14.9844;
  PUscale[5] = 2.51391;
  PUscale[6] = 1.28157;
  PUscale[7] = 1.51784;
  PUscale[8] = 2.30492;
  PUscale[9] = 2.64525;
  PUscale[10] = 2.71343;
  PUscale[11] = 2.63104;
  PUscale[12] = 2.29897;
  PUscale[13] = 1.75475;
  PUscale[14] = 1.14609;
  PUscale[15] = 0.631047;
  PUscale[16] = 0.297615;
  PUscale[17] = 0.127901;
  PUscale[18] = 0.0583678;
  PUscale[19] = 0.0324697;
  PUscale[20] = 0.0201164;
  PUscale[21] = 0.0116749;
  PUscale[22] = 0.00583743;
  PUscale[23] = 0.00250793;
  PUscale[24] = 0.000968993;
  PUscale[25] = 0.000381575;
  PUscale[26] = 0.000178639;
  PUscale[27] = 0.000105926;
  PUscale[28] = 7.49081e-05;
  PUscale[29] = 5.80434e-05;
  PUscale[30] = 4.3952e-05;
  PUscale[31] = 2.91454e-05;
  PUscale[32] = 1.59293e-05;
  PUscale[33] = 7.53183e-06;
  PUscale[34] = 3.1692e-06;
  PUscale[35] = 1.24681e-06;
  PUscale[36] = 4.55633e-07;
  PUscale[37] = 1.64214e-07;
  PUscale[38] = 5.73296e-08;
  PUscale[39] = 1.87497e-08;
  PUscale[40] = 5.9201e-09;
  PUscale[41] = 1.81867e-09;
  PUscale[42] = 5.49979e-10;
  PUscale[43] = 1.55764e-10;
  PUscale[44] = 4.45584e-11;
  PUscale[45] = 1.20037e-11;
  PUscale[46] = 7.62258e-12;
  PUscale[47] = 2.83344e-12;
  PUscale[48] = 2.39624e-12;
  PUscale[49] = 6.33247e-13;
  }
  else if( iSys==-1 ){
  PUscale[0] = 129.869;
  PUscale[1] = 153.371;
  PUscale[2] = 109.048;
  PUscale[3] = 36.945;
  PUscale[4] = 22.3864;
  PUscale[5] = 4.90528;
  PUscale[6] = 3.4667;
  PUscale[7] = 4.54737;
  PUscale[8] = 5.06312;
  PUscale[9] = 4.18434;
  PUscale[10] = 3.34869;
  PUscale[11] = 2.64698;
  PUscale[12] = 1.8101;
  PUscale[13] = 1.01708;
  PUscale[14] = 0.471125;
  PUscale[15] = 0.187396;
  PUscale[16] = 0.0743086;
  PUscale[17] = 0.0355349;
  PUscale[18] = 0.0189565;
  PUscale[19] = 0.00902975;
  PUscale[20] = 0.00350245;
  PUscale[21] = 0.00111462;
  PUscale[22] = 0.000317195;
  PUscale[23] = 9.74159e-05;
  PUscale[24] = 3.82395e-05;
  PUscale[25] = 1.92163e-05;
  PUscale[26] = 1.13268e-05;
  PUscale[27] = 7.38894e-06;
  PUscale[28] = 5.08901e-06;
  PUscale[29] = 3.51665e-06;
  PUscale[30] = 2.22272e-06;
  PUscale[31] = 1.17268e-06;
  PUscale[32] = 4.93402e-07;
  PUscale[33] = 1.75551e-07;
  PUscale[34] = 5.4654e-08;
  PUscale[35] = 1.56911e-08;
  PUscale[36] = 4.13353e-09;
  PUscale[37] = 1.06154e-09;
  PUscale[38] = 2.61117e-10;
  PUscale[39] = 5.95034e-11;
  PUscale[40] = 1.29466e-11;
  PUscale[41] = 2.71056e-12;
  PUscale[42] = 5.52478e-13;
  PUscale[43] = 1.04327e-13;
  PUscale[44] = 1.97209e-14;
  PUscale[45] = 3.46906e-15;
  PUscale[46] = 1.35347e-15;
  PUscale[47] = 0;
  PUscale[48] = 0;
  PUscale[49] = 0;
  }

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


/// algos / supporting functions:

void fillCSVhistos(TFile* fileHF, TFile* fileLF){

  for( int iSys=0; iSys<9; iSys++ ){
    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) h_csv_wgt_hf[iSys][iPt] = NULL;
    for( int iPt=0; iPt<3; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
    }
  }
  for( int iSys=0; iSys<5; iSys++ ){
    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) c_csv_wgt_hf[iSys][iPt] = NULL;
  }

  // CSV reweighting /// only care about the nominal ones
  for( int iSys=0; iSys<9; iSys++ ){
    TString syst_csv_suffix_hf = "final";
    TString syst_csv_suffix_c = "final";
    TString syst_csv_suffix_lf = "final";
    
    switch( iSys ){
    case 0:
      // this is the nominal case
      break;
    case 1:
      // JESUp
      syst_csv_suffix_hf = "final_JESUp"; syst_csv_suffix_lf = "final_JESUp";
      syst_csv_suffix_c  = "final_cErr1Up";
      break;
    case 2:
      // JESDown
      syst_csv_suffix_hf = "final_JESDown"; syst_csv_suffix_lf = "final_JESDown";
      syst_csv_suffix_c  = "final_cErr1Down";
      break;
    case 3:
      // purity up
      syst_csv_suffix_hf = "final_LFUp"; syst_csv_suffix_lf = "final_HFUp";
      syst_csv_suffix_c  = "final_cErr2Up";
      break;
    case 4:
      // purity down
      syst_csv_suffix_hf = "final_LFDown"; syst_csv_suffix_lf = "final_HFDown";
      syst_csv_suffix_c  = "final_cErr2Down";
      break;
    case 5:
      // stats1 up
      syst_csv_suffix_hf = "final_Stats1Up"; syst_csv_suffix_lf = "final_Stats1Up";
      break;
    case 6:
      // stats1 down
      syst_csv_suffix_hf = "final_Stats1Down"; syst_csv_suffix_lf = "final_Stats1Down";
      break;
    case 7:
      // stats2 up
      syst_csv_suffix_hf = "final_Stats2Up"; syst_csv_suffix_lf = "final_Stats2Up";
      break;
    case 8:
      // stats2 down
      syst_csv_suffix_hf = "final_Stats2Down"; syst_csv_suffix_lf = "final_Stats2Down";
      break;
    }

    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) h_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );

    if( iSys<5 ){
      for( int iPt=0; iPt<PtBinsHF_; iPt++ ) c_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
    }
    
    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
    }
  }

  return;
}



double get_csv_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCSV, vint jetFlavor, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF ){

  // force using no systematic here
  //iSys = 0;

  int iSysHF = 0;
  switch(iSys){
  case 7:  iSysHF=1; break; //JESUp
  case 8:  iSysHF=2; break; //JESDown
  case 9:  iSysHF=3; break; //LFUp
  case 10: iSysHF=4; break; //LFDown
  case 13: iSysHF=5; break; //Stats1Up
  case 14: iSysHF=6; break; //Stats1Down
  case 15: iSysHF=7; break; //Stats2Up
  case 16: iSysHF=8; break; //Stats2Down
  default : iSysHF = 0; break; //NoSys
  }

  int iSysC = 0;
  switch(iSys){
  case 21: iSysC=1; break;
  case 22: iSysC=2; break;
  case 23: iSysC=3; break;
  case 24: iSysC=4; break;
  default : iSysC = 0; break;
  }

  int iSysLF = 0;
  switch(iSys){
  case 7:  iSysLF=1; break; //JESUp
  case 8:  iSysLF=2; break; //JESDown
  case 11: iSysLF=3; break; //HFUp
  case 12: iSysLF=4; break; //HFDown
  case 17: iSysLF=5; break; //Stats1Up
  case 18: iSysLF=6; break; //Stats1Down
  case 19: iSysLF=7; break; //Stats2Up
  case 20: iSysLF=8; break; //Stats2Down
  default : iSysLF = 0; break; //NoSys
  }

  double csvWgthf = 1.;
  double csvWgtC  = 1.;
  double csvWgtlf = 1.;

  for( int iJet=0; iJet<int(jetPt.size()); iJet++ ){

    double csv = jetCSV[iJet];
    double pt = jetPt[iJet];
    double jetAbsEta = fabs(jetEta[iJet]);
    int flavor = jetFlavor[iJet];

    int iPt = -1; int iEta = -1;
    if (pt >=19.99 && pt<30) iPt = 0;
    else if (pt >=30 && pt<40) iPt = 1;
    else if (pt >=40 && pt<60) iPt = 2;
    else if (pt >=60 && pt<100) iPt = 3;
    else if (pt >=100)          iPt = 4;
    //else if (pt >=100 && pt<160) iPt = 4;
    //else if (pt >=160 && pt<10000) iPt = 5;

    if( PtBinsHF_ > 5 && pt >=160 )  iPt = 5;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8 ) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6 )  iEta = 1;
    else if ( jetAbsEta>=1.6 && jetAbsEta<2.41 ) iEta = 2;

    if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, pt = " << pt << ", jetAbsEta = " << jetAbsEta << std::endl;

    //std::cout << " flavor = " << flavor << ", csv = " << csv << ", iPt = " << iPt << ", iEta = " << iEta << ", iSysHF = " << iSysHF << ", iSysC = " << iSysC << ", iSysLF = " << iSysLF << std::endl;

    if (abs(flavor) == 5 ){
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;

      // if( iSysHF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      //  			     flavor, pt, jetAbsEta, csv, iCSVWgtHF );
    }
    else if( abs(flavor) == 4 ){
      int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
      double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;
      // if( iSysC==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      // 			    flavor, pt, jetAbsEta, csv, iCSVWgtC );
    }
    else {
      if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
      if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;
      // if( iSysLF==0 ) printf(" iSysLF = %d, iJet = %d,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      //  			     iSysLF, iJet, flavor, pt, jetAbsEta, csv, iCSVWgtLF );
    }
  }

  double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

  csvWgtHF = csvWgthf;
  csvWgtLF = csvWgtlf;
  csvWgtCF = csvWgtC;

  //std::cout << "\t iSys = " << iSys << " total csvWgtLF = " << csvWgtLF << std::endl;

  return csvWgtTotal;
}



void fillCSVEffhistos(TFile* file){

  h_a_jet_pt_eta_all_ = (TH2D*)file->Get("h_a_jet_pt_eta_all")->Clone("h_a_jet_pt_eta_all_temp");
  h_b_jet_pt_eta_all_ = (TH2D*)file->Get("h_b_jet_pt_eta_all")->Clone("h_b_jet_pt_eta_all_temp");
  h_c_jet_pt_eta_all_ = (TH2D*)file->Get("h_c_jet_pt_eta_all")->Clone("h_c_jet_pt_eta_all_temp");
  h_l_jet_pt_eta_all_ = (TH2D*)file->Get("h_l_jet_pt_eta_all")->Clone("h_l_jet_pt_eta_all_temp");

  h_a_jet_pt_eta_eff_ = (TH2D*)file->Get("h_a_jet_pt_eta_csvM")->Clone("h_a_jet_pt_eta_eff_temp");
  h_b_jet_pt_eta_eff_ = (TH2D*)file->Get("h_b_jet_pt_eta_csvM")->Clone("h_b_jet_pt_eta_eff_temp");
  h_c_jet_pt_eta_eff_ = (TH2D*)file->Get("h_c_jet_pt_eta_csvM")->Clone("h_c_jet_pt_eta_eff_temp");
  h_l_jet_pt_eta_eff_ = (TH2D*)file->Get("h_l_jet_pt_eta_csvM")->Clone("h_l_jet_pt_eta_eff_temp");

  h_a_jet_pt_eta_eff_->Divide(h_a_jet_pt_eta_all_);
  h_b_jet_pt_eta_eff_->Divide(h_b_jet_pt_eta_all_);
  h_c_jet_pt_eta_eff_->Divide(h_c_jet_pt_eta_all_);
  h_l_jet_pt_eta_eff_->Divide(h_l_jet_pt_eta_all_);
}


double get_btv_csv_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCSV, vint jetFlavor, TString sys_name, bool use_btv, bool verbose ){

  double mcTag = 1.;
  double mcNoTag = 1.;
  double dataTag = 1.;
  double dataNoTag = 1.;

  bool debug = verbose;

  int numJet = 0;
  int numBtag = 0;
  for( int iJet=0; iJet<int(jetPt.size()); iJet++ ){

    double csv = jetCSV[iJet];
    double pt = jetPt[iJet];
    double eta = jetEta[iJet];
    double jetAbsEta = fabs(jetEta[iJet]);
    int flavor = jetFlavor[iJet];

    if( !( (pt > 30) && (jetAbsEta < 2.4) ) ) continue;

    BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;

    if( abs(flavor)==5 )      jf = BTagEntry::FLAV_B;
    else if( abs(flavor)==4 ) jf = BTagEntry::FLAV_C;
    else                      jf = BTagEntry::FLAV_UDSG;

    double eff = get_csv_efficiency(pt,eta,flavor);

    if( !(eff > -0.000001 && eff < 1.01) ) debug = true;

    bool isTag = ( csv > 0.89 );

    numJet++;
    if( isTag ) numBtag++;

    double SF = 1.;
    if( use_btv ){
      if( abs(flavor)==5 || abs(flavor)==4 ){
	if( sys_name.Contains("LFUp") )         SF = reader_btv_mujets_up.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("LFDown") )  SF = reader_btv_mujets_down.eval(jf, eta, pt, csv);
	else                                    SF = reader_btv_mujets.eval(jf, eta, pt);
      }
      else {
	if( sys_name.Contains("HFUp") )         SF = reader_btv_comb_up.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("HFDown") )  SF = reader_btv_comb_down.eval(jf, eta, pt, csv);
	else                                    SF = reader_btv_comb.eval(jf, eta, pt);
      }
    }

    if( isTag && eff==0. ) debug = true;
    if( !isTag && eff==1. ) debug = true;

    if( isTag ){
      mcTag *= eff; 
      dataTag *= eff*SF;
    }
    else{
      mcNoTag *= (1 - eff); 
      dataNoTag *= (1 - eff*SF);
    }

    if( !sys_name.Contains("_") && debug ){
      printf("  iJet = %d, jet pt = %.1f, eta = %+.2f, csv = %.2f, flavor = %+d, eff = %.3f, SF = %.3f \n",
	     iJet, pt, eta, csv, flavor, eff, SF);
    }
  }

  double evt_wgt = (dataNoTag * dataTag ) / ( mcNoTag * mcTag );

  if( evt_wgt != evt_wgt ) debug = true;

  if( evt_wgt>10. || evt_wgt < 0. ) debug = true;

  //if( !sys_name.Contains("_JES") && debug ){
  if( !sys_name.Contains("_") &&  debug ){
    printf(" ==> Event: weight = %.3f, dataNoTag = %.3f, dataTag = %.3f, mcNoTag = %.3f, mcTag = %.3f, numJet = %d, numBtag = %d, sys = %s \n",
	   evt_wgt, dataNoTag, dataTag, mcNoTag, mcTag, numJet, numBtag, sys_name.Data());
  }

  if( evt_wgt != evt_wgt ) evt_wgt = 1.0;

  if( mcTag==0 && dataTag==0 ){
    mcTag = 1.;
    dataTag = 1.;
    evt_wgt = (dataNoTag * dataTag ) / ( mcNoTag * mcTag );
  }

  if( mcNoTag==0 ){
    mcNoTag = 1.;
    evt_wgt = (dataNoTag * dataTag ) / ( mcNoTag * mcTag );
  }

  if( !sys_name.Contains("_") && debug ){
    printf(" ==> NEW Event: weight = %.3f, dataNoTag = %.3f, dataTag = %.3f, mcNoTag = %.3f, mcTag = %.3f, numJet = %d, numBtag = %d, sys = %s \n",
	   evt_wgt, dataNoTag, dataTag, mcNoTag, mcTag, numJet, numBtag, sys_name.Data());
  }

  return evt_wgt;
}


double get_csv_efficiency( double jetPt, double jetEta, int jetFlavor ){

  if( jetPt>500. ) jetPt = 500 - 0.01;

  jetEta = fabs(jetEta);

  int useBin = h_a_jet_pt_eta_eff_->FindBin(jetPt,jetEta);

  double content = h_a_jet_pt_eta_all_->GetBinContent(useBin);
  if( abs(jetFlavor)==5 )      content = h_b_jet_pt_eta_all_->GetBinContent(useBin);
  else if( abs(jetFlavor)==4 ) content = h_c_jet_pt_eta_all_->GetBinContent(useBin);
  else                         content = h_l_jet_pt_eta_all_->GetBinContent(useBin);

  if( content < 1 ){
    while( content < 1 ){
      // printf("  get_csv_efficiency: useBin %d has content < 1: jetPt = %.1f, jetEta = %+.2f, jetFlavor = %d: content = %4.2f \n",
      // 	     useBin, jetPt, jetEta, jetFlavor, content);
      useBin -= 1;
      if( useBin<1 ) break;
      if( abs(jetFlavor)==5 )      content = h_b_jet_pt_eta_all_->GetBinContent(useBin);
      else if( abs(jetFlavor)==4 ) content = h_c_jet_pt_eta_all_->GetBinContent(useBin);
      else                         content = h_l_jet_pt_eta_all_->GetBinContent(useBin);
    }
  }

  double efficiency = 1.;
  if( abs(jetFlavor)==5 )      efficiency = h_b_jet_pt_eta_eff_->GetBinContent(useBin);
  else if( abs(jetFlavor)==4 ) efficiency = h_c_jet_pt_eta_eff_->GetBinContent(useBin);
  else                         efficiency = h_l_jet_pt_eta_eff_->GetBinContent(useBin);

  return efficiency;
}


/*





 */
