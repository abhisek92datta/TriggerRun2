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
void fillCSVhistos_old(TFile *fileHF, TFile *fileLF);
void fillCMVAhistos(TFile *fileHF, TFile *fileLF);
double get_csv_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCSV, vint jetFlavor, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );
double get_csv_wgt_old( vdouble jetPt, vdouble jetEta, vdouble jetCSV, vint jetFlavor, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );
double get_cmva_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCMVA, vint jetFlavor, int iSys, double &cmvaWgtHF, double &cmvaWgtLF, double &cmvaWgtCF );

double get_btv_csv_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCSV, vint jetFlavor, TString sys_name, bool use_csv, bool verbose );

void fillCSVEffhistos(TFile *file);
double get_csv_efficiency( double jetPt, double jetEta, int jetFlavor );


int PtBinsHF_ = 5;

// CSV reweighting
TH1D* h_csv_wgt_hf[9][5];
TH1D* c_csv_wgt_hf[9][5];
TH1D* h_csv_wgt_lf[9][4][3];

// old CSV reweighting, for comparison
TH1D* h_csv_wgt_hf_old[9][5];
TH1D* c_csv_wgt_hf_old[9][5];
TH1D* h_csv_wgt_lf_old[9][4][3];

// CSV reweighting
TH1D* h_cmva_wgt_hf[9][5];
TH1D* c_cmva_wgt_hf[9][5];
TH1D* h_cmva_wgt_lf[9][4][3];

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



// setup calibration readers
BTagCalibration calib_btv_dil_csvv2("csvv2", "CSVv2_TagCountTT.csv");

BTagCalibrationReader reader_btv_dil(&calib_btv_dil_csvv2, // calibration instance
				     BTagEntry::OP_MEDIUM, // operating point
				     "tctt", // measurement type
				     "central"); // systematics type

BTagCalibrationReader reader_btv_dil_up(&calib_btv_dil_csvv2, // calibration instance
					BTagEntry::OP_MEDIUM, // operating point
					"tctt", // measurement type
					"up"); // systematics type

BTagCalibrationReader reader_btv_dil_down(&calib_btv_dil_csvv2, // calibration instance
					  BTagEntry::OP_MEDIUM, // operating point
					  "tctt", // measurement type
					  "down"); // systematics type



// setup calibration readers
//BTagCalibration calib_csvv2("csvv2", "ttH_BTV_CSVv2_13TeV_2015D_20151122.csv");
BTagCalibration calib_csvv2("csvv2", "ttH_BTV_CSVv2_13TeV_2015D_2016_02_08.csv");
BTagCalibrationReader reader(&calib_csvv2, // calibration instance
			     BTagEntry::OP_RESHAPING, // operating point
			     "iterativefit", // measurement type
			     "central"); // systematics type

// JESUp
BTagCalibrationReader reader_JESUp(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_jes"); // systematics type
// JESDown
BTagCalibrationReader reader_JESDown(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_jes"); // systematics type

// LFUp
BTagCalibrationReader reader_LFUp(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_lf"); // systematics type
// LFDown
BTagCalibrationReader reader_LFDown(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_lf"); // systematics type

// HFUp
BTagCalibrationReader reader_HFUp(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_hf"); // systematics type
// HFDown
BTagCalibrationReader reader_HFDown(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_hf"); // systematics type

// HFStats1Up
BTagCalibrationReader reader_HFStats1Up(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_hfstats1"); // systematics type
// HFStats1Down
BTagCalibrationReader reader_HFStats1Down(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_hfstats1"); // systematics type

// HFStats2Up
BTagCalibrationReader reader_HFStats2Up(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_hfstats2"); // systematics type
// HFStats2Down
BTagCalibrationReader reader_HFStats2Down(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_hfstats2"); // systematics type

// LFStats1Up
BTagCalibrationReader reader_LFStats1Up(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_lfstats1"); // systematics type
// LFStats1Down
BTagCalibrationReader reader_LFStats1Down(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_lfstats1"); // systematics type

// LFStats2Up
BTagCalibrationReader reader_LFStats2Up(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_lfstats2"); // systematics type
// LFStats2Down
BTagCalibrationReader reader_LFStats2Down(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_lfstats2"); // systematics type

// CFErr1Up
BTagCalibrationReader reader_CFErr1Up(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_cferr1"); // systematics type
// CFErr1Down
BTagCalibrationReader reader_CFErr1Down(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_cferr1"); // systematics type

// CFErr2Up
BTagCalibrationReader reader_CFErr2Up(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_cferr2"); // systematics type
// CFErr2Down
BTagCalibrationReader reader_CFErr2Down(&calib_csvv2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_cferr2"); // systematics type





BTagCalibration calib_cmvav2("cmvav2", "ttH_BTV_CMVAv2_13TeV_2015D_2016_02_12.csv");
BTagCalibrationReader reader_cmvav2(&calib_cmvav2, // calibration instance
			     BTagEntry::OP_RESHAPING, // operating point
			     "iterativefit", // measurement type
			     "central"); // systematics type

// JESUp
BTagCalibrationReader reader_cmvav2_JESUp(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_jes"); // systematics type
// JESDown
BTagCalibrationReader reader_cmvav2_JESDown(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_jes"); // systematics type

// LFUp
BTagCalibrationReader reader_cmvav2_LFUp(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_lf"); // systematics type
// LFDown
BTagCalibrationReader reader_cmvav2_LFDown(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_lf"); // systematics type

// HFUp
BTagCalibrationReader reader_cmvav2_HFUp(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_hf"); // systematics type
// HFDown
BTagCalibrationReader reader_cmvav2_HFDown(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_hf"); // systematics type

// HFStats1Up
BTagCalibrationReader reader_cmvav2_HFStats1Up(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_hfstats1"); // systematics type
// HFStats1Down
BTagCalibrationReader reader_cmvav2_HFStats1Down(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_hfstats1"); // systematics type

// HFStats2Up
BTagCalibrationReader reader_cmvav2_HFStats2Up(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_hfstats2"); // systematics type
// HFStats2Down
BTagCalibrationReader reader_cmvav2_HFStats2Down(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_hfstats2"); // systematics type

// LFStats1Up
BTagCalibrationReader reader_cmvav2_LFStats1Up(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_lfstats1"); // systematics type
// LFStats1Down
BTagCalibrationReader reader_cmvav2_LFStats1Down(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_lfstats1"); // systematics type

// LFStats2Up
BTagCalibrationReader reader_cmvav2_LFStats2Up(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_lfstats2"); // systematics type
// LFStats2Down
BTagCalibrationReader reader_cmvav2_LFStats2Down(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_lfstats2"); // systematics type

// CFErr1Up
BTagCalibrationReader reader_cmvav2_CFErr1Up(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_cferr1"); // systematics type
// CFErr1Down
BTagCalibrationReader reader_cmvav2_CFErr1Down(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_cferr1"); // systematics type

// CFErr2Up
BTagCalibrationReader reader_cmvav2_CFErr2Up(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "up_cferr2"); // systematics type
// CFErr2Down
BTagCalibrationReader reader_cmvav2_CFErr2Down(&calib_cmvav2, // calibration instance
				   BTagEntry::OP_RESHAPING, // operating point
				   "iterativefit", // measurement type
				   "down_cferr2"); // systematics type



//*****************************************************************************

void ttHbb_data2mc_treeReader( int insample=1, int maxNentries=-1, int Njobs=1, int jobN=1, double intLumi=-1, int ttCat_=-1, bool useHTbins_=false, bool useCondor_=false ) {

  std::string inputFileHF = "data/csv_rwt_hf_IT_FlatSF_2015_07_27.root";
  std::string inputFileLF = "data/csv_rwt_lf_IT_FlatSF_2015_07_27.root";


   TFile* f_CSVwgt_HF = new TFile("/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/TriggerRun2/TriggerAnalyzer/data/csv_rwt_fit_hf_76x_2016_02_08.root");
   TFile* f_CSVwgt_LF = new TFile("/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/TriggerRun2/TriggerAnalyzer/data/csv_rwt_fit_lf_76x_2016_02_08.root");

  
  fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);

  TFile* f_CSVwgt_HF_old = new TFile("/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/TriggerRun2/TriggerAnalyzer/data/csv_rwt_fit_hf_2015_11_20.root");
  TFile* f_CSVwgt_LF_old = new TFile("/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/TriggerRun2/TriggerAnalyzer/data/csv_rwt_fit_lf_2015_11_20.root");

  fillCSVhistos_old(f_CSVwgt_HF_old, f_CSVwgt_LF_old);

  TFile* f_CMVAwgt_HF = new TFile("/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/csvReweightingRun2/csvTreeMaker/cmva_rwt_fit_hf_v2_final_2016_02_11_76x_rescaleHF_rescaleCF.root");
  TFile* f_CMVAwgt_LF = new TFile("/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/csvReweightingRun2/csvTreeMaker/cmva_rwt_fit_lf_v2_final_2016_02_11_76x_rescaleLF.root");


  
  fillCMVAhistos(f_CMVAwgt_HF, f_CMVAwgt_LF);


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
    mySample_nGen_ = 97994442;//116591749;//25357774;//25446993;
    mySample_sampleName_ = "ttbar_powheg";
    // ext
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1_triggerTree_v1/160223_114300/0000/");
    if( ttCat_>=0 ){
      if( ttCat_==0 ) mySample_sampleName_ = "ttlf_powheg";
      if( ttCat_==1 ) mySample_sampleName_ = "ttcc_powheg";
      if( ttCat_==2 ) mySample_sampleName_ = "ttb_powheg";
      if( ttCat_==3 ) mySample_sampleName_ = "tt2b_powheg";
      if( ttCat_==4 ) mySample_sampleName_ = "ttbb_powheg";
    }
  }
  else if( insample==2501 ){
    mySample_xSec_ = 831.76;//https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
    mySample_nGen_ = 10215131;//11158755;//25357774;//25446993;
    mySample_sampleName_ = "ttbar_madgraph";
    // standard
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115302/0000/");
    if( ttCat_>=0 ){
      if( ttCat_==0 ) mySample_sampleName_ = "ttlf_madgraph";
      if( ttCat_==1 ) mySample_sampleName_ = "ttcc_madgraph";
      if( ttCat_==2 ) mySample_sampleName_ = "ttb_madgraph";
      if( ttCat_==3 ) mySample_sampleName_ = "tt2b_madgraph";
      if( ttCat_==4 ) mySample_sampleName_ = "ttbb_madgraph";
    }
  }
  else if( insample==2502 ){
    mySample_xSec_ = 831.76;//https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
    mySample_nGen_ = 14188545;//42784971//25357774;//25446993; CHECK 19090400;//
    mySample_sampleName_ = "ttbar_amcatnlo";
    // standard
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/TT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115235/0000/");
    if( ttCat_>=0 ){
      if( ttCat_==0 ) mySample_sampleName_ = "ttlf_amcatnlo";
      if( ttCat_==1 ) mySample_sampleName_ = "ttcc_amcatnlo";
      if( ttCat_==2 ) mySample_sampleName_ = "ttb_amcatnlo";
      if( ttCat_==3 ) mySample_sampleName_ = "tt2b_amcatnlo";
      if( ttCat_==4 ) mySample_sampleName_ = "ttbb_amcatnlo";
    }
  }
  else if( insample==2530 ){
    mySample_xSec_ = 87.31;//https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
    mySample_nGen_ = 22533300;//42784971//25357774;//25446993;
    mySample_sampleName_ = "ttbar_powheg_dil";
    // standard
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/TTTo2L2Nu_13TeV-powheg/RunIIFall15DR76-RoWHepCloudTest-v14_MiniAOD76_v2_76X_mcRun2_asymptotic_v12_triggerTree_v1/160114_172830/0000/");
    if( ttCat_>=0 ){
      if( ttCat_==0 ) mySample_sampleName_ = "ttlf_powheg_dil";
      if( ttCat_==1 ) mySample_sampleName_ = "ttcc_powheg_dil";
      if( ttCat_==2 ) mySample_sampleName_ = "ttb_powheg_dil";
      if( ttCat_==3 ) mySample_sampleName_ = "tt2b_powheg_dil";
      if( ttCat_==4 ) mySample_sampleName_ = "ttbb_powheg_dil";
    }
  }
  else if( insample==2300 ){
    mySample_xSec_ = 6025.2; 
    mySample_nGen_ = 19310834;//59000;//28825132; CHECK 29193937;//
    mySample_sampleName_ = "ZJets_M50_amcatnlo";
    if( useHTbins_ ){
      mySample_xSec_ = 6025.2 * 9.59938e-01; // FIXME
      mySample_nGen_ = 19310834 * 9.59938e-01;// FIXME
      mySample_sampleName_ = "ZJets_M50_HT0to100_amcatnlo";
    }
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv1-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115327/0000/");
  }
  else if( insample==2301 ){
    mySample_xSec_ = 6025.2 * 3.02711e-02;//139.4 * 1.23; 
    mySample_nGen_ = 2448437;//2725655;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M50_HT100to200_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115508/0000/");
  }
  else if( insample==2302 ){
    mySample_xSec_ = 6025.2 * 8.26415e-03;//42.75 * 1.23; 
    mySample_nGen_ = 962195;//973937;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M50_HT200to400_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115536/0000/");
  }
  else if( insample==2303 ){
    mySample_xSec_ = 6025.2 * 1.12585e-03;//5.497 * 1.23; 
    mySample_nGen_ = 1069003;//1067758;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M50_HT400to600_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115604/0000/");
  }
  else if( insample==2304 ){
    mySample_xSec_ = 6025.2 * 4.01406e-04;//2.21 * 1.23; 
    mySample_nGen_ = 1031103;//998912;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M50_HT600toInf_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115629/0000/");
  }
  else if( insample==2305 ){
    mySample_xSec_ = 6025.2 * 9.59938e-01; 
    mySample_nGen_ = 9004328;//9042031;//59000;//28825132;
    mySample_sampleName_ = "ZJets_M50_madgraph";
    if( useHTbins_ ){
      mySample_xSec_ = 6025.2 * 9.59938e-01; // FIXME
      mySample_nGen_ = 9042031 * 9.59938e-01;// FIXME
      mySample_sampleName_ = "ZJets_M50_HT0to100_madgraph";
    }
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_114450/0000/");
  }
  else if( insample==2310 ){
    mySample_xSec_ = 18610;//22635.09; 
    mySample_nGen_ = 22217467;//19925500;//59000;//28825132; CHECK DEFINITELY WRONG 93034762;//
    mySample_sampleName_ = "ZJets_M10to50_amcatnlo";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_114357/0000/");
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1_triggerTree_v1/160223_114425/0000/");
  }
  else if( insample==2400 ){
    mySample_xSec_ = 61526.7;  
    mySample_nGen_ = 16518218; // CHECK 24156124;//
    mySample_sampleName_ = "WJets_amcatnlo";
    if( useHTbins_ ){
      mySample_xSec_ = 61526.7 * 9.648036e-01; // FIXME
      mySample_nGen_ = 16518218 * 9.648036e-01;// FIXME
      mySample_sampleName_ = "WJets_HT0To100_amcatnlo";
    }
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_114611/0000/");
  }
  else if( insample==2401 ){
    mySample_xSec_ = 1345 * 1.21;  
    mySample_nGen_ = 10205377;//10152718;
    mySample_sampleName_ = "WJets_HT100To200_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_114905/0000/");
  }
  else if( insample==2402 ){
    mySample_xSec_ = 359.7 * 1.21;  
    mySample_nGen_ = 4949568;//5221599;
    mySample_sampleName_ = "WJets_HT200To400_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_114932/0000/");
  }
  else if( insample==2403 ){
    mySample_xSec_ = 48.91 * 1.21;  
    mySample_nGen_ = 1943664;//1745914;
    mySample_sampleName_ = "WJets_HT400To600_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_114957/0000/");
  }
  else if( insample==2404 ){
    mySample_xSec_ = 18.77 * 1.21;  
    mySample_nGen_ = 1041358;//1039152;
    mySample_sampleName_ = "WJets_HT600ToInf_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115026/0000/");
  }
  else if( insample==2405 ){
    mySample_xSec_ = 61526.7;  
    mySample_nGen_ = 72207128;
    mySample_sampleName_ = "WJets_madgraph";
    if( useHTbins_ ){
      mySample_xSec_ = 61526.7 * 9.648036e-01; // FIXME
      mySample_nGen_ = 72207128 * 9.648036e-01;// FIXME
      mySample_sampleName_ = "WJets_HT0To100_madgraph";
    }
    /////mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/puigh/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151119_132138/0000/");
  }
  else if( insample==2524 ){
    mySample_xSec_ = 0.435;  
    mySample_nGen_ = 430330; // CHECK 831847;//
    mySample_sampleName_ = "ttW_had_amcatnlo";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115121/0000/");
  }
  else if( insample==2525 ){
    mySample_xSec_ = 0.21;  
    mySample_nGen_ = 129850; // CHECK 250307;//
    mySample_sampleName_ = "ttW_lep_amcatnlo";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115202/0000/");
  }
  else if( insample==2523 ){
    mySample_xSec_ = 0.611;  
    mySample_nGen_ = 351398; // CHECK 747000;//
    mySample_sampleName_ = "ttZ_had_amcatnlo";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115055/0000/");
  }
  else if( insample==2522 ){
    mySample_xSec_ = 0.2529;
    mySample_nGen_ = 184990;
    mySample_sampleName_ = "ttZ_lep_amcatnlo";
    //mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/puigh/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151115_164931/0000/");
  }
  else if( insample==2510 ){
    mySample_xSec_ = 136.02 * 1./3.; 
    mySample_nGen_ = 3299200;//3299800;
    mySample_sampleName_ = "st_tchan_powheg";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115421/0000/");
  }
  else if( insample==2511 ){
    mySample_xSec_ = 80.95 * 1./3.;  
    mySample_nGen_ = 1630900;//1695400;
    mySample_sampleName_ = "stbar_tchan_powheg";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115354/0000/");
  }
  else if( insample==2512 ){
    mySample_xSec_ = 35.85;  
    mySample_nGen_ = 995600;//995600;
    mySample_sampleName_ = "st_tWchan_powheg";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_114644/0000/");
  }
  else if( insample==2513 ){
    mySample_xSec_ = 35.85;  
    mySample_nGen_ = 999400;//1000000;
    mySample_sampleName_ = "stbar_tWchan_powheg";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_114718/0000/");
  }
  else if( insample==2514 ){
    mySample_xSec_ = 10.32 * 1./3.;  
    mySample_nGen_ = 613384; // CHECK 998400;//
    mySample_sampleName_ = "st_schan_amcatnlo";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2_triggerTree_v1/160223_115444/0000/");
  }
  else if( insample==2424 ){
    mySample_xSec_ = 118.7;  
    mySample_nGen_ = 988418;//994416;
    mySample_sampleName_ = "WW_pythia8";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/WW_TuneCUETP8M1_13TeV-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_114745/0000/");
  }
  else if( insample==2423 ){
    mySample_xSec_ = 44.9;  
    mySample_nGen_ = 1000000;//991232;
    mySample_sampleName_ = "WZ_pythia8";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/WZ_TuneCUETP8M1_13TeV-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_114815/0000/");
  }
  else if( insample==2323 ){
    mySample_xSec_ = 15.4;  
    mySample_nGen_ = 985600;//996168;
    mySample_sampleName_ = "ZZ_pythia8";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_114839/0000/");
  }
  else if( insample==1001 ){
    mySample_xSec_ = 27850000;  
    mySample_nGen_ = 81637494;
    mySample_sampleName_ = "QCD_HT100to200_madgraph";
    //mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/puigh/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151210_174712/0000/");
  }
  else if( insample==1002 ){
    mySample_xSec_ = 1717000;  
    mySample_nGen_ = 18718905;
    mySample_sampleName_ = "QCD_HT200to300_madgraph";
    //mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/puigh/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1_triggerTree_v1/151210_175438/0000/");
  }
  else if( insample==1003 ){
    mySample_xSec_ = 351300;  
    mySample_nGen_ = 19826197;
    mySample_sampleName_ = "QCD_HT300to500_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115832/0000/");
  }
  else if( insample==1004 ){
    mySample_xSec_ = 31630;  
    mySample_nGen_ = 19665695;//19664159;
    mySample_sampleName_ = "QCD_HT500to700_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115657/0000/");
  }
  else if( insample==1005 ){
    mySample_xSec_ = 6802;  
    mySample_nGen_ = 15356448;
    mySample_sampleName_ = "QCD_HT700to1000_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115903/0000/");
  }
  else if( insample==1006 ){
    mySample_xSec_ = 1206;  
    mySample_nGen_ = 4963895;
    mySample_sampleName_ = "QCD_HT1000to1500_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_triggerTree_QCD_HT1000to1500_13TeV_RunIIFall15MiniAODv2_2016_02_23/160223_115931/0000/");
  }
  else if( insample==1007 ){
    mySample_xSec_ = 120.4;  
    mySample_nGen_ = 3939077;//3868886;
    mySample_sampleName_ = "QCD_HT1500to2000_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115736/0000/");
  }
  else if( insample==1008 ){
    mySample_xSec_ = 25.24;  
    mySample_nGen_ = 1981228;//1912529;
    mySample_sampleName_ = "QCD_HT2000toInf_madgraph";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_115806/0000/");
  }
  else if( insample==9125 ){
    mySample_xSec_ = 0.5071 * 0.577;// YR4 * BR(all)  
    mySample_nGen_ = 3772012;//3933404;//199000;
    mySample_sampleName_ = "ttHTobb_powheg";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/ttHTobb_M125_13TeV_powheg_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_114522/0000/");
  }
  else if( insample==8125 ){
    mySample_xSec_ = 0.5071 * (1 - 0.577);// YR3 * BR(all)  
    mySample_nGen_ = 3945824;//3796398;//199000;
    mySample_sampleName_ = "ttHnonbb_powheg";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/ttHToNonbb_M125_13TeV_powheg_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_triggerTree_v1/160223_114546/0000/");
  }
  else if( insample==-11 ){
    mySample_sampleName_ = "SingleElectron_Run2015D_16Dec2015";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/SingleElectron/Run2015D-16Dec2015-v1_SilverJSON_2016_01_31_triggerTree_v1/160223_114019/0000/");
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/SingleElectron/Run2015D-16Dec2015-v1_SilverJSON_2016_01_31_triggerTree_v1/160223_114019/0001/");
  }
  else if( insample==-13 ){
    mySample_sampleName_ = "SingleMuon_Run2015D_16Dec2015";
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/SingleMuon/Run2015D-16Dec2015-v1_SilverJSON_2016_01_31_triggerTree_v1/160223_114100/0000/");
    mySample_inputDirs_.push_back("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/SingleMuon/Run2015D-16Dec2015-v1_SilverJSON_2016_01_31_triggerTree_v1/160223_114100/0001/");
  }



  std::string s_end = "histo_" + str_jobN + ".root";
  if( Njobs==1 ) s_end = "histo.root";

  std::string histofilename = Form("HistoFiles/ttHbb_data2mc_treeReader_%s_%s", mySample_sampleName_.c_str(), s_end.c_str());
  if( useCondor_ ) histofilename = Form("ttHbb_data2mc_treeReader_%s_%s", mySample_sampleName_.c_str(), s_end.c_str());


  TChain *chain = new TChain("triggeranalzyer/triggerTree");
  for( int iFile=0; iFile<int(mySample_inputDirs_.size()); iFile++ ){
    std::string treefilename = mySample_inputDirs_[iFile] + "trigger_analyzer*.root";
    std::cout << "  treefilename " << iFile << ": " << treefilename.c_str() << std::endl;
    //if( insample!=2500 ) chain->Add(treefilename.c_str());
    chain->Add(treefilename.c_str());
  }

  //chain->Add("/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/TriggerRun2/TriggerAnalyzer/trigger_analyzer.root");
  //if( insample==2500 ) chain->Add("root://cmseos.fnal.gov//store/user/tth/TriggerAnalyzer_triggerTree_76x/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1_triggerTree_v1/160201_050346/0000/trigger_analyzer_22*.root");

  std::string csvefffilename = Form("/uscms_data/d2/dpuigh/TTH/triggerRun2/CMSSW_7_4_15/src/TriggerRun2/TriggerAnalyzer/backup_2015_12_12_HistoFiles/ttHbb_data2mc_treeReader_%s_histo.root", mySample_sampleName_.c_str());

  csvefffilename = "/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/TriggerRun2/TriggerAnalyzer/HistoFiles/ttHbb_data2mc_treeReader_ttbar_histoEff.root";

  std::cout << "  csv eff filename = " << csvefffilename.c_str() << std::endl;

  TFile* f_CSVeff = new TFile(csvefffilename.c_str());

  fillCSVEffhistos(f_CSVeff);


  TFile* f_ele_id_sf = new TFile("/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/TriggerRun2/TriggerAnalyzer/data/ScaleFactor_GsfElectronToRECO_passingTrigWP80.txt.egamma_SF2D.root");
  TFile* f_ele_iso_sf = new TFile("/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/TriggerRun2/TriggerAnalyzer/data/Isolation_SF.root");
  TFile* f_ele_trig_sf = new TFile("/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/TriggerRun2/TriggerAnalyzer/data/eleTrig_SF.root");

  TFile* f_mu_id_sf = new TFile("/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/TriggerRun2/TriggerAnalyzer/data/MuonID_Z_RunCD_Reco76X_Feb15.root");
  TFile* f_mu_iso_sf = new TFile("/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/TriggerRun2/TriggerAnalyzer/data/MuonIso_Z_RunCD_Reco76X_Feb15.root");
  TFile* f_mu_trig_sf = new TFile("/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/TriggerRun2/TriggerAnalyzer/data/SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root");


  TH2D* h_ele_id_sf = (TH2D*)f_ele_id_sf->Get("EGamma_SF2D"); //max pt 200, x=sc eta, y=pt
  TH2D* h_ele_iso_sf = (TH2D*)f_ele_iso_sf->Get("IsolationSF"); //max pt 500, x=eta, y=pt
  TH2D* h_ele_trig_sf = (TH2D*)f_ele_trig_sf->Get("h_eleTrig_SF"); //max pt 200, x=pt, y=eta

  TH2D* h_mu_id_sf = (TH2D*)f_mu_id_sf->Get("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio"); //max pt 120, x=abseta, y=pt
  TH2D* h_mu_iso_sf = (TH2D*)f_mu_iso_sf->Get("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio"); //max pt 120, x=abseta, y=pt
  TH2D* h_mu_trig1_sf = (TH2D*)f_mu_trig_sf->Get("runD_IsoMu20_OR_IsoTkMu20_HLTv4p2_PtEtaBins/abseta_pt_ratio"); //max pt 120, x=abseta, y=pt
  TH2D* h_mu_trig2_sf = (TH2D*)f_mu_trig_sf->Get("runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins/abseta_pt_ratio"); //max pt 120, x=abseta, y=pt

  double mu_trig_fraction_1 = 754.394 / (754.394 + 1908.010);
  double mu_trig_fraction_2 = 1908.010 / (754.394 + 1908.010);
  
  TH2D* h_mu_trig_sf = (TH2D*)h_mu_trig1_sf->Clone("h_mu_trig_sf");
  for( int iBinX=0; iBinX < h_mu_trig_sf->GetNbinsX(); iBinX++ ){
   for( int iBinY=0; iBinY < h_mu_trig_sf->GetNbinsY(); iBinY++ ){
     double mu_hlt_sf_1 = h_mu_trig1_sf->GetBinContent(iBinX+1,iBinY+1);
     double mu_hlt_sf_2 = h_mu_trig2_sf->GetBinContent(iBinX+1,iBinY+1);

     double mu_hlt_sf = mu_trig_fraction_1 * mu_hlt_sf_1 + mu_trig_fraction_2 * mu_hlt_sf_2;
     double err_mu_hlt_sf = 0.01;
     h_mu_trig_sf->SetBinContent(iBinX+1,iBinY+1,mu_hlt_sf);
     h_mu_trig_sf->SetBinError(iBinX+1,iBinY+1,err_mu_hlt_sf);
   }
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

  TH1D* h_lheHT = new TH1D("h_lheHT",";LHE H_{T}", 2500, 0, MaxLHEHT );
  TH1D* h_lheHT_wgt = new TH1D("h_lheHT_wgt",";LHE H_{T}", 2500, 0, MaxLHEHT );
  TH1D* h_lheHT_wgt_afterHT = new TH1D("h_lheHT_wgt_afterHT",";LHE H_{T}", 2500, 0, MaxLHEHT );


  TH1D* h_numEvents = new TH1D("h_numEvents",";Number of events", 2, 0, 2 );
  TH1D* h_numEvents_wgt = new TH1D("h_numEvents_wgt",";Predicted number of events", 6, 0, 6 );

  TH1D* h_numEvents_wgtCSV = new TH1D("h_numEvents_wgtCSV",";Predicted number of events", 6, 0, 6 );
  TH1D* h_numEvents_wgtCSV_hf = new TH1D("h_numEvents_wgtCSV_hf",";Predicted number of events", 6, 0, 6 );
  TH1D* h_numEvents_wgtCSV_lf = new TH1D("h_numEvents_wgtCSV_lf",";Predicted number of events", 6, 0, 6 );
  TH1D* h_numEvents_wgtCSV_cf = new TH1D("h_numEvents_wgtCSV_cf",";Predicted number of events", 6, 0, 6 );

  TH1D* h_numEvents_wgtCSV_partonFlavourNoPU = new TH1D("h_numEvents_wgtCSV_partonFlavourNoPU",";Predicted number of events", 6, 0, 6 );
  TH1D* h_numEvents_wgtCSV_hf_partonFlavourNoPU = new TH1D("h_numEvents_wgtCSV_hf_partonFlavourNoPU",";Predicted number of events", 6, 0, 6 );
  TH1D* h_numEvents_wgtCSV_lf_partonFlavourNoPU = new TH1D("h_numEvents_wgtCSV_lf_partonFlavourNoPU",";Predicted number of events", 6, 0, 6 );
  TH1D* h_numEvents_wgtCSV_cf_partonFlavourNoPU = new TH1D("h_numEvents_wgtCSV_cf_partonFlavourNoPU",";Predicted number of events", 6, 0, 6 );

  TH1D* h_numEvents_nowgtCSV = new TH1D("h_numEvents_nowgtCSV",";Predicted number of events", 6, 0, 6 );

  TH1D* h_numEvents_lheWeight = new TH1D("h_numEvents_lheWeight",";Number of events per LHE weight", 250, 0, 250 );

  
  TH1D* h_numTruePVs = new TH1D("h_numTruePVs",";Number of True PVs", 50, 0, 50 );

  TH1D* h_numPVs = new TH1D("h_numPVs",";Number of PVs", 50, 0-0.5, 50-0.5 );


  TH1D* h_additionalJetEventId = new TH1D("h_additionalJetEventId",";additionalJetEventId", 201, -100-0.5, 101-0.5 );


  TH1D* h_deltaR_jet_lep = new TH1D("h_deltaR_jet_lep",";#DeltaR(jet,lep)", 61, 0., 6.1 );


  TH1D* h_ll_diLepMass = new TH1D("h_ll_diLepMass",";M(lep,lep)", 100, 50, 150 );
  TH1D* h_ee_diLepMass = new TH1D("h_ee_diLepMass",";M(e,e)", 100, 50, 150 );
  TH1D* h_mm_diLepMass = new TH1D("h_mm_diLepMass",";M(#mu,#mu)", 100, 50, 150 );

  TH1D* h_em_diLepMass = new TH1D("h_em_diLepMass",";M(e,#mu)", 100, 50, 150 );


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

  std::vector<TString> sys_cat_labels;
  sys_cat_labels.push_back("");               //0
  sys_cat_labels.push_back("_lepIdSFUp");     //1
  sys_cat_labels.push_back("_lepIdSFDown");   //2
  sys_cat_labels.push_back("_PUUp");          //3
  sys_cat_labels.push_back("_PUDown");        //4
  sys_cat_labels.push_back("_JERUp");         //5
  sys_cat_labels.push_back("_JERDown");       //6
  sys_cat_labels.push_back("_JESUp");         //7
  sys_cat_labels.push_back("_JESDown");       //8
  sys_cat_labels.push_back("_CSVLFUp");       //9
  sys_cat_labels.push_back("_CSVLFDown");     //10
  sys_cat_labels.push_back("_CSVHFUp");       //11
  sys_cat_labels.push_back("_CSVHFDown");     //12
  sys_cat_labels.push_back("_CSVHFStats1Up");     //13
  sys_cat_labels.push_back("_CSVHFStats1Down");   //14
  sys_cat_labels.push_back("_CSVHFStats2Up");     //15
  sys_cat_labels.push_back("_CSVHFStats2Down");   //16
  sys_cat_labels.push_back("_CSVLFStats1Up");     //17
  sys_cat_labels.push_back("_CSVLFStats1Down");   //18
  sys_cat_labels.push_back("_CSVLFStats2Up");     //19
  sys_cat_labels.push_back("_CSVLFStats2Down");   //20
  sys_cat_labels.push_back("_CSVCFErr1Up");     //21
  sys_cat_labels.push_back("_CSVCFErr1Down");   //22
  sys_cat_labels.push_back("_CSVCFErr2Up");     //23
  sys_cat_labels.push_back("_CSVCFErr2Down");   //24
  sys_cat_labels.push_back("_muFUp");           //25
  sys_cat_labels.push_back("_muFDown");         //26
  sys_cat_labels.push_back("_muRUp");           //27
  sys_cat_labels.push_back("_muRDown");         //28
  sys_cat_labels.push_back("_muRmuFUp");        //29
  sys_cat_labels.push_back("_muRmuFDown");      //30
  sys_cat_labels.push_back("_NNPDFUp");        //31
  sys_cat_labels.push_back("_NNPDFDown");      //32


  int NumSysCat = int(sys_cat_labels.size());


  TH1D* h_numEvents_perSys = new TH1D("h_numEvents_perSys",";Systematic", NumSysCat, 0-0.5, NumSysCat-0.5 );
  TH1D* h_numEvents_perSys_wgtCSV = new TH1D("h_numEvents_perSys_wgtCSV",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );
  TH1D* h_numEvents_perSys_wgtCSV_4j2t = new TH1D("h_numEvents_perSys_wgtCSV_4j2t",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );
  TH1D* h_numEvents_perSys_wgtCSV_6j4t = new TH1D("h_numEvents_perSys_wgtCSV_6j4t",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );

  TH1D* h_numEvents_perSys_wgtCSV2 = new TH1D("h_numEvents_perSys_wgtCSV2",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );
  TH1D* h_numEvents_perSys_wgtCSV2_4j2t = new TH1D("h_numEvents_perSys_wgtCSV2_4j2t",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );
  TH1D* h_numEvents_perSys_wgtCSV2_6j4t = new TH1D("h_numEvents_perSys_wgtCSV2_6j4t",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );

  TH1D* h_numEvents_perSys_wgtCMVA = new TH1D("h_numEvents_perSys_wgtCMVA",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );
  TH1D* h_numEvents_perSys_wgtCMVA_4j2t = new TH1D("h_numEvents_perSys_wgtCMVA_4j2t",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );
  TH1D* h_numEvents_perSys_wgtCMVA_6j4t = new TH1D("h_numEvents_perSys_wgtCMVA_6j4t",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );

  TH1D* h_numEvents_perSys_wgtCMVA2 = new TH1D("h_numEvents_perSys_wgtCMVA2",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );
  TH1D* h_numEvents_perSys_wgtCMVA2_4j2t = new TH1D("h_numEvents_perSys_wgtCMVA2_4j2t",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );
  TH1D* h_numEvents_perSys_wgtCMVA2_6j4t = new TH1D("h_numEvents_perSys_wgtCMVA2_6j4t",";Systematic", NumSysCat+1, 0-0.5-1, NumSysCat-0.5 );

  
  TH2D* h_diff_wgtCSV_wgtCSV2_perSys = new TH2D("h_diff_wgtCSV_wgtCSV2_perSys",";Systematic;wgtCSV2 - wgtCSV", NumSysCat+1, 0-0.5-1, NumSysCat-0.5, 100, -1, 1 );
  TH2D* h_hf_diff_wgtCSV_wgtCSV2_perSys = new TH2D("h_hf_diff_wgtCSV_wgtCSV2_perSys",";Systematic;hf wgtCSV2 - wgtCSV", NumSysCat+1, 0-0.5-1, NumSysCat-0.5, 100, -1, 1 );
  TH2D* h_lf_diff_wgtCSV_wgtCSV2_perSys = new TH2D("h_lf_diff_wgtCSV_wgtCSV2_perSys",";Systematic;lf wgtCSV2 - wgtCSV", NumSysCat+1, 0-0.5-1, NumSysCat-0.5, 100, -1, 1 );
  TH2D* h_cf_diff_wgtCSV_wgtCSV2_perSys = new TH2D("h_cf_diff_wgtCSV_wgtCSV2_perSys",";Systematic;cf wgtCSV2 - wgtCSV", NumSysCat+1, 0-0.5-1, NumSysCat-0.5, 100, -1, 1 );

  TH2D* h_hf_wgtCSV_perSys = new TH2D("h_hf_wgtCSV_perSys",";Systematic;hf wgtCSV", NumSysCat+1, 0-0.5-1, NumSysCat-0.5, 100, 0., 3. );
  TH2D* h_lf_wgtCSV_perSys = new TH2D("h_lf_wgtCSV_perSys",";Systematic;lf wgtCSV", NumSysCat+1, 0-0.5-1, NumSysCat-0.5, 100, 0., 3. );
  TH2D* h_cf_wgtCSV_perSys = new TH2D("h_cf_wgtCSV_perSys",";Systematic;cf wgtCSV", NumSysCat+1, 0-0.5-1, NumSysCat-0.5, 100, 0., 3. );

  TH2D* h_hf_wgtCSV2_perSys = new TH2D("h_hf_wgtCSV2_perSys",";Systematic;hf wgtCSV2", NumSysCat+1, 0-0.5-1, NumSysCat-0.5, 100, 0., 3. );
  TH2D* h_lf_wgtCSV2_perSys = new TH2D("h_lf_wgtCSV2_perSys",";Systematic;lf wgtCSV2", NumSysCat+1, 0-0.5-1, NumSysCat-0.5, 100, 0., 3. );
  TH2D* h_cf_wgtCSV2_perSys = new TH2D("h_cf_wgtCSV2_perSys",";Systematic;cf wgtCSV2", NumSysCat+1, 0-0.5-1, NumSysCat-0.5, 100, 0., 3. );

  
  TH1D* h_numJet_nolepreq = new TH1D("h_numJet_nolepreq",";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
  TH1D* h_numJet_nolepreq_wgtCSV = new TH1D("h_numJet_nolepreq_wgtCSV",";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );


  int NumCuts = 5;
  TH1D* h_event_selection[NumSysCat];
  TH1D* h_ele_event_selection[NumSysCat];
  TH1D* h_mu_event_selection[NumSysCat];

  TH1D* h_event_selection_wgtCSV[NumSysCat];
  TH1D* h_ele_event_selection_wgtCSV[NumSysCat];
  TH1D* h_mu_event_selection_wgtCSV[NumSysCat];

  TH1D* h_numJet[NumSysCat];
  TH1D* h_numJet_wgtCSV[NumSysCat];
  TH1D* h_numJet_wgtCSV2[NumSysCat];
  TH1D* h_numJet_wgtCSV3[NumSysCat];
  TH1D* h_numJet_wgtCSV4[NumSysCat];

  TH1D* h_numBtag[NumSysCat];
  TH1D* h_numBtag_wgtCSV[NumSysCat];
  TH1D* h_numBtag_wgtCSV2[NumSysCat];
  TH1D* h_numBtag_wgtCSV3[NumSysCat];
  TH1D* h_numBtag_wgtCSV4[NumSysCat];

  TH1D* h_numJet_4j2t[NumSysCat];
  TH1D* h_numBtag_4j2t[NumSysCat];

  TH1D* h_numJet_wgtCSV_4j2t[NumSysCat];
  TH1D* h_numJet_wgtCSV2_4j2t[NumSysCat];
  TH1D* h_numJet_wgtCSV3_4j2t[NumSysCat];
  TH1D* h_numJet_wgtCSV4_4j2t[NumSysCat];
  TH1D* h_numJet_wgtCSV5_4j2t[NumSysCat];

  TH1D* h_numBtag_wgtCSV_4j2t[NumSysCat];

  TH1D* h_numBtag_4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV_4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV5_4j[NumSysCat];

  TH1D* h_numBtag_CMVA_4j[NumSysCat];
  TH1D* h_numBtag_CMVA_wgtCSV_4j[NumSysCat];
  TH1D* h_numBtag_CMVA_wgtCMVA_4j[NumSysCat];

  TH1D* h_numBtag_eq4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV_eq4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_eq4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_eq4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_eq4j[NumSysCat];



  TH1D* h_numBtag_1e_4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV_1e_4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_1e_4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_1e_4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_1e_4j[NumSysCat];

  TH1D* h_numBtag_1e_eq4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV_1e_eq4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_1e_eq4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_1e_eq4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_1e_eq4j[NumSysCat];

  TH1D* h_numBtag_1m_4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV_1m_4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_1m_4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_1m_4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_1m_4j[NumSysCat];

  TH1D* h_numBtag_1m_eq4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV_1m_eq4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_1m_eq4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_1m_eq4j[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_1m_eq4j[NumSysCat];


  TH1D* h_numBtag_4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV_4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_4j_met30[NumSysCat];

  TH1D* h_numBtag_eq4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV_eq4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_eq4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_eq4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_eq4j_met30[NumSysCat];

  TH1D* h_numBtag_1e_4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV_1e_4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_1e_4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_1e_4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_1e_4j_met30[NumSysCat];

  TH1D* h_numBtag_1e_eq4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV_1e_eq4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_1e_eq4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_1e_eq4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_1e_eq4j_met30[NumSysCat];

  TH1D* h_numBtag_1m_4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV_1m_4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_1m_4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_1m_4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_1m_4j_met30[NumSysCat];

  TH1D* h_numBtag_1m_eq4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV_1m_eq4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_1m_eq4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_1m_eq4j_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_1m_eq4j_met30[NumSysCat];


  TH1D* h_numBtag_4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV_4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_4j_met50[NumSysCat];

  TH1D* h_numBtag_eq4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV_eq4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_eq4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_eq4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_eq4j_met50[NumSysCat];

  TH1D* h_numBtag_1e_4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV_1e_4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_1e_4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_1e_4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_1e_4j_met50[NumSysCat];

  TH1D* h_numBtag_1e_eq4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV_1e_eq4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_1e_eq4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_1e_eq4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_1e_eq4j_met50[NumSysCat];

  TH1D* h_numBtag_1m_4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV_1m_4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_1m_4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_1m_4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_1m_4j_met50[NumSysCat];

  TH1D* h_numBtag_1m_eq4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV_1m_eq4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_1m_eq4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_1m_eq4j_met50[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_1m_eq4j_met50[NumSysCat];


  /// Dilepton - Z control region

  TH1D* h_numJet_2l[NumSysCat];
  TH1D* h_numJet_wgtCSV_2l[NumSysCat];
  TH1D* h_numJet_wgtCSV2_2l[NumSysCat];
  TH1D* h_numJet_wgtCSV3_2l[NumSysCat];
  TH1D* h_numJet_wgtCSV4_2l[NumSysCat];

  TH1D* h_numBtag_2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l[NumSysCat];

  TH1D* h_numJet_2l_met30[NumSysCat];
  TH1D* h_numJet_wgtCSV_2l_met30[NumSysCat];
  TH1D* h_numJet_wgtCSV2_2l_met30[NumSysCat];
  TH1D* h_numJet_wgtCSV3_2l_met30[NumSysCat];
  TH1D* h_numJet_wgtCSV4_2l_met30[NumSysCat];

  TH1D* h_numBtag_2l_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30[NumSysCat];


  ///// leq2j
  TH1D* h_numBtag_2l_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_leq2j[NumSysCat];

  TH1D* h_numBtag_2l_met30_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_leq2j[NumSysCat];

  // 0 jet
  TH1D* h_numBtag_2l_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_leq2j_0b0c0l[NumSysCat];

  // 1 jet
  TH1D* h_numBtag_2l_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_leq2j_1b0c0l[NumSysCat];

  TH1D* h_numBtag_2l_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_leq2j_0b1c0l[NumSysCat];

  TH1D* h_numBtag_2l_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_leq2j_0b0c1l[NumSysCat];

  // 2 jet
  TH1D* h_numBtag_2l_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_leq2j_2b0c0l[NumSysCat];

  TH1D* h_numBtag_2l_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_leq2j_0b2c0l[NumSysCat];

  TH1D* h_numBtag_2l_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_leq2j_0b0c2l[NumSysCat];

  TH1D* h_numBtag_2l_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_leq2j_1b1c0l[NumSysCat];

  TH1D* h_numBtag_2l_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_leq2j_1b0c1l[NumSysCat];

  TH1D* h_numBtag_2l_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_leq2j_0b1c1l[NumSysCat];



  //// met30

  // 0 jet
  TH1D* h_numBtag_2l_met30_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_leq2j_0b0c0l[NumSysCat];

  // 1 jet
  TH1D* h_numBtag_2l_met30_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_leq2j_1b0c0l[NumSysCat];

  TH1D* h_numBtag_2l_met30_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_leq2j_0b1c0l[NumSysCat];

  TH1D* h_numBtag_2l_met30_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_leq2j_0b0c1l[NumSysCat];

  // 2 jet
  TH1D* h_numBtag_2l_met30_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_leq2j_2b0c0l[NumSysCat];

  TH1D* h_numBtag_2l_met30_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_leq2j_0b2c0l[NumSysCat];

  TH1D* h_numBtag_2l_met30_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_leq2j_0b0c2l[NumSysCat];

  TH1D* h_numBtag_2l_met30_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_leq2j_1b1c0l[NumSysCat];

  TH1D* h_numBtag_2l_met30_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_leq2j_1b0c1l[NumSysCat];

  TH1D* h_numBtag_2l_met30_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_leq2j_0b1c1l[NumSysCat];

  //// new

  ///// geq1j_leq2j
  TH1D* h_numBtag_2l_geq1j_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_geq1j_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_geq1j_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_geq1j_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_geq1j_leq2j[NumSysCat];

  TH1D* h_numBtag_2l_met30_geq1j_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_geq1j_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j[NumSysCat];

  // 0 jet
  TH1D* h_numBtag_2l_geq1j_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_geq1j_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b0c0l[NumSysCat];

  // 1 jet
  TH1D* h_numBtag_2l_geq1j_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_geq1j_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_geq1j_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_geq1j_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_geq1j_leq2j_1b0c0l[NumSysCat];

  TH1D* h_numBtag_2l_geq1j_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_geq1j_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b1c0l[NumSysCat];

  TH1D* h_numBtag_2l_geq1j_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_geq1j_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b0c1l[NumSysCat];

  // 2 jet
  TH1D* h_numBtag_2l_geq1j_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_geq1j_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_geq1j_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_geq1j_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_geq1j_leq2j_2b0c0l[NumSysCat];

  TH1D* h_numBtag_2l_geq1j_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_geq1j_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b2c0l[NumSysCat];

  TH1D* h_numBtag_2l_geq1j_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_geq1j_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b0c2l[NumSysCat];

  TH1D* h_numBtag_2l_geq1j_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_geq1j_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_geq1j_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_geq1j_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_geq1j_leq2j_1b1c0l[NumSysCat];

  TH1D* h_numBtag_2l_geq1j_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_geq1j_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_geq1j_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_geq1j_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_geq1j_leq2j_1b0c1l[NumSysCat];

  TH1D* h_numBtag_2l_geq1j_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_geq1j_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b1c1l[NumSysCat];



  //// met30

  // 0 jet
  TH1D* h_numBtag_2l_met30_geq1j_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b0c0l[NumSysCat];

  // 1 jet
  TH1D* h_numBtag_2l_met30_geq1j_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_1b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_1b0c0l[NumSysCat];

  TH1D* h_numBtag_2l_met30_geq1j_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b1c0l[NumSysCat];

  TH1D* h_numBtag_2l_met30_geq1j_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b0c1l[NumSysCat];

  // 2 jet
  TH1D* h_numBtag_2l_met30_geq1j_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_2b0c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_2b0c0l[NumSysCat];

  TH1D* h_numBtag_2l_met30_geq1j_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b2c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b2c0l[NumSysCat];

  TH1D* h_numBtag_2l_met30_geq1j_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b0c2l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b0c2l[NumSysCat];

  TH1D* h_numBtag_2l_met30_geq1j_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_1b1c0l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_1b1c0l[NumSysCat];

  TH1D* h_numBtag_2l_met30_geq1j_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_1b0c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_1b0c1l[NumSysCat];

  TH1D* h_numBtag_2l_met30_geq1j_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b1c1l[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b1c1l[NumSysCat];


  // EleMu CR
  TH1D* h_numJet_em[NumSysCat];
  TH1D* h_numJet_wgtCSV_em[NumSysCat];
  TH1D* h_numJet_wgtCSV2_em[NumSysCat];
  TH1D* h_numJet_wgtCSV3_em[NumSysCat];
  TH1D* h_numJet_wgtCSV4_em[NumSysCat];

  TH1D* h_numBtag_em_2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV_em_2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV2_em_2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV3_em_2j[NumSysCat];
  TH1D* h_numBtag_wgtCSV4_em_2j[NumSysCat];


  //// new

  TH1D* h_category_yield[NumSysCat];
  TH1D* h_category_yield_1e[NumSysCat];
  TH1D* h_category_yield_1m[NumSysCat];

  TH1D* h_category_yield_wgtCSV[NumSysCat];
  TH1D* h_category_yield_wgtCSV_1e[NumSysCat];
  TH1D* h_category_yield_wgtCSV_1m[NumSysCat];

  TH1D* h_category_yield_wgtCSV2[NumSysCat];
  TH1D* h_category_yield_wgtCSV2_1e[NumSysCat];
  TH1D* h_category_yield_wgtCSV2_1m[NumSysCat];

  TH1D* h_category_yield_wgtCSV5[NumSysCat];
  TH1D* h_category_yield_wgtCSV5_1e[NumSysCat];
  TH1D* h_category_yield_wgtCSV5_1m[NumSysCat];

  TH1D* h_category_yield_cmva[NumSysCat];
  TH1D* h_category_yield_cmva_1e[NumSysCat];
  TH1D* h_category_yield_cmva_1m[NumSysCat];

  TH1D* h_category_yield_cmva_wgtCMVA[NumSysCat];
  TH1D* h_category_yield_cmva_wgtCMVA_1e[NumSysCat];
  TH1D* h_category_yield_cmva_wgtCMVA_1m[NumSysCat];

  TH1D* h_category_yield_cmva_wgtCMVA2[NumSysCat];
  TH1D* h_category_yield_cmva_wgtCMVA2_1e[NumSysCat];
  TH1D* h_category_yield_cmva_wgtCMVA2_1m[NumSysCat];

  TH1D* h_numPVs_wgt[NumSysCat];
  TH1D* h_numPVs_noPUwgt[NumSysCat];

  TH1D* h_PV_x_wgt[NumSysCat];
  TH1D* h_PV_y_wgt[NumSysCat];
  TH1D* h_PV_z_wgt[NumSysCat];

  TH1D* h_pfMETNoHF_pt_2l[NumSysCat];
  TH1D* h_pfMETNoHF_pt_2e[NumSysCat];
  TH1D* h_pfMETNoHF_pt_2m[NumSysCat];
  TH1D* h_pfMETNoHF_pt_em[NumSysCat];

  TH1D* h_pfMETNoHF_pt_4j_1l[NumSysCat];
  TH1D* h_pfMETNoHF_pt_4j_1e[NumSysCat];
  TH1D* h_pfMETNoHF_pt_4j_1m[NumSysCat];



  //// Histograming 
  double metmax   = 500.;
  double lepPtMax = 300.;
  double jetptmax = 500.;
  double htmax    = 2000.;
  int NcsvBins = 132;

  int NmetBins   = int( metmax/20. + 0.0001 );
  int NlepPtBins = int( lepPtMax/10. + 0.0001 );
  int NjetptBins = int( jetptmax/10. + 0.0001 );
  int NhtBins    = int( htmax/50. + 0.0001 );

  for( int iSys=0; iSys<NumSysCat; iSys++ ){

    TString suffix = sys_cat_labels[iSys];

    h_mu_event_selection[iSys] = new TH1D("h_mu_event_selection" + suffix,";cut", NumCuts, 0, NumCuts );
    h_mu_event_selection[iSys]->GetXaxis()->SetBinLabel(1,"All");
    h_mu_event_selection[iSys]->GetXaxis()->SetBinLabel(2,"HLT Mu");
    h_mu_event_selection[iSys]->GetXaxis()->SetBinLabel(3,"==1 muon");
    h_mu_event_selection[iSys]->GetXaxis()->SetBinLabel(4,">=4 jets");
    h_mu_event_selection[iSys]->GetXaxis()->SetBinLabel(5,">=2 b-jets");

    h_ele_event_selection[iSys] = new TH1D("h_ele_event_selection" + suffix,";cut", NumCuts, 0, NumCuts );
    h_ele_event_selection[iSys]->GetXaxis()->SetBinLabel(1,"All");
    h_ele_event_selection[iSys]->GetXaxis()->SetBinLabel(2,"HLT Ele");
    h_ele_event_selection[iSys]->GetXaxis()->SetBinLabel(3,"==1 electron");
    h_ele_event_selection[iSys]->GetXaxis()->SetBinLabel(4,">=4 jets");
    h_ele_event_selection[iSys]->GetXaxis()->SetBinLabel(5,">=2 b-jets");

    h_event_selection[iSys]  = new TH1D("h_event_selection" + suffix,";cut", NumCuts, 0, NumCuts );
    h_event_selection[iSys]->GetXaxis()->SetBinLabel(1,"All");
    h_event_selection[iSys]->GetXaxis()->SetBinLabel(2,"HLT");
    h_event_selection[iSys]->GetXaxis()->SetBinLabel(3,"==1 lepton");
    h_event_selection[iSys]->GetXaxis()->SetBinLabel(4,">=4 jets");
    h_event_selection[iSys]->GetXaxis()->SetBinLabel(5,">=2 b-jets");

    // wgtCSV selection histograms
    h_mu_event_selection_wgtCSV[iSys] = new TH1D("h_mu_event_selection_wgtCSV" + suffix,";cut", NumCuts, 0, NumCuts );
    h_mu_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(1,"All");
    h_mu_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(2,"HLT Mu");
    h_mu_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(3,"==1 muon");
    h_mu_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(4,">=4 jets");
    h_mu_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(5,">=2 b-jets");

    h_ele_event_selection_wgtCSV[iSys] = new TH1D("h_ele_event_selection_wgtCSV" + suffix,";cut", NumCuts, 0, NumCuts );
    h_ele_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(1,"All");
    h_ele_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(2,"HLT Ele");
    h_ele_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(3,"==1 electron");
    h_ele_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(4,">=4 jets");
    h_ele_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(5,">=2 b-jets");

    h_event_selection_wgtCSV[iSys]  = new TH1D("h_event_selection_wgtCSV" + suffix,";cut", NumCuts, 0, NumCuts );
    h_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(1,"All");
    h_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(2,"HLT");
    h_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(3,"==1 lepton");
    h_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(4,">=4 jets");
    h_event_selection_wgtCSV[iSys]->GetXaxis()->SetBinLabel(5,">=2 b-jets");


    h_numJet[iSys] = new TH1D("h_numJet" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV[iSys] = new TH1D("h_numJet_wgtCSV" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV2[iSys] = new TH1D("h_numJet_wgtCSV2" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV3[iSys] = new TH1D("h_numJet_wgtCSV3" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV4[iSys] = new TH1D("h_numJet_wgtCSV4" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );

    h_numBtag[iSys] = new TH1D("h_numBtag" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV[iSys] = new TH1D("h_numBtag_wgtCSV" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2[iSys] = new TH1D("h_numBtag_wgtCSV2" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3[iSys] = new TH1D("h_numBtag_wgtCSV3" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4[iSys] = new TH1D("h_numBtag_wgtCSV4" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );



    h_numJet_4j2t[iSys] = new TH1D("h_numJet_4j2t" + suffix,";Number of Jets", MaxNjet-4, 4-0.5, MaxNjet-0.5 );
    h_numBtag_4j2t[iSys] = new TH1D("h_numBtag_4j2t" + suffix,";Number of b-tagged Jets", MaxNbtag-2, 2-0.5, MaxNbtag-0.5 );

    h_numJet_wgtCSV_4j2t[iSys] = new TH1D("h_numJet_wgtCSV_4j2t" + suffix,";Number of Jets", MaxNjet-4, 4-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV2_4j2t[iSys] = new TH1D("h_numJet_wgtCSV2_4j2t" + suffix,";Number of Jets", MaxNjet-4, 4-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV3_4j2t[iSys] = new TH1D("h_numJet_wgtCSV3_4j2t" + suffix,";Number of Jets", MaxNjet-4, 4-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV4_4j2t[iSys] = new TH1D("h_numJet_wgtCSV4_4j2t" + suffix,";Number of Jets", MaxNjet-4, 4-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV5_4j2t[iSys] = new TH1D("h_numJet_wgtCSV5_4j2t" + suffix,";Number of Jets", MaxNjet-4, 4-0.5, MaxNjet-0.5 );

    h_numBtag_wgtCSV_4j2t[iSys] = new TH1D("h_numBtag_wgtCSV_4j2t" + suffix,";Number of b-tagged Jets", MaxNbtag-2, 2-0.5, MaxNbtag-0.5 );


    h_numBtag_4j[iSys] = new TH1D("h_numBtag_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_4j[iSys] = new TH1D("h_numBtag_wgtCSV_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_4j[iSys] = new TH1D("h_numBtag_wgtCSV2_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_4j[iSys] = new TH1D("h_numBtag_wgtCSV3_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_4j[iSys] = new TH1D("h_numBtag_wgtCSV4_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV5_4j[iSys] = new TH1D("h_numBtag_wgtCSV5_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );


    h_numBtag_CMVA_4j[iSys] = new TH1D("h_numBtag_CMVA_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_CMVA_wgtCSV_4j[iSys] = new TH1D("h_numBtag_CMVA_wgtCSV_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_CMVA_wgtCMVA_4j[iSys] = new TH1D("h_numBtag_CMVA_wgtCMVA_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numBtag_eq4j[iSys] = new TH1D("h_numBtag_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_eq4j[iSys] = new TH1D("h_numBtag_wgtCSV_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_eq4j[iSys] = new TH1D("h_numBtag_wgtCSV2_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_eq4j[iSys] = new TH1D("h_numBtag_wgtCSV3_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_eq4j[iSys] = new TH1D("h_numBtag_wgtCSV4_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );


    h_numBtag_1e_4j[iSys] = new TH1D("h_numBtag_1e_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_1e_4j[iSys] = new TH1D("h_numBtag_wgtCSV_1e_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_1e_4j[iSys] = new TH1D("h_numBtag_wgtCSV2_1e_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_1e_4j[iSys] = new TH1D("h_numBtag_wgtCSV3_1e_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_1e_4j[iSys] = new TH1D("h_numBtag_wgtCSV4_1e_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numBtag_1e_eq4j[iSys] = new TH1D("h_numBtag_1e_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_1e_eq4j[iSys] = new TH1D("h_numBtag_wgtCSV_1e_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_1e_eq4j[iSys] = new TH1D("h_numBtag_wgtCSV2_1e_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_1e_eq4j[iSys] = new TH1D("h_numBtag_wgtCSV3_1e_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_1e_eq4j[iSys] = new TH1D("h_numBtag_wgtCSV4_1e_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numBtag_1m_4j[iSys] = new TH1D("h_numBtag_1m_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_1m_4j[iSys] = new TH1D("h_numBtag_wgtCSV_1m_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_1m_4j[iSys] = new TH1D("h_numBtag_wgtCSV2_1m_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_1m_4j[iSys] = new TH1D("h_numBtag_wgtCSV3_1m_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_1m_4j[iSys] = new TH1D("h_numBtag_wgtCSV4_1m_4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numBtag_1m_eq4j[iSys] = new TH1D("h_numBtag_1m_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_1m_eq4j[iSys] = new TH1D("h_numBtag_wgtCSV_1m_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_1m_eq4j[iSys] = new TH1D("h_numBtag_wgtCSV2_1m_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_1m_eq4j[iSys] = new TH1D("h_numBtag_wgtCSV3_1m_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_1m_eq4j[iSys] = new TH1D("h_numBtag_wgtCSV4_1m_eq4j" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );



    h_numBtag_4j_met30[iSys] = new TH1D("h_numBtag_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV2_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV3_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV4_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numBtag_eq4j_met30[iSys] = new TH1D("h_numBtag_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_eq4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_eq4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV2_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_eq4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV3_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_eq4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV4_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );


    h_numBtag_1e_4j_met30[iSys] = new TH1D("h_numBtag_1e_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_1e_4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV_1e_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_1e_4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV2_1e_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_1e_4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV3_1e_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_1e_4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV4_1e_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numBtag_1e_eq4j_met30[iSys] = new TH1D("h_numBtag_1e_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_1e_eq4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV_1e_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_1e_eq4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV2_1e_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_1e_eq4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV3_1e_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_1e_eq4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV4_1e_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numBtag_1m_4j_met30[iSys] = new TH1D("h_numBtag_1m_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_1m_4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV_1m_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_1m_4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV2_1m_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_1m_4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV3_1m_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_1m_4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV4_1m_4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numBtag_1m_eq4j_met30[iSys] = new TH1D("h_numBtag_1m_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_1m_eq4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV_1m_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_1m_eq4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV2_1m_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_1m_eq4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV3_1m_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_1m_eq4j_met30[iSys] = new TH1D("h_numBtag_wgtCSV4_1m_eq4j_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );


    h_numBtag_4j_met50[iSys] = new TH1D("h_numBtag_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV2_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV3_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV4_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numBtag_eq4j_met50[iSys] = new TH1D("h_numBtag_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_eq4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_eq4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV2_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_eq4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV3_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_eq4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV4_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numBtag_1e_4j_met50[iSys] = new TH1D("h_numBtag_1e_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_1e_4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV_1e_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_1e_4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV2_1e_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_1e_4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV3_1e_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_1e_4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV4_1e_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numBtag_1e_eq4j_met50[iSys] = new TH1D("h_numBtag_1e_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_1e_eq4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV_1e_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_1e_eq4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV2_1e_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_1e_eq4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV3_1e_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_1e_eq4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV4_1e_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numBtag_1m_4j_met50[iSys] = new TH1D("h_numBtag_1m_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_1m_4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV_1m_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_1m_4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV2_1m_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_1m_4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV3_1m_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_1m_4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV4_1m_4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numBtag_1m_eq4j_met50[iSys] = new TH1D("h_numBtag_1m_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_1m_eq4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV_1m_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_1m_eq4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV2_1m_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_1m_eq4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV3_1m_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_1m_eq4j_met50[iSys] = new TH1D("h_numBtag_wgtCSV4_1m_eq4j_met50" + suffix,";Number of b-tagged Jets", MaxNbtag, -0.5, MaxNbtag-0.5 );

    h_numJet_em[iSys] = new TH1D("h_numJet_em" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV_em[iSys] = new TH1D("h_numJet_wgtCSV_em" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV2_em[iSys] = new TH1D("h_numJet_wgtCSV2_em" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV3_em[iSys] = new TH1D("h_numJet_wgtCSV3_em" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV4_em[iSys] = new TH1D("h_numJet_wgtCSV4_em" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );

    h_numBtag_em_2j[iSys] = new TH1D("h_numBtag_em_2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_em_2j[iSys] = new TH1D("h_numBtag_wgtCSV_em_2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_em_2j[iSys] = new TH1D("h_numBtag_wgtCSV2_em_2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_em_2j[iSys] = new TH1D("h_numBtag_wgtCSV3_em_2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_em_2j[iSys] = new TH1D("h_numBtag_wgtCSV4_em_2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );


    h_numJet_2l[iSys] = new TH1D("h_numJet_2l" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV_2l[iSys] = new TH1D("h_numJet_wgtCSV_2l" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV2_2l[iSys] = new TH1D("h_numJet_wgtCSV2_2l" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV3_2l[iSys] = new TH1D("h_numJet_wgtCSV3_2l" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV4_2l[iSys] = new TH1D("h_numJet_wgtCSV4_2l" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );

    h_numBtag_2l[iSys] = new TH1D("h_numBtag_2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l[iSys] = new TH1D("h_numBtag_wgtCSV_2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );


    h_numJet_2l_met30[iSys] = new TH1D("h_numJet_2l_met30" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV_2l_met30[iSys] = new TH1D("h_numJet_wgtCSV_2l_met30" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV2_2l_met30[iSys] = new TH1D("h_numJet_wgtCSV2_2l_met30" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV3_2l_met30[iSys] = new TH1D("h_numJet_wgtCSV3_2l_met30" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );
    h_numJet_wgtCSV4_2l_met30[iSys] = new TH1D("h_numJet_wgtCSV4_2l_met30" + suffix,";Number of Jets", MaxNjet, 0-0.5, MaxNjet-0.5 );

    h_numBtag_2l_met30[iSys] = new TH1D("h_numBtag_2l_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    /// leq2j
    h_numBtag_2l_leq2j[iSys] = new TH1D("h_numBtag_2l_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV_2l_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_leq2j[iSys] = new TH1D("h_numBtag_2l_met30_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    // 0 jet
    h_numBtag_2l_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_2l_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    // 1 jet
    h_numBtag_2l_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_2l_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_2l_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_2l_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    // 2 jet
    h_numBtag_2l_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_2l_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_2l_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_2l_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_2l_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_2l_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_2l_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );


    /// met30

    // 0 jet
    h_numBtag_2l_met30_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_2l_met30_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    // 1 jet
    h_numBtag_2l_met30_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_2l_met30_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_2l_met30_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_2l_met30_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    // 2 jet
    h_numBtag_2l_met30_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_2l_met30_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_2l_met30_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_2l_met30_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_2l_met30_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_2l_met30_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_2l_met30_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );


    /// new
    /// geq1j_leq2j
    h_numBtag_2l_geq1j_leq2j[iSys] = new TH1D("h_numBtag_2l_geq1j_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_geq1j_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV_2l_geq1j_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_geq1j_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_geq1j_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_geq1j_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_geq1j_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_geq1j_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_geq1j_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_geq1j_leq2j[iSys] = new TH1D("h_numBtag_2l_met30_geq1j_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_geq1j_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_geq1j_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    // 0 jet
    h_numBtag_2l_geq1j_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_2l_geq1j_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_geq1j_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_geq1j_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    // 1 jet
    h_numBtag_2l_geq1j_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_2l_geq1j_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_geq1j_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_geq1j_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_geq1j_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_geq1j_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_geq1j_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_geq1j_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_geq1j_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_geq1j_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_geq1j_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_2l_geq1j_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_geq1j_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_geq1j_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_geq1j_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_2l_geq1j_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_geq1j_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_geq1j_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    // 2 jet
    h_numBtag_2l_geq1j_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_2l_geq1j_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_geq1j_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_geq1j_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_geq1j_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_geq1j_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_geq1j_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_geq1j_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_geq1j_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_geq1j_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_geq1j_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_2l_geq1j_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_geq1j_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_geq1j_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_geq1j_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_geq1j_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_geq1j_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_geq1j_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_geq1j_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_geq1j_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_geq1j_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_2l_geq1j_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_geq1j_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_geq1j_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_geq1j_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_geq1j_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_geq1j_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_geq1j_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_geq1j_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_geq1j_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_geq1j_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_2l_geq1j_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_geq1j_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_geq1j_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_geq1j_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_2l_geq1j_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_geq1j_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_geq1j_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_geq1j_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_2l_geq1j_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_geq1j_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_geq1j_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );


    /// met30

    // 0 jet
    h_numBtag_2l_met30_geq1j_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_2l_met30_geq1j_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    // 1 jet
    h_numBtag_2l_met30_geq1j_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_2l_met30_geq1j_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_1b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_1b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_geq1j_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_2l_met30_geq1j_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_geq1j_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_2l_met30_geq1j_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    // 2 jet
    h_numBtag_2l_met30_geq1j_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_2l_met30_geq1j_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_2b0c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_2b0c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_geq1j_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_2l_met30_geq1j_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_1b1c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_1b1c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_geq1j_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_2l_met30_geq1j_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_1b0c1l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_1b0c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_geq1j_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_2l_met30_geq1j_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b2c0l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b2c0l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_geq1j_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_2l_met30_geq1j_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b1c1l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b1c1l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    h_numBtag_2l_met30_geq1j_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_2l_met30_geq1j_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );
    h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b0c2l[iSys] = new TH1D("h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b0c2l" + suffix,";Number of b-tagged Jets", MaxNbtag, 0-0.5, MaxNbtag-0.5 );

    //// new


    h_category_yield[iSys] = new TH1D("h_category_yield" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_1e[iSys] = new TH1D("h_category_yield_1e" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_1m[iSys] = new TH1D("h_category_yield_1m" + suffix, ";category", NumCat, 0, NumCat );

    h_category_yield_wgtCSV[iSys] = new TH1D("h_category_yield_wgtCSV" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_wgtCSV_1e[iSys] = new TH1D("h_category_yield_wgtCSV_1e" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_wgtCSV_1m[iSys] = new TH1D("h_category_yield_wgtCSV_1m" + suffix, ";category", NumCat, 0, NumCat );

    h_category_yield_wgtCSV2[iSys] = new TH1D("h_category_yield_wgtCSV2" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_wgtCSV2_1e[iSys] = new TH1D("h_category_yield_wgtCSV2_1e" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_wgtCSV2_1m[iSys] = new TH1D("h_category_yield_wgtCSV2_1m" + suffix, ";category", NumCat, 0, NumCat );

    h_category_yield_wgtCSV5[iSys] = new TH1D("h_category_yield_wgtCSV5" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_wgtCSV5_1e[iSys] = new TH1D("h_category_yield_wgtCSV5_1e" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_wgtCSV5_1m[iSys] = new TH1D("h_category_yield_wgtCSV5_1m" + suffix, ";category", NumCat, 0, NumCat );


    

    h_category_yield_cmva[iSys] = new TH1D("h_category_yield_cmva" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_cmva_1e[iSys] = new TH1D("h_category_yield_cmva_1e" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_cmva_1m[iSys] = new TH1D("h_category_yield_cmva_1m" + suffix, ";category", NumCat, 0, NumCat );

    h_category_yield_cmva_wgtCMVA[iSys] = new TH1D("h_category_yield_cmva_wgtCMVA" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_cmva_wgtCMVA_1e[iSys] = new TH1D("h_category_yield_cmva_wgtCMVA_1e" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_cmva_wgtCMVA_1m[iSys] = new TH1D("h_category_yield_cmva_wgtCMVA_1m" + suffix, ";category", NumCat, 0, NumCat );

    h_category_yield_cmva_wgtCMVA2[iSys] = new TH1D("h_category_yield_cmva_wgtCMVA2" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_cmva_wgtCMVA2_1e[iSys] = new TH1D("h_category_yield_cmva_wgtCMVA2_1e" + suffix, ";category", NumCat, 0, NumCat );
    h_category_yield_cmva_wgtCMVA2_1m[iSys] = new TH1D("h_category_yield_cmva_wgtCMVA2_1m" + suffix, ";category", NumCat, 0, NumCat );

    for( int iCat=0; iCat<NumCat; iCat++ ){
      h_category_yield[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_1e[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_1m[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);

      h_category_yield_wgtCSV[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_wgtCSV_1e[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_wgtCSV_1m[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);

      h_category_yield_wgtCSV2[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_wgtCSV2_1e[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_wgtCSV2_1m[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);

      h_category_yield_wgtCSV5[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_wgtCSV5_1e[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_wgtCSV5_1m[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);


      h_category_yield_cmva[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_cmva_1e[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_cmva_1m[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);

      h_category_yield_cmva_wgtCMVA[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_cmva_wgtCMVA_1e[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_cmva_wgtCMVA_1m[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);

      h_category_yield_cmva_wgtCMVA2[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_cmva_wgtCMVA2_1e[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);
      h_category_yield_cmva_wgtCMVA2_1m[iSys]->GetXaxis()->SetBinLabel(iCat+1,cat_labels[iCat]);

    }


    h_numPVs_wgt[iSys] = new TH1D("h_numPVs_wgt" + suffix,";Number of PVs", 50, 0-0.5, 50-0.5 );
    h_numPVs_noPUwgt[iSys] = new TH1D("h_numPVs_noPUwgt" + suffix,";Number of PVs", 50, 0-0.5, 50-0.5 );

    h_PV_x_wgt[iSys] = new TH1D("h_PV_x_wgt" + suffix,";PV x", 100, 0.06, 0.14 );
    h_PV_y_wgt[iSys] = new TH1D("h_PV_y_wgt" + suffix,";PV y", 100, 0.14, 0.20 );
    h_PV_z_wgt[iSys] = new TH1D("h_PV_z_wgt" + suffix,";PV z", 100, -25, 25 );

    h_pfMETNoHF_pt_2l[iSys]  = new TH1D("h_pfMETNoHF_pt_2l" + suffix,";pfMETNoHF p_{T}", int(metmax+0.001), 0, metmax );
    h_pfMETNoHF_pt_2e[iSys]  = new TH1D("h_pfMETNoHF_pt_2e" + suffix,";pfMETNoHF p_{T}", int(metmax+0.001), 0, metmax );
    h_pfMETNoHF_pt_2m[iSys]  = new TH1D("h_pfMETNoHF_pt_2m" + suffix,";pfMETNoHF p_{T}", int(metmax+0.001), 0, metmax );

    h_pfMETNoHF_pt_em[iSys]  = new TH1D("h_pfMETNoHF_pt_em" + suffix,";pfMETNoHF p_{T}", int(metmax+0.001), 0, metmax );

    h_pfMETNoHF_pt_4j_1l[iSys]  = new TH1D("h_pfMETNoHF_pt_4j_1l" + suffix,";pfMETNoHF p_{T}", int(metmax+0.001), 0, metmax );
    h_pfMETNoHF_pt_4j_1e[iSys]  = new TH1D("h_pfMETNoHF_pt_4j_1e" + suffix,";pfMETNoHF p_{T}", int(metmax+0.001), 0, metmax );
    h_pfMETNoHF_pt_4j_1m[iSys]  = new TH1D("h_pfMETNoHF_pt_4j_1m" + suffix,";pfMETNoHF p_{T}", int(metmax+0.001), 0, metmax );
  } // end loop over systematics




  double HTmax = 1000.;
  int numHTbins = 1000;
  double L1HTmax = HTmax;
  int numL1HTTbins = numHTbins;





  std::vector<int> csv_systematics;
  csv_systematics.push_back(0);
  csv_systematics.push_back(7);
  csv_systematics.push_back(8);
  csv_systematics.push_back(9);
  csv_systematics.push_back(10);
  csv_systematics.push_back(11);
  csv_systematics.push_back(12);
  csv_systematics.push_back(13);
  csv_systematics.push_back(14);
  csv_systematics.push_back(15);
  csv_systematics.push_back(16);
  csv_systematics.push_back(17);
  csv_systematics.push_back(18);
  csv_systematics.push_back(19);
  csv_systematics.push_back(20);
  csv_systematics.push_back(21);
  csv_systematics.push_back(22);
  csv_systematics.push_back(23);
  csv_systematics.push_back(24);

  int numCSVSys = int(csv_systematics.size());


  TH1D* h_numEvents_Sys = new TH1D("h_numEvents_Sys",";systematic;predicted number of events", csv_systematics[numCSVSys-1]+1, 0, csv_systematics[numCSVSys-1]+1 );
  TH1D* h_numEvents_4j2t_Sys = new TH1D("h_numEvents_4j2t_Sys",";systematic;predicted number of events", csv_systematics[numCSVSys-1]+1, 0, csv_systematics[numCSVSys-1]+1 );
  TH1D* h_numEvents_4j_Sys = new TH1D("h_numEvents_4j_Sys",";systematic;predicted number of events", csv_systematics[numCSVSys-1]+1, 0, csv_systematics[numCSVSys-1]+1 );

  TH1D* h_numEvents_Sys_PU = new TH1D("h_numEvents_Sys_PU",";systematic;predicted number of events", csv_systematics[numCSVSys-1]+1, 0, csv_systematics[numCSVSys-1]+1 );

  
  TH1D* h_numEvents_jet30_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_jet30_nowgtCSV_WP[numCSVSys];

  TH1D* h_numEvents_jet20_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_jet20_nowgtCSV_WP[numCSVSys];

  TH1D* h_numEvents_jet30_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet30_nowgtCMVA_WP[numCSVSys];

  TH1D* h_numEvents_jet20_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet20_nowgtCMVA_WP[numCSVSys];


  TH1D* h_numEvents_jet30to60_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet30to60_nowgtCMVA_WP[numCSVSys];

  TH1D* h_numEvents_jet60to120_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet60to120_nowgtCMVA_WP[numCSVSys];

  TH1D* h_numEvents_jet120toInf_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet120toInf_nowgtCMVA_WP[numCSVSys];

  TH1D* h_wgt_csv_Sys[numCSVSys];
  TH1D* h_wgt_csv_hf_Sys[numCSVSys];
  TH1D* h_wgt_csv_lf_Sys[numCSVSys];
  TH1D* h_wgt_csv_cf_Sys[numCSVSys];

  TH1D* h_wgt_csv_noPU_Sys[numCSVSys];
  TH1D* h_wgt_csv_hf_noPU_Sys[numCSVSys];
  TH1D* h_wgt_csv_lf_noPU_Sys[numCSVSys];
  TH1D* h_wgt_csv_cf_noPU_Sys[numCSVSys];

  TH1D* h_wgt_csv_4j2t_Sys[numCSVSys];
  TH1D* h_wgt_csv_hf_4j2t_Sys[numCSVSys];
  TH1D* h_wgt_csv_lf_4j2t_Sys[numCSVSys];
  TH1D* h_wgt_csv_cf_4j2t_Sys[numCSVSys];

  TH1D* h_wgt_csv_4j_Sys[numCSVSys];
  TH1D* h_wgt_csv_hf_4j_Sys[numCSVSys];
  TH1D* h_wgt_csv_lf_4j_Sys[numCSVSys];
  TH1D* h_wgt_csv_cf_4j_Sys[numCSVSys];

  TH1D* h_wgt_csv_4j_noPU_Sys[numCSVSys];
  TH1D* h_wgt_csv_hf_4j_noPU_Sys[numCSVSys];
  TH1D* h_wgt_csv_lf_4j_noPU_Sys[numCSVSys];
  TH1D* h_wgt_csv_cf_4j_noPU_Sys[numCSVSys];

  TH1D* h_perjet_wgt_csv_Sys[5][3][numCSVSys];
  TH1D* h_perjet_wgt_csv_hf_Sys[5][1][numCSVSys];
  TH1D* h_perjet_wgt_csv_lf_Sys[4][3][numCSVSys];
  TH1D* h_perjet_wgt_csv_cf_Sys[5][1][numCSVSys];

  TH1D* h_perjet_wgt_cmva_Sys[5][3][numCSVSys];
  TH1D* h_perjet_wgt_cmva_hf_Sys[5][1][numCSVSys];
  TH1D* h_perjet_wgt_cmva_lf_Sys[4][3][numCSVSys];
  TH1D* h_perjet_wgt_cmva_cf_Sys[5][1][numCSVSys];

  for( int iSys=0; iSys<numCSVSys; iSys++ ){
    int useSys = csv_systematics[iSys];

    h_numEvents_jet30_wgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet30_wgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet30_nowgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet30_nowgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet20_wgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet20_wgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet20_nowgtCSV_WP[iSys] = new TH1D(Form("h_numEvents_jet20_nowgtCSV_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet30_wgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet30_wgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet30_nowgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet30_nowgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet20_wgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet20_wgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet20_nowgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet20_nowgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet30to60_wgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet30to60_wgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet30to60_nowgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet30to60_nowgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet60to120_wgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet60to120_wgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet60to120_nowgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet60to120_nowgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_numEvents_jet120toInf_wgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet120toInf_wgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );
    h_numEvents_jet120toInf_nowgtCMVA_WP[iSys] = new TH1D(Form("h_numEvents_jet120toInf_nowgtCMVA_WP_iSys_%d",useSys),";Predicted number of events", 3, 0, 3 );

    h_wgt_csv_Sys[iSys] = new TH1D(Form("h_wgt_csv_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_hf_Sys[iSys] = new TH1D(Form("h_wgt_csv_hf_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_lf_Sys[iSys] = new TH1D(Form("h_wgt_csv_lf_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_cf_Sys[iSys] = new TH1D(Form("h_wgt_csv_cf_Sys_%d",useSys),";weight", 300, 0., 3. );

    h_wgt_csv_noPU_Sys[iSys] = new TH1D(Form("h_wgt_csv_noPU_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_hf_noPU_Sys[iSys] = new TH1D(Form("h_wgt_csv_hf_noPU_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_lf_noPU_Sys[iSys] = new TH1D(Form("h_wgt_csv_lf_noPU_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_cf_noPU_Sys[iSys] = new TH1D(Form("h_wgt_csv_cf_noPU_Sys_%d",useSys),";weight", 300, 0., 3. );

    h_wgt_csv_4j2t_Sys[iSys] = new TH1D(Form("h_wgt_csv_4j2t_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_hf_4j2t_Sys[iSys] = new TH1D(Form("h_wgt_csv_hf_4j2t_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_lf_4j2t_Sys[iSys] = new TH1D(Form("h_wgt_csv_lf_4j2t_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_cf_4j2t_Sys[iSys] = new TH1D(Form("h_wgt_csv_cf_4j2t_Sys_%d",useSys),";weight", 300, 0., 3. );

    h_wgt_csv_4j_Sys[iSys] = new TH1D(Form("h_wgt_csv_4j_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_hf_4j_Sys[iSys] = new TH1D(Form("h_wgt_csv_hf_4j_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_lf_4j_Sys[iSys] = new TH1D(Form("h_wgt_csv_lf_4j_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_cf_4j_Sys[iSys] = new TH1D(Form("h_wgt_csv_cf_4j_Sys_%d",useSys),";weight", 300, 0., 3. );

    h_wgt_csv_4j_noPU_Sys[iSys] = new TH1D(Form("h_wgt_csv_4j_noPU_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_hf_4j_noPU_Sys[iSys] = new TH1D(Form("h_wgt_csv_hf_4j_noPU_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_lf_4j_noPU_Sys[iSys] = new TH1D(Form("h_wgt_csv_lf_4j_noPU_Sys_%d",useSys),";weight", 300, 0., 3. );
    h_wgt_csv_cf_4j_noPU_Sys[iSys] = new TH1D(Form("h_wgt_csv_cf_4j_noPU_Sys_%d",useSys),";weight", 300, 0., 3. );

    
    for( int iPt=0; iPt < 5; iPt++ ){
      for( int iEta=0; iEta < 3; iEta++ ){
	h_perjet_wgt_csv_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_csv_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
	if( iEta<1 ) h_perjet_wgt_csv_hf_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_csv_hf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
	if( iEta<1 ) h_perjet_wgt_csv_cf_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_csv_cf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
	if( iPt<4  ) h_perjet_wgt_csv_lf_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_csv_lf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );

	h_perjet_wgt_cmva_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_cmva_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
	if( iEta<1 ) h_perjet_wgt_cmva_hf_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_cmva_hf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
	if( iEta<1 ) h_perjet_wgt_cmva_cf_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_cmva_cf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
	if( iPt<4  ) h_perjet_wgt_cmva_lf_Sys[iPt][iEta][iSys] = new TH1D(Form("h_perjet_wgt_cmva_lf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys),";weight", 300, 0., 3. );
      }
    }
  }

  TH1D* h_event_wgt_csv = new TH1D("h_event_wgt_csv",";weight", 300, 0., 3. );
  TH1D* h_event_wgt_csv_hf = new TH1D("h_event_wgt_csv_hf",";weight", 300, 0., 3. );
  TH1D* h_event_wgt_csv_lf = new TH1D("h_event_wgt_csv_lf",";weight", 300, 0., 3. );
  TH1D* h_event_wgt_csv_cf = new TH1D("h_event_wgt_csv_cf",";weight", 300, 0., 3. );

  TH1D* h_event_wgt_csv_4j = new TH1D("h_event_wgt_csv_4j",";weight", 300, 0., 3. );
  TH1D* h_event_wgt_csv_4j_hf = new TH1D("h_event_wgt_csv_4j_hf",";weight", 300, 0., 3. );
  TH1D* h_event_wgt_csv_4j_lf = new TH1D("h_event_wgt_csv_4j_lf",";weight", 300, 0., 3. );
  TH1D* h_event_wgt_csv_4j_cf = new TH1D("h_event_wgt_csv_4j_cf",";weight", 300, 0., 3. );

  TH1D* h_event_wgt_csv_1l = new TH1D("h_event_wgt_csv_1l",";weight", 300, 0., 3. );
  TH1D* h_event_wgt_csv_1l_hf = new TH1D("h_event_wgt_csv_1l_hf",";weight", 300, 0., 3. );
  TH1D* h_event_wgt_csv_1l_lf = new TH1D("h_event_wgt_csv_1l_lf",";weight", 300, 0., 3. );
  TH1D* h_event_wgt_csv_1l_cf = new TH1D("h_event_wgt_csv_1l_cf",";weight", 300, 0., 3. );

  TH1D* h_event_wgt_csv_1l_4j = new TH1D("h_event_wgt_csv_1l_4j",";weight", 300, 0., 3. );
  TH1D* h_event_wgt_csv_1l_4j_hf = new TH1D("h_event_wgt_csv_1l_4j_hf",";weight", 300, 0., 3. );
  TH1D* h_event_wgt_csv_1l_4j_lf = new TH1D("h_event_wgt_csv_1l_4j_lf",";weight", 300, 0., 3. );
  TH1D* h_event_wgt_csv_1l_4j_cf = new TH1D("h_event_wgt_csv_1l_4j_cf",";weight", 300, 0., 3. );

  TH1D* h_perjet_wgt_csv[5][3];
  TH1D* h_perjet_wgt_csv_hf[5][1];
  TH1D* h_perjet_wgt_csv_lf[4][3];
  TH1D* h_perjet_wgt_csv_cf[5][1];

  TH1D* h_perjet_wgt_csv_lf_partonFlavourNoPU[4][3];

  TH1D* h_perjet_wgt_csv_4j[5][3];
  TH1D* h_perjet_wgt_csv_4j_hf[5][1];
  TH1D* h_perjet_wgt_csv_4j_lf[4][3];
  TH1D* h_perjet_wgt_csv_4j_cf[5][1];

  TH1D* h_perjet_wgt_csv_1l[5][3];
  TH1D* h_perjet_wgt_csv_1l_hf[5][1];
  TH1D* h_perjet_wgt_csv_1l_lf[4][3];
  TH1D* h_perjet_wgt_csv_1l_cf[5][1];

  TH1D* h_perjet_wgt_csv_1l_4j[5][3];
  TH1D* h_perjet_wgt_csv_1l_4j_hf[5][1];
  TH1D* h_perjet_wgt_csv_1l_4j_lf[4][3];
  TH1D* h_perjet_wgt_csv_1l_4j_cf[5][1];


  TH1D* h_perjet_csv_bPar[5][1];
  TH1D* h_perjet_csv_lfPar[4][3];
  TH1D* h_perjet_csv_udsgPar[4][3];
  TH1D* h_perjet_csv_oPar[4][3];
  TH1D* h_perjet_csv_cPar[5][1];

  TH1D* h_perjet_csv[5][3];
  TH1D* h_perjet_csv_hf[5][1];
  TH1D* h_perjet_csv_lf[4][3];
  TH1D* h_perjet_csv_cf[5][1];

  TH1D* h_perjet_csv_4j[5][3];
  TH1D* h_perjet_csv_4j_hf[5][1];
  TH1D* h_perjet_csv_4j_lf[4][3];
  TH1D* h_perjet_csv_4j_cf[5][1];

  TH1D* h_perjet_csv_1l[5][3];
  TH1D* h_perjet_csv_1l_hf[5][1];
  TH1D* h_perjet_csv_1l_lf[4][3];
  TH1D* h_perjet_csv_1l_cf[5][1];

  TH1D* h_perjet_csv_1l_4j[5][3];
  TH1D* h_perjet_csv_1l_4j_hf[5][1];
  TH1D* h_perjet_csv_1l_4j_lf[4][3];
  TH1D* h_perjet_csv_1l_4j_cf[5][1];

  
  TH1D* h_perjet_csv_PU[5][3];
  TH1D* h_perjet_csv_PU_hf[5][1];
  TH1D* h_perjet_csv_PU_lf[4][3];
  TH1D* h_perjet_csv_PU_cf[5][1];

  TH1D* h_perjet_csv_PU_4j[5][3];
  TH1D* h_perjet_csv_PU_4j_hf[5][1];
  TH1D* h_perjet_csv_PU_4j_lf[4][3];
  TH1D* h_perjet_csv_PU_4j_cf[5][1];

  TH1D* h_perjet_csv_PU_1l[5][3];
  TH1D* h_perjet_csv_PU_1l_hf[5][1];
  TH1D* h_perjet_csv_PU_1l_lf[4][3];
  TH1D* h_perjet_csv_PU_1l_cf[5][1];

  TH1D* h_perjet_csv_PU_1l_4j[5][3];
  TH1D* h_perjet_csv_PU_1l_4j_hf[5][1];
  TH1D* h_perjet_csv_PU_1l_4j_lf[4][3];
  TH1D* h_perjet_csv_PU_1l_4j_cf[5][1];

  for( int iPt=0; iPt < 5; iPt++ ){
    for( int iEta=0; iEta < 3; iEta++ ){
      h_perjet_wgt_csv[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );
      if( iEta<1 ) h_perjet_wgt_csv_hf[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_hf_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );
      if( iEta<1 ) h_perjet_wgt_csv_cf[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_cf_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );
      if( iPt<4  ) h_perjet_wgt_csv_lf[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_lf_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );

      if( iPt<4  ) h_perjet_wgt_csv_lf_partonFlavourNoPU[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_lf_partonFlavourNoPU_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );

      h_perjet_wgt_csv_4j[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_4j_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );
      if( iEta<1 ) h_perjet_wgt_csv_4j_hf[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_4j_hf_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );
      if( iEta<1 ) h_perjet_wgt_csv_4j_cf[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_4j_cf_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );
      if( iPt<4  ) h_perjet_wgt_csv_4j_lf[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_4j_lf_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );

      h_perjet_wgt_csv_1l[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_1l_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );
      if( iEta<1 ) h_perjet_wgt_csv_1l_hf[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_1l_hf_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );
      if( iEta<1 ) h_perjet_wgt_csv_1l_cf[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_1l_cf_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );
      if( iPt<4  ) h_perjet_wgt_csv_1l_lf[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_1l_lf_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );

      h_perjet_wgt_csv_1l_4j[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_1l_4j_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );
      if( iEta<1 ) h_perjet_wgt_csv_1l_4j_hf[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_1l_4j_hf_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );
      if( iEta<1 ) h_perjet_wgt_csv_1l_4j_cf[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_1l_4j_cf_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );
      if( iPt<4  ) h_perjet_wgt_csv_1l_4j_lf[iPt][iEta] = new TH1D(Form("h_perjet_wgt_csv_1l_4j_lf_iPt%d_iEta%d",iPt,iEta),";weight", 300, 0., 3. );

      ////
      h_perjet_csv[iPt][iEta] = new TH1D(Form("h_perjet_csv_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_hf[iPt][iEta] = new TH1D(Form("h_perjet_csv_hf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_cf[iPt][iEta] = new TH1D(Form("h_perjet_csv_cf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iPt<4  ) h_perjet_csv_lf[iPt][iEta] = new TH1D(Form("h_perjet_csv_lf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );

      if( iEta<1 ) h_perjet_csv_bPar[iPt][iEta] = new TH1D(Form("h_perjet_csv_bPar_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_cPar[iPt][iEta] = new TH1D(Form("h_perjet_csv_cPar_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iPt<4  ) h_perjet_csv_lfPar[iPt][iEta] = new TH1D(Form("h_perjet_csv_lfPar_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iPt<4  ) h_perjet_csv_udsgPar[iPt][iEta] = new TH1D(Form("h_perjet_csv_udsgPar_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iPt<4  ) h_perjet_csv_oPar[iPt][iEta] = new TH1D(Form("h_perjet_csv_oPar_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );

      h_perjet_csv_4j[iPt][iEta] = new TH1D(Form("h_perjet_csv_4j_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_4j_hf[iPt][iEta] = new TH1D(Form("h_perjet_csv_4j_hf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_4j_cf[iPt][iEta] = new TH1D(Form("h_perjet_csv_4j_cf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iPt<4  ) h_perjet_csv_4j_lf[iPt][iEta] = new TH1D(Form("h_perjet_csv_4j_lf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );

      h_perjet_csv_1l[iPt][iEta] = new TH1D(Form("h_perjet_csv_1l_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_1l_hf[iPt][iEta] = new TH1D(Form("h_perjet_csv_1l_hf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_1l_cf[iPt][iEta] = new TH1D(Form("h_perjet_csv_1l_cf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iPt<4  ) h_perjet_csv_1l_lf[iPt][iEta] = new TH1D(Form("h_perjet_csv_1l_lf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );

      h_perjet_csv_1l_4j[iPt][iEta] = new TH1D(Form("h_perjet_csv_1l_4j_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_1l_4j_hf[iPt][iEta] = new TH1D(Form("h_perjet_csv_1l_4j_hf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_1l_4j_cf[iPt][iEta] = new TH1D(Form("h_perjet_csv_1l_4j_cf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iPt<4  ) h_perjet_csv_1l_4j_lf[iPt][iEta] = new TH1D(Form("h_perjet_csv_1l_4j_lf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );

      ///
      h_perjet_csv_PU[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_PU_hf[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_hf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_PU_cf[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_cf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iPt<4  ) h_perjet_csv_PU_lf[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_lf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );

      h_perjet_csv_PU_4j[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_4j_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_PU_4j_hf[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_4j_hf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_PU_4j_cf[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_4j_cf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iPt<4  ) h_perjet_csv_PU_4j_lf[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_4j_lf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );

      h_perjet_csv_PU_1l[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_1l_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_PU_1l_hf[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_1l_hf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_PU_1l_cf[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_1l_cf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iPt<4  ) h_perjet_csv_PU_1l_lf[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_1l_lf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );

      h_perjet_csv_PU_1l_4j[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_1l_4j_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_PU_1l_4j_hf[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_1l_4j_hf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iEta<1 ) h_perjet_csv_PU_1l_4j_cf[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_1l_4j_cf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );
      if( iPt<4  ) h_perjet_csv_PU_1l_4j_lf[iPt][iEta] = new TH1D(Form("h_perjet_csv_PU_1l_4j_lf_iPt%d_iEta%d",iPt,iEta),";jet CSV", NcsvBins, -0.06, 1.01 );

    }
  }

  //
  // Lepton plots
  //
  TH1D* h_lepton_pt[NumCat][NumSysCat];
  TH1D* h_lepton_phi[NumCat][NumSysCat];
  TH1D* h_lepton_eta[NumCat][NumSysCat];
  TH1D* h_lepton_d0[NumCat][NumSysCat];
  TH1D* h_lepton_dZ[NumCat][NumSysCat];

  TH1D* h_electron_pt[NumCat][NumSysCat];
  TH1D* h_electron_phi[NumCat][NumSysCat];
  TH1D* h_electron_eta[NumCat][NumSysCat];
  TH1D* h_electron_relIso[NumCat][NumSysCat];
  TH1D* h_electron_trigMVAOutput[NumCat][NumSysCat];
  TH1D* h_electron_d0[NumCat][NumSysCat];
  TH1D* h_electron_dZ[NumCat][NumSysCat];

  TH1D* h_muon_pt[NumCat][NumSysCat];
  TH1D* h_muon_phi[NumCat][NumSysCat];
  TH1D* h_muon_eta[NumCat][NumSysCat];
  TH1D* h_muon_relIso[NumCat][NumSysCat];
  TH1D* h_muon_d0[NumCat][NumSysCat];
  TH1D* h_muon_dZ[NumCat][NumSysCat];

  TH1D* h_lepton_pt_wgtCSV[NumCat][NumSysCat];
  TH1D* h_lepton_phi_wgtCSV[NumCat][NumSysCat];
  TH1D* h_lepton_eta_wgtCSV[NumCat][NumSysCat];
  TH1D* h_lepton_d0_wgtCSV[NumCat][NumSysCat];
  TH1D* h_lepton_dZ_wgtCSV[NumCat][NumSysCat];

  TH1D* h_electron_pt_wgtCSV[NumCat][NumSysCat];
  TH1D* h_electron_phi_wgtCSV[NumCat][NumSysCat];
  TH1D* h_electron_eta_wgtCSV[NumCat][NumSysCat];
  TH1D* h_electron_relIso_wgtCSV[NumCat][NumSysCat];
  TH1D* h_electron_trigMVAOutput_wgtCSV[NumCat][NumSysCat];
  TH1D* h_electron_d0_wgtCSV[NumCat][NumSysCat];
  TH1D* h_electron_dZ_wgtCSV[NumCat][NumSysCat];

  TH1D* h_muon_pt_wgtCSV[NumCat][NumSysCat];
  TH1D* h_muon_phi_wgtCSV[NumCat][NumSysCat];
  TH1D* h_muon_eta_wgtCSV[NumCat][NumSysCat];
  TH1D* h_muon_relIso_wgtCSV[NumCat][NumSysCat];
  TH1D* h_muon_d0_wgtCSV[NumCat][NumSysCat];
  TH1D* h_muon_dZ_wgtCSV[NumCat][NumSysCat];


  // 
  // Jet plots
  //
  TH1D* h_jet_pt[NumCat][NumSysCat];
  TH1D* h_jet_eta[NumCat][NumSysCat];
  TH1D* h_jet_phi[NumCat][NumSysCat];
  TH1D* h_jet_csv[NumCat][NumSysCat];
  TH1D* h_jet_cmva[NumCat][NumSysCat];
  TH1D* h_jet_puMVA[NumCat][NumSysCat];

  TH1D* h_jet_pt_wgtCSV[NumCat][NumSysCat];
  TH1D* h_jet_eta_wgtCSV[NumCat][NumSysCat];
  TH1D* h_jet_phi_wgtCSV[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV[NumCat][NumSysCat];
  TH1D* h_jet_cmva_wgtCSV[NumCat][NumSysCat];
  TH1D* h_jet_puMVA_wgtCSV[NumCat][NumSysCat];

  TH1D* h_jet_pt_wgtCSV2[NumCat][NumSysCat];
  TH1D* h_jet_eta_wgtCSV2[NumCat][NumSysCat];
  TH1D* h_jet_phi_wgtCSV2[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV2[NumCat][NumSysCat];
  TH1D* h_jet_puMVA_wgtCSV2[NumCat][NumSysCat];

  TH1D* h_jet_pt_wgtCSV3[NumCat][NumSysCat];
  TH1D* h_jet_eta_wgtCSV3[NumCat][NumSysCat];
  TH1D* h_jet_phi_wgtCSV3[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV3[NumCat][NumSysCat];
  TH1D* h_jet_puMVA_wgtCSV3[NumCat][NumSysCat];

  TH1D* h_jet_pt_wgtCSV4[NumCat][NumSysCat];
  TH1D* h_jet_eta_wgtCSV4[NumCat][NumSysCat];
  TH1D* h_jet_phi_wgtCSV4[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV4[NumCat][NumSysCat];
  TH1D* h_jet_puMVA_wgtCSV4[NumCat][NumSysCat];

  TH1D* h_jet_pt_wgtCSV5[NumCat][NumSysCat];
  TH1D* h_jet_eta_wgtCSV5[NumCat][NumSysCat];
  TH1D* h_jet_phi_wgtCSV5[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV5[NumCat][NumSysCat];
  TH1D* h_jet_puMVA_wgtCSV5[NumCat][NumSysCat];


  TH1D* h_cmva_jet_pt[NumCat][NumSysCat];
  TH1D* h_cmva_jet_eta[NumCat][NumSysCat];
  TH1D* h_cmva_jet_phi[NumCat][NumSysCat];
  TH1D* h_cmva_jet_csv[NumCat][NumSysCat];
  TH1D* h_cmva_jet_cmva[NumCat][NumSysCat];
  TH1D* h_cmva_jet_puMVA[NumCat][NumSysCat];

  TH1D* h_cmva_jet_pt_wgtCMVA[NumCat][NumSysCat];
  TH1D* h_cmva_jet_eta_wgtCMVA[NumCat][NumSysCat];
  TH1D* h_cmva_jet_phi_wgtCMVA[NumCat][NumSysCat];
  TH1D* h_cmva_jet_csv_wgtCMVA[NumCat][NumSysCat];
  TH1D* h_cmva_jet_cmva_wgtCMVA[NumCat][NumSysCat];
  TH1D* h_cmva_jet_puMVA_wgtCMVA[NumCat][NumSysCat];

  
  TH1D* h_jet_csv_bFlav[NumCat][NumSysCat];
  TH1D* h_jet_csv_cFlav[NumCat][NumSysCat];
  TH1D* h_jet_csv_lFlav[NumCat][NumSysCat];

  TH1D* h_jet_csv_wgtCSV_bFlav[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV_cFlav[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV_lFlav[NumCat][NumSysCat];

  TH1D* h_jet_csv_wgtCSV2_bFlav[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV2_cFlav[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV2_lFlav[NumCat][NumSysCat];

  TH1D* h_jet_csv_wgtCSV3_bFlav[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV3_cFlav[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV3_lFlav[NumCat][NumSysCat];

  TH1D* h_jet_csv_wgtCSV4_bFlav[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV4_cFlav[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV4_lFlav[NumCat][NumSysCat];

  TH1D* h_jet_csv_wgtCSV5_bFlav[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV5_cFlav[NumCat][NumSysCat];
  TH1D* h_jet_csv_wgtCSV5_lFlav[NumCat][NumSysCat];

  
  TH1D* h_jet_cmva_bFlav[NumCat][NumSysCat];
  TH1D* h_jet_cmva_cFlav[NumCat][NumSysCat];
  TH1D* h_jet_cmva_lFlav[NumCat][NumSysCat];

  TH1D* h_jet_cmva_wgtCMVA_bFlav[NumCat][NumSysCat];
  TH1D* h_jet_cmva_wgtCMVA_cFlav[NumCat][NumSysCat];
  TH1D* h_jet_cmva_wgtCMVA_lFlav[NumCat][NumSysCat];


  //
  // Energy sum plots
  //
  TH1D* h_pfMET_pt[NumCat][NumSysCat];
  TH1D* h_pfMET_phi[NumCat][NumSysCat];
  TH1D* h_pfMETNoHF_pt[NumCat][NumSysCat];
  TH1D* h_pfMETNoHF_phi[NumCat][NumSysCat];
  TH1D* h_puppiMET_pt[NumCat][NumSysCat];
  TH1D* h_puppiMET_phi[NumCat][NumSysCat];
  TH1D* h_mht_pt[NumCat][NumSysCat];
  TH1D* h_mht_phi[NumCat][NumSysCat];

  TH1D* h_HT[NumCat][NumSysCat];

  TH1D* h_pfMET_pt_wgtCSV[NumCat][NumSysCat];
  TH1D* h_pfMET_phi_wgtCSV[NumCat][NumSysCat];


  std::vector<double> jet_pt_xbins;
  for( int iBin=0; iBin<200; iBin++ ) jet_pt_xbins.push_back(double(iBin));
  for( int iBin=200; iBin<400; iBin+=2 ) jet_pt_xbins.push_back(double(iBin));
  for( int iBin=400; iBin<500; iBin+=10 ) jet_pt_xbins.push_back(double(iBin));
  jet_pt_xbins.push_back(500);
  int Njetbins = int(jet_pt_xbins.size());

  double xbins_jetpt[311];
  for( int iBin=0; iBin<Njetbins; iBin++ ) xbins_jetpt[iBin] = jet_pt_xbins[iBin];

  TH2D* h_a_jet_pt_eta_csvM = new TH2D("h_a_jet_pt_eta_csvM",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );
  TH2D* h_b_jet_pt_eta_csvM = new TH2D("h_b_jet_pt_eta_csvM",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );
  TH2D* h_c_jet_pt_eta_csvM = new TH2D("h_c_jet_pt_eta_csvM",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );
  TH2D* h_l_jet_pt_eta_csvM = new TH2D("h_l_jet_pt_eta_csvM",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );

  TH2D* h_a_jet_pt_eta_all = new TH2D("h_a_jet_pt_eta_all",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );
  TH2D* h_b_jet_pt_eta_all = new TH2D("h_b_jet_pt_eta_all",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );
  TH2D* h_c_jet_pt_eta_all = new TH2D("h_c_jet_pt_eta_all",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );
  TH2D* h_l_jet_pt_eta_all = new TH2D("h_l_jet_pt_eta_all",";jet p_{T};jet |#eta|", Njetbins-1, xbins_jetpt, 25, 0, 2.5 );

  TH1D* h_jet_csv_b = new TH1D("h_jet_csv_b",";jet CSV", NcsvBins, -0.06, 1.01 );
  TH1D* h_jet_csv_c = new TH1D("h_jet_csv_c",";jet CSV", NcsvBins, -0.06, 1.01 );
  TH1D* h_jet_csv_l = new TH1D("h_jet_csv_l",";jet CSV", NcsvBins, -0.06, 1.01 );
  TH1D* h_jet_csv_a = new TH1D("h_jet_csv_a",";jet CSV", NcsvBins, -0.06, 1.01 );

  TH1D* h_jet_cmva_b = new TH1D("h_jet_cmva_b",";jet CMVA", NcsvBins, -1.01, 1.01 );
  TH1D* h_jet_cmva_c = new TH1D("h_jet_cmva_c",";jet CMVA", NcsvBins, -1.01, 1.01 );
  TH1D* h_jet_cmva_l = new TH1D("h_jet_cmva_l",";jet CMVA", NcsvBins, -1.01, 1.01 );
  TH1D* h_jet_cmva_a = new TH1D("h_jet_cmva_a",";jet CMVA", NcsvBins, -1.01, 1.01 );

  TH2D* h_b_csv_wgt = new TH2D("h_b_csv_wgt",";jet CSV;scale factor", NcsvBins, -0.06, 1.01, 200, 0, 2.0 );
  TH2D* h_b_csv_wgt2 = new TH2D("h_b_csv_wgt2",";jet CSV;scale factor", NcsvBins, -0.06, 1.01, 200, 0, 2.0 );
  TH2D* h_b_csv_diff_wgt2_wgt = new TH2D("h_b_csv_diff_wgt2_wgt",";jet CSV;scale factor difference", NcsvBins, -0.06, 1.01, 200, -1.0, 1.0 );

  TH2D* h_c_csv_wgt = new TH2D("h_c_csv_wgt",";jet CSV;scale factor", NcsvBins, -0.06, 1.01, 200, 0, 2.0 );
  TH2D* h_c_csv_wgt2 = new TH2D("h_c_csv_wgt2",";jet CSV;scale factor", NcsvBins, -0.06, 1.01, 200, 0, 2.0 );
  TH2D* h_c_csv_diff_wgt2_wgt = new TH2D("h_c_csv_diff_wgt2_wgt",";jet CSV;scale factor difference", NcsvBins, -0.06, 1.01, 200, -1.0, 1.0 );

  TH2D* h_l_csv_wgt = new TH2D("h_l_csv_wgt",";jet CSV;scale factor", NcsvBins, -0.06, 1.01, 200, 0, 2.0 );
  TH2D* h_l_csv_wgt2 = new TH2D("h_l_csv_wgt2",";jet CSV;scale factor", NcsvBins, -0.06, 1.01, 200, 0, 2.0 );
  TH2D* h_l_csv_diff_wgt2_wgt = new TH2D("h_l_csv_diff_wgt2_wgt",";jet CSV;scale factor difference", NcsvBins, -0.06, 1.01, 200, -1.0, 1.0 );


  for( int iSys=0; iSys<NumSysCat; iSys++ ){
    TString sys_suffix = sys_cat_labels[iSys];

    for( int iCat=0; iCat<NumCat; iCat++ ){
      TString cat_suffix = "_" + cat_labels[iCat];

      TString suffix = cat_suffix + sys_suffix;

      h_lepton_pt[iCat][iSys] = new TH1D("h_lepton_pt" + suffix,";lepton p_{T}", NlepPtBins, 0, lepPtMax );
      h_lepton_phi[iCat][iSys] = new TH1D("h_lepton_phi" + suffix,";lepton #phi", 34, -3.4, 3.4 );
      h_lepton_eta[iCat][iSys] = new TH1D("h_lepton_eta" + suffix,";lepton #eta", 25, -2.5, 2.5 );
      h_lepton_d0[iCat][iSys] = new TH1D("h_lepton_d0" + suffix,";lepton d0", 100, -0.05, 0.05 );
      h_lepton_dZ[iCat][iSys] = new TH1D("h_lepton_dZ" + suffix,";lepton dZ", 100, -0.05, 0.05 );

      h_electron_pt[iCat][iSys] = new TH1D("h_electron_pt" + suffix,";electron p_{T}", NlepPtBins, 0, lepPtMax );
      h_electron_phi[iCat][iSys] = new TH1D("h_electron_phi" + suffix,";electron #phi", 34, -3.4, 3.4 );
      h_electron_eta[iCat][iSys] = new TH1D("h_electron_eta" + suffix,";electron #eta", 25, -2.5, 2.5 );
      h_electron_relIso[iCat][iSys] = new TH1D("h_electron_relIso" + suffix,";electron relIso", 210, -0.01, 0.20 );
      h_electron_trigMVAOutput[iCat][iSys] = new TH1D("h_electron_trigMVAOutput" + suffix,";electron trigMVAOutput", 220, -1.1, 1.1 );
      h_electron_d0[iCat][iSys] = new TH1D("h_electron_d0" + suffix,";electron d0", 100, -0.05, 0.05 );
      h_electron_dZ[iCat][iSys] = new TH1D("h_electron_dZ" + suffix,";electron dZ", 100, -0.05, 0.05 );

      h_muon_pt[iCat][iSys] = new TH1D("h_muon_pt" + suffix,";muon p_{T}", NlepPtBins, 0, lepPtMax );
      h_muon_phi[iCat][iSys] = new TH1D("h_muon_phi" + suffix,";muon #phi", 34, -3.4, 3.4 );
      h_muon_eta[iCat][iSys] = new TH1D("h_muon_eta" + suffix,";muon #eta", 25, -2.5, 2.5 );
      h_muon_relIso[iCat][iSys] = new TH1D("h_muon_relIso" + suffix,";muon relIso", 210, -0.01, 0.20 );
      h_muon_d0[iCat][iSys] = new TH1D("h_muon_d0" + suffix,";muon d0", 100, -0.05, 0.05 );
      h_muon_dZ[iCat][iSys] = new TH1D("h_muon_dZ" + suffix,";muon dZ", 100, -0.05, 0.05 );

      
      h_lepton_pt_wgtCSV[iCat][iSys] = new TH1D("h_lepton_pt_wgtCSV" + suffix,";lepton p_{T}", NlepPtBins, 0, lepPtMax );
      h_lepton_phi_wgtCSV[iCat][iSys] = new TH1D("h_lepton_phi_wgtCSV" + suffix,";lepton #phi", 34, -3.4, 3.4 );
      h_lepton_eta_wgtCSV[iCat][iSys] = new TH1D("h_lepton_eta_wgtCSV" + suffix,";lepton #eta", 25, -2.5, 2.5 );
      h_lepton_d0_wgtCSV[iCat][iSys] = new TH1D("h_lepton_d0_wgtCSV" + suffix,";lepton d0", 100, -0.05, 0.05 );
      h_lepton_dZ_wgtCSV[iCat][iSys] = new TH1D("h_lepton_dZ_wgtCSV" + suffix,";lepton dZ", 100, -0.05, 0.05 );

      h_electron_pt_wgtCSV[iCat][iSys] = new TH1D("h_electron_pt_wgtCSV" + suffix,";electron p_{T}", NlepPtBins, 0, lepPtMax );
      h_electron_phi_wgtCSV[iCat][iSys] = new TH1D("h_electron_phi_wgtCSV" + suffix,";electron #phi", 34, -3.4, 3.4 );
      h_electron_eta_wgtCSV[iCat][iSys] = new TH1D("h_electron_eta_wgtCSV" + suffix,";electron #eta", 25, -2.5, 2.5 );
      h_electron_relIso_wgtCSV[iCat][iSys] = new TH1D("h_electron_relIso_wgtCSV" + suffix,";electron relIso", 210, -0.01, 0.20 );
      h_electron_trigMVAOutput_wgtCSV[iCat][iSys] = new TH1D("h_electron_trigMVAOutput_wgtCSV" + suffix,";electron trigMVAOutput", 220, -1.1, 1.1 );
      h_electron_d0_wgtCSV[iCat][iSys] = new TH1D("h_electron_d0_wgtCSV" + suffix,";electron d0", 100, -0.05, 0.05 );
      h_electron_dZ_wgtCSV[iCat][iSys] = new TH1D("h_electron_dZ_wgtCSV" + suffix,";electron dZ", 100, -0.05, 0.05 );

      h_muon_pt_wgtCSV[iCat][iSys] = new TH1D("h_muon_pt_wgtCSV" + suffix,";muon p_{T}", NlepPtBins, 0, lepPtMax );
      h_muon_phi_wgtCSV[iCat][iSys] = new TH1D("h_muon_phi_wgtCSV" + suffix,";muon #phi", 34, -3.4, 3.4 );
      h_muon_eta_wgtCSV[iCat][iSys] = new TH1D("h_muon_eta_wgtCSV" + suffix,";muon #eta", 25, -2.5, 2.5 );
      h_muon_relIso_wgtCSV[iCat][iSys] = new TH1D("h_muon_relIso_wgtCSV" + suffix,";muon relIso", 210, -0.01, 0.20 );
      h_muon_d0_wgtCSV[iCat][iSys] = new TH1D("h_muon_d0_wgtCSV" + suffix,";muon d0", 100, -0.05, 0.05 );
      h_muon_dZ_wgtCSV[iCat][iSys] = new TH1D("h_muon_dZ_wgtCSV" + suffix,";muon dZ", 100, -0.05, 0.05 );

      h_jet_pt[iCat][iSys] = new TH1D("h_jet_pt" + suffix,";jet p_{T}", NjetptBins, 0, jetptmax );
      h_jet_eta[iCat][iSys] = new TH1D("h_jet_eta" + suffix,";jet #eta", 25, -2.5, 2.5 );
      h_jet_phi[iCat][iSys] = new TH1D("h_jet_phi" + suffix,";jet #phi", 34, -3.4, 3.4 );
      h_jet_csv[iCat][iSys] = new TH1D("h_jet_csv" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_cmva[iCat][iSys] = new TH1D("h_jet_cmva" + suffix,";jet CMVA", NcsvBins, -1.01, 1.01 );
      h_jet_puMVA[iCat][iSys] = new TH1D("h_jet_puMVA" + suffix,";jet pileupJetId fullDiscriminant", 220, -1.1, 1.1 );

      h_jet_pt_wgtCSV[iCat][iSys] = new TH1D("h_jet_pt_wgtCSV" + suffix,";jet p_{T}", NjetptBins, 0, jetptmax );
      h_jet_eta_wgtCSV[iCat][iSys] = new TH1D("h_jet_eta_wgtCSV" + suffix,";jet #eta", 25, -2.5, 2.5 );
      h_jet_phi_wgtCSV[iCat][iSys] = new TH1D("h_jet_phi_wgtCSV" + suffix,";jet #phi", 34, -3.4, 3.4 );
      h_jet_csv_wgtCSV[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_cmva_wgtCSV[iCat][iSys] = new TH1D("h_jet_cmva_wgtCSV" + suffix,";jet CMVA", NcsvBins, -1.01, 1.01 );
      h_jet_puMVA_wgtCSV[iCat][iSys] = new TH1D("h_jet_puMVA_wgtCSV" + suffix,";jet pileupJetId fullDiscriminant", 220, -1.1, 1.1 );

      h_jet_pt_wgtCSV2[iCat][iSys] = new TH1D("h_jet_pt_wgtCSV2" + suffix,";jet p_{T}", NjetptBins, 0, jetptmax );
      h_jet_eta_wgtCSV2[iCat][iSys] = new TH1D("h_jet_eta_wgtCSV2" + suffix,";jet #eta", 25, -2.5, 2.5 );
      h_jet_phi_wgtCSV2[iCat][iSys] = new TH1D("h_jet_phi_wgtCSV2" + suffix,";jet #phi", 34, -3.4, 3.4 );
      h_jet_csv_wgtCSV2[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV2" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_puMVA_wgtCSV2[iCat][iSys] = new TH1D("h_jet_puMVA_wgtCSV2" + suffix,";jet pileupJetId fullDiscriminant", 220, -1.1, 1.1 );

      h_jet_pt_wgtCSV3[iCat][iSys] = new TH1D("h_jet_pt_wgtCSV3" + suffix,";jet p_{T}", NjetptBins, 0, jetptmax );
      h_jet_eta_wgtCSV3[iCat][iSys] = new TH1D("h_jet_eta_wgtCSV3" + suffix,";jet #eta", 25, -2.5, 2.5 );
      h_jet_phi_wgtCSV3[iCat][iSys] = new TH1D("h_jet_phi_wgtCSV3" + suffix,";jet #phi", 34, -3.4, 3.4 );
      h_jet_csv_wgtCSV3[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV3" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_puMVA_wgtCSV3[iCat][iSys] = new TH1D("h_jet_puMVA_wgtCSV3" + suffix,";jet pileupJetId fullDiscriminant", 220, -1.1, 1.1 );

      h_jet_pt_wgtCSV4[iCat][iSys] = new TH1D("h_jet_pt_wgtCSV4" + suffix,";jet p_{T}", NjetptBins, 0, jetptmax );
      h_jet_eta_wgtCSV4[iCat][iSys] = new TH1D("h_jet_eta_wgtCSV4" + suffix,";jet #eta", 25, -2.5, 2.5 );
      h_jet_phi_wgtCSV4[iCat][iSys] = new TH1D("h_jet_phi_wgtCSV4" + suffix,";jet #phi", 34, -3.4, 3.4 );
      h_jet_csv_wgtCSV4[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV4" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_puMVA_wgtCSV4[iCat][iSys] = new TH1D("h_jet_puMVA_wgtCSV4" + suffix,";jet pileupJetId fullDiscriminant", 220, -1.1, 1.1 );

      h_jet_pt_wgtCSV5[iCat][iSys] = new TH1D("h_jet_pt_wgtCSV5" + suffix,";jet p_{T}", NjetptBins, 0, jetptmax );
      h_jet_eta_wgtCSV5[iCat][iSys] = new TH1D("h_jet_eta_wgtCSV5" + suffix,";jet #eta", 25, -2.5, 2.5 );
      h_jet_phi_wgtCSV5[iCat][iSys] = new TH1D("h_jet_phi_wgtCSV5" + suffix,";jet #phi", 34, -3.4, 3.4 );
      h_jet_csv_wgtCSV5[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV5" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_puMVA_wgtCSV5[iCat][iSys] = new TH1D("h_jet_puMVA_wgtCSV5" + suffix,";jet pileupJetId fullDiscriminant", 220, -1.1, 1.1 );


      h_jet_csv_bFlav[iCat][iSys] = new TH1D("h_jet_csv_bFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_csv_cFlav[iCat][iSys] = new TH1D("h_jet_csv_cFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_csv_lFlav[iCat][iSys] = new TH1D("h_jet_csv_lFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );

      h_jet_csv_wgtCSV_bFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV_bFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_csv_wgtCSV_cFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV_cFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_csv_wgtCSV_lFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV_lFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );

      h_jet_csv_wgtCSV2_bFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV2_bFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_csv_wgtCSV2_cFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV2_cFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_csv_wgtCSV2_lFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV2_lFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );

      h_jet_csv_wgtCSV3_bFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV3_bFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_csv_wgtCSV3_cFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV3_cFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_csv_wgtCSV3_lFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV3_lFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );

      h_jet_csv_wgtCSV4_bFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV4_bFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_csv_wgtCSV4_cFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV4_cFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_csv_wgtCSV4_lFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV4_lFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );

      h_jet_csv_wgtCSV5_bFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV5_bFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_csv_wgtCSV5_cFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV5_cFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_jet_csv_wgtCSV5_lFlav[iCat][iSys] = new TH1D("h_jet_csv_wgtCSV5_lFlav" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );


      
      h_cmva_jet_pt[iCat][iSys] = new TH1D("h_cmva_jet_pt" + suffix,";jet p_{T}", NjetptBins, 0, jetptmax );
      h_cmva_jet_eta[iCat][iSys] = new TH1D("h_cmva_jet_eta" + suffix,";jet #eta", 25, -2.5, 2.5 );
      h_cmva_jet_phi[iCat][iSys] = new TH1D("h_cmva_jet_phi" + suffix,";jet #phi", 34, -3.4, 3.4 );
      h_cmva_jet_csv[iCat][iSys] = new TH1D("h_cmva_jet_csv" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_cmva_jet_cmva[iCat][iSys] = new TH1D("h_cmva_jet_cmva" + suffix,";jet CMVA", NcsvBins, -1.01, 1.01 );
      h_cmva_jet_puMVA[iCat][iSys] = new TH1D("h_cmva_jet_puMVA" + suffix,";jet pileupJetId fullDiscriminant", 220, -1.1, 1.1 );

      h_cmva_jet_pt_wgtCMVA[iCat][iSys] = new TH1D("h_cmva_jet_pt_wgtCMVA" + suffix,";jet p_{T}", NjetptBins, 0, jetptmax );
      h_cmva_jet_eta_wgtCMVA[iCat][iSys] = new TH1D("h_cmva_jet_eta_wgtCMVA" + suffix,";jet #eta", 25, -2.5, 2.5 );
      h_cmva_jet_phi_wgtCMVA[iCat][iSys] = new TH1D("h_cmva_jet_phi_wgtCMVA" + suffix,";jet #phi", 34, -3.4, 3.4 );
      h_cmva_jet_csv_wgtCMVA[iCat][iSys] = new TH1D("h_cmva_jet_csv_wgtCMVA" + suffix,";jet CSV", NcsvBins, -0.06, 1.01 );
      h_cmva_jet_cmva_wgtCMVA[iCat][iSys] = new TH1D("h_cmva_jet_cmva_wgtCMVA" + suffix,";jet CMVA", NcsvBins, -1.01, 1.01 );
      h_cmva_jet_puMVA_wgtCMVA[iCat][iSys] = new TH1D("h_cmva_jet_puMVA_wgtCMVA" + suffix,";jet pileupJetId fullDiscriminant", 220, -1.1, 1.1 );


      
      h_jet_cmva_bFlav[iCat][iSys] = new TH1D("h_jet_cmva_bFlav" + suffix,";jet CMVA", NcsvBins, -1.01, 1.01 );
      h_jet_cmva_cFlav[iCat][iSys] = new TH1D("h_jet_cmva_cFlav" + suffix,";jet CMVA", NcsvBins, -1.01, 1.01 );
      h_jet_cmva_lFlav[iCat][iSys] = new TH1D("h_jet_cmva_lFlav" + suffix,";jet CMVA", NcsvBins, -1.01, 1.01 );

      h_jet_cmva_wgtCMVA_bFlav[iCat][iSys] = new TH1D("h_jet_cmva_wgtCMVA_bFlav" + suffix,";jet CMVA", NcsvBins, -1.01, 1.01 );
      h_jet_cmva_wgtCMVA_cFlav[iCat][iSys] = new TH1D("h_jet_cmva_wgtCMVA_cFlav" + suffix,";jet CMVA", NcsvBins, -1.01, 1.01 );
      h_jet_cmva_wgtCMVA_lFlav[iCat][iSys] = new TH1D("h_jet_cmva_wgtCMVA_lFlav" + suffix,";jet CMVA", NcsvBins, -1.01, 1.01 );


      h_pfMET_pt[iCat][iSys]  = new TH1D("h_pfMET_pt" + suffix,";pfMET p_{T}", NmetBins, 0, metmax );
      h_pfMET_phi[iCat][iSys] = new TH1D("h_pfMET_phi" + suffix,";pfMET #phi", 34, -3.4, 3.4 );
      h_pfMETNoHF_pt[iCat][iSys]  = new TH1D("h_pfMETNoHF_pt" + suffix,";pfMETNoHF p_{T}", NmetBins, 0, metmax );
      h_pfMETNoHF_phi[iCat][iSys] = new TH1D("h_pfMETNoHF_phi" + suffix,";pfMETNoHF #phi", 34, -3.4, 3.4 );
      h_puppiMET_pt[iCat][iSys]  = new TH1D("h_puppiMET_pt" + suffix,";puppiMET p_{T}", NmetBins, 0, metmax );
      h_puppiMET_phi[iCat][iSys] = new TH1D("h_puppiMET_phi" + suffix,";puppiMET #phi", 34, -3.4, 3.4 );
      h_mht_pt[iCat][iSys]  = new TH1D("h_mht_pt" + suffix,";MHT p_{T}", NmetBins, 0, metmax );
      h_mht_phi[iCat][iSys] = new TH1D("h_mht_phi" + suffix,";MHT #phi", 34, -3.4, 3.4 );

      h_pfMET_pt_wgtCSV[iCat][iSys]  = new TH1D("h_pfMET_pt_wgtCSV" + suffix,";pfMET p_{T}", NmetBins, 0, metmax );
      h_pfMET_phi_wgtCSV[iCat][iSys] = new TH1D("h_pfMET_phi_wgtCSV" + suffix,";pfMET #phi", 34, -3.4, 3.4 );

      h_HT[iCat][iSys] = new TH1D("h_HT" + suffix,";H_{T} (jets)", NhtBins, 0, htmax );
    }
  }


  int NumHists = 6;
  int NumSysC = 5;
  TH1D* h_c_jet_csv[NumHists];
  TH1D* h_c_jet_csv_CharmCSVSF[NumSysC][NumHists];

  TH1D* h_c_jet_cmva[NumHists];
  TH1D* h_c_jet_cmva_CharmCMVASF[NumSysC][NumHists];

  int NfullcsvBins = 1000;
  for( int iPt=0; iPt<NumHists; iPt++ ){
    h_c_jet_csv[iPt] = new TH1D(Form("h_c_jet_csv_Pt%d_Eta%d", iPt, 0),";jet CSV", NfullcsvBins, -0.1, 1.1 );
    h_c_jet_cmva[iPt] = new TH1D(Form("h_c_jet_cmva_Pt%d_Eta%d", iPt, 0),";jet CMVA", NfullcsvBins, -1.1, 1.1 );
    for( int iSysC=0; iSysC<NumSysC; iSysC++ ){
      h_c_jet_csv_CharmCSVSF[iSysC][iPt] = new TH1D(Form("h_c_jet_csv_CharmCSVSF_Pt%d_Eta%d_SysC%d", iPt, 0, iSysC),";jet CSV", NfullcsvBins, -0.1, 1.1 );
      h_c_jet_cmva_CharmCMVASF[iSysC][iPt] = new TH1D(Form("h_c_jet_cmva_CharmCMVASF_Pt%d_Eta%d_SysC%d", iPt, 0, iSysC),";jet CMVA", NfullcsvBins, -1.1, 1.1 );
    }
  }


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

    if( verbose_ ) std::cout << " ===> test 0.1 " << std::endl;

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

      // if( keepTTbarEvent==false && ttCat_>=0 ){
      // 	std::cout << " ERROR!! keepTTbarEvent == false! additionalJetEventId = " << additionalJetEventId << " and ttCat_ = " << ttCat_ << std::endl;
      // }

      if( !keepTTbarEvent ) continue;
      //std::cout << " I have made it past the gate! " << std::endl;
    }

    h_additionalJetEventId->Fill(additionalJetEventId);

    int numTruePVs = eve->numTruePVs_;

    //double wgt_pu = ( insample < 0 ) ? 1. : reweightPU(numPVs);
    double wgt_pu = ( insample < 0 ) ? 1. : reweightPU(numTruePVs,0);
    double wgt_pu_up = ( insample < 0 ) ? 1. : reweightPU(numTruePVs,1);
    double wgt_pu_down = ( insample < 0 ) ? 1. : reweightPU(numTruePVs,-1);
    if( verbose_ ) std::cout << " ===> test 0.2 " << std::endl;

    // wgt_pu = 1.;
    // wgt_pu_up = 1.;
    // wgt_pu_down = 1.;
    
    double wgt_gen = ( insample > 0 ) ? eve->wgt_generator_ : 1;
    wgt_gen = ( wgt_gen > 0 ) ? 1. : -1.;
    //if( ievt<3 ){
    std::vector<double> lhe_event_weights = eve->LHEEvent_weights_;
    //std::cout << "\t originalXWGTUP = " << eve->originalXWGTUP_ << std::endl;
    for( unsigned int iwgt = 0; iwgt < lhe_event_weights.size(); iwgt++ ){
      h_numEvents_lheWeight->Fill(iwgt, lhe_event_weights[iwgt]/eve->originalXWGTUP_);
    	//printf(" weight %5d = %f,\t ratio = %4.5f\n", iwgt, lhe_event_weights[iwgt], lhe_event_weights[iwgt]/eve->originalXWGTUP_);
    }
      //printf(" Event %3lld: additionalJetEventId = %3d \n", ievt, additionalJetEventId);
    //}



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
	if( pt > 25 && abs(eta)<2.1 && isTight && relIsoR04 < 0.15 ) ind_mu.push_back(iLep);
	if( pt > 15 && abs(eta)<2.4 && isLoose && relIsoR04 < 0.25 ) ind_mu_loose.push_back(iLep);
      }

      if( !isMuon ){
	if( pt > 30 && abs(eta)<2.1 && isTrigMVAM && !isCrack && relIso < 0.15 ) ind_ele.push_back(iLep);
	if( pt > 15 && abs(eta)<2.4 && isTrigMVAM && !isCrack && relIso < 0.15 ) ind_ele_loose.push_back(iLep);
      }
    }

        if( verbose_ ) std::cout << " ===> test 0.3 " << std::endl;

    std::vector<double> rejetPts;
    std::vector<double> rejetEtas;
    std::vector<double> rejetCSVs;
    std::vector<int>    rejetFlavors;

    int njet = 0;
    int nbtag = 0;

    for( int iJet=0; iJet<int(eve->jet_pt_.size()); iJet++ ){
       if( verbose_ ) std::cout << " ===> test 0.3.0.1 " << std::endl;
      double pt  = eve->jet_pt_[iJet];
      double eta = eve->jet_eta_[iJet];
      int flavor = eve->jet_hadronFlavour_[iJet];
      double csv =  eve->jet_csv_[iJet];
      double cmva = eve->jet_cmva_[iJet];
       if( verbose_ ) std::cout << " ===> test 0.3.0.2 " << std::endl;

      double jetAbsEta = fabs(eta);
      
      if( csv < 0.0 ) csv = -0.05;
      if( csv > 1.0 ) csv = 1.0;

      if( !(pt>20. && fabs(eta)<2.4) ) continue;
             if( verbose_ ) std::cout << " ===> test 0.3.0.3 " << std::endl;

      int iPt = -1, iEta = -1;
      if (pt >=19.99 && pt<30) iPt = 0;
      else if (pt >=30 && pt<40) iPt = 1;
      else if (pt >=40 && pt<60) iPt = 2;
      else if (pt >=60 && pt<100) iPt = 3;
      else if (pt >=100)          iPt = 4;
      //else if (pt >=100 && pt<160) iPt = 4;
      //else if (pt >=160 && pt<10000) iPt = 5;

      if( PtBinsHF_ > 5 && pt >=160 )  iPt = 5;
       if( verbose_ ) std::cout << " ===> test 0.3.0.4 " << std::endl;

      
      if (jetAbsEta >=0 &&  jetAbsEta<0.8 ) iEta = 0;
      else if ( jetAbsEta>=0.8 && jetAbsEta<1.6 )  iEta = 1;
      else if ( jetAbsEta>=1.6 && jetAbsEta<2.41 ) iEta = 2;
        if( verbose_ ) std::cout << " ===> test 0.3.1 " << std::endl;


      if( abs(flavor)==4 && fabs(eta)<2.4 && pt>20. ){
	h_c_jet_csv[iPt]->Fill(csv);
	h_c_jet_cmva[iPt]->Fill(cmva);
	for( int iSysC=0; iSysC<NumSysC; iSysC++ ){
	  int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
	  double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);

	  int useCMVABin = c_cmva_wgt_hf[iSysC][iPt]->FindBin(cmva);
	  double iCMVAWgtC = c_cmva_wgt_hf[iSysC][iPt]->GetBinContent(useCMVABin);

	  h_c_jet_csv_CharmCSVSF[iSysC][iPt]->Fill(csv,iCSVWgtC);
	  h_c_jet_cmva_CharmCMVASF[iSysC][iPt]->Fill(cmva,iCMVAWgtC);
	}
      }
       if( verbose_ ) std::cout << " ===> test 0.3.2 " << std::endl;

      if( !(pt>30. && fabs(eta)<2.4) ) continue;

      h_jet_csv_a->Fill(csv);
      if( abs(flavor)==5 )      h_jet_csv_b->Fill(csv);
      else if( abs(flavor)==4 ) h_jet_csv_c->Fill(csv);
      else                      h_jet_csv_l->Fill(csv);

      h_jet_cmva_a->Fill(csv);
      if( abs(flavor)==5 )      h_jet_cmva_b->Fill(cmva);
      else if( abs(flavor)==4 ) h_jet_cmva_c->Fill(cmva);
      else                      h_jet_cmva_l->Fill(cmva);

      njet++;
      if( csv>0.80 ) nbtag++;

      rejetPts.push_back(pt);
      rejetEtas.push_back(eta);
      rejetCSVs.push_back(csv);
      rejetFlavors.push_back(flavor);


      pt = std::min(pt, jetptmax-0.001);
      eta = fabs(eta);
      
      h_a_jet_pt_eta_all->Fill(pt, eta);
      if( abs(flavor)==5 )      h_b_jet_pt_eta_all->Fill(pt, eta);
      else if( abs(flavor)==4 ) h_c_jet_pt_eta_all->Fill(pt, eta);
      else                      h_l_jet_pt_eta_all->Fill(pt, eta);

      if( csv>0.80 ){
	h_a_jet_pt_eta_csvM->Fill(pt, eta);
	if( abs(flavor)==5 )      h_b_jet_pt_eta_csvM->Fill(pt, eta);
	else if( abs(flavor)==4 ) h_c_jet_pt_eta_csvM->Fill(pt, eta);
	else                      h_l_jet_pt_eta_csvM->Fill(pt, eta);
      }	
    }
    if( verbose_ ) std::cout << " ===> test 0.4 " << std::endl;




    double wgt_csv_noSys_partonFlavourNoPU = 1;
    double wgt_csv_hf_noSys_partonFlavourNoPU = 1;
    double wgt_csv_lf_noSys_partonFlavourNoPU = 1;
    double wgt_csv_cf_noSys_partonFlavourNoPU = 1;

    for( int iJet=0; iJet<int(eve->jet_pt_.size()); iJet++ ){
      double pt  = eve->jet_pt_[iJet];
      double eta = eve->jet_eta_[iJet];
      int flavor = eve->jet_hadronFlavour_[iJet];
      double csv =  eve->jet_csv_[iJet];

      int partonFlavor = eve->jet_partonFlavour_[iJet];
      
      double jetAbsEta = fabs(eta);
      
      if( csv < 0.0 ) csv = -0.05;
      if( csv > 1.0 ) csv = 1.0;

      if( !(pt>20. && fabs(eta)<2.4) ) continue;
      
      int iPt = -1, iEta = -1;
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


      int iSysHF = 0;
      int iSysC  = 0;
      int iSysLF = 0;

      if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, pt = " << pt << ", jetAbsEta = " << jetAbsEta << std::endl;
      if( verbose_ ) std::cout << " ===> 1 test 1.3 jet " << iJet << std::endl;

      if( abs(partonFlavor)==5 ){
	if( iEta>=1 ) iEta=0;
	h_perjet_csv_bPar[iPt][iEta]->Fill(csv,wgt_gen);
      }
      else if( abs(partonFlavor)==4 ){
	if( iEta>=1 ) iEta=0;
	h_perjet_csv_cPar[iPt][iEta]->Fill(csv,wgt_gen);
      }
      else if( (abs(partonFlavor)>0 && abs(partonFlavor)<4) || abs(partonFlavor)==21 ){
	if (iPt >=3) iPt=3;
	h_perjet_csv_udsgPar[iPt][iEta]->Fill(csv,wgt_gen);
	h_perjet_csv_lfPar[iPt][iEta]->Fill(csv,wgt_gen);
      }
      else{
	if (iPt >=3) iPt=3;
	h_perjet_csv_oPar[iPt][iEta]->Fill(csv,wgt_gen);
	h_perjet_csv_lfPar[iPt][iEta]->Fill(csv,wgt_gen);
      }
      
      if (abs(flavor) == 5 ){
	if( iEta>=1 ) iEta=0;
	int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
	double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
	
	if( iCSVWgtHF!=0 && partonFlavor!=0 ) wgt_csv_hf_noSys_partonFlavourNoPU *= iCSVWgtHF;

	h_perjet_wgt_csv_hf[iPt][iEta]->Fill(iCSVWgtHF,wgt_gen);
	h_perjet_wgt_csv[iPt][iEta]->Fill(iCSVWgtHF,wgt_gen);
	if( njet>=4 ){
	  h_perjet_wgt_csv_4j_hf[iPt][iEta]->Fill(iCSVWgtHF,wgt_gen);
	  h_perjet_wgt_csv_4j[iPt][iEta]->Fill(iCSVWgtHF,wgt_gen);
	}

	h_perjet_csv_hf[iPt][iEta]->Fill(csv,wgt_gen);
	h_perjet_csv[iPt][iEta]->Fill(csv,wgt_gen);
	if( njet>=4 ){
	  h_perjet_csv_4j_hf[iPt][iEta]->Fill(csv,wgt_gen);
	  h_perjet_csv_4j[iPt][iEta]->Fill(csv,wgt_gen);
	}

	h_perjet_csv_PU_hf[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	h_perjet_csv_PU[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	if( njet>=4 ){
	  h_perjet_csv_PU_4j_hf[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	  h_perjet_csv_PU_4j[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	}
      }
      else if( abs(flavor) == 4 ){
	if( iEta>=1 ) iEta=0;
	int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
	double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);

	if( iCSVWgtC!=0 && partonFlavor!=0 ) wgt_csv_cf_noSys_partonFlavourNoPU *= iCSVWgtC;

	h_perjet_wgt_csv_cf[iPt][iEta]->Fill(iCSVWgtC,wgt_gen);
	h_perjet_wgt_csv[iPt][iEta]->Fill(iCSVWgtC,wgt_gen);
	if( njet>=4 ){
	  h_perjet_wgt_csv_4j_cf[iPt][iEta]->Fill(iCSVWgtC,wgt_gen);
	  h_perjet_wgt_csv_4j[iPt][iEta]->Fill(iCSVWgtC,wgt_gen);
	}
	
	h_perjet_csv_cf[iPt][iEta]->Fill(csv,wgt_gen);
	h_perjet_csv[iPt][iEta]->Fill(csv,wgt_gen);
	if( njet>=4 ){
	  h_perjet_csv_4j_cf[iPt][iEta]->Fill(csv,wgt_gen);
	  h_perjet_csv_4j[iPt][iEta]->Fill(csv,wgt_gen);
	}
	
	h_perjet_csv_PU_cf[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	h_perjet_csv_PU[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	if( njet>=4 ){
	  h_perjet_csv_PU_4j_cf[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	  h_perjet_csv_PU_4j[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	}
      }
      else {
	if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
	int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
	double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);

	if( iCSVWgtLF!=0 && partonFlavor!=0 ) wgt_csv_lf_noSys_partonFlavourNoPU *= iCSVWgtLF;

	h_perjet_wgt_csv_lf[iPt][iEta]->Fill(iCSVWgtLF,wgt_gen);
	h_perjet_wgt_csv[iPt][iEta]->Fill(iCSVWgtLF,wgt_gen);
	if( njet>=4 ){
	  h_perjet_wgt_csv_4j_lf[iPt][iEta]->Fill(iCSVWgtLF,wgt_gen);
	  h_perjet_wgt_csv_4j[iPt][iEta]->Fill(iCSVWgtLF,wgt_gen);
	}


	if( partonFlavor!=0 ){
	  h_perjet_wgt_csv_lf_partonFlavourNoPU[iPt][iEta]->Fill(iCSVWgtLF,wgt_gen);
	}
	
	h_perjet_csv_lf[iPt][iEta]->Fill(csv,wgt_gen);
	h_perjet_csv[iPt][iEta]->Fill(csv,wgt_gen);
	if( njet>=4 ){
	  h_perjet_csv_4j_lf[iPt][iEta]->Fill(csv,wgt_gen);
	  h_perjet_csv_4j[iPt][iEta]->Fill(csv,wgt_gen);
	}

	h_perjet_csv_PU_lf[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	h_perjet_csv_PU[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	if( njet>=4 ){
	  h_perjet_csv_PU_4j_lf[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	  h_perjet_csv_PU_4j[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	}
      }
    }

    wgt_csv_noSys_partonFlavourNoPU = wgt_csv_hf_noSys_partonFlavourNoPU * wgt_csv_cf_noSys_partonFlavourNoPU * wgt_csv_lf_noSys_partonFlavourNoPU;
        if( verbose_ ) std::cout << " ===> test 0.5 " << std::endl;


    // JESup
    std::vector<double> rejetPts_JESup;
    std::vector<double> rejetEtas_JESup;
    std::vector<double> rejetCSVs_JESup;
    std::vector<int>    rejetFlavors_JESup;

    int njet_JESup = 0;
    int nbtag_JESup = 0;

    for( int iJet=0; iJet<int(eve->jet_JESup_pt_.size()); iJet++ ){
      double pt  = eve->jet_JESup_pt_[iJet];
      double eta = eve->jet_JESup_eta_[iJet];
      int flavor = eve->jet_JESup_hadronFlavour_[iJet];
      double csv =  eve->jet_JESup_csv_[iJet];

      if( csv < 0.0 ) csv = -0.05;
      if( csv > 1.0 ) csv = 1.0;

      if( !(pt>30. && fabs(eta)<2.4) ) continue;

      njet_JESup++;
      if( csv>0.80 ) nbtag_JESup++;

      rejetPts_JESup.push_back(pt);
      rejetEtas_JESup.push_back(eta);
      rejetCSVs_JESup.push_back(csv);
      rejetFlavors_JESup.push_back(flavor);
    }



    // JESdown
    std::vector<double> rejetPts_JESdown;
    std::vector<double> rejetEtas_JESdown;
    std::vector<double> rejetCSVs_JESdown;
    std::vector<int>    rejetFlavors_JESdown;

    int njet_JESdown = 0;
    int nbtag_JESdown = 0;

    for( int iJet=0; iJet<int(eve->jet_JESdown_pt_.size()); iJet++ ){
      double pt  = eve->jet_JESdown_pt_[iJet];
      double eta = eve->jet_JESdown_eta_[iJet];
      int flavor = eve->jet_JESdown_hadronFlavour_[iJet];
      double csv =  eve->jet_JESdown_csv_[iJet];
      if( csv < 0.0 ) csv = -0.05;
      if( csv > 1.0 ) csv = 1.0;

      if( !(pt>30. && fabs(eta)<2.4) ) continue;

      njet_JESdown++;
      if( csv>0.80 ) nbtag_JESdown++;

      rejetPts_JESdown.push_back(pt);
      rejetEtas_JESdown.push_back(eta);
      rejetCSVs_JESdown.push_back(csv);
      rejetFlavors_JESdown.push_back(flavor);
    }

    if( verbose_ ) std::cout << " ===> test 1 " << std::endl;

    int numPVs = eve->numPVs_;

    double wgt_csv_noSys = 1;
    double wgt_csv_hf_noSys = 1;
    double wgt_csv_lf_noSys = 1;
    double wgt_csv_cf_noSys = 1;
    if( verbose_ ) std::cout << " ===> test 0.6 " << std::endl;

    for( int useSys=0; useSys<int(csv_systematics.size()); useSys++ ){
      int mySys = csv_systematics[useSys];

      double wgt_csv_hf, wgt_csv_lf, wgt_csv_cf;

      double temp_wgt_csv = 1.;//( insample<0 ) ? 1 : get_csv_wgt(rejetPts, rejetEtas, rejetCSVs, rejetFlavors, mySys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);

      if( !(insample < 0) ){
      	if( mySys==7 )      temp_wgt_csv = get_csv_wgt(rejetPts_JESup, rejetEtas_JESup, rejetCSVs_JESup, rejetFlavors_JESup, mySys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
      	else if( mySys==8 ) temp_wgt_csv = get_csv_wgt(rejetPts_JESdown, rejetEtas_JESdown, rejetCSVs_JESdown, rejetFlavors_JESdown, mySys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
      	else                temp_wgt_csv = get_csv_wgt(rejetPts, rejetEtas, rejetCSVs, rejetFlavors, mySys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
      }

      if( useSys==0 ){
	wgt_csv_noSys = temp_wgt_csv;
	wgt_csv_hf_noSys = wgt_csv_hf;
	wgt_csv_lf_noSys = wgt_csv_lf;
	wgt_csv_cf_noSys = wgt_csv_cf;


	h_numJet_nolepreq->Fill(njet,wgt_gen*wgt_pu);
	h_numJet_nolepreq_wgtCSV->Fill(njet,wgt_gen*wgt_pu*temp_wgt_csv);
      }

      
      h_numEvents_Sys->Fill(0.5+mySys,temp_wgt_csv);
      h_numEvents_Sys_PU->Fill(0.5+mySys,temp_wgt_csv*wgt_pu);

      h_wgt_csv_Sys[useSys]->Fill(temp_wgt_csv,wgt_pu);
      h_wgt_csv_hf_Sys[useSys]->Fill(wgt_csv_hf,wgt_pu);
      h_wgt_csv_lf_Sys[useSys]->Fill(wgt_csv_lf,wgt_pu);
      h_wgt_csv_cf_Sys[useSys]->Fill(wgt_csv_cf,wgt_pu);

      h_wgt_csv_noPU_Sys[useSys]->Fill(temp_wgt_csv,1);
      h_wgt_csv_hf_noPU_Sys[useSys]->Fill(wgt_csv_hf,1);
      h_wgt_csv_lf_noPU_Sys[useSys]->Fill(wgt_csv_lf,1);
      h_wgt_csv_cf_noPU_Sys[useSys]->Fill(wgt_csv_cf,1);

      int temp_njet = njet;
      int temp_nbtag = nbtag;

      if( insample < 0 ){
	if( mySys==7 ){
	  temp_njet = njet_JESup;
	  temp_nbtag = nbtag_JESup;
	}
	else if( mySys==8 ){
	  temp_njet = njet_JESdown;
	  temp_nbtag = nbtag_JESdown;
	}
	else {
	  temp_njet = njet_JESdown;
	  temp_nbtag = nbtag_JESdown;
	}
      }

      if( temp_njet>=4 && temp_nbtag>=2 ){
	h_numEvents_4j2t_Sys->Fill(0.5+mySys,temp_wgt_csv*wgt_pu);
	h_wgt_csv_4j2t_Sys[useSys]->Fill(temp_wgt_csv,wgt_pu);
	h_wgt_csv_hf_4j2t_Sys[useSys]->Fill(wgt_csv_hf,wgt_pu);
	h_wgt_csv_lf_4j2t_Sys[useSys]->Fill(wgt_csv_lf,wgt_pu);
	h_wgt_csv_cf_4j2t_Sys[useSys]->Fill(wgt_csv_cf,wgt_pu);
      }
      
      if( temp_njet>=4 ){
	h_numEvents_4j_Sys->Fill(0.5+mySys,temp_wgt_csv*wgt_pu);
	h_wgt_csv_4j_Sys[useSys]->Fill(temp_wgt_csv,wgt_pu);
	h_wgt_csv_hf_4j_Sys[useSys]->Fill(wgt_csv_hf,wgt_pu);
	h_wgt_csv_lf_4j_Sys[useSys]->Fill(wgt_csv_lf,wgt_pu);
	h_wgt_csv_cf_4j_Sys[useSys]->Fill(wgt_csv_cf,wgt_pu);

	h_wgt_csv_4j_noPU_Sys[useSys]->Fill(temp_wgt_csv,1);
	h_wgt_csv_hf_4j_noPU_Sys[useSys]->Fill(wgt_csv_hf,1);
	h_wgt_csv_lf_4j_noPU_Sys[useSys]->Fill(wgt_csv_lf,1);
	h_wgt_csv_cf_4j_noPU_Sys[useSys]->Fill(wgt_csv_cf,1);
      }



      //// NEW

      vdouble use_jet_pt = eve->jet_pt_;
      vdouble use_jet_eta = eve->jet_eta_;
      vdouble use_jet_csv = eve->jet_csv_;
      vdouble use_jet_cmva = eve->jet_cmva_;
      vint use_jet_hadronFlavour = eve->jet_hadronFlavour_;

      if( mySys==7 ){
	use_jet_pt = eve->jet_JESup_pt_;
	use_jet_eta = eve->jet_JESup_eta_;
	use_jet_csv = eve->jet_JESup_csv_;
	use_jet_cmva = eve->jet_JESup_cmva_;
	use_jet_hadronFlavour = eve->jet_JESup_hadronFlavour_;
      }
      else if( mySys==8 ){
	use_jet_pt = eve->jet_JESdown_pt_;
	use_jet_eta = eve->jet_JESdown_eta_;
	use_jet_csv = eve->jet_JESdown_csv_;
	use_jet_cmva = eve->jet_JESdown_cmva_;
	use_jet_hadronFlavour = eve->jet_JESdown_hadronFlavour_;
      }
      else {
	use_jet_pt = eve->jet_pt_;
	use_jet_eta = eve->jet_eta_;
	use_jet_csv = eve->jet_csv_;
	use_jet_cmva = eve->jet_cmva_;
	use_jet_hadronFlavour = eve->jet_hadronFlavour_;
      }


      double wgt_jet30_CSVL = 0.;
      double wgt_jet30_CSVM = 0.;
      double wgt_jet30_CSVT = 0.;

      double nowgt_jet30_CSVL = 0.;
      double nowgt_jet30_CSVM = 0.;
      double nowgt_jet30_CSVT = 0.;

      double wgt_jet20_CSVL = 0.;
      double wgt_jet20_CSVM = 0.;
      double wgt_jet20_CSVT = 0.;

      double nowgt_jet20_CSVL = 0.;
      double nowgt_jet20_CSVM = 0.;
      double nowgt_jet20_CSVT = 0.;


      double wgt_jet30_CMVAL = 0.;
      double wgt_jet30_CMVAM = 0.;
      double wgt_jet30_CMVAT = 0.;

      double nowgt_jet30_CMVAL = 0.;
      double nowgt_jet30_CMVAM = 0.;
      double nowgt_jet30_CMVAT = 0.;

      double wgt_jet20_CMVAL = 0.;
      double wgt_jet20_CMVAM = 0.;
      double wgt_jet20_CMVAT = 0.;

      double nowgt_jet20_CMVAL = 0.;
      double nowgt_jet20_CMVAM = 0.;
      double nowgt_jet20_CMVAT = 0.;


      double wgt_jet30to60_CMVAL = 0.;
      double wgt_jet30to60_CMVAM = 0.;
      double wgt_jet30to60_CMVAT = 0.;

      double nowgt_jet30to60_CMVAL = 0.;
      double nowgt_jet30to60_CMVAM = 0.;
      double nowgt_jet30to60_CMVAT = 0.;


      double wgt_jet60to120_CMVAL = 0.;
      double wgt_jet60to120_CMVAM = 0.;
      double wgt_jet60to120_CMVAT = 0.;

      double nowgt_jet60to120_CMVAL = 0.;
      double nowgt_jet60to120_CMVAM = 0.;
      double nowgt_jet60to120_CMVAT = 0.;


      double wgt_jet120toInf_CMVAL = 0.;
      double wgt_jet120toInf_CMVAM = 0.;
      double wgt_jet120toInf_CMVAT = 0.;

      double nowgt_jet120toInf_CMVAL = 0.;
      double nowgt_jet120toInf_CMVAM = 0.;
      double nowgt_jet120toInf_CMVAT = 0.;

      for( int iJet=0; iJet<int(use_jet_pt.size()); iJet++ ){
	double pt  = use_jet_pt[iJet];
	double eta = use_jet_eta[iJet];
	int flavor = use_jet_hadronFlavour[iJet];
	double csv =  use_jet_csv[iJet];
	double cmva =  use_jet_cmva[iJet];

	double jetAbsEta = fabs(eta);
      
	if( csv < 0.0 ) csv = -0.05;
	if( csv > 1.0 ) csv = 1.0;

	if( !(pt>20. && fabs(eta)<2.4) ) continue;
      
	int iPt = -1, iEta = -1;
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

	int iSysHF = 0;
	switch(mySys){
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
	switch(mySys){
	case 21: iSysC=1; break;
	case 22: iSysC=2; break;
	case 23: iSysC=3; break;
	case 24: iSysC=4; break;
	default : iSysC = 0; break;
	}

	int iSysLF = 0;
	switch(mySys){
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
  
	if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, pt = " << pt << ", jetAbsEta = " << jetAbsEta << std::endl;
	if( verbose_ ) std::cout << " ===> 2 test 1.3 jet " << iJet << std::endl;
      
	if (abs(flavor) == 5 ){
	  if( iEta>=1 ) iEta=0;
	  int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
	  double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);

	  h_perjet_wgt_csv_hf_Sys[iPt][iEta][useSys]->Fill(iCSVWgtHF,wgt_gen);
	  h_perjet_wgt_csv_Sys[iPt][iEta][useSys]->Fill(iCSVWgtHF,wgt_gen);

	  if( verbose_ ) std::cout << " ===> 2 test 1.4 jet " << iJet << std::endl;
	  int useCMVABin = h_cmva_wgt_hf[iSysHF][iPt]->FindBin(cmva);
	  double iCMVAWgtHF = h_cmva_wgt_hf[iSysHF][iPt]->GetBinContent(useCMVABin);
	  if( verbose_ ) std::cout << " ===> 2 test 1.5 jet " << iJet << std::endl;

	  h_perjet_wgt_cmva_hf_Sys[iPt][iEta][useSys]->Fill(iCMVAWgtHF,wgt_gen);
	  h_perjet_wgt_cmva_Sys[iPt][iEta][useSys]->Fill(iCMVAWgtHF,wgt_gen);
	  if( verbose_ ) std::cout << " ===> 2 test 1.6 jet " << iJet << std::endl;

	  if( csv > 0.460 && pt > 30 && iCSVWgtHF>0. ){
	    wgt_jet30_CSVL += iCSVWgtHF;
	    nowgt_jet30_CSVL += 1.0;
	  }
	  if( csv > 0.800 && pt > 30 && iCSVWgtHF>0. ){
	    wgt_jet30_CSVM += iCSVWgtHF;
	    nowgt_jet30_CSVM += 1.0;
	  }
	  if( csv > 0.935 && pt > 30 && iCSVWgtHF>0. ){
	    wgt_jet30_CSVT += iCSVWgtHF;
	    nowgt_jet30_CSVT += 1.0;
	  }

	  
	  if( csv > 0.460 && pt > 20 && iCSVWgtHF>0. ){
	    wgt_jet20_CSVL += iCSVWgtHF;
	    nowgt_jet20_CSVL += 1.0;
	  }
	  if( csv > 0.800 && pt > 20 && iCSVWgtHF>0. ){
	    wgt_jet20_CSVM += iCSVWgtHF;
	    nowgt_jet20_CSVM += 1.0;
	  }
	  if( csv > 0.935 && pt > 20 && iCSVWgtHF>0. ){
	    wgt_jet20_CSVT += iCSVWgtHF;
	    nowgt_jet20_CSVT += 1.0;
	  }



	  if( cmva > -0.715 && pt > 30 && iCMVAWgtHF>0. ){
	    wgt_jet30_CMVAL += iCMVAWgtHF;
	    nowgt_jet30_CMVAL += 1.0;
	  }
	  if( cmva > 0.185 && pt > 30 && iCMVAWgtHF>0. ){
	    wgt_jet30_CMVAM += iCMVAWgtHF;
	    nowgt_jet30_CMVAM += 1.0;
	  }
	  if( cmva > 0.875 && pt > 30 && iCMVAWgtHF>0. ){
	    wgt_jet30_CMVAT += iCMVAWgtHF;
	    nowgt_jet30_CMVAT += 1.0;
	  }

	  
	  if( cmva > -0.715 && pt > 20 && iCMVAWgtHF>0. ){
	    wgt_jet20_CMVAL += iCMVAWgtHF;
	    nowgt_jet20_CMVAL += 1.0;
	  }
	  if( cmva > 0.185 && pt > 20 && iCMVAWgtHF>0. ){
	    wgt_jet20_CMVAM += iCMVAWgtHF;
	    nowgt_jet20_CMVAM += 1.0;
	  }
	  if( cmva > 0.875 && pt > 20 && iCMVAWgtHF>0. ){
	    wgt_jet20_CMVAT += iCMVAWgtHF;
	    nowgt_jet20_CMVAT += 1.0;
	  }


	  

	  if( pt>=30 && pt<60 && iCMVAWgtHF>0. ){
	    if( cmva > -0.715 ){
	      wgt_jet30to60_CMVAL += iCMVAWgtHF;
	      nowgt_jet30to60_CMVAL += 1.0;
	    }
	    if( cmva > 0.185 ){
	      wgt_jet30to60_CMVAM += iCMVAWgtHF;
	      nowgt_jet30to60_CMVAM += 1.0;
	    }
	    if( cmva > 0.875 ){
	      wgt_jet30to60_CMVAT += iCMVAWgtHF;
	      nowgt_jet30to60_CMVAT += 1.0;
	    }
	  }

	  

	  if( pt>=60 && pt<120 && iCMVAWgtHF>0. ){
	    if( cmva > -0.715 ){
	      wgt_jet60to120_CMVAL += iCMVAWgtHF;
	      nowgt_jet60to120_CMVAL += 1.0;
	    }
	    if( cmva > 0.185 ){
	      wgt_jet60to120_CMVAM += iCMVAWgtHF;
	      nowgt_jet60to120_CMVAM += 1.0;
	    }
	    if( cmva > 0.875 ){
	      wgt_jet60to120_CMVAT += iCMVAWgtHF;
	      nowgt_jet60to120_CMVAT += 1.0;
	    }
	  }
 

	  if( pt>=120 && iCMVAWgtHF>0. ){
	    if( cmva > -0.715 ){
	      wgt_jet120toInf_CMVAL += iCMVAWgtHF;
	      nowgt_jet120toInf_CMVAL += 1.0;
	    }
	    if( cmva > 0.185 ){
	      wgt_jet120toInf_CMVAM += iCMVAWgtHF;
	      nowgt_jet120toInf_CMVAM += 1.0;
	    }
	    if( cmva > 0.875 ){
	      wgt_jet120toInf_CMVAT += iCMVAWgtHF;
	      nowgt_jet120toInf_CMVAT += 1.0;
	    }
	  }



	  
	}
	else if( abs(flavor) == 4 ){
	  if( iEta>=1 ) iEta=0;
	  int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
	  double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);

	  h_perjet_wgt_csv_cf_Sys[iPt][iEta][useSys]->Fill(iCSVWgtC,wgt_gen);
	  h_perjet_wgt_csv_Sys[iPt][iEta][useSys]->Fill(iCSVWgtC,wgt_gen);

	  
	  if( verbose_ ) std::cout << " ===> 2 test 1.7 jet " << iJet << std::endl;
	  int useCMVABin = c_cmva_wgt_hf[iSysC][iPt]->FindBin(cmva);
	  double iCMVAWgtC = c_cmva_wgt_hf[iSysC][iPt]->GetBinContent(useCMVABin);
	  if( verbose_ ) std::cout << " ===> 2 test 1.8 jet " << iJet << std::endl;

	  h_perjet_wgt_cmva_cf_Sys[iPt][iEta][useSys]->Fill(iCMVAWgtC,wgt_gen);
	  h_perjet_wgt_cmva_Sys[iPt][iEta][useSys]->Fill(iCMVAWgtC,wgt_gen);
	  if( verbose_ ) std::cout << " ===> 2 test 1.9 jet " << iJet << std::endl;
	}
	else {
	  if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
	  int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
	  double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);

	  h_perjet_wgt_csv_lf_Sys[iPt][iEta][useSys]->Fill(iCSVWgtLF,wgt_gen);
	  h_perjet_wgt_csv_Sys[iPt][iEta][useSys]->Fill(iCSVWgtLF,wgt_gen);

	  
	  if( verbose_ ) std::cout << " ===> 2 test 1.10 jet " << iJet << std::endl;
	  int useCMVABin = h_cmva_wgt_lf[iSysLF][iPt][iEta]->FindBin(cmva);
	  double iCMVAWgtLF = h_cmva_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCMVABin);
	  if( verbose_ ) std::cout << " ===> 2 test 1.11 jet " << iJet << std::endl;

	  h_perjet_wgt_cmva_lf_Sys[iPt][iEta][useSys]->Fill(iCMVAWgtLF,wgt_gen);
	  h_perjet_wgt_cmva_Sys[iPt][iEta][useSys]->Fill(iCMVAWgtLF,wgt_gen);
	  if( verbose_ ) std::cout << " ===> 2 test 1.12 jet " << iJet << std::endl;
	}
      }
    if( verbose_ ) std::cout << " ===> test 0.7 " << std::endl;

      h_numEvents_jet30_wgtCSV_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet30_CSVL);
      h_numEvents_jet30_wgtCSV_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet30_CSVM);
      h_numEvents_jet30_wgtCSV_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet30_CSVT);

      h_numEvents_jet30_nowgtCSV_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet30_CSVL);
      h_numEvents_jet30_nowgtCSV_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet30_CSVM);
      h_numEvents_jet30_nowgtCSV_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet30_CSVT);

      h_numEvents_jet20_wgtCSV_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet20_CSVL);
      h_numEvents_jet20_wgtCSV_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet20_CSVM);
      h_numEvents_jet20_wgtCSV_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet20_CSVT);

      h_numEvents_jet20_nowgtCSV_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet20_CSVL);
      h_numEvents_jet20_nowgtCSV_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet20_CSVM);
      h_numEvents_jet20_nowgtCSV_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet20_CSVT);



      h_numEvents_jet30_wgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet30_CMVAL);
      h_numEvents_jet30_wgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet30_CMVAM);
      h_numEvents_jet30_wgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet30_CMVAT);

      h_numEvents_jet30_nowgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet30_CMVAL);
      h_numEvents_jet30_nowgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet30_CMVAM);
      h_numEvents_jet30_nowgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet30_CMVAT);

      h_numEvents_jet20_wgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet20_CMVAL);
      h_numEvents_jet20_wgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet20_CMVAM);
      h_numEvents_jet20_wgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet20_CMVAT);

      h_numEvents_jet20_nowgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet20_CMVAL);
      h_numEvents_jet20_nowgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet20_CMVAM);
      h_numEvents_jet20_nowgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet20_CMVAT);


      
      h_numEvents_jet30to60_wgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet30to60_CMVAL);
      h_numEvents_jet30to60_wgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet30to60_CMVAM);
      h_numEvents_jet30to60_wgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet30to60_CMVAT);

      h_numEvents_jet30to60_nowgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet30to60_CMVAL);
      h_numEvents_jet30to60_nowgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet30to60_CMVAM);
      h_numEvents_jet30to60_nowgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet30to60_CMVAT);

      
      h_numEvents_jet60to120_wgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet60to120_CMVAL);
      h_numEvents_jet60to120_wgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet60to120_CMVAM);
      h_numEvents_jet60to120_wgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet60to120_CMVAT);

      h_numEvents_jet60to120_nowgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet60to120_CMVAL);
      h_numEvents_jet60to120_nowgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet60to120_CMVAM);
      h_numEvents_jet60to120_nowgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet60to120_CMVAT);

      
      h_numEvents_jet120toInf_wgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*wgt_jet120toInf_CMVAL);
      h_numEvents_jet120toInf_wgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*wgt_jet120toInf_CMVAM);
      h_numEvents_jet120toInf_wgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*wgt_jet120toInf_CMVAT);

      h_numEvents_jet120toInf_nowgtCMVA_WP[useSys]->Fill(0.5,wgt_gen*nowgt_jet120toInf_CMVAL);
      h_numEvents_jet120toInf_nowgtCMVA_WP[useSys]->Fill(1.5,wgt_gen*nowgt_jet120toInf_CMVAM);
      h_numEvents_jet120toInf_nowgtCMVA_WP[useSys]->Fill(2.5,wgt_gen*nowgt_jet120toInf_CMVAT);


      //// END NEW
    }
    if( verbose_ ) std::cout << " ===> test 2 " << std::endl;

    if( insample<0 ){
      wgt_csv_noSys = 1.0;
      wgt_csv_hf_noSys = 1.0;
      wgt_csv_lf_noSys = 1.0;
      wgt_csv_cf_noSys = 1.0;

      wgt_csv_noSys_partonFlavourNoPU = 1.0;
      wgt_csv_hf_noSys_partonFlavourNoPU = 1.0;
      wgt_csv_lf_noSys_partonFlavourNoPU = 1.0;
      wgt_csv_cf_noSys_partonFlavourNoPU = 1.0;
    }
    int run = eve->run_;

    //// --------- various weights: PU, topPt, triggerSF, leptonSF...
    // double  wgt_topPtSF = eve->wgt_topPt_;
    double Xsec = mySample_xSec_;//eve->wgt_xs_;
    double nGen = ( maxNentries>0 ) ? maxNentries : mySample_nGen_;//eve->wgt_nGen_;
    if( maxNentries==1000000 && insample==2300 ) nGen = 670121;
    double lumi = ( intLumi > 0 ) ? intLumi : 10000 ;

    double wgt_lumi = ( insample > 0 ) ? lumi * (Xsec/nGen) : 1;//"weight_PU*topPtWgt*osTriggerSF*lepIDAndIsoSF*"; // various weights

    h_numEvents->Fill(0.5,1.);
    h_numEvents->Fill(1.5,wgt_gen);

    h_numEvents_wgt->Fill(0.5,1.);
    h_numEvents_wgt->Fill(1.5,wgt_gen);
    h_numEvents_wgt->Fill(2.5,wgt_lumi);
    h_numEvents_wgt->Fill(3.5,wgt_gen * wgt_lumi);
    h_numEvents_wgt->Fill(4.5,wgt_gen * wgt_lumi * wgt_pu);
    h_numEvents_wgt->Fill(5.5,wgt_gen * wgt_lumi * wgt_pu * wgt_csv_noSys);


    h_numEvents_wgtCSV->Fill(0.5,wgt_gen*wgt_csv_noSys);
    h_numEvents_wgtCSV_hf->Fill(0.5,wgt_gen*wgt_csv_hf_noSys);
    h_numEvents_wgtCSV_lf->Fill(0.5,wgt_gen*wgt_csv_lf_noSys);
    h_numEvents_wgtCSV_cf->Fill(0.5,wgt_gen*wgt_csv_cf_noSys);

    h_numEvents_wgtCSV_partonFlavourNoPU->Fill(0.5,wgt_gen*wgt_csv_noSys_partonFlavourNoPU);
    h_numEvents_wgtCSV_hf_partonFlavourNoPU->Fill(0.5,wgt_gen*wgt_csv_hf_noSys_partonFlavourNoPU);
    h_numEvents_wgtCSV_lf_partonFlavourNoPU->Fill(0.5,wgt_gen*wgt_csv_lf_noSys_partonFlavourNoPU);
    h_numEvents_wgtCSV_cf_partonFlavourNoPU->Fill(0.5,wgt_gen*wgt_csv_cf_noSys_partonFlavourNoPU);

    h_numEvents_nowgtCSV->Fill(0.5,wgt_gen);

    // std::cout << " ************************************************************* " << std::endl;

    // printf("  wgt_CSVL = %.3f, nowgt_CSVL = %.3f, ratioL = %.3f \n", wgt_jet_CSVL, nowgt_jet_CSVL, wgt_jet_CSVL/nowgt_jet_CSVL );
    // printf("  wgt_CSVM = %.3f, nowgt_CSVM = %.3f, ratioM = %.3f \n", wgt_jet_CSVM, nowgt_jet_CSVM, wgt_jet_CSVM/nowgt_jet_CSVM );
    // printf("  wgt_CSVT = %.3f, nowgt_CSVT = %.3f, ratioT = %.3f \n", wgt_jet_CSVT, nowgt_jet_CSVT, wgt_jet_CSVT/nowgt_jet_CSVT );


    h_event_wgt_csv->Fill(wgt_gen*wgt_csv_noSys);
    h_event_wgt_csv_hf->Fill(wgt_gen*wgt_csv_hf_noSys);
    h_event_wgt_csv_lf->Fill(wgt_gen*wgt_csv_lf_noSys);
    h_event_wgt_csv_cf->Fill(wgt_gen*wgt_csv_cf_noSys);

    if( njet>=4 ){
      h_numEvents_wgtCSV->Fill(1.5,wgt_gen*wgt_csv_noSys);
      h_numEvents_wgtCSV_hf->Fill(1.5,wgt_gen*wgt_csv_hf_noSys);
      h_numEvents_wgtCSV_lf->Fill(1.5,wgt_gen*wgt_csv_lf_noSys);
      h_numEvents_wgtCSV_cf->Fill(1.5,wgt_gen*wgt_csv_cf_noSys);

      h_numEvents_nowgtCSV->Fill(1.5,wgt_gen);

      h_event_wgt_csv_4j->Fill(wgt_gen*wgt_csv_noSys);
      h_event_wgt_csv_4j_hf->Fill(wgt_gen*wgt_csv_hf_noSys);
      h_event_wgt_csv_4j_lf->Fill(wgt_gen*wgt_csv_lf_noSys);
      h_event_wgt_csv_4j_cf->Fill(wgt_gen*wgt_csv_cf_noSys);
    }

      
    double wgt = wgt_gen * wgt_lumi * wgt_pu;
    double wgtCSV = wgt * wgt_csv_noSys;

    numEvents_all += 1.;
    numEvents_wgt_gen += wgt_gen;
    numEvents_wgt_lumi += wgt_lumi;
    numEvents_wgt_gen_lumi += wgt_gen * wgt_lumi;
    numEvents_wgt_gen_lumi_pu += wgt_gen * wgt_lumi * wgt_pu;
    numEvents_wgt_gen_lumi_pu_csv += wgt_gen * wgt_lumi * wgt_pu * wgt_csv_noSys;


    h_numTruePVs->Fill(numTruePVs,wgt_gen);

    h_lheHT_wgt->Fill(lheHT,wgt_gen * wgt_lumi);


    // Skip HT range, if necessary
    if( useHTbins_ && (insample==2300 || insample==2305 || insample==2400 || insample==2405) ){
      if( lheHT > 100 ) continue;
    }

    h_lheHT_wgt_afterHT->Fill(lheHT,wgt_gen * wgt_lumi);

    if( verbose_ ) std::cout << " ===> test 3 " << std::endl;


    ///////////////////
    ////// selections
    ///////////////////

    h_event_selection[0]->Fill(0.5, wgt); // all
    h_mu_event_selection[0]->Fill(0.5, wgt); // all
    h_ele_event_selection[0]->Fill(0.5, wgt); // all

    h_event_selection_wgtCSV[0]->Fill(0.5, wgtCSV); // all
    h_mu_event_selection_wgtCSV[0]->Fill(0.5, wgtCSV); // all
    h_ele_event_selection_wgtCSV[0]->Fill(0.5, wgtCSV); // all

    if( verbose_ ) std::cout << " ===> test 3.0.0 " << std::endl;

    // require first PV is good
    bool goodFirstVertex = ( eve->goodFirstVertex_ );
    if( !goodFirstVertex ) continue;
    

    bool pass_trigger_ele = false;
    // if( insample<0 ) pass_trigger_ele = ( eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_==1 );
    // else             pass_trigger_ele = ( eve->pass_HLT_Ele27_WP85_Gsf_v_==1 );
    pass_trigger_ele = ( eve->pass_HLT_Ele27_eta2p1_WPLoose_Gsf_v_==1 );
    
    bool pass_trigger_mu  = false;
    // if( insample<0 ) pass_trigger_mu = ( eve->pass_HLT_IsoMu20_v_==1 );
    // else             pass_trigger_mu = ( eve->pass_HLT_IsoMu20_v_==1 );
    pass_trigger_mu = ( eve->pass_HLT_IsoMu20_v_==1 );


    if( insample<0 ){
      pass_trigger_ele = ( insample==-11 && pass_trigger_ele );
      pass_trigger_mu  = ( insample==-13 && pass_trigger_mu  );
    }

    if( pass_trigger_mu || pass_trigger_ele ) h_event_selection[0]->Fill(0.5+1, wgt); // all
    if( pass_trigger_mu )  h_mu_event_selection[0]->Fill(0.5+1, wgt); // all
    if( pass_trigger_ele ) h_ele_event_selection[0]->Fill(0.5+1, wgt); // all

    if( pass_trigger_mu || pass_trigger_ele ) h_event_selection_wgtCSV[0]->Fill(0.5+1, wgtCSV); // all
    if( pass_trigger_mu )  h_mu_event_selection_wgtCSV[0]->Fill(0.5+1, wgtCSV); // all
    if( pass_trigger_ele ) h_ele_event_selection_wgtCSV[0]->Fill(0.5+1, wgtCSV); // all

    // Pass Trigger selection
    if( !(pass_trigger_mu || pass_trigger_ele) ) continue;

    // require at least one PV
    if( numPVs<1 ) continue;


    h_numPVs->Fill(numPVs,wgt);


    h_numPVs_wgt[0]->Fill(numPVs,wgt_gen * wgt_lumi * wgt_pu);
    h_numPVs_wgt[3]->Fill(numPVs,wgt_gen * wgt_lumi * wgt_pu_up);
    h_numPVs_wgt[4]->Fill(numPVs,wgt_gen * wgt_lumi * wgt_pu_down);

    h_numPVs_noPUwgt[0]->Fill(numPVs,wgt_gen * wgt_lumi);
    h_numPVs_noPUwgt[3]->Fill(numPVs,wgt_gen * wgt_lumi);
    h_numPVs_noPUwgt[4]->Fill(numPVs,wgt_gen * wgt_lumi);

    for( int iSys=0; iSys<NumSysCat; iSys++ ){
      TString sys_name = sys_cat_labels[iSys];

      double use_wgt = wgt_gen * wgt_lumi * wgt_pu;
      if( sys_name.Contains("PUUp") )        use_wgt = wgt_gen * wgt_lumi * wgt_pu_up;
      else if( sys_name.Contains("PUDown") ) use_wgt = wgt_gen * wgt_lumi * wgt_pu_down;

      if( insample<0 ) use_wgt = 1.0;
      
      h_PV_x_wgt[iSys]->Fill(eve->PV_x_,use_wgt);
      h_PV_y_wgt[iSys]->Fill(eve->PV_y_,use_wgt);
      h_PV_z_wgt[iSys]->Fill(eve->PV_z_,use_wgt);
    }

    int numEle = int( ind_ele.size() );
    int numMu  = int( ind_mu.size() );

    int numEle_loose = int( ind_ele_loose.size() );
    int numMu_loose  = int( ind_mu_loose.size() );

    bool oneEle = ( numEle==1 && numEle_loose==1 && numMu_loose==0 );
    bool oneMu  = ( numMu==1 && numMu_loose==1 && numEle_loose==0 );

    oneEle = ( oneEle && pass_trigger_ele );
    oneMu  = ( oneMu  && pass_trigger_mu );

    bool oneLep = ( oneEle || oneMu );

    if( oneLep ) h_event_selection[0]->Fill(0.5+2, wgt); // all
    if( oneMu )  h_mu_event_selection[0]->Fill(0.5+2, wgt); // all
    if( oneEle ) h_ele_event_selection[0]->Fill(0.5+2, wgt); // all

    if( oneLep ) h_event_selection_wgtCSV[0]->Fill(0.5+2, wgtCSV); // all
    if( oneMu )  h_mu_event_selection_wgtCSV[0]->Fill(0.5+2, wgtCSV); // all
    if( oneEle ) h_ele_event_selection_wgtCSV[0]->Fill(0.5+2, wgtCSV); // all


    bool TwoEle = ( numEle>=1 && numEle_loose==2 && numMu_loose==0 );
    bool TwoMu  = ( numMu>=1 && numMu_loose==2 && numEle_loose==0 );

    TwoEle = ( TwoEle && pass_trigger_ele );
    TwoMu  = ( TwoMu  && pass_trigger_mu );

    bool TwoLep = ( TwoEle || TwoMu );

    bool EleMu = ( (numEle==1 || numMu==1) && numMu_loose==1 && numEle_loose==1 );

    EleMu = ( EleMu && ((numEle==1 && pass_trigger_ele) || (numMu==1 && pass_trigger_mu)) );
    if( verbose_ ) std::cout << " ===> test 3.0.1 " << std::endl;

    // Require at least two leptons
    if( TwoLep ){
      if( verbose_ ) std::cout << " ===> test 3.0 " << std::endl;

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


      /////
      /////  Loop over systematics
      /////
      if( verbose_ ) std::cout << " ===> test 4 " << std::endl;

      for( int iSys=0; iSys<NumSysCat; iSys++ ){
	TString sys_name = sys_cat_labels[iSys];

	vdouble use_jet_pt = eve->jet_pt_;
	vdouble use_jet_eta = eve->jet_eta_;
	vdouble use_jet_phi = eve->jet_phi_;
	vdouble use_jet_energy = eve->jet_energy_;

	vdouble use_jet_csv = eve->jet_csv_;
	vdouble use_jet_pileupJetId_fullDiscriminant = eve->jet_pileupJetId_fullDiscriminant_;

	vint use_jet_hadronFlavour = eve->jet_hadronFlavour_;


	double pfMET = eve->pfMET_pt_;
	double pfMET_phi = eve->pfMET_phi_;

	double pfMETNoHF = eve->pfMETNoHF_pt_;
	double pfMETNoHF_phi = eve->pfMETNoHF_phi_;

	double puppiMET = eve->puppiMET_pt_;
	double puppiMET_phi = eve->puppiMET_phi_;


	if( sys_name.Contains("JESUp") ){
	  use_jet_pt = eve->jet_JESup_pt_;
	  use_jet_eta = eve->jet_JESup_eta_;
	  use_jet_phi = eve->jet_JESup_phi_;
	  use_jet_energy = eve->jet_JESup_energy_;

	  use_jet_csv = eve->jet_JESup_csv_;
	  use_jet_pileupJetId_fullDiscriminant = eve->jet_JESup_pileupJetId_fullDiscriminant_;

	  use_jet_hadronFlavour = eve->jet_JESup_hadronFlavour_;


	  pfMET = eve->pfMET_pt_JESup_;
	  pfMET_phi = eve->pfMET_phi_JESup_;

	  pfMETNoHF = eve->pfMETNoHF_pt_JESup_;
	  pfMETNoHF_phi = eve->pfMETNoHF_phi_JESup_;

	  puppiMET = eve->puppiMET_pt_JESup_;
	  puppiMET_phi = eve->puppiMET_phi_JESup_;
	}
	else if( sys_name.Contains("JESDown") ){
	  use_jet_pt = eve->jet_JESdown_pt_;
	  use_jet_eta = eve->jet_JESdown_eta_;
	  use_jet_phi = eve->jet_JESdown_phi_;
	  use_jet_energy = eve->jet_JESdown_energy_;

	  use_jet_csv = eve->jet_JESdown_csv_;
	  use_jet_pileupJetId_fullDiscriminant = eve->jet_JESdown_pileupJetId_fullDiscriminant_;

	  use_jet_hadronFlavour = eve->jet_JESdown_hadronFlavour_;

	  pfMET = eve->pfMET_pt_JESdown_;
	  pfMET_phi = eve->pfMET_phi_JESdown_;

	  pfMETNoHF = eve->pfMETNoHF_pt_JESdown_;
	  pfMETNoHF_phi = eve->pfMETNoHF_phi_JESdown_;

	  puppiMET = eve->puppiMET_pt_JESdown_;
	  puppiMET_phi = eve->puppiMET_phi_JESdown_;
	}
	else if( sys_name.Contains("JERUp") ){
	  use_jet_pt = eve->jet_JERup_pt_;
	  use_jet_eta = eve->jet_JERup_eta_;
	  use_jet_phi = eve->jet_JERup_phi_;
	  use_jet_energy = eve->jet_JERup_energy_;

	  use_jet_csv = eve->jet_JERup_csv_;
	  use_jet_pileupJetId_fullDiscriminant = eve->jet_JERup_pileupJetId_fullDiscriminant_;

	  use_jet_hadronFlavour = eve->jet_JERup_hadronFlavour_;

	  pfMET = eve->pfMET_pt_JERup_;
	  pfMET_phi = eve->pfMET_phi_JERup_;

	  pfMETNoHF = eve->pfMETNoHF_pt_JERup_;
	  pfMETNoHF_phi = eve->pfMETNoHF_phi_JERup_;

	  puppiMET = eve->puppiMET_pt_JERup_;
	  puppiMET_phi = eve->puppiMET_phi_JERup_;
	}
	else if( sys_name.Contains("JERDown") ){
	  use_jet_pt = eve->jet_JERdown_pt_;
	  use_jet_eta = eve->jet_JERdown_eta_;
	  use_jet_phi = eve->jet_JERdown_phi_;
	  use_jet_energy = eve->jet_JERdown_energy_;

	  use_jet_csv = eve->jet_JERdown_csv_;
	  use_jet_pileupJetId_fullDiscriminant = eve->jet_JERdown_pileupJetId_fullDiscriminant_;

	  use_jet_hadronFlavour = eve->jet_JERdown_hadronFlavour_;

	  pfMET = eve->pfMET_pt_JERdown_;
	  pfMET_phi = eve->pfMET_phi_JERdown_;

	  pfMETNoHF = eve->pfMETNoHF_pt_JERdown_;
	  pfMETNoHF_phi = eve->pfMETNoHF_phi_JERdown_;

	  puppiMET = eve->puppiMET_pt_JERdown_;
	  puppiMET_phi = eve->puppiMET_phi_JERdown_;
	}


	// Set lepton scale factor equal to 1 for now
	double wgt_lepIdSF = 1.;
	if( sys_name.Contains("lepIdSFUp") )        wgt_lepIdSF = 1.;
	else if( sys_name.Contains("lepIdSFDown") ) wgt_lepIdSF = 1.;

	// PU systematic
	double use_wgt_pu = wgt_pu;
	if( sys_name.Contains("PUUp") )        use_wgt_pu = wgt_pu_up;
	else if( sys_name.Contains("PUDown") ) use_wgt_pu = wgt_pu_down;


    if( verbose_ ) std::cout << " ===> test 5 " << std::endl;

	// Factorization and normalization scale uncertainty
	double wgt_LHEscale = 1.;
	if( insample>-1 ){
	  double originalXWGTUP = eve->originalXWGTUP_;
	  vdouble lhe_event_weights = eve->LHEEvent_weights_;

	  if( lhe_event_weights.size()>76 && originalXWGTUP!=0 ){
	    int use_weight_index = 0;

	    if( sys_name.Contains("muFUp") )        use_weight_index = 1;
	    else if( sys_name.Contains("muFDown") ) use_weight_index = 2;
	    else if( sys_name.Contains("muRUp") )   use_weight_index = 3;
	    else if( sys_name.Contains("muRDown") ) use_weight_index = 6;
	    else if( sys_name.Contains("muRmuFUp") )   use_weight_index = 4;
	    else if( sys_name.Contains("muRmuFDown") ) use_weight_index = 8;
	    else if( sys_name.Contains("muRmuFDown") ) use_weight_index = 8;
	    else if( sys_name.Contains("NNPDFUp") ) use_weight_index = 75;
	    else if( sys_name.Contains("NNPDFDown") ) use_weight_index = 13;

	    wgt_LHEscale = lhe_event_weights[use_weight_index] / originalXWGTUP;
	  }
	}


	//wgt = wgt_gen * wgt_lumi * wgt_pu * wgt_lepIdSF * wgt_LHEscale;
	wgt = wgt_gen * wgt_lumi * use_wgt_pu * wgt_lepIdSF * wgt_LHEscale;


	double fill_pfMETNoHF = std::min( metmax-0.001, pfMETNoHF);

	// Fill MET histograms
	h_pfMETNoHF_pt_2l[iSys]->Fill(fill_pfMETNoHF,wgt);
	if( TwoEle )     h_pfMETNoHF_pt_2e[iSys]->Fill(fill_pfMETNoHF,wgt);
	else if( TwoMu ) h_pfMETNoHF_pt_2m[iSys]->Fill(fill_pfMETNoHF,wgt);


	vdouble jetPts;
	vdouble jetEtas;
	vdouble jetCSVs;
	vint    jetFlavors;

	vdouble temp_test_jetPts;
	vdouble temp_test_jetEtas;
	vdouble temp_test_jetCSVs;
	vint    temp_test_jetFlavors;

	vint ind_jet;

	int numJet = 0;
	int numBtag = 0;

	int numJet_b = 0;
	int numJet_c = 0;
	int numJet_l = 0;

	double HT30=0;
	TLorentzVector sumJet;
	bool firstJet = false;

	double wgt_csv2=1.;
	double wgt_csv2_hf=1., wgt_csv2_lf=1., wgt_csv2_cf=1.;

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

	  temp_test_jetPts.clear();
	  temp_test_jetEtas.clear();
	  temp_test_jetCSVs.clear();
	  temp_test_jetFlavors.clear();

	  temp_test_jetPts.push_back(pt);
	  temp_test_jetEtas.push_back(eta);
	  temp_test_jetCSVs.push_back(csv);
	  temp_test_jetFlavors.push_back(flavor);


	  jetPts.push_back(pt);
	  jetEtas.push_back(eta);
	  jetCSVs.push_back(csv);
	  jetFlavors.push_back(flavor);

	  numJet++;
	  if( csv > 0.800 ){
	    numBtag++;
	  }

	  if( pt > 1000 ) pt = 999.;

	  double my_jet_sf = 1.;
	  BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;

	  if( abs(flavor)==5 )      jf = BTagEntry::FLAV_B;
	  else if( abs(flavor)==4 ) jf = BTagEntry::FLAV_C;
	  else                      jf = BTagEntry::FLAV_UDSG;

	  bool isBFlav = false;
	  bool isCFlav = false;
	  bool isLFlav = false;
	  if( abs(flavor)==5 )      isBFlav = true;
	  else if( abs(flavor)==4 ) isCFlav = true;
	  else                      isLFlav = true;

	  if( abs(flavor)==5 )      numJet_b++;
	  else if( abs(flavor)==4 ) numJet_c++;
	  else                      numJet_l++;


	  if( sys_name.Contains("JESUp") )        my_jet_sf = reader_JESUp.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("JESDown") ) my_jet_sf = reader_JESDown.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("LFUp") && isBFlav )    my_jet_sf = reader_LFUp.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("LFDown") && isBFlav )  my_jet_sf = reader_LFDown.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("HFUp") && isLFlav )    my_jet_sf = reader_HFUp.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("HFDown") && isLFlav )  my_jet_sf = reader_HFDown.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVHFStats1Up") && isBFlav )   my_jet_sf = reader_HFStats1Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVHFStats1Down") && isBFlav ) my_jet_sf = reader_HFStats1Down.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVHFStats2Up") && isBFlav )   my_jet_sf = reader_HFStats2Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVHFStats2Down") && isBFlav ) my_jet_sf = reader_HFStats2Down.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVLFStats1Up") && isLFlav )   my_jet_sf = reader_LFStats1Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVLFStats1Down") && isLFlav ) my_jet_sf = reader_LFStats1Down.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVLFStats2Up") && isLFlav )   my_jet_sf = reader_LFStats2Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVLFStats2Down") && isLFlav ) my_jet_sf = reader_LFStats2Down.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVCFErr1Up") && isCFlav )   my_jet_sf = reader_CFErr1Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVCFErr1Down") && isCFlav ) my_jet_sf = reader_CFErr1Down.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVCFErr2Up") && isCFlav )   my_jet_sf = reader_CFErr2Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVCFErr2Down") && isCFlav ) my_jet_sf = reader_CFErr2Down.eval(jf, eta, pt, csv);
	  else my_jet_sf = reader.eval(jf, eta, pt, csv);


	 
	  double my_jet_sf_hf = 1.;
	  double my_jet_sf_lf = 1.;
	  double my_jet_sf_cf = 1.;
	  if( abs(flavor)==5 ){
	    if( sys_name.Contains("JESUp") )        my_jet_sf_hf = reader_JESUp.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("JESDown") ) my_jet_sf_hf = reader_JESDown.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("LFUp") )    my_jet_sf_hf = reader_LFUp.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("LFDown") )  my_jet_sf_hf = reader_LFDown.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("CSVHFStats1Up") )   my_jet_sf_hf = reader_HFStats1Up.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("CSVHFStats1Down") ) my_jet_sf_hf = reader_HFStats1Down.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("CSVHFStats2Up") )   my_jet_sf_hf = reader_HFStats2Up.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("CSVHFStats2Down") ) my_jet_sf_hf = reader_HFStats2Down.eval(jf, eta, pt, csv);
	    else my_jet_sf_hf = reader.eval(jf, eta, pt, csv);
	  }
	  else if( abs(flavor)==4 ){
	    if( sys_name.Contains("JESUp") )        my_jet_sf_cf = reader_JESUp.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("JESDown") ) my_jet_sf_cf = reader_JESDown.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("CSVCFErr1Up") )   my_jet_sf_cf = reader_CFErr1Up.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("CSVCFErr1Down") ) my_jet_sf_cf = reader_CFErr1Down.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("CSVCFErr2Up") )   my_jet_sf_cf = reader_CFErr2Up.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("CSVCFErr2Down") ) my_jet_sf_cf = reader_CFErr2Down.eval(jf, eta, pt, csv);
	    else my_jet_sf_cf = reader.eval(jf, eta, pt, csv);
	  }
	  else {
	    if( sys_name.Contains("JESUp") )        my_jet_sf_lf = reader_JESUp.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("JESDown") ) my_jet_sf_lf = reader_JESDown.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("HFUp") )    my_jet_sf_lf = reader_HFUp.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("HFDown") )  my_jet_sf_lf = reader_HFDown.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("CSVLFStats1Up") )   my_jet_sf_lf = reader_LFStats1Up.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("CSVLFStats1Down") ) my_jet_sf_lf = reader_LFStats1Down.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("CSVLFStats2Up") )   my_jet_sf_lf = reader_LFStats2Up.eval(jf, eta, pt, csv);
	    else if( sys_name.Contains("CSVLFStats2Down") ) my_jet_sf_lf = reader_LFStats2Down.eval(jf, eta, pt, csv);
	    else my_jet_sf_lf = reader.eval(jf, eta, pt, csv);
	  }


	  assert (my_jet_sf > 0.);
	  wgt_csv2 *= my_jet_sf;

	  wgt_csv2_hf *= my_jet_sf_hf;
	  wgt_csv2_cf *= my_jet_sf_cf;
	  wgt_csv2_lf *= my_jet_sf_lf;

	    
	  double temp_test_wgt_csv_hf=1, temp_test_wgt_csv_lf=1, temp_test_wgt_csv_cf=1;
	  double temp_test_wgt_csv = get_csv_wgt(temp_test_jetPts, temp_test_jetEtas, temp_test_jetCSVs, temp_test_jetFlavors, iSys, 
						 temp_test_wgt_csv_hf, temp_test_wgt_csv_lf, temp_test_wgt_csv_cf);

	  if( iSys==0 ){
	    if( abs(flavor)==5 ){
	      h_b_csv_wgt->Fill(csv,temp_test_wgt_csv);
	      h_b_csv_wgt2->Fill(csv,my_jet_sf);
	      h_b_csv_diff_wgt2_wgt->Fill(csv,my_jet_sf-temp_test_wgt_csv);
	    }
	    else if( abs(flavor)==4 ){
	      h_c_csv_wgt->Fill(csv,temp_test_wgt_csv);
	      h_c_csv_wgt2->Fill(csv,my_jet_sf);
	      h_c_csv_diff_wgt2_wgt->Fill(csv,my_jet_sf-temp_test_wgt_csv);
	    }
	    else {
	      h_l_csv_wgt->Fill(csv,temp_test_wgt_csv);
	      h_l_csv_wgt2->Fill(csv,my_jet_sf);
	      h_l_csv_diff_wgt2_wgt->Fill(csv,my_jet_sf-temp_test_wgt_csv);
	    }
	  }

	  if( sys_name.Contains("JESUp") && false ){
	    printf("  iJet = %d: pt = %4.1f, eta = %+.2f, flavor = %+d, CSVv2 = %.2f: weight csv file = %.3f, weight hist file = %.3f \n",
		   iJet, pt, eta, flavor, csv, my_jet_sf, temp_test_wgt_csv);
	  }
	}

    if( verbose_ ) std::cout << " ===> test 6 " << std::endl;



	// Calculate CSV weight
	double wgt_csv_hf, wgt_csv_lf, wgt_csv_cf;
	double wgt_csv = ( insample<0 ) ? 1 : get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);

	wgtCSV = wgt * wgt_csv;


	bool use_verbose = false;
	if( numJet>=1 && numBtag>=1 ) use_verbose = !true;
	double wgt_csv3 = ( insample<0 ) ? 1 : get_btv_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, sys_name, true, use_verbose);
	double wgt_csv4 = ( insample<0 ) ? 1 : get_btv_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, sys_name, false, false);

	if( insample<0 ) wgt_csv2 = 1;
	if( insample<0 ) wgt_csv3 = 1;
	if( insample<0 ) wgt_csv4 = 1;

	double wgtCSV2 = wgt * wgt_csv2;

	double wgtCSV3 = wgt * wgt_csv3;

	double wgtCSV4 = wgt * wgt_csv4;



	int fill_numJet  = std::min( MaxNjet-1, numJet );
	int fill_numBtag = std::min( MaxNbtag-1, numBtag );

	h_numJet_2l[iSys]->Fill(fill_numJet,wgt);
	h_numJet_wgtCSV_2l[iSys]->Fill(fill_numJet,wgtCSV);
	h_numJet_wgtCSV2_2l[iSys]->Fill(fill_numJet,wgtCSV2);
	h_numJet_wgtCSV3_2l[iSys]->Fill(fill_numJet,wgtCSV3);
	h_numJet_wgtCSV4_2l[iSys]->Fill(fill_numJet,wgtCSV4);

	h_numBtag_2l[iSys]->Fill(fill_numBtag,wgt);
	h_numBtag_wgtCSV_2l[iSys]->Fill(fill_numBtag,wgtCSV);
	h_numBtag_wgtCSV2_2l[iSys]->Fill(fill_numBtag,wgtCSV2);
	h_numBtag_wgtCSV3_2l[iSys]->Fill(fill_numBtag,wgtCSV3);
	h_numBtag_wgtCSV4_2l[iSys]->Fill(fill_numBtag,wgtCSV4);

	if( pfMETNoHF<30 ){
	  h_numJet_2l_met30[iSys]->Fill(fill_numJet,wgt);
	  h_numJet_wgtCSV_2l_met30[iSys]->Fill(fill_numJet,wgtCSV);
	  h_numJet_wgtCSV2_2l_met30[iSys]->Fill(fill_numJet,wgtCSV2);
	  h_numJet_wgtCSV3_2l_met30[iSys]->Fill(fill_numJet,wgtCSV3);
	  h_numJet_wgtCSV4_2l_met30[iSys]->Fill(fill_numJet,wgtCSV4);

	  h_numBtag_2l_met30[iSys]->Fill(fill_numBtag,wgt);
	  h_numBtag_wgtCSV_2l_met30[iSys]->Fill(fill_numBtag,wgtCSV);
	  h_numBtag_wgtCSV2_2l_met30[iSys]->Fill(fill_numBtag,wgtCSV2);
	  h_numBtag_wgtCSV3_2l_met30[iSys]->Fill(fill_numBtag,wgtCSV3);
	  h_numBtag_wgtCSV4_2l_met30[iSys]->Fill(fill_numBtag,wgtCSV4);
	}



	if( numJet<=2 ){
	  h_numBtag_2l_leq2j[iSys]->Fill(fill_numBtag,wgt);
	  h_numBtag_wgtCSV_2l_leq2j[iSys]->Fill(fill_numBtag,wgtCSV);
	  h_numBtag_wgtCSV2_2l_leq2j[iSys]->Fill(fill_numBtag,wgtCSV2);
	  h_numBtag_wgtCSV3_2l_leq2j[iSys]->Fill(fill_numBtag,wgtCSV3);
	  h_numBtag_wgtCSV4_2l_leq2j[iSys]->Fill(fill_numBtag,wgtCSV4);

	  // 0 jet
	  if( numJet_b==0 && numJet_c==0 && numJet_l==0 ){
	    h_numBtag_2l_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }

	  // 1 jet
	  if( numJet_b==1 && numJet_c==0 && numJet_l==0 ){
	    h_numBtag_2l_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==0 && numJet_c==1 && numJet_l==0 ){
	    h_numBtag_2l_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==0 && numJet_c==0 && numJet_l==1 ){
	    h_numBtag_2l_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }

	  // 2 jet
	  if( numJet_b==2 && numJet_c==0 && numJet_l==0 ){
	    h_numBtag_2l_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==0 && numJet_c==2 && numJet_l==0 ){
	    h_numBtag_2l_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==0 && numJet_c==0 && numJet_l==2 ){
	    h_numBtag_2l_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==1 && numJet_c==1 && numJet_l==0 ){
	    h_numBtag_2l_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==1 && numJet_c==0 && numJet_l==1 ){
	    h_numBtag_2l_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==0 && numJet_c==1 && numJet_l==1 ){
	    h_numBtag_2l_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }


	  if( pfMETNoHF<30 ){
	    h_numBtag_2l_met30_leq2j[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_met30_leq2j[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_met30_leq2j[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_met30_leq2j[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_met30_leq2j[iSys]->Fill(fill_numBtag,wgtCSV4);

	    // 0 jet
	    if( numJet_b==0 && numJet_c==0 && numJet_l==0 ){
	      h_numBtag_2l_met30_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }

	    // 1 jet
	    if( numJet_b==1 && numJet_c==0 && numJet_l==0 ){
	      h_numBtag_2l_met30_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==0 && numJet_c==1 && numJet_l==0 ){
	      h_numBtag_2l_met30_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==0 && numJet_c==0 && numJet_l==1 ){
	      h_numBtag_2l_met30_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }

	    // 2 jet
	    if( numJet_b==2 && numJet_c==0 && numJet_l==0 ){
	      h_numBtag_2l_met30_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==0 && numJet_c==2 && numJet_l==0 ){
	      h_numBtag_2l_met30_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==0 && numJet_c==0 && numJet_l==2 ){
	      h_numBtag_2l_met30_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==1 && numJet_c==1 && numJet_l==0 ){
	      h_numBtag_2l_met30_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==1 && numJet_c==0 && numJet_l==1 ){
	      h_numBtag_2l_met30_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==0 && numJet_c==1 && numJet_l==1 ){
	      h_numBtag_2l_met30_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	  }
	}




	if( numJet>=1 && numJet<=2 ){
	  h_numBtag_2l_geq1j_leq2j[iSys]->Fill(fill_numBtag,wgt);
	  h_numBtag_wgtCSV_2l_geq1j_leq2j[iSys]->Fill(fill_numBtag,wgtCSV);
	  h_numBtag_wgtCSV2_2l_geq1j_leq2j[iSys]->Fill(fill_numBtag,wgtCSV2);
	  h_numBtag_wgtCSV3_2l_geq1j_leq2j[iSys]->Fill(fill_numBtag,wgtCSV3);
	  h_numBtag_wgtCSV4_2l_geq1j_leq2j[iSys]->Fill(fill_numBtag,wgtCSV4);

	  // 0 jet
	  if( numJet_b==0 && numJet_c==0 && numJet_l==0 ){
	    h_numBtag_2l_geq1j_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_geq1j_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }

	  // 1 jet
	  if( numJet_b==1 && numJet_c==0 && numJet_l==0 ){
	    h_numBtag_2l_geq1j_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_geq1j_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_geq1j_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_geq1j_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_geq1j_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==0 && numJet_c==1 && numJet_l==0 ){
	    h_numBtag_2l_geq1j_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_geq1j_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==0 && numJet_c==0 && numJet_l==1 ){
	    h_numBtag_2l_geq1j_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_geq1j_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }

	  // 2 jet
	  if( numJet_b==2 && numJet_c==0 && numJet_l==0 ){
	    h_numBtag_2l_geq1j_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_geq1j_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_geq1j_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_geq1j_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_geq1j_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==0 && numJet_c==2 && numJet_l==0 ){
	    h_numBtag_2l_geq1j_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_geq1j_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==0 && numJet_c==0 && numJet_l==2 ){
	    h_numBtag_2l_geq1j_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_geq1j_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==1 && numJet_c==1 && numJet_l==0 ){
	    h_numBtag_2l_geq1j_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_geq1j_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_geq1j_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_geq1j_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_geq1j_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==1 && numJet_c==0 && numJet_l==1 ){
	    h_numBtag_2l_geq1j_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_geq1j_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_geq1j_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_geq1j_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_geq1j_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	  if( numJet_b==0 && numJet_c==1 && numJet_l==1 ){
	    h_numBtag_2l_geq1j_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_geq1j_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_geq1j_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_geq1j_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_geq1j_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }


	  if( pfMETNoHF<30 ){
	    h_numBtag_2l_met30_geq1j_leq2j[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_2l_met30_geq1j_leq2j[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j[iSys]->Fill(fill_numBtag,wgtCSV4);

	    // 0 jet
	    if( numJet_b==0 && numJet_c==0 && numJet_l==0 ){
	      h_numBtag_2l_met30_geq1j_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b0c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }

	    // 1 jet
	    if( numJet_b==1 && numJet_c==0 && numJet_l==0 ){
	      h_numBtag_2l_met30_geq1j_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_1b0c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==0 && numJet_c==1 && numJet_l==0 ){
	      h_numBtag_2l_met30_geq1j_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b1c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==0 && numJet_c==0 && numJet_l==1 ){
	      h_numBtag_2l_met30_geq1j_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b0c1l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }

	    // 2 jet
	    if( numJet_b==2 && numJet_c==0 && numJet_l==0 ){
	      h_numBtag_2l_met30_geq1j_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_2b0c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==0 && numJet_c==2 && numJet_l==0 ){
	      h_numBtag_2l_met30_geq1j_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b2c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==0 && numJet_c==0 && numJet_l==2 ){
	      h_numBtag_2l_met30_geq1j_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b0c2l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==1 && numJet_c==1 && numJet_l==0 ){
	      h_numBtag_2l_met30_geq1j_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_1b1c0l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==1 && numJet_c==0 && numJet_l==1 ){
	      h_numBtag_2l_met30_geq1j_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_1b0c1l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	    if( numJet_b==0 && numJet_c==1 && numJet_l==1 ){
	      h_numBtag_2l_met30_geq1j_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgt);
	      h_numBtag_wgtCSV_2l_met30_geq1j_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV);
	      h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV2);
	      h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV3);
	      h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j_0b1c1l[iSys]->Fill(fill_numBtag,wgtCSV4);
	    }
	  }
	}


      } // end loop over systematics
    } // end require exactly two leptons TwoLep
    if( verbose_ ) std::cout << " ===> test 3.0.2 " << std::endl;



    /// EleMu control region

    // Require at least two leptons
    if( EleMu ){
    if( verbose_ ) std::cout << " ===> test 3.1.0 " << std::endl;

      int lepInd1 = ind_ele_loose[0];
      int lepInd2 = ind_mu_loose[0];
    if( verbose_ ) std::cout << " ===> test 3.1.1 " << std::endl;

      TLorentzVector myLep1;
      myLep1.SetPtEtaPhiE( eve->lepton_pt_[lepInd1], eve->lepton_eta_[lepInd1], eve->lepton_phi_[lepInd1], eve->lepton_energy_[lepInd1] );

      TLorentzVector myLep2;
      myLep2.SetPtEtaPhiE( eve->lepton_pt_[lepInd2], eve->lepton_eta_[lepInd2], eve->lepton_phi_[lepInd2], eve->lepton_energy_[lepInd2] );

      int charge1 = eve->lepton_trkCharge_[lepInd1];
      int charge2 = eve->lepton_trkCharge_[lepInd2];
    if( verbose_ ) std::cout << " ===> test 3.1.2 " << std::endl;

      //
      // REQUIRE OPPOSITE CHARGE
      //

      if( !(charge1*charge2<0) ) continue;

      TLorentzVector diLep = myLep1 + myLep2;

      double diLepMass = diLep.M();

      //
      // REQUIRE LOOSE DILEPTON MASS
      //

      if( !(diLepMass > 20) ) continue;

      h_em_diLepMass->Fill(diLepMass,wgt);

    if( verbose_ ) std::cout << " ===> test 3.1.3 " << std::endl;

      /////
      /////  Loop over systematics
      /////

      for( int iSys=0; iSys<NumSysCat; iSys++ ){
	TString sys_name = sys_cat_labels[iSys];

	vdouble use_jet_pt = eve->jet_pt_;
	vdouble use_jet_eta = eve->jet_eta_;
	vdouble use_jet_phi = eve->jet_phi_;
	vdouble use_jet_energy = eve->jet_energy_;

	vdouble use_jet_csv = eve->jet_csv_;
	vdouble use_jet_pileupJetId_fullDiscriminant = eve->jet_pileupJetId_fullDiscriminant_;

	vint use_jet_hadronFlavour = eve->jet_hadronFlavour_;


	double pfMET = eve->pfMET_pt_;
	double pfMET_phi = eve->pfMET_phi_;

	double pfMETNoHF = eve->pfMETNoHF_pt_;
	double pfMETNoHF_phi = eve->pfMETNoHF_phi_;

	double puppiMET = eve->puppiMET_pt_;
	double puppiMET_phi = eve->puppiMET_phi_;


	if( sys_name.Contains("JESUp") ){
	  use_jet_pt = eve->jet_JESup_pt_;
	  use_jet_eta = eve->jet_JESup_eta_;
	  use_jet_phi = eve->jet_JESup_phi_;
	  use_jet_energy = eve->jet_JESup_energy_;

	  use_jet_csv = eve->jet_JESup_csv_;
	  use_jet_pileupJetId_fullDiscriminant = eve->jet_JESup_pileupJetId_fullDiscriminant_;

	  use_jet_hadronFlavour = eve->jet_JESup_hadronFlavour_;


	  pfMET = eve->pfMET_pt_JESup_;
	  pfMET_phi = eve->pfMET_phi_JESup_;

	  pfMETNoHF = eve->pfMETNoHF_pt_JESup_;
	  pfMETNoHF_phi = eve->pfMETNoHF_phi_JESup_;

	  puppiMET = eve->puppiMET_pt_JESup_;
	  puppiMET_phi = eve->puppiMET_phi_JESup_;
	}
	else if( sys_name.Contains("JESDown") ){
	  use_jet_pt = eve->jet_JESdown_pt_;
	  use_jet_eta = eve->jet_JESdown_eta_;
	  use_jet_phi = eve->jet_JESdown_phi_;
	  use_jet_energy = eve->jet_JESdown_energy_;

	  use_jet_csv = eve->jet_JESdown_csv_;
	  use_jet_pileupJetId_fullDiscriminant = eve->jet_JESdown_pileupJetId_fullDiscriminant_;

	  use_jet_hadronFlavour = eve->jet_JESdown_hadronFlavour_;

	  pfMET = eve->pfMET_pt_JESdown_;
	  pfMET_phi = eve->pfMET_phi_JESdown_;

	  pfMETNoHF = eve->pfMETNoHF_pt_JESdown_;
	  pfMETNoHF_phi = eve->pfMETNoHF_phi_JESdown_;

	  puppiMET = eve->puppiMET_pt_JESdown_;
	  puppiMET_phi = eve->puppiMET_phi_JESdown_;
	}
	else if( sys_name.Contains("JERUp") ){
	  use_jet_pt = eve->jet_JERup_pt_;
	  use_jet_eta = eve->jet_JERup_eta_;
	  use_jet_phi = eve->jet_JERup_phi_;
	  use_jet_energy = eve->jet_JERup_energy_;

	  use_jet_csv = eve->jet_JERup_csv_;
	  use_jet_pileupJetId_fullDiscriminant = eve->jet_JERup_pileupJetId_fullDiscriminant_;

	  use_jet_hadronFlavour = eve->jet_JERup_hadronFlavour_;

	  pfMET = eve->pfMET_pt_JERup_;
	  pfMET_phi = eve->pfMET_phi_JERup_;

	  pfMETNoHF = eve->pfMETNoHF_pt_JERup_;
	  pfMETNoHF_phi = eve->pfMETNoHF_phi_JERup_;

	  puppiMET = eve->puppiMET_pt_JERup_;
	  puppiMET_phi = eve->puppiMET_phi_JERup_;
	}
	else if( sys_name.Contains("JERDown") ){
	  use_jet_pt = eve->jet_JERdown_pt_;
	  use_jet_eta = eve->jet_JERdown_eta_;
	  use_jet_phi = eve->jet_JERdown_phi_;
	  use_jet_energy = eve->jet_JERdown_energy_;

	  use_jet_csv = eve->jet_JERdown_csv_;
	  use_jet_pileupJetId_fullDiscriminant = eve->jet_JERdown_pileupJetId_fullDiscriminant_;

	  use_jet_hadronFlavour = eve->jet_JERdown_hadronFlavour_;

	  pfMET = eve->pfMET_pt_JERdown_;
	  pfMET_phi = eve->pfMET_phi_JERdown_;

	  pfMETNoHF = eve->pfMETNoHF_pt_JERdown_;
	  pfMETNoHF_phi = eve->pfMETNoHF_phi_JERdown_;

	  puppiMET = eve->puppiMET_pt_JERdown_;
	  puppiMET_phi = eve->puppiMET_phi_JERdown_;
	}
    if( verbose_ ) std::cout << " ===> test 3.1.4 " << std::endl;


	// Set lepton scale factor equal to 1 for now
	double wgt_lepIdSF = 1.;
	if( sys_name.Contains("lepIdSFUp") )        wgt_lepIdSF = 1.;
	else if( sys_name.Contains("lepIdSFDown") ) wgt_lepIdSF = 1.;

	// PU systematic
	double use_wgt_pu = wgt_pu;
	if( sys_name.Contains("PUUp") )        use_wgt_pu = wgt_pu_up;
	else if( sys_name.Contains("PUDown") ) use_wgt_pu = wgt_pu_down;



	// Factorization and normalization scale uncertainty
	double wgt_LHEscale = 1.;
	if( insample>-1 ){
	  double originalXWGTUP = eve->originalXWGTUP_;
	  vdouble lhe_event_weights = eve->LHEEvent_weights_;

	  if( lhe_event_weights.size()>76 && originalXWGTUP!=0 ){
	    int use_weight_index = 0;

	    if( sys_name.Contains("muFUp") )        use_weight_index = 1;
	    else if( sys_name.Contains("muFDown") ) use_weight_index = 2;
	    else if( sys_name.Contains("muRUp") )   use_weight_index = 3;
	    else if( sys_name.Contains("muRDown") ) use_weight_index = 6;
	    else if( sys_name.Contains("muRmuFUp") )   use_weight_index = 4;
	    else if( sys_name.Contains("muRmuFDown") ) use_weight_index = 8;
	    else if( sys_name.Contains("muRmuFDown") ) use_weight_index = 8;
	    else if( sys_name.Contains("NNPDFUp") ) use_weight_index = 75;
	    else if( sys_name.Contains("NNPDFDown") ) use_weight_index = 13;

	    wgt_LHEscale = lhe_event_weights[use_weight_index] / originalXWGTUP;
	  }
	}


	//wgt = wgt_gen * wgt_lumi * wgt_pu * wgt_lepIdSF * wgt_LHEscale;
	wgt = wgt_gen * wgt_lumi * use_wgt_pu * wgt_lepIdSF * wgt_LHEscale;


	double fill_pfMETNoHF = std::min( metmax-0.001, pfMETNoHF);

	// Fill MET histograms
	h_pfMETNoHF_pt_em[iSys]->Fill(fill_pfMETNoHF,wgt);
    if( verbose_ ) std::cout << " ===> test 3.1.5 " << std::endl;

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

	double HT30=0;
	TLorentzVector sumJet;
	bool firstJet = false;

	double wgt_csv2=1.;
	double wgt_csv2_hf=1., wgt_csv2_lf=1., wgt_csv2_cf=1.;

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
	  if( csv > 0.800 ){
	    numBtag++;
	  }

	  if( pt > 1000 ) pt = 999.;

	  double my_jet_sf = 1.;
	  BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;

	  if( abs(flavor)==5 )      jf = BTagEntry::FLAV_B;
	  else if( abs(flavor)==4 ) jf = BTagEntry::FLAV_C;
	  else                      jf = BTagEntry::FLAV_UDSG;

	  bool isBFlav = false;
	  bool isCFlav = false;
	  bool isLFlav = false;
	  if( abs(flavor)==5 )      isBFlav = true;
	  else if( abs(flavor)==4 ) isCFlav = true;
	  else                      isLFlav = true;

	  if( abs(flavor)==5 )      numJet_b++;
	  else if( abs(flavor)==4 ) numJet_c++;
	  else                      numJet_l++;


	  if( sys_name.Contains("JESUp") )        my_jet_sf = reader_JESUp.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("JESDown") ) my_jet_sf = reader_JESDown.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("LFUp") && isBFlav )    my_jet_sf = reader_LFUp.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("LFDown") && isBFlav )  my_jet_sf = reader_LFDown.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("HFUp") && isLFlav )    my_jet_sf = reader_HFUp.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("HFDown") && isLFlav )  my_jet_sf = reader_HFDown.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVHFStats1Up") && isBFlav )   my_jet_sf = reader_HFStats1Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVHFStats1Down") && isBFlav ) my_jet_sf = reader_HFStats1Down.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVHFStats2Up") && isBFlav )   my_jet_sf = reader_HFStats2Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVHFStats2Down") && isBFlav ) my_jet_sf = reader_HFStats2Down.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVLFStats1Up") && isLFlav )   my_jet_sf = reader_LFStats1Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVLFStats1Down") && isLFlav ) my_jet_sf = reader_LFStats1Down.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVLFStats2Up") && isLFlav )   my_jet_sf = reader_LFStats2Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVLFStats2Down") && isLFlav ) my_jet_sf = reader_LFStats2Down.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVCFErr1Up") && isCFlav )   my_jet_sf = reader_CFErr1Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVCFErr1Down") && isCFlav ) my_jet_sf = reader_CFErr1Down.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVCFErr2Up") && isCFlav )   my_jet_sf = reader_CFErr2Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVCFErr2Down") && isCFlav ) my_jet_sf = reader_CFErr2Down.eval(jf, eta, pt, csv);
	  else my_jet_sf = reader.eval(jf, eta, pt, csv);

	  assert (my_jet_sf > 0.);
	  wgt_csv2 *= my_jet_sf;
	}


    if( verbose_ ) std::cout << " ===> test 3.1.6 " << std::endl;


	// Calculate CSV weight
	double wgt_csv_hf, wgt_csv_lf, wgt_csv_cf;
	double wgt_csv = ( insample<0 ) ? 1 : get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);

	wgtCSV = wgt * wgt_csv;


	bool use_verbose = false;
	if( numJet>=1 && numBtag>=1 ) use_verbose = !true;
	double wgt_csv3 = ( insample<0 ) ? 1 : get_btv_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, sys_name, true, use_verbose);
	double wgt_csv4 = ( insample<0 ) ? 1 : get_btv_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, sys_name, false, false);

	if( insample<0 ) wgt_csv2 = 1;
	if( insample<0 ) wgt_csv3 = 1;
	if( insample<0 ) wgt_csv4 = 1;

	double wgtCSV2 = wgt * wgt_csv2;

	double wgtCSV3 = wgt * wgt_csv3;

	double wgtCSV4 = wgt * wgt_csv4;
    if( verbose_ ) std::cout << " ===> test 3.1.7 " << std::endl;


	int fill_numJet  = std::min( MaxNjet-1, numJet );
	int fill_numBtag = std::min( MaxNbtag-1, numBtag );

	h_numJet_em[iSys]->Fill(fill_numJet,wgt);
	h_numJet_wgtCSV_em[iSys]->Fill(fill_numJet,wgtCSV);
	h_numJet_wgtCSV2_em[iSys]->Fill(fill_numJet,wgtCSV2);
	h_numJet_wgtCSV3_em[iSys]->Fill(fill_numJet,wgtCSV3);
	h_numJet_wgtCSV4_em[iSys]->Fill(fill_numJet,wgtCSV4);
    if( verbose_ ) std::cout << " ===> test 3.1.8 " << std::endl;

	if( numJet>=2 ){
	  h_numBtag_em_2j[iSys]->Fill(fill_numBtag,wgt);
	  h_numBtag_wgtCSV_em_2j[iSys]->Fill(fill_numBtag,wgtCSV);
	  h_numBtag_wgtCSV2_em_2j[iSys]->Fill(fill_numBtag,wgtCSV2);
	  h_numBtag_wgtCSV3_em_2j[iSys]->Fill(fill_numBtag,wgtCSV3);
	  h_numBtag_wgtCSV4_em_2j[iSys]->Fill(fill_numBtag,wgtCSV4);
	}
    if( verbose_ ) std::cout << " ===> test 3.1.9 " << std::endl;

      } // end loop over systematics
    } // end require exactly two leptons
  
    /// end EleMu control region
    if( verbose_ ) std::cout << " ===> test 3.0.4 " << std::endl;

      
    // Require exactly one lepton
    if( !oneLep ) continue;
    if( verbose_ ) std::cout << " ===> test 3.1 " << std::endl;

    h_numEvents_wgtCSV->Fill(2.5,wgt_gen*wgt_csv_noSys);
    h_numEvents_wgtCSV_hf->Fill(2.5,wgt_gen*wgt_csv_hf_noSys);
    h_numEvents_wgtCSV_lf->Fill(2.5,wgt_gen*wgt_csv_lf_noSys);
    h_numEvents_wgtCSV_cf->Fill(2.5,wgt_gen*wgt_csv_cf_noSys);

    h_numEvents_nowgtCSV->Fill(2.5,wgt_gen);


    h_event_wgt_csv_1l->Fill(wgt_gen*wgt_csv_noSys);
    h_event_wgt_csv_1l_hf->Fill(wgt_gen*wgt_csv_hf_noSys);
    h_event_wgt_csv_1l_lf->Fill(wgt_gen*wgt_csv_lf_noSys);
    h_event_wgt_csv_1l_cf->Fill(wgt_gen*wgt_csv_cf_noSys);

    if( njet>=4 ){
      h_numEvents_wgtCSV->Fill(3.5,wgt_gen*wgt_csv_noSys);
      h_numEvents_wgtCSV_hf->Fill(3.5,wgt_gen*wgt_csv_hf_noSys);
      h_numEvents_wgtCSV_lf->Fill(3.5,wgt_gen*wgt_csv_lf_noSys);
      h_numEvents_wgtCSV_cf->Fill(3.5,wgt_gen*wgt_csv_cf_noSys);

      h_numEvents_nowgtCSV->Fill(3.5,wgt_gen);

      h_event_wgt_csv_1l_4j->Fill(wgt_gen*wgt_csv_noSys);
      h_event_wgt_csv_1l_4j_hf->Fill(wgt_gen*wgt_csv_hf_noSys);
      h_event_wgt_csv_1l_4j_lf->Fill(wgt_gen*wgt_csv_lf_noSys);
      h_event_wgt_csv_1l_4j_cf->Fill(wgt_gen*wgt_csv_cf_noSys);
    }

    for( int iJet=0; iJet<int(eve->jet_pt_.size()); iJet++ ){
      double pt  = eve->jet_pt_[iJet];
      double eta = eve->jet_eta_[iJet];
      int flavor = eve->jet_hadronFlavour_[iJet];
      double csv =  eve->jet_csv_[iJet];

      double jetAbsEta = fabs(eta);
      
      if( csv < 0.0 ) csv = -0.05;
      if( csv > 1.0 ) csv = 1.0;

      if( !(pt>20. && fabs(eta)<2.4) ) continue;
      
      int iPt = -1, iEta = -1;
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

      int iSysHF = 0;
      int iSysC  = 0;
      int iSysLF = 0;

      if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, pt = " << pt << ", jetAbsEta = " << jetAbsEta << std::endl;
      if( verbose_ ) std::cout << " ===> 3 test 1.3 jet " << iJet << std::endl;

      if (abs(flavor) == 5 ){
	if( iEta>=1 ) iEta=0;
	int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
	double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);

	h_perjet_wgt_csv_1l_hf[iPt][iEta]->Fill(iCSVWgtHF,wgt_gen);
	h_perjet_wgt_csv_1l[iPt][iEta]->Fill(iCSVWgtHF,wgt_gen);
	if( njet>=4 ){
	  h_perjet_wgt_csv_1l_4j_hf[iPt][iEta]->Fill(iCSVWgtHF,wgt_gen);
	  h_perjet_wgt_csv_1l_4j[iPt][iEta]->Fill(iCSVWgtHF,wgt_gen);
	}

	h_perjet_csv_1l_hf[iPt][iEta]->Fill(csv,wgt_gen);
	h_perjet_csv_1l[iPt][iEta]->Fill(csv,wgt_gen);
	if( njet>=4 ){
	  h_perjet_csv_1l_4j_hf[iPt][iEta]->Fill(csv,wgt_gen);
	  h_perjet_csv_1l_4j[iPt][iEta]->Fill(csv,wgt_gen);
	}

	h_perjet_csv_PU_1l_hf[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	h_perjet_csv_PU_1l[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	if( njet>=4 ){
	  h_perjet_csv_PU_1l_4j_hf[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	  h_perjet_csv_PU_1l_4j[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	}
      }
      else if( abs(flavor) == 4 ){
	if( iEta>=1 ) iEta=0;
	int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
	double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);

	h_perjet_wgt_csv_1l_cf[iPt][iEta]->Fill(iCSVWgtC,wgt_gen);
	h_perjet_wgt_csv_1l[iPt][iEta]->Fill(iCSVWgtC,wgt_gen);
	if( njet>=4 ){
	  h_perjet_wgt_csv_1l_4j_cf[iPt][iEta]->Fill(iCSVWgtC,wgt_gen);
	  h_perjet_wgt_csv_1l_4j[iPt][iEta]->Fill(iCSVWgtC,wgt_gen);
	}
	
	h_perjet_csv_1l_cf[iPt][iEta]->Fill(csv,wgt_gen);
	h_perjet_csv_1l[iPt][iEta]->Fill(csv,wgt_gen);
	if( njet>=4 ){
	  h_perjet_csv_1l_4j_cf[iPt][iEta]->Fill(csv,wgt_gen);
	  h_perjet_csv_1l_4j[iPt][iEta]->Fill(csv,wgt_gen);
	}

      	h_perjet_csv_PU_1l_cf[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	h_perjet_csv_PU_1l[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	if( njet>=4 ){
	  h_perjet_csv_PU_1l_4j_cf[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	  h_perjet_csv_PU_1l_4j[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	}
      }
      else {
	if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
	int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
	double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);

	h_perjet_wgt_csv_1l_lf[iPt][iEta]->Fill(iCSVWgtLF,wgt_gen);
	h_perjet_wgt_csv_1l[iPt][iEta]->Fill(iCSVWgtLF,wgt_gen);
	if( njet>=4 ){
	  h_perjet_wgt_csv_1l_4j_lf[iPt][iEta]->Fill(iCSVWgtLF,wgt_gen);
	  h_perjet_wgt_csv_1l_4j[iPt][iEta]->Fill(iCSVWgtLF,wgt_gen);
	}

	h_perjet_csv_1l_lf[iPt][iEta]->Fill(csv,wgt_gen);
	h_perjet_csv_1l[iPt][iEta]->Fill(csv,wgt_gen);
	if( njet>=4 ){
	  h_perjet_csv_1l_4j_lf[iPt][iEta]->Fill(csv,wgt_gen);
	  h_perjet_csv_1l_4j[iPt][iEta]->Fill(csv,wgt_gen);
	}

	h_perjet_csv_PU_1l_lf[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	h_perjet_csv_PU_1l[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	if( njet>=4 ){
	  h_perjet_csv_PU_1l_4j_lf[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	  h_perjet_csv_PU_1l_4j[iPt][iEta]->Fill(csv,wgt_gen*wgt_pu);
	}
      }
    }

    int lepInd = ( numEle>0 ) ? ind_ele[0] : ind_mu[0];

    TLorentzVector myLep;
    myLep.SetPtEtaPhiE( eve->lepton_pt_[lepInd], eve->lepton_eta_[lepInd], eve->lepton_phi_[lepInd], eve->lepton_energy_[lepInd] );


    /////
    /////  Loop over systematics
    /////

    for( int iSys=0; iSys<NumSysCat; iSys++ ){
      if( verbose_ ) std::cout << " ===> test 3.3.17 " << std::endl;

      TString sys_name = sys_cat_labels[iSys];

      vdouble use_jet_pt = eve->jet_pt_;
      vdouble use_jet_eta = eve->jet_eta_;
      vdouble use_jet_phi = eve->jet_phi_;
      vdouble use_jet_energy = eve->jet_energy_;

      vdouble use_jet_csv = eve->jet_csv_;
      vdouble use_jet_cmva = eve->jet_cmva_;
      vdouble use_jet_pileupJetId_fullDiscriminant = eve->jet_pileupJetId_fullDiscriminant_;

      vint use_jet_hadronFlavour = eve->jet_hadronFlavour_;


      double pfMET = eve->pfMET_pt_;
      double pfMET_phi = eve->pfMET_phi_;

      double pfMETNoHF = eve->pfMETNoHF_pt_;
      double pfMETNoHF_phi = eve->pfMETNoHF_phi_;

      double puppiMET = eve->puppiMET_pt_;
      double puppiMET_phi = eve->puppiMET_phi_;


      if( sys_name.Contains("JESUp") ){
	use_jet_pt = eve->jet_JESup_pt_;
	use_jet_eta = eve->jet_JESup_eta_;
	use_jet_phi = eve->jet_JESup_phi_;
	use_jet_energy = eve->jet_JESup_energy_;

	use_jet_csv = eve->jet_JESup_csv_;
	use_jet_cmva = eve->jet_JESup_cmva_;
	use_jet_pileupJetId_fullDiscriminant = eve->jet_JESup_pileupJetId_fullDiscriminant_;

	use_jet_hadronFlavour = eve->jet_JESup_hadronFlavour_;


	pfMET = eve->pfMET_pt_JESup_;
	pfMET_phi = eve->pfMET_phi_JESup_;

	pfMETNoHF = eve->pfMETNoHF_pt_JESup_;
	pfMETNoHF_phi = eve->pfMETNoHF_phi_JESup_;

	puppiMET = eve->puppiMET_pt_JESup_;
	puppiMET_phi = eve->puppiMET_phi_JESup_;
      }
      else if( sys_name.Contains("JESDown") ){
	use_jet_pt = eve->jet_JESdown_pt_;
	use_jet_eta = eve->jet_JESdown_eta_;
	use_jet_phi = eve->jet_JESdown_phi_;
	use_jet_energy = eve->jet_JESdown_energy_;

	use_jet_csv = eve->jet_JESdown_csv_;
	use_jet_cmva = eve->jet_JESdown_cmva_;
	use_jet_pileupJetId_fullDiscriminant = eve->jet_JESdown_pileupJetId_fullDiscriminant_;

	use_jet_hadronFlavour = eve->jet_JESdown_hadronFlavour_;

	pfMET = eve->pfMET_pt_JESdown_;
	pfMET_phi = eve->pfMET_phi_JESdown_;

	pfMETNoHF = eve->pfMETNoHF_pt_JESdown_;
	pfMETNoHF_phi = eve->pfMETNoHF_phi_JESdown_;

	puppiMET = eve->puppiMET_pt_JESdown_;
	puppiMET_phi = eve->puppiMET_phi_JESdown_;
      }
      else if( sys_name.Contains("JERUp") ){
	use_jet_pt = eve->jet_JERup_pt_;
	use_jet_eta = eve->jet_JERup_eta_;
	use_jet_phi = eve->jet_JERup_phi_;
	use_jet_energy = eve->jet_JERup_energy_;

	use_jet_csv = eve->jet_JERup_csv_;
	use_jet_cmva = eve->jet_JERup_cmva_;
	use_jet_pileupJetId_fullDiscriminant = eve->jet_JERup_pileupJetId_fullDiscriminant_;

	use_jet_hadronFlavour = eve->jet_JERup_hadronFlavour_;

	pfMET = eve->pfMET_pt_JERup_;
	pfMET_phi = eve->pfMET_phi_JERup_;

	pfMETNoHF = eve->pfMETNoHF_pt_JERup_;
	pfMETNoHF_phi = eve->pfMETNoHF_phi_JERup_;

	puppiMET = eve->puppiMET_pt_JERup_;
	puppiMET_phi = eve->puppiMET_phi_JERup_;
      }
      else if( sys_name.Contains("JERDown") ){
	use_jet_pt = eve->jet_JERdown_pt_;
	use_jet_eta = eve->jet_JERdown_eta_;
	use_jet_phi = eve->jet_JERdown_phi_;
	use_jet_energy = eve->jet_JERdown_energy_;

	use_jet_csv = eve->jet_JERdown_csv_;
	use_jet_cmva = eve->jet_JERdown_cmva_;
	use_jet_pileupJetId_fullDiscriminant = eve->jet_JERdown_pileupJetId_fullDiscriminant_;

	use_jet_hadronFlavour = eve->jet_JERdown_hadronFlavour_;

	pfMET = eve->pfMET_pt_JERdown_;
	pfMET_phi = eve->pfMET_phi_JERdown_;

	pfMETNoHF = eve->pfMETNoHF_pt_JERdown_;
	pfMETNoHF_phi = eve->pfMETNoHF_phi_JERdown_;

	puppiMET = eve->puppiMET_pt_JERdown_;
	puppiMET_phi = eve->puppiMET_phi_JERdown_;
      }

     if( verbose_ ) std::cout << " ===> test 3.3.18 " << std::endl;

      // Set lepton scale factor equal to 1 for now
      double wgt_lepIdSF = 1.;
      if( sys_name.Contains("lepIdSFUp") )        wgt_lepIdSF = 1.;
      else if( sys_name.Contains("lepIdSFDown") ) wgt_lepIdSF = 1.;


      
      if( insample>=0 ){
	double use_id_pt = std::min(myLep.Pt(), 119.8);
	double use_id_eta = myLep.Eta();
	double use_id_abseta = fabs(use_id_eta);

	double id_sf = 1., iso_sf = 1., hlt_sf = 1.;
	if( oneEle ){
	  double use_id_scEta = eve->lepton_eta_[lepInd];

	  id_sf = h_ele_id_sf->GetBinContent(h_ele_id_sf->FindBin(use_id_scEta,use_id_pt));
	  iso_sf = h_ele_iso_sf->GetBinContent(h_ele_iso_sf->FindBin(use_id_eta,use_id_pt));
	  hlt_sf = h_ele_trig_sf->GetBinContent(h_ele_trig_sf->FindBin(use_id_pt,use_id_scEta));
	}
	else if( oneMu ){
	  id_sf = h_mu_id_sf->GetBinContent(h_mu_id_sf->FindBin(use_id_abseta,use_id_pt));
	  iso_sf = h_mu_iso_sf->GetBinContent(h_mu_iso_sf->FindBin(use_id_abseta,use_id_pt));
	  hlt_sf = h_mu_trig_sf->GetBinContent(h_mu_trig_sf->FindBin(use_id_abseta,use_id_pt));
	}

	wgt_lepIdSF = id_sf * iso_sf * hlt_sf;
	double hlt_err = ( oneEle ) ? 0.02 : 0.01;
	if( sys_name.Contains("lepIdSFUp") )        wgt_lepIdSF = id_sf*(1 + 0.01) * iso_sf*(1 + 0.01) * hlt_sf*(1 + hlt_err);
	else if( sys_name.Contains("lepIdSFDown") ) wgt_lepIdSF = id_sf*(1 - 0.01) * iso_sf*(1 - 0.01) * hlt_sf*(1 - hlt_err);

	// if( iSys==0 ){
	//   if( oneEle ) printf("Electron: pT = %.1f, eta = %+.2f, id_sf = %.3f, iso_sf = %.3f, hlt_sf = %.3f, wgt_lepIdSF = %.3f \n", use_id_pt, use_id_eta, id_sf, iso_sf, hlt_sf, wgt_lepIdSF);
	//   else if( oneMu ) printf("Muon: pT = %.1f, eta = %+.2f, id_sf = %.3f, iso_sf = %.3f, hlt_sf = %.3f, wgt_lepIdSF = %.3f \n", use_id_pt, use_id_eta, id_sf, iso_sf, hlt_sf, wgt_lepIdSF);
	// }
      }
 
      // PU systematic
      double use_wgt_pu = wgt_pu;
      if( sys_name.Contains("PUUp") )        use_wgt_pu = wgt_pu_up;
      else if( sys_name.Contains("PUDown") ) use_wgt_pu = wgt_pu_down;


     if( verbose_ ) std::cout << " ===> test 3.3.19 " << std::endl;

      // Factorization and normalization scale uncertainty
      double wgt_LHEscale = 1.;
      if( insample>-1 ){
	double originalXWGTUP = eve->originalXWGTUP_;
	vdouble lhe_event_weights = eve->LHEEvent_weights_;

	if( lhe_event_weights.size()>76 && originalXWGTUP!=0 ){
	  int use_weight_index = 0;

	  if( sys_name.Contains("muFUp") )        use_weight_index = 1;
	  else if( sys_name.Contains("muFDown") ) use_weight_index = 2;
	  else if( sys_name.Contains("muRUp") )   use_weight_index = 3;
	  else if( sys_name.Contains("muRDown") ) use_weight_index = 6;
	  else if( sys_name.Contains("muRmuFUp") )   use_weight_index = 4;
	  else if( sys_name.Contains("muRmuFDown") ) use_weight_index = 8;
	  else if( sys_name.Contains("muRmuFDown") ) use_weight_index = 8;
	  else if( sys_name.Contains("NNPDFUp") ) use_weight_index = 75;
	  else if( sys_name.Contains("NNPDFDown") ) use_weight_index = 13;

	  wgt_LHEscale = lhe_event_weights[use_weight_index] / originalXWGTUP;
	}
      }
     if( verbose_ ) std::cout << " ===> test 3.3.20 " << std::endl;


      //wgt = wgt_gen * wgt_lumi * wgt_pu * wgt_lepIdSF * wgt_LHEscale;
      wgt = wgt_gen * wgt_lumi * use_wgt_pu * wgt_lepIdSF * wgt_LHEscale;


      vdouble jetPts;
      vdouble jetEtas;
      vdouble jetCSVs;
      vdouble jetCMVAs;
      vint    jetFlavors;

      vdouble temp_test_jetPts;
      vdouble temp_test_jetEtas;
      vdouble temp_test_jetCSVs;
      vint    temp_test_jetFlavors;

      vint ind_jet;

      int numJet = 0;
      int numBtag = 0;
      int numBtag_CMVA = 0;

      double HT30=0;
      TLorentzVector sumJet;
      bool firstJet = false;

      double wgt_csv2=1.;
      double wgt_csv2_hf=1., wgt_csv2_lf=1., wgt_csv2_cf=1.;

      double wgt_cmva2=1.;
	    
      if( verbose_ ) std::cout << " ===> test 3.3.21 " << std::endl;

      for( int iJet=0; iJet<int(use_jet_pt.size()); iJet++ ){
	TLorentzVector myJet;
	myJet.SetPtEtaPhiE( use_jet_pt[iJet], use_jet_eta[iJet], use_jet_phi[iJet], use_jet_energy[iJet] );

	if( !firstJet ) sumJet = myJet;
	else            sumJet += myJet;

	double pt  = use_jet_pt[iJet];
	double eta = use_jet_eta[iJet];

	if( !(pt>30. && fabs(eta)<2.4) ) continue;

	double dR = myJet.DeltaR(myLep);
	if( iSys==0 ) h_deltaR_jet_lep->Fill(dR,wgt);

	if( dR < 0.4 ) continue;

	double csv = use_jet_csv[iJet];
	if( csv < 0.0 ) csv = -0.05;
	if( csv > 1.0 ) csv = 1.0;

	double cmva = use_jet_cmva[iJet];

	int flavor = use_jet_hadronFlavour[iJet];

	ind_jet.push_back(iJet);

	HT30 += pt;

	temp_test_jetPts.clear();
	temp_test_jetEtas.clear();
	temp_test_jetCSVs.clear();
	temp_test_jetFlavors.clear();

	temp_test_jetPts.push_back(pt);
	temp_test_jetEtas.push_back(eta);
	temp_test_jetCSVs.push_back(csv);
	temp_test_jetFlavors.push_back(flavor);


	jetPts.push_back(pt);
	jetEtas.push_back(eta);
	jetCSVs.push_back(csv);
	jetCMVAs.push_back(cmva);
	jetFlavors.push_back(flavor);

	numJet++;
	if( csv > 0.800 ){
	  numBtag++;
	}
	if( cmva > 0.185 ){
	  numBtag_CMVA++;
	}

	if( pt > 1000 ) pt = 999.;

	double my_jet_sf = 1.;
	BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;

	if( abs(flavor)==5 )      jf = BTagEntry::FLAV_B;
	else if( abs(flavor)==4 ) jf = BTagEntry::FLAV_C;
	else                      jf = BTagEntry::FLAV_UDSG;

	// printf(" iSys = %d, iJet = %d, pt = %.1f, eta = %.2f, csv = %.2f, flavor = %d \n",
	//        iSys, iJet, pt, eta, csv, flavor);

	//my_jet_sf = reader.eval(jf, eta, pt, csv);

	bool isBFlav = false;
	bool isCFlav = false;
	bool isLFlav = false;
	if( abs(flavor)==5 )      isBFlav = true;
	else if( abs(flavor)==4 ) isCFlav = true;
	else                      isLFlav = true;
	if( verbose_ ) std::cout << " ===> test 3.3.22 sys_name = " << sys_name << ", flav = " << flavor << ", jet" << iJet << " pt = " << pt << ", eta = " << eta << ", csv = " << csv << std::endl;

	if( sys_name.Contains("JESUp") )        my_jet_sf = reader_JESUp.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("JESDown") ) my_jet_sf = reader_JESDown.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("LFUp") && isBFlav )    my_jet_sf = reader_LFUp.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("LFDown") && isBFlav )  my_jet_sf = reader_LFDown.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("HFUp") && isLFlav )    my_jet_sf = reader_HFUp.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("HFDown") && isLFlav )  my_jet_sf = reader_HFDown.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("CSVHFStats1Up") && isBFlav )   my_jet_sf = reader_HFStats1Up.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("CSVHFStats1Down") && isBFlav ) my_jet_sf = reader_HFStats1Down.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("CSVHFStats2Up") && isBFlav )   my_jet_sf = reader_HFStats2Up.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("CSVHFStats2Down") && isBFlav ) my_jet_sf = reader_HFStats2Down.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("CSVLFStats1Up") && isLFlav )   my_jet_sf = reader_LFStats1Up.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("CSVLFStats1Down") && isLFlav ) my_jet_sf = reader_LFStats1Down.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("CSVLFStats2Up") && isLFlav )   my_jet_sf = reader_LFStats2Up.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("CSVLFStats2Down") && isLFlav ) my_jet_sf = reader_LFStats2Down.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("CSVCFErr1Up") && isCFlav )   my_jet_sf = reader_CFErr1Up.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("CSVCFErr1Down") && isCFlav ) my_jet_sf = reader_CFErr1Down.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("CSVCFErr2Up") && isCFlav )   my_jet_sf = reader_CFErr2Up.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("CSVCFErr2Down") && isCFlav ) my_jet_sf = reader_CFErr2Down.eval(jf, eta, pt, csv);
	else my_jet_sf = reader.eval(jf, eta, pt, csv);


	double my_jet_sf_cmva = 1.;
	if( sys_name.Contains("JESUp") )        my_jet_sf_cmva = reader_cmvav2_JESUp.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("JESDown") ) my_jet_sf_cmva = reader_cmvav2_JESDown.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("LFUp") && isBFlav )    my_jet_sf_cmva = reader_cmvav2_LFUp.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("LFDown") && isBFlav )  my_jet_sf_cmva = reader_cmvav2_LFDown.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("HFUp") && isLFlav )    my_jet_sf_cmva = reader_cmvav2_HFUp.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("HFDown") && isLFlav )  my_jet_sf_cmva = reader_cmvav2_HFDown.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("CSVHFStats1Up") && isBFlav )   my_jet_sf_cmva = reader_cmvav2_HFStats1Up.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("CSVHFStats1Down") && isBFlav ) my_jet_sf_cmva = reader_cmvav2_HFStats1Down.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("CSVHFStats2Up") && isBFlav )   my_jet_sf_cmva = reader_cmvav2_HFStats2Up.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("CSVHFStats2Down") && isBFlav ) my_jet_sf_cmva = reader_cmvav2_HFStats2Down.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("CSVLFStats1Up") && isLFlav )   my_jet_sf_cmva = reader_cmvav2_LFStats1Up.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("CSVLFStats1Down") && isLFlav ) my_jet_sf_cmva = reader_cmvav2_LFStats1Down.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("CSVLFStats2Up") && isLFlav )   my_jet_sf_cmva = reader_cmvav2_LFStats2Up.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("CSVLFStats2Down") && isLFlav ) my_jet_sf_cmva = reader_cmvav2_LFStats2Down.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("CSVCFErr1Up") && isCFlav )   my_jet_sf_cmva = reader_cmvav2_CFErr1Up.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("CSVCFErr1Down") && isCFlav ) my_jet_sf_cmva = reader_cmvav2_CFErr1Down.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("CSVCFErr2Up") && isCFlav )   my_jet_sf_cmva = reader_cmvav2_CFErr2Up.eval(jf, eta, pt, cmva);
	else if( sys_name.Contains("CSVCFErr2Down") && isCFlav ) my_jet_sf_cmva = reader_cmvav2_CFErr2Down.eval(jf, eta, pt, cmva);
	else my_jet_sf_cmva = reader_cmvav2.eval(jf, eta, pt, cmva);

	
	assert (my_jet_sf_cmva > 0.);
	wgt_cmva2 *= my_jet_sf_cmva;

	if( cnt<100 && iSys==0 && false ){
	  printf("  iJet = %d: pt = %4.1f, eta = %+.2f, flavor = %+d, CMVA = %.3f: weight csv file = %.3f \n",
		 iJet, pt, eta, flavor, cmva, my_jet_sf_cmva);
	}
	if( verbose_ ) std::cout << " ===> test 3.3.23 " << std::endl;

	double my_jet_sf_hf = 1.;
	double my_jet_sf_lf = 1.;
	double my_jet_sf_cf = 1.;
	if( abs(flavor)==5 ){
	  if( sys_name.Contains("JESUp") )        my_jet_sf_hf = reader_JESUp.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("JESDown") ) my_jet_sf_hf = reader_JESDown.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("LFUp") )    my_jet_sf_hf = reader_LFUp.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("LFDown") )  my_jet_sf_hf = reader_LFDown.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVHFStats1Up") )   my_jet_sf_hf = reader_HFStats1Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVHFStats1Down") ) my_jet_sf_hf = reader_HFStats1Down.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVHFStats2Up") )   my_jet_sf_hf = reader_HFStats2Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVHFStats2Down") ) my_jet_sf_hf = reader_HFStats2Down.eval(jf, eta, pt, csv);
	  else my_jet_sf_hf = reader.eval(jf, eta, pt, csv);
	}
	else if( abs(flavor)==4 ){
	  if( sys_name.Contains("JESUp") )        my_jet_sf_cf = reader_JESUp.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("JESDown") ) my_jet_sf_cf = reader_JESDown.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVCFErr1Up") )   my_jet_sf_cf = reader_CFErr1Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVCFErr1Down") ) my_jet_sf_cf = reader_CFErr1Down.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVCFErr2Up") )   my_jet_sf_cf = reader_CFErr2Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVCFErr2Down") ) my_jet_sf_cf = reader_CFErr2Down.eval(jf, eta, pt, csv);
	  else my_jet_sf_cf = reader.eval(jf, eta, pt, csv);
	}
	else {
	  if( sys_name.Contains("JESUp") )        my_jet_sf_lf = reader_JESUp.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("JESDown") ) my_jet_sf_lf = reader_JESDown.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("HFUp") )    my_jet_sf_lf = reader_HFUp.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("HFDown") )  my_jet_sf_lf = reader_HFDown.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVLFStats1Up") )   my_jet_sf_lf = reader_LFStats1Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVLFStats1Down") ) my_jet_sf_lf = reader_LFStats1Down.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVLFStats2Up") )   my_jet_sf_lf = reader_LFStats2Up.eval(jf, eta, pt, csv);
	  else if( sys_name.Contains("CSVLFStats2Down") ) my_jet_sf_lf = reader_LFStats2Down.eval(jf, eta, pt, csv);
	  else my_jet_sf_lf = reader.eval(jf, eta, pt, csv);
	}


	if( verbose_ ) std::cout << " ===> test 3.3.24 " << std::endl;

	assert (my_jet_sf > 0.);
	wgt_csv2 *= my_jet_sf;

	wgt_csv2_hf *= my_jet_sf_hf;
	wgt_csv2_cf *= my_jet_sf_cf;
	wgt_csv2_lf *= my_jet_sf_lf;


	
	double temp_test_wgt_csv_hf=1, temp_test_wgt_csv_lf=1, temp_test_wgt_csv_cf=1;
	double temp_test_wgt_csv = get_csv_wgt(temp_test_jetPts, temp_test_jetEtas, temp_test_jetCSVs, temp_test_jetFlavors, iSys, 
					       temp_test_wgt_csv_hf, temp_test_wgt_csv_lf, temp_test_wgt_csv_cf);

	if( iSys==0 ){
	  if( abs(flavor)==5 ){
	    h_b_csv_wgt->Fill(csv,temp_test_wgt_csv);
	    h_b_csv_wgt2->Fill(csv,my_jet_sf);
	    h_b_csv_diff_wgt2_wgt->Fill(csv,my_jet_sf-temp_test_wgt_csv);
	  }
	  else if( abs(flavor)==4 ){
	    h_c_csv_wgt->Fill(csv,temp_test_wgt_csv);
	    h_c_csv_wgt2->Fill(csv,my_jet_sf);
	    h_c_csv_diff_wgt2_wgt->Fill(csv,my_jet_sf-temp_test_wgt_csv);
	  }
	  else {
	    h_l_csv_wgt->Fill(csv,temp_test_wgt_csv);
	    h_l_csv_wgt2->Fill(csv,my_jet_sf);
	    h_l_csv_diff_wgt2_wgt->Fill(csv,my_jet_sf-temp_test_wgt_csv);
	  }
	}

	if( sys_name.Contains("JESUp") && false ){
	  printf("  iJet = %d: pt = %4.1f, eta = %+.2f, flavor = %+d, CSVv2 = %.2f: weight csv file = %.3f, weight hist file = %.3f \n",
		 iJet, pt, eta, flavor, csv, my_jet_sf, temp_test_wgt_csv);
	}
      }


      // Calculate CSV weight
      double wgt_csv_hf, wgt_csv_lf, wgt_csv_cf;
      double wgt_csv = ( insample<0 ) ? 1 : get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);

      wgtCSV = wgt * wgt_csv;
    if( verbose_ ) std::cout << " ===> test 3.3.1 " << std::endl;

      double wgt_csv3 = ( insample<0 ) ? 1 : get_btv_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, sys_name, true, false);
    if( verbose_ ) std::cout << " ===> test 3.3.2 " << std::endl;
      double wgt_csv4 = ( insample<0 ) ? 1 : get_btv_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, sys_name, false, false);
    if( verbose_ ) std::cout << " ===> test 3.3.3 " << std::endl;

      if( insample<0 ) wgt_csv2 = 1;
      if( insample<0 ) wgt_csv3 = 1;
      if( insample<0 ) wgt_csv4 = 1;

      double wgtCSV2 = wgt * wgt_csv2;

      double wgtCSV3 = wgt * wgt_csv3;

      double wgtCSV4 = wgt * wgt_csv4;

      // Calculate CMVA weight
      double wgt_cmva_hf, wgt_cmva_lf, wgt_cmva_cf;
      double wgt_cmva = ( insample<0 ) ? 1 : get_cmva_wgt(jetPts, jetEtas, jetCMVAs, jetFlavors, iSys, wgt_cmva_hf, wgt_cmva_lf, wgt_cmva_cf);
    if( verbose_ ) std::cout << " ===> test 3.3.4 " << std::endl;

      double wgtCMVA = wgt * wgt_cmva;

      if( insample<0 ) wgt_cmva2 = 1;
      double wgtCMVA2 = wgt * wgt_cmva2;

      
      // Calculate CSV weight
      double wgt_csv_hf_v5, wgt_csv_lf_v5, wgt_csv_cf_v5;
      double wgt_csv5 = ( insample<0 ) ? 1 : get_csv_wgt_old(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf_v5, wgt_csv_lf_v5, wgt_csv_cf_v5);
    if( verbose_ ) std::cout << " ===> test 3.3.5 " << std::endl;

      double wgtCSV5 = wgt * wgt_csv5;

      //if( iSys==0 ) printf("\t total event weight: csv file = %.3f, hist file = %.3f, btv csv file = %.3f \n", wgt_csv2, wgt_csv, wgt_csv3);

      //std::cout << "\t main: iSys = " << iSys << ",\t get_csv_wgt wgt_csv_lf = " << wgt_csv_lf << std::endl;

      h_diff_wgtCSV_wgtCSV2_perSys->Fill(iSys,wgt_csv2 - wgt_csv);
      h_hf_diff_wgtCSV_wgtCSV2_perSys->Fill(iSys,wgt_csv2_hf - wgt_csv_hf);
      h_lf_diff_wgtCSV_wgtCSV2_perSys->Fill(iSys,wgt_csv2_lf - wgt_csv_lf);
      h_cf_diff_wgtCSV_wgtCSV2_perSys->Fill(iSys,wgt_csv2_cf - wgt_csv_cf);

      h_hf_wgtCSV_perSys->Fill(iSys,wgt_csv_hf);
      h_lf_wgtCSV_perSys->Fill(iSys,wgt_csv_lf);
      h_cf_wgtCSV_perSys->Fill(iSys,wgt_csv_cf);

      h_hf_wgtCSV2_perSys->Fill(iSys,wgt_csv2_hf);
      h_lf_wgtCSV2_perSys->Fill(iSys,wgt_csv2_lf);
      h_cf_wgtCSV2_perSys->Fill(iSys,wgt_csv2_cf);


      h_numEvents_perSys->Fill(iSys,wgt);
      h_numEvents_perSys_wgtCSV->Fill(iSys,wgtCSV);
      if( iSys==0 ) h_numEvents_perSys_wgtCSV->Fill(-1,wgt);

      h_numEvents_perSys_wgtCSV2->Fill(iSys,wgtCSV2);
      if( iSys==0 ) h_numEvents_perSys_wgtCSV2->Fill(-1,wgt);

      h_numEvents_perSys_wgtCMVA->Fill(iSys,wgtCMVA);
      if( iSys==0 ) h_numEvents_perSys_wgtCMVA->Fill(-1,wgt);

      h_numEvents_perSys_wgtCMVA2->Fill(iSys,wgtCMVA2);
      if( iSys==0 ) h_numEvents_perSys_wgtCMVA2->Fill(-1,wgt);


      int fill_numJet  = std::min( MaxNjet-1, numJet );
      int fill_numBtag = std::min( MaxNbtag-1, numBtag );
      int fill_numBtag_CMVA = std::min( MaxNbtag-1, numBtag_CMVA );

      h_numJet[iSys]->Fill(fill_numJet,wgt);
      h_numJet_wgtCSV[iSys]->Fill(fill_numJet,wgtCSV);
      h_numJet_wgtCSV2[iSys]->Fill(fill_numJet,wgtCSV2);
      h_numJet_wgtCSV3[iSys]->Fill(fill_numJet,wgtCSV3);
      h_numJet_wgtCSV4[iSys]->Fill(fill_numJet,wgtCSV4);

      h_numBtag[iSys]->Fill(fill_numBtag,wgt);
      h_numBtag_wgtCSV[iSys]->Fill(fill_numBtag,wgtCSV);
      h_numBtag_wgtCSV2[iSys]->Fill(fill_numBtag,wgtCSV2);
      h_numBtag_wgtCSV3[iSys]->Fill(fill_numBtag,wgtCSV3);
      h_numBtag_wgtCSV4[iSys]->Fill(fill_numBtag,wgtCSV4);




      // Require at least four jets
      if( !(numJet>=4) ) continue;


      if( oneLep ) h_event_selection[iSys]->Fill(0.5+3, wgt); // all
      if( oneMu )  h_mu_event_selection[iSys]->Fill(0.5+3, wgt); // all
      if( oneEle ) h_ele_event_selection[iSys]->Fill(0.5+3, wgt); // all

      if( oneLep ) h_event_selection_wgtCSV[iSys]->Fill(0.5+3, wgtCSV); // all
      if( oneMu )  h_mu_event_selection_wgtCSV[iSys]->Fill(0.5+3, wgtCSV); // all
      if( oneEle ) h_ele_event_selection_wgtCSV[iSys]->Fill(0.5+3, wgtCSV); // all


      h_numBtag_4j[iSys]->Fill(fill_numBtag,wgt);
      h_numBtag_wgtCSV_4j[iSys]->Fill(fill_numBtag,wgtCSV);
      h_numBtag_wgtCSV2_4j[iSys]->Fill(fill_numBtag,wgtCSV2);
      h_numBtag_wgtCSV3_4j[iSys]->Fill(fill_numBtag,wgtCSV3);
      h_numBtag_wgtCSV4_4j[iSys]->Fill(fill_numBtag,wgtCSV4);
      h_numBtag_wgtCSV5_4j[iSys]->Fill(fill_numBtag,wgtCSV5);

      
      h_numBtag_CMVA_4j[iSys]->Fill(fill_numBtag_CMVA,wgt);
      h_numBtag_CMVA_wgtCSV_4j[iSys]->Fill(fill_numBtag_CMVA,wgtCSV);
      h_numBtag_CMVA_wgtCMVA_4j[iSys]->Fill(fill_numBtag_CMVA,wgtCMVA);

      
      if( oneEle ){
	h_numBtag_1e_4j[iSys]->Fill(fill_numBtag,wgt);
	h_numBtag_wgtCSV_1e_4j[iSys]->Fill(fill_numBtag,wgtCSV);
	h_numBtag_wgtCSV2_1e_4j[iSys]->Fill(fill_numBtag,wgtCSV2);
	h_numBtag_wgtCSV3_1e_4j[iSys]->Fill(fill_numBtag,wgtCSV3);
	h_numBtag_wgtCSV4_1e_4j[iSys]->Fill(fill_numBtag,wgtCSV4);
      }
      else if( oneMu ){
	h_numBtag_1m_4j[iSys]->Fill(fill_numBtag,wgt);
	h_numBtag_wgtCSV_1m_4j[iSys]->Fill(fill_numBtag,wgtCSV);
	h_numBtag_wgtCSV2_1m_4j[iSys]->Fill(fill_numBtag,wgtCSV2);
	h_numBtag_wgtCSV3_1m_4j[iSys]->Fill(fill_numBtag,wgtCSV3);
	h_numBtag_wgtCSV4_1m_4j[iSys]->Fill(fill_numBtag,wgtCSV4);
      }

      if( numJet==4 ){
	h_numBtag_eq4j[iSys]->Fill(fill_numBtag,wgt);
	h_numBtag_wgtCSV_eq4j[iSys]->Fill(fill_numBtag,wgtCSV);
	h_numBtag_wgtCSV2_eq4j[iSys]->Fill(fill_numBtag,wgtCSV2);
	h_numBtag_wgtCSV3_eq4j[iSys]->Fill(fill_numBtag,wgtCSV3);
	h_numBtag_wgtCSV4_eq4j[iSys]->Fill(fill_numBtag,wgtCSV4);

	if( oneEle ){
	  h_numBtag_1e_eq4j[iSys]->Fill(fill_numBtag,wgt);
	  h_numBtag_wgtCSV_1e_eq4j[iSys]->Fill(fill_numBtag,wgtCSV);
	  h_numBtag_wgtCSV2_1e_eq4j[iSys]->Fill(fill_numBtag,wgtCSV2);
	  h_numBtag_wgtCSV3_1e_eq4j[iSys]->Fill(fill_numBtag,wgtCSV3);
	  h_numBtag_wgtCSV4_1e_eq4j[iSys]->Fill(fill_numBtag,wgtCSV4);
	}
	else if( oneMu ){
	  h_numBtag_1m_eq4j[iSys]->Fill(fill_numBtag,wgt);
	  h_numBtag_wgtCSV_1m_eq4j[iSys]->Fill(fill_numBtag,wgtCSV);
	  h_numBtag_wgtCSV2_1m_eq4j[iSys]->Fill(fill_numBtag,wgtCSV2);
	  h_numBtag_wgtCSV3_1m_eq4j[iSys]->Fill(fill_numBtag,wgtCSV3);
	  h_numBtag_wgtCSV4_1m_eq4j[iSys]->Fill(fill_numBtag,wgtCSV4);
	}
      }

    if( verbose_ ) std::cout << " ===> test 3.3.6 " << std::endl;


      double fill_pfMETNoHF = std::min( metmax-0.001, pfMETNoHF);

      // Fill MET histograms
      h_pfMETNoHF_pt_4j_1l[iSys]->Fill(fill_pfMETNoHF,wgt);
      if( oneEle )     h_pfMETNoHF_pt_4j_1e[iSys]->Fill(fill_pfMETNoHF,wgt);
      else if( oneMu ) h_pfMETNoHF_pt_4j_1m[iSys]->Fill(fill_pfMETNoHF,wgt);

      if( pfMETNoHF > 30 ){
	h_numBtag_4j_met30[iSys]->Fill(fill_numBtag,wgt);
	h_numBtag_wgtCSV_4j_met30[iSys]->Fill(fill_numBtag,wgtCSV);
	h_numBtag_wgtCSV2_4j_met30[iSys]->Fill(fill_numBtag,wgtCSV2);
	h_numBtag_wgtCSV3_4j_met30[iSys]->Fill(fill_numBtag,wgtCSV3);
	h_numBtag_wgtCSV4_4j_met30[iSys]->Fill(fill_numBtag,wgtCSV4);

	if( numJet==4 ){
	  h_numBtag_eq4j_met30[iSys]->Fill(fill_numBtag,wgt);
	  h_numBtag_wgtCSV_eq4j_met30[iSys]->Fill(fill_numBtag,wgtCSV);
	  h_numBtag_wgtCSV2_eq4j_met30[iSys]->Fill(fill_numBtag,wgtCSV2);
	  h_numBtag_wgtCSV3_eq4j_met30[iSys]->Fill(fill_numBtag,wgtCSV3);
	  h_numBtag_wgtCSV4_eq4j_met30[iSys]->Fill(fill_numBtag,wgtCSV4);
	}

	if( oneEle ){
	  h_numBtag_1e_4j_met30[iSys]->Fill(fill_numBtag,wgt);
	  h_numBtag_wgtCSV_1e_4j_met30[iSys]->Fill(fill_numBtag,wgtCSV);
	  h_numBtag_wgtCSV2_1e_4j_met30[iSys]->Fill(fill_numBtag,wgtCSV2);
	  h_numBtag_wgtCSV3_1e_4j_met30[iSys]->Fill(fill_numBtag,wgtCSV3);
	  h_numBtag_wgtCSV4_1e_4j_met30[iSys]->Fill(fill_numBtag,wgtCSV4);

	  if( numJet==4 ){
	    h_numBtag_1e_eq4j_met30[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_1e_eq4j_met30[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_1e_eq4j_met30[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_1e_eq4j_met30[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_1e_eq4j_met30[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	}

	if( oneMu ){
	  h_numBtag_1m_4j_met30[iSys]->Fill(fill_numBtag,wgt);
	  h_numBtag_wgtCSV_1m_4j_met30[iSys]->Fill(fill_numBtag,wgtCSV);
	  h_numBtag_wgtCSV2_1m_4j_met30[iSys]->Fill(fill_numBtag,wgtCSV2);
	  h_numBtag_wgtCSV3_1m_4j_met30[iSys]->Fill(fill_numBtag,wgtCSV3);
	  h_numBtag_wgtCSV4_1m_4j_met30[iSys]->Fill(fill_numBtag,wgtCSV4);

	  if( numJet==4 ){
	    h_numBtag_1m_eq4j_met30[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_1m_eq4j_met30[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_1m_eq4j_met30[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_1m_eq4j_met30[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_1m_eq4j_met30[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	}
      }


      if( pfMETNoHF > 50 ){
	h_numBtag_4j_met50[iSys]->Fill(fill_numBtag,wgt);
	h_numBtag_wgtCSV_4j_met50[iSys]->Fill(fill_numBtag,wgtCSV);
	h_numBtag_wgtCSV2_4j_met50[iSys]->Fill(fill_numBtag,wgtCSV2);
	h_numBtag_wgtCSV3_4j_met50[iSys]->Fill(fill_numBtag,wgtCSV3);
	h_numBtag_wgtCSV4_4j_met50[iSys]->Fill(fill_numBtag,wgtCSV4);

	if( numJet==4 ){
	  h_numBtag_eq4j_met50[iSys]->Fill(fill_numBtag,wgt);
	  h_numBtag_wgtCSV_eq4j_met50[iSys]->Fill(fill_numBtag,wgtCSV);
	  h_numBtag_wgtCSV2_eq4j_met50[iSys]->Fill(fill_numBtag,wgtCSV2);
	  h_numBtag_wgtCSV3_eq4j_met50[iSys]->Fill(fill_numBtag,wgtCSV3);
	  h_numBtag_wgtCSV4_eq4j_met50[iSys]->Fill(fill_numBtag,wgtCSV4);
	}

	if( oneEle ){
	  h_numBtag_1e_4j_met50[iSys]->Fill(fill_numBtag,wgt);
	  h_numBtag_wgtCSV_1e_4j_met50[iSys]->Fill(fill_numBtag,wgtCSV);
	  h_numBtag_wgtCSV2_1e_4j_met50[iSys]->Fill(fill_numBtag,wgtCSV2);
	  h_numBtag_wgtCSV3_1e_4j_met50[iSys]->Fill(fill_numBtag,wgtCSV3);
	  h_numBtag_wgtCSV4_1e_4j_met50[iSys]->Fill(fill_numBtag,wgtCSV4);

	  if( numJet==4 ){
	    h_numBtag_1e_eq4j_met50[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_1e_eq4j_met50[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_1e_eq4j_met50[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_1e_eq4j_met50[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_1e_eq4j_met50[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	}

	if( oneMu ){
	  h_numBtag_1m_4j_met50[iSys]->Fill(fill_numBtag,wgt);
	  h_numBtag_wgtCSV_1m_4j_met50[iSys]->Fill(fill_numBtag,wgtCSV);
	  h_numBtag_wgtCSV2_1m_4j_met50[iSys]->Fill(fill_numBtag,wgtCSV2);
	  h_numBtag_wgtCSV3_1m_4j_met50[iSys]->Fill(fill_numBtag,wgtCSV3);
	  h_numBtag_wgtCSV4_1m_4j_met50[iSys]->Fill(fill_numBtag,wgtCSV4);

	  if( numJet==4 ){
	    h_numBtag_1m_eq4j_met50[iSys]->Fill(fill_numBtag,wgt);
	    h_numBtag_wgtCSV_1m_eq4j_met50[iSys]->Fill(fill_numBtag,wgtCSV);
	    h_numBtag_wgtCSV2_1m_eq4j_met50[iSys]->Fill(fill_numBtag,wgtCSV2);
	    h_numBtag_wgtCSV3_1m_eq4j_met50[iSys]->Fill(fill_numBtag,wgtCSV3);
	    h_numBtag_wgtCSV4_1m_eq4j_met50[iSys]->Fill(fill_numBtag,wgtCSV4);
	  }
	}
      }

    if( verbose_ ) std::cout << " ===> test 3.3.7 " << std::endl;


      if( (numBtag_CMVA>=2) ){

	
	h_numEvents_perSys_wgtCMVA_4j2t->Fill(iSys,wgtCMVA);
	if( iSys==0 ) h_numEvents_perSys_wgtCMVA_4j2t->Fill(-1,wgt);

	if( numJet>=6 && numBtag_CMVA>=4 ){
	  h_numEvents_perSys_wgtCMVA_6j4t->Fill(iSys,wgtCMVA);
	  if( iSys==0 ) h_numEvents_perSys_wgtCMVA_6j4t->Fill(-1,wgt);
	}

	h_numEvents_perSys_wgtCMVA2_4j2t->Fill(iSys,wgtCMVA2);
	if( iSys==0 ) h_numEvents_perSys_wgtCMVA2_4j2t->Fill(-1,wgt);

	if( numJet>=6 && numBtag_CMVA>=4 ){
	  h_numEvents_perSys_wgtCMVA2_6j4t->Fill(iSys,wgtCMVA2);
	  if( iSys==0 ) h_numEvents_perSys_wgtCMVA2_6j4t->Fill(-1,wgt);
	}
      
	int this_category_cmva = -1;
	if( numJet==4 && numBtag_CMVA==2) this_category_cmva=1;
	if( numJet==5 && numBtag_CMVA==2) this_category_cmva=2;
	if( numJet>=6 && numBtag_CMVA==2) this_category_cmva=3;	
	if( numJet==4 && numBtag_CMVA==3) this_category_cmva=4;
	if( numJet==5 && numBtag_CMVA==3) this_category_cmva=5;
	if( numJet>=6 && numBtag_CMVA==3) this_category_cmva=6;
	if( numJet==4 && numBtag_CMVA>=4) this_category_cmva=7;
	if( numJet==5 && numBtag_CMVA>=4) this_category_cmva=8;
	if( numJet>=6 && numBtag_CMVA>=4) this_category_cmva=9;


	h_category_yield_cmva[iSys]->Fill(0.5,wgt);
	h_category_yield_cmva[iSys]->Fill(this_category_cmva,wgt);

	if( oneEle ){
	  h_category_yield_cmva_1e[iSys]->Fill(0.5,wgt);
	  h_category_yield_cmva_1e[iSys]->Fill(this_category_cmva,wgt);
	}
	else if( oneMu ){
	  h_category_yield_cmva_1m[iSys]->Fill(0.5,wgt);
	  h_category_yield_cmva_1m[iSys]->Fill(this_category_cmva,wgt);
	}


	h_category_yield_cmva_wgtCMVA[iSys]->Fill(0.5,wgtCMVA);
	h_category_yield_cmva_wgtCMVA[iSys]->Fill(this_category_cmva,wgtCMVA);

	if( oneEle ){
	  h_category_yield_cmva_wgtCMVA_1e[iSys]->Fill(0.5,wgtCMVA);
	  h_category_yield_cmva_wgtCMVA_1e[iSys]->Fill(this_category_cmva,wgtCMVA);
	}
	else if( oneMu ){
	  h_category_yield_cmva_wgtCMVA_1m[iSys]->Fill(0.5,wgtCMVA);
	  h_category_yield_cmva_wgtCMVA_1m[iSys]->Fill(this_category_cmva,wgtCMVA);
	}


	

	h_category_yield_cmva_wgtCMVA2[iSys]->Fill(0.5,wgtCMVA2);
	h_category_yield_cmva_wgtCMVA2[iSys]->Fill(this_category_cmva,wgtCMVA2);

	if( oneEle ){
	  h_category_yield_cmva_wgtCMVA2_1e[iSys]->Fill(0.5,wgtCMVA2);
	  h_category_yield_cmva_wgtCMVA2_1e[iSys]->Fill(this_category_cmva,wgtCMVA2);
	}
	else if( oneMu ){
	  h_category_yield_cmva_wgtCMVA2_1m[iSys]->Fill(0.5,wgtCMVA2);
	  h_category_yield_cmva_wgtCMVA2_1m[iSys]->Fill(this_category_cmva,wgtCMVA2);
	}
	
          if( verbose_ ) std::cout << " ===> test 3.3.8 " << std::endl;

	for( int kJet=0; kJet<int(ind_jet.size()); kJet++ ){
	  int iJet = ind_jet[kJet];

	  TLorentzVector myJet;
	  myJet.SetPtEtaPhiE( use_jet_pt[iJet], use_jet_eta[iJet], use_jet_phi[iJet], use_jet_energy[iJet] );

	  if( !firstJet ) sumJet = myJet;
	  else            sumJet += myJet;

	  double pt  = use_jet_pt[iJet];
	  double eta = use_jet_eta[iJet];
	  double phi = use_jet_phi[iJet];
	  int flavor = use_jet_hadronFlavour[iJet];

	  if( !(pt>30. && fabs(eta)<2.4) ) continue;

	  double dR = myJet.DeltaR(myLep);
	  if( dR < 0.4 ) continue;

	  double csv =  use_jet_csv[iJet];
	  if( csv < 0.0 ) csv = -0.05;
	  if( csv > 1.0 ) csv = 1.0;

	  double cmva =  use_jet_cmva[iJet];

	  double puMVA = use_jet_pileupJetId_fullDiscriminant[iJet];
	  if( puMVA < -1.0 ) puMVA = -1.0;
	  if( puMVA >  1.0 ) puMVA = 1.0;

	  h_cmva_jet_pt[0][iSys]->Fill(pt,wgt);
	  h_cmva_jet_eta[0][iSys]->Fill(eta,wgt);
	  h_cmva_jet_phi[0][iSys]->Fill(phi,wgt);
	  h_cmva_jet_csv[0][iSys]->Fill(csv,wgt);
	  h_cmva_jet_cmva[0][iSys]->Fill(cmva,wgt);
	  h_cmva_jet_puMVA[0][iSys]->Fill(puMVA,wgt);

	  h_cmva_jet_pt[this_category_cmva][iSys]->Fill(pt,wgt);
	  h_cmva_jet_eta[this_category_cmva][iSys]->Fill(eta,wgt);
	  h_cmva_jet_phi[this_category_cmva][iSys]->Fill(phi,wgt);
	  h_cmva_jet_csv[this_category_cmva][iSys]->Fill(csv,wgt);
	  h_cmva_jet_cmva[this_category_cmva][iSys]->Fill(cmva,wgt);
	  h_cmva_jet_puMVA[this_category_cmva][iSys]->Fill(puMVA,wgt);


	  // wgtCSV
	  h_cmva_jet_pt_wgtCMVA[0][iSys]->Fill(pt,wgtCMVA);
	  h_cmva_jet_eta_wgtCMVA[0][iSys]->Fill(eta,wgtCMVA);
	  h_cmva_jet_phi_wgtCMVA[0][iSys]->Fill(phi,wgtCMVA);
	  h_cmva_jet_csv_wgtCMVA[0][iSys]->Fill(csv,wgtCMVA);
	  h_cmva_jet_cmva_wgtCMVA[0][iSys]->Fill(cmva,wgtCMVA);
	  h_cmva_jet_puMVA_wgtCMVA[0][iSys]->Fill(puMVA,wgtCMVA);

	  h_cmva_jet_pt_wgtCMVA[this_category_cmva][iSys]->Fill(pt,wgtCMVA);
	  h_cmva_jet_eta_wgtCMVA[this_category_cmva][iSys]->Fill(eta,wgtCMVA);
	  h_cmva_jet_phi_wgtCMVA[this_category_cmva][iSys]->Fill(phi,wgtCMVA);
	  h_cmva_jet_csv_wgtCMVA[this_category_cmva][iSys]->Fill(csv,wgtCMVA);
	  h_cmva_jet_cmva_wgtCMVA[this_category_cmva][iSys]->Fill(cmva,wgtCMVA);
	  h_cmva_jet_puMVA_wgtCMVA[this_category_cmva][iSys]->Fill(puMVA,wgtCMVA);

	  if( abs(flavor)==5 ){
	    h_jet_cmva_bFlav[0][iSys]->Fill(cmva,wgt);
	    h_jet_cmva_bFlav[this_category_cmva][iSys]->Fill(cmva,wgt);

	    h_jet_cmva_wgtCMVA_bFlav[0][iSys]->Fill(cmva,wgtCMVA);
	    h_jet_cmva_wgtCMVA_bFlav[this_category_cmva][iSys]->Fill(cmva,wgtCMVA);
	  }
	  else if( abs(flavor)==4 ){
	    h_jet_cmva_cFlav[0][iSys]->Fill(cmva,wgt);
	    h_jet_cmva_cFlav[this_category_cmva][iSys]->Fill(cmva,wgt);

	    h_jet_cmva_wgtCMVA_cFlav[0][iSys]->Fill(cmva,wgtCMVA);
	    h_jet_cmva_wgtCMVA_cFlav[this_category_cmva][iSys]->Fill(cmva,wgtCMVA);
	  }
	  else{
	    h_jet_cmva_lFlav[0][iSys]->Fill(cmva,wgt);
	    h_jet_cmva_lFlav[this_category_cmva][iSys]->Fill(cmva,wgt);

	    h_jet_cmva_wgtCMVA_lFlav[0][iSys]->Fill(cmva,wgtCMVA);
	    h_jet_cmva_wgtCMVA_lFlav[this_category_cmva][iSys]->Fill(cmva,wgtCMVA);
	  }
	}
      }

    if( verbose_ ) std::cout << " ===> test 3.3.9 " << std::endl;

      
      // Require at least two b-tags
      if( !(numBtag>=2) ) continue;

      if( oneLep ) h_event_selection[iSys]->Fill(0.5+4, wgt); // all
      if( oneMu )  h_mu_event_selection[iSys]->Fill(0.5+4, wgt); // all
      if( oneEle ) h_ele_event_selection[iSys]->Fill(0.5+4, wgt); // all

      if( oneLep ) h_event_selection_wgtCSV[iSys]->Fill(0.5+4, wgtCSV); // all
      if( oneMu )  h_mu_event_selection_wgtCSV[iSys]->Fill(0.5+4, wgtCSV); // all
      if( oneEle ) h_ele_event_selection_wgtCSV[iSys]->Fill(0.5+4, wgtCSV); // all


      h_numEvents_perSys_wgtCSV_4j2t->Fill(iSys,wgtCSV);
      if( iSys==0 ) h_numEvents_perSys_wgtCSV_4j2t->Fill(-1,wgt);

      if( numJet>=6 && numBtag>=4 ){
	h_numEvents_perSys_wgtCSV_6j4t->Fill(iSys,wgtCSV);
	if( iSys==0 ) h_numEvents_perSys_wgtCSV_6j4t->Fill(-1,wgt);
      }

      
      h_numEvents_perSys_wgtCSV2_4j2t->Fill(iSys,wgtCSV2);
      if( iSys==0 ) h_numEvents_perSys_wgtCSV2_4j2t->Fill(-1,wgt);

      if( numJet>=6 && numBtag>=4 ){
	h_numEvents_perSys_wgtCSV2_6j4t->Fill(iSys,wgtCSV2);
	if( iSys==0 ) h_numEvents_perSys_wgtCSV2_6j4t->Fill(-1,wgt);
      }
      
      numEvents_wgt_gen_lumi_pu_4j2t += wgt_gen * wgt_lumi * wgt_pu;
      numEvents_wgt_gen_lumi_pu_csv_4j2t += wgt_gen * wgt_lumi * wgt_pu * wgt_csv_noSys;

      h_numJet_4j2t[iSys]->Fill(fill_numJet,wgt);
      h_numBtag_4j2t[iSys]->Fill(fill_numBtag,wgt);

      h_numJet_wgtCSV_4j2t[iSys]->Fill(fill_numJet,wgtCSV);
      h_numJet_wgtCSV2_4j2t[iSys]->Fill(fill_numJet,wgtCSV2);
      h_numJet_wgtCSV3_4j2t[iSys]->Fill(fill_numJet,wgtCSV3);
      h_numJet_wgtCSV4_4j2t[iSys]->Fill(fill_numJet,wgtCSV4);
      h_numJet_wgtCSV5_4j2t[iSys]->Fill(fill_numJet,wgtCSV5);


      h_numBtag_wgtCSV_4j2t[iSys]->Fill(fill_numBtag,wgtCSV);


      ////////////

    if( verbose_ ) std::cout << " ===> test 3.3.10 " << std::endl;



      // std::sort(jetPts.begin(), jetPts.end());
      // std::reverse(jetPts.begin(), jetPts.end());

      // if( jetPts.size()>=1 ) h_jet_1_pt->Fill(jetPts[0],wgt);
      // if( jetPts.size()>=2 ) h_jet_2_pt->Fill(jetPts[1],wgt);
      // if( jetPts.size()>=3 ) h_jet_3_pt->Fill(jetPts[2],wgt);
      // if( jetPts.size()>=4 ) h_jet_4_pt->Fill(jetPts[3],wgt);



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

      h_category_yield[iSys]->Fill(0.5,wgt);
      h_category_yield[iSys]->Fill(this_category,wgt);

      if( oneEle ){
	h_category_yield_1e[iSys]->Fill(0.5,wgt);
	h_category_yield_1e[iSys]->Fill(this_category,wgt);
      }
      else if( oneMu ){
	h_category_yield_1m[iSys]->Fill(0.5,wgt);
	h_category_yield_1m[iSys]->Fill(this_category,wgt);
      }


      h_category_yield_wgtCSV[iSys]->Fill(0.5,wgtCSV);
      h_category_yield_wgtCSV[iSys]->Fill(this_category,wgtCSV);

      if( oneEle ){
	h_category_yield_wgtCSV_1e[iSys]->Fill(0.5,wgtCSV);
	h_category_yield_wgtCSV_1e[iSys]->Fill(this_category,wgtCSV);
      }
      else if( oneMu ){
	h_category_yield_wgtCSV_1m[iSys]->Fill(0.5,wgtCSV);
	h_category_yield_wgtCSV_1m[iSys]->Fill(this_category,wgtCSV);
      }



      h_category_yield_wgtCSV2[iSys]->Fill(0.5,wgtCSV2);
      h_category_yield_wgtCSV2[iSys]->Fill(this_category,wgtCSV2);

      if( oneEle ){
	h_category_yield_wgtCSV2_1e[iSys]->Fill(0.5,wgtCSV2);
	h_category_yield_wgtCSV2_1e[iSys]->Fill(this_category,wgtCSV2);
      }
      else if( oneMu ){
	h_category_yield_wgtCSV2_1m[iSys]->Fill(0.5,wgtCSV2);
	h_category_yield_wgtCSV2_1m[iSys]->Fill(this_category,wgtCSV2);
      }



      h_category_yield_wgtCSV5[iSys]->Fill(0.5,wgtCSV5);
      h_category_yield_wgtCSV5[iSys]->Fill(this_category,wgtCSV5);

      if( oneEle ){
	h_category_yield_wgtCSV5_1e[iSys]->Fill(0.5,wgtCSV5);
	h_category_yield_wgtCSV5_1e[iSys]->Fill(this_category,wgtCSV5);
      }
      else if( oneMu ){
	h_category_yield_wgtCSV5_1m[iSys]->Fill(0.5,wgtCSV5);
	h_category_yield_wgtCSV5_1m[iSys]->Fill(this_category,wgtCSV5);
      }


      double lepPt = myLep.Pt();
      double lepEta = myLep.Eta();
      double lepPhi = myLep.Phi();

      double lepRelIso = eve->lepton_relIso_[lepInd];
      double lepRelIsoR04 = eve->lepton_relIsoR04_[lepInd];

      double lepTrigMVAOutput = eve->lepton_trigMVAOutput_[lepInd];

      double lepD0 = eve->lepton_d0_[lepInd];
      double lepDZ = eve->lepton_dZ_[lepInd];

      h_lepton_pt[0][iSys]->Fill(lepPt,wgt);
      h_lepton_eta[0][iSys]->Fill(lepEta,wgt);
      h_lepton_phi[0][iSys]->Fill(lepPhi,wgt);
      h_lepton_d0[0][iSys]->Fill(lepD0,wgt);
      h_lepton_dZ[0][iSys]->Fill(lepDZ,wgt);

      h_lepton_pt[this_category][iSys]->Fill(lepPt,wgt);
      h_lepton_eta[this_category][iSys]->Fill(lepEta,wgt);
      h_lepton_phi[this_category][iSys]->Fill(lepPhi,wgt);
      h_lepton_d0[this_category][iSys]->Fill(lepD0,wgt);
      h_lepton_dZ[this_category][iSys]->Fill(lepDZ,wgt);

      h_lepton_pt_wgtCSV[0][iSys]->Fill(lepPt,wgtCSV);
      h_lepton_eta_wgtCSV[0][iSys]->Fill(lepEta,wgtCSV);
      h_lepton_phi_wgtCSV[0][iSys]->Fill(lepPhi,wgtCSV);
      h_lepton_d0_wgtCSV[0][iSys]->Fill(lepD0,wgtCSV);
      h_lepton_dZ_wgtCSV[0][iSys]->Fill(lepDZ,wgtCSV);

      h_lepton_pt_wgtCSV[this_category][iSys]->Fill(lepPt,wgtCSV);
      h_lepton_eta_wgtCSV[this_category][iSys]->Fill(lepEta,wgtCSV);
      h_lepton_phi_wgtCSV[this_category][iSys]->Fill(lepPhi,wgtCSV);
      h_lepton_d0_wgtCSV[this_category][iSys]->Fill(lepD0,wgtCSV);
      h_lepton_dZ_wgtCSV[this_category][iSys]->Fill(lepDZ,wgtCSV);

      if( oneEle ){
	h_electron_pt[0][iSys]->Fill(lepPt,wgt);
	h_electron_eta[0][iSys]->Fill(lepEta,wgt);
	h_electron_phi[0][iSys]->Fill(lepPhi,wgt);
	h_electron_relIso[0][iSys]->Fill(lepRelIso,wgt);
	h_electron_trigMVAOutput[0][iSys]->Fill(lepTrigMVAOutput,wgt);
	h_electron_d0[0][iSys]->Fill(lepD0,wgt);
	h_electron_dZ[0][iSys]->Fill(lepDZ,wgt);

	h_electron_pt[this_category][iSys]->Fill(lepPt,wgt);
	h_electron_eta[this_category][iSys]->Fill(lepEta,wgt);
	h_electron_phi[this_category][iSys]->Fill(lepPhi,wgt);
	h_electron_relIso[this_category][iSys]->Fill(lepRelIso,wgt);
	h_electron_trigMVAOutput[this_category][iSys]->Fill(lepTrigMVAOutput,wgt);
	h_electron_d0[this_category][iSys]->Fill(lepD0,wgt);
	h_electron_dZ[this_category][iSys]->Fill(lepDZ,wgt);

	h_electron_pt_wgtCSV[0][iSys]->Fill(lepPt,wgtCSV);
	h_electron_eta_wgtCSV[0][iSys]->Fill(lepEta,wgtCSV);
	h_electron_phi_wgtCSV[0][iSys]->Fill(lepPhi,wgtCSV);
	h_electron_relIso_wgtCSV[0][iSys]->Fill(lepRelIso,wgtCSV);
	h_electron_trigMVAOutput_wgtCSV[0][iSys]->Fill(lepTrigMVAOutput,wgtCSV);
	h_electron_d0_wgtCSV[0][iSys]->Fill(lepD0,wgtCSV);
	h_electron_dZ_wgtCSV[0][iSys]->Fill(lepDZ,wgtCSV);

	h_electron_pt_wgtCSV[this_category][iSys]->Fill(lepPt,wgtCSV);
	h_electron_eta_wgtCSV[this_category][iSys]->Fill(lepEta,wgtCSV);
	h_electron_phi_wgtCSV[this_category][iSys]->Fill(lepPhi,wgtCSV);
	h_electron_relIso_wgtCSV[this_category][iSys]->Fill(lepRelIso,wgtCSV);
	h_electron_trigMVAOutput_wgtCSV[this_category][iSys]->Fill(lepTrigMVAOutput,wgtCSV);
	h_electron_d0_wgtCSV[this_category][iSys]->Fill(lepD0,wgtCSV);
	h_electron_dZ_wgtCSV[this_category][iSys]->Fill(lepDZ,wgtCSV);
      }
      else if( oneMu ){
	h_muon_pt[0][iSys]->Fill(lepPt,wgt);
	h_muon_eta[0][iSys]->Fill(lepEta,wgt);
	h_muon_phi[0][iSys]->Fill(lepPhi,wgt);
	h_muon_relIso[0][iSys]->Fill(lepRelIsoR04,wgt);
	h_muon_d0[0][iSys]->Fill(lepD0,wgt);
	h_muon_dZ[0][iSys]->Fill(lepDZ,wgt);

	h_muon_pt[this_category][iSys]->Fill(lepPt,wgt);
	h_muon_eta[this_category][iSys]->Fill(lepEta,wgt);
	h_muon_phi[this_category][iSys]->Fill(lepPhi,wgt);
	h_muon_relIso[this_category][iSys]->Fill(lepRelIsoR04,wgt);
	h_muon_d0[this_category][iSys]->Fill(lepD0,wgt);
	h_muon_dZ[this_category][iSys]->Fill(lepDZ,wgt);

	h_muon_pt_wgtCSV[0][iSys]->Fill(lepPt,wgtCSV);
	h_muon_eta_wgtCSV[0][iSys]->Fill(lepEta,wgtCSV);
	h_muon_phi_wgtCSV[0][iSys]->Fill(lepPhi,wgtCSV);
	h_muon_relIso_wgtCSV[0][iSys]->Fill(lepRelIsoR04,wgtCSV);
	h_muon_d0_wgtCSV[0][iSys]->Fill(lepD0,wgtCSV);
	h_muon_dZ_wgtCSV[0][iSys]->Fill(lepDZ,wgtCSV);

	h_muon_pt_wgtCSV[this_category][iSys]->Fill(lepPt,wgtCSV);
	h_muon_eta_wgtCSV[this_category][iSys]->Fill(lepEta,wgtCSV);
	h_muon_phi_wgtCSV[this_category][iSys]->Fill(lepPhi,wgtCSV);
	h_muon_relIso_wgtCSV[this_category][iSys]->Fill(lepRelIsoR04,wgtCSV);
	h_muon_d0_wgtCSV[this_category][iSys]->Fill(lepD0,wgtCSV);
	h_muon_dZ_wgtCSV[this_category][iSys]->Fill(lepDZ,wgtCSV);
      }


    if( verbose_ ) std::cout << " ===> test 3.3.11 " << std::endl;


      //
      // Fill MET
      //

      // pfMET
      h_pfMET_pt[0][iSys]->Fill(pfMET,wgt);
      h_pfMET_pt[this_category][iSys]->Fill(pfMET,wgt);

      h_pfMET_phi[0][iSys]->Fill(pfMET_phi,wgt);
      h_pfMET_phi[this_category][iSys]->Fill(pfMET_phi,wgt);

      h_pfMET_pt_wgtCSV[0][iSys]->Fill(pfMET,wgtCSV);
      h_pfMET_pt_wgtCSV[this_category][iSys]->Fill(pfMET,wgtCSV);

      h_pfMET_phi_wgtCSV[0][iSys]->Fill(pfMET_phi,wgtCSV);
      h_pfMET_phi_wgtCSV[this_category][iSys]->Fill(pfMET_phi,wgtCSV);


      // pfMETNoHF
      h_pfMETNoHF_pt[0][iSys]->Fill(pfMETNoHF,wgt);
      h_pfMETNoHF_pt[this_category][iSys]->Fill(pfMETNoHF,wgt);

      h_pfMETNoHF_phi[0][iSys]->Fill(pfMETNoHF_phi,wgt);
      h_pfMETNoHF_phi[this_category][iSys]->Fill(pfMETNoHF_phi,wgt);

      // puppiMET
      h_puppiMET_pt[0][iSys]->Fill(puppiMET,wgt);
      h_puppiMET_pt[this_category][iSys]->Fill(puppiMET,wgt);

      h_puppiMET_phi[0][iSys]->Fill(puppiMET_phi,wgt);
      h_puppiMET_phi[this_category][iSys]->Fill(puppiMET_phi,wgt);


      //MHT
      double mht = sumJet.Pt();
      double mht_phi = sumJet.Phi();

      h_mht_pt[0][iSys]->Fill(mht,wgt);
      h_mht_pt[this_category][iSys]->Fill(mht,wgt);

      h_mht_phi[0][iSys]->Fill(mht_phi,wgt);
      h_mht_phi[this_category][iSys]->Fill(mht_phi,wgt);

      h_HT[0][iSys]->Fill(HT30,wgt);
      h_HT[this_category][iSys]->Fill(HT30,wgt);


    if( verbose_ ) std::cout << " ===> test 3.3.12 " << std::endl;

      for( int kJet=0; kJet<int(ind_jet.size()); kJet++ ){
	int iJet = ind_jet[kJet];

	TLorentzVector myJet;
	myJet.SetPtEtaPhiE( use_jet_pt[iJet], use_jet_eta[iJet], use_jet_phi[iJet], use_jet_energy[iJet] );

	if( !firstJet ) sumJet = myJet;
	else            sumJet += myJet;

	double pt  = use_jet_pt[iJet];
	double eta = use_jet_eta[iJet];
	double phi = use_jet_phi[iJet];
	int flavor = use_jet_hadronFlavour[iJet];

	if( !(pt>30. && fabs(eta)<2.4) ) continue;

	double dR = myJet.DeltaR(myLep);
	if( dR < 0.4 ) continue;

	double csv =  use_jet_csv[iJet];
	if( csv < 0.0 ) csv = -0.05;
	if( csv > 1.0 ) csv = 1.0;

	double cmva =  use_jet_cmva[iJet];

	double puMVA = use_jet_pileupJetId_fullDiscriminant[iJet];
	if( puMVA < -1.0 ) puMVA = -1.0;
	if( puMVA >  1.0 ) puMVA = 1.0;

	h_jet_pt[0][iSys]->Fill(pt,wgt);
	h_jet_eta[0][iSys]->Fill(eta,wgt);
	h_jet_phi[0][iSys]->Fill(phi,wgt);
	h_jet_csv[0][iSys]->Fill(csv,wgt);
	h_jet_cmva[0][iSys]->Fill(cmva,wgt);
	h_jet_puMVA[0][iSys]->Fill(puMVA,wgt);

	h_jet_pt[this_category][iSys]->Fill(pt,wgt);
	h_jet_eta[this_category][iSys]->Fill(eta,wgt);
	h_jet_phi[this_category][iSys]->Fill(phi,wgt);
	h_jet_csv[this_category][iSys]->Fill(csv,wgt);
	h_jet_cmva[this_category][iSys]->Fill(cmva,wgt);
	h_jet_puMVA[this_category][iSys]->Fill(puMVA,wgt);


	// wgtCSV
	h_jet_pt_wgtCSV[0][iSys]->Fill(pt,wgtCSV);
	h_jet_eta_wgtCSV[0][iSys]->Fill(eta,wgtCSV);
	h_jet_phi_wgtCSV[0][iSys]->Fill(phi,wgtCSV);
	h_jet_csv_wgtCSV[0][iSys]->Fill(csv,wgtCSV);
	h_jet_cmva_wgtCSV[0][iSys]->Fill(cmva,wgtCSV);
	h_jet_puMVA_wgtCSV[0][iSys]->Fill(puMVA,wgtCSV);

	h_jet_pt_wgtCSV[this_category][iSys]->Fill(pt,wgtCSV);
	h_jet_eta_wgtCSV[this_category][iSys]->Fill(eta,wgtCSV);
	h_jet_phi_wgtCSV[this_category][iSys]->Fill(phi,wgtCSV);
	h_jet_csv_wgtCSV[this_category][iSys]->Fill(csv,wgtCSV);
	h_jet_cmva_wgtCSV[this_category][iSys]->Fill(cmva,wgtCSV);
	h_jet_puMVA_wgtCSV[this_category][iSys]->Fill(puMVA,wgtCSV);

    if( verbose_ ) std::cout << " ===> test 3.3.14 " << std::endl;

	// wgtCSV2
	h_jet_pt_wgtCSV2[0][iSys]->Fill(pt,wgtCSV2);
	h_jet_eta_wgtCSV2[0][iSys]->Fill(eta,wgtCSV2);
	h_jet_phi_wgtCSV2[0][iSys]->Fill(phi,wgtCSV2);
	h_jet_csv_wgtCSV2[0][iSys]->Fill(csv,wgtCSV2);
	h_jet_puMVA_wgtCSV2[0][iSys]->Fill(puMVA,wgtCSV2);

	h_jet_pt_wgtCSV2[this_category][iSys]->Fill(pt,wgtCSV2);
	h_jet_eta_wgtCSV2[this_category][iSys]->Fill(eta,wgtCSV2);
	h_jet_phi_wgtCSV2[this_category][iSys]->Fill(phi,wgtCSV2);
	h_jet_csv_wgtCSV2[this_category][iSys]->Fill(csv,wgtCSV2);
	h_jet_puMVA_wgtCSV2[this_category][iSys]->Fill(puMVA,wgtCSV2);

	// wgtCSV3
	h_jet_pt_wgtCSV3[0][iSys]->Fill(pt,wgtCSV3);
	h_jet_eta_wgtCSV3[0][iSys]->Fill(eta,wgtCSV3);
	h_jet_phi_wgtCSV3[0][iSys]->Fill(phi,wgtCSV3);
	h_jet_csv_wgtCSV3[0][iSys]->Fill(csv,wgtCSV3);
	h_jet_puMVA_wgtCSV3[0][iSys]->Fill(puMVA,wgtCSV3);

	h_jet_pt_wgtCSV3[this_category][iSys]->Fill(pt,wgtCSV3);
	h_jet_eta_wgtCSV3[this_category][iSys]->Fill(eta,wgtCSV3);
	h_jet_phi_wgtCSV3[this_category][iSys]->Fill(phi,wgtCSV3);
	h_jet_csv_wgtCSV3[this_category][iSys]->Fill(csv,wgtCSV3);
	h_jet_puMVA_wgtCSV3[this_category][iSys]->Fill(puMVA,wgtCSV3);

	// wgtCSV4
	h_jet_pt_wgtCSV4[0][iSys]->Fill(pt,wgtCSV4);
	h_jet_eta_wgtCSV4[0][iSys]->Fill(eta,wgtCSV4);
	h_jet_phi_wgtCSV4[0][iSys]->Fill(phi,wgtCSV4);
	h_jet_csv_wgtCSV4[0][iSys]->Fill(csv,wgtCSV4);
	h_jet_puMVA_wgtCSV4[0][iSys]->Fill(puMVA,wgtCSV4);

	h_jet_pt_wgtCSV4[this_category][iSys]->Fill(pt,wgtCSV4);
	h_jet_eta_wgtCSV4[this_category][iSys]->Fill(eta,wgtCSV4);
	h_jet_phi_wgtCSV4[this_category][iSys]->Fill(phi,wgtCSV4);
	h_jet_csv_wgtCSV4[this_category][iSys]->Fill(csv,wgtCSV4);
	h_jet_puMVA_wgtCSV4[this_category][iSys]->Fill(puMVA,wgtCSV4);

	// wgtCSV5
	h_jet_pt_wgtCSV5[0][iSys]->Fill(pt,wgtCSV5);
	h_jet_eta_wgtCSV5[0][iSys]->Fill(eta,wgtCSV5);
	h_jet_phi_wgtCSV5[0][iSys]->Fill(phi,wgtCSV5);
	h_jet_csv_wgtCSV5[0][iSys]->Fill(csv,wgtCSV5);
	h_jet_puMVA_wgtCSV5[0][iSys]->Fill(puMVA,wgtCSV5);

	h_jet_pt_wgtCSV5[this_category][iSys]->Fill(pt,wgtCSV5);
	h_jet_eta_wgtCSV5[this_category][iSys]->Fill(eta,wgtCSV5);
	h_jet_phi_wgtCSV5[this_category][iSys]->Fill(phi,wgtCSV5);
	h_jet_csv_wgtCSV5[this_category][iSys]->Fill(csv,wgtCSV5);
	h_jet_puMVA_wgtCSV5[this_category][iSys]->Fill(puMVA,wgtCSV5);

    if( verbose_ ) std::cout << " ===> test 3.3.15 " << std::endl;


	if( abs(flavor)==5 ){
	  h_jet_csv_bFlav[0][iSys]->Fill(csv,wgt);
	  h_jet_csv_bFlav[this_category][iSys]->Fill(csv,wgt);

	  h_jet_csv_wgtCSV_bFlav[0][iSys]->Fill(csv,wgtCSV);
	  h_jet_csv_wgtCSV_bFlav[this_category][iSys]->Fill(csv,wgtCSV);

	  h_jet_csv_wgtCSV2_bFlav[0][iSys]->Fill(csv,wgtCSV2);
	  h_jet_csv_wgtCSV2_bFlav[this_category][iSys]->Fill(csv,wgtCSV2);

	  h_jet_csv_wgtCSV3_bFlav[0][iSys]->Fill(csv,wgtCSV3);
	  h_jet_csv_wgtCSV3_bFlav[this_category][iSys]->Fill(csv,wgtCSV3);

	  h_jet_csv_wgtCSV4_bFlav[0][iSys]->Fill(csv,wgtCSV4);
	  h_jet_csv_wgtCSV4_bFlav[this_category][iSys]->Fill(csv,wgtCSV4);

	  h_jet_csv_wgtCSV5_bFlav[0][iSys]->Fill(csv,wgtCSV5);
	  h_jet_csv_wgtCSV5_bFlav[this_category][iSys]->Fill(csv,wgtCSV5);
	}
	else if( abs(flavor)==4 ){
	  h_jet_csv_cFlav[0][iSys]->Fill(csv,wgt);
	  h_jet_csv_cFlav[this_category][iSys]->Fill(csv,wgt);

	  h_jet_csv_wgtCSV_cFlav[0][iSys]->Fill(csv,wgtCSV);
	  h_jet_csv_wgtCSV_cFlav[this_category][iSys]->Fill(csv,wgtCSV);

	  h_jet_csv_wgtCSV2_cFlav[0][iSys]->Fill(csv,wgtCSV2);
	  h_jet_csv_wgtCSV2_cFlav[this_category][iSys]->Fill(csv,wgtCSV2);

	  h_jet_csv_wgtCSV3_cFlav[0][iSys]->Fill(csv,wgtCSV3);
	  h_jet_csv_wgtCSV3_cFlav[this_category][iSys]->Fill(csv,wgtCSV3);

	  h_jet_csv_wgtCSV4_cFlav[0][iSys]->Fill(csv,wgtCSV4);
	  h_jet_csv_wgtCSV4_cFlav[this_category][iSys]->Fill(csv,wgtCSV4);

	  h_jet_csv_wgtCSV5_cFlav[0][iSys]->Fill(csv,wgtCSV5);
	  h_jet_csv_wgtCSV5_cFlav[this_category][iSys]->Fill(csv,wgtCSV5);
	}
	else{
	  h_jet_csv_lFlav[0][iSys]->Fill(csv,wgt);
	  h_jet_csv_lFlav[this_category][iSys]->Fill(csv,wgt);

	  h_jet_csv_wgtCSV_lFlav[0][iSys]->Fill(csv,wgtCSV);
	  h_jet_csv_wgtCSV_lFlav[this_category][iSys]->Fill(csv,wgtCSV);

	  h_jet_csv_wgtCSV2_lFlav[0][iSys]->Fill(csv,wgtCSV2);
	  h_jet_csv_wgtCSV2_lFlav[this_category][iSys]->Fill(csv,wgtCSV2);

	  h_jet_csv_wgtCSV3_lFlav[0][iSys]->Fill(csv,wgtCSV3);
	  h_jet_csv_wgtCSV3_lFlav[this_category][iSys]->Fill(csv,wgtCSV3);

	  h_jet_csv_wgtCSV4_lFlav[0][iSys]->Fill(csv,wgtCSV4);
	  h_jet_csv_wgtCSV4_lFlav[this_category][iSys]->Fill(csv,wgtCSV4);

	  h_jet_csv_wgtCSV5_lFlav[0][iSys]->Fill(csv,wgtCSV5);
	  h_jet_csv_wgtCSV5_lFlav[this_category][iSys]->Fill(csv,wgtCSV5);

    if( verbose_ ) std::cout << " ===> test 3.3.16 " << std::endl;

	}
      }

    } // End loop over systematics

    if( verbose_ ) std::cout << " ===> test 3.3.13 " << std::endl;

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

  if( nPU>49 ) nPU = 49;

  if( iSys==0 ){
  PUscale[0] = 0.544141;
  PUscale[1] = 0.683152;
  PUscale[2] = 1.0588;
  PUscale[3] = 1.36648;
  PUscale[4] = 1.62214;
  PUscale[5] = 1.93965;
  PUscale[6] = 1.45916;
  PUscale[7] = 1.29301;
  PUscale[8] = 1.39344;
  PUscale[9] = 1.37393;
  PUscale[10] = 1.26638;
  PUscale[11] = 1.1637;
  PUscale[12] = 1.05777;
  PUscale[13] = 0.901478;
  PUscale[14] = 0.695434;
  PUscale[15] = 0.48665;
  PUscale[16] = 0.323436;
  PUscale[17] = 0.233179;
  PUscale[18] = 0.204042;
  PUscale[19] = 0.179695;
  PUscale[20] = 0.117906;
  PUscale[21] = 0.0535705;
  PUscale[22] = 0.0179761;
  PUscale[23] = 0.00486845;
  PUscale[24] = 0.00120473;
  PUscale[25] = 0.000306356;
  PUscale[26] = 8.44097e-05;
  PUscale[27] = 2.35712e-05;
  PUscale[28] = 6.33309e-06;
  PUscale[29] = 1.62779e-06;
  PUscale[30] = 3.92723e-07;
  PUscale[31] = 8.94453e-08;
  PUscale[32] = 1.90178e-08;
  PUscale[33] = 3.84203e-09;
  PUscale[34] = 7.38767e-10;
  PUscale[35] = 1.28456e-10;
  PUscale[36] = 2.08504e-11;
  PUscale[37] = 3.31182e-12;
  PUscale[38] = 4.63316e-13;
  PUscale[39] = 4.72854e-14;
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
  }
  else if( iSys==1 ){
  PUscale[0] = 0.458261;
  PUscale[1] = 0.626711;
  PUscale[2] = 0.971785;
  PUscale[3] = 1.21776;
  PUscale[4] = 1.34764;
  PUscale[5] = 1.4248;
  PUscale[6] = 0.922526;
  PUscale[7] = 0.744979;
  PUscale[8] = 0.913723;
  PUscale[9] = 1.06132;
  PUscale[10] = 1.11114;
  PUscale[11] = 1.12292;
  PUscale[12] = 1.13425;
  PUscale[13] = 1.10829;
  PUscale[14] = 1.00558;
  PUscale[15] = 0.83577;
  PUscale[16] = 0.639871;
  PUscale[17] = 0.469319;
  PUscale[18] = 0.361301;
  PUscale[19] = 0.306711;
  PUscale[20] = 0.237734;
  PUscale[21] = 0.142901;
  PUscale[22] = 0.0656738;
  PUscale[23] = 0.0238299;
  PUscale[24] = 0.00722283;
  PUscale[25] = 0.00197294;
  PUscale[26] = 0.000550982;
  PUscale[27] = 0.000165108;
  PUscale[28] = 5.14307e-05;
  PUscale[29] = 1.60005e-05;
  PUscale[30] = 4.76288e-06;
  PUscale[31] = 1.35321e-06;
  PUscale[32] = 3.62213e-07;
  PUscale[33] = 9.29401e-08;
  PUscale[34] = 2.28987e-08;
  PUscale[35] = 5.14682e-09;
  PUscale[36] = 1.08953e-09;
  PUscale[37] = 2.27672e-10;
  PUscale[38] = 4.24244e-11;
  PUscale[39] = 6.6564e-12;
  PUscale[40] = 1.21726e-12;
  PUscale[41] = 1.25024e-13;
  PUscale[42] = 0;
  PUscale[43] = 0;
  PUscale[44] = 0;
  PUscale[45] = 0;
  PUscale[46] = 0;
  PUscale[47] = 0;
  PUscale[48] = 0;
  PUscale[49] = 0;
  }
  else if( iSys==-1 ){
  PUscale[0] = 0.637352;
  PUscale[1] = 0.745524;
  PUscale[2] = 1.16977;
  PUscale[3] = 1.55999;
  PUscale[4] = 2.01896;
  PUscale[5] = 2.79883;
  PUscale[6] = 2.50005;
  PUscale[7] = 2.22714;
  PUscale[8] = 2.00474;
  PUscale[9] = 1.67907;
  PUscale[10] = 1.37446;
  PUscale[11] = 1.13578;
  PUscale[12] = 0.899819;
  PUscale[13] = 0.646461;
  PUscale[14] = 0.412783;
  PUscale[15] = 0.243661;
  PUscale[16] = 0.153684;
  PUscale[17] = 0.125858;
  PUscale[18] = 0.116543;
  PUscale[19] = 0.0860169;
  PUscale[20] = 0.0414331;
  PUscale[21] = 0.0133115;
  PUscale[22] = 0.00326222;
  PUscale[23] = 0.000723524;
  PUscale[24] = 0.000169442;
  PUscale[25] = 4.20509e-05;
  PUscale[26] = 1.04015e-05;
  PUscale[27] = 2.41401e-06;
  PUscale[28] = 5.19674e-07;
  PUscale[29] = 1.0525e-07;
  PUscale[30] = 1.97832e-08;
  PUscale[31] = 3.47419e-09;
  PUscale[32] = 5.63765e-10;
  PUscale[33] = 8.60403e-11;
  PUscale[34] = 1.23701e-11;
  PUscale[35] = 1.59276e-12;
  PUscale[36] = 1.89108e-13;
  PUscale[37] = 1.98468e-14;
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
  }

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



void fillCSVhistos_old(TFile* fileHF, TFile* fileLF){

  for( int iSys=0; iSys<9; iSys++ ){
    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) h_csv_wgt_hf_old[iSys][iPt] = NULL;
    for( int iPt=0; iPt<3; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf_old[iSys][iPt][iEta] = NULL;
    }
  }
  for( int iSys=0; iSys<5; iSys++ ){
    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) c_csv_wgt_hf_old[iSys][iPt] = NULL;
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

    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) h_csv_wgt_hf_old[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );

    if( iSys<5 ){
      for( int iPt=0; iPt<PtBinsHF_; iPt++ ) c_csv_wgt_hf_old[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
    }
    
    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf_old[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
    }
  }

  return;
}


void fillCMVAhistos(TFile* fileHF, TFile* fileLF){

  for( int iSys=0; iSys<9; iSys++ ){
    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) h_cmva_wgt_hf[iSys][iPt] = NULL;
    for( int iPt=0; iPt<3; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_cmva_wgt_lf[iSys][iPt][iEta] = NULL;
    }
  }
  for( int iSys=0; iSys<5; iSys++ ){
    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) c_cmva_wgt_hf[iSys][iPt] = NULL;
  }

  // CSV reweighting /// only care about the nominal ones
  for( int iSys=0; iSys<9; iSys++ ){
    TString syst_cmva_suffix_hf = "final";
    TString syst_cmva_suffix_c = "final";
    TString syst_cmva_suffix_lf = "final";
    
    switch( iSys ){
    case 0:
      // this is the nominal case
      break;
    case 1:
      // JESUp
      syst_cmva_suffix_hf = "final_JESUp"; syst_cmva_suffix_lf = "final_JESUp";
      syst_cmva_suffix_c  = "final_cErr1Up";
      break;
    case 2:
      // JESDown
      syst_cmva_suffix_hf = "final_JESDown"; syst_cmva_suffix_lf = "final_JESDown";
      syst_cmva_suffix_c  = "final_cErr1Down";
      break;
    case 3:
      // purity up
      syst_cmva_suffix_hf = "final_LFUp"; syst_cmva_suffix_lf = "final_HFUp";
      syst_cmva_suffix_c  = "final_cErr2Up";
      break;
    case 4:
      // purity down
      syst_cmva_suffix_hf = "final_LFDown"; syst_cmva_suffix_lf = "final_HFDown";
      syst_cmva_suffix_c  = "final_cErr2Down";
      break;
    case 5:
      // stats1 up
      syst_cmva_suffix_hf = "final_Stats1Up"; syst_cmva_suffix_lf = "final_Stats1Up";
      break;
    case 6:
      // stats1 down
      syst_cmva_suffix_hf = "final_Stats1Down"; syst_cmva_suffix_lf = "final_Stats1Down";
      break;
    case 7:
      // stats2 up
      syst_cmva_suffix_hf = "final_Stats2Up"; syst_cmva_suffix_lf = "final_Stats2Up";
      break;
    case 8:
      // stats2 down
      syst_cmva_suffix_hf = "final_Stats2Down"; syst_cmva_suffix_lf = "final_Stats2Down";
      break;
    }

    for( int iPt=0; iPt<PtBinsHF_; iPt++ ) h_cmva_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_cmva_suffix_hf.Data()) );

    if( iSys<5 ){
      for( int iPt=0; iPt<PtBinsHF_; iPt++ ) c_cmva_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_cmva_suffix_c.Data()) );
    }
    
    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_cmva_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_cmva_suffix_lf.Data()) );
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



double get_csv_wgt_old( vdouble jetPt, vdouble jetEta, vdouble jetCSV, vint jetFlavor, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF ){

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
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf_old[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf_old[iSysHF][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;

      // if( iSysHF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      //  			     flavor, pt, jetAbsEta, csv, iCSVWgtHF );
    }
    else if( abs(flavor) == 4 ){
      int useCSVBin = (csv>=0.) ? c_csv_wgt_hf_old[iSysC][iPt]->FindBin(csv) : 1;
      double iCSVWgtC = c_csv_wgt_hf_old[iSysC][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;
      // if( iSysC==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      // 			    flavor, pt, jetAbsEta, csv, iCSVWgtC );
    }
    else {
      if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf_old[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf_old[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
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



double get_cmva_wgt( vdouble jetPt, vdouble jetEta, vdouble jetCMVA, vint jetFlavor, int iSys, double &cmvaWgtHF, double &cmvaWgtLF, double &cmvaWgtCF ){

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

  double cmvaWgthf = 1.;
  double cmvaWgtC  = 1.;
  double cmvaWgtlf = 1.;

  for( int iJet=0; iJet<int(jetPt.size()); iJet++ ){

    double cmva = jetCMVA[iJet];
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

    //std::cout << " flavor = " << flavor << ", cmva = " << cmva << ", iPt = " << iPt << ", iEta = " << iEta << ", iSysHF = " << iSysHF << ", iSysC = " << iSysC << ", iSysLF = " << iSysLF << std::endl;

    if (abs(flavor) == 5 ){
      int useCMVABin = h_cmva_wgt_hf[iSysHF][iPt]->FindBin(cmva);
      double iCMVAWgtHF = h_cmva_wgt_hf[iSysHF][iPt]->GetBinContent(useCMVABin);
      if( iCMVAWgtHF!=0 ) cmvaWgthf *= iCMVAWgtHF;

      // if( iSysHF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t cmva=%.3f,\t wgt=%.2f \n",
      //  			     flavor, pt, jetAbsEta, cmva, iCMVAWgtHF );
    }
    else if( abs(flavor) == 4 ){
      int useCMVABin = c_cmva_wgt_hf[iSysC][iPt]->FindBin(cmva);
      double iCMVAWgtC = c_cmva_wgt_hf[iSysC][iPt]->GetBinContent(useCMVABin);
      if( iCMVAWgtC!=0 ) cmvaWgtC *= iCMVAWgtC;
      // if( iSysC==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t cmva=%.3f,\t wgt=%.2f \n",
      // 			    flavor, pt, jetAbsEta, cmva, iCMVAWgtC );
    }
    else {
      if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
      int useCMVABin = h_cmva_wgt_lf[iSysLF][iPt][iEta]->FindBin(cmva);
      double iCMVAWgtLF = h_cmva_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCMVABin);
      if( iCMVAWgtLF!=0 ) cmvaWgtlf *= iCMVAWgtLF;
      // if( iSysLF==0 ) printf(" iSysLF = %d, iJet = %d,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t cmva=%.3f,\t wgt=%.2f \n",
      //  			     iSysLF, iJet, flavor, pt, jetAbsEta, cmva, iCMVAWgtLF );
    }
  }

  double cmvaWgtTotal = cmvaWgthf * cmvaWgtC * cmvaWgtlf;

  cmvaWgtHF = cmvaWgthf;
  cmvaWgtLF = cmvaWgtlf;
  cmvaWgtCF = cmvaWgtC;

  //std::cout << "\t iSys = " << iSys << " total cmvaWgtLF = " << cmvaWgtLF << std::endl;

  return cmvaWgtTotal;
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

    bool isTag = ( csv > 0.80 );

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
    else {
      if( abs(flavor)==5 || abs(flavor)==4 ){
	if( sys_name.Contains("LFUp") )         SF = reader_btv_dil_up.eval(jf, eta, pt, csv);
	else if( sys_name.Contains("LFDown") )  SF = reader_btv_dil_down.eval(jf, eta, pt, csv);
	else                                    SF = reader_btv_dil.eval(jf, eta, pt);
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
