#include "TFile.h"
#include "TChain.h"
#include "THStack.h"
#include "TF1.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TList.h"
#include "TLatex.h"
#include "TLine.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TKey.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>

#include "TMath.h"

//*****************************************************************************

//*****************************************************************************


void makePlots_ttHbb_data2mc( TString dirpostfix_ = "", bool printPDF_ = false, bool includeCatPlots_ = true, bool renormMC_ = false, bool includeFlavPlots_ = false, double rescaleMC_ = -1 ){

  TH1::SetDefaultSumw2();

  TString histo_dir_prefix = "HistoFiles";
  //histo_dir_prefix = "backup_2015_12_12_HistoFiles";

  int NumSamples = 12;
  TFile* file[NumSamples];
  file[0] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_SingleLepton_Run2015D_16Dec2015_histo.root");
  file[7] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttV_histo.root");
  file[6] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_singlet_histo.root");

  file[8] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_QCD_HT100toInf_histo.root");

  if( dirpostfix_.Contains("HT") ){
    file[9] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_EWK_HTbins_histo.root");
  }
  else{
    file[9] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_EWK_histo.root");
  }

  if( dirpostfix_.Contains("madgraph") ){
    file[5] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttlf_madgraph_histo.root");
    file[4] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttcc_madgraph_histo.root");
    file[3] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttb_madgraph_histo.root");
    file[2] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_tt2b_madgraph_histo.root");
    file[1] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttbb_madgraph_histo.root");
  }
  else if( dirpostfix_.Contains("amcatnlo") ){
    file[5] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttlf_amcatnlo_histo.root");
    file[4] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttcc_amcatnlo_histo.root");
    file[3] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttb_amcatnlo_histo.root");
    file[2] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_tt2b_amcatnlo_histo.root");
    file[1] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttbb_amcatnlo_histo.root");
  }
  else {
    file[5] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttlf_powheg_histo.root");
    file[4] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttcc_powheg_histo.root");
    file[3] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttb_powheg_histo.root");
    file[2] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_tt2b_powheg_histo.root");
    file[1] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttbb_powheg_histo.root");
  }

  file[10] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttHnonbb_powheg_histo.root");
  file[11] = new TFile(histo_dir_prefix + "/ttHbb_data2mc_treeReader_ttHTobb_powheg_histo.root");



  std::vector<TString> histLabels(NumSamples);
  histLabels[0] = "Data    ";
  histLabels[9] = "EWK";
  histLabels[8] = "QCD";
  histLabels[7] = "tt+W,Z";
  histLabels[6] = "single t";
  histLabels[5] = "tt+lf";
  histLabels[4] = "tt+cc";
  histLabels[3] = "tt+b";
  histLabels[2] = "tt+2b";
  histLabels[1] = "tt+bb";
  histLabels[10] = "ttHnonbb";
  histLabels[11] = "ttHbb";


  Color_t color[NumSamples];
  color[0] = kBlack;
  color[9] = kAzure+2;
  color[8] = kOrange+1;
  color[7] = kBlue-10;
  color[6] = kMagenta;
  color[5] = kRed-7;
  color[4] = kRed+1;
  color[3] = kRed-2;
  color[2] = kRed+0;
  color[1] = kRed+3;
  // color[9] = kGreen+2;
  // color[10] = kBlue;
  color[11] = kGreen+2;
  color[10] = kBlue;




  TString dirprefix = "Images/Images_2016_03_11_ttHbb_data2mc" + dirpostfix_;
  if( renormMC_ ) dirprefix += "_renormMC";

  dirprefix += "/";


  struct stat st;
  if( stat(dirprefix.Data(),&st) != 0 )  mkdir(dirprefix.Data(),0777);


 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  // sys_cat_labels.push_back("_muFUp");           //25
  // sys_cat_labels.push_back("_muFDown");         //26
  // sys_cat_labels.push_back("_muRUp");           //27
  // sys_cat_labels.push_back("_muRDown");         //28
  sys_cat_labels.push_back("_muRmuFUp");        //29
  sys_cat_labels.push_back("_muRmuFDown");      //30
  sys_cat_labels.push_back("_NNPDFUp");        //31
  sys_cat_labels.push_back("_NNPDFDown");      //32

  int NumSysCat = int(sys_cat_labels.size());


  // has cat and sys
  std::vector<TString> histoname1;

  histoname1.push_back("h_jet_csv");
  histoname1.push_back("h_jet_csv_wgtCSV");


  histoname1.push_back("h_lepton_pt");
  histoname1.push_back("h_lepton_eta");
  histoname1.push_back("h_lepton_phi");

  histoname1.push_back("h_lepton_pt_wgtCSV");
  histoname1.push_back("h_lepton_eta_wgtCSV");
  histoname1.push_back("h_lepton_phi_wgtCSV");
  histoname1.push_back("h_lepton_d0_wgtCSV");
  histoname1.push_back("h_lepton_dZ_wgtCSV");

  histoname1.push_back("h_electron_pt_wgtCSV");
  histoname1.push_back("h_electron_eta_wgtCSV");
  histoname1.push_back("h_electron_phi_wgtCSV");
  histoname1.push_back("h_electron_relIso_wgtCSV");
  histoname1.push_back("h_electron_trigMVAOutput_wgtCSV");
  histoname1.push_back("h_electron_d0_wgtCSV");
  histoname1.push_back("h_electron_dZ_wgtCSV");

  histoname1.push_back("h_muon_pt_wgtCSV");
  histoname1.push_back("h_muon_eta_wgtCSV");
  histoname1.push_back("h_muon_phi_wgtCSV");
  histoname1.push_back("h_muon_relIso_wgtCSV");
  histoname1.push_back("h_muon_d0_wgtCSV");
  histoname1.push_back("h_muon_dZ_wgtCSV");

  histoname1.push_back("h_jet_pt");
  histoname1.push_back("h_jet_eta");
  histoname1.push_back("h_jet_phi");
  histoname1.push_back("h_jet_puMVA");

  histoname1.push_back("h_jet_pt_wgtCSV");
  histoname1.push_back("h_jet_eta_wgtCSV");
  histoname1.push_back("h_jet_phi_wgtCSV");
  histoname1.push_back("h_jet_puMVA_wgtCSV");


  
  histoname1.push_back("h_pfMET_pt");
  histoname1.push_back("h_pfMET_phi");

  histoname1.push_back("h_pfMET_pt_wgtCSV");
  histoname1.push_back("h_pfMET_phi_wgtCSV");
  
  /*

  histoname1.push_back("h_jet_csv_wgtCSV2");
  histoname1.push_back("h_jet_csv_wgtCSV3");
  histoname1.push_back("h_jet_csv_wgtCSV4");

  histoname1.push_back("h_electron_pt");
  histoname1.push_back("h_electron_eta");
  histoname1.push_back("h_electron_phi");
  histoname1.push_back("h_electron_relIso");
  histoname1.push_back("h_electron_trigMVAOutput");

  histoname1.push_back("h_muon_pt");
  histoname1.push_back("h_muon_eta");
  histoname1.push_back("h_muon_phi");
  histoname1.push_back("h_muon_relIso");

  histoname1.push_back("h_jet_pt_wgtCSV2");
  histoname1.push_back("h_jet_eta_wgtCSV2");
  histoname1.push_back("h_jet_phi_wgtCSV2");
  histoname1.push_back("h_jet_puMVA_wgtCSV2");

  histoname1.push_back("h_jet_pt_wgtCSV3");
  histoname1.push_back("h_jet_eta_wgtCSV3");
  histoname1.push_back("h_jet_phi_wgtCSV3");
  histoname1.push_back("h_jet_puMVA_wgtCSV3");

  histoname1.push_back("h_jet_pt_wgtCSV4");
  histoname1.push_back("h_jet_eta_wgtCSV4");
  histoname1.push_back("h_jet_phi_wgtCSV4");
  histoname1.push_back("h_jet_puMVA_wgtCSV4");

  histoname1.push_back("h_pfMETNoHF_pt");
  histoname1.push_back("h_pfMETNoHF_phi");
  histoname1.push_back("h_puppiMET_pt");
  histoname1.push_back("h_puppiMET_phi");

  histoname1.push_back("h_mht_pt");
  histoname1.push_back("h_mht_phi");

  histoname1.push_back("h_HT");
  */

  //
  // Energy sum plots
  //


  // has sys but not cat
  std::vector<TString> histoname2;

  // histoname2.push_back("h_numBtag_4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV_4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV2_4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV3_4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV4_4j_met30");

  // histoname2.push_back("h_numBtag_1e_4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV_1e_4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV2_1e_4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV3_1e_4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV4_1e_4j_met30");

  // histoname2.push_back("h_numBtag_1m_4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV_1m_4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV2_1m_4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV3_1m_4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV4_1m_4j_met30");


  histoname2.push_back("h_numBtag_4j");
  histoname2.push_back("h_numBtag_wgtCSV_4j");
  histoname2.push_back("h_numBtag_wgtCSV2_4j");
  // histoname2.push_back("h_numBtag_wgtCSV3_4j");
  // histoname2.push_back("h_numBtag_wgtCSV4_4j");
  histoname2.push_back("h_numBtag_wgtCSV5_4j");

  histoname2.push_back("h_numBtag_1e_4j");
  histoname2.push_back("h_numBtag_wgtCSV_1e_4j");
  // histoname2.push_back("h_numBtag_wgtCSV2_1e_4j");
  // histoname2.push_back("h_numBtag_wgtCSV3_1e_4j");
  // histoname2.push_back("h_numBtag_wgtCSV4_1e_4j");

  histoname2.push_back("h_numBtag_1m_4j");
  histoname2.push_back("h_numBtag_wgtCSV_1m_4j");
  // histoname2.push_back("h_numBtag_wgtCSV2_1m_4j");
  // histoname2.push_back("h_numBtag_wgtCSV3_1m_4j");
  // histoname2.push_back("h_numBtag_wgtCSV4_1m_4j");


  histoname2.push_back("h_numBtag_CMVA_4j");
  histoname2.push_back("h_numBtag_CMVA_wgtCMVA_4j");

  histoname2.push_back("h_numJet_4j2t");
  histoname2.push_back("h_numJet_wgtCSV_4j2t");
  // histoname2.push_back("h_numJet_wgtCSV2_4j2t");
  // histoname2.push_back("h_numJet_wgtCSV3_4j2t");
  // histoname2.push_back("h_numJet_wgtCSV4_4j2t");
  histoname2.push_back("h_numJet_wgtCSV5_4j2t");

  // histoname2.push_back("h_numBtag_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV2_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV3_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV4_4j_met50");

  // histoname2.push_back("h_numBtag_1m_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV_1m_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV2_1m_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV3_1m_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV4_1m_4j_met50");

  // histoname2.push_back("h_numBtag_1e_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV_1e_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV2_1e_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV3_1e_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV4_1e_4j_met50");

  histoname2.push_back("h_category_yield");
  histoname2.push_back("h_category_yield_1e");
  histoname2.push_back("h_category_yield_1m");

  histoname2.push_back("h_category_yield_wgtCSV");
  histoname2.push_back("h_category_yield_wgtCSV_1e");
  histoname2.push_back("h_category_yield_wgtCSV_1m");

  histoname2.push_back("h_category_yield_wgtCSV5");
  histoname2.push_back("h_category_yield_wgtCSV5_1e");
  histoname2.push_back("h_category_yield_wgtCSV5_1m");

  
  histoname2.push_back("h_category_yield_cmva");
  histoname2.push_back("h_category_yield_cmva_1e");
  histoname2.push_back("h_category_yield_cmva_1m");

  histoname2.push_back("h_category_yield_cmva_wgtCMVA");
  histoname2.push_back("h_category_yield_cmva_wgtCMVA_1e");
  histoname2.push_back("h_category_yield_cmva_wgtCMVA_1m");

  // histoname2.push_back("h_event_selection");
  // histoname2.push_back("h_mu_event_selection");
  // histoname2.push_back("h_ele_event_selection");

  // histoname2.push_back("h_event_selection_wgtCSV");
  // histoname2.push_back("h_mu_event_selection_wgtCSV");
  // histoname2.push_back("h_ele_event_selection_wgtCSV");

  // histoname2.push_back("h_numJet");
  // histoname2.push_back("h_numBtag");
  // histoname2.push_back("h_numJet_wgtCSV");
  // histoname2.push_back("h_numBtag_wgtCSV");

  // histoname2.push_back("h_numBtag_4j2t");
  // histoname2.push_back("h_numBtag_wgtCSV_4j2t");

  histoname2.push_back("h_numBtag_eq4j");
  histoname2.push_back("h_numBtag_wgtCSV_eq4j");
  // histoname2.push_back("h_numBtag_wgtCSV2_eq4j");
  // histoname2.push_back("h_numBtag_wgtCSV3_eq4j");
  // histoname2.push_back("h_numBtag_wgtCSV4_eq4j");

  // histoname2.push_back("h_numBtag_2l");
  // histoname2.push_back("h_numBtag_wgtCSV_2l");
  // histoname2.push_back("h_numBtag_wgtCSV2_2l");
  // histoname2.push_back("h_numBtag_wgtCSV3_2l");
  // histoname2.push_back("h_numBtag_wgtCSV4_2l");

  // histoname2.push_back("h_numBtag_2l_met30");
  // histoname2.push_back("h_numBtag_wgtCSV_2l_met30");
  // histoname2.push_back("h_numBtag_wgtCSV2_2l_met30");
  // histoname2.push_back("h_numBtag_wgtCSV3_2l_met30");
  // histoname2.push_back("h_numBtag_wgtCSV4_2l_met30");


  // histoname2.push_back("h_numJet_2l");
  // histoname2.push_back("h_numJet_wgtCSV_2l");
  // histoname2.push_back("h_numJet_wgtCSV2_2l");
  // histoname2.push_back("h_numJet_wgtCSV3_2l");
  // histoname2.push_back("h_numJet_wgtCSV4_2l");

  // histoname2.push_back("h_numJet_2l_met30");
  // histoname2.push_back("h_numJet_wgtCSV_2l_met30");
  // histoname2.push_back("h_numJet_wgtCSV2_2l_met30");
  // histoname2.push_back("h_numJet_wgtCSV3_2l_met30");
  // histoname2.push_back("h_numJet_wgtCSV4_2l_met30");


  // histoname2.push_back("h_numBtag_2l_leq2j");
  // histoname2.push_back("h_numBtag_wgtCSV_2l_leq2j");
  // histoname2.push_back("h_numBtag_wgtCSV2_2l_leq2j");
  // histoname2.push_back("h_numBtag_wgtCSV3_2l_leq2j");
  // histoname2.push_back("h_numBtag_wgtCSV4_2l_leq2j");

  // histoname2.push_back("h_numBtag_2l_met30_leq2j");
  // histoname2.push_back("h_numBtag_wgtCSV_2l_met30_leq2j");
  // histoname2.push_back("h_numBtag_wgtCSV2_2l_met30_leq2j");
  // histoname2.push_back("h_numBtag_wgtCSV3_2l_met30_leq2j");
  // histoname2.push_back("h_numBtag_wgtCSV4_2l_met30_leq2j");


  // //histoname2.push_back("h_deltaR_jet_lep");

  histoname2.push_back("h_numPVs_wgt");
  histoname2.push_back("h_numPVs_noPUwgt");


  histoname2.push_back("h_PV_x_wgt");
  histoname2.push_back("h_PV_y_wgt");
  histoname2.push_back("h_PV_z_wgt");

  // histoname2.push_back("h_pfMETNoHF_pt_2l");
  // histoname2.push_back("h_pfMETNoHF_pt_2e");
  // histoname2.push_back("h_pfMETNoHF_pt_2m");

  // histoname2.push_back("h_ll_diLepMass");
  // histoname2.push_back("h_ee_diLepMass");
  // histoname2.push_back("h_mm_diLepMass");


  // histoname2.push_back("h_pfMETNoHF_pt_4j_1l");
  // histoname2.push_back("h_pfMETNoHF_pt_4j_1e");
  // histoname2.push_back("h_pfMETNoHF_pt_4j_1m");


  // histoname2.push_back("h_numBtag_eq4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV_eq4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV2_eq4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV3_eq4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV4_eq4j_met30");


  // histoname2.push_back("h_numBtag_1e_eq4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV_1e_eq4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV2_1e_eq4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV3_1e_eq4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV4_1e_eq4j_met30");


  // histoname2.push_back("h_numBtag_1m_eq4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV_1m_eq4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV2_1m_eq4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV3_1m_eq4j_met30");
  // histoname2.push_back("h_numBtag_wgtCSV4_1m_eq4j_met30");






  // histoname2.push_back("h_numBtag_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV2_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV3_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV4_4j_met50");

  // histoname2.push_back("h_numBtag_eq4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV_eq4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV2_eq4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV3_eq4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV4_eq4j_met50");


  // histoname2.push_back("h_numBtag_1e_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV_1e_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV2_1e_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV3_1e_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV4_1e_4j_met50");

  // histoname2.push_back("h_numBtag_1e_eq4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV_1e_eq4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV2_1e_eq4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV3_1e_eq4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV4_1e_eq4j_met50");


  // histoname2.push_back("h_numBtag_1m_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV_1m_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV2_1m_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV3_1m_4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV4_1m_4j_met50");

  // histoname2.push_back("h_numBtag_1m_eq4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV_1m_eq4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV2_1m_eq4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV3_1m_eq4j_met50");
  // histoname2.push_back("h_numBtag_wgtCSV4_1m_eq4j_met50");




  for( int iHist=0; iHist<int(histoname1.size()); iHist++ ){
    for( int iCat=0; iCat<NumCat; iCat++ ){
      TString suffix = "_" + cat_labels[iCat];
      TString new_name = histoname1[iHist] + suffix;
      if( includeCatPlots_ && (new_name.Contains("d0") || new_name.Contains("dZ")) ) histoname2.push_back(new_name);
    }
  }


  std::vector<TString> histoname3;
  //histoname3.push_back("h_category_yield");
  // histoname3.push_back("h_category_yield_1e");
  // histoname3.push_back("h_category_yield_1m");

  histoname3.push_back("h_category_yield_wgtCSV");
  histoname3.push_back("h_category_yield_cmva_wgtCMVA");
  //histoname3.push_back("h_category_yield_wgtCSV5");
  // histoname3.push_back("h_category_yield_wgtCSV_1e");
  // histoname3.push_back("h_category_yield_wgtCSV_1m");


  std::vector<TString> histoname4;

  // histoname4.push_back("h_numBtag_2l_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV_2l_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV2_2l_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV3_2l_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV4_2l_leq2j");

  // histoname4.push_back("h_numBtag_2l_met30_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV_2l_met30_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV2_2l_met30_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV3_2l_met30_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV4_2l_met30_leq2j");


  // histoname4.push_back("h_numBtag_2l_geq1j_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV_2l_geq1j_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV2_2l_geq1j_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV3_2l_geq1j_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV4_2l_geq1j_leq2j");

  // histoname4.push_back("h_numBtag_2l_met30_geq1j_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV_2l_met30_geq1j_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV2_2l_met30_geq1j_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV3_2l_met30_geq1j_leq2j");
  // histoname4.push_back("h_numBtag_wgtCSV4_2l_met30_geq1j_leq2j");


  std::vector<TString> flavor_suffix;
  //flavor_suffix.push_back("0b0c0l");
  flavor_suffix.push_back("1b0c0l");
  flavor_suffix.push_back("0b1c0l");
  flavor_suffix.push_back("0b0c1l");
  flavor_suffix.push_back("2b0c0l");
  flavor_suffix.push_back("1b1c0l");
  flavor_suffix.push_back("1b0c1l");
  flavor_suffix.push_back("0b2c0l");
  flavor_suffix.push_back("0b1c1l");
  flavor_suffix.push_back("0b0c2l");

  int NumFlavComb = int(flavor_suffix.size());

  Color_t flavcolor[NumFlavComb];
  flavcolor[0] = kAzure+2;
  flavcolor[1] = kBlue-10;
  flavcolor[2] = kMagenta;
  flavcolor[3] = kRed-7;
  flavcolor[4] = kRed+1;
  flavcolor[5] = kRed-2;
  flavcolor[6] = kRed+0;
  flavcolor[7] = kRed+3;
  flavcolor[8] = kGreen+2;
  flavcolor[9] = kBlue;



  std::vector<TString> histoname5;
  if( includeFlavPlots_ ){
    histoname5.push_back("h_jet_csv");
    histoname5.push_back("h_jet_csv_wgtCSV");
    //histoname5.push_back("h_jet_cmva");
    //histoname5.push_back("h_jet_cmva_wgtCMVA");
    // histoname5.push_back("h_jet_csv_wgtCSV2");
    // histoname5.push_back("h_jet_csv_wgtCSV3");
    // histoname5.push_back("h_jet_csv_wgtCSV4");
    //histoname5.push_back("h_jet_csv_wgtCSV5");
  }

  std::vector<TString> jet_flavor_suffix;
  jet_flavor_suffix.push_back("lFlav");
  jet_flavor_suffix.push_back("cFlav");
  jet_flavor_suffix.push_back("bFlav");

  int NumJetFlav = int(jet_flavor_suffix.size());

  Color_t jetflavcolor[NumFlavComb];
  // jetflavcolor[0] = kAzure+2;
  // jetflavcolor[1] = kGreen+2;
  // jetflavcolor[2] = kRed+1;
  jetflavcolor[0] = kBlue;
  jetflavcolor[1] = kGreen+2;
  jetflavcolor[2] = kRed;


  double lumi_err = 0.046;

 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TGaxis::SetMaxDigits(4);

  TString lumiinfo = "2.67 fb^{-1} (13 TeV)";
  TLatex LumiInfoLatex(0.70, 0.91, lumiinfo);
  LumiInfoLatex.SetNDC(); LumiInfoLatex.SetTextFont(42);
  LumiInfoLatex.SetTextSize(0.04);

  //TString cmsinfo =   "CMS Preliminary";
  TString cmsinfo =   "CMS";
  TLatex CMSInfoLatex(0.185, 0.91, cmsinfo);
  CMSInfoLatex.SetNDC(); CMSInfoLatex.SetTextFont(42);
  CMSInfoLatex.SetTextFont(61);
  CMSInfoLatex.SetTextSize(0.055); //SBOUTLE

  TString publishinfo =   "Preliminary"; //DPUIGH
  TLatex PublishInfoLatex(0.285, 0.91, publishinfo); //SBOUTLE
  PublishInfoLatex.SetNDC();
  PublishInfoLatex.SetTextFont(52);
  PublishInfoLatex.SetTextSize(0.045); //SBOUTLE


  
  TString lumiinfo2 = "2.6 fb^{-1} (13 TeV, 25 ns)";
  TLatex LumiInfoLatex2(0.63, 0.91, lumiinfo2);
  LumiInfoLatex2.SetNDC(); LumiInfoLatex2.SetTextFont(42);
  LumiInfoLatex2.SetTextSize(0.04);

  //TString cmsinfo =   "CMS Preliminary";
  //TString cmsinfo =   "CMS";
  TLatex CMSInfoLatex2(0.155, 0.83, cmsinfo);
  CMSInfoLatex2.SetNDC(); CMSInfoLatex2.SetTextFont(42);
  CMSInfoLatex2.SetTextFont(61);
  CMSInfoLatex2.SetTextSize(0.055); //SBOUTLE

  //TString publishinfo =   "Preliminary"; //DPUIGH
  TLatex PublishInfoLatex2(0.155, 0.77, publishinfo); //SBOUTLE
  PublishInfoLatex2.SetNDC();
  PublishInfoLatex2.SetTextFont(52);
  PublishInfoLatex2.SetTextSize(0.045); //SBOUTLE



  double scale_ttH = 30.;


  TString plotname;


  TCanvas* myC1 = new TCanvas("myC1", "myC1", 600,700);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  Float_t small = 1.e-5;
  myC1->Divide(1,2,small,small);
  const float padding=1e-5; const float ydivide=0.3;
  myC1->GetPad(1)->SetPad( padding, ydivide + padding , 1-padding, 1-padding);
  myC1->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
  // myC1->GetPad(1)->SetLeftMargin(.11);
  // myC1->GetPad(2)->SetLeftMargin(.11);
  myC1->GetPad(1)->SetLeftMargin(.12);
  myC1->GetPad(2)->SetLeftMargin(.12);
  myC1->GetPad(1)->SetRightMargin(.05);
  myC1->GetPad(2)->SetRightMargin(.05);
  myC1->GetPad(1)->SetBottomMargin(.3);
  myC1->GetPad(2)->SetBottomMargin(.3);
  myC1->GetPad(1)->Modified();
  myC1->GetPad(2)->Modified();
  myC1->cd(1);
  gPad->SetBottomMargin(small);
  gPad->Modified();



  //histoname2.clear();
  for( int i=0; i<int(histoname2.size()); i++ ){

    TString temp = histoname2[i];
    
    TString temp_mh = temp;
    temp_mh.ReplaceAll("h_","");

    // //TLegend *legend = new TLegend(0.1,0.91,0.9,0.99);
    // TLegend *legend = new TLegend(0.13,0.83,0.88,0.89);

    // legend->SetFillColor(kWhite);
    // legend->SetLineColor(kWhite);
    // legend->SetShadowColor(kWhite);
    // legend->SetTextFont(42);
    // legend->SetTextSize(0.04);

    // legend->SetNColumns(4);


    // //// Original vertical legend
    // TLegend *legend = new TLegend(0.76,0.50,0.84,0.89);

    // legend->SetFillColor(kWhite);
    // legend->SetLineColor(kWhite);
    // legend->SetShadowColor(kWhite);
    // legend->SetTextFont(42);
    // legend->SetTextSize(0.035);

    //// horizontal legend
    TLegend *legend = new TLegend(0.16,0.75,0.89,0.89);

    legend->SetFillColor(kWhite);
    legend->SetLineColor(kWhite);
    legend->SetShadowColor(kWhite);
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);

    legend->SetNColumns(5);


    TLegend *legend_left = new TLegend(0.15,0.50,0.32,0.89);

    legend_left->SetFillColor(kWhite);
    legend_left->SetLineColor(kWhite);
    legend_left->SetShadowColor(kWhite);
    legend_left->SetTextFont(42);
    legend_left->SetTextSize(0.035);



    int rebin = 1;

    // if( temp.Contains("_mass") ) rebin = 2;
    // if( temp.Contains("_L1HTT") ) rebin = 2;
    // if( temp.Contains("_met") ) rebin = 5;
    // if( temp.Contains("_jet_") && temp.Contains("_pt") ) rebin = 10;
    // if( temp.Contains("_jet_") && temp.Contains("_eta") ) rebin = 4;
    // if( temp.Contains("_jet_") && temp.Contains("_phi") ) rebin = 4;
    if( temp.Contains("_csv") )       rebin = 3;
    if( temp.Contains("_relIso") )    rebin = 5;
    if( temp.Contains("_jet_puMVA") ) rebin = 5;
    if( temp.Contains("electron_d0") ) rebin = 2;
    if( temp.Contains("electron_dZ") ) rebin = 2;
    if( temp.Contains("muon_d0") ) rebin = 5;
    if( temp.Contains("muon_dZ") ) rebin = 5;

    // if( temp.Contains("_diele_mass_closestZmass") ) rebin = 5;

    if( temp.Contains("_HT") && !temp.Contains("4j") ) rebin = 2;

    double total_data = 0;
    double total_bkg  = 0;

    if( renormMC_  ){
      for( int iSample=0; iSample<NumSamples-2; iSample++ ){
	TH1D* h_delete_tmp = (TH1D*)file[iSample]->Get(temp)->Clone(Form("h_delete_tmp_%d_%s",iSample,temp.Data()));
	double total = h_delete_tmp->Integral();
	if( iSample==0 ) total_data = total;
	else             total_bkg += total;
      }
    }


    TH1D* hist_sys[NumSamples][NumSysCat];
    for( int iSample=0; iSample<NumSamples; iSample++ ){
      for( int iSys=0; iSys<NumSysCat; iSys++ ){

	TString temp_sys = temp + sys_cat_labels[iSys];
	if( temp.Contains("_diLepMass") ) temp_sys = temp;

	hist_sys[iSample][iSys] = (TH1D*)file[iSample]->Get(temp_sys)->Clone(Form("%s_std",temp_sys.Data()));
	hist_sys[iSample][iSys]->Rebin(rebin);

	if( renormMC_ && total_bkg>0 && total_data>0 && (iSample>=0 && iSample<NumSamples-2) ) hist_sys[iSample][iSys]->Scale( total_data/total_bkg );
	if( rescaleMC_>0 && !renormMC_ && iSample>0 ) hist_sys[iSample][iSys]->Scale( rescaleMC_ );
	//if( iSample>=1 && iSample<=5 ) hist_sys[iSample][iSys]->Scale( 2537./2430 );
      }
    }

    TH1D* hist[NumSamples];
    TH1D* hist_sum = NULL;
    bool firstFill = true;
    for( int iSample=0; iSample<NumSamples; iSample++ ){
      hist[iSample] = (TH1D*)hist_sys[iSample][0]->Clone(Form("%s_%d_std",temp.Data(),iSample));

      if( iSample>0 ){
	hist[iSample]->SetLineColor(color[iSample]);
	if( iSample<NumSamples-2 ) hist[iSample]->SetFillColor(color[iSample]);
	else hist[iSample]->SetLineWidth(3);

	if( firstFill && iSample>0 ){
	  firstFill = false;
	  hist_sum = (TH1D*)hist[iSample]->Clone("sum_"+temp);
	}
	else if( iSample<NumSamples-2 ){
	  hist_sum->Add( hist[iSample] );
	}
      }
      else {
	hist[iSample]->SetMarkerStyle(20);
	//hist[iSample]->SetLineWidth(20);
      }

      if( iSample==NumSamples-1 || iSample==NumSamples-2 ) hist[iSample]->Scale(scale_ttH);

      if( iSample==0 )                legend->AddEntry(hist[iSample],histLabels[iSample],"pe1");
      else if( iSample<NumSamples-2 ) legend->AddEntry(hist[iSample],histLabels[iSample],"f");
      else if( iSample==NumSamples-1 ) legend->AddEntry(hist[iSample],histLabels[iSample]+Form(" (x%d)",int(scale_ttH+0.0001)),"l");

      if( iSample==0 )                legend_left->AddEntry(hist[iSample],histLabels[iSample],"pe1");
      else if( iSample<NumSamples-2 ) legend_left->AddEntry(hist[iSample],histLabels[iSample],"f");
      else if( iSample==NumSamples-1 ) legend_left->AddEntry(hist[iSample],histLabels[iSample]+Form(" (x%d)",int(scale_ttH+0.0001)),"l");
    }// end loop over samples


    double data_integral = hist[0]->Integral();
    double bkg_integral = hist_sum->Integral();
    double integral_ratio = data_integral / bkg_integral;

    if( temp.Contains("numPV") || temp.Contains("_PV_") ){
      hist_sum->Scale(integral_ratio);
      for( int iSample=1; iSample<NumSamples; iSample++ ){
	hist[iSample]->Scale(integral_ratio);
	for( int iSys=0; iSys<NumSysCat; iSys++ ){
	  hist_sys[iSample][iSys]->Scale(integral_ratio);
	}
      }
    }


    int nbins = hist[0]->GetNbinsX();

    double xmin = hist[0]->GetBinLowEdge(1);
    double xmax = hist[0]->GetBinLowEdge(nbins) + hist[0]->GetBinWidth(nbins);



    if( false && (temp=="h_numBtag_wgtCSV_4j" || temp=="h_numBtag_wgtCSV2_4j") ){
      if( temp.Contains("h_numBtag_wgtCSV_4j") ) std::cout << "h_numBtag_wgtCSV_4j" << std::endl;
      else if( temp.Contains("h_numBtag_wgtCSV2_4j") ) std::cout << "h_numBtag_wgtCSV2_4j" << std::endl;

      std::cout << "   Data " << std::endl;
      for( int iBin=0; iBin<nbins; iBin++ ){
	printf("     iBin = %d: bin content = %.1f \n", iBin, hist[0]->GetBinContent(iBin+1));
      }

      std::cout << "   Sum Bkg " << std::endl;
      for( int iBin=0; iBin<nbins; iBin++ ){
	printf("     iBin = %d: bin content = %.1f \n", iBin, hist_sum->GetBinContent(iBin+1));
      }
    }


    THStack *hs = new THStack("hs","");
    for( int iSample=NumSamples-1; iSample>-1; iSample-- ){
      if( iSample>0 && iSample<NumSamples-2 ) hs->Add(hist[iSample]);
    }

    TH1D* h_sum_bkg_sys[NumSysCat];
    for( int iSys=0; iSys<NumSysCat; iSys++ ){
      bool firstSample = true;
      for( int iSample=1; iSample<NumSamples-2; iSample++ ){
	if( firstSample ){
	  firstSample = false;
	  h_sum_bkg_sys[iSys] = (TH1D*)hist_sys[iSample][iSys]->Clone(Form("%s_sys_%d",temp.Data(),iSys));
	}
	else h_sum_bkg_sys[iSys]->Add(hist_sys[iSample][iSys]);
      }
    }


    TH1D* h_bkg_err_1sig = new TH1D("h_bkg_err_1sig_"+temp,"", nbins, hist[0]->GetBinLowEdge(1), hist[0]->GetBinLowEdge(nbins) + hist[0]->GetBinWidth(nbins) );
    h_bkg_err_1sig->SetFillColor(kGreen);
    double sum_mc_err_tot_2 = 0.;
    for( int iBin=0; iBin<nbins; iBin++ ){

      double nom = h_sum_bkg_sys[0]->GetBinContent(iBin+1);
      double mcstat = hist_sum->GetBinError(iBin+1);

      sum_mc_err_tot_2 += mcstat*mcstat;

      double up2 = 0, down2 = 0;
      double previous_diff = 0;
      for( int iSys=1; iSys<NumSysCat; iSys++ ){
	if( temp.Contains("h_numPVs_") && !sys_cat_labels[iSys].Contains("PU") ) continue;
	double diff = h_sum_bkg_sys[iSys]->GetBinContent(iBin+1) - nom;
	double diff_2 = diff*diff;
	if( iBin!=1 && iBin%2==0 && (diff*previous_diff)>0 ){
	  if( fabs(diff)>fabs(previous_diff) ) diff_2 = diff*diff - previous_diff*previous_diff;
	} 
	if( diff>0 ) up2 += diff_2;
	else if( diff<0. ) down2 += diff_2;
	previous_diff = diff;

	if( temp.Contains("h_numBtag_wgtCSV4_2l") && !temp.Contains("h_numBtag_wgtCSV_4j") && !temp.Contains("h_numBtag_wgtCSV_4j2t") ){
	  printf(" Sys = %20s, iBin = %d, nominal = %10.3f, diff = %+10.3f, abs(diff)/nominal = %.3f, abs(up2)/nom = %.3f, abs(down2)/nom = %.3f \n",
		 sys_cat_labels[iSys].Data(), iBin, nom, diff, fabs(diff)/nom, sqrt(up2)/nom, sqrt(down2)/nom );
	}
      }

      up2 += mcstat*mcstat;
      down2 += mcstat*mcstat;

      double lumi_err2 = lumi_err * lumi_err * nom * nom;

      up2 += lumi_err2;
      down2 += lumi_err2;


      double up_err   = nom + sqrt(up2);
      double down_err = nom - sqrt(down2);

      double ave = 0.5 * ( up_err + down_err );
      h_bkg_err_1sig->SetBinContent(iBin+1,ave);
      h_bkg_err_1sig->SetBinError(iBin+1,up_err-ave);
    }



    // double sum_up2 = 0, sum_down2 = 0;
    // double sum_previous_diff = 0.;
    // for( int b=1; b<NumSysCat; b++ ){
    //   double diff = h_sum_bkg_sys[b]->Integral(first_bin,last_bin) - h_sum_bkg_sys[0]->Integral(first_bin,last_bin);
    //   double diff_2 = diff*diff;
    //   if( b!=1 && b%2==0 && (diff*sum_previous_diff)>0 ){
    // 	if( fabs(diff)>fabs(sum_previous_diff) ) diff_2 = diff*diff - sum_previous_diff*sum_previous_diff;
    //   } 
    //   if( diff>0 ) sum_up2 += diff_2;
    //   else if( diff<0. ) sum_down2 += diff_2;
    //   sum_previous_diff = diff;
    // }

    // double sum_sys_err_up   = sqrt( sum_up2 );
    // double sum_sys_err_down = sqrt( sum_down2 );
    // double sum_sys_err = ( sum_sys_err_up>sum_sys_err_down ) ? sum_sys_err_up : sum_sys_err_down;



    double ratioMax = 1.6;
    double ratioMin = 0.5;
    // double ratioMax = 2.3;
    // double ratioMin = 0.0;

    TH1D* myRatio = new TH1D("ratio", "", nbins, xmin, xmax );
    TH1D* myRatio_1sig = new TH1D("ratio_1sig_"+temp, "", nbins, xmin, xmax );

    myRatio->SetStats(0);
    myRatio->Sumw2();
    myRatio->SetLineColor(kBlack);
    myRatio->SetMarkerColor(kBlack);
    myRatio->Divide(hist[0],hist_sum);



    for( int iBin=0; iBin<nbins; iBin++ ){
      double bkg  = h_bkg_err_1sig->GetBinContent(iBin+1);
      double bkg_1sig  = h_bkg_err_1sig->GetBinError(iBin+1);
      double data = hist[0]->GetBinContent(iBin+1);
      double ratio = ( bkg>0. ) ? data/bkg : 0.;
      double ratio_err = ( bkg>0. ) ? sqrt(data)/bkg : 0.;

      double bkg_noshift  = hist_sum->GetBinContent(iBin+1);
      if( bkg_noshift>0. ) ratio = data/bkg_noshift;
      if( bkg_noshift>0. ) ratio_err = ( bkg>0. ) ? sqrt(data)/bkg_noshift : 0.;

      //if( fabs(bkg-bkg_noshift)>0.001 ) std::cout << " hist = " << temp << ", bin = " << bin << ", bkg = " << bkg << ", bkg_noshift


      double up_err = bkg + bkg_1sig;
      double down_err = bkg - bkg_1sig;

      if( bkg>0. && bkg_noshift>0. ){
	myRatio->SetBinContent(iBin+1,ratio);
	myRatio->SetBinError(iBin+1,ratio_err);

	up_err *= 1./bkg_noshift;
	down_err *= 1./bkg_noshift;

	double new_ave = 0.5 * ( up_err + down_err );

	myRatio_1sig->SetBinContent(iBin+1,new_ave);
	myRatio_1sig->SetBinError(iBin+1,up_err - new_ave);

	//std::cout << " hist = " << temp << ", bin = " << bin << ", new_ave = " << new_ave << ", up_err = " << up_err << ", down_err = " << down_err << std::endl;
      }

      if( (ratio>ratioMax) && ((ratio - ratio_err)<ratioMax) ){
	double minner = ratio - ratio_err;
	myRatio->SetBinContent(iBin+1,ratioMax-0.0001);
	myRatio->SetBinError(iBin+1,ratioMax-0.0001-minner);
      }

    }






    myRatio->SetMinimum(ratioMin);
    myRatio->SetMaximum(ratioMax);
    // double ratioMax = 2.3;
    // double ratioMin = 0.0;
    //myRatio->GetYaxis()->SetNdivisions(50000+404);
    // double ratioMax = 1.6;
    // double ratioMin = 0.5;
    myRatio->GetYaxis()->SetNdivisions(50000+204);
    myRatio->GetYaxis()->SetLabelSize(0.1); //make y label bigger
    myRatio->GetXaxis()->SetLabelSize(0.12); //make y label bigger
    myRatio->GetXaxis()->SetTitleOffset(1.1);
    myRatio->GetXaxis()->SetTitle(hist[0]->GetXaxis()->GetTitle()); //make y label bigger
    if( temp.Contains("h_num") ){
      myRatio->GetXaxis()->SetLabelSize(0.16);
      myRatio->GetXaxis()->SetLabelOffset(0.034);
    }
    else {
      myRatio->GetXaxis()->SetLabelSize(0.12);
      myRatio->GetXaxis()->SetLabelOffset(0.04);
    }
    myRatio->GetXaxis()->SetTitleSize(0.12);
    myRatio->GetYaxis()->SetTitle("Data/Bkg");
    myRatio->GetYaxis()->SetTitleSize(0.09);
    myRatio->GetYaxis()->SetTitleOffset(.55);
    myC1->cd(2);
    gPad->SetTopMargin(small);
    gPad->SetTickx();
    gPad->Modified();

    myRatio->GetYaxis()->CenterTitle(kTRUE);


    if( temp.Contains("_category_yield") || temp.Contains("_selection") ){
      for( int iBin=0; iBin<nbins; iBin++ ) myRatio->GetXaxis()->SetBinLabel(iBin+1,hist[0]->GetXaxis()->GetBinLabel(iBin+1));
    }

    if( temp.Contains("h_num") && !temp.Contains("numPV") ){
      for( int iBin=0; iBin<nbins; iBin++ ) myRatio->GetXaxis()->SetBinLabel(iBin+1,Form("%d",int(0.001 + hist[0]->GetXaxis()->GetBinCenter(iBin+1))));
      for( int iBin=0; iBin<nbins; iBin++ ) hist[0]->GetXaxis()->SetBinLabel(iBin+1,Form("%d",int(0.001 + hist[0]->GetXaxis()->GetBinCenter(iBin+1))));
    }


    hist[0]->SetStats(0);

    hist[0]->GetYaxis()->SetTitleOffset(1.2);
    hist[0]->GetYaxis()->SetTitleSize(0.05);

    hist[0]->GetYaxis()->SetTitle("Number of Events");


    int max_bin_data = hist[0]->GetMaximumBin();
    double max_data = hist[0]->GetBinContent(max_bin_data) + hist[0]->GetBinError(max_bin_data);

    int max_bin_mc = h_bkg_err_1sig->GetMaximumBin();
    double max_mc = h_bkg_err_1sig->GetBinContent(max_bin_mc) + h_bkg_err_1sig->GetBinError(max_bin_mc);

    double max_content = std::max(max_data, max_mc);

    hist[0]->GetYaxis()->SetRangeUser(0.,1.2 * max_content);


    bool rerange = false;

    if( (temp.Contains("_electron_dZ")) ){ rerange=true; xmin = 0.0; xmax = 0.05; }
    if( (temp.Contains("_electron_d0")) ){ rerange=true; xmin = 0.0; xmax = 0.05; }

    if( (temp.Contains("_jet_csv")) ){ rerange=true; xmin = -0.04; xmax = 1.028; }
    if( (temp.Contains("_HT")) ){ rerange=true; xmin = 100+0.0001; xmax = 1400-0.001; }
    if( (temp.Contains("_mht_pt")) ){ rerange=true; xmin = 0+0.0001; xmax = 200-0.001; }
    if( (temp.Contains("_met_pt")) ){ rerange=true; xmin = 0+0.0001; xmax = 400-0.001; }
    //if( (temp.Contains("_category")) ){ rerange=true; xmin = 0.001; xmax = 6.99; }

    if( temp.Contains("h_numBtag_wgtCSV_4j2t") ){ rerange=true; xmin = 1.5+0.0001; xmax = 4.5-0.001; }

    if( temp.Contains("h_numBtag") && temp.Contains("4j") && !temp.Contains("4j2t") && !temp.Contains("met") ){ rerange=true; xmin = 0.5+0.0001; xmax = 4.5-0.001; }
    if( temp.Contains("h_numBtag") && temp.Contains("4j") && !temp.Contains("4j2t") && temp.Contains("met")  ){ rerange=true; xmin = -0.5+0.0001; xmax = 4.5-0.001; }
    if( temp.Contains("h_numBtag") && temp.Contains("2l") ){ rerange=true; xmin = -0.5+0.0001; xmax = 2.5-0.001; }

    if( rerange ){
      hist[0]->GetXaxis()->SetRangeUser(xmin,xmax);
      myRatio->GetXaxis()->SetRangeUser(xmin,xmax);
    }

    // if( temp.Contains("_mass") ){
    //   hist[0]->GetXaxis()->SetRangeUser(40.,140.);
    //   myRatio->GetXaxis()->SetRangeUser(40.,140.);
    // }

    // if( temp.Contains("_met") ){
    //   hist[0]->GetXaxis()->SetRangeUser(0.,150.);
    //   myRatio->GetXaxis()->SetRangeUser(0.,150.);
    // }
    // // if( temp.Contains("_eta") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);
    // // if( temp.Contains("_numGenPVs") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);
    // // if( temp.Contains("_numJets") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);

    h_bkg_err_1sig->SetFillStyle(3654);
    h_bkg_err_1sig->SetFillColor(kBlack);
    h_bkg_err_1sig->SetMarkerSize(0);


    myRatio_1sig->SetMarkerColor(kGreen);
    myRatio_1sig->SetFillColor(kGreen);
    myRatio_1sig->SetMarkerSize(0); // geoff


    TLine* myLine;
    myLine = new TLine(xmin, 1, xmax, 1);


    // Plot
    myC1->cd(1);
    hist[0]->Draw("pe1");
    hs->Draw("histsame");
    h_bkg_err_1sig->Draw("e2same");
    //hist[NumSamples-2]->Draw("histsame");
    hist[NumSamples-1]->Draw("histsame");
    hist[0]->Draw("pe1same");

    // if( temp.Contains("_jet_csv") ) legend_left->Draw();
    // else                            legend->Draw();
    legend->Draw();

    LumiInfoLatex.Draw();
    CMSInfoLatex.Draw();
    PublishInfoLatex.Draw();

    myC1->cd(2);
    myRatio->SetLineWidth(2);
    myRatio->Draw("pe1");
    myRatio_1sig->Draw("e2same");
    myRatio->Draw("pe1same");

    myLine->Draw();

    myC1->GetPad(1)->SetLogy(0);

    myC1->GetPad(1)->RedrawAxis();
    myC1->GetPad(2)->RedrawAxis();

    plotname = dirprefix + temp_mh + "_data2mc_lin.png";
    myC1->Print(plotname);

    plotname = dirprefix + temp_mh + "_data2mc_lin.pdf";
    if( printPDF_ ) myC1->Print(plotname);


    if( temp.Contains("_category_yield") || temp.Contains("h_num") ){
      // log
      hist[0]->GetYaxis()->SetRangeUser(0.4,12 * max_content);

      // if including the zero bin
      if( temp.Contains("numBtag") && temp.Contains("4j") && temp.Contains("met") ) hist[0]->GetYaxis()->SetRangeUser(0.4,40 * max_content);

      // if not including the zero bin
      if( temp.Contains("numBtag") && temp.Contains("4j") && !temp.Contains("met") ) hist[0]->GetYaxis()->SetRangeUser(0.4,12 * max_content);

      if( temp.Contains("_selection") ){
	hist[0]->GetYaxis()->SetRangeUser(1000,12 * max_content);
      }

      myC1->cd(1);
      //gPad->SetLogy(1);

      myC1->GetPad(1)->SetLogy(1);

      myC1->cd(1);
      hist[0]->Draw("pe1");
      hs->Draw("histsame");
      h_bkg_err_1sig->Draw("e2same");
      //hist[NumSamples-2]->Draw("histsame");
      hist[NumSamples-1]->Draw("histsame");
      hist[0]->Draw("pe1same");

      if( temp.Contains("_jet_csv") ) legend_left->Draw();
      else                            legend->Draw();

      LumiInfoLatex.Draw();
      CMSInfoLatex.Draw();
      PublishInfoLatex.Draw();

      myC1->cd(2);
      myRatio->SetLineWidth(2);
      myRatio->Draw("pe1");
      myRatio_1sig->Draw("e2same");
      myRatio->Draw("pe1same");
      myLine->Draw();

      myC1->GetPad(1)->RedrawAxis();
      myC1->GetPad(2)->RedrawAxis();

      plotname = dirprefix + temp_mh + "_data2mc_log.png";
      myC1->Print(plotname);

      plotname = dirprefix + temp_mh + "_data2mc_log.pdf";
      if( printPDF_ ) myC1->Print(plotname);

      //gPad->SetLogy(0);
    }


    delete myRatio;
    delete myRatio_1sig;
    delete h_bkg_err_1sig;
    delete myLine;
    delete legend;
  } // end loop on hists




  for( int i=0; i<int(histoname3.size()); i++ ){

    TString temp = histoname3[i];
    
    TString temp_mh = temp;
    temp_mh.ReplaceAll("h_","");

    TH1D* hist[NumSamples];
    TH1D* hist_sum = NULL;
    bool firstFill = true;
    for( int iSample=0; iSample<NumSamples; iSample++ ){
      hist[iSample] = (TH1D*)file[iSample]->Get(temp)->Clone(Form("%s_cat",temp.Data()));

      if( iSample>0 ){
	if( firstFill && iSample>0 ){
	  firstFill = false;
	  hist_sum = (TH1D*)hist[iSample]->Clone("sum_"+temp);
	}
	else if( iSample<NumSamples-2 ){
	  hist_sum->Add( hist[iSample] );
	}
      }
    }// end loop over samples

    int NumBins = hist[0]->GetNbinsX();

    double xmin = hist[0]->GetBinLowEdge(1);
    double xmax = hist[0]->GetBinLowEdge(NumBins) + hist[0]->GetBinWidth(NumBins);

    TH1D* myRatio = new TH1D("ratio_"+temp, "", NumBins, xmin, xmax );
    myRatio->Divide(hist[0],hist_sum);


    std::cout << " ============================================= " << std::endl;

    std::cout << " \t " << temp_mh << std::endl;

    std::cout << "\\begin{tabular}{|l|c|c|c|c|c|c|c|c|c|c|} \\hline" << std::endl; 
    std::cout << "\t";
    for( int iBin=0; iBin<NumBins; iBin++ ){
      std::cout << "&\t" << hist[0]->GetXaxis()->GetBinLabel(iBin+1);
    }
    std::cout << "\\\\ \\hline \\hline" << std::endl;

    printf("%8s", histLabels[NumSamples-1].Data());
    //std::cout << histLabels[NumSamples-1];
    for( int iBin=0; iBin<NumBins; iBin++ ){
      printf(" & %6.1f $\\pm$ %3.1f", hist[NumSamples-1]->GetBinContent(iBin+1), hist[NumSamples-1]->GetBinError(iBin+1));
    }
    std::cout << "\\\\" << std::endl;    

    printf("%8s", histLabels[NumSamples-2].Data());
    //std::cout << histLabels[NumSamples-2];
    for( int iBin=0; iBin<NumBins; iBin++ ){
      printf(" & %6.1f $\\pm$ %3.1f", hist[NumSamples-2]->GetBinContent(iBin+1), hist[NumSamples-2]->GetBinError(iBin+1));
    }
    std::cout << " " << std::endl;  
    std::cout << "\\\\ \\hline" << std::endl;

    for( int iSample=1; iSample<NumSamples-2; iSample++ ){
      //std::cout << histLabels[iSample];
      printf("%8s", histLabels[iSample].Data());
      for( int iBin=0; iBin<NumBins; iBin++ ){
	printf(" & %6.1f $\\pm$ %3.1f", hist[iSample]->GetBinContent(iBin+1), hist[iSample]->GetBinError(iBin+1));
      }
      std::cout << " \\\\" << std::endl; 
    }
    std::cout << "\\hline" << std::endl;

    //std::cout << "Data";
    printf("%8s", "Data");
    for( int iBin=0; iBin<NumBins; iBin++ ){
      printf(" & %6.0f ", hist[0]->GetBinContent(iBin+1));
    }
    std::cout << " \\\\ \\hline" << std::endl;  

    printf("%8s", "TotBkg");
    for( int iBin=0; iBin<NumBins; iBin++ ){
      printf(" & %6.1f $\\pm$ %3.1f", hist_sum->GetBinContent(iBin+1), hist_sum->GetBinError(iBin+1));
    }
    std::cout << " \\\\ \\hline" << std::endl;  


    //std::cout << "Data/Bkg";
    printf("%8s", "Data/Bkg");
    for( int iBin=0; iBin<NumBins; iBin++ ){
      printf(" & %6.2f $\\pm$ %3.2f", myRatio->GetBinContent(iBin+1), myRatio->GetBinError(iBin+1));
    }
    std::cout << " \\\\ \\hline" << std::endl;  

    printf("%8s", "S/Sqrt(B)");
    for( int iBin=0; iBin<NumBins; iBin++ ){
      double tth = hist[NumSamples-1]->GetBinContent(iBin+1);
      double sumb = hist_sum->GetBinContent(iBin+1);
      double soversqrtb = ( sumb > 0. ) ? tth / sqrt(sumb) : 0.;
      printf(" & %6.3f", soversqrtb);
    }
    std::cout << " \\\\ \\hline" << std::endl;  
    
    std::cout << " \\end{tabular} " << std::endl;


    delete myRatio;
  }
  std::cout << " ============================================= " << std::endl;





  for( int i=0; i<int(histoname4.size()); i++ ){

    TString temp = histoname4[i];
    
    TString temp_mh = temp + "_jetFlavor";
    temp_mh.ReplaceAll("h_","");

    //// horizontal legend
    TLegend *legend = new TLegend(0.16,0.75,0.89,0.89);

    legend->SetFillColor(kWhite);
    legend->SetLineColor(kWhite);
    legend->SetShadowColor(kWhite);
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);

    legend->SetNColumns(5);

    int rebin = 1;


    double total_data = 0;
    double total_bkg  = 0;

    if( renormMC_ ){
      for( int iSample=0; iSample<NumSamples-2; iSample++ ){
	TH1D* h_delete_tmp = (TH1D*)file[iSample]->Get(temp)->Clone(Form("h_delete_tmp_%d_%s_jetFlavor",iSample,temp.Data()));
	double total = h_delete_tmp->Integral();
	if( iSample==0 ) total_data = total;
	else             total_bkg += total;
      }
    }

    TH1D* hist_sys[NumSamples][NumSysCat];
    for( int iSample=0; iSample<NumSamples; iSample++ ){
      for( int iSys=0; iSys<NumSysCat; iSys++ ){

	TString temp_sys = temp + sys_cat_labels[iSys];
	if( temp.Contains("_diLepMass") ) temp_sys = temp;

	hist_sys[iSample][iSys] = (TH1D*)file[iSample]->Get(temp_sys)->Clone(Form("%s_tmp_sys_jetFlav",temp_sys.Data()));
	hist_sys[iSample][iSys]->Rebin(rebin);

	if( renormMC_ && total_bkg>0 && total_data>0 && (iSample>=0 && iSample<NumSamples-2) ) hist_sys[iSample][iSys]->Scale( total_data/total_bkg );
	if( rescaleMC_>0 && !renormMC_ && iSample>0 ) hist_sys[iSample][iSys]->Scale( rescaleMC_ );
	//if( iSample>=1 && iSample<=5 ) hist_sys[iSample][iSys]->Scale( 2537./2430 );
      }
    }

    TH1D* hist[NumSamples];
    TH1D* hist_sum = NULL;
    bool firstFill = true;
    TH1D* hist_bkg_flav[NumSamples][NumFlavComb];
    for( int iSample=0; iSample<NumSamples; iSample++ ){
      hist[iSample] = (TH1D*)hist_sys[iSample][0]->Clone(Form("%s_%d_jetFlav",temp.Data(),iSample));

      if( iSample>0 && iSample<NumSamples-2 ){
	for( int iFlav=0; iFlav<NumFlavComb; iFlav++ ){
	  TString temp_flav = temp + "_" + flavor_suffix[iFlav];
	  hist_bkg_flav[iSample][iFlav] = (TH1D*)file[iSample]->Get(temp_flav)->Clone(Form("%s_%d_jetFlav",temp_flav.Data(),iSample));
	}
      }

      if( iSample>0 ){
	hist[iSample]->SetLineColor(color[iSample]);
	if( iSample<NumSamples-2 ) hist[iSample]->SetFillColor(color[iSample]);
	else hist[iSample]->SetLineWidth(3);

	if( firstFill && iSample>0 ){
	  firstFill = false;
	  hist_sum = (TH1D*)hist[iSample]->Clone("sum_"+temp+"_jetFlav");
	}
	else if( iSample<NumSamples-2 ){
	  hist_sum->Add( hist[iSample] );
	}
      }
      else {
	hist[iSample]->SetMarkerStyle(20);
	//hist[iSample]->SetLineWidth(20);
      }

      if( iSample==NumSamples-1 || iSample==NumSamples-2 ) hist[iSample]->Scale(scale_ttH);

      // if( iSample==0 )                legend->AddEntry(hist[iSample],histLabels[iSample],"pe1");
      // else if( iSample<NumSamples-2 ) legend->AddEntry(hist[iSample],histLabels[iSample],"f");
      // else if( iSample==NumSamples-1 ) legend->AddEntry(hist[iSample],histLabels[iSample]+Form(" (x%d)",int(scale_ttH+0.0001)),"l");

    }// end loop over samples


    double data_integral = hist[0]->Integral();
    double bkg_integral = hist_sum->Integral();
    double integral_ratio = data_integral / bkg_integral;

    if( temp.Contains("numPV") ){
      hist_sum->Scale(integral_ratio);
      for( int iSample=1; iSample<NumSamples; iSample++ ){
	hist[iSample]->Scale(integral_ratio);
	for( int iSys=0; iSys<NumSysCat; iSys++ ){
	  hist_sys[iSample][iSys]->Scale(integral_ratio);
	}
      }
    }


    int nbins = hist[0]->GetNbinsX();

    double xmin = hist[0]->GetBinLowEdge(1);
    double xmax = hist[0]->GetBinLowEdge(nbins) + hist[0]->GetBinWidth(nbins);



    if( false && (temp=="h_numBtag_wgtCSV_4j" || temp=="h_numBtag_wgtCSV2_4j") ){
      if( temp.Contains("h_numBtag_wgtCSV_4j") ) std::cout << "h_numBtag_wgtCSV_4j" << std::endl;
      else if( temp.Contains("h_numBtag_wgtCSV2_4j") ) std::cout << "h_numBtag_wgtCSV2_4j" << std::endl;

      std::cout << "   Data " << std::endl;
      for( int iBin=0; iBin<nbins; iBin++ ){
	printf("     iBin = %d: bin content = %.1f \n", iBin, hist[0]->GetBinContent(iBin+1));
      }

      std::cout << "   Sum Bkg " << std::endl;
      for( int iBin=0; iBin<nbins; iBin++ ){
	printf("     iBin = %d: bin content = %.1f \n", iBin, hist_sum->GetBinContent(iBin+1));
      }
    }

    legend->AddEntry(hist[0],histLabels[0],"pe1");

    TH1D* hist_bkg_sum_flav[NumFlavComb];
    for( int iFlav=0; iFlav<NumFlavComb; iFlav++ ){      
      TString temp_flav = temp + "_" + flavor_suffix[iFlav];
      bool firstFlavFill = true;
      for( int iSample=1; iSample<NumSamples-2; iSample++ ){
	if( firstFlavFill ){
	  firstFlavFill = false;
	  hist_bkg_sum_flav[iFlav] = (TH1D*)hist_bkg_flav[iSample][iFlav]->Clone(Form("%s_bkg_sum_jetFlav",temp_flav.Data()));
	}
	else {
	  hist_bkg_sum_flav[iFlav]->Add(hist_bkg_flav[iSample][iFlav]);
	}
      }

      hist_bkg_sum_flav[iFlav]->SetFillColor(flavcolor[iFlav]);

      legend->AddEntry(hist_bkg_sum_flav[iFlav],flavor_suffix[iFlav],"f");
    }
    //legend->AddEntry(hist[NumSamples-2],histLabels[NumSamples-2]+Form(" (x%d)",int(scale_ttH+0.0001)),"l");



    THStack *hs = new THStack("hs","");
    for( int iFlav=0; iFlav<NumFlavComb; iFlav++ ){      
      hs->Add(hist_bkg_sum_flav[iFlav]);
    }


    TH1D* h_sum_bkg_sys[NumSysCat];
    for( int iSys=0; iSys<NumSysCat; iSys++ ){
      bool firstSample = true;
      for( int iSample=1; iSample<NumSamples-2; iSample++ ){
	if( firstSample ){
	  firstSample = false;
	  h_sum_bkg_sys[iSys] = (TH1D*)hist_sys[iSample][iSys]->Clone(Form("%s_sys_%d_jetFlav",temp.Data(),iSys));
	}
	else h_sum_bkg_sys[iSys]->Add(hist_sys[iSample][iSys]);
      }
    }


    TH1D* h_bkg_err_1sig = new TH1D("h_bkg_err_1sig_"+temp,"", nbins, hist[0]->GetBinLowEdge(1), hist[0]->GetBinLowEdge(nbins) + hist[0]->GetBinWidth(nbins) );
    h_bkg_err_1sig->SetFillColor(kGreen);
    double sum_mc_err_tot_2 = 0.;
    for( int iBin=0; iBin<nbins; iBin++ ){

      double nom = h_sum_bkg_sys[0]->GetBinContent(iBin+1);
      double mcstat = hist_sum->GetBinError(iBin+1);

      sum_mc_err_tot_2 += mcstat*mcstat;

      double up2 = 0, down2 = 0;
      double previous_diff = 0;
      for( int iSys=1; iSys<NumSysCat; iSys++ ){
	if( temp.Contains("h_numPVs_") && !sys_cat_labels[iSys].Contains("PU") ) continue;
	double diff = h_sum_bkg_sys[iSys]->GetBinContent(iBin+1) - nom;
	double diff_2 = diff*diff;
	if( iBin!=1 && iBin%2==0 && (diff*previous_diff)>0 ){
	  if( fabs(diff)>fabs(previous_diff) ) diff_2 = diff*diff - previous_diff*previous_diff;
	} 
	if( diff>0 ) up2 += diff_2;
	else if( diff<0. ) down2 += diff_2;
	previous_diff = diff;

	if( temp.Contains("h_numBtag_wgtCSV4_2l") && !temp.Contains("h_numBtag_wgtCSV_4j") && !temp.Contains("h_numBtag_wgtCSV_4j2t") ){
	  printf(" Sys = %20s, iBin = %d, nominal = %10.3f, diff = %+10.3f, abs(diff)/nominal = %.3f, abs(up2)/nom = %.3f, abs(down2)/nom = %.3f \n",
		 sys_cat_labels[iSys].Data(), iBin, nom, diff, fabs(diff)/nom, sqrt(up2)/nom, sqrt(down2)/nom );
	}
      }

      up2 += mcstat*mcstat;
      down2 += mcstat*mcstat;

      double lumi_err2 = lumi_err * lumi_err * nom * nom;

      up2 += lumi_err2;
      down2 += lumi_err2;


      double up_err   = nom + sqrt(up2);
      double down_err = nom - sqrt(down2);

      double ave = 0.5 * ( up_err + down_err );
      h_bkg_err_1sig->SetBinContent(iBin+1,ave);
      h_bkg_err_1sig->SetBinError(iBin+1,up_err-ave);
    }


    double ratioMax = 1.6;
    double ratioMin = 0.5;
    // double ratioMax = 2.3;
    // double ratioMin = 0.0;

    TH1D* myRatio = new TH1D("ratio", "", nbins, xmin, xmax );
    TH1D* myRatio_1sig = new TH1D("ratio_1sig_"+temp, "", nbins, xmin, xmax );

    myRatio->SetStats(0);
    myRatio->Sumw2();
    myRatio->SetLineColor(kBlack);
    myRatio->SetMarkerColor(kBlack);
    myRatio->Divide(hist[0],hist_sum);



    for( int iBin=0; iBin<nbins; iBin++ ){
      double bkg  = h_bkg_err_1sig->GetBinContent(iBin+1);
      double bkg_1sig  = h_bkg_err_1sig->GetBinError(iBin+1);
      double data = hist[0]->GetBinContent(iBin+1);
      double ratio = ( bkg>0. ) ? data/bkg : 0.;
      double ratio_err = ( bkg>0. ) ? sqrt(data)/bkg : 0.;

      double bkg_noshift  = hist_sum->GetBinContent(iBin+1);
      if( bkg_noshift>0. ) ratio = data/bkg_noshift;
      if( bkg_noshift>0. ) ratio_err = ( bkg>0. ) ? sqrt(data)/bkg_noshift : 0.;

      //if( fabs(bkg-bkg_noshift)>0.001 ) std::cout << " hist = " << temp << ", bin = " << bin << ", bkg = " << bkg << ", bkg_noshift


      double up_err = bkg + bkg_1sig;
      double down_err = bkg - bkg_1sig;

      if( bkg>0. && bkg_noshift>0. ){
	myRatio->SetBinContent(iBin+1,ratio);
	myRatio->SetBinError(iBin+1,ratio_err);

	up_err *= 1./bkg_noshift;
	down_err *= 1./bkg_noshift;

	double new_ave = 0.5 * ( up_err + down_err );

	myRatio_1sig->SetBinContent(iBin+1,new_ave);
	myRatio_1sig->SetBinError(iBin+1,up_err - new_ave);

	//std::cout << " hist = " << temp << ", bin = " << bin << ", new_ave = " << new_ave << ", up_err = " << up_err << ", down_err = " << down_err << std::endl;
      }

      if( (ratio>ratioMax) && ((ratio - ratio_err)<ratioMax) ){
	double minner = ratio - ratio_err;
	myRatio->SetBinContent(iBin+1,ratioMax-0.0001);
	myRatio->SetBinError(iBin+1,ratioMax-0.0001-minner);
      }

    }






    myRatio->SetMinimum(ratioMin);
    myRatio->SetMaximum(ratioMax);
    // double ratioMax = 2.3;
    // double ratioMin = 0.0;
    //myRatio->GetYaxis()->SetNdivisions(50000+404);
    // double ratioMax = 1.6;
    // double ratioMin = 0.5;
    myRatio->GetYaxis()->SetNdivisions(50000+204);
    myRatio->GetYaxis()->SetLabelSize(0.1); //make y label bigger
    myRatio->GetXaxis()->SetLabelSize(0.12); //make y label bigger
    myRatio->GetXaxis()->SetTitleOffset(1.1);
    myRatio->GetXaxis()->SetTitle(hist[0]->GetXaxis()->GetTitle()); //make y label bigger
    if( temp.Contains("h_num") ){
      myRatio->GetXaxis()->SetLabelSize(0.16);
      myRatio->GetXaxis()->SetLabelOffset(0.034);
    }
    else {
      myRatio->GetXaxis()->SetLabelSize(0.12);
      myRatio->GetXaxis()->SetLabelOffset(0.04);
    }
    myRatio->GetXaxis()->SetTitleSize(0.12);
    myRatio->GetYaxis()->SetTitle("Data/Bkg");
    myRatio->GetYaxis()->SetTitleSize(0.09);
    myRatio->GetYaxis()->SetTitleOffset(.55);
    myC1->cd(2);
    gPad->SetTopMargin(small);
    gPad->SetTickx();
    gPad->Modified();

    myRatio->GetYaxis()->CenterTitle(kTRUE);


    if( temp.Contains("_category_yield") || temp.Contains("_selection") ){
      for( int iBin=0; iBin<nbins; iBin++ ) myRatio->GetXaxis()->SetBinLabel(iBin+1,hist[0]->GetXaxis()->GetBinLabel(iBin+1));
    }

    if( temp.Contains("h_num") && !temp.Contains("numPV") ){
      for( int iBin=0; iBin<nbins; iBin++ ) myRatio->GetXaxis()->SetBinLabel(iBin+1,Form("%d",int(0.001 + hist[0]->GetXaxis()->GetBinCenter(iBin+1))));
      for( int iBin=0; iBin<nbins; iBin++ ) hist[0]->GetXaxis()->SetBinLabel(iBin+1,Form("%d",int(0.001 + hist[0]->GetXaxis()->GetBinCenter(iBin+1))));
    }


    hist[0]->SetStats(0);

    hist[0]->GetYaxis()->SetTitleOffset(1.2);
    hist[0]->GetYaxis()->SetTitleSize(0.05);

    hist[0]->GetYaxis()->SetTitle("Number of Events");


    int max_bin_data = hist[0]->GetMaximumBin();
    double max_data = hist[0]->GetBinContent(max_bin_data) + hist[0]->GetBinError(max_bin_data);

    int max_bin_mc = h_bkg_err_1sig->GetMaximumBin();
    double max_mc = h_bkg_err_1sig->GetBinContent(max_bin_mc) + h_bkg_err_1sig->GetBinError(max_bin_mc);

    double max_content = std::max(max_data, max_mc);

    hist[0]->GetYaxis()->SetRangeUser(0.,1.2 * max_content);
    //hist[0]->GetYaxis()->SetRangeUser(0.,1.2 * hist[0]->GetBinContent(2));


    bool rerange = false;

    if( (temp.Contains("_jet_csv")) ){ rerange=true; xmin = -0.04; xmax = 1.028; }
    if( (temp.Contains("_HT")) ){ rerange=true; xmin = 100+0.0001; xmax = 1400-0.001; }
    if( (temp.Contains("_mht_pt")) ){ rerange=true; xmin = 0+0.0001; xmax = 200-0.001; }
    if( (temp.Contains("_met_pt")) ){ rerange=true; xmin = 0+0.0001; xmax = 400-0.001; }
    //if( (temp.Contains("_category")) ){ rerange=true; xmin = 0.001; xmax = 6.99; }

    if( temp.Contains("h_numBtag_wgtCSV_4j2t") ){ rerange=true; xmin = 1.5+0.0001; xmax = 4.5-0.001; }

    if( temp.Contains("h_numBtag") && temp.Contains("4j") && !temp.Contains("4j2t") && !temp.Contains("met") ){ rerange=true; xmin = 0.5+0.0001; xmax = 4.5-0.001; }
    if( temp.Contains("h_numBtag") && temp.Contains("4j") && !temp.Contains("4j2t") && temp.Contains("met")  ){ rerange=true; xmin = -0.5+0.0001; xmax = 4.5-0.001; }
    if( temp.Contains("h_numBtag") && temp.Contains("2l") ){ rerange=true; xmin = -0.5+0.0001; xmax = 2.5-0.001; }

    if( rerange ){
      hist[0]->GetXaxis()->SetRangeUser(xmin,xmax);
      myRatio->GetXaxis()->SetRangeUser(xmin,xmax);
    }

    // if( temp.Contains("_mass") ){
    //   hist[0]->GetXaxis()->SetRangeUser(40.,140.);
    //   myRatio->GetXaxis()->SetRangeUser(40.,140.);
    // }

    // if( temp.Contains("_met") ){
    //   hist[0]->GetXaxis()->SetRangeUser(0.,150.);
    //   myRatio->GetXaxis()->SetRangeUser(0.,150.);
    // }
    // // if( temp.Contains("_eta") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);
    // // if( temp.Contains("_numGenPVs") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);
    // // if( temp.Contains("_numJets") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);

    h_bkg_err_1sig->SetFillStyle(3654);
    h_bkg_err_1sig->SetFillColor(kBlack);
    h_bkg_err_1sig->SetMarkerSize(0);


    myRatio_1sig->SetMarkerColor(kGreen);
    myRatio_1sig->SetFillColor(kGreen);
    myRatio_1sig->SetMarkerSize(0); // geoff


    TLine* myLine;
    myLine = new TLine(xmin, 1, xmax, 1);


    // Plot
    myC1->cd(1);
    hist[0]->Draw("pe1");
    hs->Draw("histsame");
    h_bkg_err_1sig->Draw("e2same");
    //hist[NumSamples-2]->Draw("histsame");
    //hist[NumSamples-1]->Draw("histsame");
    hist[0]->Draw("pe1same");

    legend->Draw();

    LumiInfoLatex.Draw();
    CMSInfoLatex.Draw();
    PublishInfoLatex.Draw();

    myC1->cd(2);
    myRatio->SetLineWidth(2);
    myRatio->Draw("pe1");
    myRatio_1sig->Draw("e2same");
    myRatio->Draw("pe1same");

    myLine->Draw();

    myC1->GetPad(1)->SetLogy(0);

    myC1->GetPad(1)->RedrawAxis();
    myC1->GetPad(2)->RedrawAxis();

    plotname = dirprefix + temp_mh + "_data2mc_lin.png";
    myC1->Print(plotname);

    plotname = dirprefix + temp_mh + "_data2mc_lin.pdf";
    if( printPDF_ ) myC1->Print(plotname);


    if( temp.Contains("_category_yield") || temp.Contains("h_num") ){
      // log
      hist[0]->GetYaxis()->SetRangeUser(0.4,40 * max_content);

      // if including the zero bin
      if( temp.Contains("numBtag") && temp.Contains("4j") && temp.Contains("met") ) hist[0]->GetYaxis()->SetRangeUser(0.4,40 * max_content);

      // if not including the zero bin
      if( temp.Contains("numBtag") && temp.Contains("4j") && !temp.Contains("met") ) hist[0]->GetYaxis()->SetRangeUser(0.4,12 * max_content);

      if( temp.Contains("_selection") ){
	hist[0]->GetYaxis()->SetRangeUser(1000,12 * max_content);
      }

      myC1->cd(1);
      //gPad->SetLogy(1);

      myC1->GetPad(1)->SetLogy(1);

      myC1->cd(1);
      hist[0]->Draw("pe1");
      hs->Draw("histsame");
      h_bkg_err_1sig->Draw("e2same");
      //hist[NumSamples-2]->Draw("histsame");
      //hist[NumSamples-1]->Draw("histsame");
      hist[0]->Draw("pe1same");

      legend->Draw();

      LumiInfoLatex.Draw();
      CMSInfoLatex.Draw();
      PublishInfoLatex.Draw();

      myC1->cd(2);
      myRatio->SetLineWidth(2);
      myRatio->Draw("pe1");
      myRatio_1sig->Draw("e2same");
      myRatio->Draw("pe1same");
      myLine->Draw();

      myC1->GetPad(1)->RedrawAxis();
      myC1->GetPad(2)->RedrawAxis();

      plotname = dirprefix + temp_mh + "_data2mc_log.png";
      myC1->Print(plotname);

      plotname = dirprefix + temp_mh + "_data2mc_log.pdf";
      if( printPDF_ ) myC1->Print(plotname);

      //gPad->SetLogy(0);
    }


    delete myRatio;
    delete myRatio_1sig;
    delete h_bkg_err_1sig;
    delete myLine;
    delete legend;
  } // end loop on hists histoname 4




  for( int i=0; i<int(histoname5.size()); i++ ){

    TString temper = histoname5[i];
    
    for( int iCat=0; iCat<NumCat; iCat++ ){

      TString temp = temper + "_" + cat_labels[iCat];

      TString temp_mh = temp + "_jetFlav";
      temp_mh.ReplaceAll("h_","");

      //// horizontal legend
      TLegend *legend = new TLegend(0.62,0.60,0.85,0.89);

      legend->SetFillColor(kWhite);
      legend->SetLineColor(kWhite);
      legend->SetShadowColor(kWhite);
      legend->SetTextFont(42);
      legend->SetTextSize(0.035);

      //legend->SetNColumns(2);

      int rebin = 1;

      if( temp.Contains("_csv") )       rebin = 3;
      if( temp.Contains("_csv") && (temp.Contains("j3t") || temp.Contains("j3t")) ) rebin = 6;
      if( temp.Contains("_cmva") )       rebin = 3;
      if( temp.Contains("_cmva") && (temp.Contains("j3t") || temp.Contains("j3t")) ) rebin = 6;

      double total_data = 0;
      double total_bkg  = 0;

      if( renormMC_ ){
	for( int iSample=0; iSample<NumSamples-2; iSample++ ){
	  TH1D* h_delete_tmp = (TH1D*)file[iSample]->Get(temp)->Clone(Form("h_delete_tmp_%d_%s_jetFlav",iSample,temp.Data()));
	  double total = h_delete_tmp->Integral();
	  if( iSample==0 ) total_data = total;
	  else             total_bkg += total;
	}
      }

      TH1D* hist_sys[NumSamples][NumSysCat];
      for( int iSample=0; iSample<NumSamples; iSample++ ){
	for( int iSys=0; iSys<NumSysCat; iSys++ ){

	  TString temp_sys = temp + sys_cat_labels[iSys];
	  if( temp.Contains("_diLepMass") ) temp_sys = temp;
	  if( temp.Contains("h_jet_cmva") ) temp_sys.ReplaceAll("h_jet_cmva","h_cmva_jet_cmva");

	  hist_sys[iSample][iSys] = (TH1D*)file[iSample]->Get(temp_sys)->Clone(Form("%s_tmp_sys_jetFlav",temp_sys.Data()));;
	  hist_sys[iSample][iSys]->Rebin(rebin);

	  if( renormMC_ && total_bkg>0 && total_data>0 && (iSample>=0 && iSample<NumSamples-2) ) hist_sys[iSample][iSys]->Scale( total_data/total_bkg );
	  if( rescaleMC_>0 && !renormMC_ && iSample>0 ) hist_sys[iSample][iSys]->Scale( rescaleMC_ );
	  //if( iSample>=1 && iSample<=5 ) hist_sys[iSample][iSys]->Scale( 2537./2430 );
	}
      }


      TH1D* hist[NumSamples];
      TH1D* hist_sum = NULL;
      bool firstFill = true;
      TH1D* hist_bkg_flav[NumSamples][NumJetFlav];
      for( int iSample=0; iSample<NumSamples; iSample++ ){
	hist[iSample] = (TH1D*)hist_sys[iSample][0]->Clone(Form("%s_%d_jetFlav",temp.Data(),iSample));

	if( iSample>0 && iSample<NumSamples-2 ){
	  for( int iFlav=0; iFlav<NumJetFlav; iFlav++ ){
	    TString temp_flav = temper + "_" + jet_flavor_suffix[iFlav] + "_" + cat_labels[iCat];
	    hist_bkg_flav[iSample][iFlav] = (TH1D*)file[iSample]->Get(temp_flav)->Clone(Form("%s_%d_jetFlav",temp_flav.Data(),iSample));
	    hist_bkg_flav[iSample][iFlav]->Rebin(rebin);
	  }
	}

	if( iSample>0 ){
	  hist[iSample]->SetLineColor(color[iSample]);
	  if( iSample<NumSamples-2 ) hist[iSample]->SetFillColor(color[iSample]);
	  else hist[iSample]->SetLineWidth(3);

	  if( firstFill && iSample>0 ){
	    firstFill = false;
	    hist_sum = (TH1D*)hist[iSample]->Clone("sum_"+temp);
	  }
	  else if( iSample<NumSamples-2 ){
	    hist_sum->Add( hist[iSample] );
	  }
	}
	else {
	  hist[iSample]->SetMarkerStyle(20);
	  //hist[iSample]->SetLineWidth(20);
	}

	if( iSample==NumSamples-1 || iSample==NumSamples-2 ) hist[iSample]->Scale(scale_ttH);

	// if( iSample==0 )                legend->AddEntry(hist[iSample],histLabels[iSample],"pe1");
	// else if( iSample<NumSamples-2 ) legend->AddEntry(hist[iSample],histLabels[iSample],"f");
	// else if( iSample==NumSamples-1 ) legend->AddEntry(hist[iSample],histLabels[iSample]+Form(" (x%d)",int(scale_ttH+0.0001)),"l");

      }// end loop over samples



      double data_integral = hist[0]->Integral();
      double bkg_integral = hist_sum->Integral();
      double integral_ratio = data_integral / bkg_integral;

      if( temp.Contains("numPV") ){
	hist_sum->Scale(integral_ratio);
	for( int iSample=1; iSample<NumSamples; iSample++ ){
	  hist[iSample]->Scale(integral_ratio);
	  for( int iSys=0; iSys<NumSysCat; iSys++ ){
	    hist_sys[iSample][iSys]->Scale(integral_ratio);
	  }
	}
      }



      int nbins = hist[0]->GetNbinsX();

      double xmin = hist[0]->GetBinLowEdge(1);
      double xmax = hist[0]->GetBinLowEdge(nbins) + hist[0]->GetBinWidth(nbins);



      if( false && (temp=="h_numBtag_wgtCSV_4j" || temp=="h_numBtag_wgtCSV2_4j") ){
	if( temp.Contains("h_numBtag_wgtCSV_4j") ) std::cout << "h_numBtag_wgtCSV_4j" << std::endl;
	else if( temp.Contains("h_numBtag_wgtCSV2_4j") ) std::cout << "h_numBtag_wgtCSV2_4j" << std::endl;

	std::cout << "   Data " << std::endl;
	for( int iBin=0; iBin<nbins; iBin++ ){
	  printf("     iBin = %d: bin content = %.1f \n", iBin, hist[0]->GetBinContent(iBin+1));
	}

	std::cout << "   Sum Bkg " << std::endl;
	for( int iBin=0; iBin<nbins; iBin++ ){
	  printf("     iBin = %d: bin content = %.1f \n", iBin, hist_sum->GetBinContent(iBin+1));
	}
      }

      legend->AddEntry(hist[0],histLabels[0],"pe1");

      TH1D* hist_bkg_sum_flav[NumJetFlav];
      for( int iFlav=0; iFlav<NumJetFlav; iFlav++ ){      
	TString temp_flav = temper + "_" + jet_flavor_suffix[iFlav] + "_" + cat_labels[iCat];
	bool firstFlavFill = true;
	for( int iSample=1; iSample<NumSamples-2; iSample++ ){
	  if( firstFlavFill ){
	    firstFlavFill = false;
	    hist_bkg_sum_flav[iFlav] = (TH1D*)hist_bkg_flav[iSample][iFlav]->Clone(Form("%s_bkg_sum_jetFlav",temp_flav.Data()));
	  }
	  else {
	    hist_bkg_sum_flav[iFlav]->Add(hist_bkg_flav[iSample][iFlav]);
	  }
	}

	hist_bkg_sum_flav[iFlav]->SetFillColor(jetflavcolor[iFlav]);

	//legend->AddEntry(hist_bkg_sum_flav[iFlav],jet_flavor_suffix[iFlav],"f");
	if( jet_flavor_suffix[iFlav]=="bFlav" ) legend->AddEntry(hist_bkg_sum_flav[iFlav],"b","f");
	if( jet_flavor_suffix[iFlav]=="cFlav" ) legend->AddEntry(hist_bkg_sum_flav[iFlav],"c","f");
	if( jet_flavor_suffix[iFlav]=="lFlav" ) legend->AddEntry(hist_bkg_sum_flav[iFlav],"udsg","f");
      }
      //legend->AddEntry(hist[NumSamples-2],histLabels[NumSamples-2]+Form(" (x%d)",int(scale_ttH+0.0001)),"l");


      THStack *hs = new THStack("hs","");
      for( int iFlav=0; iFlav<NumJetFlav; iFlav++ ){      
	hs->Add(hist_bkg_sum_flav[iFlav]);
      }


      TH1D* h_sum_bkg_sys[NumSysCat];
      for( int iSys=0; iSys<NumSysCat; iSys++ ){
	bool firstSample = true;
	for( int iSample=1; iSample<NumSamples-2; iSample++ ){
	  if( firstSample ){
	    firstSample = false;
	    h_sum_bkg_sys[iSys] = (TH1D*)hist_sys[iSample][iSys]->Clone(Form("%s_sys_%d_jetFlav",temp.Data(),iSys));
	  }
	  else h_sum_bkg_sys[iSys]->Add(hist_sys[iSample][iSys]);
	}
      }


      TH1D* h_bkg_err_1sig = new TH1D("h_bkg_err_1sig_"+temp,"", nbins, hist[0]->GetBinLowEdge(1), hist[0]->GetBinLowEdge(nbins) + hist[0]->GetBinWidth(nbins) );
      h_bkg_err_1sig->SetFillColor(kGreen);
      double sum_mc_err_tot_2 = 0.;
      for( int iBin=0; iBin<nbins; iBin++ ){

	double nom = h_sum_bkg_sys[0]->GetBinContent(iBin+1);
	double mcstat = hist_sum->GetBinError(iBin+1);

	sum_mc_err_tot_2 += mcstat*mcstat;

	double up2 = 0, down2 = 0;
	double previous_diff = 0;
	for( int iSys=1; iSys<NumSysCat; iSys++ ){
	  if( temp.Contains("h_numPVs_") && !sys_cat_labels[iSys].Contains("PU") ) continue;
	  double diff = h_sum_bkg_sys[iSys]->GetBinContent(iBin+1) - nom;
	  double diff_2 = diff*diff;
	  if( iBin!=1 && iBin%2==0 && (diff*previous_diff)>0 ){
	    if( fabs(diff)>fabs(previous_diff) ) diff_2 = diff*diff - previous_diff*previous_diff;
	  } 
	  if( diff>0 ) up2 += diff_2;
	  else if( diff<0. ) down2 += diff_2;
	  previous_diff = diff;

	  if( temp.Contains("h_numBtag_wgtCSV4_2l") && !temp.Contains("h_numBtag_wgtCSV_4j") && !temp.Contains("h_numBtag_wgtCSV_4j2t") ){
	    printf(" Sys = %20s, iBin = %d, nominal = %10.3f, diff = %+10.3f, abs(diff)/nominal = %.3f, abs(up2)/nom = %.3f, abs(down2)/nom = %.3f \n",
		   sys_cat_labels[iSys].Data(), iBin, nom, diff, fabs(diff)/nom, sqrt(up2)/nom, sqrt(down2)/nom );
	  }
	}

	up2 += mcstat*mcstat;
	down2 += mcstat*mcstat;

	double lumi_err2 = lumi_err * lumi_err * nom * nom;

	up2 += lumi_err2;
	down2 += lumi_err2;


	double up_err   = nom + sqrt(up2);
	double down_err = nom - sqrt(down2);

	double ave = 0.5 * ( up_err + down_err );
	h_bkg_err_1sig->SetBinContent(iBin+1,ave);
	h_bkg_err_1sig->SetBinError(iBin+1,up_err-ave);
      }


      double ratioMax = 1.6;
      double ratioMin = 0.5;
      // double ratioMax = 2.3;
      // double ratioMin = 0.0;

      TH1D* myRatio = new TH1D("ratio", "", nbins, xmin, xmax );
      TH1D* myRatio_1sig = new TH1D("ratio_1sig_"+temp, "", nbins, xmin, xmax );

      myRatio->SetStats(0);
      myRatio->Sumw2();
      myRatio->SetLineColor(kBlack);
      myRatio->SetMarkerColor(kBlack);
      myRatio->Divide(hist[0],hist_sum);



      for( int iBin=0; iBin<nbins; iBin++ ){
	double bkg  = h_bkg_err_1sig->GetBinContent(iBin+1);
	double bkg_1sig  = h_bkg_err_1sig->GetBinError(iBin+1);
	double data = hist[0]->GetBinContent(iBin+1);
	double ratio = ( bkg>0. ) ? data/bkg : 0.;
	double ratio_err = ( bkg>0. ) ? sqrt(data)/bkg : 0.;

	double bkg_noshift  = hist_sum->GetBinContent(iBin+1);
	if( bkg_noshift>0. ) ratio = data/bkg_noshift;
	if( bkg_noshift>0. ) ratio_err = ( bkg>0. ) ? sqrt(data)/bkg_noshift : 0.;

	//if( fabs(bkg-bkg_noshift)>0.001 ) std::cout << " hist = " << temp << ", bin = " << bin << ", bkg = " << bkg << ", bkg_noshift


	double up_err = bkg + bkg_1sig;
	double down_err = bkg - bkg_1sig;

	if( bkg>0. && bkg_noshift>0. ){
	  myRatio->SetBinContent(iBin+1,ratio);
	  myRatio->SetBinError(iBin+1,ratio_err);

	  up_err *= 1./bkg_noshift;
	  down_err *= 1./bkg_noshift;

	  double new_ave = 0.5 * ( up_err + down_err );

	  myRatio_1sig->SetBinContent(iBin+1,new_ave);
	  myRatio_1sig->SetBinError(iBin+1,up_err - new_ave);

	  //std::cout << " hist = " << temp << ", bin = " << bin << ", new_ave = " << new_ave << ", up_err = " << up_err << ", down_err = " << down_err << std::endl;
	}

	if( (ratio>ratioMax) && ((ratio - ratio_err)<ratioMax) ){
	  double minner = ratio - ratio_err;
	  myRatio->SetBinContent(iBin+1,ratioMax-0.0001);
	  myRatio->SetBinError(iBin+1,ratioMax-0.0001-minner);
	}

      }






      myRatio->SetMinimum(ratioMin);
      myRatio->SetMaximum(ratioMax);
      // double ratioMax = 2.3;
      // double ratioMin = 0.0;
      //myRatio->GetYaxis()->SetNdivisions(50000+404);
      // double ratioMax = 1.6;
      // double ratioMin = 0.5;
      myRatio->GetYaxis()->SetNdivisions(50000+204);
      myRatio->GetYaxis()->SetLabelSize(0.1); //make y label bigger
      myRatio->GetXaxis()->SetLabelSize(0.12); //make y label bigger
      myRatio->GetXaxis()->SetTitleOffset(1.1);
      myRatio->GetXaxis()->SetTitle(hist[0]->GetXaxis()->GetTitle()); //make y label bigger
      if( temp.Contains("h_num") ){
	myRatio->GetXaxis()->SetLabelSize(0.16);
	myRatio->GetXaxis()->SetLabelOffset(0.034);
      }
      else {
	myRatio->GetXaxis()->SetLabelSize(0.12);
	myRatio->GetXaxis()->SetLabelOffset(0.04);
      }
      myRatio->GetXaxis()->SetTitleSize(0.12);
      myRatio->GetYaxis()->SetTitle("Data/MC");
      myRatio->GetYaxis()->SetTitleSize(0.09);
      myRatio->GetYaxis()->SetTitleOffset(.55);
      myC1->cd(2);
      gPad->SetTopMargin(small);
      gPad->SetTickx();
      gPad->Modified();

      myRatio->GetYaxis()->CenterTitle(kTRUE);


      if( temp.Contains("_category_yield") || temp.Contains("_selection") ){
	for( int iBin=0; iBin<nbins; iBin++ ) myRatio->GetXaxis()->SetBinLabel(iBin+1,hist[0]->GetXaxis()->GetBinLabel(iBin+1));
      }

      if( temp.Contains("h_num") && !temp.Contains("numPV") ){
	for( int iBin=0; iBin<nbins; iBin++ ) myRatio->GetXaxis()->SetBinLabel(iBin+1,Form("%d",int(0.001 + hist[0]->GetXaxis()->GetBinCenter(iBin+1))));
	for( int iBin=0; iBin<nbins; iBin++ ) hist[0]->GetXaxis()->SetBinLabel(iBin+1,Form("%d",int(0.001 + hist[0]->GetXaxis()->GetBinCenter(iBin+1))));
      }


      hist[0]->SetStats(0);

      //hist[0]->GetYaxis()->SetTitleOffset(1.2);
      hist[0]->GetYaxis()->SetTitleOffset(1.0);
      hist[0]->GetYaxis()->SetTitleSize(0.05);

      hist[0]->GetYaxis()->SetTitle("Number of Events");


      int max_bin_data = hist[0]->GetMaximumBin();
      double max_data = hist[0]->GetBinContent(max_bin_data) + hist[0]->GetBinError(max_bin_data);

      int max_bin_mc = h_bkg_err_1sig->GetMaximumBin();
      double max_mc = h_bkg_err_1sig->GetBinContent(max_bin_mc) + h_bkg_err_1sig->GetBinError(max_bin_mc);

      double max_content = std::max(max_data, max_mc);

      hist[0]->GetYaxis()->SetRangeUser(0.,1.2 * max_content);
      //hist[0]->GetYaxis()->SetRangeUser(0.,1.2 * hist[0]->GetBinContent(2));


      bool rerange = false;

      if( (temp.Contains("_jet_csv")) ){ rerange=true; xmin = -0.04; xmax = 1.028; }
      if( (temp.Contains("_HT")) ){ rerange=true; xmin = 100+0.0001; xmax = 1400-0.001; }
      if( (temp.Contains("_mht_pt")) ){ rerange=true; xmin = 0+0.0001; xmax = 200-0.001; }
      if( (temp.Contains("_met_pt")) ){ rerange=true; xmin = 0+0.0001; xmax = 400-0.001; }
      //if( (temp.Contains("_category")) ){ rerange=true; xmin = 0.001; xmax = 6.99; }

      if( temp.Contains("h_numBtag_wgtCSV_4j2t") ){ rerange=true; xmin = 1.5+0.0001; xmax = 4.5-0.001; }

      if( temp.Contains("h_numBtag") && temp.Contains("4j") && !temp.Contains("4j2t") && !temp.Contains("met") ){ rerange=true; xmin = 0.5+0.0001; xmax = 4.5-0.001; }
      if( temp.Contains("h_numBtag") && temp.Contains("4j") && !temp.Contains("4j2t") && temp.Contains("met")  ){ rerange=true; xmin = -0.5+0.0001; xmax = 4.5-0.001; }
      if( temp.Contains("h_numBtag") && temp.Contains("2l") ){ rerange=true; xmin = 0.5+0.0001; xmax = 2.5-0.001; }

      if( rerange ){
	hist[0]->GetXaxis()->SetRangeUser(xmin,xmax);
	myRatio->GetXaxis()->SetRangeUser(xmin,xmax);
      }


      if( temp.Contains("csv") ){
	hist[0]->GetYaxis()->SetTitle("Jets/0.024");
	myRatio->GetXaxis()->SetTitle("CSVv2 Discriminator");
      }
      
      // if( temp.Contains("_mass") ){
      //   hist[0]->GetXaxis()->SetRangeUser(40.,140.);
      //   myRatio->GetXaxis()->SetRangeUser(40.,140.);
      // }

      // if( temp.Contains("_met") ){
      //   hist[0]->GetXaxis()->SetRangeUser(0.,150.);
      //   myRatio->GetXaxis()->SetRangeUser(0.,150.);
      // }
      // // if( temp.Contains("_eta") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);
      // // if( temp.Contains("_numGenPVs") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);
      // // if( temp.Contains("_numJets") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);

      h_bkg_err_1sig->SetFillStyle(3654);
      h_bkg_err_1sig->SetFillColor(kBlack);
      h_bkg_err_1sig->SetMarkerSize(0);


      myRatio_1sig->SetMarkerColor(kGreen);
      myRatio_1sig->SetFillColor(kGreen);
      myRatio_1sig->SetMarkerSize(0); // geoff


      TLine* myLine;
      myLine = new TLine(xmin, 1, xmax, 1);


      // Plot
      myC1->cd(1);
      hist[0]->Draw("pe1");
      hs->Draw("histsame");
      h_bkg_err_1sig->Draw("e2same");
      //hist[NumSamples-2]->Draw("histsame");
      //hist[NumSamples-1]->Draw("histsame");
      hist[0]->Draw("pe1same");

      legend->Draw();

      LumiInfoLatex2.Draw();
      CMSInfoLatex2.Draw();
      PublishInfoLatex2.Draw();

      myC1->cd(2);
      myRatio->SetLineWidth(2);
      myRatio->Draw("pe1");
      myRatio_1sig->Draw("e2same");
      myRatio->Draw("pe1same");

      myLine->Draw();

      myC1->GetPad(1)->SetLogy(0);

      myC1->GetPad(1)->RedrawAxis();
      myC1->GetPad(2)->RedrawAxis();

      plotname = dirprefix + temp_mh + "_data2mc_lin.png";
      myC1->Print(plotname);

      plotname = dirprefix + temp_mh + "_data2mc_lin.pdf";
      if( printPDF_ ) myC1->Print(plotname);


      if( temp.Contains("_category_yield") || temp.Contains("h_num") || temp.Contains("csv") ){
	// log
	hist[0]->GetYaxis()->SetRangeUser(0.4,40 * max_content);

	// if including the zero bin
	if( temp.Contains("numBtag") && temp.Contains("4j") && temp.Contains("met") ) hist[0]->GetYaxis()->SetRangeUser(0.4,40 * max_content);

	// if not including the zero bin
	if( temp.Contains("numBtag") && temp.Contains("4j") && !temp.Contains("met") ) hist[0]->GetYaxis()->SetRangeUser(0.4,12 * max_content);

	if( temp.Contains("_selection") ){
	  hist[0]->GetYaxis()->SetRangeUser(1000,12 * max_content);
	}

	myC1->cd(1);
	//gPad->SetLogy(1);

	myC1->GetPad(1)->SetLogy(1);

	myC1->cd(1);
	hist[0]->Draw("pe1");
	hs->Draw("histsame");
	h_bkg_err_1sig->Draw("e2same");
	//hist[NumSamples-2]->Draw("histsame");
	//hist[NumSamples-1]->Draw("histsame");
	hist[0]->Draw("pe1same");

	legend->Draw();

	LumiInfoLatex.Draw();
	CMSInfoLatex.Draw();
	PublishInfoLatex.Draw();

	myC1->cd(2);
	myRatio->SetLineWidth(2);
	myRatio->Draw("pe1");
	myRatio_1sig->Draw("e2same");
	myRatio->Draw("pe1same");
	myLine->Draw();

	myC1->GetPad(1)->RedrawAxis();
	myC1->GetPad(2)->RedrawAxis();

	plotname = dirprefix + temp_mh + "_data2mc_log.png";
	myC1->Print(plotname);

	plotname = dirprefix + temp_mh + "_data2mc_log.pdf";
	if( printPDF_ ) myC1->Print(plotname);

	//gPad->SetLogy(0);
      }


      delete myRatio;
      delete myRatio_1sig;
      delete h_bkg_err_1sig;
      delete myLine;
      delete legend;
    } // end loop on categories
  }// end loop on hists histoname5




  // Close the files
  std::cout << " Closing all files..." << std::endl;
  for( int iFile=0; iFile<NumSamples; iFile++ ) file[iFile]->Close();
  std::cout << " Done! " << std::endl;
}


/*


  std::vector<TString> sys_cat_labels;
  sys_cat_labels.push_back("No CSV");               //0
  sys_cat_labels.push_back("No Sys");               //0
  sys_cat_labels.push_back("lepIdSFUp");     //1
  sys_cat_labels.push_back("lepIdSFDown");   //2
  sys_cat_labels.push_back("PUUp");          //3
  sys_cat_labels.push_back("PUDown");        //4
  sys_cat_labels.push_back("JERUp");         //5
  sys_cat_labels.push_back("JERDown");       //6
  sys_cat_labels.push_back("JESUp");         //7
  sys_cat_labels.push_back("JESDown");       //8
  sys_cat_labels.push_back("CSVLFUp");       //9
  sys_cat_labels.push_back("CSVLFDown");     //10
  sys_cat_labels.push_back("CSVHFUp");       //11
  sys_cat_labels.push_back("CSVHFDown");     //12
  sys_cat_labels.push_back("CSVHFStats1Up");     //13
  sys_cat_labels.push_back("CSVHFStats1Down");   //14
  sys_cat_labels.push_back("CSVHFStats2Up");     //15
  sys_cat_labels.push_back("CSVHFStats2Down");   //16
  sys_cat_labels.push_back("CSVLFStats1Up");     //17
  sys_cat_labels.push_back("CSVLFStats1Down");   //18
  sys_cat_labels.push_back("CSVLFStats2Up");     //19
  sys_cat_labels.push_back("CSVLFStats2Down");   //20
  sys_cat_labels.push_back("CSVCFErr1Up");     //21
  sys_cat_labels.push_back("CSVCFErr1Down");   //22
  sys_cat_labels.push_back("CSVCFErr2Up");     //23
  sys_cat_labels.push_back("CSVCFErr2Down");   //24
  sys_cat_labels.push_back("muFUp");           //25
  sys_cat_labels.push_back("muFDown");         //26
  sys_cat_labels.push_back("muRUp");           //27
  sys_cat_labels.push_back("muRDown");         //28
  sys_cat_labels.push_back("muRmuFUp");        //29
  sys_cat_labels.push_back("muRmuFDown");      //30

  int NumSysCat = int(sys_cat_labels.size());

  for( int iBin=0; iBin<NumSysCat; iBin++ ) h_numEvents_perSys_wgtCSV->GetXaxis()->SetBinLabel(iBin+1, sys_cat_labels[iBin] );



  sys_cat_labels.push_back("_CSVCFErr1Up");     //21
  sys_cat_labels.push_back("_CSVCFErr1Down");   //22
  sys_cat_labels.push_back("_CSVCFErr2Up");     //23
  sys_cat_labels.push_back("_CSVCFErr2Down");   //24



h_numBtag_wgtCSV2_4j_CSVHFUp->Divide(h_numBtag_wgtCSV_4j_CSVHFUp);
h_numBtag_wgtCSV2_4j_CSVHFUp->Draw();

h_numBtag_wgtCSV2_4j_CSVHFDown->Divide(h_numBtag_wgtCSV_4j_CSVHFDown);
h_numBtag_wgtCSV2_4j_CSVHFDown->Draw();


h_numBtag_wgtCSV2_4j_CSVHFStats1Up->Divide(h_numBtag_wgtCSV_4j_CSVHFStats1Up);
h_numBtag_wgtCSV2_4j_CSVHFStats1Up->Draw();

h_numBtag_wgtCSV2_4j_CSVHFStats1Down->Divide(h_numBtag_wgtCSV_4j_CSVHFStats1Down);
h_numBtag_wgtCSV2_4j_CSVHFStats1Down->Draw();


h_numBtag_wgtCSV2_4j_CSVHFStats2Up->Divide(h_numBtag_wgtCSV_4j_CSVHFStats2Up);
h_numBtag_wgtCSV2_4j_CSVHFStats2Up->Draw();

h_numBtag_wgtCSV2_4j_CSVHFStats2Down->Divide(h_numBtag_wgtCSV_4j_CSVHFStats2Down);
h_numBtag_wgtCSV2_4j_CSVHFStats2Down->Draw();


h_numBtag_wgtCSV2_4j_CSVLFUp->Divide(h_numBtag_wgtCSV_4j_CSVLFUp);
h_numBtag_wgtCSV2_4j_CSVLFUp->Draw();

h_numBtag_wgtCSV2_4j_CSVLFDown->Divide(h_numBtag_wgtCSV_4j_CSVLFDown);
h_numBtag_wgtCSV2_4j_CSVLFDown->Draw();



h_numBtag_wgtCSV2_4j_CSVLFStats1Up->Divide(h_numBtag_wgtCSV_4j_CSVLFStats1Up);
h_numBtag_wgtCSV2_4j_CSVLFStats1Up->Draw();

h_numBtag_wgtCSV2_4j_CSVLFStats1Down->Divide(h_numBtag_wgtCSV_4j_CSVLFStats1Down);
h_numBtag_wgtCSV2_4j_CSVLFStats1Down->Draw();


h_numBtag_wgtCSV2_4j_CSVLFStats2Up->Divide(h_numBtag_wgtCSV_4j_CSVLFStats2Up);
h_numBtag_wgtCSV2_4j_CSVLFStats2Up->Draw();

h_numBtag_wgtCSV2_4j_CSVLFStats2Down->Divide(h_numBtag_wgtCSV_4j_CSVLFStats2Down);
h_numBtag_wgtCSV2_4j_CSVLFStats2Down->Draw();



h_numBtag_wgtCSV2_4j_JESUp->Divide(h_numBtag_wgtCSV_4j_JESUp);
h_numBtag_wgtCSV2_4j_JESUp->Draw();

h_numBtag_wgtCSV2_4j_JESDown->Divide(h_numBtag_wgtCSV_4j_JESDown);
h_numBtag_wgtCSV2_4j_JESDown->Draw();



h_numBtag_wgtCSV2_4j_CSVCFErr1Up->Divide(h_numBtag_wgtCSV_4j_CSVCFErr1Up);
h_numBtag_wgtCSV2_4j_CSVCFErr1Up->Draw();

h_numBtag_wgtCSV2_4j_CSVCFErr1Down->Divide(h_numBtag_wgtCSV_4j_CSVCFErr1Down);
h_numBtag_wgtCSV2_4j_CSVCFErr1Down->Draw();


h_numBtag_wgtCSV2_4j_CSVCFErr2Up->Divide(h_numBtag_wgtCSV_4j_CSVCFErr2Up);
h_numBtag_wgtCSV2_4j_CSVCFErr2Up->Draw();

h_numBtag_wgtCSV2_4j_CSVCFErr2Down->Divide(h_numBtag_wgtCSV_4j_CSVCFErr2Down);
h_numBtag_wgtCSV2_4j_CSVCFErr2Down->Draw();


 */
