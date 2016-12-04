#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH3.h"
#include "TH2F.h"
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
#include "TEfficiency.h"
#include "THStack.h"

void combine_control_region() {
	
	TH1::SetDefaultSumw2();
	
	double L = 26.149; // 1/fb
	
	TFile *f1 = new TFile("Control_Region_mc_dy.root");
	TH1D *dy_pt_loose_27 = (TH1D*)f1->Get("Pt_WPLoose_27");
	TH1D *dy_pt2_loose_27 = (TH1D*)f1->Get("Pt2_WPLoose_27");
	TH1D *dy_eta_loose_27 = (TH1D*)f1->Get("Eta_WPLoose_27");
	TH1D *dy_phi_loose_27 = (TH1D*)f1->Get("Phi_WPLoose_27");
	TH1D *dy_ht_loose_27 = (TH1D*)f1->Get("HT_WPLoose_27");
	TH1D *dy_npv_loose_27 = (TH1D*)f1->Get("NumPV_WPLoose_27");
	TH1D *dy_jet1pt_loose_27 = (TH1D*)f1->Get("Jet1_pt_WPLoose_27");
	TH1D *dy_jet1csv_loose_27 = (TH1D*)f1->Get("Jet1_csv_WPLoose_27");
	TH1D *dy_njets_loose_27 = (TH1D*)f1->Get("Njets_WPLoose_27");
	TH1D *dy_nbtags_loose_27 = (TH1D*)f1->Get("Nbtags_WPLoose_27");
	TH1D *dy_pt_tight_27 = (TH1D*)f1->Get("Pt_WPTight_27");
	TH1D *dy_pt2_tight_27 = (TH1D*)f1->Get("Pt2_WPTight_27");
	TH1D *dy_eta_tight_27 = (TH1D*)f1->Get("Eta_WPTight_27");
	TH1D *dy_phi_tight_27 = (TH1D*)f1->Get("Phi_WPTight_27");
	TH1D *dy_ht_tight_27 = (TH1D*)f1->Get("HT_WPTight_27");
	TH1D *dy_npv_tight_27 = (TH1D*)f1->Get("NumPV_WPTight_27");
	TH1D *dy_jet1pt_tight_27 = (TH1D*)f1->Get("Jet1_pt_WPTight_27");
	TH1D *dy_jet1csv_tight_27 = (TH1D*)f1->Get("Jet1_csv_WPTight_27");
	TH1D *dy_njets_tight_27 = (TH1D*)f1->Get("Njets_WPTight_27");
	TH1D *dy_nbtags_tight_27 = (TH1D*)f1->Get("Nbtags_WPTight_27");
	TH1D *dy_pt_tight_32 = (TH1D*)f1->Get("Pt_WPTight_32");
	TH1D *dy_pt2_tight_32 = (TH1D*)f1->Get("Pt2_WPTight_32");
	TH1D *dy_eta_tight_32 = (TH1D*)f1->Get("Eta_WPTight_32");
	TH1D *dy_phi_tight_32 = (TH1D*)f1->Get("Phi_WPTight_32");
	TH1D *dy_ht_tight_32 = (TH1D*)f1->Get("HT_WPTight_32");
	TH1D *dy_npv_tight_32 = (TH1D*)f1->Get("NumPV_WPTight_32");
	TH1D *dy_jet1pt_tight_32 = (TH1D*)f1->Get("Jet1_pt_WPTight_32");
	TH1D *dy_jet1csv_tight_32 = (TH1D*)f1->Get("Jet1_csv_WPTight_32");
	TH1D *dy_njets_tight_32 = (TH1D*)f1->Get("Njets_WPTight_32");
	TH1D *dy_nbtags_tight_32 = (TH1D*)f1->Get("Nbtags_WPTight_32");
	TH1D *dy_pt_loose_27_ht200 = (TH1D*)f1->Get("Pt_WPLoose_27_HT200");
	TH1D *dy_pt2_loose_27_ht200 = (TH1D*)f1->Get("Pt2_WPLoose_27_HT200");
	TH1D *dy_eta_loose_27_ht200 = (TH1D*)f1->Get("Eta_WPLoose_27_HT200");
	TH1D *dy_phi_loose_27_ht200 = (TH1D*)f1->Get("Phi_WPLoose_27_HT200");
	TH1D *dy_ht_loose_27_ht200 = (TH1D*)f1->Get("HT_WPLoose_27_HT200");
	TH1D *dy_npv_loose_27_ht200 = (TH1D*)f1->Get("NumPV_WPLoose_27_HT200");
	TH1D *dy_jet1pt_loose_27_ht200 = (TH1D*)f1->Get("Jet1_pt_WPLoose_27_HT200");
	TH1D *dy_jet1csv_loose_27_ht200 = (TH1D*)f1->Get("Jet1_csv_WPLoose_27_HT200");
	TH1D *dy_njets_loose_27_ht200 = (TH1D*)f1->Get("Njets_WPLoose_27_HT200");
	TH1D *dy_nbtags_loose_27_ht200 = (TH1D*)f1->Get("Nbtags_WPLoose_27_HT200");
	TH1D *dy_pt_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *dy_pt2_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Pt2_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *dy_eta_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Eta_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *dy_phi_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Phi_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *dy_ht_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("HT_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *dy_npv_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("NumPV_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *dy_jet1pt_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Jet1_pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *dy_jet1csv_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Jet1_csv_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *dy_njets_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Njets_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *dy_nbtags_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Nbtags_WPTight_27_OR_WPLoose_27_HT200");
	
	dy_pt_loose_27->SetName("dy_pt_WPLoose_27");
	dy_pt_loose_27->SetTitle("MC dy : Tag Electron pT Distribution for WPLoose_27");
	dy_pt2_loose_27->SetName("dy_pt2_WPLoose_27");
	dy_pt2_loose_27->SetTitle("MC dy : Probe Electron pT Distribution for WPLoose_27");
	dy_eta_loose_27->SetName("dy_eta_WPLoose_27");
	dy_eta_loose_27->SetTitle("MC dy : eta Distribution for WPLoose_27");
	dy_phi_loose_27->SetName("dy_phi_WPLoose_27");
	dy_phi_loose_27->SetTitle("MC dy : phi Distribution for WPLoose_27");
	dy_ht_loose_27->SetName("dy_ht_WPLoose_27");
	dy_ht_loose_27->SetTitle("MC dy : HT Distribution for WPLoose_27");
	dy_npv_loose_27->SetName("dy_npv_WPLoose_27");
	dy_npv_loose_27->SetTitle("MC dy : NPV Distribution for WPLoose_27");
	dy_jet1pt_loose_27->SetName("dy_jet1pt_WPLoose_27");
	dy_jet1pt_loose_27->SetTitle("MC dy : Jet1 pT Distribution for WPLoose_27");
	dy_jet1csv_loose_27->SetName("dy_jet1csv_WPLoose_27");
	dy_jet1csv_loose_27->SetTitle("MC dy : Jet1 CSV Distribution for WPLoose_27");
	dy_njets_loose_27->SetName("dy_njets_WPLoose_27");
	dy_njets_loose_27->SetTitle("MC dy : Nr. of Jets Distribution for WPLoose_27");
	dy_nbtags_loose_27->SetName("dy_nbtags_WPLoose_27");
	dy_nbtags_loose_27->SetTitle("MC dy : Nr. of Btags Distribution for WPLoose_27");
	dy_pt_tight_27->SetName("dy_pt_WPTight_27");
	dy_pt_tight_27->SetTitle("MC dy : Tag ELectron pT Distribution for WPTight_27");
	dy_pt2_tight_27->SetName("dy_pt2_WPTight_27");
	dy_pt2_tight_27->SetTitle("MC dy : Probe Electron pT Distribution for WPTight_27");
	dy_eta_tight_27->SetName("dy_eta_WPTight_27");
	dy_eta_tight_27->SetTitle("MC dy : eta Distribution for WPTight_27");
	dy_phi_tight_27->SetName("dy_phi_WPTight_27");
	dy_phi_tight_27->SetTitle("MC dy : phi Distribution for WPTight_27");
	dy_ht_tight_27->SetName("dy_ht_WPTight_27");
	dy_ht_tight_27->SetTitle("MC dy : HT Distribution for WPTight_27");
	dy_npv_tight_27->SetName("dy_npv_WPTight_27");
	dy_npv_tight_27->SetTitle("MC dy : NPV Distribution for WPTight_27");
	dy_jet1pt_tight_27->SetName("dy_jet1pt_WPTight_27");
	dy_jet1pt_tight_27->SetTitle("MC dy : Jet1 pT Distribution for WPTight_27");
	dy_jet1csv_tight_27->SetName("dy_jet1csv_WPTight_27");
	dy_jet1csv_tight_27->SetTitle("MC dy : Jet1 CSV Distribution for WPTight_27");
	dy_njets_tight_27->SetName("dy_njets_WPTight_27");
	dy_njets_tight_27->SetTitle("MC dy : Nr. of Jets Distribution for WPTight_27");
	dy_nbtags_tight_27->SetName("dy_nbtags_WPTight_27");
	dy_nbtags_tight_27->SetTitle("MC dy : Nr. of Btags Distribution for WPTight_27");
	dy_pt_tight_32->SetName("dy_pt_WPTight_32");
	dy_pt_tight_32->SetTitle("MC dy : Tag Electron pT Distribution for WPTight_32");
	dy_pt2_tight_32->SetName("dy_pt2_WPTight_32");
	dy_pt2_tight_32->SetTitle("MC dy : Probe Electron pT Distribution for WPTight_32");
	dy_eta_tight_32->SetName("dy_eta_WPTight_32");
	dy_eta_tight_32->SetTitle("MC dy : eta Distribution for WPTight_32");
	dy_phi_tight_32->SetName("dy_phi_WPTight_32");
	dy_phi_tight_32->SetTitle("MC dy : phi Distribution for WPTight_32");
	dy_ht_tight_32->SetName("dy_ht_WPTight_32");
	dy_ht_tight_32->SetTitle("MC dy : HT Distribution for WPTight_32");
	dy_npv_tight_32->SetName("dy_npv_WPTight_32");
	dy_npv_tight_32->SetTitle("MC dy : NPV Distribution for WPTight_32");
	dy_jet1pt_tight_32->SetName("dy_jet1pt_WPTight_32");
	dy_jet1pt_tight_32->SetTitle("MC dy : Jet1 pT Distribution for WPTight_32");
	dy_jet1csv_tight_32->SetName("dy_jet1csv_WPTight_32");
	dy_jet1csv_tight_32->SetTitle("MC dy : Jet1 CSV Distribution for WPTight_32");
	dy_njets_tight_32->SetName("dy_njets_WPTight_32");
	dy_njets_tight_32->SetTitle("MC dy : Nr. of Jets Distribution for WPTight_32");
	dy_nbtags_tight_32->SetName("dy_nbtags_WPTight_32");
	dy_nbtags_tight_32->SetTitle("MC dy : Nr. of Btags Distribution for WPTight_32");
	dy_pt_loose_27_ht200->SetName("dy_pt_WPloose_27_ht200");
	dy_pt_loose_27_ht200->SetTitle("MC dy : Tag Electron pT Distribution for WPloose_27_ht200");
	dy_pt2_loose_27_ht200->SetName("dy_pt2_WPloose_27_ht200");
	dy_pt2_loose_27_ht200->SetTitle("MC dy : Probe Electron pT Distribution for WPloose_27_ht200");
	dy_eta_loose_27_ht200->SetName("dy_eta_WPloose_27_ht200");
	dy_eta_loose_27_ht200->SetTitle("MC dy : eta Distribution for WPloose_27_ht200");
	dy_phi_loose_27_ht200->SetName("dy_phi_WPloose_27_ht200");
	dy_phi_loose_27_ht200->SetTitle("MC dy : phi Distribution for WPloose_27_ht200");
	dy_ht_loose_27_ht200->SetName("dy_ht_WPloose_27_ht200");
	dy_ht_loose_27_ht200->SetTitle("MC dy : HT Distribution for WPloose_27_ht200");
	dy_npv_loose_27_ht200->SetName("dy_npv_WPloose_27_ht200");
	dy_npv_loose_27_ht200->SetTitle("MC dy : NPV Distribution for WPloose_27_ht200");
	dy_jet1pt_loose_27_ht200->SetName("dy_jet1pt_WPloose_27_ht200");
	dy_jet1pt_loose_27_ht200->SetTitle("MC dy : Jet1 pT Distribution for WPloose_27_ht200");
	dy_jet1csv_loose_27_ht200->SetName("dy_jet1csv_WPloose_27_ht200");
	dy_jet1csv_loose_27_ht200->SetTitle("MC dy : Jet1 CSV Distribution for WPloose_27_ht200");
	dy_njets_loose_27_ht200->SetName("dy_njets_WPloose_27_ht200");
	dy_njets_loose_27_ht200->SetTitle("MC dy : Nr. of Jets Distribution for WPloose_27_ht200");
	dy_nbtags_loose_27_ht200->SetName("dy_nbtags_WPloose_27_ht200");
	dy_nbtags_loose_27_ht200->SetTitle("MC dy : Nr. of Btags Distribution for WPloose_27_ht200");
	dy_pt_tight_27_or_loose_27_ht200->SetName("dy_pt_WPtight_27_or_loose_27_ht200");
	dy_pt_tight_27_or_loose_27_ht200->SetTitle("MC dy : Tag Electron pT Distribution for WPtight_27_or_loose_27_ht200");
	dy_pt2_tight_27_or_loose_27_ht200->SetName("dy_pt2_WPtight_27_or_loose_27_ht200");
	dy_pt2_tight_27_or_loose_27_ht200->SetTitle("MC dy : Probe Electron pT Distribution for WPtight_27_or_loose_27_ht200");
	dy_eta_tight_27_or_loose_27_ht200->SetName("dy_eta_WPtight_27_or_loose_27_ht200");
	dy_eta_tight_27_or_loose_27_ht200->SetTitle("MC dy : eta Distribution for WPtight_27_or_loose_27_ht200");
	dy_phi_tight_27_or_loose_27_ht200->SetName("dy_phi_WPtight_27_or_loose_27_ht200");
	dy_phi_tight_27_or_loose_27_ht200->SetTitle("MC dy : phi Distribution for WPtight_27_or_loose_27_ht200");
	dy_ht_tight_27_or_loose_27_ht200->SetName("dy_ht_WPtight_27_or_loose_27_ht200");
	dy_ht_tight_27_or_loose_27_ht200->SetTitle("MC dy : HT Distribution for WPtight_27_or_loose_27_ht200");
	dy_npv_tight_27_or_loose_27_ht200->SetName("dy_npv_WPtight_27_or_loose_27_ht200");
	dy_npv_tight_27_or_loose_27_ht200->SetTitle("MC dy : NPV Distribution for WPtight_27_or_loose_27_ht200");
	dy_jet1pt_tight_27_or_loose_27_ht200->SetName("dy_jet1pt_WPtight_27_or_loose_27_ht200");
	dy_jet1pt_tight_27_or_loose_27_ht200->SetTitle("MC dy : Jet1 pT Distribution for WPtight_27_or_loose_27_ht200");
	dy_jet1csv_tight_27_or_loose_27_ht200->SetName("dy_jet1csv_WPtight_27_or_loose_27_ht200");
	dy_jet1csv_tight_27_or_loose_27_ht200->SetTitle("MC dy : Jet1 CSV Distribution for WPtight_27_or_loose_27_ht200");
	dy_njets_tight_27_or_loose_27_ht200->SetName("dy_njets_WPtight_27_or_loose_27_ht200");
	dy_njets_tight_27_or_loose_27_ht200->SetTitle("MC dy : Nr. of Jets Distribution for WPtight_27_or_loose_27_ht200");
	dy_nbtags_tight_27_or_loose_27_ht200->SetName("dy_nbtags_WPtight_27_or_loose_27_ht200");
	dy_nbtags_tight_27_or_loose_27_ht200->SetTitle("MC dy : Nr. of Btags Distribution for WPtight_27_or_loose_27_ht200");
	
	
	TFile *f5 = new TFile("Control_Region_data_2016.root");
	TH1D *data_pt_loose_27 = (TH1D*)f5->Get("Pt_WPLoose_27");
	TH1D *data_pt2_loose_27 = (TH1D*)f5->Get("Pt2_WPLoose_27");
	TH1D *data_eta_loose_27 = (TH1D*)f5->Get("Eta_WPLoose_27");
	TH1D *data_phi_loose_27 = (TH1D*)f5->Get("Phi_WPLoose_27");
	TH1D *data_ht_loose_27 = (TH1D*)f5->Get("HT_WPLoose_27");
	TH1D *data_npv_loose_27 = (TH1D*)f5->Get("NumPV_WPLoose_27");
	TH1D *data_jet1pt_loose_27 = (TH1D*)f5->Get("Jet1_pt_WPLoose_27");
	TH1D *data_jet1csv_loose_27 = (TH1D*)f5->Get("Jet1_csv_WPLoose_27");
	TH1D *data_njets_loose_27 = (TH1D*)f5->Get("Njets_WPLoose_27");
	TH1D *data_nbtags_loose_27 = (TH1D*)f5->Get("Nbtags_WPLoose_27");
	TH1D *data_pt_tight_27 = (TH1D*)f5->Get("Pt_WPTight_27");
	TH1D *data_pt2_tight_27 = (TH1D*)f5->Get("Pt2_WPTight_27");
	TH1D *data_eta_tight_27 = (TH1D*)f5->Get("Eta_WPTight_27");
	TH1D *data_phi_tight_27 = (TH1D*)f5->Get("Phi_WPTight_27");
	TH1D *data_ht_tight_27 = (TH1D*)f5->Get("HT_WPTight_27");
	TH1D *data_npv_tight_27 = (TH1D*)f5->Get("NumPV_WPTight_27");
	TH1D *data_jet1pt_tight_27 = (TH1D*)f5->Get("Jet1_pt_WPTight_27");
	TH1D *data_jet1csv_tight_27 = (TH1D*)f5->Get("Jet1_csv_WPTight_27");
	TH1D *data_njets_tight_27 = (TH1D*)f5->Get("Njets_WPTight_27");
	TH1D *data_nbtags_tight_27 = (TH1D*)f5->Get("Nbtags_WPTight_27");
	TH1D *data_pt_tight_32 = (TH1D*)f5->Get("Pt_WPTight_32");
	TH1D *data_pt2_tight_32 = (TH1D*)f5->Get("Pt2_WPTight_32");
	TH1D *data_eta_tight_32 = (TH1D*)f5->Get("Eta_WPTight_32");
	TH1D *data_phi_tight_32 = (TH1D*)f5->Get("Phi_WPTight_32");
	TH1D *data_ht_tight_32 = (TH1D*)f5->Get("HT_WPTight_32");
	TH1D *data_npv_tight_32 = (TH1D*)f5->Get("NumPV_WPTight_32");
	TH1D *data_jet1pt_tight_32 = (TH1D*)f5->Get("Jet1_pt_WPTight_32");
	TH1D *data_jet1csv_tight_32 = (TH1D*)f5->Get("Jet1_csv_WPTight_32");
	TH1D *data_njets_tight_32 = (TH1D*)f5->Get("Njets_WPTight_32");
	TH1D *data_nbtags_tight_32 = (TH1D*)f5->Get("Nbtags_WPTight_32");
	TH1D *data_pt_loose_27_ht200 = (TH1D*)f5->Get("Pt_WPLoose_27_HT200");
	TH1D *data_pt2_loose_27_ht200 = (TH1D*)f5->Get("Pt2_WPLoose_27_HT200");
	TH1D *data_eta_loose_27_ht200 = (TH1D*)f5->Get("Eta_WPLoose_27_HT200");
	TH1D *data_phi_loose_27_ht200 = (TH1D*)f5->Get("Phi_WPLoose_27_HT200");
	TH1D *data_ht_loose_27_ht200 = (TH1D*)f5->Get("HT_WPLoose_27_HT200");
	TH1D *data_npv_loose_27_ht200 = (TH1D*)f5->Get("NumPV_WPLoose_27_HT200");
	TH1D *data_jet1pt_loose_27_ht200 = (TH1D*)f5->Get("Jet1_pt_WPLoose_27_HT200");
	TH1D *data_jet1csv_loose_27_ht200 = (TH1D*)f5->Get("Jet1_csv_WPLoose_27_HT200");
	TH1D *data_njets_loose_27_ht200 = (TH1D*)f5->Get("Njets_WPLoose_27_HT200");
	TH1D *data_nbtags_loose_27_ht200 = (TH1D*)f5->Get("Nbtags_WPLoose_27_HT200");
	TH1D *data_pt_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_pt2_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Pt2_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_eta_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Eta_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_phi_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Phi_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_ht_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("HT_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_npv_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("NumPV_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_jet1pt_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Jet1_pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_jet1csv_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Jet1_csv_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_njets_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Njets_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_nbtags_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Nbtags_WPTight_27_OR_WPLoose_27_HT200");
	
	data_pt_loose_27->SetName("data_pt_WPLoose_27");
	data_pt_loose_27->SetTitle("DATA : Tag Electron pT Distribution for WPLoose_27");
	data_pt2_loose_27->SetName("data_pt2_WPLoose_27");
	data_pt2_loose_27->SetTitle("DATA : Probe ELectron pT Distribution for WPLoose_27");
	data_eta_loose_27->SetName("data_eta_WPLoose_27");
	data_eta_loose_27->SetTitle("DATA : eta Distribution for WPLoose_27");
	data_phi_loose_27->SetName("data_phi_WPLoose_27");
	data_phi_loose_27->SetTitle("DATA : phi Distribution for WPLoose_27");
	data_ht_loose_27->SetName("data_ht_WPLoose_27");
	data_ht_loose_27->SetTitle("DATA : HT Distribution for WPLoose_27");
	data_npv_loose_27->SetName("data_npv_WPLoose_27");
	data_npv_loose_27->SetTitle("DATA : NPV Distribution for WPLoose_27");
	data_jet1pt_loose_27->SetName("data_jet1pt_WPLoose_27");
	data_jet1pt_loose_27->SetTitle("DATA : Jet1 pT Distribution for WPLoose_27");
	data_jet1csv_loose_27->SetName("data_jet1csv_WPLoose_27");
	data_jet1csv_loose_27->SetTitle("DATA : Jet1 CSV Distribution for WPLoose_27");
	data_njets_loose_27->SetName("data_njets_WPLoose_27");
	data_njets_loose_27->SetTitle("DATA : Nr. of Jets Distribution for WPLoose_27");
	data_nbtags_loose_27->SetName("data_nbtags_WPLoose_27");
	data_nbtags_loose_27->SetTitle("DATA : Nr. of Btags Distribution for WPLoose_27");	
	data_pt_tight_27->SetName("data_pt_WPTight_27");
	data_pt_tight_27->SetTitle("DATA : Tag Electron pT Distribution for WPTight_27");
	data_pt2_tight_27->SetName("data_pt2_WPTight_27");
	data_pt2_tight_27->SetTitle("DATA : Probe Electron pT Distribution for WPTight_27");
	data_eta_tight_27->SetName("data_eta_WPTight_27");
	data_eta_tight_27->SetTitle("DATA : eta Distribution for WPTight_27");
	data_phi_tight_27->SetName("data_phi_WPTight_27");
	data_phi_tight_27->SetTitle("DATA : phi Distribution for WPTight_27");
	data_ht_tight_27->SetName("data_ht_WPTight_27");
	data_ht_tight_27->SetTitle("DATA : HT Distribution for WPTight_27");
	data_npv_tight_27->SetName("data_npv_WPTight_27");
	data_npv_tight_27->SetTitle("DATA : NPV Distribution for WPTight_27");
	data_jet1pt_tight_27->SetName("data_jet1pt_WPTight_27");
	data_jet1pt_tight_27->SetTitle("DATA : Jet1 pT Distribution for WPTight_27");
	data_jet1csv_tight_27->SetName("data_jet1csv_WPTight_27");
	data_jet1csv_tight_27->SetTitle("DATA : Jet1 CSV Distribution for WPTight_27");
	data_njets_tight_27->SetName("data_njets_WPTight_27");
	data_njets_tight_27->SetTitle("DATA : Nr. of Jets Distribution for WPTight_27");
	data_nbtags_tight_27->SetName("data_nbtags_WPTight_27");
	data_nbtags_tight_27->SetTitle("DATA : Nr. of Btags Distribution for WPTight_27");
	data_pt_tight_32->SetName("data_pt_WPTight_32");
	data_pt_tight_32->SetTitle("DATA : Tag Electron pT Distribution for WPTight_32");
	data_pt2_tight_32->SetName("data_pt2_WPTight_32");
	data_pt2_tight_32->SetTitle("DATA : Probe Electron pT Distribution for WPTight_32");
	data_eta_tight_32->SetName("data_eta_WPTight_32");
	data_eta_tight_32->SetTitle("DATA : eta Distribution for WPTight_32");
	data_phi_tight_32->SetName("data_phi_WPTight_32");
	data_phi_tight_32->SetTitle("DATA : phi Distribution for WPTight_32");
	data_ht_tight_32->SetName("data_ht_WPTight_32");
	data_ht_tight_32->SetTitle("DATA : HT Distribution for WPTight_32");
	data_npv_tight_32->SetName("data_npv_WPTight_32");
	data_npv_tight_32->SetTitle("DATA : NPV Distribution for WPTight_32");
	data_jet1pt_tight_32->SetName("data_jet1pt_WPTight_32");
	data_jet1pt_tight_32->SetTitle("DATA : Jet1 pT Distribution for WPTight_32");
	data_jet1csv_tight_32->SetName("data_jet1csv_WPTight_32");
	data_jet1csv_tight_32->SetTitle("DATA : Jet1 CSV Distribution for WPTight_32");
	data_njets_tight_32->SetName("data_njets_WPTight_32");
	data_njets_tight_32->SetTitle("DATA : Nr. of Jets Distribution for WPTight_32");
	data_nbtags_tight_32->SetName("data_nbtags_WPTight_32");
	data_nbtags_tight_32->SetTitle("DATA : Nr. of Btags Distribution for WPTight_32");
	data_pt_loose_27_ht200->SetName("data_pt_WPloose_27_ht200");
	data_pt_loose_27_ht200->SetTitle("DATA : Tag Electron pT Distribution for WPloose_27_ht200");
	data_pt2_loose_27_ht200->SetName("data_pt2_WPloose_27_ht200");
	data_pt2_loose_27_ht200->SetTitle("DATA : Probe ELectron pT Distribution for WPloose_27_ht200");
	data_eta_loose_27_ht200->SetName("data_eta_WPloose_27_ht200");
	data_eta_loose_27_ht200->SetTitle("DATA : eta Distribution for WPloose_27_ht200");
	data_phi_loose_27_ht200->SetName("data_phi_WPloose_27_ht200");
	data_phi_loose_27_ht200->SetTitle("DATA : phi Distribution for WPloose_27_ht200");
	data_ht_loose_27_ht200->SetName("data_ht_WPloose_27_ht200");
	data_ht_loose_27_ht200->SetTitle("DATA : HT Distribution for WPloose_27_ht200");
	data_npv_loose_27_ht200->SetName("data_npv_WPloose_27_ht200");
	data_npv_loose_27_ht200->SetTitle("DATA : NPV Distribution for WPloose_27_ht200");
	data_jet1pt_loose_27_ht200->SetName("data_jet1pt_WPloose_27_ht200");
	data_jet1pt_loose_27_ht200->SetTitle("DATA : Jet1 pT Distribution for WPloose_27_ht200");
	data_jet1csv_loose_27_ht200->SetName("data_jet1csv_WPloose_27_ht200");
	data_jet1csv_loose_27_ht200->SetTitle("DATA : Jet1 CSV Distribution for WPloose_27_ht200");
	data_njets_loose_27_ht200->SetName("data_njets_WPloose_27_ht200");
	data_njets_loose_27_ht200->SetTitle("DATA : Nr. of Jets Distribution for WPloose_27_ht200");
	data_nbtags_loose_27_ht200->SetName("data_nbtags_WPloose_27_ht200");
	data_nbtags_loose_27_ht200->SetTitle("DATA : Nr. of Btags Distribution for WPloose_27_ht200");
	data_pt_tight_27_or_loose_27_ht200->SetName("data_pt_WPtight_27_or_loose_27_ht200");
	data_pt_tight_27_or_loose_27_ht200->SetTitle("DATA : Tag Electron pT Distribution for WPtight_27_or_loose_27_ht200");
	data_pt2_tight_27_or_loose_27_ht200->SetName("data_pt2_WPtight_27_or_loose_27_ht200");
	data_pt2_tight_27_or_loose_27_ht200->SetTitle("DATA : Probe Electron pT Distribution for WPtight_27_or_loose_27_ht200");
	data_eta_tight_27_or_loose_27_ht200->SetName("data_eta_WPtight_27_or_loose_27_ht200");
	data_eta_tight_27_or_loose_27_ht200->SetTitle("DATA : eta Distribution for WPtight_27_or_loose_27_ht200");
	data_phi_tight_27_or_loose_27_ht200->SetName("data_phi_WPtight_27_or_loose_27_ht200");
	data_phi_tight_27_or_loose_27_ht200->SetTitle("DATA : phi Distribution for WPtight_27_or_loose_27_ht200");
	data_ht_tight_27_or_loose_27_ht200->SetName("data_ht_WPtight_27_or_loose_27_ht200");
	data_ht_tight_27_or_loose_27_ht200->SetTitle("DATA : HT Distribution for WPtight_27_or_loose_27_ht200");
	data_npv_tight_27_or_loose_27_ht200->SetName("data_npv_WPtight_27_or_loose_27_ht200");
	data_npv_tight_27_or_loose_27_ht200->SetTitle("DATA : NPV Distribution for WPtight_27_or_loose_27_ht200");
	data_jet1pt_tight_27_or_loose_27_ht200->SetName("data_jet1pt_WPtight_27_or_loose_27_ht200");
	data_jet1pt_tight_27_or_loose_27_ht200->SetTitle("DATA : Jet1 pT Distribution for WPtight_27_or_loose_27_ht200");
	data_jet1csv_tight_27_or_loose_27_ht200->SetName("data_jet1csv_WPtight_27_or_loose_27_ht200");
	data_jet1csv_tight_27_or_loose_27_ht200->SetTitle("DATA : Jet1 CSV Distribution for WPtight_27_or_loose_27_ht200");
	data_njets_tight_27_or_loose_27_ht200->SetName("data_njets_WPtight_27_or_loose_27_ht200");
	data_njets_tight_27_or_loose_27_ht200->SetTitle("DATA : Nr. of Jets Distribution for WPtight_27_or_loose_27_ht200");
	data_nbtags_tight_27_or_loose_27_ht200->SetName("data_nbtags_WPtight_27_or_loose_27_ht200");
	data_nbtags_tight_27_or_loose_27_ht200->SetTitle("DATA : Nr. of Btags Distribution for WPtight_27_or_loose_27_ht200");
	
	
	///////////////////////////////////////////////////////////
	// Scaling and Normlization of mc signal : dy
	///////////////////////////////////////////////////////////
	
	double factor_dy = 1;
	//double sigma_dy = 3*2008.4; // pb
	double sigma_dy = 3*1921.8; // pb
	double N_total_dy = 104113466;
	double norm_dy = (L*sigma_dy*1000)/N_total_dy;
	double scale_dy = factor_dy*norm_dy;
	
	dy_pt_loose_27->Scale(scale_dy);
	dy_pt2_loose_27->Scale(scale_dy);
	dy_eta_loose_27->Scale(scale_dy);
	dy_phi_loose_27->Scale(scale_dy);
	dy_ht_loose_27->Scale(scale_dy);
	dy_npv_loose_27->Scale(scale_dy);
	dy_jet1pt_loose_27->Scale(scale_dy);
	dy_jet1csv_loose_27->Scale(scale_dy);
	dy_njets_loose_27->Scale(scale_dy);
	dy_nbtags_loose_27->Scale(scale_dy);
	dy_pt_tight_27->Scale(scale_dy);
	dy_eta_tight_27->Scale(scale_dy);
	dy_phi_tight_27->Scale(scale_dy);
	dy_ht_tight_27->Scale(scale_dy);
	dy_npv_tight_27->Scale(scale_dy);
	dy_jet1pt_tight_27->Scale(scale_dy);
	dy_jet1csv_tight_27->Scale(scale_dy);
	dy_njets_tight_27->Scale(scale_dy);
	dy_nbtags_tight_27->Scale(scale_dy);
	dy_pt_tight_32->Scale(scale_dy);
	dy_eta_tight_32->Scale(scale_dy);
	dy_phi_tight_32->Scale(scale_dy);
	dy_ht_tight_32->Scale(scale_dy);
	dy_npv_tight_32->Scale(scale_dy);
	dy_jet1pt_tight_32->Scale(scale_dy);
	dy_jet1csv_tight_32->Scale(scale_dy);
	dy_njets_tight_32->Scale(scale_dy);
	dy_nbtags_tight_32->Scale(scale_dy);
	dy_pt_loose_27_ht200->Scale(scale_dy);
	dy_eta_loose_27_ht200->Scale(scale_dy);
	dy_phi_loose_27_ht200->Scale(scale_dy);
	dy_ht_loose_27_ht200->Scale(scale_dy);
	dy_npv_loose_27_ht200->Scale(scale_dy);
	dy_jet1pt_loose_27_ht200->Scale(scale_dy);
	dy_jet1csv_loose_27_ht200->Scale(scale_dy);
	dy_njets_loose_27_ht200->Scale(scale_dy);
	dy_nbtags_loose_27_ht200->Scale(scale_dy);
	dy_pt_tight_27_or_loose_27_ht200->Scale(scale_dy);
	dy_eta_tight_27_or_loose_27_ht200->Scale(scale_dy);
	dy_phi_tight_27_or_loose_27_ht200->Scale(scale_dy);
	dy_ht_tight_27_or_loose_27_ht200->Scale(scale_dy);
	dy_npv_tight_27_or_loose_27_ht200->Scale(scale_dy);
	dy_jet1pt_tight_27_or_loose_27_ht200->Scale(scale_dy);
	dy_jet1csv_tight_27_or_loose_27_ht200->Scale(scale_dy);
	dy_njets_tight_27_or_loose_27_ht200->Scale(scale_dy);
	dy_nbtags_tight_27_or_loose_27_ht200->Scale(scale_dy);
	
	////////////////////////////////////////////////////////
	//Plotting
	////////////////////////////////////////////////////////
	
	TCanvas *c1 = new TCanvas("c1","test",1100,650);
	c1->DrawFrame(0,0,300,300000,"Tight Electron pT (for WPTight_27) ; pT (Electron) [GeV] ; Nr. of Events");
	TLegend* leg1 = new TLegend(0.65,0.70,0.85,0.85);
    leg1->SetFillColor(kWhite);
    leg1->SetFillStyle(1001);
    leg1->AddEntry(dy_pt_tight_27,"MC : DY","L");
    leg1->AddEntry(data_pt_tight_27,"Data","L");	
	dy_pt_tight_27->SetLineColor(kRed);
	dy_pt_tight_27->SetLineWidth(2);
	//dy_pt_tight_27->SetFillStyle(3004);
	//dy_pt_tight_27->SetFillColor(kOrange); 
	dy_pt_tight_27->Draw("same");
	data_pt_tight_27->SetLineColor(kBlue);
	data_pt_tight_27->SetLineWidth(2);
	data_pt_tight_27->Draw("same");
	leg1->Draw("same");
	c1->Print("tightele_pt_tight_27_control.png");
	delete c1;
	delete leg1;
	
	
	TCanvas *c2= new TCanvas("c2","test",1100,650);
	c2->DrawFrame(0,0,300,300000,"Tight Electron pT (for WPTight_32) ; pT (Electron) [GeV] ; Nr. of Events");
	TLegend* leg2 = new TLegend(0.65,0.70,0.85,0.85);
    leg2->SetFillColor(kWhite);
    leg2->SetFillStyle(1001);
    leg2->AddEntry(dy_pt_tight_32,"MC : DY","L");
    leg2->AddEntry(data_pt_tight_32,"Data","L");	
	dy_pt_tight_32->SetLineColor(kRed);
	dy_pt_tight_32->SetLineWidth(2);
	dy_pt_tight_32->Draw("same");
	data_pt_tight_32->SetLineColor(kBlue);
	data_pt_tight_32->SetLineWidth(2);
	data_pt_tight_32->Draw("same");
	leg2->Draw("same");
	c2->Print("tightele_pt_tight_32_control.png");
	delete c2;
	delete leg2;
	
	TCanvas *c3 = new TCanvas("c3","test",1100,650);
	c3->DrawFrame(0,0,300,11000,"Tight Electron pT (for WPLoose_27_HT200) ; pT (Electron) [GeV] ; Nr. of Events");
	TLegend* leg3 = new TLegend(0.65,0.70,0.85,0.85);
    leg3->SetFillColor(kWhite);
    leg3->SetFillStyle(1001);
    leg3->AddEntry(dy_pt_loose_27_ht200,"MC : DY","L");
    leg3->AddEntry(data_pt_loose_27_ht200,"Data","L");	
	dy_pt_loose_27_ht200->SetLineColor(kRed);
	dy_pt_loose_27_ht200->SetLineWidth(2);
	dy_pt_loose_27_ht200->Draw("same");
	data_pt_loose_27_ht200->SetLineColor(kBlue);
	data_pt_loose_27_ht200->SetLineWidth(2);
	data_pt_loose_27_ht200->Draw("same");
	leg3->Draw("same");
	c3->Print("tightele_pt_loose_27_ht200_control.png");
	delete c3;
	delete leg3;
	
	TCanvas *c4= new TCanvas("c4","test",1100,650);
	c4->DrawFrame(0,0,300,300000,"Tight Electron pT (for WPTight_27_OR_WPLoose_27_HT200) ; pT (Electron) [GeV] ; Nr. of Events");
	TLegend* leg4 = new TLegend(0.65,0.70,0.85,0.85);
    leg4->SetFillColor(kWhite);
    leg4->SetFillStyle(1001);
    leg4->AddEntry(dy_pt_tight_27_or_loose_27_ht200,"MC : DY","L");
    leg4->AddEntry(data_pt_tight_27_or_loose_27_ht200,"Data","L");	
	dy_pt_tight_27_or_loose_27_ht200->SetLineColor(kRed);
	dy_pt_tight_27_or_loose_27_ht200->SetLineWidth(2);
	dy_pt_tight_27_or_loose_27_ht200->Draw("same");
	data_pt_tight_27_or_loose_27_ht200->SetLineColor(kBlue);
	data_pt_tight_27_or_loose_27_ht200->SetLineWidth(2);
	data_pt_tight_27_or_loose_27_ht200->Draw("same");
	leg4->Draw("same");
	c4->Print("tightele_pt_tight_27_or_loose_27_ht200_control.png");
	delete c4;
	delete leg4;
	
	TCanvas *c5= new TCanvas("c5","test",1100,650);
	c5->DrawFrame(0,0,300,300000,"Tight Electron pT (for WPLoose_27) ; pT (Electron) [GeV] ; Nr. of Events");
	TLegend* leg5 = new TLegend(0.65,0.70,0.85,0.85);
    leg5->SetFillColor(kWhite);
    leg5->SetFillStyle(1001);
    leg5->AddEntry(dy_pt_loose_27,"MC : DY","L");
    leg5->AddEntry(data_pt_loose_27,"Data","L");	
	dy_pt_loose_27->SetLineColor(kRed);
	dy_pt_loose_27->SetLineWidth(2);
	dy_pt_loose_27->Draw("same");
	data_pt_loose_27->SetLineColor(kBlue);
	data_pt_loose_27->SetLineWidth(2);
	data_pt_loose_27->Draw("same");
	leg5->Draw("same");
	c5->Print("tightele_pt_loose_27_control.png");
	delete c5;
	delete leg5;
	
	
	
	
	
	
	
	
	
	
	
	
	

	delete f1;
	delete f5;
	return;
}


