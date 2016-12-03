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

void combine_data_mc() {
	
	TH1::SetDefaultSumw2();
	
	double L = 26.149; // 1/fb
	
	TFile *f1 = new TFile("Distribution_mc_tthbb.root");
	TH1D *tthbb_pt_tight_27 = (TH1D*)f1->Get("Pt_WPTight_27");
	TH1D *tthbb_eta_tight_27 = (TH1D*)f1->Get("Eta_WPTight_27");
	TH1D *tthbb_phi_tight_27 = (TH1D*)f1->Get("Phi_WPTight_27");
	TH1D *tthbb_ht_tight_27 = (TH1D*)f1->Get("HT_WPTight_27");
	TH1D *tthbb_npv_tight_27 = (TH1D*)f1->Get("NumPV_WPTight_27");
	TH1D *tthbb_jet1pt_tight_27 = (TH1D*)f1->Get("Jet1_pt_WPTight_27");
	TH1D *tthbb_jet1csv_tight_27 = (TH1D*)f1->Get("Jet1_csv_WPTight_27");
	TH1D *tthbb_njets_tight_27 = (TH1D*)f1->Get("Njets_WPTight_27");
	TH1D *tthbb_nbtags_tight_27 = (TH1D*)f1->Get("Nbtags_WPTight_27");
	TH1D *tthbb_pt_tight_32 = (TH1D*)f1->Get("Pt_WPTight_32");
	TH1D *tthbb_eta_tight_32 = (TH1D*)f1->Get("Eta_WPTight_32");
	TH1D *tthbb_phi_tight_32 = (TH1D*)f1->Get("Phi_WPTight_32");
	TH1D *tthbb_ht_tight_32 = (TH1D*)f1->Get("HT_WPTight_32");
	TH1D *tthbb_npv_tight_32 = (TH1D*)f1->Get("NumPV_WPTight_32");
	TH1D *tthbb_jet1pt_tight_32 = (TH1D*)f1->Get("Jet1_pt_WPTight_32");
	TH1D *tthbb_jet1csv_tight_32 = (TH1D*)f1->Get("Jet1_csv_WPTight_32");
	TH1D *tthbb_njets_tight_32 = (TH1D*)f1->Get("Njets_WPTight_32");
	TH1D *tthbb_nbtags_tight_32 = (TH1D*)f1->Get("Nbtags_WPTight_32");
	TH1D *tthbb_pt_loose_27_ht200 = (TH1D*)f1->Get("Pt_WPLoose_27_HT200");
	TH1D *tthbb_eta_loose_27_ht200 = (TH1D*)f1->Get("Eta_WPLoose_27_HT200");
	TH1D *tthbb_phi_loose_27_ht200 = (TH1D*)f1->Get("Phi_WPLoose_27_HT200");
	TH1D *tthbb_ht_loose_27_ht200 = (TH1D*)f1->Get("HT_WPLoose_27_HT200");
	TH1D *tthbb_npv_loose_27_ht200 = (TH1D*)f1->Get("NumPV_WPLoose_27_HT200");
	TH1D *tthbb_jet1pt_loose_27_ht200 = (TH1D*)f1->Get("Jet1_pt_WPLoose_27_HT200");
	TH1D *tthbb_jet1csv_loose_27_ht200 = (TH1D*)f1->Get("Jet1_csv_WPLoose_27_HT200");
	TH1D *tthbb_njets_loose_27_ht200 = (TH1D*)f1->Get("Njets_WPLoose_27_HT200");
	TH1D *tthbb_nbtags_loose_27_ht200 = (TH1D*)f1->Get("Nbtags_WPLoose_27_HT200");
	TH1D *tthbb_pt_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthbb_eta_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Eta_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthbb_phi_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Phi_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthbb_ht_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("HT_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthbb_npv_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("NumPV_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthbb_jet1pt_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Jet1_pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthbb_jet1csv_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Jet1_csv_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthbb_njets_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Njets_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthbb_nbtags_tight_27_or_loose_27_ht200 = (TH1D*)f1->Get("Nbtags_WPTight_27_OR_WPLoose_27_HT200");
	
	tthbb_pt_tight_27->SetName("ttHbb_pt_WPTight_27");
	tthbb_pt_tight_27->SetTitle("MC ttHbb : pT Distribution for WPTight_27");
	tthbb_eta_tight_27->SetName("ttHbb_eta_WPTight_27");
	tthbb_eta_tight_27->SetTitle("MC ttHbb : eta Distribution for WPTight_27");
	tthbb_phi_tight_27->SetName("ttHbb_phi_WPTight_27");
	tthbb_phi_tight_27->SetTitle("MC ttHbb : phi Distribution for WPTight_27");
	tthbb_ht_tight_27->SetName("ttHbb_ht_WPTight_27");
	tthbb_ht_tight_27->SetTitle("MC ttHbb : HT Distribution for WPTight_27");
	tthbb_npv_tight_27->SetName("ttHbb_npv_WPTight_27");
	tthbb_npv_tight_27->SetTitle("MC ttHbb : NPV Distribution for WPTight_27");
	tthbb_jet1pt_tight_27->SetName("ttHbb_jet1pt_WPTight_27");
	tthbb_jet1pt_tight_27->SetTitle("MC ttHbb : Jet1 pT Distribution for WPTight_27");
	tthbb_jet1csv_tight_27->SetName("ttHbb_jet1csv_WPTight_27");
	tthbb_jet1csv_tight_27->SetTitle("MC ttHbb : Jet1 CSV Distribution for WPTight_27");
	tthbb_njets_tight_27->SetName("ttHbb_njets_WPTight_27");
	tthbb_njets_tight_27->SetTitle("MC ttHbb : Nr. of Jets Distribution for WPTight_27");
	tthbb_nbtags_tight_27->SetName("ttHbb_nbtags_WPTight_27");
	tthbb_nbtags_tight_27->SetTitle("MC ttHbb : Nr. of Btags Distribution for WPTight_27");
	tthbb_pt_tight_32->SetName("ttHbb_pt_WPTight_32");
	tthbb_pt_tight_32->SetTitle("MC ttHbb : pT Distribution for WPTight_32");
	tthbb_eta_tight_32->SetName("ttHbb_eta_WPTight_32");
	tthbb_eta_tight_32->SetTitle("MC ttHbb : eta Distribution for WPTight_32");
	tthbb_phi_tight_32->SetName("ttHbb_phi_WPTight_32");
	tthbb_phi_tight_32->SetTitle("MC ttHbb : phi Distribution for WPTight_32");
	tthbb_ht_tight_32->SetName("ttHbb_ht_WPTight_32");
	tthbb_ht_tight_32->SetTitle("MC ttHbb : HT Distribution for WPTight_32");
	tthbb_npv_tight_32->SetName("ttHbb_npv_WPTight_32");
	tthbb_npv_tight_32->SetTitle("MC ttHbb : NPV Distribution for WPTight_32");
	tthbb_jet1pt_tight_32->SetName("ttHbb_jet1pt_WPTight_32");
	tthbb_jet1pt_tight_32->SetTitle("MC ttHbb : Jet1 pT Distribution for WPTight_32");
	tthbb_jet1csv_tight_32->SetName("ttHbb_jet1csv_WPTight_32");
	tthbb_jet1csv_tight_32->SetTitle("MC ttHbb : Jet1 CSV Distribution for WPTight_32");
	tthbb_njets_tight_32->SetName("ttHbb_njets_WPTight_32");
	tthbb_njets_tight_32->SetTitle("MC ttHbb : Nr. of Jets Distribution for WPTight_32");
	tthbb_nbtags_tight_32->SetName("ttHbb_nbtags_WPTight_32");
	tthbb_nbtags_tight_32->SetTitle("MC ttHbb : Nr. of Btags Distribution for WPTight_32");
	tthbb_pt_loose_27_ht200->SetName("ttHbb_pt_WPloose_27_ht200");
	tthbb_pt_loose_27_ht200->SetTitle("MC ttHbb : pT Distribution for WPloose_27_ht200");
	tthbb_eta_loose_27_ht200->SetName("ttHbb_eta_WPloose_27_ht200");
	tthbb_eta_loose_27_ht200->SetTitle("MC ttHbb : eta Distribution for WPloose_27_ht200");
	tthbb_phi_loose_27_ht200->SetName("ttHbb_phi_WPloose_27_ht200");
	tthbb_phi_loose_27_ht200->SetTitle("MC ttHbb : phi Distribution for WPloose_27_ht200");
	tthbb_ht_loose_27_ht200->SetName("ttHbb_ht_WPloose_27_ht200");
	tthbb_ht_loose_27_ht200->SetTitle("MC ttHbb : HT Distribution for WPloose_27_ht200");
	tthbb_npv_loose_27_ht200->SetName("ttHbb_npv_WPloose_27_ht200");
	tthbb_npv_loose_27_ht200->SetTitle("MC ttHbb : NPV Distribution for WPloose_27_ht200");
	tthbb_jet1pt_loose_27_ht200->SetName("ttHbb_jet1pt_WPloose_27_ht200");
	tthbb_jet1pt_loose_27_ht200->SetTitle("MC ttHbb : Jet1 pT Distribution for WPloose_27_ht200");
	tthbb_jet1csv_loose_27_ht200->SetName("ttHbb_jet1csv_WPloose_27_ht200");
	tthbb_jet1csv_loose_27_ht200->SetTitle("MC ttHbb : Jet1 CSV Distribution for WPloose_27_ht200");
	tthbb_njets_loose_27_ht200->SetName("ttHbb_njets_WPloose_27_ht200");
	tthbb_njets_loose_27_ht200->SetTitle("MC ttHbb : Nr. of Jets Distribution for WPloose_27_ht200");
	tthbb_nbtags_loose_27_ht200->SetName("ttHbb_nbtags_WPloose_27_ht200");
	tthbb_nbtags_loose_27_ht200->SetTitle("MC ttHbb : Nr. of Btags Distribution for WPloose_27_ht200");
	tthbb_pt_tight_27_or_loose_27_ht200->SetName("ttHbb_pt_WPtight_27_or_loose_27_ht200");
	tthbb_pt_tight_27_or_loose_27_ht200->SetTitle("MC ttHbb : pT Distribution for WPtight_27_or_loose_27_ht200");
	tthbb_eta_tight_27_or_loose_27_ht200->SetName("ttHbb_eta_WPtight_27_or_loose_27_ht200");
	tthbb_eta_tight_27_or_loose_27_ht200->SetTitle("MC ttHbb : eta Distribution for WPtight_27_or_loose_27_ht200");
	tthbb_phi_tight_27_or_loose_27_ht200->SetName("ttHbb_phi_WPtight_27_or_loose_27_ht200");
	tthbb_phi_tight_27_or_loose_27_ht200->SetTitle("MC ttHbb : phi Distribution for WPtight_27_or_loose_27_ht200");
	tthbb_ht_tight_27_or_loose_27_ht200->SetName("ttHbb_ht_WPtight_27_or_loose_27_ht200");
	tthbb_ht_tight_27_or_loose_27_ht200->SetTitle("MC ttHbb : HT Distribution for WPtight_27_or_loose_27_ht200");
	tthbb_npv_tight_27_or_loose_27_ht200->SetName("ttHbb_npv_WPtight_27_or_loose_27_ht200");
	tthbb_npv_tight_27_or_loose_27_ht200->SetTitle("MC ttHbb : NPV Distribution for WPtight_27_or_loose_27_ht200");
	tthbb_jet1pt_tight_27_or_loose_27_ht200->SetName("ttHbb_jet1pt_WPtight_27_or_loose_27_ht200");
	tthbb_jet1pt_tight_27_or_loose_27_ht200->SetTitle("MC ttHbb : Jet1 pT Distribution for WPtight_27_or_loose_27_ht200");
	tthbb_jet1csv_tight_27_or_loose_27_ht200->SetName("ttHbb_jet1csv_WPtight_27_or_loose_27_ht200");
	tthbb_jet1csv_tight_27_or_loose_27_ht200->SetTitle("MC ttHbb : Jet1 CSV Distribution for WPtight_27_or_loose_27_ht200");
	tthbb_njets_tight_27_or_loose_27_ht200->SetName("ttHbb_njets_WPtight_27_or_loose_27_ht200");
	tthbb_njets_tight_27_or_loose_27_ht200->SetTitle("MC ttHbb : Nr. of Jets Distribution for WPtight_27_or_loose_27_ht200");
	tthbb_nbtags_tight_27_or_loose_27_ht200->SetName("ttHbb_nbtags_WPtight_27_or_loose_27_ht200");
	tthbb_nbtags_tight_27_or_loose_27_ht200->SetTitle("MC ttHbb : Nr. of Btags Distribution for WPtight_27_or_loose_27_ht200");
	
	
	TFile *f2 = new TFile("Distribution_mc_tthnonbb.root");
	TH1D *tthnonbb_pt_tight_27 = (TH1D*)f2->Get("Pt_WPTight_27");
	TH1D *tthnonbb_eta_tight_27 = (TH1D*)f2->Get("Eta_WPTight_27");
	TH1D *tthnonbb_phi_tight_27 = (TH1D*)f2->Get("Phi_WPTight_27");
	TH1D *tthnonbb_ht_tight_27 = (TH1D*)f2->Get("HT_WPTight_27");
	TH1D *tthnonbb_npv_tight_27 = (TH1D*)f2->Get("NumPV_WPTight_27");
	TH1D *tthnonbb_jet1pt_tight_27 = (TH1D*)f2->Get("Jet1_pt_WPTight_27");
	TH1D *tthnonbb_jet1csv_tight_27 = (TH1D*)f2->Get("Jet1_csv_WPTight_27");
	TH1D *tthnonbb_njets_tight_27 = (TH1D*)f2->Get("Njets_WPTight_27");
	TH1D *tthnonbb_nbtags_tight_27 = (TH1D*)f2->Get("Nbtags_WPTight_27");
	TH1D *tthnonbb_pt_tight_32 = (TH1D*)f2->Get("Pt_WPTight_32");
	TH1D *tthnonbb_eta_tight_32 = (TH1D*)f2->Get("Eta_WPTight_32");
	TH1D *tthnonbb_phi_tight_32 = (TH1D*)f2->Get("Phi_WPTight_32");
	TH1D *tthnonbb_ht_tight_32 = (TH1D*)f2->Get("HT_WPTight_32");
	TH1D *tthnonbb_npv_tight_32 = (TH1D*)f2->Get("NumPV_WPTight_32");
	TH1D *tthnonbb_jet1pt_tight_32 = (TH1D*)f2->Get("Jet1_pt_WPTight_32");
	TH1D *tthnonbb_jet1csv_tight_32 = (TH1D*)f2->Get("Jet1_csv_WPTight_32");
	TH1D *tthnonbb_njets_tight_32 = (TH1D*)f2->Get("Njets_WPTight_32");
	TH1D *tthnonbb_nbtags_tight_32 = (TH1D*)f2->Get("Nbtags_WPTight_32");
	TH1D *tthnonbb_pt_loose_27_ht200 = (TH1D*)f2->Get("Pt_WPLoose_27_HT200");
	TH1D *tthnonbb_eta_loose_27_ht200 = (TH1D*)f2->Get("Eta_WPLoose_27_HT200");
	TH1D *tthnonbb_phi_loose_27_ht200 = (TH1D*)f2->Get("Phi_WPLoose_27_HT200");
	TH1D *tthnonbb_ht_loose_27_ht200 = (TH1D*)f2->Get("HT_WPLoose_27_HT200");
	TH1D *tthnonbb_npv_loose_27_ht200 = (TH1D*)f2->Get("NumPV_WPLoose_27_HT200");
	TH1D *tthnonbb_jet1pt_loose_27_ht200 = (TH1D*)f2->Get("Jet1_pt_WPLoose_27_HT200");
	TH1D *tthnonbb_jet1csv_loose_27_ht200 = (TH1D*)f2->Get("Jet1_csv_WPLoose_27_HT200");
	TH1D *tthnonbb_njets_loose_27_ht200 = (TH1D*)f2->Get("Njets_WPLoose_27_HT200");
	TH1D *tthnonbb_nbtags_loose_27_ht200 = (TH1D*)f2->Get("Nbtags_WPLoose_27_HT200");
	TH1D *tthnonbb_pt_tight_27_or_loose_27_ht200 = (TH1D*)f2->Get("Pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthnonbb_eta_tight_27_or_loose_27_ht200 = (TH1D*)f2->Get("Eta_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthnonbb_phi_tight_27_or_loose_27_ht200 = (TH1D*)f2->Get("Phi_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthnonbb_ht_tight_27_or_loose_27_ht200 = (TH1D*)f2->Get("HT_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthnonbb_npv_tight_27_or_loose_27_ht200 = (TH1D*)f2->Get("NumPV_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthnonbb_jet1pt_tight_27_or_loose_27_ht200 = (TH1D*)f2->Get("Jet1_pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthnonbb_jet1csv_tight_27_or_loose_27_ht200 = (TH1D*)f2->Get("Jet1_csv_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthnonbb_njets_tight_27_or_loose_27_ht200 = (TH1D*)f2->Get("Njets_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *tthnonbb_nbtags_tight_27_or_loose_27_ht200 = (TH1D*)f2->Get("Nbtags_WPTight_27_OR_WPLoose_27_HT200");
	
	tthnonbb_pt_tight_27->SetName("tthnonbb_pt_WPTight_27");
	tthnonbb_pt_tight_27->SetTitle("MC tthnonbb : pT Distribution for WPTight_27");
	tthnonbb_eta_tight_27->SetName("tthnonbb_eta_WPTight_27");
	tthnonbb_eta_tight_27->SetTitle("MC tthnonbb : eta Distribution for WPTight_27");
	tthnonbb_phi_tight_27->SetName("tthnonbb_phi_WPTight_27");
	tthnonbb_phi_tight_27->SetTitle("MC tthnonbb : phi Distribution for WPTight_27");
	tthnonbb_ht_tight_27->SetName("tthnonbb_ht_WPTight_27");
	tthnonbb_ht_tight_27->SetTitle("MC tthnonbb : HT Distribution for WPTight_27");
	tthnonbb_npv_tight_27->SetName("tthnonbb_npv_WPTight_27");
	tthnonbb_npv_tight_27->SetTitle("MC tthnonbb : NPV Distribution for WPTight_27");
	tthnonbb_jet1pt_tight_27->SetName("tthnonbb_jet1pt_WPTight_27");
	tthnonbb_jet1pt_tight_27->SetTitle("MC tthnonbb : Jet1 pT Distribution for WPTight_27");
	tthnonbb_jet1csv_tight_27->SetName("tthnonbb_jet1csv_WPTight_27");
	tthnonbb_jet1csv_tight_27->SetTitle("MC tthnonbb : Jet1 CSV Distribution for WPTight_27");
	tthnonbb_njets_tight_27->SetName("tthnonbb_njets_WPTight_27");
	tthnonbb_njets_tight_27->SetTitle("MC tthnonbb : Nr. of Jets Distribution for WPTight_27");
	tthnonbb_nbtags_tight_27->SetName("tthnonbb_nbtags_WPTight_27");
	tthnonbb_nbtags_tight_27->SetTitle("MC tthnonbb : Nr. of Btags Distribution for WPTight_27");
	tthnonbb_pt_tight_32->SetName("tthnonbb_pt_WPTight_32");
	tthnonbb_pt_tight_32->SetTitle("MC tthnonbb : pT Distribution for WPTight_32");
	tthnonbb_eta_tight_32->SetName("tthnonbb_eta_WPTight_32");
	tthnonbb_eta_tight_32->SetTitle("MC tthnonbb : eta Distribution for WPTight_32");
	tthnonbb_phi_tight_32->SetName("tthnonbb_phi_WPTight_32");
	tthnonbb_phi_tight_32->SetTitle("MC tthnonbb : phi Distribution for WPTight_32");
	tthnonbb_ht_tight_32->SetName("tthnonbb_ht_WPTight_32");
	tthnonbb_ht_tight_32->SetTitle("MC tthnonbb : HT Distribution for WPTight_32");
	tthnonbb_npv_tight_32->SetName("tthnonbb_npv_WPTight_32");
	tthnonbb_npv_tight_32->SetTitle("MC tthnonbb : NPV Distribution for WPTight_32");
	tthnonbb_jet1pt_tight_32->SetName("tthnonbb_jet1pt_WPTight_32");
	tthnonbb_jet1pt_tight_32->SetTitle("MC tthnonbb : Jet1 pT Distribution for WPTight_32");
	tthnonbb_jet1csv_tight_32->SetName("tthnonbb_jet1csv_WPTight_32");
	tthnonbb_jet1csv_tight_32->SetTitle("MC tthnonbb : Jet1 CSV Distribution for WPTight_32");
	tthnonbb_njets_tight_32->SetName("tthnonbb_njets_WPTight_32");
	tthnonbb_njets_tight_32->SetTitle("MC tthnonbb : Nr. of Jets Distribution for WPTight_32");
	tthnonbb_nbtags_tight_32->SetName("tthnonbb_nbtags_WPTight_32");
	tthnonbb_nbtags_tight_32->SetTitle("MC tthnonbb : Nr. of Btags Distribution for WPTight_32");
	tthnonbb_pt_loose_27_ht200->SetName("tthnonbb_pt_WPloose_27_ht200");
	tthnonbb_pt_loose_27_ht200->SetTitle("MC tthnonbb : pT Distribution for WPloose_27_ht200");
	tthnonbb_eta_loose_27_ht200->SetName("tthnonbb_eta_WPloose_27_ht200");
	tthnonbb_eta_loose_27_ht200->SetTitle("MC tthnonbb : eta Distribution for WPloose_27_ht200");
	tthnonbb_phi_loose_27_ht200->SetName("tthnonbb_phi_WPloose_27_ht200");
	tthnonbb_phi_loose_27_ht200->SetTitle("MC tthnonbb : phi Distribution for WPloose_27_ht200");
	tthnonbb_ht_loose_27_ht200->SetName("tthnonbb_ht_WPloose_27_ht200");
	tthnonbb_ht_loose_27_ht200->SetTitle("MC tthnonbb : HT Distribution for WPloose_27_ht200");
	tthnonbb_npv_loose_27_ht200->SetName("tthnonbb_npv_WPloose_27_ht200");
	tthnonbb_npv_loose_27_ht200->SetTitle("MC tthnonbb : NPV Distribution for WPloose_27_ht200");
	tthnonbb_jet1pt_loose_27_ht200->SetName("tthnonbb_jet1pt_WPloose_27_ht200");
	tthnonbb_jet1pt_loose_27_ht200->SetTitle("MC tthnonbb : Jet1 pT Distribution for WPloose_27_ht200");
	tthnonbb_jet1csv_loose_27_ht200->SetName("tthnonbb_jet1csv_WPloose_27_ht200");
	tthnonbb_jet1csv_loose_27_ht200->SetTitle("MC tthnonbb : Jet1 CSV Distribution for WPloose_27_ht200");
	tthnonbb_njets_loose_27_ht200->SetName("tthnonbb_njets_WPloose_27_ht200");
	tthnonbb_njets_loose_27_ht200->SetTitle("MC tthnonbb : Nr. of Jets Distribution for WPloose_27_ht200");
	tthnonbb_nbtags_loose_27_ht200->SetName("tthnonbb_nbtags_WPloose_27_ht200");
	tthnonbb_nbtags_loose_27_ht200->SetTitle("MC tthnonbb : Nr. of Btags Distribution for WPloose_27_ht200");
	tthnonbb_pt_tight_27_or_loose_27_ht200->SetName("tthnonbb_pt_WPtight_27_or_loose_27_ht200");
	tthnonbb_pt_tight_27_or_loose_27_ht200->SetTitle("MC tthnonbb : pT Distribution for WPtight_27_or_loose_27_ht200");
	tthnonbb_eta_tight_27_or_loose_27_ht200->SetName("tthnonbb_eta_WPtight_27_or_loose_27_ht200");
	tthnonbb_eta_tight_27_or_loose_27_ht200->SetTitle("MC tthnonbb : eta Distribution for WPtight_27_or_loose_27_ht200");
	tthnonbb_phi_tight_27_or_loose_27_ht200->SetName("tthnonbb_phi_WPtight_27_or_loose_27_ht200");
	tthnonbb_phi_tight_27_or_loose_27_ht200->SetTitle("MC tthnonbb : phi Distribution for WPtight_27_or_loose_27_ht200");
	tthnonbb_ht_tight_27_or_loose_27_ht200->SetName("tthnonbb_ht_WPtight_27_or_loose_27_ht200");
	tthnonbb_ht_tight_27_or_loose_27_ht200->SetTitle("MC tthnonbb : HT Distribution for WPtight_27_or_loose_27_ht200");
	tthnonbb_npv_tight_27_or_loose_27_ht200->SetName("tthnonbb_npv_WPtight_27_or_loose_27_ht200");
	tthnonbb_npv_tight_27_or_loose_27_ht200->SetTitle("MC tthnonbb : NPV Distribution for WPtight_27_or_loose_27_ht200");
	tthnonbb_jet1pt_tight_27_or_loose_27_ht200->SetName("tthnonbb_jet1pt_WPtight_27_or_loose_27_ht200");
	tthnonbb_jet1pt_tight_27_or_loose_27_ht200->SetTitle("MC tthnonbb : Jet1 pT Distribution for WPtight_27_or_loose_27_ht200");
	tthnonbb_jet1csv_tight_27_or_loose_27_ht200->SetName("tthnonbb_jet1csv_WPtight_27_or_loose_27_ht200");
	tthnonbb_jet1csv_tight_27_or_loose_27_ht200->SetTitle("MC tthnonbb : Jet1 CSV Distribution for WPtight_27_or_loose_27_ht200");
	tthnonbb_njets_tight_27_or_loose_27_ht200->SetName("tthnonbb_njets_WPtight_27_or_loose_27_ht200");
	tthnonbb_njets_tight_27_or_loose_27_ht200->SetTitle("MC tthnonbb : Nr. of Jets Distribution for WPtight_27_or_loose_27_ht200");
	tthnonbb_nbtags_tight_27_or_loose_27_ht200->SetName("tthnonbb_nbtags_WPtight_27_or_loose_27_ht200");
	tthnonbb_nbtags_tight_27_or_loose_27_ht200->SetTitle("MC tthnonbb : Nr. of Btags Distribution for WPtight_27_or_loose_27_ht200");
	
	TFile *f3 = new TFile("Distribution_mc_ttjets_sl.root");
	TH1D *ttjets_sl_pt_tight_27 = (TH1D*)f3->Get("Pt_WPTight_27");
	TH1D *ttjets_sl_eta_tight_27 = (TH1D*)f3->Get("Eta_WPTight_27");
	TH1D *ttjets_sl_phi_tight_27 = (TH1D*)f3->Get("Phi_WPTight_27");
	TH1D *ttjets_sl_ht_tight_27 = (TH1D*)f3->Get("HT_WPTight_27");
	TH1D *ttjets_sl_npv_tight_27 = (TH1D*)f3->Get("NumPV_WPTight_27");
	TH1D *ttjets_sl_jet1pt_tight_27 = (TH1D*)f3->Get("Jet1_pt_WPTight_27");
	TH1D *ttjets_sl_jet1csv_tight_27 = (TH1D*)f3->Get("Jet1_csv_WPTight_27");
	TH1D *ttjets_sl_njets_tight_27 = (TH1D*)f3->Get("Njets_WPTight_27");
	TH1D *ttjets_sl_nbtags_tight_27 = (TH1D*)f3->Get("Nbtags_WPTight_27");
	TH1D *ttjets_sl_pt_tight_32 = (TH1D*)f3->Get("Pt_WPTight_32");
	TH1D *ttjets_sl_eta_tight_32 = (TH1D*)f3->Get("Eta_WPTight_32");
	TH1D *ttjets_sl_phi_tight_32 = (TH1D*)f3->Get("Phi_WPTight_32");
	TH1D *ttjets_sl_ht_tight_32 = (TH1D*)f3->Get("HT_WPTight_32");
	TH1D *ttjets_sl_npv_tight_32 = (TH1D*)f3->Get("NumPV_WPTight_32");
	TH1D *ttjets_sl_jet1pt_tight_32 = (TH1D*)f3->Get("Jet1_pt_WPTight_32");
	TH1D *ttjets_sl_jet1csv_tight_32 = (TH1D*)f3->Get("Jet1_csv_WPTight_32");
	TH1D *ttjets_sl_njets_tight_32 = (TH1D*)f3->Get("Njets_WPTight_32");
	TH1D *ttjets_sl_nbtags_tight_32 = (TH1D*)f3->Get("Nbtags_WPTight_32");
	TH1D *ttjets_sl_pt_loose_27_ht200 = (TH1D*)f3->Get("Pt_WPLoose_27_HT200");
	TH1D *ttjets_sl_eta_loose_27_ht200 = (TH1D*)f3->Get("Eta_WPLoose_27_HT200");
	TH1D *ttjets_sl_phi_loose_27_ht200 = (TH1D*)f3->Get("Phi_WPLoose_27_HT200");
	TH1D *ttjets_sl_ht_loose_27_ht200 = (TH1D*)f3->Get("HT_WPLoose_27_HT200");
	TH1D *ttjets_sl_npv_loose_27_ht200 = (TH1D*)f3->Get("NumPV_WPLoose_27_HT200");
	TH1D *ttjets_sl_jet1pt_loose_27_ht200 = (TH1D*)f3->Get("Jet1_pt_WPLoose_27_HT200");
	TH1D *ttjets_sl_jet1csv_loose_27_ht200 = (TH1D*)f3->Get("Jet1_csv_WPLoose_27_HT200");
	TH1D *ttjets_sl_njets_loose_27_ht200 = (TH1D*)f3->Get("Njets_WPLoose_27_HT200");
	TH1D *ttjets_sl_nbtags_loose_27_ht200 = (TH1D*)f3->Get("Nbtags_WPLoose_27_HT200");
	TH1D *ttjets_sl_pt_tight_27_or_loose_27_ht200 = (TH1D*)f3->Get("Pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_sl_eta_tight_27_or_loose_27_ht200 = (TH1D*)f3->Get("Eta_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_sl_phi_tight_27_or_loose_27_ht200 = (TH1D*)f3->Get("Phi_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_sl_ht_tight_27_or_loose_27_ht200 = (TH1D*)f3->Get("HT_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_sl_npv_tight_27_or_loose_27_ht200 = (TH1D*)f3->Get("NumPV_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_sl_jet1pt_tight_27_or_loose_27_ht200 = (TH1D*)f3->Get("Jet1_pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_sl_jet1csv_tight_27_or_loose_27_ht200 = (TH1D*)f3->Get("Jet1_csv_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_sl_njets_tight_27_or_loose_27_ht200 = (TH1D*)f3->Get("Njets_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_sl_nbtags_tight_27_or_loose_27_ht200 = (TH1D*)f3->Get("Nbtags_WPTight_27_OR_WPLoose_27_HT200");
	
	ttjets_sl_pt_tight_27->SetName("ttjets_sl_pt_WPTight_27");
	ttjets_sl_pt_tight_27->SetTitle("MC ttjets_sl : pT Distribution for WPTight_27");
	ttjets_sl_eta_tight_27->SetName("ttjets_sl_eta_WPTight_27");
	ttjets_sl_eta_tight_27->SetTitle("MC ttjets_sl : eta Distribution for WPTight_27");
	ttjets_sl_phi_tight_27->SetName("ttjets_sl_phi_WPTight_27");
	ttjets_sl_phi_tight_27->SetTitle("MC ttjets_sl : phi Distribution for WPTight_27");
	ttjets_sl_ht_tight_27->SetName("ttjets_sl_ht_WPTight_27");
	ttjets_sl_ht_tight_27->SetTitle("MC ttjets_sl : HT Distribution for WPTight_27");
	ttjets_sl_npv_tight_27->SetName("ttjets_sl_npv_WPTight_27");
	ttjets_sl_npv_tight_27->SetTitle("MC ttjets_sl : NPV Distribution for WPTight_27");
	ttjets_sl_jet1pt_tight_27->SetName("ttjets_sl_jet1pt_WPTight_27");
	ttjets_sl_jet1pt_tight_27->SetTitle("MC ttjets_sl : Jet1 pT Distribution for WPTight_27");
	ttjets_sl_jet1csv_tight_27->SetName("ttjets_sl_jet1csv_WPTight_27");
	ttjets_sl_jet1csv_tight_27->SetTitle("MC ttjets_sl : Jet1 CSV Distribution for WPTight_27");
	ttjets_sl_njets_tight_27->SetName("ttjets_sl_njets_WPTight_27");
	ttjets_sl_njets_tight_27->SetTitle("MC ttjets_sl : Nr. of Jets Distribution for WPTight_27");
	ttjets_sl_nbtags_tight_27->SetName("ttjets_sl_nbtags_WPTight_27");
	ttjets_sl_nbtags_tight_27->SetTitle("MC ttjets_sl : Nr. of Btags Distribution for WPTight_27");
	ttjets_sl_pt_tight_32->SetName("ttjets_sl_pt_WPTight_32");
	ttjets_sl_pt_tight_32->SetTitle("MC ttjets_sl : pT Distribution for WPTight_32");
	ttjets_sl_eta_tight_32->SetName("ttjets_sl_eta_WPTight_32");
	ttjets_sl_eta_tight_32->SetTitle("MC ttjets_sl : eta Distribution for WPTight_32");
	ttjets_sl_phi_tight_32->SetName("ttjets_sl_phi_WPTight_32");
	ttjets_sl_phi_tight_32->SetTitle("MC ttjets_sl : phi Distribution for WPTight_32");
	ttjets_sl_ht_tight_32->SetName("ttjets_sl_ht_WPTight_32");
	ttjets_sl_ht_tight_32->SetTitle("MC ttjets_sl : HT Distribution for WPTight_32");
	ttjets_sl_npv_tight_32->SetName("ttjets_sl_npv_WPTight_32");
	ttjets_sl_npv_tight_32->SetTitle("MC ttjets_sl : NPV Distribution for WPTight_32");
	ttjets_sl_jet1pt_tight_32->SetName("ttjets_sl_jet1pt_WPTight_32");
	ttjets_sl_jet1pt_tight_32->SetTitle("MC ttjets_sl : Jet1 pT Distribution for WPTight_32");
	ttjets_sl_jet1csv_tight_32->SetName("ttjets_sl_jet1csv_WPTight_32");
	ttjets_sl_jet1csv_tight_32->SetTitle("MC ttjets_sl : Jet1 CSV Distribution for WPTight_32");
	ttjets_sl_njets_tight_32->SetName("ttjets_sl_njets_WPTight_32");
	ttjets_sl_njets_tight_32->SetTitle("MC ttjets_sl : Nr. of Jets Distribution for WPTight_32");
	ttjets_sl_nbtags_tight_32->SetName("ttjets_sl_nbtags_WPTight_32");
	ttjets_sl_nbtags_tight_32->SetTitle("MC ttjets_sl : Nr. of Btags Distribution for WPTight_32");
	ttjets_sl_pt_loose_27_ht200->SetName("ttjets_sl_pt_WPloose_27_ht200");
	ttjets_sl_pt_loose_27_ht200->SetTitle("MC ttjets_sl : pT Distribution for WPloose_27_ht200");
	ttjets_sl_eta_loose_27_ht200->SetName("ttjets_sl_eta_WPloose_27_ht200");
	ttjets_sl_eta_loose_27_ht200->SetTitle("MC ttjets_sl : eta Distribution for WPloose_27_ht200");
	ttjets_sl_phi_loose_27_ht200->SetName("ttjets_sl_phi_WPloose_27_ht200");
	ttjets_sl_phi_loose_27_ht200->SetTitle("MC ttjets_sl : phi Distribution for WPloose_27_ht200");
	ttjets_sl_ht_loose_27_ht200->SetName("ttjets_sl_ht_WPloose_27_ht200");
	ttjets_sl_ht_loose_27_ht200->SetTitle("MC ttjets_sl : HT Distribution for WPloose_27_ht200");
	ttjets_sl_npv_loose_27_ht200->SetName("ttjets_sl_npv_WPloose_27_ht200");
	ttjets_sl_npv_loose_27_ht200->SetTitle("MC ttjets_sl : NPV Distribution for WPloose_27_ht200");
	ttjets_sl_jet1pt_loose_27_ht200->SetName("ttjets_sl_jet1pt_WPloose_27_ht200");
	ttjets_sl_jet1pt_loose_27_ht200->SetTitle("MC ttjets_sl : Jet1 pT Distribution for WPloose_27_ht200");
	ttjets_sl_jet1csv_loose_27_ht200->SetName("ttjets_sl_jet1csv_WPloose_27_ht200");
	ttjets_sl_jet1csv_loose_27_ht200->SetTitle("MC ttjets_sl : Jet1 CSV Distribution for WPloose_27_ht200");
	ttjets_sl_njets_loose_27_ht200->SetName("ttjets_sl_njets_WPloose_27_ht200");
	ttjets_sl_njets_loose_27_ht200->SetTitle("MC ttjets_sl : Nr. of Jets Distribution for WPloose_27_ht200");
	ttjets_sl_nbtags_loose_27_ht200->SetName("ttjets_sl_nbtags_WPloose_27_ht200");
	ttjets_sl_nbtags_loose_27_ht200->SetTitle("MC ttjets_sl : Nr. of Btags Distribution for WPloose_27_ht200");
	ttjets_sl_pt_tight_27_or_loose_27_ht200->SetName("ttjets_sl_pt_WPtight_27_or_loose_27_ht200");
	ttjets_sl_pt_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_sl : pT Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_sl_eta_tight_27_or_loose_27_ht200->SetName("ttjets_sl_eta_WPtight_27_or_loose_27_ht200");
	ttjets_sl_eta_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_sl : eta Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_sl_phi_tight_27_or_loose_27_ht200->SetName("ttjets_sl_phi_WPtight_27_or_loose_27_ht200");
	ttjets_sl_phi_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_sl : phi Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_sl_ht_tight_27_or_loose_27_ht200->SetName("ttjets_sl_ht_WPtight_27_or_loose_27_ht200");
	ttjets_sl_ht_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_sl : HT Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_sl_npv_tight_27_or_loose_27_ht200->SetName("ttjets_sl_npv_WPtight_27_or_loose_27_ht200");
	ttjets_sl_npv_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_sl : NPV Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_sl_jet1pt_tight_27_or_loose_27_ht200->SetName("ttjets_sl_jet1pt_WPtight_27_or_loose_27_ht200");
	ttjets_sl_jet1pt_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_sl : Jet1 pT Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_sl_jet1csv_tight_27_or_loose_27_ht200->SetName("ttjets_sl_jet1csv_WPtight_27_or_loose_27_ht200");
	ttjets_sl_jet1csv_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_sl : Jet1 CSV Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_sl_njets_tight_27_or_loose_27_ht200->SetName("ttjets_sl_njets_WPtight_27_or_loose_27_ht200");
	ttjets_sl_njets_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_sl : Nr. of Jets Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_sl_nbtags_tight_27_or_loose_27_ht200->SetName("ttjets_sl_nbtags_WPtight_27_or_loose_27_ht200");
	ttjets_sl_nbtags_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_sl : Nr. of Btags Distribution for WPtight_27_or_loose_27_ht200");
	
	TFile *f4 = new TFile("Distribution_mc_ttjets_dl.root");
	TH1D *ttjets_dl_pt_tight_27 = (TH1D*)f4->Get("Pt_WPTight_27");
	TH1D *ttjets_dl_eta_tight_27 = (TH1D*)f4->Get("Eta_WPTight_27");
	TH1D *ttjets_dl_phi_tight_27 = (TH1D*)f4->Get("Phi_WPTight_27");
	TH1D *ttjets_dl_ht_tight_27 = (TH1D*)f4->Get("HT_WPTight_27");
	TH1D *ttjets_dl_npv_tight_27 = (TH1D*)f4->Get("NumPV_WPTight_27");
	TH1D *ttjets_dl_jet1pt_tight_27 = (TH1D*)f4->Get("Jet1_pt_WPTight_27");
	TH1D *ttjets_dl_jet1csv_tight_27 = (TH1D*)f4->Get("Jet1_csv_WPTight_27");
	TH1D *ttjets_dl_njets_tight_27 = (TH1D*)f4->Get("Njets_WPTight_27");
	TH1D *ttjets_dl_nbtags_tight_27 = (TH1D*)f4->Get("Nbtags_WPTight_27");
	TH1D *ttjets_dl_pt_tight_32 = (TH1D*)f4->Get("Pt_WPTight_32");
	TH1D *ttjets_dl_eta_tight_32 = (TH1D*)f4->Get("Eta_WPTight_32");
	TH1D *ttjets_dl_phi_tight_32 = (TH1D*)f4->Get("Phi_WPTight_32");
	TH1D *ttjets_dl_ht_tight_32 = (TH1D*)f4->Get("HT_WPTight_32");
	TH1D *ttjets_dl_npv_tight_32 = (TH1D*)f4->Get("NumPV_WPTight_32");
	TH1D *ttjets_dl_jet1pt_tight_32 = (TH1D*)f4->Get("Jet1_pt_WPTight_32");
	TH1D *ttjets_dl_jet1csv_tight_32 = (TH1D*)f4->Get("Jet1_csv_WPTight_32");
	TH1D *ttjets_dl_njets_tight_32 = (TH1D*)f4->Get("Njets_WPTight_32");
	TH1D *ttjets_dl_nbtags_tight_32 = (TH1D*)f4->Get("Nbtags_WPTight_32");
	TH1D *ttjets_dl_pt_loose_27_ht200 = (TH1D*)f4->Get("Pt_WPLoose_27_HT200");
	TH1D *ttjets_dl_eta_loose_27_ht200 = (TH1D*)f4->Get("Eta_WPLoose_27_HT200");
	TH1D *ttjets_dl_phi_loose_27_ht200 = (TH1D*)f4->Get("Phi_WPLoose_27_HT200");
	TH1D *ttjets_dl_ht_loose_27_ht200 = (TH1D*)f4->Get("HT_WPLoose_27_HT200");
	TH1D *ttjets_dl_npv_loose_27_ht200 = (TH1D*)f4->Get("NumPV_WPLoose_27_HT200");
	TH1D *ttjets_dl_jet1pt_loose_27_ht200 = (TH1D*)f4->Get("Jet1_pt_WPLoose_27_HT200");
	TH1D *ttjets_dl_jet1csv_loose_27_ht200 = (TH1D*)f4->Get("Jet1_csv_WPLoose_27_HT200");
	TH1D *ttjets_dl_njets_loose_27_ht200 = (TH1D*)f4->Get("Njets_WPLoose_27_HT200");
	TH1D *ttjets_dl_nbtags_loose_27_ht200 = (TH1D*)f4->Get("Nbtags_WPLoose_27_HT200");
	TH1D *ttjets_dl_pt_tight_27_or_loose_27_ht200 = (TH1D*)f4->Get("Pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_dl_eta_tight_27_or_loose_27_ht200 = (TH1D*)f4->Get("Eta_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_dl_phi_tight_27_or_loose_27_ht200 = (TH1D*)f4->Get("Phi_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_dl_ht_tight_27_or_loose_27_ht200 = (TH1D*)f4->Get("HT_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_dl_npv_tight_27_or_loose_27_ht200 = (TH1D*)f4->Get("NumPV_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_dl_jet1pt_tight_27_or_loose_27_ht200 = (TH1D*)f4->Get("Jet1_pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_dl_jet1csv_tight_27_or_loose_27_ht200 = (TH1D*)f4->Get("Jet1_csv_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_dl_njets_tight_27_or_loose_27_ht200 = (TH1D*)f4->Get("Njets_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *ttjets_dl_nbtags_tight_27_or_loose_27_ht200 = (TH1D*)f4->Get("Nbtags_WPTight_27_OR_WPLoose_27_HT200");
	
	ttjets_dl_pt_tight_27->SetName("ttjets_dl_pt_WPTight_27");
	ttjets_dl_pt_tight_27->SetTitle("MC ttjets_dl : pT Distribution for WPTight_27");
	ttjets_dl_eta_tight_27->SetName("ttjets_dl_eta_WPTight_27");
	ttjets_dl_eta_tight_27->SetTitle("MC ttjets_dl : eta Distribution for WPTight_27");
	ttjets_dl_phi_tight_27->SetName("ttjets_dl_phi_WPTight_27");
	ttjets_dl_phi_tight_27->SetTitle("MC ttjets_dl : phi Distribution for WPTight_27");
	ttjets_dl_ht_tight_27->SetName("ttjets_dl_ht_WPTight_27");
	ttjets_dl_ht_tight_27->SetTitle("MC ttjets_dl : HT Distribution for WPTight_27");
	ttjets_dl_npv_tight_27->SetName("ttjets_dl_npv_WPTight_27");
	ttjets_dl_npv_tight_27->SetTitle("MC ttjets_dl : NPV Distribution for WPTight_27");
	ttjets_dl_jet1pt_tight_27->SetName("ttjets_dl_jet1pt_WPTight_27");
	ttjets_dl_jet1pt_tight_27->SetTitle("MC ttjets_dl : Jet1 pT Distribution for WPTight_27");
	ttjets_dl_jet1csv_tight_27->SetName("ttjets_dl_jet1csv_WPTight_27");
	ttjets_dl_jet1csv_tight_27->SetTitle("MC ttjets_dl : Jet1 CSV Distribution for WPTight_27");
	ttjets_dl_njets_tight_27->SetName("ttjets_dl_njets_WPTight_27");
	ttjets_dl_njets_tight_27->SetTitle("MC ttjets_dl : Nr. of Jets Distribution for WPTight_27");
	ttjets_dl_nbtags_tight_27->SetName("ttjets_dl_nbtags_WPTight_27");
	ttjets_dl_nbtags_tight_27->SetTitle("MC ttjets_dl : Nr. of Btags Distribution for WPTight_27");
	ttjets_dl_pt_tight_32->SetName("ttjets_dl_pt_WPTight_32");
	ttjets_dl_pt_tight_32->SetTitle("MC ttjets_dl : pT Distribution for WPTight_32");
	ttjets_dl_eta_tight_32->SetName("ttjets_dl_eta_WPTight_32");
	ttjets_dl_eta_tight_32->SetTitle("MC ttjets_dl : eta Distribution for WPTight_32");
	ttjets_dl_phi_tight_32->SetName("ttjets_dl_phi_WPTight_32");
	ttjets_dl_phi_tight_32->SetTitle("MC ttjets_dl : phi Distribution for WPTight_32");
	ttjets_dl_ht_tight_32->SetName("ttjets_dl_ht_WPTight_32");
	ttjets_dl_ht_tight_32->SetTitle("MC ttjets_dl : HT Distribution for WPTight_32");
	ttjets_dl_npv_tight_32->SetName("ttjets_dl_npv_WPTight_32");
	ttjets_dl_npv_tight_32->SetTitle("MC ttjets_dl : NPV Distribution for WPTight_32");
	ttjets_dl_jet1pt_tight_32->SetName("ttjets_dl_jet1pt_WPTight_32");
	ttjets_dl_jet1pt_tight_32->SetTitle("MC ttjets_dl : Jet1 pT Distribution for WPTight_32");
	ttjets_dl_jet1csv_tight_32->SetName("ttjets_dl_jet1csv_WPTight_32");
	ttjets_dl_jet1csv_tight_32->SetTitle("MC ttjets_dl : Jet1 CSV Distribution for WPTight_32");
	ttjets_dl_njets_tight_32->SetName("ttjets_dl_njets_WPTight_32");
	ttjets_dl_njets_tight_32->SetTitle("MC ttjets_dl : Nr. of Jets Distribution for WPTight_32");
	ttjets_dl_nbtags_tight_32->SetName("ttjets_dl_nbtags_WPTight_32");
	ttjets_dl_nbtags_tight_32->SetTitle("MC ttjets_dl : Nr. of Btags Distribution for WPTight_32");
	ttjets_dl_pt_loose_27_ht200->SetName("ttjets_dl_pt_WPloose_27_ht200");
	ttjets_dl_pt_loose_27_ht200->SetTitle("MC ttjets_dl : pT Distribution for WPloose_27_ht200");
	ttjets_dl_eta_loose_27_ht200->SetName("ttjets_dl_eta_WPloose_27_ht200");
	ttjets_dl_eta_loose_27_ht200->SetTitle("MC ttjets_dl : eta Distribution for WPloose_27_ht200");
	ttjets_dl_phi_loose_27_ht200->SetName("ttjets_dl_phi_WPloose_27_ht200");
	ttjets_dl_phi_loose_27_ht200->SetTitle("MC ttjets_dl : phi Distribution for WPloose_27_ht200");
	ttjets_dl_ht_loose_27_ht200->SetName("ttjets_dl_ht_WPloose_27_ht200");
	ttjets_dl_ht_loose_27_ht200->SetTitle("MC ttjets_dl : HT Distribution for WPloose_27_ht200");
	ttjets_dl_npv_loose_27_ht200->SetName("ttjets_dl_npv_WPloose_27_ht200");
	ttjets_dl_npv_loose_27_ht200->SetTitle("MC ttjets_dl : NPV Distribution for WPloose_27_ht200");
	ttjets_dl_jet1pt_loose_27_ht200->SetName("ttjets_dl_jet1pt_WPloose_27_ht200");
	ttjets_dl_jet1pt_loose_27_ht200->SetTitle("MC ttjets_dl : Jet1 pT Distribution for WPloose_27_ht200");
	ttjets_dl_jet1csv_loose_27_ht200->SetName("ttjets_dl_jet1csv_WPloose_27_ht200");
	ttjets_dl_jet1csv_loose_27_ht200->SetTitle("MC ttjets_dl : Jet1 CSV Distribution for WPloose_27_ht200");
	ttjets_dl_njets_loose_27_ht200->SetName("ttjets_dl_njets_WPloose_27_ht200");
	ttjets_dl_njets_loose_27_ht200->SetTitle("MC ttjets_dl : Nr. of Jets Distribution for WPloose_27_ht200");
	ttjets_dl_nbtags_loose_27_ht200->SetName("ttjets_dl_nbtags_WPloose_27_ht200");
	ttjets_dl_nbtags_loose_27_ht200->SetTitle("MC ttjets_dl : Nr. of Btags Distribution for WPloose_27_ht200");
	ttjets_dl_pt_tight_27_or_loose_27_ht200->SetName("ttjets_dl_pt_WPtight_27_or_loose_27_ht200");
	ttjets_dl_pt_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_dl : pT Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_dl_eta_tight_27_or_loose_27_ht200->SetName("ttjets_dl_eta_WPtight_27_or_loose_27_ht200");
	ttjets_dl_eta_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_dl : eta Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_dl_phi_tight_27_or_loose_27_ht200->SetName("ttjets_dl_phi_WPtight_27_or_loose_27_ht200");
	ttjets_dl_phi_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_dl : phi Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_dl_ht_tight_27_or_loose_27_ht200->SetName("ttjets_dl_ht_WPtight_27_or_loose_27_ht200");
	ttjets_dl_ht_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_dl : HT Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_dl_npv_tight_27_or_loose_27_ht200->SetName("ttjets_dl_npv_WPtight_27_or_loose_27_ht200");
	ttjets_dl_npv_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_dl : NPV Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_dl_jet1pt_tight_27_or_loose_27_ht200->SetName("ttjets_dl_jet1pt_WPtight_27_or_loose_27_ht200");
	ttjets_dl_jet1pt_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_dl : Jet1 pT Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_dl_jet1csv_tight_27_or_loose_27_ht200->SetName("ttjets_dl_jet1csv_WPtight_27_or_loose_27_ht200");
	ttjets_dl_jet1csv_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_dl : Jet1 CSV Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_dl_njets_tight_27_or_loose_27_ht200->SetName("ttjets_dl_njets_WPtight_27_or_loose_27_ht200");
	ttjets_dl_njets_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_dl : Nr. of Jets Distribution for WPtight_27_or_loose_27_ht200");
	ttjets_dl_nbtags_tight_27_or_loose_27_ht200->SetName("ttjets_dl_nbtags_WPtight_27_or_loose_27_ht200");
	ttjets_dl_nbtags_tight_27_or_loose_27_ht200->SetTitle("MC ttjets_dl : Nr. of Btags Distribution for WPtight_27_or_loose_27_ht200");
	
	
	TFile *f5 = new TFile("Distribution_data_2016.root");
	TH1D *data_pt_tight_27 = (TH1D*)f5->Get("Pt_WPTight_27");
	TH1D *data_eta_tight_27 = (TH1D*)f5->Get("Eta_WPTight_27");
	TH1D *data_phi_tight_27 = (TH1D*)f5->Get("Phi_WPTight_27");
	TH1D *data_ht_tight_27 = (TH1D*)f5->Get("HT_WPTight_27");
	TH1D *data_npv_tight_27 = (TH1D*)f5->Get("NumPV_WPTight_27");
	TH1D *data_jet1pt_tight_27 = (TH1D*)f5->Get("Jet1_pt_WPTight_27");
	TH1D *data_jet1csv_tight_27 = (TH1D*)f5->Get("Jet1_csv_WPTight_27");
	TH1D *data_njets_tight_27 = (TH1D*)f5->Get("Njets_WPTight_27");
	TH1D *data_nbtags_tight_27 = (TH1D*)f5->Get("Nbtags_WPTight_27");
	TH1D *data_pt_tight_32 = (TH1D*)f5->Get("Pt_WPTight_32");
	TH1D *data_eta_tight_32 = (TH1D*)f5->Get("Eta_WPTight_32");
	TH1D *data_phi_tight_32 = (TH1D*)f5->Get("Phi_WPTight_32");
	TH1D *data_ht_tight_32 = (TH1D*)f5->Get("HT_WPTight_32");
	TH1D *data_npv_tight_32 = (TH1D*)f5->Get("NumPV_WPTight_32");
	TH1D *data_jet1pt_tight_32 = (TH1D*)f5->Get("Jet1_pt_WPTight_32");
	TH1D *data_jet1csv_tight_32 = (TH1D*)f5->Get("Jet1_csv_WPTight_32");
	TH1D *data_njets_tight_32 = (TH1D*)f5->Get("Njets_WPTight_32");
	TH1D *data_nbtags_tight_32 = (TH1D*)f5->Get("Nbtags_WPTight_32");
	TH1D *data_pt_loose_27_ht200 = (TH1D*)f5->Get("Pt_WPLoose_27_HT200");
	TH1D *data_eta_loose_27_ht200 = (TH1D*)f5->Get("Eta_WPLoose_27_HT200");
	TH1D *data_phi_loose_27_ht200 = (TH1D*)f5->Get("Phi_WPLoose_27_HT200");
	TH1D *data_ht_loose_27_ht200 = (TH1D*)f5->Get("HT_WPLoose_27_HT200");
	TH1D *data_npv_loose_27_ht200 = (TH1D*)f5->Get("NumPV_WPLoose_27_HT200");
	TH1D *data_jet1pt_loose_27_ht200 = (TH1D*)f5->Get("Jet1_pt_WPLoose_27_HT200");
	TH1D *data_jet1csv_loose_27_ht200 = (TH1D*)f5->Get("Jet1_csv_WPLoose_27_HT200");
	TH1D *data_njets_loose_27_ht200 = (TH1D*)f5->Get("Njets_WPLoose_27_HT200");
	TH1D *data_nbtags_loose_27_ht200 = (TH1D*)f5->Get("Nbtags_WPLoose_27_HT200");
	TH1D *data_pt_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_eta_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Eta_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_phi_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Phi_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_ht_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("HT_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_npv_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("NumPV_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_jet1pt_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Jet1_pt_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_jet1csv_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Jet1_csv_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_njets_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Njets_WPTight_27_OR_WPLoose_27_HT200");
	TH1D *data_nbtags_tight_27_or_loose_27_ht200 = (TH1D*)f5->Get("Nbtags_WPTight_27_OR_WPLoose_27_HT200");
	
	data_pt_tight_27->SetName("data_pt_WPTight_27");
	data_pt_tight_27->SetTitle("DATA : pT Distribution for WPTight_27");
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
	data_pt_tight_32->SetTitle("DATA : pT Distribution for WPTight_32");
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
	data_pt_loose_27_ht200->SetTitle("DATA : pT Distribution for WPloose_27_ht200");
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
	data_pt_tight_27_or_loose_27_ht200->SetTitle("DATA : pT Distribution for WPtight_27_or_loose_27_ht200");
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
	// Scaling and Normlization of mc signal : tthbb
	///////////////////////////////////////////////////////////
	
	double factor_tthbb = 600.0;
	double sigma_tthbb = 0.5824*0.5071; // pb
	double N_total_tthbb = 3993304;
	double norm_tthbb = (L*sigma_tthbb*1000)/N_total_tthbb;
	double scale_tthbb = factor_tthbb*norm_tthbb;
	
	tthbb_pt_tight_27->Scale(scale_tthbb);
	tthbb_eta_tight_27->Scale(scale_tthbb);
	tthbb_phi_tight_27->Scale(scale_tthbb);
	tthbb_ht_tight_27->Scale(scale_tthbb);
	tthbb_npv_tight_27->Scale(scale_tthbb);
	tthbb_jet1pt_tight_27->Scale(scale_tthbb);
	tthbb_jet1csv_tight_27->Scale(scale_tthbb);
	tthbb_njets_tight_27->Scale(scale_tthbb);
	tthbb_nbtags_tight_27->Scale(scale_tthbb);
	tthbb_pt_tight_32->Scale(scale_tthbb);
	tthbb_eta_tight_32->Scale(scale_tthbb);
	tthbb_phi_tight_32->Scale(scale_tthbb);
	tthbb_ht_tight_32->Scale(scale_tthbb);
	tthbb_npv_tight_32->Scale(scale_tthbb);
	tthbb_jet1pt_tight_32->Scale(scale_tthbb);
	tthbb_jet1csv_tight_32->Scale(scale_tthbb);
	tthbb_njets_tight_32->Scale(scale_tthbb);
	tthbb_nbtags_tight_32->Scale(scale_tthbb);
	tthbb_pt_loose_27_ht200->Scale(scale_tthbb);
	tthbb_eta_loose_27_ht200->Scale(scale_tthbb);
	tthbb_phi_loose_27_ht200->Scale(scale_tthbb);
	tthbb_ht_loose_27_ht200->Scale(scale_tthbb);
	tthbb_npv_loose_27_ht200->Scale(scale_tthbb);
	tthbb_jet1pt_loose_27_ht200->Scale(scale_tthbb);
	tthbb_jet1csv_loose_27_ht200->Scale(scale_tthbb);
	tthbb_njets_loose_27_ht200->Scale(scale_tthbb);
	tthbb_nbtags_loose_27_ht200->Scale(scale_tthbb);
	tthbb_pt_tight_27_or_loose_27_ht200->Scale(scale_tthbb);
	tthbb_eta_tight_27_or_loose_27_ht200->Scale(scale_tthbb);
	tthbb_phi_tight_27_or_loose_27_ht200->Scale(scale_tthbb);
	tthbb_ht_tight_27_or_loose_27_ht200->Scale(scale_tthbb);
	tthbb_npv_tight_27_or_loose_27_ht200->Scale(scale_tthbb);
	tthbb_jet1pt_tight_27_or_loose_27_ht200->Scale(scale_tthbb);
	tthbb_jet1csv_tight_27_or_loose_27_ht200->Scale(scale_tthbb);
	tthbb_njets_tight_27_or_loose_27_ht200->Scale(scale_tthbb);
	tthbb_nbtags_tight_27_or_loose_27_ht200->Scale(scale_tthbb);
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	//Adding and Normalization of mc background : ttjets_SL + ttjets_DL + tthnonbb
	///////////////////////////////////////////////////////////////////////////////////////////////
	
	double factor_tthnonbb = 0;
	double sigma_tthnonbb = 0.4176*0.5071; // pb
	double N_total_tthnonbb = 3975854;
	double norm_tthnonbb = (L*sigma_tthnonbb*1000)/N_total_tthnonbb;
	double scale_tthnonbb = factor_tthnonbb*norm_tthnonbb;
	
	double factor_ttjets_sl = 1;
	double sigma_ttjets_sl = 364.3; // pb
	double N_total_ttjets_sl = 157387260;
	double norm_ttjets_sl = (L*sigma_ttjets_sl*1000)/N_total_ttjets_sl;
	double scale_ttjets_sl = factor_ttjets_sl*norm_ttjets_sl;
	
	double factor_ttjets_dl = 0;
	double sigma_ttjets_dl = 87.3; // pb
	double N_total_ttjets_dl = 75383000;
	double norm_ttjets_dl = (L*sigma_ttjets_dl*1000)/N_total_ttjets_dl;
	double scale_ttjets_dl = factor_ttjets_dl*norm_ttjets_dl;
	
	TH1D* bkgd_pt_tight_27 = new TH1D("Pt_WPTight_27","Pt Distribution for WPTight_27;pT (GeV);Nr. of Events",150,0,300);
    TH1D* bkgd_eta_tight_27 = new TH1D("Eta_WPTight_27","Eta Distribution for WPTight_27;#eta;Nr. of Events",30,-3,3);
    TH1D* bkgd_phi_tight_27 = new TH1D("Phi_WPTight_27","Phi Distribution for WPTight_27;#phi;Nr. of Events",30,-3,3);
    TH1D* bkgd_ht_tight_27 = new TH1D("HT_WPTight_27","HT Distribution for WPTight_27;HT (GeV);Nr. of Events",200,0,1000);	
    TH1D* bkgd_npv_tight_27 = new TH1D("NumPV_WPTight_27","NumPV Distribution for WPTight_27;numPV;Nr. of Events",50,0,50);		
    TH1D* bkgd_jet1pt_tight_27 = new TH1D("Jet1_pt_WPTight_27","Jet1 pT Distribution for WPTight_27;pT (GeV);Nr. of Events",250,0,500);
    TH1D* bkgd_jet1csv_tight_27 = new TH1D("Jet1_csv_WPTight_27","Jet1 csv Distribution for WPTight_27;csv;Nr. of Events",50,0,1);
    TH1D* bkgd_njets_tight_27 = new TH1D("Njets_WPTight_27","N_jets Distribution for WPTight_27;Nr. of jets;Nr. of Events",9,4,13);
    TH1D* bkgd_nbtags_tight_27 = new TH1D("Nbtags_WPTight_27","N_btags Distribution for WPTight_27;Nr. of btags;Nr. of Events",7,2,9);	
    TH1D* bkgd_pt_tight_32 = new TH1D("Pt_WPTight_32","Pt Distribution for WPTight_32;pT (GeV);Nr. of Events",150,0,300);
    TH1D* bkgd_eta_tight_32 = new TH1D("Eta_WPTight_32","Eta Distribution for WPTight_32;#eta;Nr. of Events",30,-3,3);
    TH1D* bkgd_phi_tight_32 = new TH1D("Phi_WPTight_32","Phi Distribution for WPTight_32;#phi;Nr. of Events",30,-3,3);
    TH1D* bkgd_ht_tight_32 = new TH1D("HT_WPTight_32","HT Distribution for WPTight_32;HT (GeV);Nr. of Events",200,0,1000);	
    TH1D* bkgd_npv_tight_32 = new TH1D("NumPV_WPTight_32","NumPV Distribution for WPTight_32;numPV;Nr. of Events",50,0,50);
    TH1D* bkgd_jet1pt_tight_32 = new TH1D("Jet1_pt_WPTight_32","Jet1 pT Distribution for WPTight_32;pT (GeV);Nr. of Events",250,0,500);
    TH1D* bkgd_jet1csv_tight_32 = new TH1D("Jet1_csv_WPTight_32","Jet1 csv Distribution for WPTight_32;csv;Nr. of Events",50,0,1);	
    TH1D* bkgd_njets_tight_32 = new TH1D("Njets_WPTight_32","N_jets Distribution for WPTight_32;Nr. of jets;Nr. of Events",9,4,13);
    TH1D* bkgd_nbtags_tight_32 = new TH1D("Nbtags_WPTight_32","N_btags Distribution for WPTight_32;Nr. of btags;Nr. of Events",7,2,9);	
    TH1D* bkgd_pt_loose_27_ht200 = new TH1D("Pt_WPLoose_27_HT200","Pt Distribution for WPLoose_27_HT200;pT (GeV);Nr. of Events",150,0,300);
    TH1D* bkgd_eta_loose_27_ht200 = new TH1D("Eta_WPLoose_27_HT200","Eta Distribution for WPLoose_27_HT200;#eta;Nr. of Events",30,-3,3);
    TH1D* bkgd_phi_loose_27_ht200 = new TH1D("Phi_WPLoose_27_HT200","Phi Distribution for WPLoose_27_HT200;#phi;Nr. of Events",30,-3,3);
    TH1D* bkgd_ht_loose_27_ht200 = new TH1D("HT_WPLoose_27_HT200","HT Distribution for WPLoose_27_HT200;HT (GeV);Nr. of Events",200,0,1000);	
    TH1D* bkgd_npv_loose_27_ht200 = new TH1D("NumPV_WPLoose_27_HT200","NumPV Distribution for WPLoose_27_HT200;numPV;Nr. of Events",50,0,50);
    TH1D* bkgd_jet1pt_loose_27_ht200 = new TH1D("Jet1_pt_WPLoose_27_HT200","Jet1 pT Distribution for WPLoose_27_HT200;pT (GeV);Nr. of Events",250,0,500);
    TH1D* bkgd_jet1csv_loose_27_ht200 = new TH1D("Jet1_csv_WPLoose_27_HT200","Jet1 csv Distribution for WPLoose_27_HT200;csv;Nr. of Events",50,0,1);		
    TH1D* bkgd_njets_loose_27_ht200 = new TH1D("Njets_WPLoose_27_HT200","N_jets Distribution for WPLoose_27_HT200;Nr. of jets;Nr. of Events",9,4,13);
    TH1D* bkgd_nbtags_loose_27_ht200 = new TH1D("Nbtags_WPLoose_27_HT200","N_btags Distribution for WPLoose_27_HT200;Nr. of btags;Nr. of Events",7,2,9);	
    TH1D* bkgd_pt_tight_27_or_loose_27_ht200 = new TH1D("Pt_WPTight_27_OR_WPLoose_27_HT200","Pt Distribution for WPTight_27_OR_WPLoose_27_HT200;pT (GeV);Nr. of Events",150,0,300);
    TH1D* bkgd_eta_tight_27_or_loose_27_ht200 = new TH1D("Eta_WPTight_27_OR_WPLoose_27_HT200","Eta Distribution for WPTight_27_OR_WPLoose_27_HT200;#eta;Nr. of Events",30,-3,3);
    TH1D* bkgd_phi_tight_27_or_loose_27_ht200 = new TH1D("Phi_WPTight_27_OR_WPLoose_27_HT200","Phi Distribution for WPTight_27_OR_WPLoose_27_HT200;#phi;Nr. of Events",30,-3,3);
    TH1D* bkgd_ht_tight_27_or_loose_27_ht200 = new TH1D("HT_WPTight_27_OR_WPLoose_27_HT200","HT Distribution for WPTight_27_OR_WPLoose_27_HT200;HT (GeV);Nr. of Events",200,0,1000);	
    TH1D* bkgd_npv_tight_27_or_loose_27_ht200 = new TH1D("NumPV_WPTight_27_OR_WPLoose_27_HT200","NumPV Distribution for WPTight_27_OR_WPLoose_27_HT200;numPV;Nr. of Events",50,0,50);
    TH1D* bkgd_jet1pt_tight_27_or_loose_27_ht200 = new TH1D("Jet1_pt_WPTight_27_OR_WPLoose_27_HT200","Jet1 pT Distribution for WPTight_27_OR_WPLoose_27_HT200;pT (GeV);Nr. of Events",250,0,500);
    TH1D* bkgd_jet1csv_tight_27_or_loose_27_ht200 = new TH1D("Jet1_csv_WPTight_27_OR_WPLoose_27_HT200","Jet1 csv Distribution for WPTight_27_OR_WPLoose_27_HT200;csv;Nr. of Events",50,0,1);			
    TH1D* bkgd_njets_tight_27_or_loose_27_ht200 = new TH1D("Njets_WPTight_27_OR_WPLoose_27_HT200","N_jets Distribution for WPTight_27_OR_WPLoose_27_HT200;Nr. of jets;Nr. of Events",9,4,13);
    TH1D* bkgd_nbtags_tight_27_or_loose_27_ht200 = new TH1D("Nbtags_WPTight_27_OR_WPLoose_27_HT200","N_btags Distribution for WPTight_27_OR_WPLoose_27_HT200;Nr. of btags;Nr. of Events",7,2,9);
	
	bkgd_pt_tight_27->Add(ttjets_sl_pt_tight_27, ttjets_dl_pt_tight_27, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_eta_tight_27->Add(ttjets_sl_eta_tight_27, ttjets_dl_eta_tight_27, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_phi_tight_27->Add(ttjets_sl_phi_tight_27, ttjets_dl_phi_tight_27, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_ht_tight_27->Add(ttjets_sl_ht_tight_27, ttjets_dl_ht_tight_27, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_npv_tight_27->Add(ttjets_sl_npv_tight_27, ttjets_dl_npv_tight_27, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_jet1pt_tight_27->Add(ttjets_sl_jet1pt_tight_27, ttjets_dl_jet1pt_tight_27, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_jet1csv_tight_27->Add(ttjets_sl_jet1csv_tight_27, ttjets_dl_jet1csv_tight_27, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_njets_tight_27->Add(ttjets_sl_njets_tight_27, ttjets_dl_njets_tight_27, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_nbtags_tight_27->Add(ttjets_sl_nbtags_tight_27, ttjets_dl_nbtags_tight_27, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_pt_tight_32->Add(ttjets_sl_pt_tight_32, ttjets_dl_pt_tight_32, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_eta_tight_32->Add(ttjets_sl_eta_tight_32, ttjets_dl_eta_tight_32, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_phi_tight_32->Add(ttjets_sl_phi_tight_32, ttjets_dl_phi_tight_32, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_ht_tight_32->Add(ttjets_sl_ht_tight_32, ttjets_dl_ht_tight_32, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_npv_tight_32->Add(ttjets_sl_npv_tight_32, ttjets_dl_npv_tight_32, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_jet1pt_tight_32->Add(ttjets_sl_jet1pt_tight_32, ttjets_dl_jet1pt_tight_32, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_jet1csv_tight_32->Add(ttjets_sl_jet1csv_tight_32, ttjets_dl_jet1csv_tight_32, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_njets_tight_32->Add(ttjets_sl_njets_tight_32, ttjets_dl_njets_tight_32, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_nbtags_tight_32->Add(ttjets_sl_nbtags_tight_32, ttjets_dl_nbtags_tight_32, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_pt_loose_27_ht200->Add(ttjets_sl_pt_loose_27_ht200, ttjets_dl_pt_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_eta_loose_27_ht200->Add(ttjets_sl_eta_loose_27_ht200, ttjets_dl_eta_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_phi_loose_27_ht200->Add(ttjets_sl_phi_loose_27_ht200, ttjets_dl_phi_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_ht_loose_27_ht200->Add(ttjets_sl_ht_loose_27_ht200, ttjets_dl_ht_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_npv_loose_27_ht200->Add(ttjets_sl_npv_loose_27_ht200, ttjets_dl_npv_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_jet1pt_loose_27_ht200->Add(ttjets_sl_jet1pt_loose_27_ht200, ttjets_dl_jet1pt_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_jet1csv_loose_27_ht200->Add(ttjets_sl_jet1csv_loose_27_ht200, ttjets_dl_jet1csv_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_njets_loose_27_ht200->Add(ttjets_sl_njets_loose_27_ht200, ttjets_dl_njets_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_nbtags_loose_27_ht200->Add(ttjets_sl_nbtags_loose_27_ht200, ttjets_dl_nbtags_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_pt_tight_27_or_loose_27_ht200->Add(ttjets_sl_pt_tight_27_or_loose_27_ht200, ttjets_dl_pt_tight_27_or_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_eta_tight_27_or_loose_27_ht200->Add(ttjets_sl_eta_tight_27_or_loose_27_ht200, ttjets_dl_eta_tight_27_or_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_phi_tight_27_or_loose_27_ht200->Add(ttjets_sl_phi_tight_27_or_loose_27_ht200, ttjets_dl_phi_tight_27_or_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_ht_tight_27_or_loose_27_ht200->Add(ttjets_sl_ht_tight_27_or_loose_27_ht200, ttjets_dl_ht_tight_27_or_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_npv_tight_27_or_loose_27_ht200->Add(ttjets_sl_npv_tight_27_or_loose_27_ht200, ttjets_dl_npv_tight_27_or_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_jet1pt_tight_27_or_loose_27_ht200->Add(ttjets_sl_jet1pt_tight_27_or_loose_27_ht200, ttjets_dl_jet1pt_tight_27_or_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_jet1csv_tight_27_or_loose_27_ht200->Add(ttjets_sl_jet1csv_tight_27_or_loose_27_ht200, ttjets_dl_jet1csv_tight_27_or_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_njets_tight_27_or_loose_27_ht200->Add(ttjets_sl_njets_tight_27_or_loose_27_ht200, ttjets_dl_njets_tight_27_or_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	bkgd_nbtags_tight_27_or_loose_27_ht200->Add(ttjets_sl_nbtags_tight_27_or_loose_27_ht200, ttjets_dl_nbtags_tight_27_or_loose_27_ht200, scale_ttjets_sl, scale_ttjets_dl);
	
	bkgd_pt_tight_27->Add(tthnonbb_pt_tight_27, scale_tthnonbb);
	bkgd_eta_tight_27->Add(tthnonbb_eta_tight_27, scale_tthnonbb );
	bkgd_phi_tight_27->Add(tthnonbb_phi_tight_27, scale_tthnonbb);
	bkgd_ht_tight_27->Add(tthnonbb_ht_tight_27, scale_tthnonbb);
	bkgd_npv_tight_27->Add(tthnonbb_npv_tight_27, scale_tthnonbb );
	bkgd_jet1pt_tight_27->Add(tthnonbb_jet1pt_tight_27, scale_tthnonbb );
	bkgd_jet1csv_tight_27->Add(tthnonbb_jet1csv_tight_27, scale_tthnonbb );
	bkgd_njets_tight_27->Add(tthnonbb_njets_tight_27, scale_tthnonbb);
	bkgd_nbtags_tight_27->Add(tthnonbb_nbtags_tight_27, scale_tthnonbb);
	bkgd_pt_tight_32->Add(tthnonbb_pt_tight_32, scale_tthnonbb);
	bkgd_eta_tight_32->Add(tthnonbb_eta_tight_32, scale_tthnonbb);
	bkgd_phi_tight_32->Add(tthnonbb_phi_tight_32, scale_tthnonbb);
	bkgd_ht_tight_32->Add(tthnonbb_ht_tight_32, scale_tthnonbb);
	bkgd_npv_tight_32->Add(tthnonbb_npv_tight_32, scale_tthnonbb);
	bkgd_jet1pt_tight_32->Add(tthnonbb_jet1pt_tight_32, scale_tthnonbb);
	bkgd_jet1csv_tight_32->Add(tthnonbb_jet1csv_tight_32, scale_tthnonbb);
	bkgd_njets_tight_32->Add(tthnonbb_njets_tight_32, scale_tthnonbb);
	bkgd_nbtags_tight_32->Add(tthnonbb_nbtags_tight_32, scale_tthnonbb);
	bkgd_pt_loose_27_ht200->Add(tthnonbb_pt_loose_27_ht200, scale_tthnonbb);
	bkgd_eta_loose_27_ht200->Add(tthnonbb_eta_loose_27_ht200, scale_tthnonbb);
	bkgd_phi_loose_27_ht200->Add(tthnonbb_phi_loose_27_ht200, scale_tthnonbb);
	bkgd_ht_loose_27_ht200->Add(tthnonbb_ht_loose_27_ht200, scale_tthnonbb);
	bkgd_npv_loose_27_ht200->Add(tthnonbb_npv_loose_27_ht200, scale_tthnonbb);
	bkgd_jet1pt_loose_27_ht200->Add(tthnonbb_jet1pt_loose_27_ht200, scale_tthnonbb);
	bkgd_jet1csv_loose_27_ht200->Add(tthnonbb_jet1csv_loose_27_ht200, scale_tthnonbb);
	bkgd_njets_loose_27_ht200->Add(tthnonbb_njets_loose_27_ht200, scale_tthnonbb);
	bkgd_nbtags_loose_27_ht200->Add(tthnonbb_nbtags_loose_27_ht200, scale_tthnonbb);
	bkgd_pt_tight_27_or_loose_27_ht200->Add(tthnonbb_pt_tight_27_or_loose_27_ht200, scale_tthnonbb);
	bkgd_eta_tight_27_or_loose_27_ht200->Add(tthnonbb_eta_tight_27_or_loose_27_ht200, scale_tthnonbb);
	bkgd_phi_tight_27_or_loose_27_ht200->Add(tthnonbb_phi_tight_27_or_loose_27_ht200, scale_tthnonbb);
	bkgd_ht_tight_27_or_loose_27_ht200->Add(tthnonbb_ht_tight_27_or_loose_27_ht200, scale_tthnonbb);
	bkgd_npv_tight_27_or_loose_27_ht200->Add(tthnonbb_npv_tight_27_or_loose_27_ht200, scale_tthnonbb);
	bkgd_jet1pt_tight_27_or_loose_27_ht200->Add(tthnonbb_jet1pt_tight_27_or_loose_27_ht200, scale_tthnonbb);
	bkgd_jet1csv_tight_27_or_loose_27_ht200->Add(tthnonbb_jet1csv_tight_27_or_loose_27_ht200, scale_tthnonbb);
	bkgd_njets_tight_27_or_loose_27_ht200->Add(tthnonbb_njets_tight_27_or_loose_27_ht200, scale_tthnonbb);
	bkgd_nbtags_tight_27_or_loose_27_ht200->Add(tthnonbb_nbtags_tight_27_or_loose_27_ht200, scale_tthnonbb);
    
	
	////////////////////////////////////////////////////////
	//Plotting
	////////////////////////////////////////////////////////
	
	TCanvas *c1 = new TCanvas("c1","test",1100,650);
	c1->DrawFrame(0,0,300,11000,"Electron pT (for WPTight_27) ; pT (Electron) [GeV] ; Nr. of Events");
	TLegend* leg1 = new TLegend(0.65,0.70,0.85,0.85);
    leg1->SetFillColor(kWhite);
    leg1->SetFillStyle(1001);
    leg1->AddEntry(tthbb_pt_tight_27,"MC : ttH x 600","L");
    leg1->AddEntry(bkgd_pt_tight_27,"MC : Background","L");
    leg1->AddEntry(data_pt_tight_27,"Data","L");	
	//tthbb_pt_tight_27->SetLineColor(kRed);
	//tthbb_pt_tight_27->SetLineWidth(2);
	tthbb_pt_tight_27->SetFillStyle(3004);
	tthbb_pt_tight_27->SetFillColor(kOrange); 
	tthbb_pt_tight_27->Draw("same");
	//bkgd_pt_tight_27->SetLineColor(kGreen+3);
	//bkgd_pt_tight_27->SetLineWidth(2);
	bkgd_pt_tight_27->SetFillStyle(3004);
	bkgd_pt_tight_27->SetFillColor(kOrange); 
	//bkgd_pt_tight_27->Draw("same");
	data_pt_tight_27->SetLineColor(kBlue);
	data_pt_tight_27->SetLineWidth(2);
	//data_pt_tight_27->Draw("same");
	leg1->Draw("same");
	c1->Print("pt_tight_27_data_mc.png");
	delete c1;
	delete leg1;
	
	/*
	TCanvas *c2= new TCanvas("c2","test",1100,650);
	c2->DrawFrame(0,0,300,11000,"Electron pT (for WPTight_32) ; pT (Electron) [GeV] ; Nr. of Events");
	TLegend* leg2 = new TLegend(0.65,0.70,0.85,0.85);
    leg2->SetFillColor(kWhite);
    leg2->SetFillStyle(1001);
    leg2->AddEntry(tthbb_pt_tight_32,"MC : ttH x 600","L");
    leg2->AddEntry(bkgd_pt_tight_32,"MC : Background","L");
    leg2->AddEntry(data_pt_tight_32,"Data","L");	
	tthbb_pt_tight_32->SetLineColor(kRed);
	tthbb_pt_tight_32->SetLineWidth(2);
	tthbb_pt_tight_32->Draw("same");
	//bkgd_pt_tight_32->SetLineColor(kGreen+3);
	//bkgd_pt_tight_32->SetLineWidth(2);
	bkgd_pt_tight_32->SetFillStyle(1001);
	bkgd_pt_tight_32->SetFillColor(kOrange); 
	bkgd_pt_tight_32->Draw("same");
	data_pt_tight_32->SetLineColor(kBlue);
	data_pt_tight_32->SetLineWidth(2);
	data_pt_tight_32->Draw("same");
	leg2->Draw("same");
	c2->Print("pt_tight_32_data_mc.png");
	delete c2;
	delete leg2;
	
	TCanvas *c3 = new TCanvas("c3","test",1100,650);
	c3->DrawFrame(0,0,300,11000,"Electron pT (for WPLoose_27_HT200) ; pT (Electron) [GeV] ; Nr. of Events");
	TLegend* leg3 = new TLegend(0.65,0.70,0.85,0.85);
    leg3->SetFillColor(kWhite);
    leg3->SetFillStyle(1001);
    leg3->AddEntry(tthbb_pt_loose_27_ht200,"MC : ttH x 600","L");
    leg3->AddEntry(bkgd_pt_loose_27_ht200,"MC : Background","L");
    leg3->AddEntry(data_pt_loose_27_ht200,"Data","L");	
	tthbb_pt_loose_27_ht200->SetLineColor(kRed);
	tthbb_pt_loose_27_ht200->SetLineWidth(2);
	tthbb_pt_loose_27_ht200->Draw("same");
	//bkgd_pt_loose_27_ht200->SetLineColor(kGreen+3);
	//bkgd_pt_loose_27_ht200->SetLineWidth(2);
	bkgd_pt_loose_27_ht200->SetFillStyle(1001);
	bkgd_pt_loose_27_ht200->SetFillColor(kOrange); 
	bkgd_pt_loose_27_ht200->Draw("same");
	data_pt_loose_27_ht200->SetLineColor(kBlue);
	data_pt_loose_27_ht200->SetLineWidth(2);
	data_pt_loose_27_ht200->Draw("same");
	leg3->Draw("same");
	c3->Print("pt_loose_27_ht200_data_mc.png");
	delete c3;
	delete leg3;
	
	TCanvas *c4= new TCanvas("c4","test",1100,650);
	c4->DrawFrame(0,0,300,11000,"Electron pT (for WPTight_27_OR_WPLoose_27_HT200) ; pT (Electron) [GeV] ; Nr. of Events");
	TLegend* leg4 = new TLegend(0.65,0.70,0.85,0.85);
    leg4->SetFillColor(kWhite);
    leg4->SetFillStyle(1001);
    leg4->AddEntry(tthbb_pt_tight_27_or_loose_27_ht200,"MC : ttH x 600","L");
    leg4->AddEntry(bkgd_pt_tight_27_or_loose_27_ht200,"MC : Background","L");
    leg4->AddEntry(data_pt_tight_27_or_loose_27_ht200,"Data","L");	
	tthbb_pt_tight_27_or_loose_27_ht200->SetLineColor(kRed);
	tthbb_pt_tight_27_or_loose_27_ht200->SetLineWidth(2);
	tthbb_pt_tight_27_or_loose_27_ht200->Draw("same");
	//bkgd_pt_tight_27_or_loose_27_ht200->SetLineColor(kGreen+3);
	//bkgd_pt_tight_27_or_loose_27_ht200->SetLineWidth(2);
	bkgd_pt_tight_27_or_loose_27_ht200->SetFillStyle(1001);
	bkgd_pt_tight_27_or_loose_27_ht200->SetFillColor(kOrange); 
	bkgd_pt_tight_27_or_loose_27_ht200->Draw("same");
	data_pt_tight_27_or_loose_27_ht200->SetLineColor(kBlue);
	data_pt_tight_27_or_loose_27_ht200->SetLineWidth(2);
	data_pt_tight_27_or_loose_27_ht200->Draw("same");
	leg4->Draw("same");
	c4->Print("pt_tight_27_or_loose_27_ht200_data_mc.png");
	delete c4;
	delete leg4;
	*/
	
	
	
	
	
	
	
	
	
	
	

	delete f1;
	delete f2;
	delete f3;
	delete f4;
	delete f5;
	return;
}


