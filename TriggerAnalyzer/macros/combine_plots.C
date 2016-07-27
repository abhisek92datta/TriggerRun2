#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
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


void combine_plots() {
	
	 //TFile histofile("Single_Electron_Efficiency_plots.root","recreate");
   	 //histofile.cd();

	TFile *f1 = new TFile("Trigger_efficiency_single_ele_tth.root");
	TEfficiency *tth_h1 = (TEfficiency*)f1->Get("Eff_pt_sing");
	TEfficiency *tth_h2 = (TEfficiency*)f1->Get("Eff_eta_sing");
	TEfficiency *tth_h3 = (TEfficiency*)f1->Get("Eff_pt_cross");
	TEfficiency *tth_h4 = (TEfficiency*)f1->Get("Eff_eta_cross");
	delete f1;
	
	TFile *f2 = new TFile("Trigger_efficiency_single_ele_ttjets.root");
	TEfficiency *ttjets_h1 = (TEfficiency*)f2->Get("Eff_pt_sing");
	TEfficiency *ttjets_h2 = (TEfficiency*)f2->Get("Eff_eta_sing");
	TEfficiency *ttjets_h3 = (TEfficiency*)f2->Get("Eff_pt_cross");
	TEfficiency *ttjets_h4 = (TEfficiency*)f2->Get("Eff_eta_cross");
	delete f2;
	
	TFile *f3 = new TFile("Trigger_efficiency_single_ele_data.root");
	TEfficiency *data_h1 = (TEfficiency*)f3->Get("Eff_pt_sing");
	TEfficiency *data_h2 = (TEfficiency*)f3->Get("Eff_eta_sing");
	TEfficiency *data_h3 = (TEfficiency*)f3->Get("Eff_pt_cross");
	TEfficiency *data_h4 = (TEfficiency*)f3->Get("Eff_eta_cross");
	delete f3;
	
	TFile *f = new TFile("Efficiency_plots_single_electron_triggers.root","recreate");
    //TTree *T = (TTree*)f->Get("T");
	
	//TCanvas *c1 = new TCanvas("c1","test",10,10,600,1000);
	
	//T->SetBranchAddress("tth_single_eff_pt",&tth_h1);
	tth_h1->SetName("tth_single_eff_pt");
	tth_h1->SetTitle("MC ttH : Efficiency vs pT for single trigger;pT (GeV); #epsilon");
	tth_h1->Write();
	
	tth_h2->SetName("tth_single_eff_eta");
	tth_h2->SetTitle("MC ttH : Efficiency vs eta for single trigger;#eta; #epsilon");
	tth_h2->Write();
	
	tth_h3->SetName("tth_cross_eff_pt");
	tth_h3->SetTitle("MC ttH : Efficiency vs pT for cross trigger;pT (GeV); #epsilon");
	tth_h3->Write();
	
	tth_h4->SetName("tth_cross_eff_eta");
	tth_h4->SetTitle("MC ttH : Efficiency vs eta for cross trigger;#eta; #epsilon");
	tth_h4->Write();
	
	ttjets_h1->SetName("ttjets_single_eff_pt");
	ttjets_h1->SetTitle("MC ttjets : Efficiency vs pT for single trigger;pT (GeV); #epsilon");
	ttjets_h1->Write();
	
	ttjets_h2->SetName("ttjets_single_eff_eta");
	ttjets_h2->SetTitle("MC ttjets : Efficiency vs eta for single trigger;#eta; #epsilon");
	ttjets_h2->Write();
	
	ttjets_h3->SetName("ttjets_cross_eff_pt");
	ttjets_h3->SetTitle("MC ttjets : Efficiency vs pT for cross trigger;pT (GeV); #epsilon");
	ttjets_h3->Write();
	
	ttjets_h4->SetName("ttjets_cross_eff_eta");
	ttjets_h4->SetTitle("MC ttjets : Efficiency vs eta for cross trigger;#eta; #epsilon");
	ttjets_h4->Write();
	
	data_h1->SetName("data_single_eff_pt");
	data_h1->SetTitle("DATA : Efficiency vs pT for single trigger;pT (GeV); #epsilon");
	data_h1->Write();
	
	data_h2->SetName("data_single_eff_eta");
	data_h2->SetTitle("DATA : Efficiency vs eta for single trigger;#eta; #epsilon");
	data_h2->Write();
	
	data_h3->SetName("data_cross_eff_pt");
	data_h3->SetTitle("DATA : Efficiency vs pT for cross trigger;pT (GeV); #epsilon");
	data_h3->Write();
	
	data_h4->SetName("data_cross_eff_eta");
	data_h4->SetTitle("DATA : Efficiency vs eta for cross trigger;#eta; #epsilon");
	data_h4->Write();
	
	//T->Fill();
	//tth_h2->Draw();
	
	
	//histofile.Write();
	//histofile.Close();
	f->Close();
	delete f;
	return;
}


