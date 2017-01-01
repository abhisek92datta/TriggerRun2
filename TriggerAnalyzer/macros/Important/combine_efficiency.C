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


void combine_efficiency() {

	TH1::SetDefaultSumw2();
	
	TFile *f1 = new TFile("Efficiency_mc_tthbb.root");
	TEfficiency *tth_h1 = (TEfficiency*)f1->Get("Eff_pt_WPTight_27");
	TEfficiency *tth_h2 = (TEfficiency*)f1->Get("Eff_eta_WPTight_27");
	TEfficiency *tth_h3 = (TEfficiency*)f1->Get("Eff_phi_WPTight_27");
	TEfficiency *tth_h51 = (TEfficiency*)f1->Get("Eff_HT_WPTight_27");
	TEfficiency *tth_h4 = (TEfficiency*)f1->Get("Eff_pt_WPTight_32");
	TEfficiency *tth_h5 = (TEfficiency*)f1->Get("Eff_eta_WPTight_32");
	TEfficiency *tth_h6 = (TEfficiency*)f1->Get("Eff_phi_WPTight_32");
	TEfficiency *tth_h52 = (TEfficiency*)f1->Get("Eff_HT_WPTight_32");
	TEfficiency *tth_h7 = (TEfficiency*)f1->Get("Eff_pt_WPLoose_27_HT200");
	TEfficiency *tth_h8 = (TEfficiency*)f1->Get("Eff_eta_WPLoose_27_HT200");
	TEfficiency *tth_h9 = (TEfficiency*)f1->Get("Eff_phi_WPLoose_27_HT200");
	TEfficiency *tth_h53 = (TEfficiency*)f1->Get("Eff_HT_WPLoose_27_HT200");
	TEfficiency *tth_ho1 = (TEfficiency*)f1->Get("Eff_pt_WPTight_27_OR_WPLoose_27_HT200");
	TEfficiency *tth_ho2 = (TEfficiency*)f1->Get("Eff_eta_WPTight_27_OR_WPLoose_27_HT200");
	TEfficiency *tth_ho3 = (TEfficiency*)f1->Get("Eff_phi_WPTight_27_OR_WPLoose_27_HT200");
	TEfficiency *tth_ho4 = (TEfficiency*)f1->Get("Eff_HT_WPTight_27_OR_WPLoose_27_HT200");
	delete f1;
	
	TFile *f2 = new TFile("Efficiency_mc_ttjets_sl.root");
	TEfficiency *tth_h11 = (TEfficiency*)f2->Get("Eff_pt_WPTight_27");
	TEfficiency *tth_h12 = (TEfficiency*)f2->Get("Eff_eta_WPTight_27");
	TEfficiency *tth_h13 = (TEfficiency*)f2->Get("Eff_phi_WPTight_27");
	TEfficiency *tth_h61 = (TEfficiency*)f2->Get("Eff_HT_WPTight_27");
	TEfficiency *tth_h14 = (TEfficiency*)f2->Get("Eff_pt_WPTight_32");
	TEfficiency *tth_h15 = (TEfficiency*)f2->Get("Eff_eta_WPTight_32");
	TEfficiency *tth_h16 = (TEfficiency*)f2->Get("Eff_phi_WPTight_32");
	TEfficiency *tth_h62 = (TEfficiency*)f2->Get("Eff_HT_WPTight_32");
	TEfficiency *tth_h17 = (TEfficiency*)f2->Get("Eff_pt_WPLoose_27_HT200");
	TEfficiency *tth_h18 = (TEfficiency*)f2->Get("Eff_eta_WPLoose_27_HT200");
	TEfficiency *tth_h19 = (TEfficiency*)f2->Get("Eff_phi_WPLoose_27_HT200");
	TEfficiency *tth_h63 = (TEfficiency*)f2->Get("Eff_HT_WPLoose_27_HT200");
	TEfficiency *tth_ho11 = (TEfficiency*)f2->Get("Eff_pt_WPTight_27_OR_WPLoose_27_HT200");
	TEfficiency *tth_ho12 = (TEfficiency*)f2->Get("Eff_eta_WPTight_27_OR_WPLoose_27_HT200");
	TEfficiency *tth_ho13 = (TEfficiency*)f2->Get("Eff_phi_WPTight_27_OR_WPLoose_27_HT200");
	TEfficiency *tth_ho14 = (TEfficiency*)f2->Get("Eff_HT_WPTight_27_OR_WPLoose_27_HT200");
	delete f2;
	
	TFile *f3 = new TFile("Efficiency_data_2016.root");
	TEfficiency *tth_h21 = (TEfficiency*)f3->Get("Eff_pt_WPTight_27");
	TEfficiency *tth_h22 = (TEfficiency*)f3->Get("Eff_eta_WPTight_27");
	TEfficiency *tth_h23 = (TEfficiency*)f3->Get("Eff_phi_WPTight_27");
	TEfficiency *tth_h71 = (TEfficiency*)f3->Get("Eff_HT_WPTight_27");
	TEfficiency *tth_h24 = (TEfficiency*)f3->Get("Eff_pt_WPTight_32");
	TEfficiency *tth_h25 = (TEfficiency*)f3->Get("Eff_eta_WPTight_32");
	TEfficiency *tth_h26 = (TEfficiency*)f3->Get("Eff_phi_WPTight_32");
	TEfficiency *tth_h72 = (TEfficiency*)f3->Get("Eff_HT_WPTight_32");
	TEfficiency *tth_h27 = (TEfficiency*)f3->Get("Eff_pt_WPLoose_27_HT200");
	TEfficiency *tth_h28 = (TEfficiency*)f3->Get("Eff_eta_WPLoose_27_HT200");
	TEfficiency *tth_h29 = (TEfficiency*)f3->Get("Eff_phi_WPLoose_27_HT200");
	TEfficiency *tth_h73 = (TEfficiency*)f3->Get("Eff_HT_WPLoose_27_HT200");
	TEfficiency *tth_ho21 = (TEfficiency*)f3->Get("Eff_pt_WPTight_27_OR_WPLoose_27_HT200");
	TEfficiency *tth_ho22 = (TEfficiency*)f3->Get("Eff_eta_WPTight_27_OR_WPLoose_27_HT200");
	TEfficiency *tth_ho23 = (TEfficiency*)f3->Get("Eff_phi_WPTight_27_OR_WPLoose_27_HT200");
	TEfficiency *tth_ho24 = (TEfficiency*)f3->Get("Eff_HT_WPTight_27_OR_WPLoose_27_HT200");
	delete f3;

	
	//ttH
	
	tth_h1->SetName("ttH_Eff_pt_WPTight_27");
	tth_h1->SetTitle("MC ttH : Efficiency vs pT for WPTight_27");
		
	tth_h2->SetName("ttH_Eff_eta_WPTight_27");
	tth_h2->SetTitle("MC ttH : Efficiency vs eta for WPTight_27");
	
	tth_h3->SetName("ttH_Eff_phi_WPTight_27");
	tth_h3->SetTitle("MC ttH : Efficiency vs phi for WPTight_27");
	
	tth_h51->SetName("ttH_Eff_HT_WPTight_27");
	tth_h51->SetTitle("MC ttH : Efficiency vs HT for WPTight_27");
	
	tth_h4->SetName("ttH_Eff_pt_WPTight_32");
	tth_h4->SetTitle("MC ttH : Efficiency vs pT for WPTight_32");
	
	tth_h5->SetName("ttH_Eff_eta_WPTight_32");
	tth_h5->SetTitle("MC ttH : Efficiency vs eta for WPTight_32");
	
	tth_h6->SetName("ttH_Eff_phi_WPTight_32");
	tth_h6->SetTitle("MC ttH : Efficiency vs phi for WPTight_32");
	
	tth_h52->SetName("ttH_Eff_HT_WPTight_32");
	tth_h52->SetTitle("MC ttH : Efficiency vs HT for WPTight_32");

	tth_h7->SetName("ttH_Eff_pt_WPLoose_27_HT200");
	tth_h7->SetTitle("MC ttH : Efficiency vs pT for WPLoose_27_HT200");
	
	tth_h8->SetName("ttH_Eff_eta_WPLoose_27_HT200");
	tth_h8->SetTitle("MC ttH : Efficiency vs eta for WPLoose_27_HT200");
	
	tth_h9->SetName("ttH_Eff_phi_WPLoose_27_HT200");
	tth_h9->SetTitle("MC ttH : Efficiency vs phi for WPLoose_27_HT200");
	
	tth_h53->SetName("ttH_Eff_HT_WPLoose_27_HT200");
	tth_h53->SetTitle("MC ttH : Efficiency vs HT for WPLoose_27_HT200");
	
	tth_ho1->SetName("ttH_Eff_pt_WPTight_27_OR_WPLoose_27_HT200");
	tth_ho1->SetTitle("MC ttH : Efficiency vs pT for WPTight_27_OR_WPLoose_27_HT200");
	
	tth_ho2->SetName("ttH_Eff_eta_WPTight_27_OR_WPLoose_27_HT200");
	tth_ho2->SetTitle("MC ttH : Efficiency vs eta for WPTight_27_OR_WPLoose_27_HT200");
	
	tth_ho3->SetName("ttH_Eff_phi_WPTight_27_OR_WPLoose_27_HT200");
	tth_ho3->SetTitle("MC ttH : Efficiency vs phi for WPTight_27_OR_WPLoose_27_HT200");
	
	tth_ho4->SetName("ttH_Eff_HT_WPTight_27_OR_WPLoose_27_HT200");
	tth_ho4->SetTitle("MC ttH : Efficiency vs HT for WPTight_27_OR_WPLoose_27_HT200");
	
	TCanvas *c1 = new TCanvas("c1","test",1100,650);
	//c1->DrawFrame(0,0.85,330,1.03,"MC ttH : Efficiency vs Electron pT; pT (Electron) [GeV] ; (Trig_pass + Cont_trig + Eve_sel)/(Cont_trig + Eve_sel)");
	c1->DrawFrame(0,0.85,330,1.03,"MC ttH : Efficiency vs Electron pT; pT (Electron) [GeV] ; Trigger_pass Efficiency");
	TLegend* leg1 = new TLegend(0.50,0.15,0.70,0.3);
	TLine *line1 = new TLine(0,1,330,1);
  	line1->SetLineStyle(kDashed);
  	line1->SetLineWidth(2);
    leg1->SetFillColor(kWhite);
    leg1->SetFillStyle(1001);
    leg1->AddEntry(tth_h1,"WPTight_27","L");
    leg1->AddEntry(tth_h4,"WPTight_32","L");
    leg1->AddEntry(tth_h7,"WPLoose_27_HT200","L");	
    leg1->AddEntry(tth_ho1,"WPTight_27_OR_WPLoose_27_HT200","L");	
	tth_h1->SetLineColor(kRed);
	tth_h1->SetLineWidth(2);
	tth_h1->Draw("same");
	tth_h7->SetLineColor(kGreen+3);
	tth_h7->SetLineWidth(2);
	tth_h7->Draw("same");
	tth_h4->SetLineColor(kBlue);
	tth_h4->SetLineWidth(2);
	tth_h4->Draw("same");
	tth_ho1->SetLineColor(kOrange);
	tth_ho1->SetLineWidth(2);
	tth_ho1->Draw("same");
	leg1->Draw("same");
	line1->Draw("same");
	c1->Print("tth_eff_pt_all_triggers.png");
	delete c1;
	delete leg1;
	delete line1;
	
	TCanvas *c2 = new TCanvas("c2","test",1100,650);
	//c2->DrawFrame(-3,0.85,3,1.03,"MC ttH : Efficiency vs Electron eta; #eta (Electron) ; (Trigger_pass + Event_sel)/Event_sel");
	c2->DrawFrame(-3,0.85,3,1.03,"MC ttH : Efficiency vs Electron eta; #eta (Electron) ; Trigger_pass Efficiency");
	TLegend* leg2 = new TLegend(0.50,0.15,0.7,0.30);
	TLine *line2 = new TLine(-3,1,3,1);
  	line2->SetLineStyle(kDashed);
  	line2->SetLineWidth(2);
    leg2->SetFillColor(kWhite);
    leg2->SetFillStyle(1001);
    leg2->AddEntry(tth_h2,"WPTight_27","L");
    leg2->AddEntry(tth_h5,"WPTight_32","L");
    leg2->AddEntry(tth_h8,"WPLoose_27_HT200","L");	
    leg2->AddEntry(tth_ho2,"WPTight_27_OR_WPLoose_27_HT200","L");
	tth_h2->SetLineColor(kRed);
	tth_h2->SetLineWidth(2);
	tth_h2->Draw("same");
	tth_h8->SetLineColor(kGreen+3);
	tth_h8->SetLineWidth(2);
	tth_h8->Draw("same");
	tth_h5->SetLineColor(kBlue);
	tth_h5->SetLineWidth(2);
	tth_h5->Draw("same");
	tth_ho2->SetLineColor(kOrange);
	tth_ho2->SetLineWidth(2);
	tth_ho2->Draw("same");
	leg2->Draw("same");
	line2->Draw("same");
	c2->Print("tth_eff_eta_all_triggers.png");
	delete c2;
	delete leg2;
	delete line2;
	
	TCanvas *c3 = new TCanvas("c3","test",1100,650);
	c3->DrawFrame(-3,0.9,3,1.03,"MC ttH : Efficiency vs Electron phi; #phi (Electron) ; Trigger_pass Efficiency");
	TLegend* leg3 = new TLegend(0.65,0.15,0.85,0.30);
	TLine *line3 = new TLine(-3,1,3,1);
  	line3->SetLineStyle(kDashed);
  	line3->SetLineWidth(2);
    leg3->SetFillColor(kWhite);
    leg3->SetFillStyle(1001);
    leg3->AddEntry(tth_h3,"WPTight_27","L");
    leg3->AddEntry(tth_h6,"WPTight_32","L");
    leg3->AddEntry(tth_h9,"WPLoose_27_HT200","L");	
    leg3->AddEntry(tth_ho3,"WPTight_27_OR_WPLoose_27_HT200","L");
	tth_h3->SetLineColor(kRed);
	tth_h3->SetLineWidth(2);
	tth_h3->Draw("same");
	tth_h9->SetLineColor(kGreen+3);
	tth_h9->SetLineWidth(2);
	tth_h9->Draw("same");
	tth_h6->SetLineColor(kBlue);
	tth_h6->SetLineWidth(2);
	tth_h6->Draw("same");
	tth_ho3->SetLineColor(kOrange);
	tth_ho3->SetLineWidth(2);
	tth_ho3->Draw("same");
	leg3->Draw("same");
	line3->Draw("same");
	c3->Print("tth_eff_phi_all_triggers.png");
	delete c3;
	delete leg3;
	delete line3;
	
	TCanvas *c4 = new TCanvas("c4","test",1100,650);
	c4->DrawFrame(100,0.7,1030,1.03,"MC ttH : Efficiency vs HT; HT [GeV] ; Trigger_pass Efficiency");
	TLegend* leg4 = new TLegend(0.65,0.15,0.85,0.3);
	TLine *line4 = new TLine(100,1,1030,1);
  	line4->SetLineStyle(kDashed);
  	line4->SetLineWidth(2);
    leg4->SetFillColor(kWhite);
    leg4->SetFillStyle(1001);
    leg4->AddEntry(tth_h51,"WPTight_27","L");
    leg4->AddEntry(tth_h52,"WPTight_32","L");
    leg4->AddEntry(tth_h53,"WPLoose_27_HT200","L");	
    leg4->AddEntry(tth_ho4,"WPTight_27_OR_WPLoose_27_HT200","L");
	tth_h51->SetLineColor(kRed);
	tth_h51->SetLineWidth(2);
	tth_h51->Draw("same");
	tth_h53->SetLineColor(kGreen+3);
	tth_h53->SetLineWidth(2);
	tth_h53->Draw("same");
	tth_h52->SetLineColor(kBlue);
	tth_h52->SetLineWidth(2);
	tth_h52->Draw("same");
	tth_ho4->SetLineColor(kOrange);
	tth_ho4->SetLineWidth(2);
	tth_ho4->Draw("same");
	leg4->Draw("same");
	line4->Draw("same");
	c4->Print("tth_eff_HT_all_triggers.png");
	delete c4;
	delete leg4;
	delete line4;
	
	
	//background
	
	tth_h11->SetName("bkgd_Eff_pt_WPTight_27");
	tth_h11->SetTitle("MC bkgd : Efficiency vs pT for WPTight_27");
	
	tth_h12->SetName("bkgd_Eff_eta_WPTight_27");
	tth_h12->SetTitle("MC bkgd : Efficiency vs eta for WPTight_27");
	
	tth_h13->SetName("bkgd_Eff_phi_WPTight_27");
	tth_h13->SetTitle("MC bkgd : Efficiency vs phi for WPTight_27");
	
	tth_h61->SetName("bkgd_Eff_HT_WPTight_27");
	tth_h61->SetTitle("MC bkgd : Efficiency vs HT for WPTight_27");
	
	tth_h14->SetName("bkgd_Eff_pt_WPTight_32");
	tth_h14->SetTitle("MC bkgd : Efficiency vs pT for WPTight_32");
	
	tth_h15->SetName("bkgd_Eff_eta_WPTight_32");
	tth_h15->SetTitle("MC bkgd : Efficiency vs eta for WPTight_32");
	
	tth_h16->SetName("bkgd_Eff_phi_WPTight_32");
	tth_h16->SetTitle("MC bkgd : Efficiency vs phi for WPTight_32");
	
	tth_h62->SetName("bkgd_Eff_HT_WPTight_32");
	tth_h62->SetTitle("MC bkgd : Efficiency vs HT for WPTight_32");
	
	tth_h17->SetName("bkgd_Eff_pt_WPLoose_27_HT200");
	tth_h17->SetTitle("MC bkgd : Efficiency vs pT for WPLoose_27_HT200");
	
	tth_h18->SetName("bkgd_Eff_eta_WPLoose_27_HT200");
	tth_h18->SetTitle("MC bkgd : Efficiency vs eta for WPLoose_27_HT200");
	
	tth_h19->SetName("bkgd_Eff_phi_WPLoose_27_HT200");
	tth_h19->SetTitle("MC bkgd : Efficiency vs phi for WPLoose_27_HT200");
	
	tth_h63->SetName("bkgd_Eff_HT_WPLoose_27_HT200");
	tth_h63->SetTitle("MC bkgd : Efficiency vs HT for WPLoose_27_HT200");
	
	tth_ho11->SetName("bkgd_Eff_pt_WPTight_27_OR_WPLoose_27_HT200");
	tth_ho11->SetTitle("MC bkgd : Efficiency vs pT for WPTight_27_OR_WPLoose_27_HT200");
	
	tth_ho12->SetName("bkgd_Eff_eta_WPTight_27_OR_WPLoose_27_HT200");
	tth_ho12->SetTitle("MC bkgd : Efficiency vs eta for WPTight_27_OR_WPLoose_27_HT200");
	
	tth_ho13->SetName("bkgd_Eff_phi_WPTight_27_OR_WPLoose_27_HT200");
	tth_ho13->SetTitle("MC bkgd : Efficiency vs phi for WPTight_27_OR_WPLoose_27_HT200");
	
	tth_ho14->SetName("bkgd_Eff_HT_WPTight_27_OR_WPLoose_27_HT200");
	tth_ho14->SetTitle("MC bkgd : Efficiency vs HT for WPTight_27_OR_WPLoose_27_HT200");
	
	TCanvas *c11 = new TCanvas("c11","test",1100,650);
	c11->DrawFrame(0,0.85,330,1.03,"MC background : Efficiency vs Electron pT; pT (Electron) [GeV] ; Trigger_pass Efficiency");
	TLegend* leg11 = new TLegend(0.65,0.15,0.85,0.3);
	TLine *line11 = new TLine(0,1,330,1);
  	line11->SetLineStyle(kDashed);
  	line11->SetLineWidth(2);
    leg11->SetFillColor(kWhite);
    leg11->SetFillStyle(1001);
    leg11->AddEntry(tth_h11,"WPTight_27","L");
    leg11->AddEntry(tth_h14,"WPTight_32","L");
    leg11->AddEntry(tth_h17,"WPLoose_27_HT200","L");	
    leg11->AddEntry(tth_ho11,"WPTight_27_OR_WPLoose_27_HT200","L");	
	tth_h11->SetLineColor(kRed);
	tth_h11->SetLineWidth(2);
	tth_h11->Draw("same");
	tth_h17->SetLineColor(kGreen+3);
	tth_h17->SetLineWidth(2);
	tth_h17->Draw("same");
	tth_h14->SetLineColor(kBlue);
	tth_h14->SetLineWidth(2);
	tth_h14->Draw("same");
	tth_ho11->SetLineColor(kOrange);
	tth_ho11->SetLineWidth(2);
	tth_ho11->Draw("same");
	leg11->Draw("same");
	line11->Draw("same");
	c11->Print("bkgd_eff_pt_all_triggers.png");
	delete c11;
	delete leg11;
	delete line11;
	
	
	TCanvas *c12 = new TCanvas("c12","test",1100,650);
	c12->DrawFrame(-3,0.85,3,1.03,"MC background : Efficiency vs Electron eta; #eta (Electron) ; Trigger_pass Efficiency");
	TLegend* leg12 = new TLegend(0.50,0.15,0.70,0.3);
	TLine *line12 = new TLine(-3,1,3,1);
  	line12->SetLineStyle(kDashed);
  	line12->SetLineWidth(2);
    leg12->SetFillColor(kWhite);
    leg12->SetFillStyle(1001);
    leg12->AddEntry(tth_h12,"WPTight_27","L");
    leg12->AddEntry(tth_h15,"WPTight_32","L");
    leg12->AddEntry(tth_h18,"WPLoose_27_HT200","L");
    leg12->AddEntry(tth_ho12,"WPTight_27_OR_WPLoose_27_HT200","L");		
	tth_h12->SetLineColor(kRed);
	tth_h12->SetLineWidth(2);
	tth_h12->Draw("same");
	tth_h18->SetLineColor(kGreen+3);
	tth_h18->SetLineWidth(2);
	tth_h18->Draw("same");
	tth_h15->SetLineColor(kBlue);
	tth_h15->SetLineWidth(2);
	tth_h15->Draw("same");
	tth_ho12->SetLineColor(kOrange);
	tth_ho12->SetLineWidth(2);
	tth_ho12->Draw("same");
	leg12->Draw("same");
	line12->Draw("same");
	c12->Print("bkgd_eff_eta_all_triggers.png");
	delete c12;
	delete leg12;
	delete line12;
	
	TCanvas *c13 = new TCanvas("c13","test",1100,650);
	c13->DrawFrame(-3,0.9,3,1.03,"MC background : Efficiency vs Electron phi; #phi (Electron) ; Trigger_pass Efficiency");
	TLegend* leg13 = new TLegend(0.65,0.15,0.85,0.3);
	TLine *line13 = new TLine(-3,1,3,1);
  	line13->SetLineStyle(kDashed);
  	line13->SetLineWidth(2);
    leg13->SetFillColor(kWhite);
    leg13->SetFillStyle(1001);
    leg13->AddEntry(tth_h13,"WPTight_27","L");
    leg13->AddEntry(tth_h16,"WPTight_32","L");
    leg13->AddEntry(tth_h19,"WPLoose_27_HT200","L");	
    leg13->AddEntry(tth_ho13,"WPTight_27_OR_WPLoose_27_HT200","L");	
	tth_h13->SetLineColor(kRed);
	tth_h13->SetLineWidth(2);
	tth_h13->Draw("same");
	tth_h19->SetLineColor(kGreen+3);
	tth_h19->SetLineWidth(2);
	tth_h19->Draw("same");
	tth_h16->SetLineColor(kBlue);
	tth_h16->SetLineWidth(2);
	tth_h16->Draw("same");
	tth_ho13->SetLineColor(kOrange);
	tth_ho13->SetLineWidth(2);
	tth_ho13->Draw("same");
	leg13->Draw("same");
	line13->Draw("same");
	c13->Print("bkgd_eff_phi_all_triggers.png");
	delete c13;
	delete leg13;
	delete line13;
	
	TCanvas *c14 = new TCanvas("c14","test",1100,650);
	c14->DrawFrame(100,0.75,1030,1.03,"MC background : Efficiency vs HT; HT [GeV] ; Trigger_pass Efficiency");
	TLegend* leg14 = new TLegend(0.65,0.15,0.85,0.3);
	TLine *line14 = new TLine(100,1,1030,1);
  	line14->SetLineStyle(kDashed);
  	line14->SetLineWidth(2);
    leg14->SetFillColor(kWhite);
    leg14->SetFillStyle(1001);
    leg14->AddEntry(tth_h61,"WPTight_27","L");
    leg14->AddEntry(tth_h62,"WPTight_32","L");
    leg14->AddEntry(tth_h63,"WPLoose_27_HT200","L");
    leg14->AddEntry(tth_ho14,"WPTight_27_OR_WPLoose_27_HT200","L");		
	tth_h61->SetLineColor(kRed);
	tth_h61->SetLineWidth(2);
	tth_h61->Draw("same");
	tth_h63->SetLineColor(kGreen+3);
	tth_h63->SetLineWidth(2);
	tth_h63->Draw("same");
	tth_h62->SetLineColor(kBlue);
	tth_h62->SetLineWidth(2);
	tth_h62->Draw("same");
	tth_ho14->SetLineColor(kOrange);
	tth_ho14->SetLineWidth(2);
	tth_ho14->Draw("same");
	leg14->Draw("same");
	line14->Draw("same");
	c14->Print("bkgd_eff_HT_all_triggers.png");
	delete c14;
	delete leg14;
	delete line14;
	
	//Data
	
	tth_h21->SetName("data_Eff_pt_WPTight_27");
	tth_h21->SetTitle("DATA : Efficiency vs pT for WPTight_27");
	
	tth_h22->SetName("data_Eff_eta_WPTight_27");
	tth_h22->SetTitle("DATA : Efficiency vs eta for WPTight_27");
	
	tth_h23->SetName("data_Eff_phi_WPTight_27");
	tth_h23->SetTitle("DATA : Efficiency vs phi for WPTight_27");
	
	tth_h71->SetName("data_Eff_HT_WPTight_27");
	tth_h71->SetTitle("DATA : Efficiency vs HT for WPTight_27");
	
	tth_h24->SetName("data_Eff_pt_WPTight_32");
	tth_h24->SetTitle("DATA : Efficiency vs pT for WPTight_32");
	
	tth_h25->SetName("data_Eff_eta_WPTight_32");
	tth_h25->SetTitle("DATA : Efficiency vs eta for WPTight_32");
	
	tth_h72->SetName("data_Eff_HT_WPTight_32");
	tth_h72->SetTitle("DATA : Efficiency vs HT for WPTight_32");
	
	tth_h26->SetName("data_Eff_phi_WPTight_32");
	tth_h26->SetTitle("DATA : Efficiency vs phi for WPTight_32");
	
	tth_h27->SetName("data_Eff_pt_WPLoose_27_HT200");
	tth_h27->SetTitle("DATA : Efficiency vs pT for WPLoose_27_HT200");
	
	tth_h28->SetName("data_Eff_eta_WPLoose_27_HT200");
	tth_h28->SetTitle("DATA : Efficiency vs eta for WPLoose_27_HT200");
	
	tth_h29->SetName("data_Eff_phi_WPLoose_27_HT200");
	tth_h29->SetTitle("DATA : Efficiency vs phi for WPLoose_27_HT200");
	
	tth_h73->SetName("data_Eff_HT_WPLoose_27_HT200");
	tth_h73->SetTitle("DATA : Efficiency vs HT for WPLoose_27_HT200");
	
	tth_ho21->SetName("data_Eff_pt_WPTight_27_OR_WPLoose_27_HT200");
	tth_ho21->SetTitle("DATA : Efficiency vs pT for WPTight_27_OR_WPLoose_27_HT200");
	
	tth_ho22->SetName("data_Eff_eta_WPTight_27_OR_WPLoose_27_HT200");
	tth_ho22->SetTitle("DATA : Efficiency vs eta for WPTight_27_OR_WPLoose_27_HT200");
	
	tth_ho23->SetName("data_Eff_phi_WPTight_27_OR_WPLoose_27_HT200");
	tth_ho23->SetTitle("DATA : Efficiency vs phi for WPTight_27_OR_WPLoose_27_HT200");
	
	tth_ho24->SetName("data_Eff_HT_WPTight_27_OR_WPLoose_27_HT200");
	tth_ho24->SetTitle("DATA : Efficiency vs HT for WPTight_27_OR_WPLoose_27_HT200");
	
	TCanvas *c21 = new TCanvas("c21","test",1100,650);
	c21->DrawFrame(0,0.75,330,1.03,"DATA : Efficiency vs Electron pT; pT (Electron) [GeV] ; Trigger_pass Efficiency");
	TLegend* leg21 = new TLegend(0.3,0.15,0.5,0.3);
	TLine *line21 = new TLine(0,1,330,1);
  	line21->SetLineStyle(kDashed);
  	line21->SetLineWidth(2);
    leg21->SetFillColor(kWhite);
    leg21->SetFillStyle(1001);
    leg21->AddEntry(tth_h21,"WPTight_27","L");
    leg21->AddEntry(tth_h24,"WPTight_32","L");
    leg21->AddEntry(tth_h27,"WPLoose_27_HT200","L");	
    leg21->AddEntry(tth_ho21,"WPTight_27_OR_WPLoose_27_HT200","L");	
	tth_h21->SetLineColor(kRed);
	tth_h21->SetLineWidth(2);
	tth_h21->Draw("same");
	tth_h27->SetLineColor(kGreen+3);
	tth_h27->SetLineWidth(2);
	tth_h27->Draw("same");
	tth_h24->SetLineColor(kBlue);
	tth_h24->SetLineWidth(2);
	tth_h24->Draw("same");
	tth_ho21->SetLineColor(kOrange);
	tth_ho21->SetLineWidth(2);
	tth_ho21->Draw("same");
	leg21->Draw("same");
	line21->Draw("same");
	c21->Print("data_eff_pt_all_triggers.png");
	delete c21;
	delete leg21;
	delete line21;
	
	
	TCanvas *c22 = new TCanvas("c22","test",1100,650);
	c22->DrawFrame(-3,0.8,3,1.03,"DATA : Efficiency vs Electron eta; #eta (Electron) ; Trigger_pass Efficiency");
	TLegend* leg22 = new TLegend(0.45,0.15,0.65,0.3);
	TLine *line22 = new TLine(-3,1,3,1);
  	line22->SetLineStyle(kDashed);
  	line22->SetLineWidth(2);
    leg22->SetFillColor(kWhite);
    leg22->SetFillStyle(1001);
    leg22->AddEntry(tth_h22,"WPTight_27","L");
    leg22->AddEntry(tth_h25,"WPTight_32","L");
    leg22->AddEntry(tth_h28,"WPLoose_27_HT200","L");
    leg22->AddEntry(tth_ho22,"WPTight_27_OR_WPLoose_27_HT200","L");		
	tth_h22->SetLineColor(kRed);
	tth_h22->SetLineWidth(2);
	tth_h22->Draw("same");
	tth_h28->SetLineColor(kGreen+3);
	tth_h28->SetLineWidth(2);
	tth_h28->Draw("same");
	tth_h25->SetLineColor(kBlue);
	tth_h25->SetLineWidth(2);
	tth_h25->Draw("same");
	tth_ho22->SetLineColor(kOrange);
	tth_ho22->SetLineWidth(2);
	tth_ho22->Draw("same");
	leg22->Draw("same");
	line22->Draw("same");
	c22->Print("data_eff_eta_all_triggers.png");
	delete c22;
	delete leg22;
	delete line22;
	
	TCanvas *c23 = new TCanvas("c23","test",1100,650);
	c23->DrawFrame(-3,0.85,3,1.03,"DATA : Efficiency vs Electron phi; #phi (Electron) ; Trigger_pass Efficiency");
	TLegend* leg23 = new TLegend(0.65,0.15,0.85,0.3);
	TLine *line23 = new TLine(-3,1,3,1);
  	line23->SetLineStyle(kDashed);
  	line23->SetLineWidth(2);
    leg23->SetFillColor(kWhite);
    leg23->SetFillStyle(1001);
    leg23->AddEntry(tth_h23,"WPTight_27","L");
    leg23->AddEntry(tth_h26,"WPTight_32","L");
    leg23->AddEntry(tth_h29,"WPLoose_27_HT200","L");	
    leg23->AddEntry(tth_ho23,"WPTight_27_OR_WPLoose_27_HT200","L");	
	tth_h23->SetLineColor(kRed);
	tth_h23->SetLineWidth(2);
	tth_h23->Draw("same");
	tth_h29->SetLineColor(kGreen+3);
	tth_h29->SetLineWidth(2);
	tth_h29->Draw("same");
	tth_h26->SetLineColor(kBlue);
	tth_h26->SetLineWidth(2);
	tth_h26->Draw("same");
	tth_ho23->SetLineColor(kOrange);
	tth_ho23->SetLineWidth(2);
	tth_ho23->Draw("same");
	leg23->Draw("same");
	line23->Draw("same");
	c23->Print("data_eff_phi_all_triggers.png");
	delete c23;
	delete leg23;
	delete line23;
	
	TCanvas *c24 = new TCanvas("c24","test",1100,650);
	c24->DrawFrame(100,0.6,1030,1.03,"DATA : Efficiency vs HT; HT [GeV] ; Trigger_pass Efficiency");
	TLegend* leg24 = new TLegend(0.3,0.15,0.5,0.3);
	TLine *line24 = new TLine(100,1,1030,1);
  	line24->SetLineStyle(kDashed);
  	line24->SetLineWidth(2);
    leg24->SetFillColor(kWhite);
    leg24->SetFillStyle(1001);
    leg24->AddEntry(tth_h71,"WPTight_27","L");
    leg24->AddEntry(tth_h72,"WPTight_32","L");
    leg24->AddEntry(tth_h73,"WPLoose_27_HT200","L");	
    leg24->AddEntry(tth_ho24,"WPTight_27_OR_WPLoose_27_HT200","L");	
	tth_h71->SetLineColor(kRed);
	tth_h71->SetLineWidth(2);
	tth_h71->Draw("same");
	tth_h73->SetLineColor(kGreen+3);
	tth_h73->SetLineWidth(2);
	tth_h73->Draw("same");
	tth_h72->SetLineColor(kBlue);
	tth_h72->SetLineWidth(2);
	tth_h72->Draw("same");
	tth_ho24->SetLineColor(kOrange);
	tth_ho24->SetLineWidth(2);
	tth_ho24->Draw("same");
	leg24->Draw("same");
	line24->Draw("same");
	c24->Print("data_eff_HT_all_triggers.png");
	delete c24;
	delete leg24;
	delete line24;
	
	return;
}


