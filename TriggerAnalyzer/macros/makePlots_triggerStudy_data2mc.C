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

TH1D* divide_ratio_plots( TH1D* h_numerator_numerator, TH1D* h_numerator_denominator, TH1D* h_denominator_numerator, TH1D* h_denominator_denominator );

TH2D* divide_ratio_plots_2d( TH2D* h_numerator_numerator, TH2D* h_numerator_denominator, TH2D* h_denominator_numerator, TH2D* h_denominator_denominator );

//*****************************************************************************


void makePlots_triggerStudy_data2mc( TString dirpostfix_ = "", bool printPDF_ = false, bool plotEleEff_ = false, bool plot2dSF_ = false, bool plot1dDist_ = false, double sfMC_ = 1., bool renormMC_ = false ){

  TH1::SetDefaultSumw2();


  int NumSamples = 3;
  TFile* file[NumSamples];
  file[0] = new TFile("HistoFiles/triggerStudy_data2mc_treeReader_SingleElectron_Run2015D_05Oct2015_PromptRecov4_histo.root");
  file[1] = new TFile("HistoFiles/triggerStudy_data2mc_treeReader_ZJets_M50_histo.root");
  file[2] = new TFile("HistoFiles/triggerStudy_data2mc_treeReader_ttbar_histo.root");

  int file_data = 0;
  int file_dy = 1;
  int file_tt = 2;

  std::vector<TString> histLabels(NumSamples);
  histLabels[file_data] = "Data";
  histLabels[file_dy] = "Z+jets";
  histLabels[file_tt] = "ttbar";

  Color_t color[8];
  color[file_data] = kBlack;
  color[file_dy] = kAzure-2;
  color[file_tt] = kRed-4;



  TString dirprefix = "Images/Images_2015_12_18_triggerStudy_data2mc" + dirpostfix_;
  if( renormMC_ ) dirprefix += "_renormMC";

  dirprefix += "/";


  struct stat st;
  if( stat(dirprefix.Data(),&st) != 0 )  mkdir(dirprefix.Data(),0777);


 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<TString> histoname1;
  histoname1.push_back("h_probe_pt");
  histoname1.push_back("h_probe_eta");
  histoname1.push_back("h_probe_phi");
  histoname1.push_back("h_probe_numPVs");
  histoname1.push_back("h_probe_HT30");

  histoname1.push_back("h_2e_HT30");
  histoname1.push_back("h_1e4j2t_HT30");
  histoname1.push_back("h_2e_HT30nocc");
  histoname1.push_back("h_1e4j2t_HT30nocc");


  std::vector<TString> histoname2;
  histoname2.push_back("h_ee_diLepMass");
  histoname2.push_back("h_numPVs_wgt");
  histoname2.push_back("h_numPVs_noPUwgt");
  histoname2.push_back("h_2e_L1HTT");
  histoname2.push_back("h_2e_HT30");
  histoname2.push_back("h_2e_HT30_L1HTT100_passHLTEle27HT200");
  histoname2.push_back("h_2e_numJet");
  histoname2.push_back("h_2e_numBtag");



 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TGaxis::SetMaxDigits(4);

  TString lumiinfo = "2.6 fb^{-1} (13 TeV)";
  TLatex LumiInfoLatex(0.70, 0.91, lumiinfo);
  LumiInfoLatex.SetNDC(); LumiInfoLatex.SetTextFont(42);
  LumiInfoLatex.SetTextSize(0.04);

  TString cmsinfo =   "CMS";
  TLatex CMSInfoLatex(0.185, 0.91, cmsinfo);
  CMSInfoLatex.SetNDC(); CMSInfoLatex.SetTextFont(42);
  CMSInfoLatex.SetTextFont(61);
  CMSInfoLatex.SetTextSize(0.055);

  TString publishinfo =   "Preliminary";
  TLatex PublishInfoLatex(0.285, 0.91, publishinfo);
  PublishInfoLatex.SetNDC();
  PublishInfoLatex.SetTextFont(52);
  PublishInfoLatex.SetTextSize(0.045);


  TString plotname;

  TCanvas* myC1 = new TCanvas("myC1", "myC1", 600,700);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  Float_t small = 1.e-5;
  myC1->Divide(1,2,small,small);
  const float padding=1e-5; const float ydivide=0.3;
  myC1->GetPad(1)->SetPad( padding, ydivide + padding , 1-padding, 1-padding);
  myC1->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
  myC1->GetPad(1)->SetLeftMargin(.11);
  myC1->GetPad(2)->SetLeftMargin(.11);
  myC1->GetPad(1)->SetRightMargin(.05);
  myC1->GetPad(2)->SetRightMargin(.05);
  myC1->GetPad(1)->SetBottomMargin(.3);
  myC1->GetPad(2)->SetBottomMargin(.3);
  myC1->GetPad(1)->Modified();
  myC1->GetPad(2)->Modified();
  myC1->cd(1);
  gPad->SetBottomMargin(small);
  gPad->Modified();


  if( plotEleEff_ ){
    for( int i=0; i<int(histoname1.size()); i++ ){

      TString temp = histoname1[i];
    
      TString temp_mh = temp;
      temp_mh.ReplaceAll("h_","");

      // TLegend *legend = new TLegend(0.1,0.91,0.9,0.99);
      //TLegend *legend = new TLegend(0.2,0.91,0.9,0.99);
      //TLegend *legend = new TLegend(0.1,0.91,0.9,0.99);
      TLegend *legend = new TLegend(0.2,0.83,0.9,0.89);

      legend->SetFillColor(kWhite);
      legend->SetLineColor(kWhite);
      legend->SetShadowColor(kWhite);
      legend->SetTextFont(42);
      legend->SetTextSize(0.04);

      legend->SetNColumns(2);

      int rebin = 1;

      //if( temp.Contains("_HT30") ) rebin = 10;

      TString s_all = temp+"_all";
      TString s_l1t = temp+"_pL1T";
      TString s_hlt = temp+"_pHLT";

      if( temp.Contains("h_2e_HT30") || temp.Contains("h_1e4j2t_HT30") ){
	s_all = temp+"";
	s_l1t = temp+"_L1HTT125";
	s_hlt = temp+"_L1HTT125_passHLTEle27HT200";
      }
      if( temp.Contains("h_HT30_v2") ){
	s_all = temp+"";
	s_l1t = temp+"_L1HTT";
	s_hlt = temp+"_L1HTT_HLTHT";
      }

      TH1D* h_data_all = (TH1D*)file[file_data]->Get(s_all)->Clone(s_all+"_data");
      TH1D* h_data_l1t = (TH1D*)file[file_data]->Get(s_l1t)->Clone(s_l1t+"_data");
      TH1D* h_data_hlt = (TH1D*)file[file_data]->Get(s_hlt)->Clone(s_hlt+"_data");

      TH1D* h_mc_tt_all = (TH1D*)file[file_tt]->Get(s_all)->Clone(s_all+"_mc_tt");
      TH1D* h_mc_tt_l1t = (TH1D*)file[file_tt]->Get(s_l1t)->Clone(s_l1t+"_mc_tt");
      TH1D* h_mc_tt_hlt = (TH1D*)file[file_tt]->Get(s_hlt)->Clone(s_hlt+"_mc_tt");

      TH1D* h_mc_dy_all = (TH1D*)file[file_dy]->Get(s_all)->Clone(s_all+"_mc_dy");
      TH1D* h_mc_dy_l1t = (TH1D*)file[file_dy]->Get(s_l1t)->Clone(s_l1t+"_mc_dy");
      TH1D* h_mc_dy_hlt = (TH1D*)file[file_dy]->Get(s_hlt)->Clone(s_hlt+"_mc_dy");

      TH1D* h_mc_all = (TH1D*)h_mc_tt_all->Clone(s_all+"_mc");
      TH1D* h_mc_l1t = (TH1D*)h_mc_tt_l1t->Clone(s_l1t+"_mc");
      TH1D* h_mc_hlt = (TH1D*)h_mc_tt_hlt->Clone(s_hlt+"_mc");

      h_mc_all->Add(h_mc_dy_all);
      h_mc_l1t->Add(h_mc_dy_l1t);
      h_mc_hlt->Add(h_mc_dy_hlt);


      h_data_all->Rebin(rebin);
      h_data_l1t->Rebin(rebin);
      h_data_hlt->Rebin(rebin);

      h_mc_all->Rebin(rebin);
      h_mc_l1t->Rebin(rebin);
      h_mc_hlt->Rebin(rebin);

      double eps = 0.001;
      int nbins = h_data_all->GetNbinsX();

      for( int iBin=0; iBin<nbins+1; iBin++ ){
	double content1 = h_mc_hlt->GetBinContent(iBin+1);
	double content2 = h_mc_l1t->GetBinContent(iBin+1);
	double content3 = h_mc_all->GetBinContent(iBin+1);

	double diff12 = content1 - content2;
	if( diff12>0 ) content2 += diff12 + eps;

	double diff23 = content2 - content3;
	if( diff23>0 ) content3 += diff23 + eps;

	h_mc_l1t->SetBinContent(iBin+1,content2);
	h_mc_all->SetBinContent(iBin+1,content3);
      }


      TEfficiency* h_eff_data_l1t_all = new TEfficiency(*h_data_l1t, *h_data_all);
      TEfficiency* h_eff_data_hlt_all = new TEfficiency(*h_data_hlt, *h_data_all);
      TEfficiency* h_eff_data_hlt_l1t = new TEfficiency(*h_data_hlt, *h_data_l1t);

      TEfficiency* h_eff_mc_l1t_all = new TEfficiency(*h_mc_l1t, *h_mc_all);
      TEfficiency* h_eff_mc_hlt_all = new TEfficiency(*h_mc_hlt, *h_mc_all);
      TEfficiency* h_eff_mc_hlt_l1t = new TEfficiency(*h_mc_hlt, *h_mc_l1t);

      h_eff_data_l1t_all->SetLineColor(color[0]);
      h_eff_data_l1t_all->SetMarkerColor(color[0]);
      h_eff_data_l1t_all->SetMarkerStyle(20);

      h_eff_data_hlt_all->SetLineColor(color[0]);
      h_eff_data_hlt_all->SetMarkerColor(color[0]);
      h_eff_data_hlt_all->SetMarkerStyle(20);

      h_eff_data_hlt_l1t->SetLineColor(color[0]);
      h_eff_data_hlt_l1t->SetMarkerColor(color[0]);
      h_eff_data_hlt_l1t->SetMarkerStyle(20);

      h_eff_mc_l1t_all->SetLineColor(color[2]);
      h_eff_mc_l1t_all->SetMarkerColor(color[2]);
      h_eff_mc_l1t_all->SetMarkerStyle(20);

      h_eff_mc_hlt_all->SetLineColor(color[2]);
      h_eff_mc_hlt_all->SetMarkerColor(color[2]);
      h_eff_mc_hlt_all->SetMarkerStyle(20);

      h_eff_mc_hlt_l1t->SetLineColor(color[2]);
      h_eff_mc_hlt_l1t->SetMarkerColor(color[2]);
      h_eff_mc_hlt_l1t->SetMarkerStyle(20);


      legend->AddEntry(h_eff_data_l1t_all,histLabels[file_data],"p");
      legend->AddEntry(h_eff_mc_l1t_all,"MC","p");


      TH1D* h_ratio_eff_l1t_all = divide_ratio_plots(h_data_l1t, h_data_all, h_mc_l1t, h_mc_all);
      TH1D* h_ratio_eff_hlt_all = divide_ratio_plots(h_data_hlt, h_data_all, h_mc_hlt, h_mc_all);
      TH1D* h_ratio_eff_hlt_l1t = divide_ratio_plots(h_data_hlt, h_data_l1t, h_mc_hlt, h_mc_l1t);

      h_ratio_eff_l1t_all->SetMarkerStyle(20);
      h_ratio_eff_hlt_all->SetMarkerStyle(20);
      h_ratio_eff_hlt_l1t->SetMarkerStyle(20);


      h_data_all->SetStats(0);

      //c3->SetTopMargin(.05);
      //c1->SetRightMargin(.05);


      double xmin = h_data_all->GetBinLowEdge(1);
      double xmax = h_data_all->GetBinLowEdge(nbins) + h_data_all->GetBinWidth(nbins);


      TH1D* myRatio = (TH1D*)h_data_all->Clone("ratio_"+temp);
      //new TH1D("ratio", "", nbins, xmin, xmax );

      myRatio->SetStats(0);
      myRatio->Sumw2();
      myRatio->SetLineColor(kBlack);
      myRatio->SetMarkerColor(kBlack);
      //myRatio->Divide(hist[bin_one],hist[bin_two]);

      //myRatio->GetYaxis()->SetNdivisions(50000+404);
      // double ratioMax = 1.6;
      // double ratioMin = 0.5;
      // myRatio->GetYaxis()->SetNdivisions(50000+204);
      double ratioMax = 1.11;
      double ratioMin = 0.90;
      myRatio->GetYaxis()->SetNdivisions(50000+204);
      myRatio->GetYaxis()->SetLabelSize(0.1); //make y label bigger
      myRatio->GetXaxis()->SetLabelSize(0.1); //make y label bigger
      myRatio->GetXaxis()->SetTitleOffset(1.1);
      myRatio->GetXaxis()->SetTitle(h_data_all->GetXaxis()->GetTitle()); //make y label bigger
      myRatio->GetXaxis()->SetLabelSize(0.12);
      myRatio->GetXaxis()->SetLabelOffset(0.04);
      myRatio->GetXaxis()->SetTitleSize(0.12);
      myRatio->GetYaxis()->SetTitle("Data/MC");
      myRatio->GetYaxis()->SetTitleSize(0.09);
      myRatio->GetYaxis()->SetTitleOffset(.55);
      myC1->cd(2);
      gPad->SetTopMargin(small);
      gPad->SetTickx();
      gPad->Modified();

      myRatio->SetMinimum(ratioMin);
      myRatio->SetMaximum(ratioMax);


      myRatio->GetYaxis()->CenterTitle(kTRUE);


      double maxPt = 499.9999;
      double minPt = 0.000001;

      h_data_all->GetYaxis()->SetTitleOffset(1.0);
      h_data_all->GetYaxis()->SetTitleSize(0.05);

      h_data_all->GetYaxis()->SetTitle("Efficiency");

      h_data_all->GetYaxis()->SetRangeUser(0.,1.15);


      bool rerange = false;
      //if( (temp.Contains("_pt")) ){ rerange=true; xmin = 10.0001; xmax = 199.9; }
      if( (temp.Contains("_pt")) ){ rerange=true; xmin = 30+0.0001; xmax = 100-0.00001; }
      //if( (temp.Contains("_HT30")) ){ rerange=true; xmin = 0+0.0001; xmax = 1000-0.00001; }

      if( rerange ){
	h_data_all->GetXaxis()->SetRangeUser(xmin,xmax);
	myRatio->GetXaxis()->SetRangeUser(xmin,xmax);
      }

      if( temp.Contains("_pt") ){
	h_data_all->GetYaxis()->SetRangeUser(0.5,1.05);
      }

      if( temp.Contains("_numPVs") ){
	h_data_all->GetYaxis()->SetRangeUser(0.45,1.05);
      }

      // if( temp.Contains("_eta") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);
      // if( temp.Contains("_numGenPVs") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);
      // if( temp.Contains("_numJets") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);


      TLine* myLine;
      myLine = new TLine(xmin, 1, xmax, 1);


      //// l1t_all
      myC1->cd(1);
      h_data_all->Draw("axis");
      h_eff_data_l1t_all->Draw("pe1same");
      h_eff_mc_l1t_all->Draw("pe1same");
      legend->Draw();
      LumiInfoLatex.Draw();
      CMSInfoLatex.Draw();
      PublishInfoLatex.Draw();

      myC1->cd(2);
      myRatio->SetLineWidth(2);
      myRatio->Draw("axis");
      myLine->Draw();

      h_ratio_eff_l1t_all->Draw("pe1same");

      plotname = dirprefix + temp_mh + "_data2mc_eff_l1t_all_lin.png";
      myC1->Print(plotname);

      plotname = dirprefix + temp_mh + "_data2mc_eff_l1t_all_lin.pdf";
      if( printPDF_ ) myC1->Print(plotname);


      //// hlt_all
      myC1->cd(1);
      h_data_all->Draw("axis");
      h_eff_data_hlt_all->Draw("pe1same");
      h_eff_mc_hlt_all->Draw("pe1same");
      legend->Draw();
      LumiInfoLatex.Draw();
      CMSInfoLatex.Draw();
      PublishInfoLatex.Draw();

      myC1->cd(2);
      myRatio->SetLineWidth(2);
      myRatio->Draw("axis");
      myLine->Draw();

      h_ratio_eff_hlt_all->Draw("pe1same");

      plotname = dirprefix + temp_mh + "_data2mc_eff_hlt_all_lin.png";
      myC1->Print(plotname);

      plotname = dirprefix + temp_mh + "_data2mc_eff_hlt_all_lin.pdf";
      if( printPDF_ ) myC1->Print(plotname);


      //// hlt_l1t
      myC1->cd(1);
      h_data_all->Draw("axis");
      h_eff_data_hlt_l1t->Draw("pe1same");
      h_eff_mc_hlt_l1t->Draw("pe1same");
      legend->Draw();
      LumiInfoLatex.Draw();
      CMSInfoLatex.Draw();
      PublishInfoLatex.Draw();

      myC1->cd(2);
      myRatio->SetLineWidth(2);
      myRatio->Draw("axis");
      myLine->Draw();

      h_ratio_eff_hlt_l1t->Draw("pe1same");

      plotname = dirprefix + temp_mh + "_data2mc_eff_hlt_l1t_lin.png";
      myC1->Print(plotname);

      plotname = dirprefix + temp_mh + "_data2mc_eff_hlt_l1t_lin.pdf";
      if( printPDF_ ) myC1->Print(plotname);



      delete myRatio;
      delete myLine;
      delete legend;
    } // end loop on hists
  }// end if plotEleEff_

  ////////////////////////////////////////////////////////////////////////////


  if( plot1dDist_ ){
    for( int i=0; i<int(histoname2.size()); i++ ){

      TString temp = histoname2[i];
    
      TString temp_mh = temp;
      temp_mh.ReplaceAll("h_","");

      //TLegend *legend = new TLegend(0.1,0.91,0.9,0.99);
      TLegend *legend = new TLegend(0.2,0.83,0.88,0.89);

      legend->SetFillColor(kWhite);
      legend->SetLineColor(kWhite);
      legend->SetShadowColor(kWhite);
      legend->SetTextFont(42);
      legend->SetTextSize(0.04);

      legend->SetNColumns(3);

      int rebin = 1;

      if( temp.Contains("_mass") ) rebin = 2;
      if( temp.Contains("_L1HTT") ) rebin = 2;
      if( temp.Contains("_met") ) rebin = 5;
      if( temp.Contains("_jet_") && temp.Contains("_pt") ) rebin = 10;
      if( temp.Contains("_jet_") && temp.Contains("_eta") ) rebin = 4;
      if( temp.Contains("_jet_") && temp.Contains("_phi") ) rebin = 4;
      if( temp.Contains("_csv") ) rebin = 8;

      if( temp.Contains("_diele_mass_closestZmass") ) rebin = 5;


      TH1D* h_data = (TH1D*)file[file_data]->Get(temp)->Clone(temp+"_data");


      TH1D* h_mc_tt = (TH1D*)file[file_tt]->Get(temp)->Clone(temp+"_mc_tt");
      TH1D* h_mc_dy = (TH1D*)file[file_dy]->Get(temp)->Clone(temp+"_mc_dy");

      h_data->Rebin(rebin);
      h_mc_tt->Rebin(rebin);
      h_mc_dy->Rebin(rebin);

      h_mc_tt->SetFillColor(color[file_tt]);
      h_mc_dy->SetFillColor(color[file_dy]);


      h_mc_tt->Scale(sfMC_);
      h_mc_dy->Scale(sfMC_);

      TH1D* h_mc = (TH1D*)h_mc_tt->Clone(temp+"_mc");
      h_mc->Add(h_mc_dy);


      if( temp.Contains("_numPV") || renormMC_ ){
	h_mc_dy->Scale( h_data->Integral() / h_mc->Integral() );
	h_mc_tt->Scale( h_data->Integral() / h_mc->Integral() );
      }

      THStack *hs = new THStack("hs","");
      hs->Add(h_mc_tt);
      hs->Add(h_mc_dy);


      h_data->SetLineColor(color[0]);

      h_data->SetMarkerColor(color[0]);

      h_data->SetMarkerStyle(20);

      legend->AddEntry(h_data,histLabels[0],"pe1");
      legend->AddEntry(h_mc_dy,histLabels[file_dy],"f");
      legend->AddEntry(h_mc_tt,histLabels[file_tt],"f");

      //c3->SetTopMargin(.05);
      //c1->SetRightMargin(.05);

      double ratioMax = 1.6;
      double ratioMin = 0.5;

      int nbins = h_data->GetNbinsX();

      double xmin = h_data->GetBinLowEdge(1);
      double xmax = h_data->GetBinLowEdge(nbins) + h_data->GetBinWidth(nbins);

      TH1D* myRatio = new TH1D("ratio", "", nbins, xmin, xmax );

      myRatio->SetStats(0);
      myRatio->Sumw2();
      myRatio->SetLineColor(kBlack);
      myRatio->SetMarkerColor(kBlack);
      myRatio->Divide(h_data,h_mc);

      myRatio->SetMinimum(ratioMin);
      myRatio->SetMaximum(ratioMax);
      //myRatio->GetYaxis()->SetNdivisions(50000+404);
      myRatio->GetYaxis()->SetNdivisions(50000+204);
      myRatio->GetYaxis()->SetLabelSize(0.1); //make y label bigger
      myRatio->GetXaxis()->SetLabelSize(0.1); //make y label bigger
      myRatio->GetXaxis()->SetTitleOffset(1.1);
      myRatio->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle()); //make y label bigger
      myRatio->GetXaxis()->SetLabelSize(0.12);
      myRatio->GetXaxis()->SetLabelOffset(0.04);
      myRatio->GetXaxis()->SetTitleSize(0.12);
      myRatio->GetYaxis()->SetTitle("Data/MC");
      myRatio->GetYaxis()->SetTitleSize(0.09);
      myRatio->GetYaxis()->SetTitleOffset(.55);
      myC1->cd(2);
      gPad->SetTopMargin(small);
      gPad->SetTickx();
      gPad->Modified();

      myRatio->GetYaxis()->CenterTitle(kTRUE);


      if( temp.Contains("_selection") ){
	for( int iBin=0; iBin<nbins; iBin++ ) myRatio->GetXaxis()->SetBinLabel(iBin+1,h_data->GetXaxis()->GetBinLabel(iBin+1));
	myRatio->GetXaxis()->SetBinLabel(2,"HLT Ele27WPL");
	myRatio->GetXaxis()->SetBinLabel(8,"HLT HT200");
	myC1->GetPad(2)->SetBottomMargin(.4);
	myC1->GetPad(1)->SetRightMargin(.10);
	myC1->GetPad(2)->SetRightMargin(.10);
	myRatio->GetXaxis()->SetTitle("");
      }
      else{
	myC1->GetPad(2)->SetBottomMargin(.3);
	myC1->GetPad(1)->SetRightMargin(.05);
	myC1->GetPad(2)->SetRightMargin(.05);
      }

      h_data->SetStats(0);

      h_data->GetYaxis()->SetTitleOffset(1.0);
      h_data->GetYaxis()->SetTitleSize(0.05);

      h_data->GetYaxis()->SetTitle("Number of Events");


      int max_bin_data = h_data->GetMaximumBin();
      double max_data = h_data->GetBinContent(max_bin_data) + h_data->GetBinError(max_bin_data);

      int max_bin_mc = h_mc->GetMaximumBin();
      double max_mc = h_mc->GetBinContent(max_bin_mc) + h_mc->GetBinError(max_bin_mc);

      double max_content = std::max(max_data, max_mc);

      h_data->GetYaxis()->SetRangeUser(0.,1.2 * max_content);

      if( temp.Contains("_mass") ){
	h_data->GetXaxis()->SetRangeUser(40.,140.);
	myRatio->GetXaxis()->SetRangeUser(40.,140.);
      }

      if( temp.Contains("_met") ){
	h_data->GetXaxis()->SetRangeUser(0.,150.);
	myRatio->GetXaxis()->SetRangeUser(0.,150.);
      }
      // if( temp.Contains("_eta") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);
      // if( temp.Contains("_numGenPVs") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);
      // if( temp.Contains("_numJets") ) h_all->GetYaxis()->SetRangeUser(0.4,1.1);


      TLine* myLine;
      if( temp.Contains("_mass") )     myLine = new TLine(40, 1, 140, 1);
      else if( temp.Contains("_met") ) myLine = new TLine(0, 1, 150, 1);
      else                             myLine = new TLine(h_data->GetXaxis()->GetXmin(), 1, h_data->GetXaxis()->GetXmax(), 1);


      myC1->cd(1);
      h_data->Draw("pe1");
      hs->Draw("histsame");
      h_data->Draw("pe1same");
      legend->Draw();
      LumiInfoLatex.Draw();
      CMSInfoLatex.Draw();
      PublishInfoLatex.Draw();
      myC1->GetPad(1)->RedrawAxis();

      myC1->cd(2);
      myRatio->SetLineWidth(2);
      myRatio->Draw("pe1");
      myLine->Draw();
      myC1->GetPad(2)->RedrawAxis();

      plotname = dirprefix + temp_mh + "_data2mc_lin.png";
      myC1->Print(plotname);

      plotname = dirprefix + temp_mh + "_data2mc_lin.pdf";
      if( printPDF_ ) myC1->Print(plotname);


      // log
      h_data->GetYaxis()->SetRangeUser(0.4,12 * max_content);
      if( temp.Contains("_selection") ){
	h_data->GetYaxis()->SetRangeUser(1000,12 * max_content);
      }

      myC1->cd(1);
      gPad->SetLogy(1);
      myC1->GetPad(1)->RedrawAxis();

      plotname = dirprefix + temp_mh + "_data2mc_log.png";
      myC1->Print(plotname);

      plotname = dirprefix + temp_mh + "_data2mc_log.pdf";
      if( printPDF_ ) myC1->Print(plotname);

      gPad->SetLogy(0);



      delete myRatio;
      delete myLine;
      delete legend;
    } // end loop on hists
  }// end if plot1dDist_


  if( plot2dSF_ ){
    TString name_numerator = "h_probe_pt_eta_pHLTall";
    TString name_denominator = "h_probe_pt_eta_all";

    TH2D* h_data_num = (TH2D*)file[file_data]->Get(name_numerator)->Clone(name_numerator+"_data");
    TH2D* h_data_den = (TH2D*)file[file_data]->Get(name_denominator)->Clone(name_denominator+"_data");

    TH2D* h_mc_num = (TH2D*)file[file_dy]->Get(name_numerator)->Clone(name_numerator+"_mc");
    TH2D* h_mc_den = (TH2D*)file[file_dy]->Get(name_denominator)->Clone(name_denominator+"_mc");

    int NbinsX = h_mc_num->GetNbinsX();
    int NbinsY = h_mc_num->GetNbinsY();

    double eps = 0.001;

    for( int iBinX=0; iBinX<NbinsX; iBinX++ ){
      for( int iBinY=0; iBinY<NbinsY; iBinY++ ){
	double content1 = h_mc_num->GetBinContent(iBinX+1,iBinY+1);
	double content2 = h_mc_den->GetBinContent(iBinX+1,iBinY+1);

	double diff = content1 - content2;
	if( diff>0 ) content2 += diff + eps;

	h_mc_den->SetBinContent(iBinX+1,iBinY+1,content2);
      }
    }

    TH2D* h_ratio_2d = divide_ratio_plots_2d(h_data_num,h_data_den,h_mc_num,h_mc_den);

    TCanvas* c1 = new TCanvas("c2","",900,600);

    h_ratio_2d->SetStats(0);
    h_ratio_2d->SetMaximum(1.0);

    h_ratio_2d->Draw("colz");


    plotname = dirprefix + "eff_probe_pt_eta_data2mc.png";
    c1->Print(plotname);

    plotname = dirprefix + "eff_probe_pt_eta_data2mc.pdf";
    if( printPDF_ ) c1->Print(plotname);


    h_ratio_2d->GetXaxis()->SetRangeUser(30.001,100-0.0001);
    h_ratio_2d->Draw("colztexte1");

    gStyle->SetPaintTextFormat( "4.3f" );
    plotname = dirprefix + "eff_probe_pt_eta_data2mc_lim.png";
    c1->Print(plotname);

    plotname = dirprefix + "eff_probe_pt_eta_data2mc_lim.pdf";
    if( printPDF_ ) c1->Print(plotname);

  }

  // Close the file
  file[file_data]->Close();
  file[file_dy]->Close();
  file[file_tt]->Close();
}


TH1D* divide_ratio_plots( TH1D* h_numerator_numerator, TH1D* h_numerator_denominator, TH1D* h_denominator_numerator, TH1D* h_denominator_denominator ){
 
  TString title11 = h_numerator_numerator->GetName();
  TString title12 = h_numerator_denominator->GetName();
  TString title21 = h_denominator_numerator->GetName();
  TString title22 = h_denominator_denominator->GetName();

  TString title = title11 + "_" + title12 + "_" + title21 + "_" + title22;

  TH1D* h_ratio = (TH1D*)h_numerator_numerator->Clone(title+"_numerator");
  TH1D* h_ratio_denominator = (TH1D*)h_numerator_numerator->Clone(title+"_denominator");

  h_ratio->Reset();
  h_ratio_denominator->Reset();

  h_ratio->Divide( h_numerator_numerator, h_numerator_denominator, 1, 1, "B" );
  h_ratio_denominator->Divide( h_denominator_numerator, h_denominator_denominator, 1, 1, "B" );

  // h_ratio->Divide( h_numerator_denominator );
  //h_ratio_denominator->Divide( h_denominator_denominator );

  h_ratio->Divide( h_ratio_denominator );

  return h_ratio;
  //delete h_ratio_denominator;
}


TH2D* divide_ratio_plots_2d( TH2D* h_numerator_numerator, TH2D* h_numerator_denominator, TH2D* h_denominator_numerator, TH2D* h_denominator_denominator ){
 
  TString title11 = h_numerator_numerator->GetName();
  TString title12 = h_numerator_denominator->GetName();
  TString title21 = h_denominator_numerator->GetName();
  TString title22 = h_denominator_denominator->GetName();

  TString title = title11 + "_" + title12 + "_" + title21 + "_" + title22;

  TEfficiency* eff_numerator = new TEfficiency(*h_numerator_numerator,*h_numerator_denominator);
  TEfficiency* eff_denominator = new TEfficiency(*h_denominator_numerator,*h_denominator_denominator);

  TH2D* h_ratio = (TH2D*)h_numerator_numerator->Clone(title+"_numerator");
  TH2D* h_ratio_denominator = (TH2D*)h_numerator_numerator->Clone(title+"_denominator");

  h_ratio->Reset();
  h_ratio_denominator->Reset();

  h_ratio->Divide( h_numerator_numerator, h_numerator_denominator, 1, 1, "B" );
  h_ratio_denominator->Divide( h_denominator_numerator, h_denominator_denominator, 1, 1, "B" );

  // h_ratio->Divide( h_numerator_denominator );
  //h_ratio_denominator->Divide( h_denominator_denominator );

  h_ratio->Divide( h_ratio_denominator );


  int NbinsX = h_ratio->GetNbinsX();
  int NbinsY = h_ratio->GetNbinsY();

  for( int iBinX=0; iBinX<NbinsX; iBinX++ ){
    for( int iBinY=0; iBinY<NbinsY; iBinY++ ){

      int useBin = eff_numerator->GetGlobalBin(iBinX+1,iBinY+1);

      double numerator_eff = eff_numerator->GetEfficiency(useBin);
      double numerator_eff_up = eff_numerator->GetEfficiencyErrorUp(useBin);
      double numerator_eff_down = eff_numerator->GetEfficiencyErrorLow(useBin);

      double denominator_eff = eff_denominator->GetEfficiency(useBin);
      double denominator_eff_up = eff_denominator->GetEfficiencyErrorUp(useBin);
      double denominator_eff_down = eff_denominator->GetEfficiencyErrorLow(useBin);

      double sf = (denominator_eff>0) ? numerator_eff/denominator_eff : 0.;

      double err_num_eff_up = (numerator_eff>0) ? numerator_eff_up/numerator_eff : 0.;
      double err_den_eff_up = (denominator_eff>0) ? denominator_eff_up/denominator_eff : 0.;

      double err_num_eff_down = (numerator_eff>0) ? numerator_eff_down/numerator_eff : 0.;
      double err_den_eff_down = (denominator_eff>0) ? denominator_eff_down/denominator_eff : 0.;

      double sf_up = sf * sqrt( (err_num_eff_up*err_num_eff_up) + (err_den_eff_up*err_den_eff_up) );
      double sf_down = sf * sqrt( (err_num_eff_down*err_num_eff_down) + (err_den_eff_down*err_den_eff_down) );

      double sf_err = 0.5*(sf_up + sf_down);

      // Additional systematic uncertainty due to bad online BS
      //https://hypernews.cern.ch/HyperNews/CMS/get/susy/2115.html
      sf_err += 0.02;

      double xmin = h_numerator_numerator->GetXaxis()->GetBinLowEdge(iBinX+1);
      double xmax = h_numerator_numerator->GetXaxis()->GetBinUpEdge(iBinX+1);

      double ymin = h_numerator_numerator->GetYaxis()->GetBinLowEdge(iBinY+1);
      double ymax = h_numerator_numerator->GetYaxis()->GetBinUpEdge(iBinY+1);

      double ratio_content = h_ratio->GetBinContent(iBinX+1,iBinY+1);
      double ratio_err = h_ratio->GetBinError(iBinX+1,iBinY+1);

      h_ratio->SetBinError(iBinX+1,iBinY+1,sf_err);

      if( sf>0. ){
	printf(" useBin = %3d, bX = %2d, bY = %2d, xmin = %3.0f, xmax = %3.0f, ymin = %+.2f, ymax = %+.2f: ratio = %.3f, ratio_err = %.3f, sf = %.3f, sf_u = %.3f, sf_d = %.3f, sf_err = %.3f\n", useBin, iBinX, iBinY, xmin, xmax, ymin, ymax, ratio_content, ratio_err, sf, sf_up, sf_down, sf_err);

	// printf(" useBin = %3d, bX = %2d, bY = %2d, xmin = %.0f, xmax = %.0f, ymin = %.2f, ymax = %.2f: ratio = %.3f, sf = %.3f, sf_u = %.3f, sf_d = %.3f, num_eff = %.3f, den_eff = %.3f, num_eff_u = %.3f, num_eff_d = %.3f, den_eff_u = %.3f, den_eff_d = %.3f \n", useBin, iBinX, iBinY, xmin, xmax, ymin, ymax, ratio_content, sf, sf_up, sf_down, numerator_eff, denominator_eff, numerator_eff_up, numerator_eff_down, denominator_eff_up, denominator_eff_down);
      }
    }
  }

  return h_ratio;
  //delete h_ratio_denominator;
}
