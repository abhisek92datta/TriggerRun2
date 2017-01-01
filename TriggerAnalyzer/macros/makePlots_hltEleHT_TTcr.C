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

TEfficiency* makeEfficiency( TH1D* h_numerator, TH1D* h_denominator );

//*****************************************************************************


void makePlots_hltEleHT_TTcr( bool printPDF_ = false, int useSample_ = 0 ){

  TH1::SetDefaultSumw2();

  int NumSamples = 3;
  TFile* file[NumSamples];
  file[0] = new TFile("HistoFiles/hltEleHT_treeReader_TTcr_TT_13TeV_Spring15_Asympt25ns_histo.root");
  file[1] = new TFile("HistoFiles/hltEleHT_treeReader_TTcr_ttHTobb_M125_13TeV_powheg_pythia8_Spring15_Asympt25ns_histo.root");
  file[2] = new TFile("HistoFiles/hltEleHT_treeReader_TTcr_SingleElectron_Run2015D_PromptReco_254231_258158_histo.root");

  std::vector<TString> histLabels(NumSamples);
  histLabels[0] = "TTJets";
  histLabels[1] = "ttHTobb";
  histLabels[2] = "Data";

  Color_t color[5];
  color[0] = kBlack;
  color[1] = kBlue;
  color[2] = kRed;
  color[3] = kGreen+1;

  // color[2] = kBlack;
  // color[3] = kGreen+1;
  // color[4] = kMagenta+2;
  // color[5] = kRed+1;
  // color[6] = kGreen-5;
  // color[7] = kRed+3;

  std::vector<std::vector<int> > projection_start_bins;
  std::vector<std::vector<int> > projection_end_bins;
  std::vector<std::vector<TString> > projection_labels;

  // elePt
  std::vector<int> proj_elePt_start_bins;
  proj_elePt_start_bins.push_back(1);
  proj_elePt_start_bins.push_back(4);
  proj_elePt_start_bins.push_back(6);
  proj_elePt_start_bins.push_back(7);

  std::vector<int> proj_elePt_end_bins;
  proj_elePt_end_bins.push_back(12);
  proj_elePt_end_bins.push_back(5);
  proj_elePt_end_bins.push_back(8);
  proj_elePt_end_bins.push_back(12);

  std::vector<TString> proj_elePt_labels;
  proj_elePt_labels.push_back("elePt30toInf");
  proj_elePt_labels.push_back("elePt30to50");
  proj_elePt_labels.push_back("elePt50to80");
  proj_elePt_labels.push_back("elePt80toInf");

  // numJet
  std::vector<int> proj_numJet_start_bins;
  proj_numJet_start_bins.push_back(1);
  proj_numJet_start_bins.push_back(1);
  proj_numJet_start_bins.push_back(4);
  proj_numJet_start_bins.push_back(5);

  std::vector<int> proj_numJet_end_bins;
  proj_numJet_end_bins.push_back(8);
  proj_numJet_end_bins.push_back(3);
  proj_numJet_end_bins.push_back(4);
  proj_numJet_end_bins.push_back(8);

  std::vector<TString> proj_numJet_labels;
  proj_numJet_labels.push_back("numJet0toInf");
  proj_numJet_labels.push_back("numJet0to2");
  proj_numJet_labels.push_back("numJet3to3");
  proj_numJet_labels.push_back("numJet4toInf");


  projection_start_bins.push_back(proj_elePt_start_bins);
  projection_end_bins.push_back(proj_elePt_end_bins);
  projection_labels.push_back(proj_elePt_labels);

  projection_start_bins.push_back(proj_numJet_start_bins);
  projection_end_bins.push_back(proj_numJet_end_bins);
  projection_labels.push_back(proj_numJet_labels);

  projection_start_bins.push_back(proj_elePt_start_bins);
  projection_end_bins.push_back(proj_elePt_end_bins);
  projection_labels.push_back(proj_elePt_labels);

  projection_start_bins.push_back(proj_elePt_start_bins);
  projection_end_bins.push_back(proj_elePt_end_bins);
  projection_labels.push_back(proj_elePt_labels);

  std::vector<int> NumProjBins;
  NumProjBins.push_back( int(proj_elePt_labels.size()) );
  NumProjBins.push_back( int(proj_numJet_labels.size()) );
  NumProjBins.push_back( int(proj_elePt_labels.size()) );
  NumProjBins.push_back( int(proj_elePt_labels.size()) );



  TString dirprefix = "Images/Images_2015_10_13_hltEleHT_TTcr_" + histLabels[useSample_] + "/";

  struct stat st;
  if( stat(dirprefix.Data(),&st) != 0 )  mkdir(dirprefix.Data(),0777);


 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<std::string> histoname1;
  std::vector<std::string> histoname2;
  std::vector<std::string> histoname3;

  histoname1.push_back("h_HT30");
  histoname1.push_back("h_HT30er");
  histoname1.push_back("h_HT30_4j");

  histoname2.push_back("elePt");
  histoname2.push_back("numJet");
  histoname2.push_back("eleEBPt");
  histoname2.push_back("eleEEPt");


  histoname3.push_back("h_event_selection");


 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TGaxis::SetMaxDigits(3);

  TString lumiinfo = "553 pb^{-1} (13 TeV)";
  TLatex LumiInfoLatex(0.65, 0.94, lumiinfo);
  LumiInfoLatex.SetNDC(); LumiInfoLatex.SetTextFont(42);
  LumiInfoLatex.SetTextSize(0.04);

  //TString cmsinfo =   "CMS Preliminary";
  TString cmsinfo =   "CMS";
  TLatex CMSInfoLatex(0.13, 0.94, cmsinfo);
  CMSInfoLatex.SetNDC(); CMSInfoLatex.SetTextFont(42);
  CMSInfoLatex.SetTextFont(61);
  CMSInfoLatex.SetTextSize(0.055); //SBOUTLE

  std::string publishinfo =   "Preliminary"; //DPUIGH
  TLatex PublishInfoLatex(0.26, 0.94, publishinfo.c_str()); //SBOUTLE
  PublishInfoLatex.SetNDC();
  PublishInfoLatex.SetTextFont(52);
  PublishInfoLatex.SetTextSize(0.045); //SBOUTLE


  TString plotname;

  TCanvas* c1 = new TCanvas("c1", "c1", 600,700);


  TH2D* h_L1HTT_elePt = (TH2D*)file[useSample_]->Get("h_L1HTT_elePt");
  TH2D* h_HT30_HT30er = (TH2D*)file[useSample_]->Get("h_HT30_HT30er");
  TH2D* h_HT30_L1HTT = (TH2D*)file[useSample_]->Get("h_HT30_L1HTT");

  TProfile* p_L1HTT_elePt = (TProfile*)h_L1HTT_elePt->ProfileX("p_L1HTT_elePt");
  // TProfile* p_HT30_HT30er = (TProfile*)h_HT30_HT30er->ProfileX("p_HT30_HT30er");
  // TProfile* p_HT30_L1HTT = (TProfile*)h_HT30_L1HTT->ProfileX("p_HT30_L1HTT");

  p_L1HTT_elePt->SetMarkerStyle(20);

  h_L1HTT_elePt->Draw("colz");
  p_L1HTT_elePt->Draw("pe1same");
  plotname = dirprefix + "h_L1HTT_elePt" + "_2D_colz.png";
  c1->Print(plotname);

  h_HT30_HT30er->Draw("colz");
  //p_HT30_HT30er->Draw("pe1same");
  plotname = dirprefix + "h_HT30_HT30er" + "_2D_colz.png";
  c1->Print(plotname);

  h_HT30_L1HTT->Draw("colz");
  //p_HT30_L1HTT->Draw("pe1same");
  plotname = dirprefix + "h_HT30_L1HTT" + "_2D_colz.png";
  c1->Print(plotname);


  for( int i=0; i<int(histoname1.size()); i++ ){

    for( int j=0; j<int(histoname2.size()); j++ ){

      for( int k=0; k<2; k++ ){

	TString temp = histoname1[i];

	TString suffix = histoname2[j];

	TString l1suffix = ( k==0 ) ? "125" : "100";

	if( temp!="h_HT30" && suffix!="elePt" ) continue;

	TString temp_L1 = temp + "_L1HTT" + l1suffix + "_" + suffix;
	TString temp_HLT = temp + "_L1HTT" + l1suffix + "_passHLTEle27HT200_" + suffix;

	TString temp_mh = temp;
	temp_mh.ReplaceAll("h_","");


	temp = temp + "_" + suffix;

	//TLegend *legend = new TLegend(0.2,0.83,0.9,0.89);
	TLegend *legend = new TLegend(0.2,0.85,0.9,0.91);

	legend->SetFillColor(kWhite);
	legend->SetLineColor(kWhite);
	legend->SetShadowColor(kWhite);
	legend->SetTextFont(42);
	legend->SetTextSize(0.04);

	legend->SetNColumns(2);

	int rebin = ( useSample_==2 ) ? 25 : 10;


	TH2D* h_all = (TH2D*)file[useSample_]->Get(temp)->Clone(temp+"_"+suffix+"_"+l1suffix);
	TH2D* h_l1t = (TH2D*)file[useSample_]->Get(temp_L1)->Clone(temp_L1+"_"+suffix+"_"+l1suffix);
	TH2D* h_hlt = (TH2D*)file[useSample_]->Get(temp_HLT)->Clone(temp_HLT+"_"+suffix+"_"+l1suffix);


	TH2D* h_temp_hlt_all = (TH2D*)h_all->Clone("h_temp_hlt_all_"+suffix+"_"+l1suffix);
	TH2D* h_temp_l1t_all = (TH2D*)h_all->Clone("h_temp_l1t_all_"+suffix+"_"+l1suffix);
	TH2D* h_temp_hlt_l1t = (TH2D*)h_l1t->Clone("h_temp_hlt_l1t_"+suffix+"_"+l1suffix);

	TH2D* h_ratio_hlt_all = (TH2D*)h_hlt->Clone("h_ratio_hlt_all_"+suffix+"_"+l1suffix);
	TH2D* h_ratio_l1t_all = (TH2D*)h_l1t->Clone("h_ratio_l1t_all_"+suffix+"_"+l1suffix);
	TH2D* h_ratio_hlt_l1t = (TH2D*)h_hlt->Clone("h_ratio_hlt_l1t_"+suffix+"_"+l1suffix);

	h_temp_hlt_all->RebinY(rebin);
	h_temp_l1t_all->RebinY(rebin);
	h_temp_hlt_l1t->RebinY(rebin);

	h_ratio_hlt_all->RebinY(rebin);
	h_ratio_l1t_all->RebinY(rebin);
	h_ratio_hlt_l1t->RebinY(rebin);

	h_ratio_hlt_all->Divide(h_temp_hlt_all);
	h_ratio_l1t_all->Divide(h_temp_l1t_all);
	h_ratio_hlt_l1t->Divide(h_temp_hlt_l1t);


	TProfile* p_all = (TProfile*)h_all->ProfileX("p_all");
	TProfile* p_l1t = (TProfile*)h_l1t->ProfileX("p_l1t");
	TProfile* p_hlt = (TProfile*)h_hlt->ProfileX("p_hlt");

	p_all->SetMarkerStyle(20);
	p_l1t->SetMarkerStyle(20);
	p_hlt->SetMarkerStyle(20);

	h_all->Draw("colz");
	p_all->Draw("pe1same");
	plotname = dirprefix + temp + "_2D_colz.png";
	c1->Print(plotname);

	h_l1t->Draw("colz");
	p_l1t->Draw("pe1same");
	plotname = dirprefix + temp_L1 + "_2D_colz.png";
	c1->Print(plotname);

	h_hlt->Draw("colz");
	p_hlt->Draw("pe1same");
	plotname = dirprefix + temp_HLT + "_2D_colz.png";
	c1->Print(plotname);



	h_ratio_hlt_all->SetStats(0);
	h_ratio_hlt_all->GetYaxis()->SetRangeUser(0.,400.);
	h_ratio_hlt_all->Draw("colz");
	plotname = dirprefix + "h_ratio_hlt_all" + "_2D_colz.png";
	c1->Print(plotname);

	h_ratio_l1t_all->SetStats(0);
	h_ratio_l1t_all->GetYaxis()->SetRangeUser(0.,400.);
	h_ratio_l1t_all->Draw("colz");
	plotname = dirprefix + "h_ratio_l1t_all" + "_2D_colz.png";
	c1->Print(plotname);

	h_ratio_hlt_l1t->SetStats(0);
	h_ratio_hlt_l1t->GetYaxis()->SetRangeUser(0.,400.);
	h_ratio_hlt_l1t->Draw("colz");
	plotname = dirprefix + "h_ratio_hlt_l1t" + "_2D_colz.png";
	c1->Print(plotname);

	int numProjBin = NumProjBins[j];

	TH1D* h_py_all[numProjBin];
	TH1D* h_py_l1t[numProjBin];
	TH1D* h_py_hlt[numProjBin];

	TEfficiency* eff_hlt_all[numProjBin];
	TEfficiency* eff_l1t_all[numProjBin];
	TEfficiency* eff_hlt_l1t[numProjBin];

	for( int iBin=0; iBin<numProjBin; iBin++ ){

	  int bin_start = projection_start_bins[j][iBin];
	  int bin_end   = projection_end_bins[j][iBin];
	  h_py_all[iBin] = (TH1D*)h_all->ProjectionY(Form("h_py_all_%d",iBin),bin_start,bin_end);
	  h_py_l1t[iBin] = (TH1D*)h_l1t->ProjectionY(Form("h_py_l1t_%d",iBin),bin_start,bin_end);
	  h_py_hlt[iBin] = (TH1D*)h_hlt->ProjectionY(Form("h_py_hlt_%d",iBin),bin_start,bin_end);

	  h_py_all[iBin]->Rebin(rebin);
	  h_py_l1t[iBin]->Rebin(rebin);
	  h_py_hlt[iBin]->Rebin(rebin);

	  eff_l1t_all[iBin] = makeEfficiency(h_py_l1t[iBin], h_py_all[iBin]);
	  eff_hlt_l1t[iBin] = makeEfficiency(h_py_hlt[iBin], h_py_l1t[iBin]);
	  eff_hlt_all[iBin] = makeEfficiency(h_py_hlt[iBin], h_py_all[iBin]);

	  eff_l1t_all[iBin]->SetLineColor(color[iBin]);
	  eff_l1t_all[iBin]->SetMarkerColor(color[iBin]);
	  eff_l1t_all[iBin]->SetMarkerStyle(20);

	  eff_hlt_l1t[iBin]->SetLineColor(color[iBin]);
	  eff_hlt_l1t[iBin]->SetMarkerColor(color[iBin]);
	  eff_hlt_l1t[iBin]->SetMarkerStyle(20);

	  eff_hlt_all[iBin]->SetLineColor(color[iBin]);
	  eff_hlt_all[iBin]->SetMarkerColor(color[iBin]);
	  eff_hlt_all[iBin]->SetMarkerStyle(20);

	  legend->AddEntry(eff_l1t_all[iBin], projection_labels[j][iBin],"pe1");
	}



	h_py_all[0]->SetStats(0);
	h_py_all[0]->GetYaxis()->SetTitle("Efficiency");

	h_py_all[0]->GetYaxis()->SetRangeUser(0.,1.15);
	h_py_all[0]->GetXaxis()->SetRangeUser(0.,600.);
	// h_data_all->SetStats(0);

	c1->SetTopMargin(.07);
	c1->SetRightMargin(.05);

	// h_data_all->GetYaxis()->SetTitleOffset(1.0);
	// h_data_all->GetYaxis()->SetTitleSize(0.05);

	// h_data_all->GetYaxis()->SetRangeUser(0.,1.15);
	// if( temp.Contains("_pt") ){
	//   h_data_all->GetXaxis()->SetRangeUser(0.,150.);
	//   myRatio->GetXaxis()->SetRangeUser(0.,150.);
	// }


	//// l1t_all
	h_py_all[0]->Draw("axis");
	for( int iBin=0; iBin<numProjBin; iBin++ ) eff_l1t_all[iBin]->Draw("pe1same");
	legend->Draw();
	LumiInfoLatex.Draw();
	CMSInfoLatex.Draw();
	PublishInfoLatex.Draw();

	plotname = dirprefix + temp_mh + "_"+suffix+"_l1t" + l1suffix + "_all_lin.png";
	c1->Print(plotname);

	plotname = dirprefix + temp_mh + "_"+suffix+"_l1t" + l1suffix + "_all_lin.pdf";
	if( printPDF_ ) c1->Print(plotname);


	//// hlt_all
	h_py_all[0]->Draw("axis");
	for( int iBin=0; iBin<numProjBin; iBin++ ) eff_hlt_all[iBin]->Draw("pe1same");
	legend->Draw();
	LumiInfoLatex.Draw();
	CMSInfoLatex.Draw();
	PublishInfoLatex.Draw();

	plotname = dirprefix + temp_mh + "_"+suffix+"_hlt_l1t" + l1suffix + "_all_lin.png";
	c1->Print(plotname);

	plotname = dirprefix + temp_mh + "_"+suffix+"_hlt_l1t" + l1suffix + "_all_lin.pdf";
	if( printPDF_ ) c1->Print(plotname);


	//// hlt_l1t
	h_py_all[0]->Draw("axis");
	for( int iBin=0; iBin<numProjBin; iBin++ ) eff_hlt_l1t[iBin]->Draw("pe1same");
	legend->Draw();
	LumiInfoLatex.Draw();
	CMSInfoLatex.Draw();
	PublishInfoLatex.Draw();

	plotname = dirprefix + temp_mh + "_"+suffix+"_hlt_l1t" + l1suffix + "_lin.png";
	c1->Print(plotname);

	plotname = dirprefix + temp_mh + "_"+suffix+"_hlt_l1t" + l1suffix + "_lin.pdf";
	if( printPDF_ ) c1->Print(plotname);




	h_py_all[0]->GetYaxis()->SetRangeUser(0.85,1.15);



	//// l1t_all
	h_py_all[0]->Draw("axis");
	for( int iBin=0; iBin<numProjBin; iBin++ ) eff_l1t_all[iBin]->Draw("pe1same");
	legend->Draw();
	LumiInfoLatex.Draw();
	CMSInfoLatex.Draw();
	PublishInfoLatex.Draw();

	plotname = dirprefix + temp_mh + "_"+suffix+"_l1t" + l1suffix + "_all_lin_zoom.png";
	c1->Print(plotname);

	plotname = dirprefix + temp_mh + "_"+suffix+"_l1t" + l1suffix + "_all_lin_zoom.pdf";
	if( printPDF_ ) c1->Print(plotname);


	//// hlt_all
	h_py_all[0]->Draw("axis");
	for( int iBin=0; iBin<numProjBin; iBin++ ) eff_hlt_all[iBin]->Draw("pe1same");
	legend->Draw();
	LumiInfoLatex.Draw();
	CMSInfoLatex.Draw();
	PublishInfoLatex.Draw();

	plotname = dirprefix + temp_mh + "_"+suffix+"_hlt_l1t" + l1suffix + "_all_lin_zoom.png";
	c1->Print(plotname);

	plotname = dirprefix + temp_mh + "_"+suffix+"_hlt_l1t" + l1suffix + "_all_lin_zoom.png";
	if( printPDF_ ) c1->Print(plotname);


	//// hlt_l1t
	h_py_all[0]->Draw("axis");
	for( int iBin=0; iBin<numProjBin; iBin++ ) eff_hlt_l1t[iBin]->Draw("pe1same");
	legend->Draw();
	LumiInfoLatex.Draw();
	CMSInfoLatex.Draw();
	PublishInfoLatex.Draw();

	plotname = dirprefix + temp_mh + "_"+suffix+"_hlt_l1t" + l1suffix + "_lin_zoom.png";
	c1->Print(plotname);

	plotname = dirprefix + temp_mh + "_"+suffix+"_hlt_l1t" + l1suffix + "_lin_zoom.pdf";
	if( printPDF_ ) c1->Print(plotname);


	delete legend;
      } // end loop on hists
    }
  }


  if( true ){

    TString temp = "h_HT30";

    TString temp_L1 = temp + "_L1HTT125";
    TString temp_HLT = temp + "_L1HTT125_passHLTEle27HT200";

    TString temp_mh = temp;
    temp_mh.ReplaceAll("h_","");

    //TLegend *legend = new TLegend(0.2,0.83,0.9,0.89);
    TLegend *legend = new TLegend(0.2,0.87,0.85,0.91);

    legend->SetFillColor(kWhite);
    legend->SetLineColor(kWhite);
    legend->SetShadowColor(kWhite);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);

    legend->SetNColumns(3);

    int rebin = 25;

    TH1D* h_all[NumSamples];
    TH1D* h_l1t[NumSamples];
    TH1D* h_hlt[NumSamples];

    TEfficiency* eff_hlt_all[NumSamples];
    TEfficiency* eff_l1t_all[NumSamples];
    TEfficiency* eff_hlt_l1t[NumSamples];

    for( int iSample=0; iSample<NumSamples; iSample++ ){
      h_all[iSample] = (TH1D*)file[iSample]->Get(temp.Data())->Clone(Form("%s_%s",temp.Data(),histLabels[iSample].Data()));
      h_l1t[iSample] = (TH1D*)file[iSample]->Get(temp_L1.Data())->Clone(Form("%s_%s",temp.Data(),histLabels[iSample].Data()));
      h_hlt[iSample] = (TH1D*)file[iSample]->Get(temp_HLT.Data())->Clone(Form("%s_%s",temp.Data(),histLabels[iSample].Data()));

      h_all[iSample]->Rebin(rebin);
      h_l1t[iSample]->Rebin(rebin);
      h_hlt[iSample]->Rebin(rebin);

      eff_l1t_all[iSample] = makeEfficiency(h_l1t[iSample], h_all[iSample]);
      eff_hlt_l1t[iSample] = makeEfficiency(h_hlt[iSample], h_l1t[iSample]);
      eff_hlt_all[iSample] = makeEfficiency(h_hlt[iSample], h_all[iSample]);

      eff_l1t_all[iSample]->SetLineColor(color[iSample]);
      eff_l1t_all[iSample]->SetMarkerColor(color[iSample]);
      eff_l1t_all[iSample]->SetMarkerStyle(20);

      eff_hlt_l1t[iSample]->SetLineColor(color[iSample]);
      eff_hlt_l1t[iSample]->SetMarkerColor(color[iSample]);
      eff_hlt_l1t[iSample]->SetMarkerStyle(20);

      eff_hlt_all[iSample]->SetLineColor(color[iSample]);
      eff_hlt_all[iSample]->SetMarkerColor(color[iSample]);
      eff_hlt_all[iSample]->SetMarkerStyle(20);

      legend->AddEntry(eff_l1t_all[iSample], histLabels[iSample],"pl");
    }



    h_all[0]->SetStats(0);
    h_all[0]->GetYaxis()->SetTitle("Efficiency");

    h_all[0]->GetYaxis()->SetRangeUser(0.,1.15);
    h_all[0]->GetXaxis()->SetRangeUser(0.,600.);


    //// l1t_all
    h_all[0]->Draw("axis");
    for( int iSample=0; iSample<NumSamples; iSample++ ) eff_l1t_all[iSample]->Draw("pe1same");
    legend->Draw();
    LumiInfoLatex.Draw();
    CMSInfoLatex.Draw();
    PublishInfoLatex.Draw();

    plotname = dirprefix + temp_mh + "_compareSamples_l1t_all_lin.png";
    c1->Print(plotname);

    plotname = dirprefix + temp_mh + "_compareSamples_l1t_all_lin.pdf";
    if( printPDF_ ) c1->Print(plotname);


    //// hlt_all
    h_all[0]->Draw("axis");
    for( int iSample=0; iSample<NumSamples; iSample++ ) eff_hlt_all[iSample]->Draw("pe1same");
    legend->Draw();
    LumiInfoLatex.Draw();
    CMSInfoLatex.Draw();
    PublishInfoLatex.Draw();

    plotname = dirprefix + temp_mh + "_compareSamples_hlt_all_lin.png";
    c1->Print(plotname);

    plotname = dirprefix + temp_mh + "_compareSamples_hlt_all_lin.pdf";
    if( printPDF_ ) c1->Print(plotname);


    //// hlt_l1t
    h_all[0]->Draw("axis");
    for( int iSample=0; iSample<NumSamples; iSample++ ) eff_hlt_l1t[iSample]->Draw("pe1same");
    legend->Draw();
    LumiInfoLatex.Draw();
    CMSInfoLatex.Draw();
    PublishInfoLatex.Draw();

    plotname = dirprefix + temp_mh + "_compareSamples_hlt_l1t_lin.png";
    c1->Print(plotname);

    plotname = dirprefix + temp_mh + "_compareSamples_hlt_l1t_lin.pdf";
    if( printPDF_ ) c1->Print(plotname);




    h_all[0]->GetYaxis()->SetRangeUser(0.85,1.15);



    //// l1t_all
    h_all[0]->Draw("axis");
    for( int iSample=0; iSample<NumSamples; iSample++ ) eff_l1t_all[iSample]->Draw("pe1same");
    legend->Draw();
    LumiInfoLatex.Draw();
    CMSInfoLatex.Draw();
    PublishInfoLatex.Draw();

    plotname = dirprefix + temp_mh + "_compareSamples_l1t_all_lin_zoom.png";
    c1->Print(plotname);

    plotname = dirprefix + temp_mh + "_compareSamples_l1t_all_lin_zoom.pdf";
    if( printPDF_ ) c1->Print(plotname);


    //// hlt_all
    h_all[0]->Draw("axis");
    for( int iSample=0; iSample<NumSamples; iSample++ ) eff_hlt_all[iSample]->Draw("pe1same");
    legend->Draw();
    LumiInfoLatex.Draw();
    CMSInfoLatex.Draw();
    PublishInfoLatex.Draw();

    plotname = dirprefix + temp_mh + "_compareSamples_hlt_all_lin_zoom.png";
    c1->Print(plotname);

    plotname = dirprefix + temp_mh + "_compareSamples_hlt_all_lin_zoom.pdf";
    if( printPDF_ ) c1->Print(plotname);


    //// hlt_l1t
    h_all[0]->Draw("axis");
    for( int iSample=0; iSample<NumSamples; iSample++ ) eff_hlt_l1t[iSample]->Draw("pe1same");
    legend->Draw();
    LumiInfoLatex.Draw();
    CMSInfoLatex.Draw();
    PublishInfoLatex.Draw();

    plotname = dirprefix + temp_mh + "_compareSamples_hlt_l1t_lin_zoom.png";
    c1->Print(plotname);

    plotname = dirprefix + temp_mh + "_compareSamples_hlt_l1t_lin_zoom.pdf";
    if( printPDF_ ) c1->Print(plotname);


    delete legend;
  }


  if( true ){

    TString temp = "h_HT30";

    std::vector<TString> useL1Name;
    useL1Name.push_back(temp + "_L1HTT125");
    useL1Name.push_back(temp + "_L1HTT100");

    std::vector<TString> labelL1Name;
    labelL1Name.push_back("L1_HTT125");
    labelL1Name.push_back("L1_HTT100");

    std::vector<Color_t> useColor;
    useColor.push_back(kBlack);
    useColor.push_back(kRed);

    TString temp_mh = temp;
    temp_mh.ReplaceAll("h_","");

    //TLegend *legend = new TLegend(0.2,0.83,0.9,0.89);
    TLegend *legend = new TLegend(0.15,0.83,0.85,0.89);

    legend->SetFillColor(kWhite);
    legend->SetLineColor(kWhite);
    legend->SetShadowColor(kWhite);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);

    legend->SetNColumns(2);

    int rebin = 25;

    int NumL1Samples = int(useL1Name.size());

    TH1D* h_all[NumL1Samples];
    TH1D* h_l1t[NumL1Samples];

    TEfficiency* eff_l1t_all[NumL1Samples];

    for( int iSample=0; iSample<NumL1Samples; iSample++ ){
      h_all[iSample] = (TH1D*)file[useSample_]->Get(temp)->Clone(Form("%s_%s",temp.Data(),labelL1Name[iSample].Data()));
      h_l1t[iSample] = (TH1D*)file[useSample_]->Get(useL1Name[iSample])->Clone(Form("%s_%s",useL1Name[iSample].Data(),labelL1Name[iSample].Data()));

      h_all[iSample]->Rebin(rebin);
      h_l1t[iSample]->Rebin(rebin);

      eff_l1t_all[iSample] = makeEfficiency(h_l1t[iSample], h_all[iSample]);

      eff_l1t_all[iSample]->SetLineColor(useColor[iSample]);
      eff_l1t_all[iSample]->SetMarkerColor(useColor[iSample]);
      eff_l1t_all[iSample]->SetMarkerStyle(20);

      legend->AddEntry(eff_l1t_all[iSample], labelL1Name[iSample],"pl");
    }



    h_all[0]->SetStats(0);
    h_all[0]->GetYaxis()->SetTitle("Efficiency");

    h_all[0]->GetYaxis()->SetRangeUser(0.,1.15);
    h_all[0]->GetXaxis()->SetRangeUser(0.,600.);



    //// l1t_all
    h_all[0]->Draw("axis");
    for( int iSample=0; iSample<NumL1Samples; iSample++ ) eff_l1t_all[iSample]->Draw("pe1same");
    legend->Draw();
    LumiInfoLatex.Draw();
    CMSInfoLatex.Draw();
    PublishInfoLatex.Draw();

    plotname = dirprefix + temp_mh + "_compareL1_l1t_all_lin.png";
    c1->Print(plotname);

    plotname = dirprefix + temp_mh + "_compareL1_l1t_all_lin.pdf";
    if( printPDF_ ) c1->Print(plotname);



    h_all[0]->GetYaxis()->SetRangeUser(0.85,1.15);



    //// l1t_all
    h_all[0]->Draw("axis");
    for( int iSample=0; iSample<NumL1Samples; iSample++ ) eff_l1t_all[iSample]->Draw("pe1same");
    legend->Draw();
    LumiInfoLatex.Draw();
    CMSInfoLatex.Draw();
    PublishInfoLatex.Draw();

    plotname = dirprefix + temp_mh + "_compareL1_l1t_all_lin_zoom.png";
    c1->Print(plotname);

    plotname = dirprefix + temp_mh + "_compareL1_l1t_all_lin_zoom.pdf";
    if( printPDF_ ) c1->Print(plotname);


    delete legend;
  }




  ////////////////////////////////////////////////////////////////////////////

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


  for( int i=0; i<int(histoname3.size()); i++ ){

    TString temp = histoname3[i];
    
    TString temp_mh = temp;
    temp_mh.ReplaceAll("h_","");

    //TLegend *legend = new TLegend(0.1,0.91,0.9,0.99);
    TLegend *legend = new TLegend(0.2,0.83,0.88,0.89);

    legend->SetFillColor(kWhite);
    legend->SetLineColor(kWhite);
    legend->SetShadowColor(kWhite);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);

    legend->SetNColumns(2);

    int rebin = 1;

    if( temp.Contains("_mass") ) rebin = 2;
    if( temp.Contains("_L1HTT") ) rebin = 2;
    if( temp.Contains("_met") ) rebin = 5;
    if( temp.Contains("_jet_") && temp.Contains("_pt") ) rebin = 10;
    if( temp.Contains("_jet_") && temp.Contains("_eta") ) rebin = 4;
    if( temp.Contains("_jet_") && temp.Contains("_phi") ) rebin = 4;
    if( temp.Contains("_csv") ) rebin = 8;

    if( temp.Contains("_diele_mass_closestZmass") ) rebin = 5;


    TH1D* h_data = (TH1D*)file[2]->Get(temp.Data());
    TH1D* h_mc   = (TH1D*)file[0]->Get(temp.Data());

    h_data->Rebin(rebin);
    h_mc->Rebin(rebin);

    if( temp.Contains("_numPV") ){
      h_mc->Scale( h_data->Integral() / h_mc->Integral() );
    }


    h_data->SetLineColor(color[0]);
    h_mc->SetLineColor(color[2]);

    h_data->SetMarkerColor(color[0]);
    h_mc->SetMarkerColor(color[2]);

    h_data->SetMarkerStyle(20);
    h_mc->SetMarkerStyle(20);

    legend->AddEntry(h_data,histLabels[2],"p");
    legend->AddEntry(h_mc,histLabels[0],"p");

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


    // HLT
    myC1->cd(1);
    h_data->Draw("pe1");
    h_mc->Draw("pe1same");
    legend->Draw();
    LumiInfoLatex.Draw();
    CMSInfoLatex.Draw();
    PublishInfoLatex.Draw();

    myC1->cd(2);
    myRatio->SetLineWidth(2);
    myRatio->Draw("pe1");
    myLine->Draw();

    plotname = dirprefix + temp_mh + "_data2mc_lin.png";
    myC1->Print(plotname);

    plotname = dirprefix + temp_mh + "_data2mc_lin.pdf";
    if( printPDF_ ) myC1->Print(plotname);


    // log
    h_data->GetYaxis()->SetRangeUser(0.4,12 * max_content);
    if( temp.Contains("_selection") ){
      h_data->GetYaxis()->SetRangeUser(400,12 * max_content);
    }

    myC1->cd(1);
    gPad->SetLogy(1);

    plotname = dirprefix + temp_mh + "_data2mc_log.png";
    myC1->Print(plotname);

    plotname = dirprefix + temp_mh + "_data2mc_log.pdf";
    if( printPDF_ ) myC1->Print(plotname);

    gPad->SetLogy(0);



    delete myRatio;
    delete myLine;
    delete legend;
  } // end loop on hists

  ////////////////////////////////////////////////////////////////////////////


  // Close the file
  for( int iFile=0; iFile<NumSamples; iFile++ ) file[iFile]->Close();
  std::cout << " Done! " << std::endl;
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


TEfficiency* makeEfficiency( TH1D* h_numerator, TH1D* h_denominator ){

  int nbins = h_numerator->GetNbinsX();
  for( int iBin=0; iBin<nbins; iBin++ ){
    double numerator = h_numerator->GetBinContent(iBin+1);
    double denominator = h_denominator->GetBinContent(iBin+1);

    if( numerator>denominator ) h_numerator->SetBinContent(iBin+1, denominator-0.0000001);
  }

  TEfficiency* result = new TEfficiency(*h_numerator, *h_denominator);

  return result;
 }

/*

h_HT30_L1HTT100_passHLTEle27HT200->Rebin(25);
h_HT30->Rebin(25);

TH1D* h_ratio = (TH1D*)h_HT30_L1HTT100_passHLTEle27HT200->Clone("h_ratio");

h_ratio->Divide(h_HT30_L1HTT100_passHLTEle27HT200,h_HT30,1,1,"B");


TF1 *f1 = new TF1("f1", "[0]*(1./(1. + [2]*exp(-[1]*x)))", 100, 500);
f1->SetParameters(1,1,1);

h_ratio->Fit("f1","R");



h_HT30_L1HTT125_passHLTEle27HT200->Rebin(25);
h_HT30->Rebin(25);

TH1D* h_ratio = (TH1D*)h_HT30_L1HTT125_passHLTEle27HT200->Clone("h_ratio");

h_ratio->Divide(h_HT30_L1HTT125_passHLTEle27HT200,h_HT30,1,1,"B");


TF1 *f1 = new TF1("f1", "[0]*(1./(1. + [2]*exp(-[1]*x)))", 100, 500);
f1->SetParameters(1,1,1);

h_ratio->Fit("f1","R");

 */
