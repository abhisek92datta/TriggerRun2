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

//*****************************************************************************


void makePlots_hltEleHT_Zcr( bool printPDF_ = false, double sfMC_ = 1., bool renorm_ = false ){

  TH1::SetDefaultSumw2();

  TFile* file_data = new TFile("HistoFiles/hltEleHT_treeReader_Zcr_SingleElectron_Run2015D_PromptReco_254231_258158_histo.root");
  TFile* file_mc   = new TFile("HistoFiles/hltEleHT_treeReader_Zcr_DYJets_M50_13TeV_Spring15_Asympt25ns_histo.root");

  std::vector<TString> histLabels(2);
  histLabels[0] = "SingleElectron";
  histLabels[1] = "DYJetsToLL_M50";

  Color_t color[8];
  color[0] = kBlue;
  color[1] = kRed;

  // color[2] = kBlack;
  // color[3] = kGreen+1;
  // color[4] = kMagenta+2;
  // color[5] = kRed+1;
  // color[6] = kGreen-5;
  // color[7] = kRed+3;


  TString dirprefix = "Images/Images_2015_11_06_hltEleHT_Zcr/";
  if( renorm_ ) dirprefix = "Images/Images_2015_11_06_hltEleHT_Zcr_renorm/";

  struct stat st;
  if( stat(dirprefix.Data(),&st) != 0 )  mkdir(dirprefix.Data(),0777);


 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<std::string> histoname1;
  std::vector<std::string> histoname2;

  histoname1.push_back("h_probe_pt");
  histoname1.push_back("h_probe_eta");
  histoname1.push_back("h_probe_phi");
  histoname1.push_back("h_probe_numPVs");
  histoname1.push_back("h_HT30");
  histoname1.push_back("h_HT30_v2");

  histoname2.push_back("h_numPVs");
  histoname2.push_back("h_numPVs_PUwgt");
  histoname2.push_back("h_event_selection");
  histoname2.push_back("h_diele_mass");
  histoname2.push_back("h_diele_mass_closestZmass");
  histoname2.push_back("h_diele_deltaR");
  histoname2.push_back("h_L1HTT");
  histoname2.push_back("h_HT30");
  histoname2.push_back("h_HT30_L1HTT125");
  histoname2.push_back("h_HT30_L1HTT125_passHLTEle27HT200");

  histoname2.push_back("h_met_pt");

  histoname2.push_back("h_jet_1_pt");
  histoname2.push_back("h_jet_2_pt");
  histoname2.push_back("h_jet_3_pt");
  histoname2.push_back("h_jet_4_pt");

  histoname2.push_back("h_jet_pt");
  histoname2.push_back("h_jet_eta");
  histoname2.push_back("h_jet_phi");
  histoname2.push_back("h_jet_csv");

  histoname2.push_back("h_numJet");
  histoname2.push_back("h_numBtag");

  histoname2.push_back("h_deltaR_jet_ele");

  histoname2.push_back("h_numL1EG25");
  histoname2.push_back("h_l1EG25_pt");
  histoname2.push_back("h_l1EG25_eta");
  histoname2.push_back("h_l1EG25_phi");
  histoname2.push_back("h_l1EG25_1_pt");
  histoname2.push_back("h_l1EG25_2_pt");
  histoname2.push_back("h_l1EG25_1_eta");
  histoname2.push_back("h_l1EG25_2_eta");


 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TGaxis::SetMaxDigits(3);

  TString lumiinfo = "553 pb^{-1} (13 TeV)";
  TLatex LumiInfoLatex(0.70, 0.91, lumiinfo);
  LumiInfoLatex.SetNDC(); LumiInfoLatex.SetTextFont(42);
  LumiInfoLatex.SetTextSize(0.04);

  //TString cmsinfo =   "CMS Preliminary";
  TString cmsinfo =   "CMS";
  TLatex CMSInfoLatex(0.18, 0.91, cmsinfo);
  CMSInfoLatex.SetNDC(); CMSInfoLatex.SetTextFont(42);
  CMSInfoLatex.SetTextFont(61);
  CMSInfoLatex.SetTextSize(0.055); //SBOUTLE

  std::string publishinfo =   "Preliminary"; //DPUIGH
  TLatex PublishInfoLatex(0.28, 0.91, publishinfo.c_str()); //SBOUTLE
  PublishInfoLatex.SetNDC();
  PublishInfoLatex.SetTextFont(52);
  PublishInfoLatex.SetTextSize(0.045); //SBOUTLE


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

    if( temp.Contains("_HT30") ) rebin = 10;

    TString s_all = temp+"_all";
    TString s_l1t = temp+"_pL1T";
    TString s_hlt = temp+"_pHLT";

    if( temp.Contains("_HT30") ){
      s_all = temp+"";
      s_l1t = temp+"_L1HTT125";
      s_hlt = temp+"_L1HTT125_passHLTEle27HT200";
    }
    if( temp.Contains("_HT30_v2") ){
      s_all = temp+"";
      s_l1t = temp+"_L1HTT";
      s_hlt = temp+"_L1HTT_HLTHT";
    }

    TH1D* h_data_all = (TH1D*)file_data->Get(s_all);
    TH1D* h_data_l1t = (TH1D*)file_data->Get(s_l1t);
    TH1D* h_data_hlt = (TH1D*)file_data->Get(s_hlt);

    TH1D* h_mc_all = (TH1D*)file_mc->Get(s_all);
    TH1D* h_mc_l1t = (TH1D*)file_mc->Get(s_l1t);
    TH1D* h_mc_hlt = (TH1D*)file_mc->Get(s_hlt);


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

    h_eff_mc_l1t_all->SetLineColor(color[1]);
    h_eff_mc_l1t_all->SetMarkerColor(color[1]);
    h_eff_mc_l1t_all->SetMarkerStyle(20);

    h_eff_mc_hlt_all->SetLineColor(color[1]);
    h_eff_mc_hlt_all->SetMarkerColor(color[1]);
    h_eff_mc_hlt_all->SetMarkerStyle(20);

    h_eff_mc_hlt_l1t->SetLineColor(color[1]);
    h_eff_mc_hlt_l1t->SetMarkerColor(color[1]);
    h_eff_mc_hlt_l1t->SetMarkerStyle(20);


    legend->AddEntry(h_eff_data_l1t_all,histLabels[0],"p");
    legend->AddEntry(h_eff_mc_l1t_all,histLabels[1],"p");


    TH1D* h_ratio_eff_l1t_all = divide_ratio_plots(h_data_l1t, h_data_all, h_mc_l1t, h_mc_all);
    TH1D* h_ratio_eff_hlt_all = divide_ratio_plots(h_data_hlt, h_data_all, h_mc_hlt, h_mc_all);
    TH1D* h_ratio_eff_hlt_l1t = divide_ratio_plots(h_data_hlt, h_data_l1t, h_mc_hlt, h_mc_l1t);

    h_ratio_eff_l1t_all->SetMarkerStyle(20);
    h_ratio_eff_hlt_all->SetMarkerStyle(20);
    h_ratio_eff_hlt_l1t->SetMarkerStyle(20);


    h_data_all->SetStats(0);

    //c3->SetTopMargin(.05);
    //c1->SetRightMargin(.05);

    double ratioMax = 1.6;
    double ratioMin = 0.5;

    double xmin = h_data_all->GetBinLowEdge(1);
    double xmax = h_data_all->GetBinLowEdge(nbins) + h_data_all->GetBinWidth(nbins);


    TH1D* myRatio = (TH1D*)h_data_all->Clone("ratio_"+temp);
    //new TH1D("ratio", "", nbins, xmin, xmax );

    myRatio->SetStats(0);
    myRatio->Sumw2();
    myRatio->SetLineColor(kBlack);
    myRatio->SetMarkerColor(kBlack);
    //myRatio->Divide(hist[bin_one],hist[bin_two]);

    myRatio->SetMinimum(ratioMin);
    myRatio->SetMaximum(ratioMax);
    //myRatio->GetYaxis()->SetNdivisions(50000+404);
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


    myRatio->GetYaxis()->CenterTitle(kTRUE);


    double maxPt = 499.9999;
    double minPt = 0.000001;

    h_data_all->GetYaxis()->SetTitleOffset(1.0);
    h_data_all->GetYaxis()->SetTitleSize(0.05);

    h_data_all->GetYaxis()->SetTitle("Efficiency");

    h_data_all->GetYaxis()->SetRangeUser(0.,1.15);


    bool rerange = false;
    //if( (temp.Contains("_pt")) ){ rerange=true; xmin = 10.0001; xmax = 199.9; }
    if( (temp.Contains("_pt")) ){ rerange=true; xmin = 30.001; xmax = 199.9; }

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


  ////////////////////////////////////////////////////////////////////////////


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


    TH1D* h_data = (TH1D*)file_data->Get(temp.Data());
    TH1D* h_mc   = (TH1D*)file_mc->Get(temp.Data());

    h_data->Rebin(rebin);
    h_mc->Rebin(rebin);

    h_mc->Scale(sfMC_);

    if( temp.Contains("_numPV") || renorm_ ){
      h_mc->Scale( h_data->Integral() / h_mc->Integral() );
    }


    h_data->SetLineColor(color[0]);
    h_mc->SetLineColor(color[1]);

    h_data->SetMarkerColor(color[0]);
    h_mc->SetMarkerColor(color[1]);

    h_data->SetMarkerStyle(20);
    h_mc->SetMarkerStyle(20);

    legend->AddEntry(h_data,histLabels[0],"pe1");
    legend->AddEntry(h_mc,histLabels[1],"pe1");

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
      h_data->GetYaxis()->SetRangeUser(1000,12 * max_content);
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



  // Close the file
  file_data->Close();
  file_mc->Close();
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
