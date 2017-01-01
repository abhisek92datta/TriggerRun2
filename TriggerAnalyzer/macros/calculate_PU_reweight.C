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


void calculate_PU_reweight( TString dataFileName = "", TString mcFileName = "", TString dataFileName_up = "", TString dataFileName_down = "" ){

  TH1::SetDefaultSumw2();

  TFile* file_data = new TFile(dataFileName);
  TFile* file_data_up = new TFile(dataFileName_up);
  TFile* file_data_down = new TFile(dataFileName_down);
  TFile* file_mc   = new TFile(mcFileName);

  TH1D* h_data = (TH1D*)file_data->Get("pileup");
  TH1D* h_data_up = (TH1D*)file_data_up->Get("pileup");
  TH1D* h_data_down = (TH1D*)file_data_down->Get("pileup");
  TH1D* h_mc = (TH1D*)file_mc->Get("h_numTruePVs");
 
  h_data->Scale(1. / h_data->Integral());
  h_data_up->Scale(1. / h_data_up->Integral());
  h_data_down->Scale(1. / h_data_down->Integral());
  h_mc->Scale(1. / h_mc->Integral());
 
  printf("===> NOMINAL \n");
  TH1D* h_ratio = (TH1D*)h_data->Clone("h_ratio");
  h_ratio->Divide(h_mc);
 
  for( int iBin=0; iBin<h_ratio->GetNbinsX(); iBin++ ) std::cout << "  PUscale[" << iBin << "] = " << h_ratio->GetBinContent(iBin+1) << ";" << std::endl;
 

  printf("===> PU Up \n");

  TH1D* h_ratio_up = (TH1D*)h_data_up->Clone("h_ratio_up");
  h_ratio_up->Divide(h_mc);
 
  for( int iBin=0; iBin<h_ratio->GetNbinsX(); iBin++ ) std::cout << "  PUscale[" << iBin << "] = " << h_ratio_up->GetBinContent(iBin+1) << ";" << std::endl;
 


  printf("===> PU Down \n");

  TH1D* h_ratio_down = (TH1D*)h_data_down->Clone("h_ratio_down");
  h_ratio_down->Divide(h_mc);
 
  for( int iBin=0; iBin<h_ratio->GetNbinsX(); iBin++ ) std::cout << "  PUscale[" << iBin << "] = " << h_ratio_down->GetBinContent(iBin+1) << ";" << std::endl;
 

  // Close the file
  file_data->Close();
  file_data_up->Close();
  file_data_down->Close();
  file_mc->Close();
}
