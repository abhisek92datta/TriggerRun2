//Root includes                                   
#include "TROOT.h"
#include "Riostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TAxis.h"
#include "TKey.h"
#include "TList.h"

void dump_cumulative_btagSFs(bool isCSV = true, TString inputFileName  = "file.root") {

  TFile *histFile = TFile::Open(inputFileName);

  TString btagger = ( isCSV ) ? "CSV" : "CMVA";

  
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


  
  
  TH1D* h_numEvents_jet30_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_jet30_nowgtCSV_WP[numCSVSys];

  TH1D* h_numEvents_jet20_wgtCSV_WP[numCSVSys];
  TH1D* h_numEvents_jet20_nowgtCSV_WP[numCSVSys];

  TH1D* h_numEvents_jet30to60_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet30to60_nowgtCMVA_WP[numCSVSys];

  TH1D* h_numEvents_jet60to120_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet60to120_nowgtCMVA_WP[numCSVSys];

  TH1D* h_numEvents_jet120toInf_wgtCMVA_WP[numCSVSys];
  TH1D* h_numEvents_jet120toInf_nowgtCMVA_WP[numCSVSys];

  double sf_jet30_csvL[numCSVSys];
  double sf_jet30_csvM[numCSVSys];
  double sf_jet30_csvT[numCSVSys];
  
  double sf_jet20_csvL[numCSVSys];
  double sf_jet20_csvM[numCSVSys];
  double sf_jet20_csvT[numCSVSys];

  double sf_jet30to60_cmvaL[numCSVSys];
  double sf_jet30to60_cmvaM[numCSVSys];
  double sf_jet30to60_cmvaT[numCSVSys];
  
  double sf_jet60to120_cmvaL[numCSVSys];
  double sf_jet60to120_cmvaM[numCSVSys];
  double sf_jet60to120_cmvaT[numCSVSys];
  
  double sf_jet120toInf_cmvaL[numCSVSys];
  double sf_jet120toInf_cmvaM[numCSVSys];
  double sf_jet120toInf_cmvaT[numCSVSys];

  double sf_jet30_csvL_err_up2 = 0;
  double sf_jet30_csvM_err_up2 = 0;
  double sf_jet30_csvT_err_up2 = 0;

  double sf_jet30_csvL_err_down2 = 0;
  double sf_jet30_csvM_err_down2 = 0;
  double sf_jet30_csvT_err_down2 = 0;


  double sf_jet20_csvL_err_up2 = 0;
  double sf_jet20_csvM_err_up2 = 0;
  double sf_jet20_csvT_err_up2 = 0;
  
  double sf_jet20_csvL_err_down2 = 0;
  double sf_jet20_csvM_err_down2 = 0;
  double sf_jet20_csvT_err_down2 = 0;




  double sf_jet30to60_cmvaL_err_up2 = 0;
  double sf_jet30to60_cmvaM_err_up2 = 0;
  double sf_jet30to60_cmvaT_err_up2 = 0;

  double sf_jet30to60_cmvaL_err_down2 = 0;
  double sf_jet30to60_cmvaM_err_down2 = 0;
  double sf_jet30to60_cmvaT_err_down2 = 0;


  double sf_jet60to120_cmvaL_err_up2 = 0;
  double sf_jet60to120_cmvaM_err_up2 = 0;
  double sf_jet60to120_cmvaT_err_up2 = 0;

  double sf_jet60to120_cmvaL_err_down2 = 0;
  double sf_jet60to120_cmvaM_err_down2 = 0;
  double sf_jet60to120_cmvaT_err_down2 = 0;


  double sf_jet120toInf_cmvaL_err_up2 = 0;
  double sf_jet120toInf_cmvaM_err_up2 = 0;
  double sf_jet120toInf_cmvaT_err_up2 = 0;

  double sf_jet120toInf_cmvaL_err_down2 = 0;
  double sf_jet120toInf_cmvaM_err_down2 = 0;
  double sf_jet120toInf_cmvaT_err_down2 = 0;

  
  for( int iSys=0; iSys<numCSVSys; iSys++ ){
    int useSys = csv_systematics[iSys];

    h_numEvents_jet30_wgtCSV_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet30_wgt%s_WP_iSys_%d",btagger.Data(),useSys));
    h_numEvents_jet30_nowgtCSV_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet30_nowgt%s_WP_iSys_%d",btagger.Data(),useSys));

    h_numEvents_jet20_wgtCSV_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet20_wgt%s_WP_iSys_%d",btagger.Data(),useSys));
    h_numEvents_jet20_nowgtCSV_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet20_nowgt%s_WP_iSys_%d",btagger.Data(),useSys));


    h_numEvents_jet30to60_wgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet30to60_wgt%s_WP_iSys_%d","CMVA",useSys));
    h_numEvents_jet30to60_nowgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet30to60_nowgt%s_WP_iSys_%d","CMVA",useSys));

    h_numEvents_jet60to120_wgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet60to120_wgt%s_WP_iSys_%d","CMVA",useSys));
    h_numEvents_jet60to120_nowgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet60to120_nowgt%s_WP_iSys_%d","CMVA",useSys));

    h_numEvents_jet120toInf_wgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet120toInf_wgt%s_WP_iSys_%d","CMVA",useSys));
    h_numEvents_jet120toInf_nowgtCMVA_WP[iSys] = (TH1D*)histFile->Get(Form("h_numEvents_jet120toInf_nowgt%s_WP_iSys_%d","CMVA",useSys));




    h_numEvents_jet30_wgtCSV_WP[iSys]->Divide(h_numEvents_jet30_nowgtCSV_WP[iSys]);
    h_numEvents_jet20_wgtCSV_WP[iSys]->Divide(h_numEvents_jet20_nowgtCSV_WP[iSys]);

    
    h_numEvents_jet30to60_wgtCMVA_WP[iSys]->Divide(h_numEvents_jet30to60_nowgtCMVA_WP[iSys]);
    h_numEvents_jet60to120_wgtCMVA_WP[iSys]->Divide(h_numEvents_jet60to120_nowgtCMVA_WP[iSys]);
    h_numEvents_jet120toInf_wgtCMVA_WP[iSys]->Divide(h_numEvents_jet120toInf_nowgtCMVA_WP[iSys]);

    if( iSys==0 ){
      sf_jet30_csvL[iSys] = h_numEvents_jet30_wgtCSV_WP[iSys]->GetBinContent(1);
      sf_jet30_csvM[iSys] = h_numEvents_jet30_wgtCSV_WP[iSys]->GetBinContent(2);
      sf_jet30_csvT[iSys] = h_numEvents_jet30_wgtCSV_WP[iSys]->GetBinContent(3);

      sf_jet20_csvL[iSys] = h_numEvents_jet20_wgtCSV_WP[iSys]->GetBinContent(1);
      sf_jet20_csvM[iSys] = h_numEvents_jet20_wgtCSV_WP[iSys]->GetBinContent(2);
      sf_jet20_csvT[iSys] = h_numEvents_jet20_wgtCSV_WP[iSys]->GetBinContent(3);


      sf_jet30to60_cmvaL[iSys] = h_numEvents_jet30to60_wgtCMVA_WP[iSys]->GetBinContent(1);
      sf_jet30to60_cmvaM[iSys] = h_numEvents_jet30to60_wgtCMVA_WP[iSys]->GetBinContent(2);
      sf_jet30to60_cmvaT[iSys] = h_numEvents_jet30to60_wgtCMVA_WP[iSys]->GetBinContent(3);

      sf_jet60to120_cmvaL[iSys] = h_numEvents_jet60to120_wgtCMVA_WP[iSys]->GetBinContent(1);
      sf_jet60to120_cmvaM[iSys] = h_numEvents_jet60to120_wgtCMVA_WP[iSys]->GetBinContent(2);
      sf_jet60to120_cmvaT[iSys] = h_numEvents_jet60to120_wgtCMVA_WP[iSys]->GetBinContent(3);

      sf_jet120toInf_cmvaL[iSys] = h_numEvents_jet120toInf_wgtCMVA_WP[iSys]->GetBinContent(1);
      sf_jet120toInf_cmvaM[iSys] = h_numEvents_jet120toInf_wgtCMVA_WP[iSys]->GetBinContent(2);
      sf_jet120toInf_cmvaT[iSys] = h_numEvents_jet120toInf_wgtCMVA_WP[iSys]->GetBinContent(3);

    }
    if( iSys>0 ){
      sf_jet30_csvL[iSys] = h_numEvents_jet30_wgtCSV_WP[iSys]->GetBinContent(1) - sf_jet30_csvL[0];
      sf_jet30_csvM[iSys] = h_numEvents_jet30_wgtCSV_WP[iSys]->GetBinContent(2) - sf_jet30_csvM[0];
      sf_jet30_csvT[iSys] = h_numEvents_jet30_wgtCSV_WP[iSys]->GetBinContent(3) - sf_jet30_csvT[0];

      sf_jet20_csvL[iSys] = h_numEvents_jet20_wgtCSV_WP[iSys]->GetBinContent(1) - sf_jet20_csvL[0];
      sf_jet20_csvM[iSys] = h_numEvents_jet20_wgtCSV_WP[iSys]->GetBinContent(2) - sf_jet20_csvM[0];
      sf_jet20_csvT[iSys] = h_numEvents_jet20_wgtCSV_WP[iSys]->GetBinContent(3) - sf_jet20_csvT[0];



      sf_jet30to60_cmvaL[iSys] = h_numEvents_jet30to60_wgtCMVA_WP[iSys]->GetBinContent(1) - sf_jet30to60_cmvaL[0];
      sf_jet30to60_cmvaM[iSys] = h_numEvents_jet30to60_wgtCMVA_WP[iSys]->GetBinContent(2) - sf_jet30to60_cmvaM[0];
      sf_jet30to60_cmvaT[iSys] = h_numEvents_jet30to60_wgtCMVA_WP[iSys]->GetBinContent(3) - sf_jet30to60_cmvaT[0];

      sf_jet60to120_cmvaL[iSys] = h_numEvents_jet60to120_wgtCMVA_WP[iSys]->GetBinContent(1) - sf_jet60to120_cmvaL[0];
      sf_jet60to120_cmvaM[iSys] = h_numEvents_jet60to120_wgtCMVA_WP[iSys]->GetBinContent(2) - sf_jet60to120_cmvaM[0];
      sf_jet60to120_cmvaT[iSys] = h_numEvents_jet60to120_wgtCMVA_WP[iSys]->GetBinContent(3) - sf_jet60to120_cmvaT[0];

      sf_jet120toInf_cmvaL[iSys] = h_numEvents_jet120toInf_wgtCMVA_WP[iSys]->GetBinContent(1) - sf_jet120toInf_cmvaL[0];
      sf_jet120toInf_cmvaM[iSys] = h_numEvents_jet120toInf_wgtCMVA_WP[iSys]->GetBinContent(2) - sf_jet120toInf_cmvaM[0];
      sf_jet120toInf_cmvaT[iSys] = h_numEvents_jet120toInf_wgtCMVA_WP[iSys]->GetBinContent(3) - sf_jet120toInf_cmvaT[0];

      
      if( sf_jet30_csvL[iSys] > 0 ) sf_jet30_csvL_err_up2   += sf_jet30_csvL[iSys] * sf_jet30_csvL[iSys];
      else                          sf_jet30_csvL_err_down2 += sf_jet30_csvL[iSys] * sf_jet30_csvL[iSys];

      if( sf_jet30_csvM[iSys] > 0 ) sf_jet30_csvM_err_up2   += sf_jet30_csvM[iSys] * sf_jet30_csvM[iSys];
      else                          sf_jet30_csvM_err_down2 += sf_jet30_csvM[iSys] * sf_jet30_csvM[iSys];

      if( sf_jet30_csvT[iSys] > 0 ) sf_jet30_csvT_err_up2   += sf_jet30_csvT[iSys] * sf_jet30_csvT[iSys];
      else                          sf_jet30_csvT_err_down2 += sf_jet30_csvT[iSys] * sf_jet30_csvT[iSys];

      
      if( sf_jet20_csvL[iSys] > 0 ) sf_jet20_csvL_err_up2   += sf_jet20_csvL[iSys] * sf_jet20_csvL[iSys];
      else                          sf_jet20_csvL_err_down2 += sf_jet20_csvL[iSys] * sf_jet20_csvL[iSys];

      if( sf_jet20_csvM[iSys] > 0 ) sf_jet20_csvM_err_up2   += sf_jet20_csvM[iSys] * sf_jet20_csvM[iSys];
      else                          sf_jet20_csvM_err_down2 += sf_jet20_csvM[iSys] * sf_jet20_csvM[iSys];

      if( sf_jet20_csvT[iSys] > 0 ) sf_jet20_csvT_err_up2   += sf_jet20_csvT[iSys] * sf_jet20_csvT[iSys];
      else                          sf_jet20_csvT_err_down2 += sf_jet20_csvT[iSys] * sf_jet20_csvT[iSys];

      

      
      if( sf_jet30to60_cmvaL[iSys] > 0 ) sf_jet30to60_cmvaL_err_up2   += sf_jet30to60_cmvaL[iSys] * sf_jet30to60_cmvaL[iSys];
      else                               sf_jet30to60_cmvaL_err_down2 += sf_jet30to60_cmvaL[iSys] * sf_jet30to60_cmvaL[iSys];

      if( sf_jet30to60_cmvaM[iSys] > 0 ) sf_jet30to60_cmvaM_err_up2   += sf_jet30to60_cmvaM[iSys] * sf_jet30to60_cmvaM[iSys];
      else                               sf_jet30to60_cmvaM_err_down2 += sf_jet30to60_cmvaM[iSys] * sf_jet30to60_cmvaM[iSys];

      if( sf_jet30to60_cmvaT[iSys] > 0 ) sf_jet30to60_cmvaT_err_up2   += sf_jet30to60_cmvaT[iSys] * sf_jet30to60_cmvaT[iSys];
      else                               sf_jet30to60_cmvaT_err_down2 += sf_jet30to60_cmvaT[iSys] * sf_jet30to60_cmvaT[iSys];

      
      if( sf_jet60to120_cmvaL[iSys] > 0 ) sf_jet60to120_cmvaL_err_up2   += sf_jet60to120_cmvaL[iSys] * sf_jet60to120_cmvaL[iSys];
      else                                sf_jet60to120_cmvaL_err_down2 += sf_jet60to120_cmvaL[iSys] * sf_jet60to120_cmvaL[iSys];

      if( sf_jet60to120_cmvaM[iSys] > 0 ) sf_jet60to120_cmvaM_err_up2   += sf_jet60to120_cmvaM[iSys] * sf_jet60to120_cmvaM[iSys];
      else                                sf_jet60to120_cmvaM_err_down2 += sf_jet60to120_cmvaM[iSys] * sf_jet60to120_cmvaM[iSys];

      if( sf_jet60to120_cmvaT[iSys] > 0 ) sf_jet60to120_cmvaT_err_up2   += sf_jet60to120_cmvaT[iSys] * sf_jet60to120_cmvaT[iSys];
      else                                sf_jet60to120_cmvaT_err_down2 += sf_jet60to120_cmvaT[iSys] * sf_jet60to120_cmvaT[iSys];

      
      if( sf_jet120toInf_cmvaL[iSys] > 0 ) sf_jet120toInf_cmvaL_err_up2   += sf_jet120toInf_cmvaL[iSys] * sf_jet120toInf_cmvaL[iSys];
      else                                 sf_jet120toInf_cmvaL_err_down2 += sf_jet120toInf_cmvaL[iSys] * sf_jet120toInf_cmvaL[iSys];

      if( sf_jet120toInf_cmvaM[iSys] > 0 ) sf_jet120toInf_cmvaM_err_up2   += sf_jet120toInf_cmvaM[iSys] * sf_jet120toInf_cmvaM[iSys];
      else                                 sf_jet120toInf_cmvaM_err_down2 += sf_jet120toInf_cmvaM[iSys] * sf_jet120toInf_cmvaM[iSys];

      if( sf_jet120toInf_cmvaT[iSys] > 0 ) sf_jet120toInf_cmvaT_err_up2   += sf_jet120toInf_cmvaT[iSys] * sf_jet120toInf_cmvaT[iSys];
      else                                 sf_jet120toInf_cmvaT_err_down2 += sf_jet120toInf_cmvaT[iSys] * sf_jet120toInf_cmvaT[iSys];

    }
  }


  
  std::cout << "***********************************************************" << std::endl;
  std::cout << "  Cumulative " << btagger << "v2 scale factors for jet pT > 30 " << std::endl;
  std::cout << " " << std::endl;

  printf("  %sv2L SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet30_csvL[0], sqrt(sf_jet30_csvL_err_up2), sqrt(sf_jet30_csvL_err_down2), sf_jet30_csvL[1], sf_jet30_csvL[2], sf_jet30_csvL[3], sf_jet30_csvL[4],
	 sf_jet30_csvL[7], sf_jet30_csvL[8], sf_jet30_csvL[9], sf_jet30_csvL[10] );

  printf("  %sv2M SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet30_csvM[0], sqrt(sf_jet30_csvM_err_up2), sqrt(sf_jet30_csvM_err_down2), sf_jet30_csvM[1], sf_jet30_csvM[2], sf_jet30_csvM[3], sf_jet30_csvM[4],
	 sf_jet30_csvM[7], sf_jet30_csvM[8], sf_jet30_csvM[9], sf_jet30_csvM[10] );

  printf("  %sv2T SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet30_csvT[0], sqrt(sf_jet30_csvT_err_up2), sqrt(sf_jet30_csvT_err_down2), sf_jet30_csvT[1], sf_jet30_csvT[2], sf_jet30_csvT[3], sf_jet30_csvT[4],
	 sf_jet30_csvT[7], sf_jet30_csvT[8], sf_jet30_csvT[9], sf_jet30_csvT[10] );

  std::cout << " " << std::endl;
  std::cout << "***********************************************************" << std::endl;
  std::cout << "  Cumulative " << btagger << "v2 scale factors for jet pT > 20 " << std::endl;
  std::cout << " " << std::endl;

  printf("  %sv2L SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet20_csvL[0], sqrt(sf_jet20_csvL_err_up2), sqrt(sf_jet20_csvL_err_down2), sf_jet20_csvL[1], sf_jet20_csvL[2], sf_jet20_csvL[3], sf_jet20_csvL[4],
	 sf_jet20_csvL[7], sf_jet20_csvL[8], sf_jet20_csvL[9], sf_jet20_csvL[10] );

  printf("  %sv2M SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet20_csvM[0], sqrt(sf_jet20_csvM_err_up2), sqrt(sf_jet20_csvM_err_down2), sf_jet20_csvM[1], sf_jet20_csvM[2], sf_jet20_csvM[3], sf_jet20_csvM[4],
	 sf_jet20_csvM[7], sf_jet20_csvM[8], sf_jet20_csvM[9], sf_jet20_csvM[10] );

  printf("  %sv2T SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet20_csvT[0], sqrt(sf_jet20_csvT_err_up2), sqrt(sf_jet20_csvT_err_down2), sf_jet20_csvT[1], sf_jet20_csvT[2], sf_jet20_csvT[3], sf_jet20_csvT[4],
	 sf_jet20_csvT[7], sf_jet20_csvT[8], sf_jet20_csvT[9], sf_jet20_csvT[10] );

  std::cout << " " << std::endl;
  std::cout << "***********************************************************" << std::endl;
  if( !isCSV ){
  std::cout << "  Cumulative " << btagger << "v2 scale factors for jet pT [30, 60) " << std::endl;
  std::cout << " " << std::endl;

  printf("  %sv2L SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet30to60_cmvaL[0], sqrt(sf_jet30to60_cmvaL_err_up2), sqrt(sf_jet30to60_cmvaL_err_down2), sf_jet30to60_cmvaL[1], sf_jet30to60_cmvaL[2], sf_jet30to60_cmvaL[3], sf_jet30to60_cmvaL[4],
	 sf_jet30to60_cmvaL[7], sf_jet30to60_cmvaL[8], sf_jet30to60_cmvaL[9], sf_jet30to60_cmvaL[10] );

  printf("  %sv2M SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet30to60_cmvaM[0], sqrt(sf_jet30to60_cmvaM_err_up2), sqrt(sf_jet30to60_cmvaM_err_down2), sf_jet30to60_cmvaM[1], sf_jet30to60_cmvaM[2], sf_jet30to60_cmvaM[3], sf_jet30to60_cmvaM[4],
	 sf_jet30to60_cmvaM[7], sf_jet30to60_cmvaM[8], sf_jet30to60_cmvaM[9], sf_jet30to60_cmvaM[10] );

  printf("  %sv2T SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet30to60_cmvaT[0], sqrt(sf_jet30to60_cmvaT_err_up2), sqrt(sf_jet30to60_cmvaT_err_down2), sf_jet30to60_cmvaT[1], sf_jet30to60_cmvaT[2], sf_jet30to60_cmvaT[3], sf_jet30to60_cmvaT[4],
	 sf_jet30to60_cmvaT[7], sf_jet30to60_cmvaT[8], sf_jet30to60_cmvaT[9], sf_jet30to60_cmvaT[10] );

  std::cout << " " << std::endl;
  std::cout << "***********************************************************" << std::endl;
  std::cout << "  Cumulative " << btagger << "v2 scale factors for jet pT [60, 120) " << std::endl;
  std::cout << " " << std::endl;

  printf("  %sv2L SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet60to120_cmvaL[0], sqrt(sf_jet60to120_cmvaL_err_up2), sqrt(sf_jet60to120_cmvaL_err_down2), sf_jet60to120_cmvaL[1], sf_jet60to120_cmvaL[2], sf_jet60to120_cmvaL[3], sf_jet60to120_cmvaL[4],
	 sf_jet60to120_cmvaL[7], sf_jet60to120_cmvaL[8], sf_jet60to120_cmvaL[9], sf_jet60to120_cmvaL[10] );

  printf("  %sv2M SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet60to120_cmvaM[0], sqrt(sf_jet60to120_cmvaM_err_up2), sqrt(sf_jet60to120_cmvaM_err_down2), sf_jet60to120_cmvaM[1], sf_jet60to120_cmvaM[2], sf_jet60to120_cmvaM[3], sf_jet60to120_cmvaM[4],
	 sf_jet60to120_cmvaM[7], sf_jet60to120_cmvaM[8], sf_jet60to120_cmvaM[9], sf_jet60to120_cmvaM[10] );

  printf("  %sv2T SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet60to120_cmvaT[0], sqrt(sf_jet60to120_cmvaT_err_up2), sqrt(sf_jet60to120_cmvaT_err_down2), sf_jet60to120_cmvaT[1], sf_jet60to120_cmvaT[2], sf_jet60to120_cmvaT[3], sf_jet60to120_cmvaT[4],
	 sf_jet60to120_cmvaT[7], sf_jet60to120_cmvaT[8], sf_jet60to120_cmvaT[9], sf_jet60to120_cmvaT[10] );

  std::cout << " " << std::endl;
  std::cout << "***********************************************************" << std::endl;
  std::cout << "  Cumulative " << btagger << "v2 scale factors for jet pT [120, inf) " << std::endl;
  std::cout << " " << std::endl;

  printf("  %sv2L SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet120toInf_cmvaL[0], sqrt(sf_jet120toInf_cmvaL_err_up2), sqrt(sf_jet120toInf_cmvaL_err_down2), sf_jet120toInf_cmvaL[1], sf_jet120toInf_cmvaL[2], sf_jet120toInf_cmvaL[3], sf_jet120toInf_cmvaL[4],
	 sf_jet120toInf_cmvaL[7], sf_jet120toInf_cmvaL[8], sf_jet120toInf_cmvaL[9], sf_jet120toInf_cmvaL[10] );

  printf("  %sv2M SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet120toInf_cmvaM[0], sqrt(sf_jet120toInf_cmvaM_err_up2), sqrt(sf_jet120toInf_cmvaM_err_down2), sf_jet120toInf_cmvaM[1], sf_jet120toInf_cmvaM[2], sf_jet120toInf_cmvaM[3], sf_jet120toInf_cmvaM[4],
	 sf_jet120toInf_cmvaM[7], sf_jet120toInf_cmvaM[8], sf_jet120toInf_cmvaM[9], sf_jet120toInf_cmvaM[10] );

  printf("  %sv2T SF = %.4f +%.4f -%.4f : JES (%.4f, %.4f) : LF (%.4f, %.4f) : HFStats1 (%.4f, %.4f) : HFStats2 (%.4f, %.4f) \n",
	 btagger.Data(),
	 sf_jet120toInf_cmvaT[0], sqrt(sf_jet120toInf_cmvaT_err_up2), sqrt(sf_jet120toInf_cmvaT_err_down2), sf_jet120toInf_cmvaT[1], sf_jet120toInf_cmvaT[2], sf_jet120toInf_cmvaT[3], sf_jet120toInf_cmvaT[4],
	 sf_jet120toInf_cmvaT[7], sf_jet120toInf_cmvaT[8], sf_jet120toInf_cmvaT[9], sf_jet120toInf_cmvaT[10] );

  std::cout << " " << std::endl;
  std::cout << "***********************************************************" << std::endl;
  }

      
  std::cout << "Done." << std::endl;

}
