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

void dump_CMVA_correction_factors(TString inputFileName  = "file.root") {

  TFile *histFile = TFile::Open(inputFileName);


  
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


  
  TH1D* h_perjet_wgt_cmva_Sys[5][3][numCSVSys];
  TH1D* h_perjet_wgt_cmva_hf_Sys[5][1][numCSVSys];
  TH1D* h_perjet_wgt_cmva_lf_Sys[4][3][numCSVSys];
  TH1D* h_perjet_wgt_cmva_cf_Sys[5][1][numCSVSys];

  
  for( int iSys=0; iSys<numCSVSys; iSys++ ){
    int useSys = csv_systematics[iSys];

    for( int iPt=0; iPt < 5; iPt++ ){
      for( int iEta=0; iEta < 3; iEta++ ){
	h_perjet_wgt_cmva_Sys[iPt][iEta][iSys] = (TH1D*)histFile->Get(Form("h_perjet_wgt_cmva_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys));
	if( iEta<1 ) h_perjet_wgt_cmva_hf_Sys[iPt][iEta][iSys] = (TH1D*)histFile->Get(Form("h_perjet_wgt_cmva_hf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys));
	if( iEta<1 ) h_perjet_wgt_cmva_cf_Sys[iPt][iEta][iSys] = (TH1D*)histFile->Get(Form("h_perjet_wgt_cmva_cf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys));
	if( iPt<4  ) h_perjet_wgt_cmva_lf_Sys[iPt][iEta][iSys] = (TH1D*)histFile->Get(Form("h_perjet_wgt_cmva_lf_iPt%d_iEta%d_Sys%d",iPt,iEta,useSys));
      }
    }
  }

  
  std::cout << "***********************************************************" << std::endl;
  std::cout << "HF correction factors!" << std::endl;
  std::cout << " " << std::endl;

  for( int iPt=0; iPt < 5; iPt++ ){	   
    printf("    if( hist_name[iHist]==\"csv_ratio_Pt%d_Eta%d\" ){ \n", iPt, 0);
    printf("         h_csv_ratio_final[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][0]->GetMean() );
    printf("         h_csv_ratio_final_JESUp[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][1]->GetMean() );
    printf("         h_csv_ratio_final_JESDown[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][2]->GetMean() );
    printf("         h_csv_ratio_final_LFUp[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][3]->GetMean() );
    printf("         h_csv_ratio_final_LFDown[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][4]->GetMean() );
    printf("         h_csv_ratio_final_Stats1Up[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][7]->GetMean() );
    printf("         h_csv_ratio_final_Stats1Down[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][8]->GetMean() );
    printf("         h_csv_ratio_final_Stats2Up[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][9]->GetMean() );
    printf("         h_csv_ratio_final_Stats2Down[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][10]->GetMean() );
    std::cout << " " << std::endl;
    printf("         h_csv_ratio[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][0]->GetMean() );
    printf("         h_csv_ratio_JESUp[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][1]->GetMean() );
    printf("         h_csv_ratio_JESDown[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][2]->GetMean() );
    printf("         h_csv_ratio_LFUp[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][3]->GetMean() );
    printf("         h_csv_ratio_LFDown[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][4]->GetMean() );
    printf("         h_csv_ratio_Stats1Up[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][7]->GetMean() );
    printf("         h_csv_ratio_Stats1Down[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][8]->GetMean() );
    printf("         h_csv_ratio_Stats2Up[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][9]->GetMean() );
    printf("         h_csv_ratio_Stats2Down[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_hf_Sys[iPt][0][10]->GetMean() );
    printf("    } \n");
  }

  std::cout << " " << std::endl;
  std::cout << "***********************************************************" << std::endl;
  std::cout << "LF correction factors!" << std::endl;
  std::cout << " " << std::endl;

  for( int iPt=0; iPt < 4; iPt++ ){	   
    for( int iEta=0; iEta < 3; iEta++ ){
      printf("    if( hist_name[iHist]==\"csv_ratio_Pt%d_Eta%d\" ){ \n", iPt, iEta);
      printf("         h_csv_ratio_final[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][0]->GetMean() );
      printf("         h_csv_ratio_final_JESUp[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][1]->GetMean() );
      printf("         h_csv_ratio_final_JESDown[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][2]->GetMean() );
      printf("         h_csv_ratio_final_HFUp[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][5]->GetMean() );
      printf("         h_csv_ratio_final_HFDown[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][6]->GetMean() );
      printf("         h_csv_ratio_final_Stats1Up[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][11]->GetMean() );
      printf("         h_csv_ratio_final_Stats1Down[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][12]->GetMean() );
      printf("         h_csv_ratio_final_Stats2Up[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][13]->GetMean() );
      printf("         h_csv_ratio_final_Stats2Down[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][14]->GetMean() );
      printf("    } \n");
      std::cout << " " << std::endl;
      printf("    if( hist_name[iHist]==\"csv_ratio_Pt%d_Eta%d\" ){ \n", iPt, iEta);
      printf("         h_csv_ratio[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][0]->GetMean() );
      printf("         h_csv_ratio_JESUp[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][1]->GetMean() );
      printf("         h_csv_ratio_JESDown[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][2]->GetMean() );
      printf("         h_csv_ratio_HFUp[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][5]->GetMean() );
      printf("         h_csv_ratio_HFDown[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][6]->GetMean() );
      printf("         h_csv_ratio_Stats1Up[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][11]->GetMean() );
      printf("         h_csv_ratio_Stats1Down[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][12]->GetMean() );
      printf("         h_csv_ratio_Stats2Up[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][13]->GetMean() );
      printf("         h_csv_ratio_Stats2Down[iHist]->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_lf_Sys[iPt][iEta][14]->GetMean() );
      printf("    } \n");
    }
  }


   
  std::cout << "***********************************************************" << std::endl;
  std::cout << "CF correction factors!" << std::endl;
  std::cout << " " << std::endl;

  for( int iPt=0; iPt < 5; iPt++ ){	   
    printf("    if( hist_name[iHist]==\"csv_ratio_Pt%d_Eta%d\" ){ \n", iPt, 0);
    printf("         h_cErr1Up->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_cf_Sys[iPt][0][15]->GetMean() );
    printf("         h_cErr1Down->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_cf_Sys[iPt][0][16]->GetMean() );
    printf("         h_cErr2Up->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_cf_Sys[iPt][0][17]->GetMean() );
    printf("         h_cErr2Down->Scale( %.6f ); \n", 1./h_perjet_wgt_cmva_cf_Sys[iPt][0][18]->GetMean() );
    printf("    } \n");
  }

  std::cout << " " << std::endl;
  std::cout << "***********************************************************" << std::endl;


      
  std::cout << "Done." << std::endl;

}
