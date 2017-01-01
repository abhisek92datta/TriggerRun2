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

void get_histograms(TString inputFileName = "file.root",
		    TString outputFileName = "file.root"
		    ) {

  TFile *outFile = TFile::Open(outputFileName,"RECREATE");
  TFile *inFile = TFile::Open(inputFileName);

  //Now copy all other histograms in the file as is (preserving names...)
  TList *keys = inFile->GetListOfKeys();

  TIter nextKey(keys);
  TKey *key = 0;

  while ((key = (TKey *)nextKey())) {

    TString name = key->GetName();

    if( !name.Contains("h_numTruePVs") ) continue;

    TH1 *hist = 0;
    hist = (TH1 *)inFile->Get(name);

    outFile->cd();
    hist->Write(name);
  }

  outFile->Close();

  std::cout << "Done." << std::endl;

}
