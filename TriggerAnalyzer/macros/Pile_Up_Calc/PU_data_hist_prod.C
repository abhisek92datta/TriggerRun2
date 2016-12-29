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


void PU_data_hist_prod()
{
	TFile *f = new TFile("output.root");
	TH1F *h = (TH1F*)f->Get("pileup");
	ofstream fout;
	fout.open("PU_Data.txt");

	for(int i=1; i<=75; i++)
	{
		fout<<h->GetBinContent(i)<<"\n";
	}
	fout.close();
	delete f;
	return;
}
