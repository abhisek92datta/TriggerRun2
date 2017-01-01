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
#include "Math/Interpolator.h"
#include "TEfficiency.h"
#include <algorithm>

#include <time.h>

#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif

#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "TriggerStudyEventVars.h"

#endif

typedef std::vector< TLorentzVector >          vecTLorentzVector;
typedef std::vector<int>                       vint;
typedef std::vector<double>                    vdouble;
typedef std::vector<std::vector<double> >      vvdouble;
typedef std::vector<std::vector<int> >         vvint;

void Luminosity_calc( int maxNentries=-1, int Njobs=1, int jobN=1 ) {

  clock_t t;
  t = clock();
  
  int N_total = 0;
  
  std::cout << "   ===> load the root files! " << std::endl;
  
  double wgt=1;
  
  std::string str_jobN;
  std::stringstream stream;
  stream << jobN;
  str_jobN = stream.str();

  std::string treefilename = "/eos/uscms/store/user/adatta/Trigger_Analysis/SingleElectron/triggerTree_SingleElectron_Run2016G/161124_192828/0001/trigger_analyzer_*.root";
  //std::string treefilename2 = "/eos/uscms/store/user/adatta/Trigger_Analysis/SingleElectron/triggerTree_SingleElectron_Run2016B/160723_151923/0001/trigger_analyzer_*.root";
  //std::string treefilename3 = "/eos/uscms/store/user/adatta/Trigger_Analysis/SingleElectron/triggerTree_SingleElectron_Run2016B/160723_151923/0002/trigger_analyzer_*.root";
  
  TChain *chain = new TChain("triggeranalyzer/triggerTree");
 
  chain->Add(treefilename.c_str());
  //chain->Add(treefilename2.c_str());
  //chain->Add(treefilename3.c_str());
   
  triggerStudyEventVars *eve=0;
  chain->SetBranchAddress("eve.", &eve );

  int nentries = chain->GetEntries();
  std::cout << "\n\t Number of entries = " << nentries << std::endl;
  std::cout << "\t Max number of entries = " << maxNentries << std::endl;
  std::cout << "\n" << std::endl;

  int use_nentries = std::max( maxNentries, nentries);

  int NeventsPerJob = int( double(use_nentries)/double(Njobs) + 0.000001 ) + 1;

  int firstEvent = (jobN-1)*NeventsPerJob + 1;
  int lastEvent  = firstEvent + NeventsPerJob;
  if( jobN==Njobs ) lastEvent = -1;
  if( jobN==1 ) firstEvent = 0;
 
  int lumi_old1=0;
  int lumi_list1[4000];
  int j1=0;
  int lumi_old2=0;
  int lumi_list2[4000];
  int j2=0;
  //int check_run = 278308;
  int check_run1 = 279841;
  int check_run2 = 279931;
 
  std::cout << "========  Starting Event Loop  ========" << std::endl;
  for (Long64_t ievt=0; ievt<chain->GetEntries();ievt++) {    
  
  		N_total++;
  		
	    if( ievt<firstEvent ) continue;
	    if( ievt==lastEvent ) break;

      	if( ievt==0 )        std::cout << "     Event " << ievt+1 << std::endl;
      	if( ievt%10000==0 && ievt!=1) std::cout << "           " << ievt << "\t" 
  					     << int(double(ievt-1)/double(nentries)*100) << "% done" << std::endl;
	  
	 	if( ievt==(maxNentries+1) && ievt!=0 ) break;
	  
	   
	  	chain->GetEntry(ievt);
	    
		int run = eve->run_;
		int lumi = eve->lumi_;
	
		//if(ievt<10)
		//if(ievt>(nentries-10))
		//	std::cout<<"Run : "<<run<<"  Lumi : "<<lumi<<"\n";	
			
		if(run==check_run1) {
			if (lumi!=lumi_old1)	{
				lumi_list1[j1++] = lumi;
				lumi_old1=lumi;
			}
		}
		
		if(run==check_run2) {
			if (lumi!=lumi_old2)	{
				lumi_list2[j2++] = lumi;
				lumi_old2=lumi;
			}
		}		
		
   }
   
   
   sort(lumi_list1, lumi_list1+j1);
   sort(lumi_list2, lumi_list2+j2);
   
   std::cout<<"Run : "<<check_run1<<"\n\n";
   //for(int i=0; i<j; i++) {
   //		std::cout<<lumi_list[i]<<"\n";
   
   // finding missing lumis
    int i=75;
    int flag;
    int count=0;
    while(i<=2122)
	{	
		if(i==399) i=407;
			
		flag=0;
		for(int k=0; k<j1; k++)
		{
			if(i==lumi_list1[k]) {
				flag=1;
				break;
			}
		}
		if(flag==0) {
			std::cout<<i<<"\n";
			count++;
		}
		i++;
	}
   	std::cout<<"\nNo of missing lumis : "<<count<<"\n\n\n";
   	
   	
  	std::cout<<"Run : "<<check_run2<<"\n\n";
    i=84;
    count=0;
    while(i<=3022)
	{	
		if(i==629) i=630;
		else if (i==744) i=746;
		else if (i==802) i=803;
		else if (i==1044) i=1045;
	
		flag=0;
		for(int k=0; k<j2; k++)
		{
			if(i==lumi_list2[k]) {
				flag=1;
				break;
			}
		}
		if(flag==0) {
			std::cout<<i<<"\n";
			count++;
		}
		i++;
	}
   	std::cout<<"\nNo of missing lumis : "<<count<<"\n\n";
   	
    
  t=clock()-t;
	
  std::cout << " Done! " <<((float)t)/CLOCKS_PER_SEC<< std::endl;
  std::cout<<"**********************************************************************************************\n";
  std::cout<<"Total No. of events : "<<N_total<<"\n";
  
}


